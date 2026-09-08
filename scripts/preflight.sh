#!/usr/bin/env bash
# Preflight: verify everything the pipeline needs BEFORE queuing days of compute.
#
# Stage 02 alone is 3-7 days. Discovering a missing database or an unwritable
# path after that wait is the expensive failure mode this exists to prevent.
# Read-only; changes nothing.
#
#   ./scripts/preflight.sh
#
# Exit 0 = ready to submit. Exit 1 = at least one blocker.
set -uo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

PASS=0; WARN=0; FAIL=0
ok()   { printf '  \033[32mOK\033[0m    %-38s %s\n' "$1" "${2:-}"; PASS=$((PASS+1)); }
warn() { printf '  \033[33mWARN\033[0m  %-38s %s\n' "$1" "${2:-}"; WARN=$((WARN+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %-38s %s\n' "$1" "${2:-}"; FAIL=$((FAIL+1)); }

hr() { printf '\n\033[1m%s\033[0m\n' "$1"; }

hr "Scheduler"
if command -v sbatch >/dev/null; then ok "sbatch" "$(command -v sbatch)"; else bad "sbatch" "not found"; fi
if sacctmgr -nP show assoc user="$USER" format=Account 2>/dev/null | grep -qx "$SLURM_ACCOUNT"; then
  ok "account $SLURM_ACCOUNT" "association exists"
else
  bad "account $SLURM_ACCOUNT" "no association for $USER"
fi
for p in "$PART_CPU" "$PART_GPU"; do
  [ -z "$p" ] && continue
  if sinfo -h -p "$p" -o %P 2>/dev/null | grep -q .; then ok "partition $p" "exists"
  else bad "partition $p" "not found"; fi
done

hr "Storage"
printf '  %-44s %s\n' "SCRATCH (work + databases)" "$SCRATCH"
printf '  %-44s %s\n' "PROJ (release, logs, checkpoints)" "$PROJ"
case "$SCRATCH" in
  "$PROJ"*) warn "scratch tier" "falling back to PROJ - counts against group quota" ;;
  *)        ok   "scratch tier" "$SCRATCH" ;;
esac
for d in "$SCRATCH" "$PROJ"; do
  if [ -w "$d" ] 2>/dev/null || mkdir -p "$d" 2>/dev/null; then
    avail=$(df -BG --output=avail "$d" 2>/dev/null | tail -1 | tr -dc '0-9')
    fstype=$(df -PT "$d" 2>/dev/null | tail -1 | awk '{print $2}')
    if [ -n "$avail" ] && [ "$avail" -lt 350 ]; then
      warn "space on $d" "${avail} GB free; ~300 GB needed"
    else
      ok "space on $d" "${avail:-?} GB free ($fstype)"
    fi
  else
    bad "writable $d" "cannot create or write"
  fi
done

hr "Software"
if [ -x "$MICROMAMBA_BIN" ]; then
  ok "package backend" "micromamba $("$MICROMAMBA_BIN" --version 2>/dev/null)"
elif command -v conda >/dev/null 2>&1; then
  ok "package backend" "site conda ($(command -v conda))"
elif [ -x "$STAGED/micromamba" ]; then
  warn "package backend" "micromamba staged but not installed - run 00_setup.sh"
else
  bad "package backend" "no micromamba and no conda; try: module load shared anaconda3"
fi
# Both back ends place envs at $MAMBA_ROOT_PREFIX/envs/<name>, so one test serves.
for e in prep repeats rnaseq annot ncrna func qc; do
  if [ -d "$MAMBA_ROOT_PREFIX/envs/$e" ]; then ok "env $e"
  else bad "env $e" "missing - run scripts/00_setup.sh"; fi
done
if command -v "$CONTAINER_CMD" >/dev/null; then ok "$CONTAINER_CMD"; else bad "$CONTAINER_CMD" "needed for BRAKER3/Helixer"; fi
[ -s "$SIF_DIR/braker3.sif" ] && ok "braker3.sif" "$(du -h "$SIF_DIR/braker3.sif" | cut -f1)" \
                              || bad "braker3.sif" "missing - run scripts/00_setup.sh"
if [ "${USE_HELIXER:-1}" = "1" ]; then
  [ -s "$SIF_DIR/helixer.sif" ] && ok "helixer.sif" "$(du -h "$SIF_DIR/helixer.sif" | cut -f1)" \
                                || warn "helixer.sif" "missing - stage 07 will fail; set USE_HELIXER=0 to skip"
fi

hr "Genome"
[ -s "$WORK/genome/genome.fna" ] && ok "assembly" "$(grep -c '^>' "$WORK/genome/genome.fna" 2>/dev/null) sequences" \
  || { [ -s "$CKPT/genome.fna" ] && warn "assembly" "only in checkpoints - will be restored" \
       || bad "assembly" "run scripts/01_prepare_genome.sh"; }
[ -s "$WORK/genome/contig_name_map.tsv" ] && ok "contig name map" || warn "contig name map" "produced by stage 01"

hr "Databases"
db_check() {  # db_check <label> <path> <required|optional> [note]
  if [ -e "$2" ]; then ok "$1"
  elif [ "$3" = required ]; then bad "$1" "${4:-missing}"
  else warn "$1" "${4:-missing - that stage degrades}"; fi
}
db_check "OrthoDB Vertebrata"   "$DB/Vertebrata.fa.gz"        required "mirror_fetch.sh + rsync to \$STAGED"
db_check "NCBI taxonomy"        "$DB/taxonomy/nodes.dmp"      required "stage 03"
db_check "Dfam $DFAM_RELEASE"   "$DB/dfam/dfam40.0.h5"        required "mirror_fetch.sh"
db_check "Rfam models"          "$DB/rfam/Rfam.cm"            optional "stage 11 loses Rfam ncRNA"
db_check "BUSCO lineages"       "$DB/busco"                   optional "stage 13 falls back to online (blocked here)"
db_check "Liftoff reference"    "$DB/liftoff_ref/ref.gff"     optional "stage 10 loses the congener axis"
db_check "Swiss-Prot"           "$DB/uniprot_sprot.fasta.gz"  optional "stage 12 loses gene symbols"
db_check "eggNOG"               "$DB/eggnog/eggnog.db"        optional "stage 12 loses GO/KEGG"
if ls -d "$DB"/interproscan-5* >/dev/null 2>&1; then ok "InterProScan"
else warn "InterProScan" "stage 12 loses protein domains"; fi
db_check "RNA-seq reads"        "$WORK/rnaseq/raw"            optional "BRAKER falls back to protein-only (the v1 flaw)"

hr "Staged mirror"
if [ -d "$STAGED" ] && [ -n "$(ls -A "$STAGED" 2>/dev/null)" ]; then
  ok "staged payload" "$(du -sh "$STAGED" 2>/dev/null | cut -f1)"
  if [ -s "$STAGED/MANIFEST.md5" ]; then
    if ( cd "$STAGED" && md5sum -c MANIFEST.md5 >/dev/null 2>&1 ); then ok "staged checksums"
    else bad "staged checksums" "MANIFEST.md5 does not verify - re-rsync"; fi
  else warn "staged checksums" "no MANIFEST.md5"; fi
else
  warn "staged payload" "$STAGED empty; firewalled databases must be rsynced in"
fi

hr "Summary"
printf '  %d passed, %d warnings, %d blockers\n\n' "$PASS" "$WARN" "$FAIL"
if [ "$FAIL" -gt 0 ]; then
  echo "  NOT ready. Resolve the blockers above, then re-run."
  exit 1
fi
if [ "$WARN" -gt 0 ]; then
  echo "  Ready, but some stages will run degraded. Each warning above names"
  echo "  exactly what is lost. Proceed with:  ./scripts/run_all.sh"
else
  echo "  Ready. Submit with:  ./scripts/run_all.sh"
fi
exit 0
