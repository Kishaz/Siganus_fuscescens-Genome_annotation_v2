#!/usr/bin/env bash
# Stage 03 - fetch all external evidence. Login node (network).
# Skips anything already staged by scripts/mirror_fetch.sh.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh
mkdir -p "$DB" "$STAGED" "$WORK/rnaseq/raw"

have() { [ -s "$1" ]; }
from_staged() {  # from_staged <filename> <dest> - use the mirrored copy if present
  if have "$STAGED/$1"; then echo "[03] using staged $1"; cp -n "$STAGED/$1" "$2"; return 0; fi
  return 1
}
get() {  # get <url> <dest>
  have "$2" && { echo "[03] have $(basename "$2")"; return; }
  from_staged "$(basename "$2")" "$2" && return
  echo "[03] downloading $(basename "$2")"
  curl -fsSL --retry 5 --retry-delay 5 -C - -o "$2.part" "$1" && mv "$2.part" "$2"
}

# ---- 1. RNA-seq (conspecific only, per project decision) --------------------
cd "$WORK/rnaseq/raw"
for acc in $RNASEQ_RUNS; do
  if ! have "${acc}_1.fastq.gz"; then
    echo "[03] prefetch+dump $acc"
    mm rnaseq prefetch --max-size 100G -O . "$acc"
    mm rnaseq fasterq-dump --split-3 --threads 8 -O . "$acc"
    mm rnaseq pigz -p 8 "${acc}"_*.fastq
  fi
done
ls -la

# ---- 2. Reference proteomes (NCBI - reachable) ------------------------------
P="$DB/proteomes"; mkdir -p "$P"; cd "$P"
python3 "$LIB/resolve_proteomes.py" --out ref_proteomes.tsv
tail -n +2 ref_proteomes.tsv | while IFS=$'\t' read -r sp acc name lvl annot url; do
  short="$(echo "$sp" | awk '{print substr($1,1,1)"_"$2}')"
  get "$url" "${short}.faa.gz"
done

# ---- 3. OrthoDB v12 Vertebrata (BLOCKED host - expect staged copy) ----------
get "$ORTHODB_URL" "$DB/Vertebrata.fa.gz"

# ---- 4. Dfam 4.0 (BLOCKED host - expect staged copy) -----------------------
D="$DB/dfam"; mkdir -p "$D"
for f in $DFAM_FILES; do get "$DFAM_BASE/$f" "$D/$f"; done
for f in $DFAM_FILES; do [ -s "$D/${f%.gz}" ] || gzip -dk "$D/$f"; done

# ---- 4b. NCBI taxonomy (to cut a clade partition out of OrthoDB) ------------
if [ "${ORTHODB_CLADE_TAXID:-0}" != "0" ]; then
  TX="$DB/taxonomy"; mkdir -p "$TX"
  get "$TAXDUMP_URL" "$TX/taxdump.tar.gz"
  [ -s "$TX/nodes.dmp" ] || tar xzf "$TX/taxdump.tar.gz" -C "$TX" nodes.dmp names.dmp
fi

# ---- 4c. Liftoff reference genome + annotation ------------------------------
# Fetched here rather than in stage 08. Compute nodes on this cluster cannot
# reach ftp.ncbi.nlm.nih.gov (conda.anaconda.org resolves, NCBI does not), so
# any NCBI download from inside a job fails. Everything network-bound lives on
# the login node.
LR="$DB/liftoff_ref"; mkdir -p "$LR"
if [ ! -s "$LR/ref.fna" ] || [ ! -s "$LR/ref.gff" ]; then
  echo "[03] fetching Liftoff reference $LIFTOFF_REF_ACC"
  ( cd "$LR"
    mm prep datasets download genome accession "$LIFTOFF_REF_ACC" \
       --include genome,gff3 --filename ref.zip
    ( unzip -o -q ref.zip -d refdl || mm prep unzip -o -q ref.zip -d refdl )
    find refdl -name '*_genomic.fna' -exec cp {} ref.fna \;
    find refdl -name 'genomic.gff'   -exec cp {} ref.gff \; )
fi
[ -s "$LR/ref.gff" ] && echo "[03] Liftoff reference: $(awk -F'\t' '$3=="gene"' "$LR/ref.gff" | wc -l) genes"

# ---- 4d. Swiss-Prot (stage 12 runs on a compute node and cannot fetch it) ---
get "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" \
    "$DB/uniprot_sprot.fasta.gz" || \
  echo "[03] WARNING: Swiss-Prot unavailable - stage 12 will skip gene-symbol assignment" >&2

# ---- 5. Rfam, BUSCO lineage -------------------------------------------------
mkdir -p "$DB/rfam" "$DB/busco"
get "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"     "$DB/rfam/Rfam.cm.gz"
get "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"    "$DB/rfam/Rfam.clanin"
[ -s "$DB/rfam/Rfam.cm" ] || gzip -dk "$DB/rfam/Rfam.cm.gz"

# ---- 6. InterProScan (~50 GB) ----------------------------------------------
# The bioconda build is pinned to 5.59_91.0 (2022) with correspondingly stale
# member databases, so the current release comes straight from EBI. Login node
# only - stage 12 runs on a compute node and cannot reach ftp.ebi.ac.uk.
# ftp.ebi.ac.uk is unreachable from this cluster (login AND compute), so unlike
# the other databases there is no download fallback: the tarball must have been
# staged by scripts/mirror_fetch.sh and rsynced into $STAGED.
if [ "${SKIP_INTERPROSCAN:-0}" != "1" ] && [ -z "$(ls -d "$DB"/interproscan-5* 2>/dev/null)" ]; then
  IPR_TGZ="$(ls "$STAGED"/interproscan-*-64-bit.tar.gz 2>/dev/null | tail -1)"
  if [ -n "$IPR_TGZ" ]; then
    IPR_VER="$(basename "$IPR_TGZ" | sed -E 's/interproscan-(.*)-64-bit\.tar\.gz/\1/')"
    echo "[03] unpacking staged InterProScan $IPR_VER (~50 GB unpacked)"
    if [ -s "${IPR_TGZ}.md5" ]; then
      want=$(awk '{print $1}' "${IPR_TGZ}.md5")
      got=$(md5sum "$IPR_TGZ" | awk '{print $1}')
      [ "$want" = "$got" ] || { echo "[03] FATAL: InterProScan md5 mismatch" >&2; exit 1; }
      echo "[03] InterProScan md5 OK"
    fi
    ( cd "$DB" && tar -xzf "$IPR_TGZ" \
      && cd "interproscan-${IPR_VER}" \
      && python3 setup.py -f interproscan.properties ) \
      || echo "[03] WARNING: InterProScan setup failed - stage 12 will skip domains" >&2
  else
    echo "[03] NOTE: no InterProScan tarball in $STAGED." >&2
    echo "     ftp.ebi.ac.uk is firewalled here, so it cannot be fetched." >&2
    echo "     Run scripts/mirror_fetch.sh elsewhere and rsync it in," >&2
    echo "     or set SKIP_INTERPROSCAN=1 to run without protein domains." >&2
  fi
fi

# ---- 7. eggNOG database (~50 GB) -------------------------------------------
if [ "${SKIP_EGGNOG:-0}" != "1" ] && [ ! -s "$DB/eggnog/eggnog.db" ]; then
  mkdir -p "$DB/eggnog"
  # Prefer a mirrored copy. The emapper database is served from
  # eggnog5.embl.de - eggnog6 answers but hosts no emapperdb, so a reachability
  # check against eggnog6 tells you nothing useful.
  staged_egg=0
  for f in eggnog.db.gz eggnog_proteins.dmnd.gz eggnog.taxa.tar.gz; do
    [ -s "$STAGED/eggnog/$f" ] && { cp -n "$STAGED/eggnog/$f" "$DB/eggnog/"; staged_egg=1; }
  done
  if [ "$staged_egg" = 1 ]; then
    echo "[03] unpacking staged eggNOG database"
    ( cd "$DB/eggnog"
      [ -s eggnog.db ]             || gzip -dk eggnog.db.gz
      [ -s eggnog_proteins.dmnd ]  || gzip -dk eggnog_proteins.dmnd.gz
      [ -d taxa ] || tar -xzf eggnog.taxa.tar.gz 2>/dev/null || true )
  else
    echo "[03] eggNOG database (~11 GB from eggnog5.embl.de)"
    mm func download_eggnog_data.py -y --data_dir "$DB/eggnog" \
      || { echo "[03] WARNING: eggNOG download failed - stage 12 loses GO/KEGG." >&2
           echo "     Mirror it with MIRROR_EGGNOG=1 scripts/mirror_fetch.sh" >&2; }
  fi
fi

echo "[03] evidence summary"; du -sh "$DB"/* "$WORK/rnaseq/raw" 2>/dev/null
echo
echo "[03] readiness for the compute stages:"
chk() { printf '  %-34s %s\n' "$1" "$([ -e "$2" ] && echo present || echo MISSING)"; }
chk "OrthoDB Vertebrata"        "$DB/Vertebrata.fa.gz"
chk "NCBI taxonomy (clade cut)" "$DB/taxonomy/nodes.dmp"
chk "Dfam ${DFAM_RELEASE}"      "$DB/dfam/dfam40.0.h5"
chk "Rfam covariance models"    "$DB/rfam/Rfam.cm"
chk "Liftoff reference"         "$DB/liftoff_ref/ref.gff"
chk "Swiss-Prot"                "$DB/uniprot_sprot.fasta.gz"
chk "eggNOG database"           "$DB/eggnog/eggnog.db"
chk "RNA-seq reads"             "$WORK/rnaseq/raw"
ls -d "$DB"/interproscan-* >/dev/null 2>&1 \
  && echo "  InterProScan                       present" \
  || echo "  InterProScan                       MISSING"

done_stamp 03_fetch_evidence
