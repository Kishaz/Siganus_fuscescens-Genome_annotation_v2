#!/usr/bin/env bash
#SBATCH --job-name=sfus_combine
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --nodes=1
# Stage 09 - consolidate BRAKER3 and Helixer into one gene set.
#
# TSEBRA was written for exactly this: it scores transcripts against the
# extrinsic evidence BRAKER already generated (hintsfile.gff) and picks between
# competing models at a locus, rather than naively unioning them - which would
# inflate the gene count and produce overlapping duplicates.
#
# BRAKER3 is the anchor. Helixer contributes where BRAKER has nothing, which
# with a single RNA-seq library is a real gap. If Helixer was skipped, this
# stage passes BRAKER through unchanged so downstream stages do not branch.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

BRK="$WORK/braker/run/braker.gff3"
HLX="$WORK/helixer/helixer.gff3"
need "$BRK"
C="$WORK/combined"; mkdir -p "$C"; cd "$C"

if [ "${USE_HELIXER:-1}" != "1" ] || [ ! -s "$HLX" ]; then
  echo "[09] no Helixer input - passing BRAKER3 through unchanged"
  cp "$BRK" combined.gff3
  done_stamp 09_combine
  echo "[09] done -> $C/combined.gff3"
  exit 0
fi

HINTS="$WORK/braker/run/hintsfile.gff"
[ -s "$HINTS" ] || HINTS="$(find "$WORK/braker/run" -name 'hintsfile.gff' | head -1)"

# TSEBRA consumes GTF. AGAT round-trips both inputs so the two agree in format.
[ -s braker.gtf  ] || mm annot agat_convert_sp_gff2gtf.pl -g "$BRK" -o braker.gtf  >/dev/null
[ -s helixer.gtf ] || mm annot agat_convert_sp_gff2gtf.pl -g "$HLX" -o helixer.gtf >/dev/null

if [ -s "$HINTS" ]; then
  # pref_braker1.cfg keeps BRAKER models unless Helixer is clearly better
  # supported; that is the conservative direction for an evidence-poor genome.
  CFG="$(mm annot python3 -c 'import os,shutil;p=shutil.which("tsebra");print(os.path.join(os.path.dirname(os.path.realpath(p)),"..","config","pref_braker1.cfg"))' 2>/dev/null || true)"
  TSEBRA_ARGS=(-g braker.gtf,helixer.gtf -e "$HINTS" -o tsebra.gtf)
  [ -s "$CFG" ] && TSEBRA_ARGS+=(-c "$CFG")
  if mm annot tsebra "${TSEBRA_ARGS[@]}" 2>&1 | tee "tsebra.${SLURM_JOB_ID:-local}.log"; then
    mm annot agat_convert_sp_gxf2gxf.pl -g tsebra.gtf -o combined.gff3 >/dev/null
  else
    echo "[09] WARNING: TSEBRA failed; falling back to BRAKER3 alone" >&2
    cp "$BRK" combined.gff3
  fi
else
  echo "[09] WARNING: no BRAKER hintsfile found, so TSEBRA cannot score models." >&2
  echo "     Falling back to BRAKER3 alone rather than blindly unioning." >&2
  cp "$BRK" combined.gff3
fi

need combined.gff3
b=$(awk -F'\t' '$3=="gene"' "$BRK" | wc -l)
h=$(awk -F'\t' '$3=="gene"' "$HLX" | wc -l)
c=$(awk -F'\t' '$3=="gene"' combined.gff3 | wc -l)
printf '[09] genes: BRAKER3 %s | Helixer %s | combined %s\n' "$b" "$h" "$c"

# How much did each source actually contribute? gffcompare against BRAKER shows
# what Helixer added rather than leaving the merge as a black box.
mm annot gffcompare -r "$BRK" -o vs_braker combined.gff3 >/dev/null 2>&1 || true
[ -s vs_braker.stats ] && sed -n '1,25p' vs_braker.stats

done_stamp 09_combine
echo "[09] done -> $C/combined.gff3"
