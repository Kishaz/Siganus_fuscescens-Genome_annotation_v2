#!/usr/bin/env bash
#SBATCH --job-name=sfus_braker3
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
# Stage 06 - BRAKER3 (GeneMark-ETP + AUGUSTUS + TSEBRA) with BOTH evidence types.
#
# This is the substantive upgrade over the published annotation, which ran
# protein-only. Runs from the official container: the conda build of BRAKER3
# routinely breaks on AUGUSTUS script paths and GeneMark perl deps.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

restore_checkpoint "$WORK/repeats/final" genome.softmasked.fna
GM="$WORK/repeats/final/genome.softmasked.fna"
PROT="$WORK/evidence/proteins.fa"
need "$GM" "$PROT"
T="${SLURM_CPUS_PER_TASK:-48}"

B="$WORK/braker"; mkdir -p "$B"
SIF="$SIF_DIR/braker3.sif"
[ -s "$SIF" ] || { echo "FATAL: $SIF missing - run scripts/00_setup.sh" >&2; exit 1; }

# AUGUSTUS needs a writable config tree; the one inside the image is read-only.
ACFG="$B/augustus_config"
if [ ! -d "$ACFG/species" ]; then
  echo "[06] extracting AUGUSTUS config from container"
  "$CONTAINER_CMD" exec "$SIF" bash -c 'cp -r /opt/Augustus/config /tmp/acfg && echo ok' >/dev/null
  "$CONTAINER_CMD" exec --bind "$B:$B" "$SIF" cp -r /opt/Augustus/config "$ACFG"
fi

BAMS="$(paste -sd, "$WORK/rnaseq/bam_list.txt" 2>/dev/null || true)"
RNAOPT=()
if [ -n "$BAMS" ]; then
  RNAOPT=(--bam="$BAMS")
  echo "[06] RNA-seq evidence: $BAMS"
else
  echo "[06] WARNING: no BAMs found - falling back to protein-only, which is" >&2
  echo "     exactly what the 2024 annotation did. Check stage 05 first."   >&2
fi

# --species must be unique per run or AUGUSTUS refuses to retrain.
SPP="${SPECIES}_${ANNOT_VERSION//./_}"

"$CONTAINER_CMD" exec \
  --bind "$PROJ:$PROJ" --bind "$FASTTMP:$FASTTMP" \
  --env AUGUSTUS_CONFIG_PATH="$ACFG" \
  "$SIF" braker.pl \
    --genome="$GM" \
    --prot_seq="$PROT" \
    "${RNAOPT[@]}" \
    --species="$SPP" \
    --workingdir="$B/run" \
    --threads="$T" \
    --gff3 \
    --softmasking \
    --AUGUSTUS_CONFIG_PATH="$ACFG" \
  2>&1 | tee "$B/braker.${SLURM_JOB_ID:-local}.log"

need "$B/run/braker.gff3"
echo "[06] raw models: $(awk '$3=="gene"' "$B/run/braker.gff3" | wc -l) genes, "\
     "$(awk '$3=="mRNA"' "$B/run/braker.gff3" | wc -l) transcripts"

checkpoint "$B/run/braker.gff3" "$B/run/hintsfile.gff"

done_stamp 06_braker3
echo "[06] done -> $B/run/braker.gff3"
