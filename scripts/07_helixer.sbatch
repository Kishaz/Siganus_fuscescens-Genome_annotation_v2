#!/usr/bin/env bash
#SBATCH --job-name=sfus_helixer
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#PIPE-PARTITION: gpu
# Stage 07 - Helixer: deep-learning gene prediction, GPU.
#
# Why this exists in v2 and not v1: Helixer uses NO extrinsic evidence at all.
# It predicts from sequence alone. That makes it genuinely orthogonal to
# BRAKER3, which is driven entirely by RNA-seq and protein hints. Where the two
# independently agree on a model, the agreement means something; where BRAKER
# has no evidence to work with - and with a single RNA-seq library there will
# be many such loci - Helixer still produces a call.
#
# It is NOT a replacement for BRAKER3. Helixer tends to miss lineage-specific
# and rapidly-evolving genes and does not model alternative isoforms. Stage 09
# merges the two with TSEBRA, keeping BRAKER3 as the backbone.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

if [ "${USE_HELIXER:-1}" != "1" ]; then
  echo "[07] USE_HELIXER=0 - skipping"; done_stamp 07_helixer; exit 0
fi

GM="$WORK/repeats/final/genome.softmasked.fna"; need "$GM"
H="$WORK/helixer"; mkdir -p "$H"
SIF="$SIF_DIR/helixer.sif"
[ -s "$SIF" ] || { echo "FATAL: $SIF missing - run scripts/00_setup.sh" >&2; exit 1; }

echo "[07] GPU visible to this job:"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>/dev/null \
  || echo "  WARNING: nvidia-smi unavailable; Helixer will fall back to CPU and be very slow"

# Helixer wants UPPERCASE sequence - it does its own repeat handling internally
# and soft-masked input measurably degrades its predictions.
if [ ! -s "$H/genome.upper.fna" ]; then
  mm prep seqkit seq -u -w 60 "$GM" > "$H/genome.upper.fna"
fi

# --subsequence-length 213840 is the value Helixer's own docs give for
# vertebrates; the default is tuned for much more compact genomes.
if [ ! -s "$H/helixer.gff3" ]; then
  "$CONTAINER_CMD" exec --nv \
    --bind "$PROJ:$PROJ" --bind "$FASTTMP:$FASTTMP" \
    "$SIF" Helixer.py \
      --lineage "$HELIXER_LINEAGE" \
      --fasta-path "$H/genome.upper.fna" \
      --species "$SPECIES" \
      --gff-output-path "$H/helixer.gff3" \
      --subsequence-length 213840 \
      --temporary-dir "$FASTTMP" \
    2>&1 | tee "$H/helixer.${SLURM_JOB_ID:-local}.log"
fi

need "$H/helixer.gff3"
echo "[07] Helixer models: $(awk -F'\t' '$3=="gene"' "$H/helixer.gff3" | wc -l) genes, "\
     "$(awk -F'\t' '$3=="mRNA"' "$H/helixer.gff3" | wc -l) transcripts"

done_stamp 07_helixer
echo "[07] done -> $H/helixer.gff3"
