#!/usr/bin/env bash
#SBATCH --job-name=sfus_liftoff
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
# Stage 08 - lift the congener annotation across as an independent track.
#
# GCA_053572195.1 (Siganus canaliculatus, chromosome-level, annotated) did not
# exist when the 2024 paper was written. It is the closest annotated genome to
# ours by a wide margin - same genus, and the two species are so close that much
# of the East Asian literature conflates them.
#
# This is NOT used to define gene models: a cross-species lift inherits the
# reference's errors and cannot see genes the reference lacks. It is used as
#   1. an independent evidence axis when scoring our own models (stage 10), and
#   2. a recall check - loci Liftoff places but BRAKER and Helixer both miss are
#      exactly where to look for systematic failure.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

restore_checkpoint "$WORK/repeats/final" genome.softmasked.fna
GM="$WORK/repeats/final/genome.softmasked.fna"; need "$GM"
T="${SLURM_CPUS_PER_TASK:-32}"
L="$WORK/liftoff"; mkdir -p "$L"; cd "$L"

# ---- reference genome + annotation ------------------------------------------
# NOT downloaded here. Compute nodes on this cluster cannot reach
# ftp.ncbi.nlm.nih.gov, so stage 03 stages the reference on the login node.
LR="$DB/liftoff_ref"
if [ ! -s "$LR/ref.fna" ] || [ ! -s "$LR/ref.gff" ]; then
  echo "[08] WARNING: $LIFTOFF_REF_ACC not staged in $LR - skipping Liftoff." >&2
  echo "     Run scripts/03_fetch_evidence.sh on the login node first." >&2
  echo "     Stage 10 will score without the congener axis." >&2
  : > "$L/liftoff.gff3"
  done_stamp 08_liftoff
  exit 0
fi
ln -sf "$LR/ref.fna" ref.fna
ln -sf "$LR/ref.gff" ref.gff
echo "[08] reference: $(awk -F'\t' '$3=="gene"' ref.gff | wc -l) genes"

# ---- lift -------------------------------------------------------------------
# -exclude_partial keeps fragments out; -copies finds duplications, which is
# relevant given the moderate duplication the 2024 paper reported. -flank helps
# across the ~5-8% divergence between these congeners.
if [ ! -s liftoff.gff3 ]; then
  mm annot liftoff \
    -g ref.gff \
    -o liftoff.gff3 \
    -u unmapped_features.txt \
    -dir "$FASTTMP/liftoff_intermediate" \
    -p "$T" \
    -copies -sc 0.95 -flank 0.1 -exclude_partial \
    "$GM" ref.fna \
    2>&1 | tee "liftoff.${SLURM_JOB_ID:-local}.log"
fi

need liftoff.gff3
LIFTED=$(awk -F'\t' '$3=="gene"' liftoff.gff3 | wc -l)
REFN=$(awk -F'\t' '$3=="gene"' ref.gff | wc -l)
UNMAPPED=$(wc -l < unmapped_features.txt 2>/dev/null || echo 0)
echo "[08] lifted $LIFTED / $REFN reference genes ($UNMAPPED unmapped)"
awk -v a="$LIFTED" -v b="$REFN" 'BEGIN{
  r = (b>0 ? 100*a/b : 0);
  printf "  transfer rate %.1f%%\n", r;
  if (r < 60) print "  WARNING: low transfer rate - check assembly/reference compatibility"
}'

# BED of lifted CDS, for the stage-10 evidence scoring.
awk 'BEGIN{OFS="\t"} !/^#/ && $3=="CDS"{print $1,$4-1,$5}' liftoff.gff3 \
  | sort -k1,1 -k2,2n | mm annot bedtools merge -i - > liftoff_cds.bed

done_stamp 08_liftoff
echo "[08] done -> $L/liftoff.gff3"
