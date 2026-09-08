#!/usr/bin/env bash
#SBATCH --job-name=sfus_rnaseq
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
# Stage 04 - transcript evidence.
#
# HISAT2 rather than STAR: STAR would want ~35 GB for a 500 Mb genome and gains
# nothing here, while HISAT2's spliced alignment is more than adequate for hint
# generation and leaves memory for the rest of the node.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

GM="$WORK/repeats/final/genome.softmasked.fna"; need "$GM"
T="${SLURM_CPUS_PER_TASK:-32}"
RAW="$WORK/rnaseq/raw"; QC="$WORK/rnaseq/trimmed"; AL="$WORK/rnaseq/align"
mkdir -p "$QC" "$AL" "$WORK/rnaseq/stringtie"

# ---- trim -------------------------------------------------------------------
for acc in $RNASEQ_RUNS; do
  [ -s "$QC/${acc}_1.fq.gz" ] && continue
  echo "[04] fastp $acc"
  mm rnaseq fastp \
    -i "$RAW/${acc}_1.fastq.gz" -I "$RAW/${acc}_2.fastq.gz" \
    -o "$QC/${acc}_1.fq.gz"     -O "$QC/${acc}_2.fq.gz" \
    --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 40 \
    --thread 16 --json "$QC/${acc}.fastp.json" --html "$QC/${acc}.fastp.html"
done

# ---- index ------------------------------------------------------------------
IDX="$WORK/rnaseq/idx/genome"
if [ ! -s "${IDX}.1.ht2" ]; then
  mkdir -p "$(dirname "$IDX")"
  echo "[04] hisat2-build"
  mm rnaseq hisat2-build -p "$T" "$GM" "$IDX"
fi

# ---- align ------------------------------------------------------------------
for acc in $RNASEQ_RUNS; do
  BAM="$AL/${acc}.sorted.bam"
  [ -s "$BAM" ] && continue
  echo "[04] hisat2 $acc"
  mm rnaseq hisat2 -p "$T" -x "$IDX" \
      -1 "$QC/${acc}_1.fq.gz" -2 "$QC/${acc}_2.fq.gz" \
      --dta --new-summary --summary-file "$AL/${acc}.hisat2.summary" \
    | mm rnaseq samtools sort -@ 8 -m 2G -T "$FASTTMP/${acc}" -o "$BAM" -
  mm rnaseq samtools index -@ 8 "$BAM"
  mm rnaseq samtools flagstat -@ 8 "$BAM" > "$AL/${acc}.flagstat"
done

# ---- mapping-rate gate ------------------------------------------------------
# Conspecific data should map well; anything under 60% means the library is not
# what its metadata claims, and feeding it to BRAKER would inject noise.
echo "[04] mapping rates"; : > "$AL/mapping_rates.tsv"
for acc in $RNASEQ_RUNS; do
  r=$(grep -oE '^[0-9.]+% overall' "$AL/${acc}.hisat2.summary" 2>/dev/null | grep -oE '^[0-9.]+' || \
      awk '/overall alignment rate/{gsub("%","");print $1}' "$AL/${acc}.hisat2.summary")
  printf '%s\t%s\n' "$acc" "${r:-NA}" | tee -a "$AL/mapping_rates.tsv"
  awk -v a="$acc" -v x="${r:-0}" 'BEGIN{if (x+0 < 60) print "  WARNING: "a" maps at "x"% - inspect before trusting its hints"}'
done

# ---- transcript assembly ----------------------------------------------------
BAMS=""
for acc in $RNASEQ_RUNS; do
  ST="$WORK/rnaseq/stringtie/${acc}.gtf"
  [ -s "$ST" ] || mm rnaseq stringtie "$AL/${acc}.sorted.bam" -p "$T" -o "$ST" -l "$acc"
  BAMS="$BAMS $AL/${acc}.sorted.bam"
done

cd "$WORK/rnaseq/stringtie"
ls -1 *.gtf > merge_list.txt
[ -s merged.gtf ] || mm rnaseq stringtie --merge -p "$T" -o merged.gtf merge_list.txt
echo "[04] merged transcripts: $(awk '$3=="transcript"' merged.gtf | wc -l)"

echo "$BAMS" | tr ' ' '\n' | grep -v '^$' > "$WORK/rnaseq/bam_list.txt"
done_stamp 04_rnaseq
echo "[04] done"
