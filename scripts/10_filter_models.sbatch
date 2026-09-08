#!/usr/bin/env bash
#SBATCH --job-name=sfus_filter
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
# Stage 10 - evidence-scored filtering, UTRs, and stable IDs.
#
# The published annotation reported 28,351 genes; teleost genomes typically hold
# 20-26k. Rather than trimming to hit a target number, every model is scored on
# independent axes and only those failing ALL of them are dropped:
#   - protein homology  (DIAMOND vs the curated relative panel)
#   - expression        (StringTie -e TPM from the conspecific library)
#   - congener support  (overlap with the Liftoff transfer from S. canaliculatus)
#   - repeat overlap    (fraction of CDS inside a repeat feature)
#
# v2 adds the congener axis. With one RNA-seq library, homology and expression
# are correlated failures for the same loci - a poorly conserved gene that is
# not expressed in that tissue fails both for one underlying reason. The lift
# from a chromosome-level congener is an independent way for such a locus to
# prove it is real.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

GM="$WORK/repeats/final/genome.softmasked.fna"
# v2 input is the TSEBRA-consolidated BRAKER3 + Helixer set, not BRAKER alone.
RAW="$WORK/combined/combined.gff3"
[ -s "$RAW" ] || RAW="$WORK/braker/run/braker.gff3"
REP="$WORK/repeats/final/${ASM_NAME}.repeats.gff3"
need "$GM" "$RAW" "$REP"
echo "[10] input models: $RAW"
T="${SLURM_CPUS_PER_TASK:-32}"
F="$WORK/filter"; mkdir -p "$F"; cd "$F"

# ---- 0. normalise ------------------------------------------------------------
[ -s braker.clean.gff3 ] || mm annot agat_sp_fix_features_locations_duplicated.pl \
    --gff "$RAW" -o braker.fixed.gff3 >/dev/null 2>&1 || cp "$RAW" braker.fixed.gff3
[ -s braker.clean.gff3 ] || mm annot agat_convert_sp_gxf2gxf.pl \
    -g braker.fixed.gff3 -o braker.clean.gff3 >/dev/null

[ -s proteins.faa ] || mm annot agat_sp_extract_sequences.pl \
    -g braker.clean.gff3 -f "$GM" -p -o proteins.faa >/dev/null

# ---- 1. protein homology -----------------------------------------------------
if [ ! -s homology.tsv ]; then
  cat "$DB/proteomes"/*.faa.gz > rel.faa.gz
  [ -s rel.dmnd ] || { zcat rel.faa.gz | mm annot diamond makedb --in - -d rel -p "$T"; }
  mm annot diamond blastp -q proteins.faa -d rel -p "$T" \
      --very-sensitive --max-target-seqs 1 --evalue 1e-5 \
      --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp \
      -o homology.tsv
fi

# ---- 2. expression -----------------------------------------------------------
if [ ! -s expression.tsv ] && [ -s "$WORK/rnaseq/bam_list.txt" ]; then
  mm annot agat_convert_sp_gff2gtf.pl -g braker.clean.gff3 -o braker.clean.gtf >/dev/null
  : > expression.tsv
  while read -r bam; do
    [ -s "$bam" ] || continue
    s=$(basename "$bam" .sorted.bam)
    mm rnaseq stringtie -e -G braker.clean.gtf -p "$T" -o "expr.${s}.gtf" "$bam"
    awk -v s="$s" '$3=="transcript"{
        match($0,/transcript_id "[^"]+"/); tid=substr($0,RSTART+15,RLENGTH-16);
        match($0,/TPM "[^"]+"/);           tpm=substr($0,RSTART+5,RLENGTH-6);
        print tid"\t"s"\t"tpm }' "expr.${s}.gtf" >> expression.tsv
  done < "$WORK/rnaseq/bam_list.txt"
fi
[ -s expression.tsv ] || : > expression.tsv

# ---- 3. repeat overlap -------------------------------------------------------
if [ ! -s repeat_overlap.tsv ]; then
  awk 'BEGIN{OFS="\t"} $3=="CDS"{
        match($0,/Parent=[^;]+/); p=substr($0,RSTART+7,RLENGTH-7);
        print $1,$4-1,$5,p }' braker.clean.gff3 | sort -k1,1 -k2,2n > cds.bed
  awk 'BEGIN{OFS="\t"} !/^#/{print $1,$4-1,$5}' "$REP" | sort -k1,1 -k2,2n \
    | mm annot bedtools merge -i - > rep.merged.bed
  mm annot bedtools intersect -a cds.bed -b rep.merged.bed -wao \
    | awk 'BEGIN{OFS="\t"}{cds[$4]+=$3-$2; ov[$4]+=$NF}
           END{for(k in cds) print k, cds[k], ov[k], (cds[k]?ov[k]/cds[k]:0)}' \
    > repeat_overlap.tsv
fi

# ---- 3b. congener support (Liftoff from S. canaliculatus) -------------------
# Fraction of each model's CDS that a lifted congener CDS also covers.
if [ ! -s liftoff_overlap.tsv ]; then
  if [ -s "$WORK/liftoff/liftoff_cds.bed" ]; then
    mm annot bedtools intersect -a cds.bed -b "$WORK/liftoff/liftoff_cds.bed" -wao \
      | awk 'BEGIN{OFS="\t"}{cds[$4]+=$3-$2; ov[$4]+=$NF}
             END{for(k in cds) print k, cds[k], ov[k], (cds[k]?ov[k]/cds[k]:0)}' \
      > liftoff_overlap.tsv
  else
    echo "[10] no Liftoff track - scoring without the congener axis"
    : > liftoff_overlap.tsv
  fi
fi

# ---- 4. score and filter -----------------------------------------------------
python3 "$LIB/filter_models.py" \
  --gff3 braker.clean.gff3 \
  --proteins proteins.faa \
  --homology homology.tsv \
  --expression expression.tsv \
  --repeat-overlap repeat_overlap.tsv \
  --liftoff-overlap liftoff_overlap.tsv \
  --keep-gff3 filtered.gff3 \
  --drop-gff3 dropped.gff3 \
  --report evidence_support.tsv \
  --summary filter_summary.txt
cat filter_summary.txt

# ---- 5. UTRs from the assembled transcripts ---------------------------------
if [ -s "$WORK/rnaseq/stringtie/merged.gtf" ]; then
  python3 "$LIB/add_utrs.py" \
    --gff3 filtered.gff3 \
    --transcripts "$WORK/rnaseq/stringtie/merged.gtf" \
    --out with_utr.gff3 --report utr_report.tsv
else
  cp filtered.gff3 with_utr.gff3
fi

# ---- 6. stable IDs -----------------------------------------------------------
python3 "$LIB/assign_ids.py" \
  --gff3 with_utr.gff3 --prefix "$GENE_PREFIX" \
  --out "$WORK/filter/${ASM_NAME}.annotation.${ANNOT_VERSION}.gff3" \
  --map id_map.tsv

mm annot agat_sp_statistics.pl --gff "$WORK/filter/${ASM_NAME}.annotation.${ANNOT_VERSION}.gff3" \
  -o agat_stats.txt >/dev/null 2>&1 || true

done_stamp 10_filter_models
echo "[10] done"
