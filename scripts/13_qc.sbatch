#!/usr/bin/env bash
#SBATCH --job-name=sfus_qc
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
# Stage 13 - validation. Nothing ships until this passes.
#
# Acceptance criteria, stated in advance so the result cannot be rationalised:
#   BUSCO protein mode C >= 92%   (assembly itself is 98.5%, so the annotation
#                                  should recover nearly all of it)
#   gene count 20,000-27,000      (teleost range; 2024 release was 28,351)
#   GFF3 validates                (gt gff3validator, zero errors)
#   0 internal stop codons
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

GFF="$WORK/functional/${ASM_NAME}.annotation.${ANNOT_VERSION}.functional.gff3"
[ -s "$GFF" ] || GFF="$WORK/filter/${ASM_NAME}.annotation.${ANNOT_VERSION}.gff3"
restore_checkpoint "$WORK/repeats/final" genome.softmasked.fna
GM="$WORK/repeats/final/genome.softmasked.fna"
need "$GFF" "$GM"
T="${SLURM_CPUS_PER_TASK:-32}"
Q="$WORK/qc"; mkdir -p "$Q"; cd "$Q"

# ---- sequence products -------------------------------------------------------
[ -s proteins.faa    ] || mm annot agat_sp_extract_sequences.pl -g "$GFF" -f "$GM" -p        -o proteins.faa    >/dev/null
[ -s cds.fna         ] || mm annot agat_sp_extract_sequences.pl -g "$GFF" -f "$GM" -t cds    -o cds.fna         >/dev/null
[ -s transcripts.fna ] || mm annot agat_sp_extract_sequences.pl -g "$GFF" -f "$GM" -t exon --merge -o transcripts.fna >/dev/null

# ---- structural validity ----------------------------------------------------
echo "=== gt gff3validator ==="
mm annot gt gff3validator "$GFF" > gff3validator.log 2>&1 \
  && echo "  GFF3 valid" \
  || { echo "  GFF3 INVALID - see $Q/gff3validator.log"; tail -20 gff3validator.log; }

echo "=== internal stop codons ==="
ISTOP=$(mm annot seqkit fx2tab proteins.faa | awk -F'\t' '{s=$2; sub(/\*$/,"",s); if (index(s,"*")) n++} END{print n+0}')
echo "  $ISTOP"

# ---- BUSCO in protein mode --------------------------------------------------
# Run against every configured lineage. odb10 is the gated one because the 2024
# paper reported odb10 (98.5%, genome mode) and the comparison only means
# something on the same lineage; odb12 is reported for current-practice readers.
for L in ${BUSCO_LINEAGES:-$BUSCO_LINEAGE}; do
  OUTDIR="busco_prot_${L}"
  if [ ! -d "$OUTDIR/run_${L}" ]; then
    echo "=== BUSCO $L ==="
    mm qc busco -i proteins.faa -m proteins -l "$L" \
      -o "$OUTDIR" --out_path "$Q" -c "$T" --offline \
      --download_path "$DB/busco" 2>&1 | tail -20 \
    || mm qc busco -i proteins.faa -m proteins -l "$L" \
      -o "$OUTDIR" --out_path "$Q" -c "$T" 2>&1 | tail -20 \
    || echo "  WARNING: BUSCO failed for lineage $L"
  fi
  grep -hE '^\s+C:' "$OUTDIR"/short_summary.specific.*.txt 2>/dev/null \
    | sed "s/^/  $L /" || true
done

# The gate reads the primary lineage only.
BUSCO_SUMMARY="$(ls "busco_prot_${BUSCO_LINEAGE}"/short_summary.specific.*.txt 2>/dev/null | head -1)"

# ---- compleasm --------------------------------------------------------------
# A miniprot-based reimplementation of the same lineages. BUSCO's protein mode
# is known to report inflated fragmentation on predicted proteomes; compleasm
# is a genuinely different aligner, so a large gap between the two is itself
# diagnostic rather than something to average away. Genome mode also lets us
# state directly how much of the assembly's own BUSCO content the annotation
# actually recovered - which is the number that matters here.
if [ ! -s "compleasm_genome/summary.txt" ]; then
  echo "=== compleasm (genome mode, $BUSCO_LINEAGE) ==="
  mm qc compleasm run -a "$GM" -o compleasm_genome -l "$BUSCO_LINEAGE" \
     -t "$T" --library_path "$DB/busco" 2>&1 | tail -12 \
  || echo "  WARNING: compleasm genome run failed"
fi
[ -s compleasm_genome/summary.txt ] && sed 's/^/  genome  /' compleasm_genome/summary.txt

python3 "$LIB/compare_completeness.py" \
  --busco "$BUSCO_SUMMARY" \
  --compleasm compleasm_genome/summary.txt \
  --assembly-busco-complete 98.5 \
  --out completeness_comparison.tsv --report completeness_report.txt \
  2>/dev/null || true
[ -s completeness_report.txt ] && cat completeness_report.txt

# ---- OMArk: proteome consistency, fragmentation, contamination --------------
if command -v omark >/dev/null 2>&1 || mm qc omark --version >/dev/null 2>&1; then
  [ -s "$DB/omark/LUCA.h5" ] || echo "  (OMArk database absent - skipping)"
  if [ -s "$DB/omark/LUCA.h5" ]; then
    mm qc omamer search --db "$DB/omark/LUCA.h5" --query proteins.faa --out omamer.out
    mm qc omark -f omamer.out -d "$DB/omark/LUCA.h5" -o omark_out
    cat omark_out/*.sum 2>/dev/null | head -30
  fi
fi

# ---- comparative statistics -------------------------------------------------
python3 "$LIB/annotation_stats.py" \
  --gff3 "$GFF" --proteins proteins.faa \
  --busco "$BUSCO_SUMMARY" \
  --evidence "$WORK/filter/evidence_support.tsv" \
  --out annotation_stats.tsv --report annotation_report.txt \
  --check-thresholds
cat annotation_report.txt

done_stamp 13_qc
echo "[13] done -> $Q"
