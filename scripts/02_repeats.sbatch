#!/usr/bin/env bash
#SBATCH --job-name=sfus_repeats
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
# Stage 02 - repeat discovery, curation and softmasking. The long pole: 3-7 days.
#
# Default route is EarlGrey 7.3.1, which is what changed between v1 and v2.
# The 2024 annotation used bare RepeatModeler2 + RepeatMasker and left a large
# "Unknown" fraction; unclassified families make the repeat track hard to reuse
# and make it impossible to judge which gene models are TE-derived. EarlGrey
# adds BLAST-extend-extend boundary curation (RepeatModeler consensi are
# routinely truncated), classification, and a divergence landscape.
#
# REPEAT_METHOD=manual falls back to the v1 route.
#
# QC anchor: NCBI's own softmasking of this assembly is 21.15%. Landing far
# below means a pass silently failed; far above means over-masking.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

T="${SLURM_CPUS_PER_TASK:-$THREADS_BIG}"
R="$WORK/repeats"; OUT="$R/final"; mkdir -p "$R" "$OUT" "$WORK/genome"

# Scratch may have been purged between submission and start, or between a
# failed attempt and a resubmission. Pull anything durable back first.
restore_checkpoint "$WORK/genome" genome.fna
restore_checkpoint "$OUT" genome.softmasked.fna
restore_checkpoint "$R" "${SPECIES}-families.curated.fa"

G="$WORK/genome/genome.fna"
if [ ! -s "$G" ]; then
  echo "FATAL: $G missing and no checkpoint available." >&2
  echo "       Re-run scripts/01_prepare_genome.sh on the login node." >&2
  exit 1
fi

# A completed mask is days of compute; never redo it just because scratch moved.
if [ -s "$OUT/genome.softmasked.fna" ] && have_stamp 02_repeats; then
  echo "[02] softmasked genome already present - nothing to do"
  exit 0
fi

# Dfam 4.0 has to be visible to RepeatMasker before either route runs.
if [ -s "$DB/dfam/dfam40.0.h5" ]; then
  RM_LIBDIR="$("$MICROMAMBA_BIN" run -r "$MAMBA_ROOT_PREFIX" -n repeats \
    python3 -c 'import os,shutil;p=shutil.which("RepeatMasker");print(os.path.join(os.path.dirname(os.path.realpath(p)),"Libraries"))' 2>/dev/null || true)"
  if [ -n "$RM_LIBDIR" ] && [ -d "$RM_LIBDIR" ]; then
    mkdir -p "$RM_LIBDIR/famdb"
    for f in "$DB"/dfam/*.h5; do
      ln -sf "$f" "$RM_LIBDIR/famdb/$(basename "$f")" 2>/dev/null || true
    done
    echo "[02] linked Dfam $DFAM_RELEASE into $RM_LIBDIR/famdb"
  fi
fi

case "${REPEAT_METHOD:-earlgrey}" in

# ------------------------------------------------------------------ EarlGrey --
earlgrey)
  EG="$R/earlgrey"; mkdir -p "$EG"
  if [ ! -s "$EG/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.softmasked.fasta" ]; then
    echo "[02] EarlGrey on $(basename "$G")"
    mm repeats earlGrey \
      -g "$G" \
      -s "$SPECIES" \
      -o "$EG" \
      -t "$T" \
      -d yes \
      2>&1 | tee "$R/earlgrey.${SLURM_JOB_ID:-local}.log"
  fi
  S="$EG/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles"
  need "$S/${SPECIES}.softmasked.fasta"
  cp "$S/${SPECIES}.softmasked.fasta" "$OUT/genome.softmasked.fna"
  # EarlGrey emits a GFF directly; keep it and derive our summary from the .out
  cp "$S/${SPECIES}.filteredRepeats.gff" "$OUT/${ASM_NAME}.repeats.gff3" 2>/dev/null || true
  cp "$S"/*.familyLevelCount.txt "$OUT/" 2>/dev/null || true
  cp "$S"/*highLevelCount.txt    "$OUT/" 2>/dev/null || true
  cp "$S"/*.pdf                  "$OUT/" 2>/dev/null || true
  cp "$EG/${SPECIES}_EarlGrey/${SPECIES}_RepeatMasker_Against_Custom_Library/"*.out \
     "$OUT/full_mask.out" 2>/dev/null || true
  cp "$EG/${SPECIES}_EarlGrey/${SPECIES}_Curated_Library/"*.fasta \
     "$R/${SPECIES}-families.curated.fa" 2>/dev/null || true
  ;;

# ------------------------------------------------------------------- manual --
manual)
  DBNAME="$SPECIES"
  RM="$R/modeler"; mkdir -p "$RM"; cd "$RM"
  [ -s "${DBNAME}.translation" ] || mm repeats BuildDatabase -name "$DBNAME" -engine ncbi "$G"
  [ -s "${DBNAME}-families.fa" ] || mm repeats RepeatModeler \
      -database "$DBNAME" -threads "$T" -engine ncbi -LTRStruct \
      2>&1 | tee "$R/RepeatModeler.${SLURM_JOB_ID:-local}.log"
  need "${DBNAME}-families.fa"

  mm repeats seqkit fx2tab "${DBNAME}-families.fa" \
    | awk 'BEGIN{FS=OFS="\t"}{$1="Sfus_"$1; print}' \
    | mm repeats seqkit tab2fx -w 60 > "$R/${SPECIES}-families.fa"

  CL="$R/classified"; mkdir -p "$CL"; cd "$CL"
  mm teclass TEsorter "$R/${SPECIES}-families.fa" -db rexdb-metazoa -p "$T" || true
  python3 "$LIB/merge_te_classification.py" \
    --library "$R/${SPECIES}-families.fa" \
    --tesorter "$(basename "$R/${SPECIES}-families.fa").rexdb-metazoa.cls.tsv" \
    --out "${SPECIES}-families.classified.fa" \
    --report "${SPECIES}-classification_report.tsv"
  LIBUSE="$CL/${SPECIES}-families.classified.fa"; [ -s "$LIBUSE" ] || LIBUSE="$R/${SPECIES}-families.fa"

  PA=$(( T / 4 )); [ "$PA" -lt 1 ] && PA=1
  P1="$R/mask_pass1"; P2="$R/mask_pass2"; mkdir -p "$P1" "$P2"
  [ -s "$P1/genome.fna.masked" ] || mm repeats RepeatMasker \
      -lib "$LIBUSE" -xsmall -gff -poly -a -pa "$PA" -dir "$P1" "$G"
  [ -s "$P2/genome.fna.masked.masked" ] || mm repeats RepeatMasker \
      -species Teleostei -xsmall -gff -poly -a -pa "$PA" -dir "$P2" "$P1/genome.fna.masked"

  cp "$P2"/*.masked "$OUT/genome.softmasked.fna"
  cat "$P1"/*.out <(tail -n +4 "$P2"/*.out) > "$OUT/full_mask.out"
  ;;

*) echo "FATAL: REPEAT_METHOD must be 'earlgrey' or 'manual'" >&2; exit 2 ;;
esac

# ---------------------------------------------------------------- normalise ---
cd "$REPO"
need "$OUT/genome.softmasked.fna"
mm prep samtools faidx "$OUT/genome.softmasked.fna"

# Always regenerate our own GFF3 + summary from the .out table: EarlGrey's GFF
# is fine but its attribute set differs from the manual route, and the release
# must look identical whichever route produced it.
if [ -s "$OUT/full_mask.out" ]; then
  python3 "$LIB/repeatmasker_to_gff3.py" \
    --out "$OUT/full_mask.out" \
    --gff3 "$OUT/${ASM_NAME}.repeats.gff3" \
    --summary "$OUT/repeat_summary.tsv" \
    --source "EarlGrey7/RepeatMasker/Dfam${DFAM_RELEASE}"
fi

# ------------------------------------------------------------------- QC ------
echo "=== our softmask ==="
python3 "$LIB/fasta_tools.py" stats --fasta "$OUT/genome.softmasked.fna" \
  | tee "$OUT/stats.softmasked.txt"
echo "=== NCBI softmask of the same assembly (independent anchor: 21.15%) ==="
cat "$WORK/genome/stats.ncbi_masked.txt" 2>/dev/null || true

# Masking changes case, never bases. If this fails, something rewrote sequence.
python3 "$LIB/fasta_tools.py" verify \
  --a "$OUT/genome.softmasked.fna" --b "$WORK/genome/genome.fna" \
  || { echo "FATAL: masking altered sequence content" >&2; exit 1; }

PCT=$(awk '/^softmasked_pct/{print $2}' "$OUT/stats.softmasked.txt")
awk -v p="$PCT" 'BEGIN{
  if (p+0 < 12) print "WARNING: only " p "% masked - well below the 21.15% NCBI anchor; a pass may have failed";
  else if (p+0 > 35) print "WARNING: " p "% masked - well above the NCBI anchor; check for over-masking";
  else print "softmask " p "% - consistent with the NCBI anchor" }'

# Days of compute - copy to durable storage before the job exits.
checkpoint \
  "$OUT/genome.softmasked.fna" \
  "$OUT/${ASM_NAME}.repeats.gff3" \
  "$OUT/repeat_summary.tsv" \
  "$OUT/full_mask.out" \
  "$R/${SPECIES}-families.curated.fa" \
  "$R/${SPECIES}-families.fa"

done_stamp 02_repeats
echo "[02] done -> $OUT"
