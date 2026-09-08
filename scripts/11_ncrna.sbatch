#!/usr/bin/env bash
#SBATCH --job-name=sfus_ncrna
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
# Stage 13 - non-coding RNA. The 2024 release shipped none (14-tRNA-scan.sh was
# an empty file), so this is new content rather than a refinement.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

GM="$WORK/repeats/final/genome.softmasked.fna"; need "$GM"
T="${SLURM_CPUS_PER_TASK:-32}"
N="$WORK/ncrna"; mkdir -p "$N"; cd "$N"

# ---- tRNA -------------------------------------------------------------------
if [ ! -s trnascan.out ]; then
  mm ncrna tRNAscan-SE -E --thread "$T" \
    -o trnascan.out -f trnascan.struct -m trnascan.stats -a trnascan.fa "$GM"
fi
# EukHighConfidenceFilter removes the pseudogene-like hits that otherwise
# dominate a fish tRNA set.
if [ ! -s filtered/trnascan.hc.out ]; then
  mkdir -p filtered
  mm ncrna EukHighConfidenceFilter -i trnascan.out -s trnascan.struct \
     -o filtered -p trnascan.hc || cp trnascan.out filtered/trnascan.hc.out
fi

# ---- rRNA -------------------------------------------------------------------
[ -s rrna.gff3 ] || mm ncrna barrnap --kingdom euk --threads "$T" \
    --outseq rrna.fa "$GM" > rrna.gff3

# ---- Rfam / Infernal --------------------------------------------------------
if [ ! -s rfam.tblout ]; then
  [ -s "$DB/rfam/Rfam.cm.i1f" ] || mm ncrna cmpress "$DB/rfam/Rfam.cm"
  # -Z is the doubled genome size in Mb, required for correct E-values.
  ZS=$(awk '/^total_bp/{printf "%.2f", $2*2/1000000}' "$WORK/repeats/final/stats.softmasked.txt" \
       2>/dev/null || echo 996.34)
  mm ncrna cmscan --cpu "$T" --rfam --cut_ga --nohmmonly --noali \
      -Z "$ZS" --fmt 2 --clanin "$DB/rfam/Rfam.clanin" \
      --tblout rfam.tblout -o rfam.cmscan "$DB/rfam/Rfam.cm" "$GM"
fi

# ---- miRNA ------------------------------------------------------------------
# Rfam finds miRNA families but MirMachine uses family-specific covariance
# models built from MirGeneDB, which is markedly more precise for metazoan
# miRNAs. Its calls take precedence where the two overlap.
if [ ! -s mirmachine/mirmachine.gff ]; then
  mkdir -p mirmachine
  mm ncrna MirMachine.py -n Actinopterygii -s "$SPECIES" \
     -g "$GM" -m proto --cpu "$T" -o mirmachine \
     2>&1 | tail -20 || echo "[11] WARNING: MirMachine failed - Rfam miRNAs still stand"
  find mirmachine -name '*.PRE.gff' -o -name '*.gff' 2>/dev/null | head -1 \
    | xargs -r -I{} cp {} mirmachine/mirmachine.gff
fi

python3 "$LIB/ncrna_to_gff3.py" \
  --trnascan filtered/trnascan.hc.out --rrna rrna.gff3 --rfam rfam.tblout \
  --mirmachine mirmachine/mirmachine.gff \
  --out "$N/${ASM_NAME}.ncRNA.gff3" --summary "$N/ncrna_summary.tsv"
cat "$N/ncrna_summary.tsv"

done_stamp 11_ncrna
echo "[11] done"
