#!/usr/bin/env bash
#SBATCH --job-name=sfus_proteins
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
# Stage 05 - protein evidence set + an independent miniprot alignment track.
#
# With only one RNA-seq library, proteins carry most of the loci that have no
# transcript support, so the panel is deliberately percomorph-weighted rather
# than relying on OrthoDB Vertebrata alone (which is what the 2024 run used).
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

restore_checkpoint "$WORK/repeats/final" genome.softmasked.fna
GM="$WORK/repeats/final/genome.softmasked.fna"; need "$GM"
T="${SLURM_CPUS_PER_TASK:-32}"
P="$DB/proteomes"; E="$WORK/evidence"; mkdir -p "$E"

# ---- assemble the protein set -----------------------------------------------
# Cut the ray-finned partition out of OrthoDB Vertebrata first. Keeping all
# 19.4 M vertebrate proteins would make stage 07 the bottleneck for hints that
# are mostly tetrapod and therefore weak donors for fish splice sites.
ODB="$DB/Vertebrata.fa.gz"
if [ "${ORTHODB_CLADE_TAXID:-0}" != "0" ] && [ -s "$DB/Vertebrata.fa.gz" ]; then
  CLADE="$DB/orthodb_clade_${ORTHODB_CLADE_TAXID}.fa.gz"
  if [ ! -s "$CLADE" ]; then
    python3 "$LIB/filter_orthodb.py" \
      --fasta "$DB/Vertebrata.fa.gz" \
      --nodes "$DB/taxonomy/nodes.dmp" --names "$DB/taxonomy/names.dmp" \
      --taxon "$ORTHODB_CLADE_TAXID" \
      --out "$CLADE" --report "$E/orthodb_partition_report.tsv"
  fi
  [ -s "$CLADE" ] && ODB="$CLADE"
fi
echo "[05] OrthoDB source: $ODB"

if [ ! -s "$E/proteins.fa" ]; then
  : > "$E/proteins.raw.fa"
  if [ -s "$ODB" ]; then
    zcat "$ODB" >> "$E/proteins.raw.fa"
  else
    echo "WARNING: OrthoDB missing - run mirror_fetch.sh + stage 04" >&2
  fi
  # Curated relatives - short evolutionary distance, better splice sites
  for f in "$P"/*.faa.gz; do
    [ -s "$f" ] || continue
    tag="$(basename "$f" .faa.gz)"
    zcat "$f" | mm annot seqkit replace -p '^' -r "${tag}_" >> "$E/proteins.raw.fa"
  done

  # Drop exact duplicates and fragments; BRAKER slows badly on a redundant set.
  mm annot seqkit rmdup -s -o "$E/proteins.dedup.fa" "$E/proteins.raw.fa"
  mm annot seqkit seq -m 50 -g "$E/proteins.dedup.fa" > "$E/proteins.fa"
  rm -f "$E/proteins.raw.fa" "$E/proteins.dedup.fa"
fi
echo "[05] protein evidence: $(grep -c '>' "$E/proteins.fa") sequences"

# ---- independent spliced-protein alignment track ----------------------------
# Not fed to BRAKER; used later to rescue loci BRAKER misses and to score
# how well-supported each final model is.
if [ ! -s "$E/miniprot.gff3" ]; then
  [ -s "$E/genome.mpi" ] || mm annot miniprot -t "$T" -d "$E/genome.mpi" "$GM"
  # Relatives only - OrthoDB Vertebrata is too large and too diverse for this.
  cat "$P"/*.faa.gz > "$E/relatives.faa.gz" 2>/dev/null || true
  zcat "$E/relatives.faa.gz" > "$E/relatives.faa"
  mm annot miniprot -t "$T" --gff --outs=0.95 -I "$E/genome.mpi" "$E/relatives.faa" \
    > "$E/miniprot.gff3"
fi
echo "[05] miniprot loci: $(awk '$3=="mRNA"' "$E/miniprot.gff3" | wc -l)"

done_stamp 05_protein_evidence
echo "[05] done -> $E"
