#!/usr/bin/env bash
#SBATCH --job-name=sfus_functional
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
# Stage 12 - functional annotation. Entirely new: the 2024 release shipped no
# GO terms, KEGG assignments, domains, or gene symbols.
set -euo pipefail
cd "${REPO:?}"; source ./config.sh

GFF="$WORK/filter/${ASM_NAME}.annotation.${ANNOT_VERSION}.gff3"
restore_checkpoint "$WORK/repeats/final" genome.softmasked.fna
GM="$WORK/repeats/final/genome.softmasked.fna"
need "$GFF" "$GM"
T="${SLURM_CPUS_PER_TASK:-48}"
F="$WORK/functional"; mkdir -p "$F"; cd "$F"

[ -s proteins.faa ] || mm annot agat_sp_extract_sequences.pl -g "$GFF" -f "$GM" -p -o proteins.faa >/dev/null
# InterProScan rejects '*'; strip terminal stops rather than dropping proteins.
[ -s proteins.clean.faa ] || sed 's/\*$//' proteins.faa | mm annot seqkit seq -g > proteins.clean.faa
echo "[12] $(grep -c '>' proteins.clean.faa) proteins"

# ---- InterProScan (domains, GO, pathways) -----------------------------------
if [ ! -s interproscan.tsv ]; then
  IPR="$(ls -d "$DB"/interproscan-* 2>/dev/null | tail -1)"
  if [ -n "$IPR" ] && [ -x "$IPR/interproscan.sh" ]; then
    "$IPR/interproscan.sh" -i proteins.clean.faa -f TSV,GFF3 -b interproscan \
      -goterms -pa -iprlookup -cpu "$T" -T "$FASTTMP/ipr"
  else
    echo "WARNING: InterProScan not installed under $DB - skipping" >&2
    : > interproscan.tsv
  fi
fi

# ---- eggNOG-mapper (orthology -> GO, KEGG, COG) -----------------------------
if [ ! -s eggnog.emapper.annotations ]; then
  if [ -s "$DB/eggnog/eggnog.db" ]; then
    mm func emapper.py -i proteins.clean.faa -o eggnog --output_dir . \
      --data_dir "$DB/eggnog" -m diamond --cpu "$T" --override \
      --tax_scope 7898 --go_evidence non-electronic
  else
    echo "WARNING: eggNOG database missing under $DB/eggnog - skipping" >&2
    : > eggnog.emapper.annotations
  fi
fi

# ---- Swiss-Prot best hits (gene symbols and product names) ------------------
if [ ! -s swissprot.tsv ]; then
  # No download here: compute nodes cannot reach ftp.ebi.ac.uk / NCBI on this
  # cluster. Stage 03 stages Swiss-Prot on the login node.
  if [ ! -s "$DB/uniprot_sprot.fasta.gz" ]; then
    echo "WARNING: $DB/uniprot_sprot.fasta.gz not staged - gene symbols and" >&2
    echo "         product names will fall back to eggNOG/InterPro only." >&2
    echo "         Run scripts/03_fetch_evidence.sh on the login node." >&2
  fi
  if [ -s "$DB/uniprot_sprot.fasta.gz" ]; then
    [ -s sprot.dmnd ] || zcat "$DB/uniprot_sprot.fasta.gz" | mm func diamond makedb --in - -d sprot -p "$T"
    mm func diamond blastp -q proteins.clean.faa -d sprot -p "$T" \
      --very-sensitive --max-target-seqs 1 --evalue 1e-10 \
      --outfmt 6 qseqid sseqid stitle pident length evalue bitscore qcovhsp \
      -o swissprot.tsv
  else
    : > swissprot.tsv
  fi
fi

# ---- merge into the GFF3 and a flat table -----------------------------------
python3 "$LIB/merge_functional.py" \
  --gff3 "$GFF" \
  --interproscan interproscan.tsv \
  --eggnog eggnog.emapper.annotations \
  --swissprot swissprot.tsv \
  --out-gff3 "$F/${ASM_NAME}.annotation.${ANNOT_VERSION}.functional.gff3" \
  --out-tsv  "$F/${ASM_NAME}.functional_annotation.tsv" \
  --summary  "$F/functional_summary.txt"
cat "$F/functional_summary.txt"

done_stamp 12_functional
echo "[12] done"
