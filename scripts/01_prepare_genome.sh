#!/usr/bin/env bash
# Stage 01 - genome preparation. Login node; needs ftp.ncbi.nlm.nih.gov only.
#
# Pulls the assembly from GenBank rather than from a local copy, so the release
# is provably built on the same bytes the public downloads. Produces:
#   $WORK/genome/genome.fna              uppercase, GenBank accession headers
#   $WORK/genome/genome.ncbi_masked.fna  NCBI's own softmasking (21.15%) - QC ref
#   $WORK/genome/contig_name_map.tsv     local <-> contig_NNN <-> GenBank accn
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

G="$WORK/genome"; mkdir -p "$G"; cd "$G"
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/042/847/805/${ASM_ACC}_${ASM_NAME}"

fetch() {  # fetch <url> <dest>  - resumable, skips if already complete
  [ -s "$2" ] && { echo "[01] have $2"; return; }
  echo "[01] downloading $(basename "$2")"
  curl -fsSL --retry 5 --retry-delay 5 -C - -o "$2.part" "$1"
  mv "$2.part" "$2"
}

fetch "$BASE/${ASM_ACC}_${ASM_NAME}_genomic.fna.gz"      genomic.fna.gz
fetch "$BASE/${ASM_ACC}_${ASM_NAME}_assembly_report.txt"  assembly_report.txt

# --- verify the download against NCBI's own checksums -----------------------
if curl -fsSL "$BASE/md5checksums.txt" -o md5checksums.txt 2>/dev/null; then
  grep -E 'genomic\.fna\.gz$' md5checksums.txt | sed 's#\./##' > md5.want
  if md5sum -c md5.want 2>/dev/null | grep -q OK; then
    echo "[01] md5 OK against NCBI"
  else
    echo "[01] FATAL: md5 mismatch on genomic.fna.gz" >&2; exit 1
  fi
fi

# --- unmasked working copy (RepeatModeler input) ----------------------------
# Softmasking is case only, so uppercasing losslessly recovers the raw sequence.
if [ ! -s genome.fna ]; then
  mm prep seqkit seq -u -w 60 genomic.fna.gz \
    | mm prep seqkit replace -p '\s.*$' -r '' > genome.fna
fi
mm prep samtools faidx genome.fna

# --- keep NCBI's softmasked version as an independent QC reference ----------
if [ ! -s genome.ncbi_masked.fna ]; then
  mm prep seqkit replace -p '\s.*$' -r '' -w 60 genomic.fna.gz > genome.ncbi_masked.fna
fi

# --- name map + integrity proof ---------------------------------------------
python3 "$LIB/fasta_tools.py" map \
  --fasta genome.fna --report assembly_report.txt --out contig_name_map.tsv

python3 "$LIB/fasta_tools.py" verify \
  --a genome.fna --b genome.ncbi_masked.fna \
  --map contig_name_map.tsv --from genbank_accn --to genbank_accn

echo "--- unmasked working genome ---"
python3 "$LIB/fasta_tools.py" stats --fasta genome.fna | tee stats.genome.txt
echo "--- NCBI softmasked reference ---"
python3 "$LIB/fasta_tools.py" stats --fasta genome.ncbi_masked.fna | tee stats.ncbi_masked.txt

done_stamp 01_prepare_genome
echo "[01] done -> $G"
