#!/usr/bin/env bash
# Stage 03 - fetch all external evidence. Login node (network).
# Skips anything already staged by scripts/mirror_fetch.sh.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh
mkdir -p "$DB" "$STAGED" "$WORK/rnaseq/raw"

have() { [ -s "$1" ]; }
from_staged() {  # from_staged <filename> <dest> - use the mirrored copy if present
  if have "$STAGED/$1"; then echo "[03] using staged $1"; cp -n "$STAGED/$1" "$2"; return 0; fi
  return 1
}
get() {  # get <url> <dest>
  have "$2" && { echo "[03] have $(basename "$2")"; return; }
  from_staged "$(basename "$2")" "$2" && return
  echo "[03] downloading $(basename "$2")"
  curl -fsSL --retry 5 --retry-delay 5 -C - -o "$2.part" "$1" && mv "$2.part" "$2"
}

# ---- 1. RNA-seq (conspecific only, per project decision) --------------------
cd "$WORK/rnaseq/raw"
for acc in $RNASEQ_RUNS; do
  if ! have "${acc}_1.fastq.gz"; then
    echo "[03] prefetch+dump $acc"
    mm rnaseq prefetch --max-size 100G -O . "$acc"
    mm rnaseq fasterq-dump --split-3 --threads 8 -O . "$acc"
    mm rnaseq pigz -p 8 "${acc}"_*.fastq
  fi
done
ls -la

# ---- 2. Reference proteomes (NCBI - reachable) ------------------------------
P="$DB/proteomes"; mkdir -p "$P"; cd "$P"
python3 "$LIB/resolve_proteomes.py" --out ref_proteomes.tsv
tail -n +2 ref_proteomes.tsv | while IFS=$'\t' read -r sp acc name lvl annot url; do
  short="$(echo "$sp" | awk '{print substr($1,1,1)"_"$2}')"
  get "$url" "${short}.faa.gz"
done

# ---- 3. OrthoDB v12 Vertebrata (BLOCKED host - expect staged copy) ----------
get "$ORTHODB_URL" "$DB/Vertebrata.fa.gz"

# ---- 4. Dfam 4.0 (BLOCKED host - expect staged copy) -----------------------
D="$DB/dfam"; mkdir -p "$D"
for f in $DFAM_FILES; do get "$DFAM_BASE/$f" "$D/$f"; done
for f in $DFAM_FILES; do [ -s "$D/${f%.gz}" ] || gzip -dk "$D/$f"; done

# ---- 4b. NCBI taxonomy (to cut a clade partition out of OrthoDB) ------------
if [ "${ORTHODB_CLADE_TAXID:-0}" != "0" ]; then
  TX="$DB/taxonomy"; mkdir -p "$TX"
  get "$TAXDUMP_URL" "$TX/taxdump.tar.gz"
  [ -s "$TX/nodes.dmp" ] || tar xzf "$TX/taxdump.tar.gz" -C "$TX" nodes.dmp names.dmp
fi

# ---- 5. Rfam, BUSCO lineage -------------------------------------------------
mkdir -p "$DB/rfam" "$DB/busco"
get "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"     "$DB/rfam/Rfam.cm.gz"
get "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"    "$DB/rfam/Rfam.clanin"
[ -s "$DB/rfam/Rfam.cm" ] || gzip -dk "$DB/rfam/Rfam.cm.gz"

echo "[03] evidence summary"; du -sh "$DB"/* "$WORK/rnaseq/raw" 2>/dev/null
done_stamp 03_fetch_evidence
