#!/usr/bin/env bash
# Fetch resources the cluster firewall blocks. Run on ANY internet-connected
# machine (laptop, workstation, WSL), then rsync the output to $DB/staged/.
#
#   ./scripts/mirror_fetch.sh /some/local/staging
#   rsync -avP /some/local/staging/ cpu201:/projects/.../aigo_genome/db/staged/
#
# Total payload ~6 GB, dominated by OrthoDB Vertebrata.
set -euo pipefail
DEST="${1:?usage: mirror_fetch.sh <staging-dir>}"
mkdir -p "$DEST"; cd "$DEST"

ORTHODB_URL="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Vertebrata.fa.gz"
DFAM_BASE="https://www.dfam.org/releases/current/families/FamDB"
DFAM_FILES="dfam40.0.h5.gz dfam40.curated.consensus.0.h5.gz"
BUSCO_URL="https://busco-data.ezlab.org/v5/data/lineages/"
RFAM_BASE="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT"

get() {
  local url="$1" dst="$2"
  [ -s "$dst" ] && { echo "  have $dst"; return; }
  echo "  fetching $dst"
  curl -fL --retry 5 --retry-delay 5 -C - --progress-bar -o "$dst.part" "$url"
  mv "$dst.part" "$dst"
}

echo "[mirror] OrthoDB v12 Vertebrata (~5.9 GB)"
get "$ORTHODB_URL" Vertebrata.fa.gz

echo "[mirror] Dfam 4.0 root + curated consensus (~90 MB)"
for f in $DFAM_FILES; do
  get "$DFAM_BASE/$f" "$f"
  get "$DFAM_BASE/$f.md5" "$f.md5" || true
done
for f in $DFAM_FILES; do
  if [ -s "$f.md5" ]; then
    md5sum -c "$f.md5" || echo "  WARN: md5 mismatch on $f"
  fi
done

# busco-data.ezlab.org does NOT serve directory listings (200 with an empty
# body), so the dated tarball name has to come from the version index. That
# index also carries the md5, which is worth checking on a 200 MB download.
# odb10 keeps the new BUSCO numbers directly comparable with the 2024 paper;
# odb12 is the current lineage, so both are staged.
echo "[mirror] BUSCO lineages"
if curl -fsSL --retry 3 "https://busco-data.ezlab.org/v5/data/file_versions.tsv" \
     -o file_versions.tsv; then
  for lineage in actinopterygii_odb10 actinopterygii_odb12; do
    read -r d sum < <(awk -F'\t' -v L="$lineage" '$1==L {print $2, $3; exit}' file_versions.tsv)
    if [ -z "${d:-}" ]; then
      echo "  WARN: $lineage not present in the version index - skipping"
      continue
    fi
    tarball="${lineage}.${d}.tar.gz"
    get "${BUSCO_URL}${tarball}" "$tarball"      # -L: the path 301-redirects
    if [ -n "${sum:-}" ] && [ -s "$tarball" ]; then
      got=$(md5sum "$tarball" | awk '{print $1}')
      if [ "$got" = "$sum" ]; then
        echo "  md5 OK   $tarball"
      else
        echo "  WARN: md5 mismatch on $tarball (want $sum got $got)"
      fi
    fi
  done
else
  echo "  WARN: could not reach the BUSCO version index - skipping lineages"
fi

echo "[mirror] Rfam covariance models"
get "$RFAM_BASE/Rfam.cm.gz"  Rfam.cm.gz
get "$RFAM_BASE/Rfam.clanin" Rfam.clanin

echo; echo "[mirror] staged payload:"; du -sh .; ls -la
cat <<TXT

Next:
  rsync -avP --partial "$DEST"/ <user>@cpu201:/projects/yates_lab_hpc/sam/smwambu1/aigo_genome/db/staged/
TXT
