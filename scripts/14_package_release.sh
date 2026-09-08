#!/usr/bin/env bash
# Stage 14 - assemble the public release bundle. Login node; fast.
#
# Produces coordinate-identical files under three naming schemes so nobody has
# to rename anything: GenBank accessions (matches the NCBI download), NCBI
# sequence names, and the original assembly contig names.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

V="$ANNOT_VERSION"
OUT="$REL/${ASM_NAME}.annotation.${V}"
mkdir -p "$OUT"
MAP="$WORK/genome/contig_name_map.tsv"
GFF="$WORK/functional/${ASM_NAME}.annotation.${V}.functional.gff3"
[ -s "$GFF" ] || GFF="$WORK/filter/${ASM_NAME}.annotation.${V}.gff3"
need "$GFF" "$MAP"

cp "$GFF" "$OUT/${ASM_NAME}.annotation.${V}.gff3"
cd "$OUT"

# ---- GTF (what the original request asked for) ------------------------------
mm annot agat_convert_sp_gff2gtf.pl -g "${ASM_NAME}.annotation.${V}.gff3" \
  -o "${ASM_NAME}.annotation.${V}.gtf" >/dev/null
# Round-trip check: the GTF must describe the same CDS as the GFF3.
A=$(awk -F'\t' '$3=="CDS"' "${ASM_NAME}.annotation.${V}.gff3" | wc -l)
B=$(awk -F'\t' '$3=="CDS"' "${ASM_NAME}.annotation.${V}.gtf"  | wc -l)
[ "$A" -eq "$B" ] || { echo "FATAL: GFF3 has $A CDS but GTF has $B" >&2; exit 1; }
echo "[14] GTF round-trip OK ($A CDS features)"

# ---- alternative sequence naming --------------------------------------------
for scheme in ncbi_seq_name local_name; do
  python3 "$LIB/rename_gff.py" \
    --gff "${ASM_NAME}.annotation.${V}.gff3" --map "$MAP" \
    --from genbank_accn --to "$scheme" \
    --out "alt_names/${ASM_NAME}.annotation.${V}.${scheme}.gff3"
done

# ---- sequence products -------------------------------------------------------
for f in proteins.faa cds.fna transcripts.fna; do
  [ -s "$WORK/qc/$f" ] && cp "$WORK/qc/$f" "${ASM_NAME}.${f}"
done
cp "$WORK/ncrna/${ASM_NAME}.ncRNA.gff3"            . 2>/dev/null || true
cp "$WORK/repeats/final/${ASM_NAME}.repeats.gff3"  . 2>/dev/null || true
cp "$WORK/repeats/final/repeat_summary.tsv"        "${ASM_NAME}.repeat_summary.tsv" 2>/dev/null || true
cp "$WORK/functional/${ASM_NAME}.functional_annotation.tsv" . 2>/dev/null || true
cp "$WORK/filter/evidence_support.tsv"             "${ASM_NAME}.evidence_support.tsv" 2>/dev/null || true
cp "$WORK/filter/id_map.tsv"                       "${ASM_NAME}.id_map.tsv" 2>/dev/null || true
cp "$MAP" contig_name_map.tsv
cp "$WORK/repeats/final/genome.softmasked.fna"     "${ASM_NAME}.genomic.softmasked.fna" 2>/dev/null || true

# ---- provenance --------------------------------------------------------------
python3 "$LIB/write_methods.py" --config ../../config.sh --work "$WORK" \
  --out METHODS.md --versions software_versions.tsv 2>/dev/null || true
cp "$REPO/docs/RELEASE_README.md" README.md 2>/dev/null || true
cp "$REPO/LICENSE" . 2>/dev/null || true

# ---- compress and checksum ---------------------------------------------------
for f in *.gff3 *.gtf *.faa *.fna *.tsv; do
  [ -s "$f" ] && [ "${f##*.}" != "gz" ] && mm prep pigz -f -k "$f" 2>/dev/null || true
done
find . -type f ! -name MANIFEST.md5 -exec md5sum {} + | sort -k2 > MANIFEST.md5

echo; echo "[14] release bundle:"; ls -la; du -sh .
echo; echo "Verify anywhere with:  md5sum -c MANIFEST.md5"
done_stamp 14_package_release
