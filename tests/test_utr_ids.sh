#!/usr/bin/env bash
# Tests for lib/add_utrs.py (conservative UTR transfer) and lib/assign_ids.py.
set -euo pipefail
cd "$(dirname "$0")"
LIB="../lib"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
fail=0
check() { if [ "$2" = "$3" ]; then echo "  PASS  $1"; else echo "  FAIL  $1: want '$2' got '$3'"; fail=1; fi; }

# gA  + strand, 2 CDS exons, StringTie transcript extends both ends -> UTRs added
# gB  - strand, StringTie transcript extends -> 5'UTR must land DOWNSTREAM (higher coords)
# gC  + strand, StringTie transcript disagrees on splicing -> must be refused
# gD  + strand, no StringTie transcript at all -> passes through untouched
cat > "$T/in.gff3" <<'GFF'
##gff-version 3
ctg1	BRAKER	gene	1000	2000	.	+	.	ID=gA
ctg1	BRAKER	mRNA	1000	2000	.	+	.	ID=gA.t1;Parent=gA
ctg1	BRAKER	exon	1000	1200	.	+	.	ID=gA.t1.e1;Parent=gA.t1
ctg1	BRAKER	exon	1800	2000	.	+	.	ID=gA.t1.e2;Parent=gA.t1
ctg1	BRAKER	CDS	1000	1200	.	+	0	ID=gA.t1.c1;Parent=gA.t1
ctg1	BRAKER	CDS	1800	2000	.	+	0	ID=gA.t1.c2;Parent=gA.t1
ctg1	BRAKER	gene	5000	5500	.	-	.	ID=gB
ctg1	BRAKER	mRNA	5000	5500	.	-	.	ID=gB.t1;Parent=gB
ctg1	BRAKER	exon	5000	5500	.	-	.	ID=gB.t1.e1;Parent=gB.t1
ctg1	BRAKER	CDS	5000	5500	.	-	0	ID=gB.t1.c1;Parent=gB.t1
ctg1	BRAKER	gene	9000	10000	.	+	.	ID=gC
ctg1	BRAKER	mRNA	9000	10000	.	+	.	ID=gC.t1;Parent=gC
ctg1	BRAKER	exon	9000	9300	.	+	.	ID=gC.t1.e1;Parent=gC.t1
ctg1	BRAKER	exon	9700	10000	.	+	.	ID=gC.t1.e2;Parent=gC.t1
ctg1	BRAKER	CDS	9000	9300	.	+	0	ID=gC.t1.c1;Parent=gC.t1
ctg1	BRAKER	CDS	9700	10000	.	+	0	ID=gC.t1.c2;Parent=gC.t1
ctg1	BRAKER	gene	20000	20500	.	+	.	ID=gD
ctg1	BRAKER	mRNA	20000	20500	.	+	.	ID=gD.t1;Parent=gD
ctg1	BRAKER	exon	20000	20500	.	+	.	ID=gD.t1.e1;Parent=gD.t1
ctg1	BRAKER	CDS	20000	20500	.	+	0	ID=gD.t1.c1;Parent=gD.t1
GFF

cat > "$T/st.gtf" <<'GTF'
ctg1	StringTie	transcript	800	2300	1000	+	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.1";
ctg1	StringTie	exon	800	1200	1000	+	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.1";
ctg1	StringTie	exon	1800	2300	1000	+	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.1";
ctg1	StringTie	transcript	4800	5900	1000	-	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1";
ctg1	StringTie	exon	4800	5900	1000	-	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1";
ctg1	StringTie	transcript	8800	10200	1000	+	.	gene_id "MSTRG.3"; transcript_id "MSTRG.3.1";
ctg1	StringTie	exon	8800	9300	1000	+	.	gene_id "MSTRG.3"; transcript_id "MSTRG.3.1";
ctg1	StringTie	exon	9500	10200	1000	+	.	gene_id "MSTRG.3"; transcript_id "MSTRG.3.1";
GTF

echo "=== add_utrs ==="
python3 "$LIB/add_utrs.py" --gff3 "$T/in.gff3" --transcripts "$T/st.gtf" \
  --out "$T/utr.gff3" --report "$T/utr.tsv"
echo; column -t "$T/utr.tsv"

st() { awk -F'\t' -v t="$1" '$1==t{print $5}' "$T/utr.tsv"; }
u5() { awk -F'\t' -v t="$1" '$1==t{print $3}' "$T/utr.tsv"; }
u3() { awk -F'\t' -v t="$1" '$1==t{print $4}' "$T/utr.tsv"; }

echo; echo "--- assertions: UTR transfer ---"
check "gA accepted"                  "ok" "$(st gA.t1)"
check "gA 5'UTR = 200 bp (800-999)"  "200" "$(u5 gA.t1)"
check "gA 3'UTR = 300 bp (2001-2300)" "300" "$(u3 gA.t1)"
check "gB accepted"                  "ok" "$(st gB.t1)"
# minus strand: 5'UTR is the HIGH-coordinate side (5501-5900 = 400 bp)
check "gB 5'UTR = 400 bp (high side)" "400" "$(u5 gB.t1)"
check "gB 3'UTR = 200 bp (low side)"  "200" "$(u3 gB.t1)"
check "gC refused (splice mismatch)" "no_qualifying_transcript" "$(st gC.t1)"
check "gD refused (no transcript)"   "no_qualifying_transcript" "$(st gD.t1)"
check "gA gene start extended to 800" "800" \
  "$(awk -F'\t' '$3=="gene" && $9 ~ /ID=gA/ {print $4}' "$T/utr.gff3")"
check "gC unchanged, no UTR feature"  "0" \
  "$(awk -F'\t' '$9 ~ /Parent=gC.t1/ && $3 ~ /UTR/' "$T/utr.gff3" | wc -l)"
check "CDS count preserved (6)"       "6" \
  "$(awk -F'\t' '$3=="CDS"' "$T/utr.gff3" | wc -l)"

echo; echo "=== assign_ids ==="
python3 "$LIB/assign_ids.py" --gff3 "$T/utr.gff3" --prefix SFUS \
  --out "$T/final.gff3" --map "$T/idmap.tsv" --source BRAKER3
echo; head -14 "$T/final.gff3"

echo; echo "--- assertions: IDs ---"
check "first gene is SFUS_G000001" "SFUS_G000001" \
  "$(awk -F'\t' '$3=="gene"{match($9,/ID=[^;]+/); print substr($9,RSTART+3,RLENGTH-3); exit}' "$T/final.gff3")"
check "4 genes numbered"           "4" "$(awk -F'\t' '$3=="gene"' "$T/final.gff3" | wc -l)"
check "transcript is .1 suffixed"  "SFUS_G000001.1" \
  "$(awk -F'\t' '$3=="mRNA"{match($9,/ID=[^;]+/); print substr($9,RSTART+3,RLENGTH-3); exit}' "$T/final.gff3")"
check "old gene id preserved"      "1" \
  "$(grep -cE 'old_locus_tag=gA(;|$)' "$T/final.gff3")"
check "old transcript id preserved" "1" \
  "$(grep -cE 'old_locus_tag=gA\.t1(;|$)' "$T/final.gff3")"
check "ID is first attribute"      "0" \
  "$(awk -F'\t' '!/^#/ && $9 !~ /^ID=/' "$T/final.gff3" | wc -l)"
check "id map has header + rows"   "9" "$(wc -l < "$T/idmap.tsv")"
check "source column rewritten"    "0" "$(awk -F'\t' '$2=="BRAKER"' "$T/final.gff3" | wc -l)"
check "no orphan Parent refs"      "0" \
  "$(python3 - "$T/final.gff3" <<'PY'
import sys, re
ids=set(); bad=0; rows=[]
for line in open(sys.argv[1]):
    if line.startswith('#'): continue
    f=line.rstrip('\n').split('\t')
    if len(f)!=9: continue
    m=re.search(r'ID=([^;]+)',f[8]);  ids.add(m.group(1)) if m else None
    rows.append(f[8])
for a in rows:
    m=re.search(r'Parent=([^;]+)',a)
    if m and any(p not in ids for p in m.group(1).split(',')): bad+=1
print(bad)
PY
)"

echo
[ "$fail" -eq 0 ] && echo "ALL TESTS PASSED" || { echo "TESTS FAILED"; exit 1; }
