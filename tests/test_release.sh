#!/usr/bin/env bash
# Covers rename_gff.py, ncrna_to_gff3.py, merge_functional.py, annotation_stats.py.
set -euo pipefail
cd "$(dirname "$0")"
LIB="../lib"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
fail=0
check() { if [ "$2" = "$3" ]; then echo "  PASS  $1"; else echo "  FAIL  $1: want '$2' got '$3'"; fail=1; fi; }

# ---------------------------------------------------------------- rename ----
cat > "$T/map.tsv" <<'M'
local_name	ncbi_seq_name	genbank_accn	length
ctg000580_np1212	contig_001	BAAFLC010000001.1	17652855
ctg001050_np1212	contig_002	BAAFLC010000002.1	16605418
M
cat > "$T/ann.gff3" <<'G'
##gff-version 3
##sequence-region BAAFLC010000001.1 1 17652855
BAAFLC010000001.1	BRAKER3	gene	1000	3000	.	+	.	ID=SFUS_G000001
BAAFLC010000001.1	BRAKER3	mRNA	1000	3000	.	+	.	ID=SFUS_G000001.1;Parent=SFUS_G000001
BAAFLC010000001.1	BRAKER3	exon	1000	1600	.	+	.	ID=SFUS_G000001.1.exon1;Parent=SFUS_G000001.1
BAAFLC010000001.1	BRAKER3	exon	2400	3000	.	+	.	ID=SFUS_G000001.1.exon2;Parent=SFUS_G000001.1
BAAFLC010000001.1	BRAKER3	CDS	1000	1600	.	+	0	ID=SFUS_G000001.1.cds1;Parent=SFUS_G000001.1
BAAFLC010000001.1	BRAKER3	CDS	2400	3000	.	+	0	ID=SFUS_G000001.1.cds2;Parent=SFUS_G000001.1
BAAFLC010000002.1	BRAKER3	gene	5000	9000	.	-	.	ID=SFUS_G000002
BAAFLC010000002.1	BRAKER3	mRNA	5000	9000	.	-	.	ID=SFUS_G000002.1;Parent=SFUS_G000002
BAAFLC010000002.1	BRAKER3	exon	5000	6000	.	-	.	ID=SFUS_G000002.1.exon1;Parent=SFUS_G000002.1
BAAFLC010000002.1	BRAKER3	exon	8000	9000	.	-	.	ID=SFUS_G000002.1.exon2;Parent=SFUS_G000002.1
BAAFLC010000002.1	BRAKER3	CDS	5000	6000	.	-	0	ID=SFUS_G000002.1.cds1;Parent=SFUS_G000002.1
BAAFLC010000002.1	BRAKER3	CDS	8000	9000	.	-	0	ID=SFUS_G000002.1.cds2;Parent=SFUS_G000002.1
G

echo "=== rename_gff ==="
python3 "$LIB/rename_gff.py" --gff "$T/ann.gff3" --map "$T/map.tsv" \
  --from genbank_accn --to local_name --out "$T/alt/ann.local.gff3"
check "seqid renamed to original contig" "ctg000580_np1212" \
  "$(awk -F'\t' '!/^#/{print $1; exit}' "$T/alt/ann.local.gff3")"
check "sequence-region also renamed" "1" \
  "$(grep -c '##sequence-region ctg000580_np1212' "$T/alt/ann.local.gff3")"
check "coordinates unchanged" "1000" \
  "$(awk -F'\t' '$3=="gene"{print $4; exit}' "$T/alt/ann.local.gff3")"
check "line count preserved" \
  "$(wc -l < "$T/ann.gff3")" "$(wc -l < "$T/alt/ann.local.gff3")"
echo "  (unmapped seqid should abort)"
printf 'ctgXXX\tBRAKER3\tgene\t1\t9\t.\t+\t.\tID=x\n' > "$T/bad.gff3"
if python3 "$LIB/rename_gff.py" --gff "$T/bad.gff3" --map "$T/map.tsv" \
     --from genbank_accn --to local_name --out "$T/bad.out.gff3" 2>/dev/null; then
  echo "  FAIL  unmapped seqid was not rejected"; fail=1
else
  echo "  PASS  unmapped seqid rejected"
fi

# ----------------------------------------------------------------- ncRNA ----
cat > "$T/trna.out" <<'R'
Sequence		tRNA	Bounds	tRNA	Anti	Intron Bounds	Inf
Name		tRNA #	Begin	End	Type	Codon	Begin	End	Score
--------	------	-----	------	----	-----	-----	----	------
BAAFLC010000001.1	1	5000	5072	Ala	AGC	0	0	78.2
BAAFLC010000001.1	2	9000	8928	Gly	GCC	0	0	65.4
R
cat > "$T/rrna.gff3" <<'R'
##gff-version 3
BAAFLC010000001.1	barrnap:0.9	rRNA	20000	21800	1e-30	+	.	Name=18S_rRNA;product=18S ribosomal RNA
R
cat > "$T/rfam.tblout" <<'R'
#idx target        accession  query              accession clan  mdl mdl-from mdl-to seq-from seq-to strand trunc pass gc bias score E-value inc olp
1    tRNA          RF00005    BAAFLC010000001.1  -         CL00001 cm 1 72 5000 5072 + no 1 0.55 0.0 60.1 1e-12 !  ^
2    U2            RF00004    BAAFLC010000001.1  -         CL00001 cm 1 190 40000 40190 + no 1 0.48 0.0 95.3 1e-25 ! ^
3    mir-21        RF00651    BAAFLC010000002.1  -         -       cm 1 72 12000 12072 - no 1 0.51 0.0 55.2 1e-10 ! ^
4    SSU_rRNA_euk  RF01960    BAAFLC010000001.1  -         CL00111 cm 1 1800 20050 21700 + no 1 0.52 0.0 900.0 0 ! ^
R
echo
echo "=== ncrna_to_gff3 ==="
python3 "$LIB/ncrna_to_gff3.py" --trnascan "$T/trna.out" --rrna "$T/rrna.gff3" \
  --rfam "$T/rfam.tblout" --out "$T/ncrna.gff3" --summary "$T/ncrna.tsv"
echo; cat "$T/ncrna.tsv"
check "tRNA count = 2"  "2" "$(awk -F'\t' '$3=="tRNA"'  "$T/ncrna.gff3" | wc -l)"
check "rRNA count = 1 (Rfam dup removed)" "1" "$(awk -F'\t' '$3=="rRNA"' "$T/ncrna.gff3" | wc -l)"
check "snRNA kept"      "1" "$(awk -F'\t' '$3=="snRNA"' "$T/ncrna.gff3" | wc -l)"
check "miRNA kept"      "1" "$(awk -F'\t' '$3=="miRNA"' "$T/ncrna.gff3" | wc -l)"
check "minus-strand tRNA normalised" "8928" \
  "$(awk -F'\t' '$3=="tRNA" && $7=="-"{print $4}' "$T/ncrna.gff3")"
check "Rfam Dbxref present" "2" "$(grep -c 'Dbxref=RFAM:' "$T/ncrna.gff3")"

# MirMachine is the specialist for miRNA, so its call must displace the
# overlapping generic Rfam hit rather than both being emitted.
cat > "$T/mirmachine.gff" <<'M'
##gff-version 3
BAAFLC010000002.1	MirMachine	miRNA	11990	12080	95.4	-	.	Name=Mir-21-P1;ID=mir21
M
python3 "$LIB/ncrna_to_gff3.py" --trnascan "$T/trna.out" --rrna "$T/rrna.gff3" \
  --rfam "$T/rfam.tblout" --mirmachine "$T/mirmachine.gff" \
  --out "$T/ncrna2.gff3" --summary "$T/ncrna2.tsv"
check "miRNA still exactly one"        "1" "$(awk -F'\t' '$3=="miRNA"' "$T/ncrna2.gff3" | wc -l)"
check "MirMachine call won over Rfam"  "1" "$(grep -c 'Name=Mir-21-P1' "$T/ncrna2.gff3")"
check "Rfam mir-21 suppressed"         "0" "$(grep -c 'RFAM:RF00651' "$T/ncrna2.gff3")"
check "other classes unaffected"       "1" "$(awk -F'\t' '$3=="snRNA"' "$T/ncrna2.gff3" | wc -l)"

# ------------------------------------------------------------- functional ---
printf 'SFUS_G000001.1\tmd5\t400\tPfam\tPF00089\tTrypsin\t20\t240\t1e-40\tT\t01-01-2026\tIPR001254\tSerine protease, trypsin domain\tGO:0006508|GO:0004252\tKEGG: 00230\n' > "$T/ipr.tsv"
{
  printf '#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\n'
  printf 'SFUS_G000001.1\t7955.X\t1e-99\t400\t28IPW@33208,COG0265@1\t33208\tO\tSerine protease\tprss1\tGO:0006508,GO:0008233\t3.4.21.4\tko:K01312\tko04974\n'
  printf 'SFUS_G000002.1\t7955.Y\t1e-50\t200\t2ABCD@33208\t33208\tS\t-\t-\t-\t-\t-\t-\n'
} > "$T/egg.tsv"
printf 'SFUS_G000001.1\tsp|P07477|TRY1_HUMAN\tsp|P07477|TRY1_HUMAN Trypsin-1 OS=Homo sapiens OX=9606 GN=PRSS1 PE=1 SV=1\t72.5\t240\t1e-88\t320.0\t95.0\n' > "$T/sprot.tsv"

echo
echo "=== merge_functional ==="
python3 "$LIB/merge_functional.py" --gff3 "$T/ann.gff3" \
  --interproscan "$T/ipr.tsv" --eggnog "$T/egg.tsv" --swissprot "$T/sprot.tsv" \
  --out-gff3 "$T/func.gff3" --out-tsv "$T/func.tsv" --summary "$T/func.txt"
echo; awk -F'\t' '$3=="mRNA"{print $9}' "$T/func.gff3"
check "Swiss-Prot description wins"  "1" "$(grep -c 'product=Trypsin-1' "$T/func.gff3")"
check "gene symbol from Swiss-Prot"  "1" "$(grep -c 'gene_symbol=PRSS1' "$T/func.gff3")"
check "GO terms merged and unique"   "GO:0004252,GO:0006508,GO:0008233" \
  "$(awk -F'\t' '$3=="mRNA"' "$T/func.gff3" | grep -o 'Ontology_term=[^;]*' | head -1 | cut -d= -f2)"
check "InterPro + UniProt + KEGG xrefs" "1" \
  "$(awk -F'\t' '$3=="mRNA"' "$T/func.gff3" | grep -c 'Dbxref=InterPro:IPR001254,UniProtKB:P07477,KEGG:K01312')"
check "unannotated model marked hypothetical" "1" \
  "$(grep -c 'product=hypothetical protein' "$T/func.gff3")"
check "CDS lines untouched" "4" "$(awk -F'\t' '$3=="CDS"' "$T/func.gff3" | wc -l)"
check "TSV has 2 data rows" "2" "$(($(wc -l < "$T/func.tsv") - 1))"

# -------------------------------------------------------------- statistics --
printf '>SFUS_G000001.1\nMKAILVVLLYTFATANAD*\n>SFUS_G000002.1\nMKAILVVLLYTFATANADT\n' > "$T/prot.faa"
echo
echo "=== annotation_stats (gates should FAIL on a 2-gene toy set) ==="
if python3 "$LIB/annotation_stats.py" --gff3 "$T/func.gff3" --proteins "$T/prot.faa" \
     --out "$T/stats.tsv" --report "$T/stats.txt" --check-thresholds > "$T/stats.log" 2>&1; then
  echo "  FAIL  gates passed on a toy set (should not)"; fail=1
else
  echo "  PASS  gates correctly failed"
fi
grep -E '^  (PASS|FAIL)' "$T/stats.txt" | sed 's/^/  /'
check "gene count metric = 2" "2" "$(awk -F'\t' '$1=="genes"{print $2}' "$T/stats.tsv")"
check "mean exons metric = 2.0" "2" "$(awk -F'\t' '$1=="mean_exons_per_transcript"{print $2}' "$T/stats.tsv")"
check "no internal stops counted" "0" \
  "$(awk -F'\t' '$1=="internal_stop_proteins"{print $2}' "$T/stats.tsv")"

echo
[ "$fail" -eq 0 ] && echo "ALL TESTS PASSED" || { echo "TESTS FAILED"; exit 1; }
