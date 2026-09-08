#!/usr/bin/env bash
# Unit test for lib/filter_models.py - exercises every drop branch plus the
# partial case where one transcript of a multi-transcript gene is removed.
set -euo pipefail
cd "$(dirname "$0")"
LIB="../lib"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT

cat > "$T/in.gff3" <<'GFF'
##gff-version 3
##sequence-region BAAFLC010000001.1 1 17652855
BAAFLC010000001.1	AUGUSTUS	gene	1000	3000	.	+	.	ID=g1
BAAFLC010000001.1	AUGUSTUS	mRNA	1000	3000	.	+	.	ID=g1.t1;Parent=g1
BAAFLC010000001.1	AUGUSTUS	exon	1000	1600	.	+	.	ID=g1.t1.e1;Parent=g1.t1
BAAFLC010000001.1	AUGUSTUS	exon	2400	3000	.	+	.	ID=g1.t1.e2;Parent=g1.t1
BAAFLC010000001.1	AUGUSTUS	CDS	1000	1600	.	+	0	ID=g1.t1.c1;Parent=g1.t1
BAAFLC010000001.1	AUGUSTUS	CDS	2400	3000	.	+	0	ID=g1.t1.c2;Parent=g1.t1
BAAFLC010000001.1	AUGUSTUS	gene	5000	6000	.	-	.	ID=g2
BAAFLC010000001.1	AUGUSTUS	mRNA	5000	6000	.	-	.	ID=g2.t1;Parent=g2
BAAFLC010000001.1	AUGUSTUS	exon	5000	6000	.	-	.	ID=g2.t1.e1;Parent=g2.t1
BAAFLC010000001.1	AUGUSTUS	CDS	5000	6000	.	-	0	ID=g2.t1.c1;Parent=g2.t1
BAAFLC010000001.1	AUGUSTUS	gene	8000	8600	.	+	.	ID=g3
BAAFLC010000001.1	AUGUSTUS	mRNA	8000	8600	.	+	.	ID=g3.t1;Parent=g3
BAAFLC010000001.1	AUGUSTUS	exon	8000	8600	.	+	.	ID=g3.t1.e1;Parent=g3.t1
BAAFLC010000001.1	AUGUSTUS	CDS	8000	8600	.	+	0	ID=g3.t1.c1;Parent=g3.t1
BAAFLC010000002.1	AUGUSTUS	gene	100	200	.	+	.	ID=g4
BAAFLC010000002.1	AUGUSTUS	mRNA	100	200	.	+	.	ID=g4.t1;Parent=g4
BAAFLC010000002.1	AUGUSTUS	exon	100	200	.	+	.	ID=g4.t1.e1;Parent=g4.t1
BAAFLC010000002.1	AUGUSTUS	CDS	100	200	.	+	0	ID=g4.t1.c1;Parent=g4.t1
BAAFLC010000002.1	AUGUSTUS	gene	9000	12000	.	+	.	ID=g5
BAAFLC010000002.1	AUGUSTUS	mRNA	9000	12000	.	+	.	ID=g5.t1;Parent=g5
BAAFLC010000002.1	AUGUSTUS	exon	9000	9800	.	+	.	ID=g5.t1.e1;Parent=g5.t1
BAAFLC010000002.1	AUGUSTUS	exon	11000	12000	.	+	.	ID=g5.t1.e2;Parent=g5.t1
BAAFLC010000002.1	AUGUSTUS	CDS	9000	9800	.	+	0	ID=g5.t1.c1;Parent=g5.t1
BAAFLC010000002.1	AUGUSTUS	CDS	11000	12000	.	+	0	ID=g5.t1.c2;Parent=g5.t1
BAAFLC010000002.1	AUGUSTUS	mRNA	9000	12000	.	+	.	ID=g5.t2;Parent=g5
BAAFLC010000002.1	AUGUSTUS	exon	9000	9800	.	+	.	ID=g5.t2.e1;Parent=g5.t2
BAAFLC010000002.1	AUGUSTUS	CDS	9000	9800	.	+	0	ID=g5.t2.c1;Parent=g5.t2
BAAFLC010000002.1	AUGUSTUS	gene	20000	23000	.	+	.	ID=g6
BAAFLC010000002.1	AUGUSTUS	mRNA	20000	23000	.	+	.	ID=g6.t1;Parent=g6
BAAFLC010000002.1	AUGUSTUS	exon	20000	20900	.	+	.	ID=g6.t1.e1;Parent=g6.t1
BAAFLC010000002.1	AUGUSTUS	exon	22000	23000	.	+	.	ID=g6.t1.e2;Parent=g6.t1
BAAFLC010000002.1	AUGUSTUS	CDS	20000	20900	.	+	0	ID=g6.t1.c1;Parent=g6.t1
BAAFLC010000002.1	AUGUSTUS	CDS	22000	23000	.	+	0	ID=g6.t1.c2;Parent=g6.t1
GFF

# g1 strong homology + expression                 -> keep
# g2 no homology, no expression, repeaty, no lift -> drop repeat_derived_unsupported
# g3 internal stop in protein                     -> drop internal_stop
# g4 CDS 101 bp (< 150)                           -> drop cds_too_short
# g5 t1 expressed -> keep; t2 unsupported -> drop, gene survives
# g6 NO homology, NO expression, but lifted from the congener -> keep.
#    This is the v2 behaviour change: under v1 rules g6 would have been dropped
#    as "no_evidence". It is the case the congener axis exists to rescue.

printf '>g1.t1\nMKAILVVLLYTFATANADTLLILGDSLSAGYRMSASAAWPALLNDKWQSKTSVVNASISGDTSQQGLARLPALLKQHQPRWVLVELGGNDGLRGFQPQQTEQTLRQILQDVKAANAEPLLMQIRLPANYGRRYNEAFSAIYPKLAKEFDVPLLPFFMEEVYLKPQWMQDDGIHPNRDAQPFIADWMAKQLQPLVNHDS*\n' > "$T/proteins.faa"
printf '>g2.t1\nMAAAAAAAAAAGGGGGGGGGGKKKKKKKKKKPPPPPPPPPPWWWWWWWWWWAAAAAAAAAAGGGGGGGGGGKKKKKKKKKKPPPPPPPPPPWWWWWWWWWWAAAAAAAAAAGGGGGGGGGGKKKKKKKKKKPPPPPPPPPPWWWWWWWWWWAAAAAAAAAAGGGGGGGGGGKKKKKKKKKKPPPPPPPPPPWWWWWWWWWW*\n' >> "$T/proteins.faa"
printf '>g3.t1\nMKAILVVLLYT*ATANADTLLILGDSLSAGYRMSASAAWPALLNDKWQSKTSVVNASISGDTSQQGLARLPALLKQHQPRWVLVELGGNDGLRGFQPQQTEQTLRQILQDVKAANAEPLLMQIRLPANYGRRYNEAFSAIYPKLAKEFDVPLLPFFMEEVYLKPQWMQDD*\n' >> "$T/proteins.faa"
printf '>g4.t1\nMKAILVVLLYTFATANADTLLILGDSLSAG*\n' >> "$T/proteins.faa"
printf '>g5.t1\nMKAILVVLLYTFATANADTLLILGDSLSAGYRMSASAAWPALLNDKWQSKTSVVNASISGDTSQQGLARLPALLKQHQPRWVLVELGGNDGLRGFQPQQTEQTLRQILQDVKAANAEPLLMQIRLPANYGRRYNEAFSAIYPKLAKEFD*\n' >> "$T/proteins.faa"
printf '>g5.t2\nMAAAAAAAAAAGGGGGGGGGGKKKKKKKKKKPPPPPPPPPPWWWWWWWWWWAAAAAAAAAAGGGGGGGGGGKKKKKKKKKKPPPPPPPPPPWWWWWWWWWWAAAAAAAAAAGGGGGGGGGGKKKKKKKKKKPPPPPPPPPPWWWWWWWWWWAAAAAAAAAAGGGGGGGGGG*\n' >> "$T/proteins.faa"

# qseqid sseqid pident length evalue bitscore qcovhsp scovhsp
printf 'g1.t1\tD_rerio_XP_001\t78.2\t400\t1e-120\t410.0\t92.0\t88.0\n' > "$T/hom.tsv"
printf 'g3.t1\tD_rerio_XP_009\t65.0\t150\t1e-40\t180.0\t80.0\t70.0\n' >> "$T/hom.tsv"
printf 'g4.t1\tD_rerio_XP_010\t70.0\t30\t1e-08\t60.0\t90.0\t20.0\n'  >> "$T/hom.tsv"
printf 'g5.t1\tS_canaliculatus_XP_1\t95.0\t149\t1e-90\t300.0\t99.0\t97.0\n' >> "$T/hom.tsv"

printf 'g1.t1\tSRR36642802\t14.220\n' > "$T/expr.tsv"
printf 'g5.t1\tSRR36642802\t8.400\n'  >> "$T/expr.tsv"
printf 'g2.t1\tSRR36642802\t0.000\n'  >> "$T/expr.tsv"

# transcript_id, cds_bp, overlap_bp, fraction
printf 'g1.t1\t1208\t0\t0.0000\n'      > "$T/rep.tsv"
printf 'g2.t1\t1001\t950\t0.9491\n'   >> "$T/rep.tsv"
printf 'g3.t1\t601\t0\t0.0000\n'      >> "$T/rep.tsv"
printf 'g4.t1\t101\t0\t0.0000\n'      >> "$T/rep.tsv"
printf 'g5.t1\t1802\t12\t0.0067\n'    >> "$T/rep.tsv"
printf 'g5.t2\t801\t40\t0.0499\n'     >> "$T/rep.tsv"
printf 'g6.t1\t1902\t0\t0.0000\n'     >> "$T/rep.tsv"

# Liftoff overlap: only g6 is covered by a lifted congener CDS.
printf 'g1.t1\t1208\t0\t0.0000\n'      > "$T/lift.tsv"
printf 'g2.t1\t1001\t0\t0.0000\n'     >> "$T/lift.tsv"
printf 'g5.t1\t1802\t1750\t0.9711\n'  >> "$T/lift.tsv"
printf 'g5.t2\t801\t0\t0.0000\n'      >> "$T/lift.tsv"
printf 'g6.t1\t1902\t1810\t0.9516\n'  >> "$T/lift.tsv"

python3 "$LIB/filter_models.py" \
  --gff3 "$T/in.gff3" --proteins "$T/proteins.faa" \
  --homology "$T/hom.tsv" --expression "$T/expr.tsv" \
  --repeat-overlap "$T/rep.tsv" --liftoff-overlap "$T/lift.tsv" \
  --keep-gff3 "$T/keep.gff3" --drop-gff3 "$T/drop.gff3" \
  --report "$T/report.tsv" --summary "$T/summary.txt"

echo
echo "--- decisions ---"
cut -f2,12,13,15,16 "$T/report.tsv" | column -t

fail=0
check() {  # check <desc> <expected> <actual>
  if [ "$2" = "$3" ]; then echo "  PASS  $1"; else echo "  FAIL  $1: want '$2' got '$3'"; fail=1; fi
}
dec() { awk -F'\t' -v t="$1" '$2==t{print $15"/"$16}' "$T/report.tsv"; }

echo
echo "--- assertions ---"
check "g1 kept (homology+expression)"     "keep/"                            "$(dec g1.t1)"
check "g2 dropped as repeat-derived"      "drop/repeat_derived_unsupported"  "$(dec g2.t1)"
check "g3 dropped for internal stop"      "drop/internal_stop"               "$(dec g3.t1)"
check "g4 dropped for short CDS"          "drop/cds_too_short"               "$(dec g4.t1)"
check "g5.t1 kept"                        "keep/"                            "$(dec g5.t1)"
check "g5.t2 dropped, no evidence"        "drop/no_evidence"                 "$(dec g5.t2)"
check "g6 RESCUED by congener axis alone" "keep/"                            "$(dec g6.t1)"
check "kept genes = 3"                    "3"  "$(grep -cP '\tgene\t' "$T/keep.gff3")"
check "dropped genes = 3"                 "3"  "$(grep -cP '\tgene\t' "$T/drop.gff3")"
check "g6 counted as congener-supported"  "1"  \
  "$(awk -F'\t' 'NR>1 && $2=="g6.t1" && $13+0 > 0' "$T/report.tsv" | wc -l)"
# Under the v1 rule set (no congener axis) g6 must fall out again.
python3 "$LIB/filter_models.py" \
  --gff3 "$T/in.gff3" --proteins "$T/proteins.faa" \
  --homology "$T/hom.tsv" --expression "$T/expr.tsv" --repeat-overlap "$T/rep.tsv" \
  --keep-gff3 "$T/keep_v1.gff3" --drop-gff3 "$T/drop_v1.gff3" \
  --report "$T/report_v1.tsv" --summary "$T/summary_v1.txt" >/dev/null
check "without congener axis g6 drops"    "drop/no_evidence" \
  "$(awk -F'\t' -v t=g6.t1 '$2==t{print $15"/"$16}' "$T/report_v1.tsv")"
check "g5 keeps only 1 transcript"        "1"  "$(awk -F'\t' '$3=="mRNA" && $9 ~ /Parent=g5/' "$T/keep.gff3" | wc -l)"
check "kept GFF3 has no g5.t2"            "0"  "$(grep -c 'ID=g5.t2' "$T/keep.gff3" || true)"

echo
[ "$fail" -eq 0 ] && echo "ALL TESTS PASSED" || { echo "TESTS FAILED"; exit 1; }
