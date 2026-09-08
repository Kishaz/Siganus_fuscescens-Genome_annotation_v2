#!/usr/bin/env python3
"""Score BRAKER gene models on independent evidence axes, then filter.

Philosophy: never drop a model merely for being "extra". A model is removed only
when it is structurally invalid, or when it fails EVERY independent line of
evidence at once. Everything kept carries its evidence in the report, so a user
can apply a stricter threshold themselves without re-running the pipeline.

Evidence axes (deliberately independent of one another):
  homology    DIAMOND blastp against the curated relative panel
  expression  StringTie -e TPM from the conspecific RNA-seq library
  congener    overlap with the Liftoff transfer from S. canaliculatus
  repeats     fraction of CDS bases inside a repeat feature

The congener axis matters more than it looks. With a single RNA-seq library,
homology and expression fail together for one underlying reason at poorly
conserved, tissue-restricted loci - so treating them as two independent votes
overstates the evidence. A lift from a chromosome-level congener gives such a
locus a genuinely separate way to prove it is real.
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gff3  # noqa: E402


def read_proteins(path):
    seqs = {}
    if not path or not os.path.exists(path):
        return seqs
    name = None
    buf = []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def read_homology(path):
    """DIAMOND outfmt6: qseqid sseqid pident length evalue bitscore qcovhsp scovhsp"""
    best = {}
    if not path or not os.path.exists(path):
        return best
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        try:
            bits, qcov = float(f[5]), float(f[6])
        except ValueError:
            continue
        if f[0] not in best or bits > best[f[0]][0]:
            best[f[0]] = (bits, qcov, f[1])
    return best


def read_expression(path):
    expr = defaultdict(float)
    if not path or not os.path.exists(path):
        return expr
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            expr[f[0]] = max(expr[f[0]], float(f[2]))
        except ValueError:
            continue
    return expr


def read_repeat_overlap(path):
    rep = {}
    if not path or not os.path.exists(path):
        return rep
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        try:
            rep[f[0]] = float(f[3])
        except ValueError:
            continue
    return rep


def classify(tx, args, hom, expr, rep, prot, lift):
    """Return (drop_reason_or_None, evidence_dict) for one transcript."""
    tid = tx.id or ""
    cl = gff3.cds_len(tx)
    ne = gff3.n_exons(tx)
    bits, qcov, hit = hom.get(tid, (0.0, 0.0, ""))
    tpm = expr.get(tid, 0.0)
    rfrac = rep.get(tid, 0.0)
    lfrac = lift.get(tid, 0.0)
    seq = prot.get(tid, "")
    internal_stop = ("*" in seq[:-1]) if seq else False

    has_hom = bits >= args.min_bitscore and qcov >= args.min_qcov
    expressed = tpm >= args.min_tpm
    lifted = lfrac >= args.min_liftoff_frac
    repeaty = rfrac >= args.max_repeat_frac
    supported = has_hom or expressed or lifted

    why = None
    if internal_stop:
        why = "internal_stop"
    elif cl < args.min_cds_len:
        why = "cds_too_short"
    elif repeaty and not supported:
        why = "repeat_derived_unsupported"
    elif not supported:
        why = "no_evidence"
    elif ne == 1 and cl < 300 and not supported:
        why = "short_single_exon_unsupported"

    ev = dict(tid=tid, exons=ne, cds_len=cl, bits=bits, qcov=qcov,
              hit=hit, tpm=tpm, rfrac=rfrac, lfrac=lfrac)
    return why, ev


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff3", required=True)
    p.add_argument("--proteins")
    p.add_argument("--homology")
    p.add_argument("--expression")
    p.add_argument("--repeat-overlap", dest="repeat_overlap")
    p.add_argument("--liftoff-overlap", dest="liftoff_overlap",
                   help="TSV: transcript_id, cds_bp, overlap_bp, fraction")
    p.add_argument("--keep-gff3", dest="keep_gff3", required=True)
    p.add_argument("--drop-gff3", dest="drop_gff3", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--summary", required=True)
    # Thresholds are surfaced as flags and echoed into the summary so the
    # release README can state exactly what produced the shipped file.
    p.add_argument("--min-bitscore", type=float, default=50.0)
    p.add_argument("--min-qcov", type=float, default=30.0)
    p.add_argument("--min-tpm", type=float, default=0.5)
    p.add_argument("--min-liftoff-frac", type=float, default=0.50)
    p.add_argument("--max-repeat-frac", type=float, default=0.60)
    p.add_argument("--min-cds-len", type=int, default=150)
    args = p.parse_args()

    genes, headers = gff3.read(args.gff3)
    hom = read_homology(args.homology)
    expr = read_expression(args.expression)
    rep = read_repeat_overlap(args.repeat_overlap)
    lift = read_repeat_overlap(args.liftoff_overlap)   # same 4-column shape
    prot = read_proteins(args.proteins)

    keep, drop = [], []
    reasons = defaultdict(int)
    rows = []

    for g in genes:
        txs = gff3.transcripts(g)
        if not txs:
            drop.append(g)
            reasons["no_transcript"] += 1
            continue

        dropped_tx = []
        kept_any = False
        for t in txs:
            why, ev = classify(t, args, hom, expr, rep, prot, lift)
            rows.append([
                g.id or "", ev["tid"], t.seqid, t.start, t.end, t.strand,
                ev["exons"], ev["cds_len"], "%.1f" % ev["bits"],
                "%.1f" % ev["qcov"], ev["hit"], "%.3f" % ev["tpm"],
                "%.3f" % ev["lfrac"], "%.3f" % ev["rfrac"],
                "keep" if why is None else "drop", why or "",
            ])
            if why is None:
                kept_any = True
            else:
                dropped_tx.append(t)
                reasons[why] += 1

        if kept_any:
            if dropped_tx:
                g.children = [c for c in g.children if c not in dropped_tx]
            keep.append(g)
        else:
            drop.append(g)

    gff3.write(args.keep_gff3, keep, headers)
    gff3.write(args.drop_gff3, drop, headers)

    with open(args.report, "w") as r:
        r.write("gene_id\ttranscript_id\tseqid\tstart\tend\tstrand\texons\tcds_len\t"
                "best_bitscore\tquery_cov\tbest_hit\tmax_tpm\tliftoff_frac\t"
                "repeat_frac\tdecision\tdrop_reason\n")
        for row in rows:
            r.write("\t".join(str(x) for x in row) + "\n")

    kept_rows = [x for x in rows if x[14] == "keep"]
    n_keep_tx = len(kept_rows)
    lines = []
    lines.append("=== model filtering ===")
    lines.append("thresholds: bitscore>=%g qcov>=%g%% tpm>=%g liftoff_frac>=%g "
                 "repeat_frac<%g cds>=%dbp"
                 % (args.min_bitscore, args.min_qcov, args.min_tpm,
                    args.min_liftoff_frac, args.max_repeat_frac, args.min_cds_len))
    lines.append("input genes            %d" % len(genes))
    lines.append("kept genes             %d" % len(keep))
    lines.append("dropped genes          %d" % len(drop))
    lines.append("input transcripts      %d" % len(rows))
    lines.append("kept transcripts       %d" % n_keep_tx)
    lines.append("--- drop reasons ---")
    for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
        lines.append("  %-32s %d" % (k, v))
    lines.append("--- evidence on kept models ---")
    if n_keep_tx:
        pc = lambda k: 100.0 * k / n_keep_tx
        with_hom = sum(1 for x in kept_rows if float(x[8]) > 0)
        with_tpm = sum(1 for x in kept_rows if float(x[11]) > 0)
        with_lift = sum(1 for x in kept_rows if float(x[12]) > 0)
        all3 = sum(1 for x in kept_rows
                   if float(x[8]) > 0 and float(x[11]) > 0 and float(x[12]) > 0)
        lines.append("  homology-supported            %d/%d (%.1f%%)"
                     % (with_hom, n_keep_tx, pc(with_hom)))
        lines.append("  transcript-backed (TPM>0)     %d/%d (%.1f%%)"
                     % (with_tpm, n_keep_tx, pc(with_tpm)))
        lines.append("  congener-supported (Liftoff)  %d/%d (%.1f%%)"
                     % (with_lift, n_keep_tx, pc(with_lift)))
        lines.append("  all three independent axes    %d/%d (%.1f%%)"
                     % (all3, n_keep_tx, pc(all3)))

    text = "\n".join(lines)
    open(args.summary, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
