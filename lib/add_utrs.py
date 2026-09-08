#!/usr/bin/env python3
"""Transfer UTRs from assembled transcripts onto CDS-only gene models.

BRAKER models are CDS-bounded, so the published annotation has no UTRs at all.
StringTie assemblies do have them, but naively grafting transcript ends onto a
prediction is how annotations acquire chimeric genes.

This is deliberately conservative. A StringTie transcript is accepted only if:
  1. same sequence and strand;
  2. every CDS base of the model lies inside the transcript's exons;
  3. the model's internal intron chain is a subset of the transcript's, so the
     two agree on splicing wherever they overlap;
  4. the extension does not run into a neighbouring gene's CDS.
When several candidates qualify, the one adding the most UTR wins. Models with
no qualifying transcript pass through untouched rather than being guessed at.
"""
import argparse
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gff3  # noqa: E402

TID = re.compile(r'transcript_id "([^"]+)"')


def read_stringtie(path):
    """GTF -> {tid: dict(seqid, strand, exons=[(s,e)...])}"""
    tx = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = TID.search(f[8])
        if not m:
            continue
        t = tx.setdefault(m.group(1),
                          {"seqid": f[0], "strand": f[6], "exons": []})
        t["exons"].append((int(f[3]), int(f[4])))
    for t in tx.values():
        t["exons"].sort()
        t["introns"] = introns_of(t["exons"])
        t["start"] = t["exons"][0][0]
        t["end"] = t["exons"][-1][1]
    return tx


def introns_of(exons):
    return frozenset((exons[i][1] + 1, exons[i + 1][0] - 1)
                     for i in range(len(exons) - 1))


def covered(blocks, s, e):
    """True if [s,e] is fully inside the union of blocks."""
    pos = s
    for bs, be in blocks:
        if be < pos:
            continue
        if bs > pos:
            return False
        pos = max(pos, be + 1)
        if pos > e:
            return True
    return pos > e


def subtract(exons, cds_blocks):
    """exon union minus CDS union -> list of (s,e) UTR blocks."""
    cut = sorted(cds_blocks)
    out = []
    for es, ee in exons:
        cur = [(es, ee)]
        for cs, ce in cut:
            nxt = []
            for s, e in cur:
                if ce < s or cs > e:
                    nxt.append((s, e))
                    continue
                if s < cs:
                    nxt.append((s, min(e, cs - 1)))
                if e > ce:
                    nxt.append((max(s, ce + 1), e))
            cur = nxt
        out.extend(cur)
    return [b for b in out if b[0] <= b[1]]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff3", required=True)
    p.add_argument("--transcripts", required=True, help="StringTie merged GTF")
    p.add_argument("--out", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--max-extension", type=int, default=10000,
                   help="refuse a UTR longer than this on one side")
    args = p.parse_args()

    genes, headers = gff3.read(args.gff3)
    stx = read_stringtie(args.transcripts)

    by_key = defaultdict(list)
    for tid, t in stx.items():
        by_key[(t["seqid"], t["strand"])].append((tid, t))

    # Every CDS block in the annotation, so an extension can be checked against
    # neighbouring genes rather than only against its own.
    cds_index = defaultdict(list)
    for g in genes:
        for t in gff3.transcripts(g):
            for c in t.children:
                if c.type == "CDS":
                    cds_index[(c.seqid, c.strand)].append((c.start, c.end, t.id))

    rows = []
    n_with_utr = 0
    for g in genes:
        gmax = None
        for t in gff3.transcripts(g):
            cds = sorted((c.start, c.end) for c in t.children if c.type == "CDS")
            if not cds:
                rows.append((t.id, "", 0, 0, "no_cds"))
                continue
            cds_s, cds_e = cds[0][0], cds[-1][1]
            model_introns = introns_of([(c.start, c.end) for c in t.children
                                        if c.type == "exon"] or cds)

            best = None
            for tid, st in by_key.get((t.seqid, t.strand), ()):
                if st["start"] > cds_s or st["end"] < cds_e:
                    continue                       # does not span the CDS
                if not all(covered(st["exons"], s, e) for s, e in cds):
                    continue                       # CDS not inside its exons
                if not model_introns <= st["introns"]:
                    continue                       # splicing disagrees
                gain = (cds_s - st["start"]) + (st["end"] - cds_e)
                if gain <= 0 or max(cds_s - st["start"],
                                    st["end"] - cds_e) > args.max_extension:
                    continue
                # Refuse an extension that would reach over another gene's CDS -
                # this is the failure mode that produces chimeric models.
                clash = False
                for o_s, o_e, o_tid in cds_index[(t.seqid, t.strand)]:
                    if o_tid == t.id:
                        continue
                    if o_e < st["start"] or o_s > st["end"]:
                        continue                   # outside the proposed span
                    if o_s < cds_s or o_e > cds_e:
                        clash = True               # foreign CDS in the extension
                        break
                if clash:
                    continue
                if best is None or gain > best[0]:
                    best = (gain, tid, st)

            if best is None:
                rows.append((t.id, "", 0, 0, "no_qualifying_transcript"))
                continue

            _gain, tid, st = best
            new_exons = st["exons"]
            utr_blocks = subtract(new_exons, cds)
            five = [b for b in utr_blocks
                    if (b[1] < cds_s) == (t.strand == "+")] if t.strand in "+-" else []
            three = [b for b in utr_blocks if b not in five]

            t.children = [c for c in t.children
                          if c.type not in ("exon", "five_prime_UTR",
                                            "three_prime_UTR")]
            for s, e in new_exons:
                f = gff3.Feature([t.seqid, t.source, "exon", str(s), str(e),
                                  ".", t.strand, ".", "Parent=%s" % t.id])
                t.children.append(f)
            for s, e in five:
                t.children.append(gff3.Feature(
                    [t.seqid, t.source, "five_prime_UTR", str(s), str(e),
                     ".", t.strand, ".", "Parent=%s" % t.id]))
            for s, e in three:
                t.children.append(gff3.Feature(
                    [t.seqid, t.source, "three_prime_UTR", str(s), str(e),
                     ".", t.strand, ".", "Parent=%s" % t.id]))

            t.start, t.end = st["start"], st["end"]
            gmax = (min(gmax[0], t.start), max(gmax[1], t.end)) if gmax else (t.start, t.end)
            n_with_utr += 1
            rows.append((t.id, tid,
                         sum(e - s + 1 for s, e in five),
                         sum(e - s + 1 for s, e in three), "ok"))

        if gmax:
            g.start, g.end = min(g.start, gmax[0]), max(g.end, gmax[1])

    gff3.write(args.out, genes, headers)
    with open(args.report, "w") as r:
        r.write("transcript_id\tsource_transcript\tutr5_bp\tutr3_bp\tstatus\n")
        for row in rows:
            r.write("\t".join(str(x) for x in row) + "\n")

    print("transcripts given UTRs: %d / %d" % (n_with_utr, len(rows)))
    print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
