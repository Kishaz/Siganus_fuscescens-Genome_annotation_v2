#!/usr/bin/env python3
"""Assign stable, sortable, versioned identifiers to a final annotation.

BRAKER emits names like g1234.t1 that carry no species information and change
completely between runs. A public release needs IDs that are stable, obviously
ours, and traceable back to the run that produced them - so the old identifier
is preserved on every feature and a full old->new map is written alongside.

  gene        SFUS_G000001
  transcript  SFUS_G000001.1
  exon/CDS    SFUS_G000001.1.exon1 / .cds1
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gff3  # noqa: E402

SUFFIX = {"exon": "exon", "CDS": "cds",
          "five_prime_UTR": "utr5", "three_prime_UTR": "utr3",
          "start_codon": "start", "stop_codon": "stop",
          "intron": "intron"}

# Canonical attribute order for a readable public file.
ATTR_ORDER = ["ID", "Parent", "Name", "gene_biotype", "product",
              "old_locus_tag", "Note", "Dbxref", "Ontology_term"]


def reorder(feat):
    """Put ID/Parent first, then known keys, then anything else alphabetically."""
    a = feat.attrs
    rest = sorted(k for k in a if k not in ATTR_ORDER)
    feat.attrs = type(a)((k, a[k]) for k in ATTR_ORDER + rest if k in a)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff3", required=True)
    p.add_argument("--prefix", default="SFUS")
    p.add_argument("--out", required=True)
    p.add_argument("--map", required=True)
    p.add_argument("--pad", type=int, default=6)
    p.add_argument("--step", type=int, default=1,
                   help="gene numbering increment; >1 leaves room for later insertions")
    p.add_argument("--source", default=None,
                   help="override the GFF3 source column, e.g. BRAKER3")
    args = p.parse_args()

    genes, headers = gff3.read(args.gff3)
    # Sort by position so IDs increase along the assembly - reviewers and
    # browsers both expect this, and it makes diffs between versions readable.
    genes.sort(key=lambda g: (g.seqid, g.start, g.end))

    rows = []
    n = 0
    for g in genes:
        n += args.step
        gid = "%s_G%0*d" % (args.prefix, args.pad, n)
        old_gid = g.id or ""
        rows.append(("gene", old_gid, gid))
        g.attrs["ID"] = gid
        if old_gid:
            g.attrs["old_locus_tag"] = old_gid
        g.attrs.pop("Name", None)
        if args.source:
            g.source = args.source

        for ti, t in enumerate(sorted(gff3.transcripts(g),
                                      key=lambda x: (x.start, x.end)), 1):
            tid = "%s.%d" % (gid, ti)
            old_tid = t.id or ""
            rows.append((t.type, old_tid, tid))
            t.attrs["ID"] = tid
            t.attrs["Parent"] = gid
            if old_tid:
                t.attrs["old_locus_tag"] = old_tid
            if args.source:
                t.source = args.source

            counters = {}
            for c in sorted(t.children, key=lambda x: (x.type, x.start)):
                sfx = SUFFIX.get(c.type, c.type.lower())
                counters[sfx] = counters.get(sfx, 0) + 1
                c.attrs["ID"] = "%s.%s%d" % (tid, sfx, counters[sfx])
                c.attrs["Parent"] = tid
                c.attrs.pop("old_locus_tag", None)
                if args.source:
                    c.source = args.source
                reorder(c)
            reorder(t)
        reorder(g)

    gff3.write(args.out, genes, headers)
    with open(args.map, "w") as m:
        m.write("feature_type\told_id\tnew_id\n")
        for r in rows:
            m.write("\t".join(r) + "\n")

    n_tx = sum(1 for r in rows if r[0] in ("mRNA", "transcript"))
    print("assigned %d genes / %d transcripts -> %s" % (len(genes), n_tx, args.out))
    print("id map -> %s" % args.map)


if __name__ == "__main__":
    main()
