#!/usr/bin/env python3
"""Rewrite the seqid column of a GFF/GTF through the contig name map.

Ships the same annotation under GenBank accessions, NCBI sequence names, and the
original assembly contig names, so a user never has to rename anything to match
whichever FASTA they already have. Coordinates are untouched.
"""
import argparse
import os
import sys


def load_map(path, frm, to):
    hdr = None
    m = {}
    for i, line in enumerate(open(path)):
        f = line.rstrip("\n").split("\t")
        if i == 0:
            hdr = f
            if frm not in hdr or to not in hdr:
                sys.exit("ERROR: map lacks column '%s' or '%s' (has: %s)"
                         % (frm, to, ",".join(hdr)))
            continue
        d = dict(zip(hdr, f))
        m[d[frm]] = d[to]
    return m


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True)
    p.add_argument("--map", required=True)
    p.add_argument("--from", dest="frm", required=True)
    p.add_argument("--to", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    m = load_map(a.map, a.frm, a.to)
    d = os.path.dirname(a.out)
    if d:
        os.makedirs(d, exist_ok=True)

    n = 0
    unmapped = set()
    with open(a.out, "w") as o:
        for line in open(a.gff):
            if line.startswith("#"):
                # ##sequence-region lines carry a seqid too
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) > 1 and parts[1] in m:
                        parts[1] = m[parts[1]]
                        o.write(" ".join(parts) + "\n")
                        continue
                o.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                o.write(line)
                continue
            if f[0] in m:
                f[0] = m[f[0]]
                n += 1
            else:
                unmapped.add(f[0])
            o.write("\t".join(f) + "\n")

    if unmapped:
        sys.exit("ERROR: %d sequence name(s) absent from the map: %s"
                 % (len(unmapped), ", ".join(sorted(unmapped)[:5])))
    print("renamed %d feature lines (%s -> %s) -> %s" % (n, a.frm, a.to, a.out))


if __name__ == "__main__":
    main()
