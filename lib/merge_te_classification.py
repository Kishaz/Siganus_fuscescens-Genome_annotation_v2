#!/usr/bin/env python3
"""Fold TEsorter superfamily calls back into a RepeatModeler family library.

RepeatModeler leaves a large fraction of families as '#Unknown'. Those entries
make a repeat track hard to reuse and, worse, make it impossible to judge which
gene models are TE-derived. TEsorter assigns classes from protein domains; this
rewrites only the still-Unknown headers and reports what changed.
"""
import argparse, csv, os, re, sys


def read_fa(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].rstrip(), []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--library", required=True)
    p.add_argument("--tesorter", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--report", required=True)
    a = p.parse_args()

    calls = {}
    if os.path.exists(a.tesorter):
        with open(a.tesorter) as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd, None)
            for row in rd:
                if len(row) < 4:
                    continue
                seqid, order, superfam = row[0], row[1], row[2]
                if order and order.lower() not in ("unknown", "na", ""):
                    sf = superfam if superfam and superfam.lower() != "unknown" else None
                    calls[seqid.split("#")[0]] = f"{order}/{sf}" if sf else order
    else:
        print(f"WARNING: {a.tesorter} not found; leaving library unchanged",
              file=sys.stderr)

    n = reclassified = still_unknown = 0
    with open(a.out, "w") as o, open(a.report, "w") as rep:
        rep.write("family\told_class\tnew_class\tsource\n")
        for hdr, seq in read_fa(a.library):
            n += 1
            base = hdr.split()[0]
            fam, _, cls = base.partition("#")
            old = cls or "Unknown"
            new, src = old, "RepeatModeler"
            if old.split("/")[0] in ("Unknown", "Unspecified", ""):
                hit = calls.get(fam) or calls.get(base)
                if hit:
                    new, src = hit, "TEsorter"
                    reclassified += 1
                else:
                    still_unknown += 1
            rep.write(f"{fam}\t{old}\t{new}\t{src}\n")
            o.write(f">{fam}#{new}\n")
            for i in range(0, len(seq), 60):
                o.write(seq[i:i + 60] + "\n")

    print(f"families={n}  reclassified_by_TEsorter={reclassified}  "
          f"still_unknown={still_unknown}  -> {a.out}")


if __name__ == "__main__":
    main()
