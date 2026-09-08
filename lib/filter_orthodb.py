#!/usr/bin/env python3
"""Carve a clade-specific protein partition out of an OrthoDB partition file.

BRAKER's guidance is to use the smallest OrthoDB clade that contains the target
species. OrthoDB v12 ships no Actinopterygii partition, so the smallest
available option containing Siganus is Vertebrata - 19.4 M proteins, most of
them tetrapod and therefore poor hint donors for a fish, while dominating
BRAKER's runtime.

OrthoDB headers carry the NCBI taxid directly (`>7955_0:001234`), so the subset
can be cut exactly: resolve every descendant of the target clade from the NCBI
taxonomy dump, then keep only those records.

  python3 filter_orthodb.py --fasta Vertebrata.fa.gz --nodes nodes.dmp \\
      --taxon 7898 --out Actinopterygii.fa.gz --report partition_report.tsv
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict


def read_nodes(path):
    """nodes.dmp -> {parent_taxid: [child_taxid, ...]}"""
    children = defaultdict(list)
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            f = line.split("\t|\t")
            if len(f) < 2:
                continue
            try:
                tid, parent = int(f[0]), int(f[1])
            except ValueError:
                continue
            if tid != parent:
                children[parent].append(tid)
    return children


def descendants(children, root):
    """Iterative DFS - the vertebrate subtree is deep enough to blow recursion."""
    seen = {root}
    stack = [root]
    while stack:
        n = stack.pop()
        for c in children.get(n, ()):
            if c not in seen:
                seen.add(c)
                stack.append(c)
    return seen


def read_names(path, wanted):
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "scientific name" not in line:
                continue
            f = line.split("\t|\t")
            try:
                tid = int(f[0])
            except (ValueError, IndexError):
                continue
            if tid in wanted:
                out[tid] = f[1].strip()
    return out


def opener(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True)
    p.add_argument("--nodes", required=True, help="NCBI taxdump nodes.dmp")
    p.add_argument("--names", help="NCBI taxdump names.dmp (for the report)")
    p.add_argument("--taxon", type=int, default=7898,
                   help="clade taxid to keep (default 7898 = Actinopterygii)")
    p.add_argument("--out", required=True)
    p.add_argument("--report")
    p.add_argument("--min-len", type=int, default=50)
    a = p.parse_args()

    sys.stderr.write("reading taxonomy...\n")
    keep = descendants(read_nodes(a.nodes), a.taxon)
    sys.stderr.write("clade %d has %d descendant taxa\n" % (a.taxon, len(keep)))

    per_taxon = defaultdict(int)
    kept = total = short = 0
    name = None
    buf = []
    keeping = False
    cur_tax = None

    def flush(out):
        nonlocal kept, short
        if not keeping or name is None:
            return
        seq = "".join(buf)
        if len(seq) < a.min_len:
            short += 1
            return
        out.write(">%s\n" % name)
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + "\n")
        kept += 1
        per_taxon[cur_tax] += 1

    with opener(a.fasta) as fh, opener(a.out, "wt") as out:
        for line in fh:
            if line.startswith(">"):
                flush(out)
                total += 1
                name = line[1:].rstrip("\n")
                head = name.split()[0]
                try:
                    cur_tax = int(head.split("_", 1)[0])
                except ValueError:
                    cur_tax = None
                keeping = cur_tax in keep
                buf = []
            else:
                if keeping:
                    buf.append(line.strip())
        flush(out)

    sys.stderr.write("proteins: %d in -> %d kept (%.2f%%), %d dropped as short\n"
                     % (total, kept, 100.0 * kept / total if total else 0, short))
    sys.stderr.write("species represented: %d\n" % len(per_taxon))

    if a.report:
        names = read_names(a.names, set(per_taxon)) if a.names else {}
        with open(a.report, "w") as r:
            r.write("taxid\tspecies\tproteins\n")
            for t, n in sorted(per_taxon.items(), key=lambda x: -x[1]):
                r.write("%s\t%s\t%d\n" % (t, names.get(t, ""), n))
            r.write("TOTAL\t\t%d\n" % kept)
        sys.stderr.write("wrote %s\n" % a.report)

    if kept == 0:
        sys.exit("ERROR: nothing kept - is the taxid correct for this file?")


if __name__ == "__main__":
    main()
