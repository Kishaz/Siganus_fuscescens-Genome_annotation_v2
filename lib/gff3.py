#!/usr/bin/env python3
"""Minimal, dependency-free GFF3 reader/writer preserving gene->mRNA->child structure.

Deliberately not using a library: this has to run inside whatever environment a
downstream user happens to have, and the structures we touch are simple. It
preserves attribute order so round-tripping does not gratuitously reshuffle files.
"""
from collections import OrderedDict


class Feature:
    __slots__ = ("seqid", "source", "type", "start", "end", "score",
                 "strand", "phase", "attrs", "children")

    def __init__(self, cols):
        self.seqid, self.source, self.type = cols[0], cols[1], cols[2]
        self.start, self.end = int(cols[3]), int(cols[4])
        self.score, self.strand, self.phase = cols[5], cols[6], cols[7]
        self.attrs = OrderedDict()
        for kv in cols[8].rstrip(";").split(";"):
            if not kv:
                continue
            k, _, v = kv.partition("=")
            self.attrs[k.strip()] = v
        self.children = []

    @property
    def id(self):
        return self.attrs.get("ID")

    @property
    def parent(self):
        p = self.attrs.get("Parent")
        return p.split(",") if p else []

    def to_row(self):
        attr = ";".join("%s=%s" % (k, v)
                        for k, v in self.attrs.items() if v is not None)
        return "\t".join([self.seqid, self.source, self.type,
                          str(self.start), str(self.end), self.score,
                          self.strand, self.phase, attr])


def read(path):
    """Return (roots, headers). Features whose Parent is absent become roots."""
    by_id = {}
    roots = []
    headers = []
    pending = []
    for line in open(path):
        if line.startswith("#"):
            if line.startswith("##gff-version") or line.startswith("##sequence-region"):
                headers.append(line.rstrip("\n"))
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 9:
            continue
        f = Feature(cols)
        if f.id:
            by_id[f.id] = f
        if f.parent:
            pending.append(f)
        else:
            roots.append(f)
    for f in pending:
        attached = False
        for pid in f.parent:
            if pid in by_id:
                by_id[pid].children.append(f)
                attached = True
                break
        if not attached:
            roots.append(f)
    return roots, headers


CHILD_ORDER = {"mRNA": 0, "transcript": 0, "exon": 1, "five_prime_UTR": 2,
               "CDS": 3, "three_prime_UTR": 4}


def _emit(out, f):
    out.write(f.to_row() + "\n")
    for c in sorted(f.children, key=lambda x: (CHILD_ORDER.get(x.type, 9), x.start)):
        _emit(out, c)


def write(path, roots, headers=None):
    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for h in (headers or []):
            if h.startswith("##sequence-region"):
                out.write(h + "\n")
        for r in sorted(roots, key=lambda x: (x.seqid, x.start, x.end)):
            _emit(out, r)


def transcripts(gene):
    return [c for c in gene.children
            if c.type in ("mRNA", "transcript", "ncRNA", "pseudogenic_transcript")]


def cds_len(tx):
    return sum(c.end - c.start + 1 for c in tx.children if c.type == "CDS")


def n_exons(tx):
    n = sum(1 for c in tx.children if c.type == "exon")
    return n or sum(1 for c in tx.children if c.type == "CDS")
