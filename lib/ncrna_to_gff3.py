#!/usr/bin/env python3
"""Merge tRNAscan-SE, barrnap and Infernal/Rfam output into one ncRNA GFF3.

Rfam hits are deduplicated against the tRNA and rRNA calls, which are more
specific for their own classes; without that step the file double-counts every
tRNA once as a tRNA and again as an Rfam family.
"""
import argparse
import os
import re
from collections import defaultdict

RFAM_TYPE = [
    (re.compile(r"^tRNA", re.I), "tRNA"),
    (re.compile(r"(^|_)rRNA|^LSU|^SSU|^5S|^5_8S", re.I), "rRNA"),
    (re.compile(r"^snoRNA|SNOR|^U3$|^U8$|^U13$", re.I), "snoRNA"),
    (re.compile(r"^mir-|^MIR|microRNA", re.I), "miRNA"),
    (re.compile(r"^U[0-9]|snRNA|spliceosom", re.I), "snRNA"),
    (re.compile(r"^SRP|7SL", re.I), "SRP_RNA"),
    (re.compile(r"RNase_P|RNaseP", re.I), "RNase_P_RNA"),
    (re.compile(r"telomerase|TERC", re.I), "telomerase_RNA"),
    (re.compile(r"riboswitch|ribozyme", re.I), "ribozyme"),
    (re.compile(r"^Y_RNA|^Vault", re.I), "ncRNA"),
]


def rfam_so(name):
    for rx, so in RFAM_TYPE:
        if rx.search(name):
            return so
    return "ncRNA"


def esc(s):
    for a, b in [("%", "%25"), (";", "%3B"), ("=", "%3D"),
                 ("&", "%26"), (",", "%2C")]:
        s = s.replace(a, b)
    return s


def read_trnascan(path):
    """tRNAscan-SE tabular: seq, idx, begin, end, type, anticodon, ..., score"""
    out = []
    if not path or not os.path.exists(path):
        return out
    for line in open(path):
        f = line.split()
        if len(f) < 9 or not f[2].isdigit():
            continue
        s, e = int(f[2]), int(f[3])
        strand = "+" if s <= e else "-"
        out.append(dict(seqid=f[0], start=min(s, e), end=max(s, e), strand=strand,
                        so="tRNA", name="tRNA-%s" % f[4], anticodon=f[5],
                        score=f[8] if len(f) > 8 else "."))
    return out


def read_barrnap(path):
    out = []
    if not path or not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        m = re.search(r"Name=([^;]+)", f[8])
        out.append(dict(seqid=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                        so="rRNA", name=m.group(1) if m else "rRNA", score=f[5]))
    return out


def read_mirmachine(path):
    """MirMachine GFF of precursor miRNA calls."""
    out = []
    if not path or not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        m = re.search(r"(?:Name|ID)=([^;]+)", f[8])
        try:
            s, e = int(f[3]), int(f[4])
        except ValueError:
            continue
        out.append(dict(seqid=f[0], start=min(s, e), end=max(s, e),
                        strand=f[6] if f[6] in "+-" else "+",
                        so="miRNA", name=m.group(1) if m else "miRNA",
                        score=f[5]))
    return out


def read_rfam(path):
    """Infernal `cmscan --fmt 2 --tblout` reader.

    Column order (0-based) is fixed by Infernal:
        0 idx  1 target  2 accession  3 query  4 accession  5 clan  6 mdl
        7 mdl-from  8 mdl-to  9 seq-from  10 seq-to  11 strand  12 trunc
        13 pass  14 gc  15 bias  16 score  17 E-value  18 inc  19 olp
    The trailing description field may contain spaces, so only these fixed
    leading indices are trusted.
    """
    out = []
    if not path or not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 20:
            continue
        name, acc = f[1], f[2]
        try:
            s, e = int(f[9]), int(f[10])
        except ValueError:
            continue
        if f[18] != "!":
            continue                      # below the reporting threshold
        if f[19] not in ("^", "="):
            continue                      # a lower-scoring member of a clan
        strand = f[11]
        out.append(dict(seqid=f[3], start=min(s, e), end=max(s, e),
                        strand=strand if strand in "+-" else "+",
                        so=rfam_so(name), name=name, rfam=acc, score=f[16]))
    return out


def overlaps(a, b, frac=0.5):
    if a["seqid"] != b["seqid"]:
        return False
    lo, hi = max(a["start"], b["start"]), min(a["end"], b["end"])
    if hi < lo:
        return False
    ov = hi - lo + 1
    return ov >= frac * min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--trnascan")
    p.add_argument("--rrna")
    p.add_argument("--rfam")
    p.add_argument("--mirmachine")
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--source", default="ncRNA-pipeline")
    a = p.parse_args()

    trna = read_trnascan(a.trnascan)
    rrna = read_barrnap(a.rrna)
    mirs = read_mirmachine(a.mirmachine)
    rfam = read_rfam(a.rfam)

    # Each specialist tool wins for its own class: tRNAscan-SE for tRNA, barrnap
    # for rRNA, MirMachine for miRNA. Rfam is the generalist and fills the rest,
    # so its hits are dropped where a specialist already called the same locus.
    specific = defaultdict(list)
    for r in trna + rrna + mirs:
        specific[r["seqid"]].append(r)
    specialist_classes = {"tRNA", "rRNA", "miRNA"}
    kept_rfam = []
    for r in rfam:
        if r["so"] in specialist_classes and \
           any(overlaps(r, s) for s in specific.get(r["seqid"], ())):
            continue
        kept_rfam.append(r)

    feats = trna + rrna + mirs + kept_rfam
    feats.sort(key=lambda x: (x["seqid"], x["start"], x["end"]))

    counters = defaultdict(int)
    with open(a.out, "w") as o:
        o.write("##gff-version 3\n")
        for f in feats:
            counters[f["so"]] += 1
            fid = "%s_%05d" % (f["so"], counters[f["so"]])
            attrs = ["ID=%s" % fid, "Name=%s" % esc(f["name"]),
                     "gene_biotype=%s" % f["so"]]
            if f.get("anticodon"):
                attrs.append("anticodon=%s" % esc(f["anticodon"]))
            if f.get("rfam"):
                attrs.append("Dbxref=RFAM:%s" % f["rfam"])
            o.write("\t".join([f["seqid"], a.source, f["so"],
                               str(f["start"]), str(f["end"]),
                               str(f.get("score", ".")), f["strand"], ".",
                               ";".join(attrs)]) + "\n")

    with open(a.summary, "w") as s:
        s.write("ncRNA_class\tcount\ttotal_bp\n")
        agg = defaultdict(lambda: [0, 0])
        for f in feats:
            agg[f["so"]][0] += 1
            agg[f["so"]][1] += f["end"] - f["start"] + 1
        for k, (n, bp) in sorted(agg.items(), key=lambda x: -x[1][0]):
            s.write("%s\t%d\t%d\n" % (k, n, bp))
        s.write("TOTAL\t%d\t%d\n" % (len(feats), sum(v[1] for v in agg.values())))

    print("ncRNA features: %d (tRNAscan %d, barrnap %d, MirMachine %d, "
          "Rfam %d kept of %d)"
          % (len(feats), len(trna), len(rrna), len(mirs), len(kept_rfam), len(rfam)))


if __name__ == "__main__":
    main()
