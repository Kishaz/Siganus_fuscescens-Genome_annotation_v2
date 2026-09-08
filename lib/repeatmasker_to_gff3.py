#!/usr/bin/env python3
"""Convert a RepeatMasker .out table into a valid, sorted GFF3 + class summary.

RepeatMasker's own -gff output is minimal and loses the class/family split and
divergence, which is exactly what makes a repeat track reusable. This keeps them.
"""
import argparse, re, sys
from collections import defaultdict

# Map RepeatMasker class/family strings onto Sequence Ontology terms so the GFF3
# is interpretable by genome browsers rather than being all "repeat_region".
SO = {
    "SINE": "SINE_element", "LINE": "LINE_element",
    "LTR": "LTR_retrotransposon", "DNA": "terminal_inverted_repeat_element",
    "RC": "repeat_region", "Retroposon": "retrotransposon",
    "Satellite": "satellite_DNA", "Simple_repeat": "microsatellite",
    "Low_complexity": "low_complexity_region", "rRNA": "rRNA_gene",
    "tRNA": "tRNA_gene", "snRNA": "snRNA_gene", "scRNA": "scRNA_gene",
    "srpRNA": "SRP_RNA_gene", "Unknown": "repeat_region",
    "Unspecified": "repeat_region", "ARTEFACT": "repeat_region",
}


def so_term(cls):
    return SO.get(cls.split("/")[0], "repeat_region")


def parse(path):
    """Yield parsed .out rows, skipping headers and concatenation artefacts."""
    num = re.compile(r"^\s*\d+\s")
    for line in open(path, errors="replace"):
        if not num.match(line):
            continue                      # header / blank / re-concatenated header
        f = line.split()
        if len(f) < 15:
            continue
        try:
            sw, div = int(f[0]), float(f[1])
            seqid, start, end = f[4], int(f[5]), int(f[6])
        except ValueError:
            continue
        strand = "-" if f[8] in ("C", "-") else "+"
        name, clsfam = f[9], f[10]
        # For '-' strand RepeatMasker prints (left) first; the two real
        # coordinates are always the middle pair of the trailing triplet.
        trip = [f[11], f[12], f[13]]
        vals = [int(x.strip("()")) for x in trip]
        rstart, rend = (vals[1], vals[2]) if strand == "-" else (vals[0], vals[1])
        yield dict(seqid=seqid, start=start, end=end, sw=sw, div=div,
                   strand=strand, name=name, clsfam=clsfam,
                   rstart=min(rstart, rend), rend=max(rstart, rend))


def esc(s):
    for a, b in [("%", "%25"), (";", "%3B"), ("=", "%3D"),
                 ("&", "%26"), (",", "%2C")]:
        s = s.replace(a, b)
    return s


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True, help="RepeatMasker .out")
    p.add_argument("--gff3", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--source", default="RepeatMasker")
    p.add_argument("--min-len", type=int, default=1)
    a = p.parse_args()

    rows = [r for r in parse(a.out) if r["end"] - r["start"] + 1 >= a.min_len]
    if not rows:
        sys.exit(f"ERROR: no usable rows parsed from {a.out}")
    rows.sort(key=lambda r: (r["seqid"], r["start"], r["end"]))

    stats = defaultdict(lambda: [0, 0])          # class -> [count, bp]
    seen = defaultdict(int)
    with open(a.gff3, "w") as g:
        g.write("##gff-version 3\n")
        for r in rows:
            cls = r["clsfam"]
            k = cls.split("/")[0]
            stats[cls][0] += 1
            stats[cls][1] += r["end"] - r["start"] + 1
            seen[k] += 1
            rid = f"repeat_{k}_{seen[k]:07d}"
            attrs = (f"ID={rid};Name={esc(r['name'])};"
                     f"class={esc(k)};family={esc(cls)};"
                     f"Target={esc(r['name'])} {r['rstart']} {r['rend']};"
                     f"perc_divergence={r['div']:.1f}")
            g.write("\t".join([r["seqid"], a.source, so_term(cls),
                               str(r["start"]), str(r["end"]), str(r["sw"]),
                               r["strand"], ".", attrs]) + "\n")

    total_bp = sum(v[1] for v in stats.values())
    with open(a.summary, "w") as s:
        s.write("class_family\telements\tbp\tpct_of_repeat_bp\n")
        for k, (n, bp) in sorted(stats.items(), key=lambda x: -x[1][1]):
            s.write(f"{k}\t{n}\t{bp}\t{bp / total_bp * 100:.3f}\n")
        s.write(f"TOTAL\t{len(rows)}\t{total_bp}\t100.000\n")

    print(f"wrote {a.gff3} ({len(rows)} features) and {a.summary} "
          f"({total_bp:,} repeat bp before overlap merging)")


if __name__ == "__main__":
    main()
