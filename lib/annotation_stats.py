#!/usr/bin/env python3
"""Annotation statistics plus explicit pass/fail gating.

The acceptance thresholds are declared here, in code, rather than being judged
after the fact. They are set from published teleost annotations and from this
assembly's own genome-mode BUSCO (98.5%), so a run that quietly degrades is
caught instead of being written up as a result.
"""
import argparse
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gff3  # noqa: E402

# Reference points for the comparison table (published annotations).
PEERS = [
    ("Siganus fuscescens (2024, this genome)", 28351, None,  None),
    ("Siganus canaliculatus GCA_053572195.1",  24000, None,  None),
    ("Dicentrarchus labrax",                   26719, 9.6,   None),
    ("Larimichthys crocea",                    25401, 9.8,   None),
    ("Oryzias latipes",                        24141, 9.9,   None),
    ("Danio rerio",                            25592, 10.5,  None),
]

THRESHOLDS = {
    "busco_complete_min": 92.0,
    "gene_count_min": 20000,
    "gene_count_max": 27000,
    "mean_exons_min": 6.0,
    "internal_stops_max": 0,
    "mono_exonic_frac_max": 0.35,
}


def parse_busco(path):
    if not path or not os.path.exists(path):
        return {}
    txt = open(path).read()
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%",
                  txt)
    if not m:
        return {}
    return dict(complete=float(m.group(1)), single=float(m.group(2)),
                dup=float(m.group(3)), frag=float(m.group(4)),
                missing=float(m.group(5)))


def read_proteins(path):
    seqs, name, buf = {}, None, []
    if not path or not os.path.exists(path):
        return seqs
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


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff3", required=True)
    p.add_argument("--proteins")
    p.add_argument("--busco")
    p.add_argument("--evidence")
    p.add_argument("--out", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--check-thresholds", action="store_true")
    a = p.parse_args()

    genes, _ = gff3.read(a.gff3)
    genes = [g for g in genes if g.type == "gene"]

    n_tx = 0
    exon_counts, cds_lens, gene_lens, intron_lens, utr5, utr3 = [], [], [], [], 0, 0
    mono = 0
    for g in genes:
        gene_lens.append(g.end - g.start + 1)
        for t in gff3.transcripts(g):
            n_tx += 1
            ne = gff3.n_exons(t)
            exon_counts.append(ne)
            if ne == 1:
                mono += 1
            cds_lens.append(gff3.cds_len(t))
            ex = sorted((c.start, c.end) for c in t.children if c.type == "exon")
            for i in range(len(ex) - 1):
                intron_lens.append(ex[i + 1][0] - ex[i][1] - 1)
            for c in t.children:
                if c.type == "five_prime_UTR":
                    utr5 += c.end - c.start + 1
                elif c.type == "three_prime_UTR":
                    utr3 += c.end - c.start + 1

    prot = read_proteins(a.proteins)
    istop = sum(1 for s in prot.values() if "*" in s.rstrip("*"))
    busco = parse_busco(a.busco)

    n_with_utr = sum(1 for g in genes for t in gff3.transcripts(g)
                     if any(c.type in ("five_prime_UTR", "three_prime_UTR")
                            for c in t.children))

    # evidence_support.tsv column order is set by lib/filter_models.py; the
    # header is read rather than assumed so the two cannot drift apart.
    ev_tpm = ev_hom = ev_lift = ev_total = 0
    if a.evidence and os.path.exists(a.evidence):
        hdr = None
        for line in open(a.evidence):
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = {name: i for i, name in enumerate(f)}
                continue
            col = lambda k: f[hdr[k]] if k in hdr and hdr[k] < len(f) else ""
            if col("decision") != "keep":
                continue
            ev_total += 1
            num = lambda k: float(col(k) or 0)
            try:
                if num("best_bitscore") > 0:
                    ev_hom += 1
                if num("max_tpm") > 0:
                    ev_tpm += 1
                if num("liftoff_frac") > 0:
                    ev_lift += 1
            except ValueError:
                pass

    med = statistics.median
    stats = [
        ("genes", len(genes)),
        ("transcripts", n_tx),
        ("transcripts_per_gene", round(n_tx / len(genes), 2) if genes else 0),
        ("mean_exons_per_transcript", round(statistics.mean(exon_counts), 2) if exon_counts else 0),
        ("mono_exonic_transcripts", mono),
        ("mono_exonic_fraction", round(mono / n_tx, 3) if n_tx else 0),
        ("median_cds_len_bp", int(med(cds_lens)) if cds_lens else 0),
        ("mean_cds_len_bp", int(statistics.mean(cds_lens)) if cds_lens else 0),
        ("median_gene_len_bp", int(med(gene_lens)) if gene_lens else 0),
        ("median_intron_len_bp", int(med(intron_lens)) if intron_lens else 0),
        ("introns", len(intron_lens)),
        ("transcripts_with_utr", n_with_utr),
        ("utr5_total_bp", utr5),
        ("utr3_total_bp", utr3),
        ("proteins", len(prot)),
        ("internal_stop_proteins", istop),
    ]
    for k, v in busco.items():
        stats.append(("busco_" + k, v))
    if ev_total:
        stats += [
            ("kept_models_scored", ev_total),
            ("homology_supported", ev_hom),
            ("homology_supported_pct", round(100.0 * ev_hom / ev_total, 1)),
            ("transcript_backed", ev_tpm),
            ("transcript_backed_pct", round(100.0 * ev_tpm / ev_total, 1)),
            ("congener_supported", ev_lift),
            ("congener_supported_pct", round(100.0 * ev_lift / ev_total, 1)),
        ]

    with open(a.out, "w") as o:
        o.write("metric\tvalue\n")
        for k, v in stats:
            o.write("%s\t%s\n" % (k, v))

    L = []
    L.append("=== annotation statistics: %s ===" % os.path.basename(a.gff3))
    for k, v in stats:
        L.append("  %-32s %s" % (k, v))

    L.append("")
    L.append("=== comparison with published teleost annotations ===")
    L.append("  %-42s %8s" % ("annotation", "genes"))
    L.append("  %-42s %8d  <== this release" % ("Siganus fuscescens (this release)", len(genes)))
    for name, ng, _me, _x in PEERS:
        L.append("  %-42s %8d" % (name, ng))

    failures = []
    if a.check_thresholds:
        L.append("")
        L.append("=== acceptance gates ===")

        def gate(name, ok, detail):
            L.append("  %-4s %s (%s)" % ("PASS" if ok else "FAIL", name, detail))
            if not ok:
                failures.append(name)

        if busco:
            gate("BUSCO completeness", busco["complete"] >= THRESHOLDS["busco_complete_min"],
                 "C=%.1f%% need >=%.1f%%" % (busco["complete"], THRESHOLDS["busco_complete_min"]))
        else:
            L.append("  SKIP BUSCO completeness (no summary supplied)")
        gate("gene count in teleost range",
             THRESHOLDS["gene_count_min"] <= len(genes) <= THRESHOLDS["gene_count_max"],
             "%d, expect %d-%d" % (len(genes), THRESHOLDS["gene_count_min"],
                                   THRESHOLDS["gene_count_max"]))
        me = statistics.mean(exon_counts) if exon_counts else 0
        gate("mean exons per transcript", me >= THRESHOLDS["mean_exons_min"],
             "%.2f need >=%.1f" % (me, THRESHOLDS["mean_exons_min"]))
        gate("no internal stop codons", istop <= THRESHOLDS["internal_stops_max"],
             "%d" % istop)
        mf = mono / n_tx if n_tx else 0
        gate("mono-exonic fraction", mf <= THRESHOLDS["mono_exonic_frac_max"],
             "%.3f need <=%.2f" % (mf, THRESHOLDS["mono_exonic_frac_max"]))

        L.append("")
        L.append("RESULT: %s" % ("ALL GATES PASSED" if not failures
                                 else "FAILED: " + ", ".join(failures)))

    text = "\n".join(L)
    open(a.report, "w").write(text + "\n")
    print(text)
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
