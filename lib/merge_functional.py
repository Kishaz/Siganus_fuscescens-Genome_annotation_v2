#!/usr/bin/env python3
"""Fold InterProScan, eggNOG-mapper and Swiss-Prot results into the annotation.

Emits two things from one pass so they can never disagree:
  - a GFF3 with product/Dbxref/Ontology_term attributes on every mRNA
  - a flat TSV for people who would rather open a spreadsheet than parse GFF3

Product-name precedence is deliberate: Swiss-Prot description first (curated and
human-readable), then eggNOG's free-text description, then an InterPro domain
name, and finally "hypothetical protein". Nothing is invented - a model with no
evidence is labelled as such rather than being given a plausible-sounding name.
"""
import argparse
import os
import re
from collections import defaultdict

SPROT_DESC = re.compile(r"^\S+\s+(.*?)\s+OS=")
SPROT_GN = re.compile(r"\bGN=(\S+)")


def esc(s):
    for a, b in [("%", "%25"), (";", "%3B"), ("=", "%3D"),
                 ("&", "%26"), (",", "%2C")]:
        s = s.replace(a, b)
    return s


def read_interproscan(path):
    """IPR TSV: prot, md5, len, analysis, sig_acc, sig_desc, start, end, score,
    status, date, [ipr_acc, ipr_desc, [GO, [pathways]]]"""
    dom = defaultdict(list)
    ipr = defaultdict(set)
    go = defaultdict(set)
    pw = defaultdict(set)
    if not path or not os.path.exists(path):
        return dom, ipr, go, pw
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        q = f[0]
        if f[3] and f[4]:
            dom[q].append("%s:%s" % (f[3], f[4]))
        if len(f) > 11 and f[11].startswith("IPR"):
            ipr[q].add(f[11])
            if len(f) > 12 and f[12]:
                dom[q].append(f[12])
        if len(f) > 13 and f[13]:
            for g in re.findall(r"GO:\d{7}", f[13]):
                go[q].add(g)
        if len(f) > 14 and f[14]:
            for t in f[14].split("|"):
                if t.strip():
                    pw[q].add(t.strip())
    return dom, ipr, go, pw


def read_eggnog(path):
    """emapper .annotations; column names come from the '#query' header line."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    hdr = None
    for line in open(path):
        if line.startswith("#"):
            if line.startswith("#query"):
                hdr = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if not hdr:
            continue
        f = line.rstrip("\n").split("\t")
        d = dict(zip(hdr, f))
        q = d.get("query") or f[0]
        clean = lambda k: (d.get(k, "") or "").strip()
        val = lambda k: "" if clean(k) in ("-", "") else clean(k)
        out[q] = dict(
            og=val("eggNOG_OGs").split(",")[0] if val("eggNOG_OGs") else "",
            desc=val("Description"),
            pname=val("Preferred_name"),
            go=set(x for x in val("GOs").split(",") if x.startswith("GO:")),
            ko=set(x for x in val("KEGG_ko").split(",") if x),
            path=set(x for x in val("KEGG_Pathway").split(",") if x),
            cog=val("COG_category"),
        )
    return out


def read_swissprot(path):
    """DIAMOND outfmt6 with stitle: q, s, stitle, pident, len, evalue, bits, qcov"""
    best = {}
    if not path or not os.path.exists(path):
        return best
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        try:
            bits = float(f[6])
        except ValueError:
            continue
        if f[0] in best and bits <= best[f[0]]["bits"]:
            continue
        title = f[2]
        m = SPROT_DESC.match(title)
        g = SPROT_GN.search(title)
        acc = f[1].split("|")[1] if "|" in f[1] else f[1]
        best[f[0]] = dict(acc=acc, desc=m.group(1) if m else "",
                          symbol=g.group(1) if g else "", bits=bits,
                          pident=f[3], evalue=f[5],
                          qcov=f[7] if len(f) > 7 else "")
    return best


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff3", required=True)
    p.add_argument("--interproscan")
    p.add_argument("--eggnog")
    p.add_argument("--swissprot")
    p.add_argument("--out-gff3", dest="out_gff3", required=True)
    p.add_argument("--out-tsv", dest="out_tsv", required=True)
    p.add_argument("--summary", required=True)
    a = p.parse_args()

    dom, ipr, ipr_go, pw = read_interproscan(a.interproscan)
    egg = read_eggnog(a.eggnog)
    sp = read_swissprot(a.swissprot)

    rows = []
    n_named = n_go = n_kegg = n_ipr = 0

    with open(a.out_gff3, "w") as o:
        for line in open(a.gff3):
            if line.startswith("#"):
                o.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] not in ("mRNA", "transcript"):
                o.write(line)
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            tid = m.group(1) if m else ""
            e = egg.get(tid, {})
            s = sp.get(tid, {})

            product = (s.get("desc") or e.get("desc")
                       or (dom[tid][0] if dom.get(tid) else "")
                       or "hypothetical protein")
            symbol = s.get("symbol") or e.get("pname") or ""
            gos = sorted(ipr_go.get(tid, set()) | e.get("go", set()))
            kos = sorted(e.get("ko", set()))
            iprs = sorted(ipr.get(tid, set()))
            paths = sorted(pw.get(tid, set()) | e.get("path", set()))

            if product != "hypothetical protein":
                n_named += 1
            if gos:
                n_go += 1
            if kos:
                n_kegg += 1
            if iprs:
                n_ipr += 1

            attrs = [f[8].rstrip(";")]
            attrs.append("product=%s" % esc(product))
            if symbol:
                attrs.append("gene_symbol=%s" % esc(symbol))
            xref = (["InterPro:%s" % x for x in iprs]
                    + (["UniProtKB:%s" % s["acc"]] if s.get("acc") else [])
                    + ["KEGG:%s" % k.replace("ko:", "") for k in kos])
            if xref:
                attrs.append("Dbxref=%s" % ",".join(xref))
            if gos:
                attrs.append("Ontology_term=%s" % ",".join(gos))
            f[8] = ";".join(attrs)
            o.write("\t".join(f) + "\n")

            rows.append([tid, symbol, product,
                         s.get("acc", ""), s.get("pident", ""), s.get("evalue", ""),
                         e.get("og", ""), e.get("cog", ""),
                         ";".join(iprs), ";".join(gos),
                         ";".join(kos), ";".join(paths),
                         ";".join(sorted(set(dom.get(tid, []))))])

    with open(a.out_tsv, "w") as t:
        t.write("\t".join([
            "transcript_id", "gene_symbol", "product", "swissprot_acc",
            "swissprot_pident", "swissprot_evalue", "eggnog_og", "cog_category",
            "interpro", "go_terms", "kegg_ko", "kegg_pathway", "domains"]) + "\n")
        for r in rows:
            t.write("\t".join(str(x) for x in r) + "\n")

    n = len(rows)
    pct = lambda x: (100.0 * x / n) if n else 0.0
    L = ["=== functional annotation coverage ===",
         "transcripts                    %d" % n,
         "with a product name            %d (%.1f%%)" % (n_named, pct(n_named)),
         "with GO terms                  %d (%.1f%%)" % (n_go, pct(n_go)),
         "with KEGG KO                   %d (%.1f%%)" % (n_kegg, pct(n_kegg)),
         "with InterPro domains          %d (%.1f%%)" % (n_ipr, pct(n_ipr)),
         "hypothetical only              %d (%.1f%%)" % (n - n_named, pct(n - n_named))]
    text = "\n".join(L)
    open(a.summary, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
