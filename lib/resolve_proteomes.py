#!/usr/bin/env python3
"""Resolve reference fish proteomes to concrete NCBI FTP protein-FASTA URLs.

Accessions in published pipelines rot. This asks NCBI Datasets for the current
annotated assembly of each species instead of hard-coding, and emits a TSV that
the download step consumes. Only ftp.ncbi.nlm.nih.gov and api.ncbi.nlm.nih.gov
are contacted - both are reachable from the cluster.
"""
import argparse, json, sys, time, urllib.error, urllib.parse, urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

# Percomorph-weighted panel: close relatives first, models for coverage.
DEFAULT_SPECIES = [
    "Siganus canaliculatus",     # congener - closest available annotated genome
    "Dicentrarchus labrax",
    "Sparus aurata",
    "Lates calcarifer",
    "Larimichthys crocea",
    "Gasterosteus aculeatus",
    "Takifugu rubripes",
    "Oryzias latipes",
    "Danio rerio",               # model, deepest functional annotation
]


def get(url, tries=5):
    """GET with backoff - the Datasets API rate-limits unauthenticated clients."""
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "sfus-annotation/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(2 ** i)
                continue
            raise
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise last


def ftp_dir(acc):
    """GCF_000002035.6 -> .../GCF/000/002/035/"""
    pre, num = acc.split("_")
    num = num.split(".")[0]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{pre}/{num[0:3]}/{num[3:6]}/{num[6:9]}"


def pick(species):
    """Prefer RefSeq + annotated + reference; fall back to any annotated assembly."""
    q = urllib.parse.quote(species)
    url = (f"{API}/genome/taxon/{q}/dataset_report"
           f"?filters.has_annotation=true&page_size=60")
    try:
        d = get(url)
    except Exception as e:
        print(f"  ! query failed for {species}: {e}", file=sys.stderr)
        return None
    reps = d.get("reports", [])
    if not reps:
        return None

    def score(r):
        ai = r.get("assembly_info", {})
        s = 0
        if r["accession"].startswith("GCF_"):        s += 8
        if ai.get("refseq_category") == "reference genome": s += 4
        lvl = ai.get("assembly_level", "")
        s += {"Chromosome": 3, "Complete Genome": 3, "Scaffold": 1}.get(lvl, 0)
        return s

    best = max(reps, key=score)
    acc = best["accession"]
    name = best.get("assembly_info", {}).get("assembly_name", "").replace(" ", "_")
    return {
        "species": species,
        "accession": acc,
        "assembly_name": name,
        "level": best.get("assembly_info", {}).get("assembly_level"),
        "annotation": best.get("annotation_info", {}).get("name", ""),
        "url": f"{ftp_dir(acc)}/{acc}_{name}/{acc}_{name}_protein.faa.gz",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--species", nargs="*", default=DEFAULT_SPECIES)
    p.add_argument("--check", action="store_true",
                   help="HEAD each URL and report size (slower, needs network)")
    a = p.parse_args()

    rows = []
    for sp in a.species:
        r = pick(sp)
        if not r:
            print(f"  ! no annotated assembly found for {sp}", file=sys.stderr)
            continue
        if a.check:
            try:
                req = urllib.request.Request(r["url"], method="HEAD")
                with urllib.request.urlopen(req, timeout=60) as h:
                    r["size_mb"] = round(int(h.headers.get("Content-Length", 0)) / 1048576, 1)
            except Exception:
                r["size_mb"] = "UNREACHABLE"
        time.sleep(0.5)   # stay under the API rate limit
        rows.append(r)
        print(f"  {sp:26s} {r['accession']:18s} {r.get('level',''):16s} "
              f"{r.get('size_mb','')}", file=sys.stderr)

    cols = ["species", "accession", "assembly_name", "level", "annotation", "url"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"wrote {a.out} ({len(rows)} proteomes)", file=sys.stderr)


if __name__ == "__main__":
    main()
