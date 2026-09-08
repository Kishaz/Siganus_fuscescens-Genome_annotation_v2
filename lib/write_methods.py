#!/usr/bin/env python3
"""Emit METHODS.md and a software-version table from what the run actually used.

Version numbers written by hand drift away from what was executed. This reads
them out of the conda environments and run logs at packaging time, so the
methods text describes the run that produced the files sitting beside it.
"""
import argparse
import json
import os
import re
import subprocess
import sys
from datetime import date

KEY_PKGS = {"repeatmodeler", "repeatmasker", "tesorter", "hisat2", "stringtie",
            "fastp", "miniprot", "diamond", "agat", "busco", "omark",
            "trnascan-se", "infernal", "barrnap", "eggnog-mapper",
            "genometools-genometools", "samtools", "bedtools", "seqkit",
            "sra-tools", "gffread"}


def env_versions(mamba_bin, root, envs):
    rows = []
    for e in envs:
        try:
            out = subprocess.run([mamba_bin, "list", "-r", root, "-n", e, "--json"],
                                 capture_output=True, text=True, timeout=180)
            for pkg in json.loads(out.stdout or "[]"):
                rows.append((e, pkg.get("name", ""), pkg.get("version", ""),
                             pkg.get("channel", "")))
        except Exception as exc:
            print("  (could not read env %s: %s)" % (e, exc), file=sys.stderr)
    return rows


def cfg(path, key, default=""):
    """Pull a value out of config.sh without sourcing it."""
    try:
        txt = open(path).read()
    except OSError:
        return default
    m = re.search(r'^export %s="\$\{%s:-([^"}]*)\}?"' % (key, key), txt, re.M)
    if m:
        return m.group(1)
    m = re.search(r'^export %s="([^"]*)"' % key, txt, re.M)
    return m.group(1) if m else default


def embed(lines, title, path):
    if not os.path.exists(path):
        return
    lines.append("## " + title)
    lines.append("")
    lines.append("```")
    lines.append(open(path).read().rstrip())
    lines.append("```")
    lines.append("")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--versions", required=True)
    a = p.parse_args()

    proj = cfg(a.config, "PROJ")
    mamba = cfg(a.config, "MICROMAMBA_BIN") or os.path.join(proj, "bin", "micromamba")
    root = cfg(a.config, "MAMBA_ROOT_PREFIX") or os.path.join(proj, "micromamba")
    envs = ["prep", "repeats", "teclass", "rnaseq", "annot", "ncrna", "func", "qc"]

    rows = env_versions(mamba, root, envs) if os.path.exists(mamba) else []
    with open(a.versions, "w") as v:
        v.write("environment\tpackage\tversion\tchannel\n")
        for r in rows:
            v.write("\t".join(r) + "\n")

    key = sorted({(n, ver) for _e, n, ver, _c in rows if n.lower() in KEY_PKGS})
    asm = cfg(a.config, "ASM_NAME", "Sfus_1.0")
    ver = cfg(a.config, "ANNOT_VERSION", "v2.0")
    acc = cfg(a.config, "ASM_ACC", "GCA_042847805.1")

    L = []
    L.append("# Methods - %s annotation %s" % (asm, ver))
    L.append("")
    L.append("Generated %s from the run that produced the files in this directory."
             % date.today())
    L.append("")
    L.append("## Assembly")
    L.append("")
    L.append("%s (%s), Siganus fuscescens isolate JG2022: 164 contigs, "
             "498,171,360 bp, N50 7.16 Mb, GC 42.92%%. Downloaded from GenBank "
             "and verified against the md5 checksums NCBI publishes alongside it."
             % (acc, asm))
    L.append("")
    L.append("## Software versions")
    L.append("")
    if key:
        L.append("| tool | version |")
        L.append("|---|---|")
        for n, v in key:
            L.append("| %s | %s |" % (n, v))
    else:
        L.append("_(environments unreadable at packaging time; "
                 "see software_versions.tsv)_")
    L.append("")
    L.append("BRAKER3 was run from the official `teambraker/braker3` container.")
    L.append("Full package inventory for every environment: `software_versions.tsv`.")
    L.append("")
    L.append("## Databases")
    L.append("")
    L.append("- Dfam %s, curated consensus partition" % cfg(a.config, "DFAM_RELEASE", "4.0"))
    L.append("- OrthoDB v12, Vertebrata partition")
    L.append("- Rfam, current release at run time")
    L.append("- UniProtKB/Swiss-Prot")
    L.append("- BUSCO lineage %s" % cfg(a.config, "BUSCO_LINEAGE", "actinopterygii_odb10"))
    L.append("- Reference proteomes resolved at run time from NCBI Datasets; the "
             "exact accessions used are in `ref_proteomes.tsv`.")
    L.append("")
    L.append("## Transcript evidence")
    L.append("")
    L.append("Conspecific RNA-seq only: %s. Cross-species data from congeners "
             "(notably the 80 public S. canaliculatus RNA-seq runs) was "
             "deliberately excluded so that provenance is unambiguous. Per-model "
             "transcript support is reported in `evidence_support.tsv`."
             % cfg(a.config, "RNASEQ_RUNS", "SRR36642802"))
    L.append("")

    embed(L, "Repeat landscape", os.path.join(a.work, "repeats", "final", "repeat_summary.tsv"))
    embed(L, "Filtering outcome", os.path.join(a.work, "filter", "filter_summary.txt"))
    embed(L, "Functional coverage", os.path.join(a.work, "functional", "functional_summary.txt"))
    embed(L, "Validation", os.path.join(a.work, "qc", "annotation_report.txt"))

    open(a.out, "w").write("\n".join(L) + "\n")
    print("wrote %s and %s (%d packages)" % (a.out, a.versions, len(rows)))


if __name__ == "__main__":
    main()
