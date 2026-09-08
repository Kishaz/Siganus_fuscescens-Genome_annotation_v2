#!/usr/bin/env python3
"""Compare annotation completeness against the assembly's own BUSCO ceiling.

A BUSCO score on a predicted proteome is easy to misread. The number that
actually matters is not "how complete is the annotation" in the abstract, but
"how much of what is demonstrably present in the assembly did the annotation
manage to capture". The assembly scored 98.5% complete in genome mode
(Mwamburi et al. 2024); a proteome scoring 92% has therefore recovered 93.4% of
what was there to find, and lost the rest.

Reporting recovery rather than the raw score makes a shortfall visible instead
of letting a respectable-looking absolute number hide it.
"""
import argparse
import os
import re
import sys

BUSCO_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%")


def parse_busco(path):
    if not path or not os.path.exists(path):
        return {}
    m = BUSCO_RE.search(open(path).read())
    if not m:
        return {}
    return dict(complete=float(m.group(1)), single=float(m.group(2)),
                duplicated=float(m.group(3)), fragmented=float(m.group(4)),
                missing=float(m.group(5)))


def parse_compleasm(path):
    """compleasm summary.txt: lines like 'S:95.12%, 3462' / 'D:0.55%, 20'."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    key = {"S": "single", "D": "duplicated", "F": "fragmented",
           "I": "incomplete", "M": "missing"}
    for line in open(path):
        m = re.match(r"\s*([SDFIM]):([\d.]+)%", line)
        if m:
            out[key[m.group(1)]] = float(m.group(2))
    if "single" in out:
        out["complete"] = out.get("single", 0.0) + out.get("duplicated", 0.0)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--busco")
    p.add_argument("--compleasm")
    p.add_argument("--assembly-busco-complete", type=float, default=98.5,
                   help="genome-mode completeness of the assembly itself")
    p.add_argument("--out", required=True)
    p.add_argument("--report", required=True)
    a = p.parse_args()

    b = parse_busco(a.busco)
    c = parse_compleasm(a.compleasm)
    if not b and not c:
        sys.exit("no completeness input could be parsed")

    rows = []
    L = ["=== completeness ===", ""]
    L.append("  %-14s %9s %9s %9s %9s" % ("source", "complete", "single",
                                          "duplicated", "missing"))
    for name, d in (("BUSCO", b), ("compleasm", c)):
        if not d:
            L.append("  %-14s %9s" % (name, "(absent)"))
            continue
        L.append("  %-14s %8.1f%% %8.1f%% %8.1f%% %8.1f%%"
                 % (name, d.get("complete", 0), d.get("single", 0),
                    d.get("duplicated", 0), d.get("missing", 0)))
        for k, v in d.items():
            rows.append(("%s_%s" % (name.lower(), k), v))

    asm = a.assembly_busco_complete
    rows.append(("assembly_busco_complete", asm))
    L.append("")
    L.append("  assembly ceiling (genome mode, 2024 paper): %.1f%%" % asm)

    best = max([d.get("complete", 0.0) for d in (b, c) if d] or [0.0])
    if asm > 0:
        rec = 100.0 * best / asm
        rows.append(("recovery_of_assembly_buscos_pct", round(rec, 1)))
        L.append("  annotation recovers %.1f%% of the assembly's BUSCO content" % rec)
        if rec < 95.0:
            L.append("  NOTE: >5%% of genes demonstrably present in the assembly are")
            L.append("        not represented in the annotation. Check the stage 10")
            L.append("        drop log before treating this release as final.")

    if b and c:
        gap = abs(b.get("complete", 0) - c.get("complete", 0))
        rows.append(("busco_compleasm_gap_pct", round(gap, 1)))
        if gap > 3.0:
            L.append("  NOTE: BUSCO and compleasm disagree by %.1f points. That is a"
                     % gap)
            L.append("        methodological artefact, not biology - report both.")

    with open(a.out, "w") as o:
        o.write("metric\tvalue\n")
        for k, v in rows:
            o.write("%s\t%s\n" % (k, v))
    text = "\n".join(L)
    open(a.report, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
