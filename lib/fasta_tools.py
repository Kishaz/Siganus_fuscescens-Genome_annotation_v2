#!/usr/bin/env python3
"""FASTA utilities for the Sfus_1.0 annotation pipeline.

Subcommands
-----------
map      build local<->NCBI<->GenBank name map from an NCBI assembly_report
rename   rewrite FASTA headers using a two-column mapping
verify   md5-compare two FASTAs (case-insensitive) through a mapping
stats    assembly statistics incl. softmasked fraction
"""
import argparse, gzip, hashlib, sys
from collections import OrderedDict


def opener(path):
    return gzip.open if str(path).endswith('.gz') else open


def read_fa(path):
    """Yield (name, description, sequence_bytes)."""
    op = opener(path)
    name = desc = None
    buf = []
    with op(path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, desc, b''.join(buf)
                parts = line[1:].rstrip().split(None, 1)
                name = parts[0].decode()
                desc = parts[1].decode() if len(parts) > 1 else ''
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, desc, b''.join(buf)


def read_report(path):
    """Parse an NCBI *_assembly_report.txt -> [(seq_name, genbank_accn, length)]."""
    out = []
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        out.append((f[0], f[4], int(f[8])))
    return out


def cmd_map(a):
    """Match a local FASTA to an NCBI assembly report.

    Primary key is sequence length; it is only accepted when every length in the
    assembly is unique, so the mapping cannot be silently wrong. Falls back to
    positional order (with a loud warning) otherwise.
    """
    loc = [(n, len(s)) for n, _, s in read_fa(a.fasta)]
    rep = read_report(a.report)
    if len(loc) != len(rep):
        sys.exit(f'ERROR: {len(loc)} local sequences vs {len(rep)} in report')

    lens_unique = len({L for _, L in loc}) == len(loc) and \
                  len({L for _, _, L in rep}) == len(rep)
    if lens_unique:
        by_len = {L: (n, acc) for n, acc, L in rep}
        pairs = [(n, *by_len[L]) for n, L in loc]
        method = 'unique-length'
    else:
        print('WARNING: sequence lengths are not unique; falling back to file '
              'order. Verify with `fasta_tools.py verify`.', file=sys.stderr)
        pairs = [(l[0], r[0], r[1]) for l, r in zip(loc, rep)]
        method = 'positional'

    lens = dict(loc)
    with open(a.out, 'w') as o:
        o.write('local_name\tncbi_seq_name\tgenbank_accn\tlength\n')
        for ln, sn, acc in pairs:
            o.write(f'{ln}\t{sn}\t{acc}\t{lens[ln]}\n')
    print(f'wrote {a.out} ({len(pairs)} sequences, method={method})')


def load_map(path, frm, to):
    hdr = None
    m = OrderedDict()
    for i, line in enumerate(open(path)):
        f = line.rstrip('\n').split('\t')
        if i == 0:
            hdr = f
            continue
        d = dict(zip(hdr, f))
        m[d[frm]] = d[to]
    return m


def cmd_rename(a):
    m = load_map(a.map, a.frm, a.to)
    n = miss = 0
    with open(a.out, 'w') as o:
        for name, _desc, seq in read_fa(a.fasta):
            if name not in m:
                miss += 1
                if a.strict:
                    sys.exit(f'ERROR: {name} not in mapping')
                new = name
            else:
                new = m[name]
            o.write(f'>{new}\n')
            for i in range(0, len(seq), 60):
                o.write(seq[i:i + 60].decode() + '\n')
            n += 1
    print(f'renamed {n} sequences -> {a.out}' + (f' ({miss} unmapped)' if miss else ''))


def cmd_verify(a):
    md5 = lambda s: hashlib.md5(s.upper()).hexdigest()
    A = {n: md5(s) for n, _, s in read_fa(a.a)}
    B = {n: md5(s) for n, _, s in read_fa(a.b)}
    m = load_map(a.map, a.frm, a.to) if a.map else {k: k for k in A}
    ok = bad = miss = 0
    bad_names = []
    for x, y in m.items():
        if x not in A or y not in B:
            miss += 1
            continue
        if A[x] == B[y]:
            ok += 1
        else:
            bad += 1
            bad_names.append((x, y))
    print(f'identical={ok}  mismatch={bad}  missing={miss}  (n={len(m)})')
    if bad_names:
        print('first mismatches:', bad_names[:10])
    sys.exit(0 if (bad == 0 and miss == 0 and ok == len(m)) else 1)


def cmd_stats(a):
    lens, soft, ns, gc, at = [], 0, 0, 0, 0
    for _n, _d, s in read_fa(a.fasta):
        lens.append(len(s))
        soft += sum(1 for c in s if 97 <= c <= 122)
        u = s.upper()
        ns += u.count(b'N')
        gc += u.count(b'G') + u.count(b'C')
        at += u.count(b'A') + u.count(b'T')
    tot = sum(lens)
    lens.sort(reverse=True)
    c = 0
    n50 = l50 = 0
    for i, L in enumerate(lens, 1):
        c += L
        if c >= tot / 2:
            n50, l50 = L, i
            break
    print(f'sequences      {len(lens)}')
    print(f'total_bp       {tot}')
    print(f'longest_bp     {lens[0]}')
    print(f'N50_bp         {n50}')
    print(f'L50            {l50}')
    print(f'N_bases        {ns}')
    print(f'GC_percent     {gc / (gc + at) * 100:.2f}')
    print(f'softmasked_bp  {soft}')
    print(f'softmasked_pct {soft / tot * 100:.2f}')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    m = sub.add_parser('map');    m.add_argument('--fasta', required=True)
    m.add_argument('--report', required=True); m.add_argument('--out', required=True)
    m.set_defaults(func=cmd_map)

    r = sub.add_parser('rename'); r.add_argument('--fasta', required=True)
    r.add_argument('--map', required=True); r.add_argument('--out', required=True)
    r.add_argument('--from', dest='frm', default='local_name')
    r.add_argument('--to', default='genbank_accn')
    r.add_argument('--strict', action='store_true')
    r.set_defaults(func=cmd_rename)

    v = sub.add_parser('verify'); v.add_argument('--a', required=True)
    v.add_argument('--b', required=True); v.add_argument('--map')
    v.add_argument('--from', dest='frm', default='local_name')
    v.add_argument('--to', default='genbank_accn')
    v.set_defaults(func=cmd_verify)

    s = sub.add_parser('stats');  s.add_argument('--fasta', required=True)
    s.set_defaults(func=cmd_stats)

    a = p.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
