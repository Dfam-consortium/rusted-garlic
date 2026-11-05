#!/usr/bin/env python3
"""
Compare inserted-element statistics: Perl GARLIC `.inserts` files vs. a
rusted-garlic `--out_bed` file.  Reports per-Mb counts and length distributions
by repeat class, and (Perl only) requested fragment length vs. inserted length.

Usage:
  compare_inserts.py --perl perl_*.inserts --rust rust.bed --mb_per_seq 5
"""
import argparse
from collections import defaultdict

import numpy as np

CLASSES = ["SINE/Alu", "SINE/MIR", "LINE/L1", "LINE/L2", "LTR", "DNA", "Simple_repeat", "other"]


def classify(family):
    if family == "SIMPLE" or family == "Simple_repeat":
        return "Simple_repeat"
    for c in CLASSES[:4]:
        if family.startswith(c):
            return c
    for c in ("LTR", "DNA"):
        if family.startswith(c):
            return c
    return "other"


def read_perl(paths):
    """Yield (class, length, requested_len or None)."""
    for p in paths:
        for line in open(p):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            ini, end, spec = int(f[0]), int(f[1]), f[3]
            spec = spec.split("[", 1)[0]
            parts = spec.split(":")
            if parts[0] == "SIMPLE":
                yield "Simple_repeat", end - ini, None
            else:
                req = parts[6] if len(parts) > 6 else None
                req = int(req) if req and req.isdigit() else None
                yield classify(parts[1]), end - ini, req


def read_rust(path):
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        fam = f[4] if len(f) > 4 and f[4] else "Simple_repeat"
        yield classify(fam), int(f[2]) - int(f[1]), None


def summarize(records, total_mb):
    by = defaultdict(list)
    for c, ln, _ in records:
        by[c].append(ln)
        by["ALL"].append(ln)
    out = {}
    for c, v in by.items():
        a = np.array(v)
        out[c] = (len(a) / total_mb, a.mean(), np.median(a), np.percentile(a, 90), a.sum() / total_mb / 1e4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perl", nargs="+", required=True)
    ap.add_argument("--rust", required=True)
    ap.add_argument("--mb_per_seq", type=float, default=5.0)
    args = ap.parse_args()

    perl = list(read_perl(args.perl))
    rust = list(read_rust(args.rust))
    n_rust_seqs = len({l.split("\t")[0] for l in open(args.rust)})
    sp = summarize(perl, len(args.perl) * args.mb_per_seq)
    sr = summarize(rust, n_rust_seqs * args.mb_per_seq)

    print("class\tperl_per_Mb\trust_per_Mb\tperl_mean\trust_mean\tperl_median\trust_median"
          "\tperl_p90\trust_p90\tperl_cov%\trust_cov%")
    for c in ["ALL"] + CLASSES:
        p = sp.get(c, (0, 0, 0, 0, 0))
        r = sr.get(c, (0, 0, 0, 0, 0))
        print(f"{c}\t{p[0]:.0f}\t{r[0]:.0f}\t{p[1]:.0f}\t{r[1]:.0f}\t{p[2]:.0f}\t{r[2]:.0f}"
              f"\t{p[3]:.0f}\t{r[3]:.0f}\t{p[4]:.2f}\t{r[4]:.2f}")

    req = np.array([(r, ln) for c, ln, r in perl if r is not None and c != "Simple_repeat"])
    if len(req):
        print(f"\n# Perl TE: requested fragment length mean {req[:, 0].mean():.0f}, "
              f"inserted length mean {req[:, 1].mean():.0f} "
              f"(inserted includes nested inserts and indels)")


if __name__ == "__main__":
    main()
