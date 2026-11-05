#!/usr/bin/env python3
"""
Tier 2 validation: compare statistical properties of sequences emitted by
GARLIC 1.4 (Perl) and rusted-garlic against the real sequence the models were
built from, and against each other.

Each FASTA record of an emitted file is treated as one replicate.

Metrics (per replicate, 1 kb non-overlapping windows):
  - window GC distribution: mean, sd, JSD of 1%-bin histogram vs. reference
  - GC autocorrelation at lag 1 and lag 10 windows (captures GC transitions)
  - 4-mer spectra per GC decile bin (deciles from the reference), JSD vs.
    reference: all bases, and background only (uppercase / unmasked bases)
  - CpG observed/expected
  - lowercase (inserted repeat) fraction and run-length summary
  - rusted-garlic `complexity` measures: KS distance of window distributions

Usage:
  compare_emitted.py --ref chr22.fa --ref_masked chr22.fa.masked \
      --group perl perl_*.fasta --group rust rust.fa \
      [--complexity path/to/rusted-garlic] [--target_te 44.3]
"""
import argparse
import subprocess
import sys
from collections import defaultdict

import numpy as np

WIN = 1000
K = 4
MIN_ACGT = 0.7

ENC = np.full(256, -1, dtype=np.int8)
for i, b in enumerate(b"ACGT"):
    ENC[b] = i
    ENC[b + 32] = i  # lowercase


def read_fasta(path):
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)


def jsd(p, q):
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(p > 0, p * np.log2(p / m), 0.0).sum()
        b = np.where(q > 0, q * np.log2(q / m), 0.0).sum()
    return 0.5 * (a + b)


def ks(x, y):
    x, y = np.sort(x), np.sort(y)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    allv = np.concatenate([x, y])
    cx = np.searchsorted(x, allv, side="right") / len(x)
    cy = np.searchsorted(y, allv, side="right") / len(y)
    return float(np.max(np.abs(cx - cy)))


class SeqStats:
    """Per-record statistics. `kmer_seq` overrides the bases used for the
    background 4-mer spectrum (used for the hard-masked reference)."""

    def __init__(self, raw, kmer_raw=None):
        arr = np.frombuffer(raw, dtype=np.uint8)
        code = ENC[arr]
        nwin = len(arr) // WIN
        n = nwin * WIN
        code, arr = code[:n], arr[:n]
        cw = code.reshape(nwin, WIN)
        acgt = (cw >= 0).sum(1)
        gcn = ((cw == 1) | (cw == 2)).sum(1)
        self.valid = acgt >= MIN_ACGT * WIN
        self.gc = np.where(self.valid, 100.0 * gcn / np.maximum(acgt, 1), np.nan)
        self.code = code
        self.nwin = nwin
        is_lower = (arr >= 97) & (arr <= 122)
        self.n_acgt = int((code >= 0).sum())
        self.lower_frac = float((is_lower & (code >= 0)).sum() / max(self.n_acgt, 1))
        # lowercase runs
        d = np.diff(np.concatenate([[0], is_lower.astype(np.int8), [0]]))
        starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
        self.lower_runs = ends - starts
        # kmer index per position (-1 if any base invalid)
        self.kidx = self._kmer_index(code)
        if kmer_raw is not None:
            bg_code = ENC[np.frombuffer(kmer_raw, dtype=np.uint8)][:n]
            self.bg_kidx = self._kmer_index(bg_code)
        else:
            bg_code = np.where(is_lower, -1, code).astype(np.int8)
            self.bg_kidx = self._kmer_index(bg_code)
        # CpG o/e
        c = code
        ncg = int(((c[:-1] == 1) & (c[1:] == 2)).sum())
        nc, ng = int((c == 1).sum()), int((c == 2).sum())
        self.cpg_oe = ncg * self.n_acgt / max(nc * ng, 1)
        lc = is_lower & (c >= 0)
        lcg = int(((c[:-1] == 1) & (c[1:] == 2) & lc[:-1] & lc[1:]).sum())
        lnc, lng = int(((c == 1) & lc).sum()), int(((c == 2) & lc).sum())
        self.cpg_oe_rep = lcg * int(lc.sum()) / max(lnc * lng, 1)

    @staticmethod
    def _kmer_index(code):
        n = len(code)
        idx = np.zeros(n - K + 1, dtype=np.int32)
        ok = np.ones(n - K + 1, dtype=bool)
        for j in range(K):
            cj = code[j:n - K + 1 + j]
            idx = idx * 4 + np.maximum(cj, 0)
            ok &= cj >= 0
        return np.where(ok, idx, -1)

    def gc_hist(self):
        g = self.gc[self.valid]
        return np.histogram(g, bins=100, range=(0, 100))[0].astype(float) + 1e-9

    def autocorr(self, lag):
        x, y = self.gc[:-lag], self.gc[lag:]
        ok = ~np.isnan(x) & ~np.isnan(y)
        if ok.sum() < 3:
            return float("nan")
        return float(np.corrcoef(x[ok], y[ok])[0, 1])

    def kmer_by_bin(self, edges, background=False):
        """Return (nbins, 4**K) count matrix of k-mers, binned by window GC."""
        kidx = self.bg_kidx if background else self.kidx
        nb = len(edges) + 1
        wbin = np.where(self.valid, np.searchsorted(edges, self.gc, side="right"), -1)
        pos_bin = np.repeat(wbin, WIN)[: len(kidx)]
        ok = (kidx >= 0) & (pos_bin >= 0)
        comb = pos_bin[ok].astype(np.int64) * 4**K + kidx[ok]
        return np.bincount(comb, minlength=nb * 4**K).reshape(nb, 4**K).astype(float)


def complexity_table(binary, path):
    """Run `rusted-garlic complexity -w 1000` and return {seq: {measure: array}}."""
    out = subprocess.run([binary, "--log-level", "error", "complexity", "-w", str(WIN), path],
                         check=True, capture_output=True, text=True).stdout.splitlines()
    hdr = out[0].split("\t")
    meas = hdr[6:]
    data = defaultdict(lambda: defaultdict(list))
    for line in out[1:]:
        f = line.split("\t")
        for m, v in zip(meas, f[6:]):
            data[f[1]][m].append(float(v))
    return meas, {s: {m: np.array(v) for m, v in d.items()} for s, d in data.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="real sequence the model was built from")
    ap.add_argument("--ref_masked", help="hard-masked (N) reference for background 4-mers")
    ap.add_argument("--group", nargs="+", action="append", required=True, metavar=("LABEL", "FASTA"),
                    help="label followed by one or more emitted FASTA files")
    ap.add_argument("--complexity", help="rusted-garlic binary; enables complexity comparison")
    ap.add_argument("--target_te", type=float, help="requested TE %% (for reference in report)")
    ap.add_argument("--ndec", type=int, default=10, help="number of GC quantile bins [10]")
    args = ap.parse_args()

    print(f"# Reference: {args.ref}", flush=True)
    ref_recs = list(read_fasta(args.ref))
    mask_recs = dict(read_fasta(args.ref_masked)) if args.ref_masked else {}
    ref_raw = b"".join(r for _, r in ref_recs)
    ref_mask = b"".join(mask_recs[n] for n, _ in ref_recs) if mask_recs else None
    ref = SeqStats(ref_raw, ref_mask)
    if ref_mask is None:
        ref.bg_kidx = ref.kidx
    edges = np.nanpercentile(ref.gc[ref.valid], np.linspace(0, 100, args.ndec + 1)[1:-1])
    ref_hist = ref.gc_hist()
    ref_k = ref.kmer_by_bin(edges) + 0.5
    ref_kbg = ref.kmer_by_bin(edges, background=True) + 0.5
    print(f"# GC bin edges (ref deciles): {np.round(edges, 1).tolist()}")

    groups = {}
    for g in args.group:
        label, files = g[0], g[1:]
        reps = []
        for f in files:
            for name, raw in read_fasta(f):
                reps.append((f, name, SeqStats(raw)))
        groups[label] = reps
        print(f"# Group {label}: {len(reps)} replicate(s) from {len(files)} file(s)", flush=True)

    # ---- per-replicate table ----
    cols = ["len_Mb", "lower%", "runs/Mb", "run_mean", "gc_mean", "gc_sd", "gcJSD",
            "ac1", "ac10", "cpg_oe", "cpg_oe_rep", "k4JSD_mean", "k4JSD_max", "bgk4JSD_mean", "bgk4JSD_max"]

    def row(s):
        k = s.kmer_by_bin(edges) + 0.5
        kb = s.kmer_by_bin(edges, background=True) + 0.5
        kj = [jsd(k[i], ref_k[i]) for i in range(len(k))]
        kbj = [jsd(kb[i], ref_kbg[i]) for i in range(len(kb))]
        g = s.gc[s.valid]
        return [s.n_acgt / 1e6, 100 * s.lower_frac,
                len(s.lower_runs) / (s.n_acgt / 1e6), float(np.mean(s.lower_runs)) if len(s.lower_runs) else 0,
                float(np.mean(g)), float(np.std(g)), jsd(s.gc_hist(), ref_hist),
                s.autocorr(1), s.autocorr(10), s.cpg_oe, s.cpg_oe_rep,
                float(np.mean(kj)), float(np.max(kj)), float(np.mean(kbj)), float(np.max(kbj))]

    print("\n## Per-replicate metrics (JSD in bits, vs. reference)")
    if args.target_te is not None:
        print(f"# requested TE fraction: {args.target_te}%")
    print("\t".join(["group", "replicate"] + cols))
    ref_row = row(ref)
    ref_row[2] = ref_row[3] = ref_row[10] = float("nan")  # reference is not soft-masked
    print("\t".join(["ref", "-"] + [f"{v:.4g}" for v in ref_row]))
    summary = {}
    for label, reps in groups.items():
        rows = []
        for f, name, s in reps:
            r = row(s)
            rows.append(r)
            print("\t".join([label, name] + [f"{v:.4g}" for v in r]))
        summary[label] = np.array(rows)

    print("\n## Group summary (mean ± sd over replicates)")
    print("\t".join(["metric", "ref"] + list(summary)))
    for j, c in enumerate(cols):
        cells = [f"{ref_row[j]:.4g}"]
        for label, a in summary.items():
            cells.append(f"{a[:, j].mean():.4g} ± {a[:, j].std():.2g}")
        print("\t".join([c] + cells))

    # ---- group-vs-group equivalence ----
    # Distributions: mean JSD over cross-group replicate pairs divided by the mean
    # JSD over within-group pairs.  ~1 means the groups differ no more than two
    # runs of the same program do (both tools are stochastic); >1 is a real
    # difference.  Scalars: difference of group means in standard errors (z).
    labels = list(groups)
    if len(labels) >= 2:
        dists = {}
        for label, reps in groups.items():
            dists[label] = [(s.gc_hist(), s.kmer_by_bin(edges) + 0.5,
                             s.kmer_by_bin(edges, background=True) + 0.5) for _, _, s in reps]

        def dist_jsd(x, y, which):
            if which == 0:
                return jsd(x[0], y[0])
            return float(np.mean([jsd(x[which][i], y[which][i]) for i in range(len(x[which]))]))

        names = ["gc_hist", "k4_all", "k4_background"]
        print("\n## Equivalence: cross-group vs within-group replicate distance")
        print("# ratio = mean JSD(group A run, group B run) / mean JSD(two runs of the same group)")
        print("# ~1: indistinguishable from run-to-run noise.  z = |mean_A - mean_B| / SE")
        print("pair\tmetric\tcross\twithin\tratio\tverdict")
        for ia in range(len(labels)):
            for ib in range(ia + 1, len(labels)):
                A, B = dists[labels[ia]], dists[labels[ib]]
                pair = f"{labels[ia]}_vs_{labels[ib]}"
                for w, nm in enumerate(names):
                    cross = [dist_jsd(x, y, w) for x in A for y in B]
                    within = [dist_jsd(G[i], G[j], w) for G in (A, B)
                              for i in range(len(G)) for j in range(i + 1, len(G))]
                    if not within:
                        continue
                    r = np.mean(cross) / np.mean(within)
                    print(f"{pair}\t{nm}\t{np.mean(cross):.3g}\t{np.mean(within):.3g}\t{r:.2f}"
                          f"\t{'OK' if r < 1.5 else 'DIFF'}")
                sa, sb = summary[labels[ia]], summary[labels[ib]]
                for j, c in enumerate(cols):
                    if c in ("len_Mb", "gcJSD") or "JSD" in c:
                        continue
                    se = np.sqrt(sa[:, j].var(ddof=1) / len(sa) + sb[:, j].var(ddof=1) / len(sb))
                    d = abs(sa[:, j].mean() - sb[:, j].mean())
                    z = d / se if se > 0 else float("inf") if d > 0 else 0.0
                    print(f"{pair}\t{c}\t{sa[:, j].mean():.4g}\t{sb[:, j].mean():.4g}\tz={z:.1f}"
                          f"\t{'OK' if z < 3 else 'DIFF'}")

        pooled = {}
        for label, reps in groups.items():
            pooled[label] = (None, sum(d[1] for d in dists[label]), None)

        print("\n## Per-GC-bin 4-mer JSD vs. reference (pooled, all bases)")
        print("bin\tGC_range\t" + "\t".join(labels))
        lo = [0.0] + list(edges)
        hi = list(edges) + [100.0]
        for i in range(len(lo)):
            cells = [f"{jsd(pooled[l][1][i], ref_k[i]):.4g}" for l in labels]
            print(f"{i}\t{lo[i]:.1f}-{hi[i]:.1f}\t" + "\t".join(cells))

    # ---- complexity ----
    if args.complexity:
        print("\n## Complexity measures (1 kb windows): mean and KS distance", flush=True)
        meas, rt = complexity_table(args.complexity, args.ref)
        ref_c = {m: np.concatenate([d[m] for d in rt.values()]) for m in meas}
        grp_c = {}
        for label, reps in groups.items():
            acc = defaultdict(list)
            for f in sorted({f for f, _, _ in reps}):
                _, t = complexity_table(args.complexity, f)
                for d in t.values():
                    for m in meas:
                        acc[m].append(d[m])
            grp_c[label] = {m: np.concatenate(v) for m, v in acc.items()}
        hdr = ["measure", "ref_mean"] + [f"{l}_mean" for l in labels] + [f"KS_{l}_ref" for l in labels]
        if len(labels) >= 2:
            hdr.append(f"KS_{labels[0]}_{labels[1]}")
        print("\t".join(hdr))
        for m in meas:
            cells = [m, f"{ref_c[m].mean():.4g}"]
            cells += [f"{grp_c[l][m].mean():.4g}" for l in labels]
            cells += [f"{ks(grp_c[l][m], ref_c[m]):.3f}" for l in labels]
            if len(labels) >= 2:
                cells.append(f"{ks(grp_c[labels[0]][m], grp_c[labels[1]][m]):.3f}")
            print("\t".join(cells))


if __name__ == "__main__":
    sys.exit(main())
