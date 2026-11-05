"""Compare simulated repeats with real TEs.

Simulated repeats are the lower-case runs of emitted FASTA files; real TEs are the
intervals of a TE annotation TSV (rusted-garlic --te_annot format) on the reference.
Reports, per group: GC %, CpG density (CG per kb), CpG o/e, frequency of CG-containing
4-mers, and the Jensen-Shannon divergence (bits) of the 4-mer spectrum vs real TEs.

Usage:
  python3 compare_repeats.py REF.fa TE_ANNOT.tsv LABEL=GLOB [LABEL=GLOB ...]
e.g.
  python3 compare_repeats.py chr22.fa chr22.fa.te_annot.tsv \\
      perl='emit/perl_*.fasta' rust=emit/rust.fa

Notes: lower-case runs include simple repeats (~4% of repeat bases on chr22).
CpG o/e depends on C and G content; CpG density and CG 4-mers do not, so use those
to judge CpG decay when repeat GC differs from real.
"""
import glob, re, sys
import numpy as np

IDX = {c: i for i, c in enumerate("ACGT")}


def fasta(path):
    buf = []
    for line in open(path):
        if line.startswith(">"):
            if buf:
                yield "".join(buf)
            buf = []
        else:
            buf.append(line.strip())
    if buf:
        yield "".join(buf)


def kmer4(seqs):
    h = np.zeros(256)
    for s in seqs:
        a = np.frombuffer(s.encode(), dtype=np.uint8)
        code = np.full(len(a), -1)
        for ch, i in IDX.items():
            code[a == ord(ch)] = i
        if len(a) < 4:
            continue
        ok = (code[:-3] >= 0) & (code[1:-2] >= 0) & (code[2:-1] >= 0) & (code[3:] >= 0)
        k = (code[:-3] * 64 + code[1:-2] * 16 + code[2:-1] * 4 + code[3:])[ok]
        h += np.bincount(k, minlength=256)
    return h / h.sum()


def jsd(p, q):
    m = (p + q) / 2
    kl = lambda x, y: np.sum(x[x > 0] * np.log2(x[x > 0] / y[x > 0]))
    return (kl(p, m) + kl(q, m)) / 2


CG_KMERS = [i for i in range(256)
            if "CG" in "".join("ACGT"[(i >> s) & 3] for s in (6, 4, 2, 0))]


def stats(regions, real=None):
    s = "".join(regions)
    n = sum(s.count(b) for b in "ACGT")
    c, g = s.count("C"), s.count("G")
    cg = sum(r.count("CG") for r in regions)
    k = kmer4(regions)
    out = dict(gc=100 * (c + g) / n, cg_kb=1000 * cg / n, oe=cg / (c * g / n),
               cg4=k[CG_KMERS].sum(), k4=k)
    out["jsd"] = jsd(k, real) if real is not None else 0.0
    return out


def main():
    ref = "".join(fasta(sys.argv[1])).upper()
    real_regions = []
    for line in open(sys.argv[2]):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        real_regions.append(ref[int(f[1]):int(f[2])])
    real = stats(real_regions)
    rows = [("real TEs", real)]
    for arg in sys.argv[3:]:
        label, pat = arg.split("=", 1)
        regions = [r.upper() for p in sorted(glob.glob(pat)) for s in fasta(p)
                   for r in re.findall(r"[acgt]+", s)]
        rows.append((label, stats(regions, real["k4"])))
    print("group\tGC%\tCG_per_kb\tCpG_oe\tCG_4mer_freq\t4mer_JSD_vs_real")
    for label, r in rows:
        print(f"{label}\t{r['gc']:.1f}\t{r['cg_kb']:.2f}\t{r['oe']:.3f}\t{r['cg4']:.4f}\t{r['jsd']:.5f}")


if __name__ == "__main__":
    main()
