"""Compare per-record complexity values: Rust (rusted-garlic complexity TSV) vs the
SeqComplex.pm reference (tests/data/complexity/seqcomplex_reference.tsv).

Usage: python3 compare_complexity.py rust.tsv reference.tsv
Prints, per measure, the max |Rust - Perl| over records, how many records differ by
more than 1e-6 (Rust writes 6 decimals), the worst record, and how many reference
values are NA.
"""
import csv, sys
rust = {r['seq']: r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')}
perl = {r['seq']: r for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}
print(f"records: rust {len(rust)}, reference {len(perl)}; "
      f"missing in reference: {sorted(set(rust) - set(perl))}; missing in rust: {sorted(set(perl) - set(rust))}")
measures = [c for c in next(iter(perl.values())) if c != 'seq']
print(f"{'measure':8s} {'max|R-P|':>10s} {'n>1e-6':>7s}  worst")
for m in measures:
    diffs, na = [], 0
    for k, p in perl.items():
        if k not in rust: continue
        if p[m] in ('NA', ''): na += 1; continue
        diffs.append((abs(float(rust[k][m]) - float(p[m])), k))
    mx = max(diffs) if diffs else (float('nan'), '')
    n = sum(1 for d, _ in diffs if d > 1e-6)
    print(f"{m:8s} {mx[0]:10.6f} {n:7d}  {mx[1] if n else ''}{'  NA:' + str(na) if na else ''}")
