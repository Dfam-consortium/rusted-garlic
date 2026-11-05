#!/usr/bin/env python3
"""
compare_models.py — compare a Perl GARLIC model directory with a rusted_garlic model.

Usage:
  python3 compare_models.py <perl_model_prefix> <rust_model_prefix>

Example:
  python3 compare_models.py out/chr22/perl/data/chr22_perl/chr22_perl out/chr22/models/rust_compat

Files read:
  Perl:
    <perl_prefix>.kmer.K4.W1000.data
    <perl_prefix>.GCt.W1000.data
    <perl_prefix>.repeats.W1000.data
    <perl_prefix>.inserts.W1000.data
  Rust:
    <rust_prefix>.model          (JSON — repeat_model embedded within)
"""

import sys
import json
import math
from collections import defaultdict

BASES = ['A', 'C', 'G', 'T']

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def section(title):
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)

def subsection(title):
    print(f"\n  --- {title} ---")

def ok(msg):    print(f"    [OK]  {msg}")
def warn(msg):  print(f"  [WARN]  {msg}")
def err(msg):   print(f" [ERROR]  {msg}")

def rel_diff(a, b):
    """Relative difference, safe against both-zero."""
    if a == 0 and b == 0:
        return 0.0
    denom = max(abs(a), abs(b))
    return abs(a - b) / denom

# ─────────────────────────────────────────────────────────────────────────────
# Perl file parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_perl_gct(path):
    """
    Parse .GCt.W1000.data → dict of (from_label, to_label) → (freq, count).
    Also returns ordered list of bin labels.
    """
    trans = {}
    bins = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            from_b, to_b, freq = parts[0], parts[1], float(parts[2])
            count = int(parts[3]) if len(parts) > 3 else 0
            trans[(from_b, to_b)] = (freq, count)
            if from_b not in bins:
                bins.append(from_b)
    return bins, trans

def parse_perl_kmer(path):
    """
    Parse .kmer.K4.W1000.data
    Format:  #GC=<range>
             AAAA  0.39...  0.02...  2681
             ...
    Returns: dict of gc_label → dict of kmer → conditional_freq
    """
    data = {}
    current_bin = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('#GC='):
                current_bin = line[4:]
                data[current_bin] = {}
            elif current_bin is not None:
                parts = line.split()
                if len(parts) >= 2:
                    kmer = parts[0]
                    cond_freq = float(parts[1])
                    data[current_bin][kmer] = cond_freq
    return data

def parse_repeats_file(path):
    """
    Parse .repeats.W1000.data (both Perl and Rust have the same format)
    Returns: dict of gc_label → list of record dicts
    """
    data = {}
    current_bin = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('#GC='):
                current_bin = line[4:]
                data[current_bin] = []
            elif current_bin is not None:
                parts = line.split(':')
                if len(parts) >= 6 and parts[0] == 'SIMPLE':
                    data[current_bin].append({
                        'type': 'SIMPLE',
                        'unit': parts[1],
                        'strand': parts[2],
                        'length': float(parts[3]),
                        'div': float(parts[4]),
                        'ins_del': float(parts[5]),
                    })
                elif len(parts) >= 7:
                    data[current_bin].append({
                        'type': 'TE',
                        'name': parts[0],
                        'family': parts[1],
                        'strand': parts[2],
                        'div': float(parts[3]),
                        'ins': float(parts[4]),
                        'del': float(parts[5]),
                        'length': int(parts[6]),
                    })
    return data

def parse_inserts_file(path):
    """
    Parse .inserts.W1000.data
    Returns: dict of gc_label → list of (host, inserted, freq)
    """
    data = {}
    current_bin = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('#GC='):
                current_bin = line[4:]
                data[current_bin] = []
            elif current_bin is not None:
                parts = line.split('\t')
                if len(parts) >= 3:
                    data[current_bin].append((parts[0], parts[1], float(parts[2])))
    return data

# ─────────────────────────────────────────────────────────────────────────────
# Rust JSON model parser
# ─────────────────────────────────────────────────────────────────────────────

def load_rust_model(path):
    with open(path) as f:
        return json.load(f)

def rust_bin_label(r, precision=1):
    """Convert a gc_range dict {start, end} to a display label like '34.2-36.4'."""
    return f"{r['start']/10:.{precision}f}-{r['end']/10:.{precision}f}"

def rust_repeats_from_json(rust_model):
    """Extract repeat profile from the JSON model's repeat_model field.

    Returns the same format as parse_repeats_file():
        dict of gc_label → list of record dicts with 'type' in ('TE', 'SIMPLE').

    Handles both compact-array serialization (current code) and named-dict
    serialization (older model files).
    """
    result = {}
    rm = rust_model.get('repeat_model')
    if rm is None:
        return result

    gc_ranges    = rust_model['gc_ranges']
    te_profile   = rm.get('te_profile', [])
    simple_profile = rm.get('simple_profile', [])

    for i, r in enumerate(gc_ranges):
        label = rust_bin_label(r)
        records = []

        if i < len(te_profile):
            for rec in te_profile[i]:
                if isinstance(rec, list):
                    # current: [name, family, dir, div, ins, del, length, nfrg]
                    records.append({
                        'type': 'TE',
                        'name':   rec[0],
                        'family': rec[1],
                        'strand': rec[2],
                        'div':    float(rec[3]),
                        'ins':    float(rec[4]),
                        'del':    float(rec[5]),
                        'length': int(rec[6]),
                    })
                elif isinstance(rec, dict):
                    # older: named fields
                    records.append({
                        'type': 'TE',
                        'name':   rec['name'],
                        'family': rec['family'],
                        'strand': rec.get('dir', '+'),
                        'div':    float(rec['div']),
                        'ins':    float(rec['ins']),
                        'del':    float(rec['del']),
                        'length': int(rec['length']),
                    })

        if i < len(simple_profile):
            for rec in simple_profile[i]:
                if isinstance(rec, list):
                    # current: [consensus, dir, period, div, indel]
                    records.append({
                        'type':   'SIMPLE',
                        'unit':   rec[0],
                        'strand': rec[1],
                        'length': float(rec[2]),
                        'div':    float(rec[3]),
                        'ins_del': float(rec[4]),
                    })
                elif isinstance(rec, dict):
                    # older: named fields
                    records.append({
                        'type':   'SIMPLE',
                        'unit':   rec.get('consensus', rec.get('unit', '')),
                        'strand': rec.get('dir', '+'),
                        'length': float(rec.get('period', 0)),
                        'div':    float(rec.get('div', 0)),
                        'ins_del': float(rec.get('indel', 0)),
                    })

        result[label] = records

    return result


def rust_inserts_from_json(rust_model):
    """Extract insert table from the JSON model's repeat_model field.

    Returns the same format as parse_inserts_file():
        dict of gc_label → list of (host, inserted, freq) tuples.
    """
    result = {}
    rm = rust_model.get('repeat_model')
    if rm is None:
        return result

    gc_ranges    = rust_model['gc_ranges']
    insert_table = rm.get('insert_table', [])

    for i, r in enumerate(gc_ranges):
        label = rust_bin_label(r)
        pairs = []
        if i < len(insert_table):
            for host_key, inner in insert_table[i].items():
                for ins_key, freq in inner.items():
                    pairs.append((host_key, ins_key, float(freq)))
        result[label] = pairs

    return result


def rust_kmer_freqs(model):
    """
    Extract kmer frequencies from Rust JSON model.
    Returns: list of dicts (one per bin), each dict maps full kmer (e.g. 'AAAA') → cond_freq
    The JSON stores prefix → [P(A|prefix), P(C|prefix), P(G|prefix), P(T|prefix)].
    Reconstruct full kmer entries for comparison with the Perl format.
    """
    result = []
    for bin_map in model['gc_binned_kmer_flank_freqs']:
        kmer_dict = {}
        for prefix, freqs in bin_map.items():
            for i, base in enumerate(BASES):
                kmer = prefix + base
                kmer_dict[kmer] = freqs[i]
        result.append(kmer_dict)
    return result

# ─────────────────────────────────────────────────────────────────────────────
# Comparison functions
# ─────────────────────────────────────────────────────────────────────────────

def compare_gc_bins(perl_bins, rust_ranges):
    """Compare GC bin boundaries between Perl and Rust."""
    subsection("GC Bin Boundaries")

    n_perl = len(perl_bins)
    n_rust = len(rust_ranges)
    if n_perl != n_rust:
        warn(f"Bin count: Perl={n_perl}, Rust={n_rust}")
    else:
        ok(f"Bin count matches: {n_perl}")

    n = min(n_perl, n_rust)
    bad = 0
    for i in range(n):
        pb = perl_bins[i]          # e.g. "34.2-36.5"
        rb = rust_bin_label(rust_ranges[i])  # e.g. "34.2-36.4"
        try:
            p_start, p_end = [float(x) for x in pb.split('-')]
        except ValueError:
            # Handle e.g. "0-34.2"
            parts = pb.split('-')
            p_start, p_end = float(parts[0]), float(parts[1])
        r_start = rust_ranges[i]['start'] / 10.0
        r_end   = rust_ranges[i]['end']   / 10.0

        start_diff = abs(p_start - r_start)
        end_diff   = abs(p_end   - r_end)
        if start_diff > 0.15 or end_diff > 0.15:
            warn(f"  bin{i:02d}: Perl={pb:15s}  Rust={rb:15s}  "
                 f"(start_Δ={start_diff:.2f}, end_Δ={end_diff:.2f})")
            bad += 1
    if bad == 0:
        ok(f"All {n} bin boundaries agree within 0.15%")
    else:
        warn(f"{bad}/{n} bins have boundary differences > 0.15%")


def compare_gc_transitions(perl_gct, perl_bins, rust_model):
    """Compare the GC transition matrix."""
    subsection("GC Transition Matrix")

    n = len(perl_bins)
    rust_trans = rust_model['gc_trans']

    max_diff = 0.0
    big_diffs = []
    for i, from_b in enumerate(perl_bins):
        for j, to_b in enumerate(perl_bins):
            if i >= len(rust_trans) or j >= len(rust_trans[i]):
                continue
            perl_freq = perl_gct.get((from_b, to_b), (0.0, 0))[0]
            rust_freq = rust_trans[i][j]
            d = abs(perl_freq - rust_freq)
            if d > max_diff:
                max_diff = d
            if d > 0.01:
                big_diffs.append((i, j, from_b, to_b, perl_freq, rust_freq, d))

    if not big_diffs:
        ok(f"Transition matrix: max abs diff = {max_diff:.6f} (all within 0.01)")
    else:
        warn(f"Transition matrix: {len(big_diffs)} cell(s) differ by >0.01 "
             f"(max diff = {max_diff:.4f})")
        for i, j, fb, tb, pf, rf, d in big_diffs[:10]:
            print(f"      [{i},{j}] {fb:15s} → {tb:15s}  Perl={pf:.5f}  Rust={rf:.5f}  Δ={d:.5f}")
        if len(big_diffs) > 10:
            print(f"      ... and {len(big_diffs)-10} more")


def compare_kmer_freqs(perl_kmer_data, rust_model, perl_bins):
    """Compare per-bin k-mer conditional frequencies."""
    subsection("K-mer Conditional Frequencies")

    rust_kmer_list = rust_kmer_freqs(rust_model)
    n = min(len(perl_bins), len(rust_kmer_list))

    total_comparisons = 0
    diffs_01  = 0   # > 0.01
    diffs_001 = 0   # > 0.001
    max_diff  = 0.0
    max_info  = None

    for i in range(n):
        p_label = perl_bins[i]
        perl_bin = perl_kmer_data.get(p_label, {})
        rust_bin = rust_kmer_list[i]

        all_kmers = set(perl_bin) | set(rust_bin)
        for kmer in all_kmers:
            pf = perl_bin.get(kmer, 0.0)
            rf = rust_bin.get(kmer, 0.0)
            d = abs(pf - rf)
            total_comparisons += 1
            if d > max_diff:
                max_diff = d
                max_info = (i, p_label, kmer, pf, rf)
            if d > 0.01:
                diffs_01 += 1
            elif d > 0.001:
                diffs_001 += 1

    ok(f"Total k-mer comparisons: {total_comparisons:,}")
    ok(f"Max abs difference: {max_diff:.6f}")
    if max_info:
        i, lbl, kmer, pf, rf = max_info
        ok(f"  (bin{i:02d} '{lbl}' kmer '{kmer}': Perl={pf:.6f} Rust={rf:.6f})")
    if diffs_01 == 0 and diffs_001 == 0:
        ok("All k-mer frequencies agree within 0.001")
    else:
        if diffs_01:
            warn(f"{diffs_01:,} k-mer entries differ by >0.01")
        if diffs_001:
            warn(f"{diffs_001:,} k-mer entries differ by >0.001")

    # Summarise which bins have the most disagreement
    bin_maxdiff = []
    for i in range(n):
        p_label = perl_bins[i]
        perl_bin = perl_kmer_data.get(p_label, {})
        rust_bin = rust_kmer_list[i]
        local_max = 0.0
        for kmer in set(perl_bin) | set(rust_bin):
            d = abs(perl_bin.get(kmer, 0.0) - rust_bin.get(kmer, 0.0))
            if d > local_max:
                local_max = d
        bin_maxdiff.append((local_max, i, p_label))
    bin_maxdiff.sort(reverse=True)
    print(f"\n    Worst 5 bins by max k-mer diff:")
    for d, i, lbl in bin_maxdiff[:5]:
        print(f"      bin{i:02d} '{lbl}':  max_diff={d:.6f}")


def compare_repeat_profiles(perl_reps, rust_reps, perl_bins, rust_ranges):
    """Compare repeat (TE + SIMPLE) record counts and aggregate statistics."""
    subsection("Repeat Profile")

    print(f"\n    {'Bin':18s}  {'Perl_TE':>8s}  {'Rust_TE':>8s}  "
          f"{'Perl_SIM':>9s}  {'Rust_SIM':>9s}  {'TE_ratio':>9s}")
    print("    " + "-" * 68)

    total_p_te = total_r_te = total_p_sim = total_r_sim = 0
    for i, pb in enumerate(perl_bins):
        if i >= len(rust_ranges):
            break
        rb = rust_bin_label(rust_ranges[i])

        p_records = perl_reps.get(pb, [])
        r_records = rust_reps.get(rb, [])

        p_te  = sum(1 for r in p_records if r['type'] == 'TE')
        r_te  = sum(1 for r in r_records if r['type'] == 'TE')
        p_sim = sum(1 for r in p_records if r['type'] == 'SIMPLE')
        r_sim = sum(1 for r in r_records if r['type'] == 'SIMPLE')

        total_p_te  += p_te
        total_r_te  += r_te
        total_p_sim += p_sim
        total_r_sim += r_sim

        ratio = r_te / p_te if p_te else float('nan')
        ratio_str = f"{ratio:.2f}" if not math.isnan(ratio) else "  N/A"
        print(f"    {pb:18s}  {p_te:8d}  {r_te:8d}  "
              f"{p_sim:9d}  {r_sim:9d}  {ratio_str:>9s}")

    print("    " + "-" * 68)
    print(f"    {'TOTAL':18s}  {total_p_te:8d}  {total_r_te:8d}  "
          f"{total_p_sim:9d}  {total_r_sim:9d}")

    te_ratio  = total_r_te  / total_p_te  if total_p_te  else float('nan')
    sim_ratio = total_r_sim / total_p_sim if total_p_sim else float('nan')
    print(f"\n    Overall TE  coverage: Rust/Perl = {te_ratio:.3f}")
    print(f"    Overall SIM coverage: Rust/Perl = {sim_ratio:.3f}")

    # Aggregate TE divergence/ins/del distributions
    subsection("TE divergence / indel rate distributions (Perl vs Rust, all bins)")
    for field in ('div', 'ins', 'del'):
        p_vals = [r[field] for recs in perl_reps.values()
                  for r in recs if r['type'] == 'TE']
        r_vals = [r[field] for recs in rust_reps.values()
                  for r in recs if r['type'] == 'TE']
        if p_vals and r_vals:
            p_mean = sum(p_vals) / len(p_vals)
            r_mean = sum(r_vals) / len(r_vals)
            print(f"    TE {field:4s}:  Perl mean={p_mean:.3f}  Rust mean={r_mean:.3f}  "
                  f"Δ={abs(p_mean-r_mean):.3f}")


def compare_insert_tables(perl_ins, rust_ins, perl_bins, rust_ranges):
    """Compare insert (nested insertion) tables."""
    subsection("Nested Insert Table")

    total_p = sum(len(v) for v in perl_ins.values())
    total_r = sum(len(v) for v in rust_ins.values())
    print(f"    Perl total insert pairs: {total_p:,}")
    print(f"    Rust total insert pairs: {total_r:,}")

    # Check overlap of (host, inserted) pairs in each bin
    matched = missing_in_rust = extra_in_rust = 0
    for i, pb in enumerate(perl_bins):
        if i >= len(rust_ranges):
            break
        rb = rust_bin_label(rust_ranges[i])
        p_pairs = {(h, ins) for h, ins, _ in perl_ins.get(pb, [])}
        r_pairs = {(h, ins) for h, ins, _ in rust_ins.get(rb, [])}
        matched        += len(p_pairs & r_pairs)
        missing_in_rust += len(p_pairs - r_pairs)
        extra_in_rust   += len(r_pairs - p_pairs)

    total_perl_pairs = matched + missing_in_rust
    if total_perl_pairs > 0:
        pct = 100.0 * matched / total_perl_pairs
        print(f"    Pairs in both:         {matched:,} ({pct:.1f}% of Perl pairs)")
        if missing_in_rust:
            warn(f"Pairs in Perl but not Rust: {missing_in_rust:,}")
        if extra_in_rust:
            print(f"    Pairs in Rust but not Perl: {extra_in_rust:,}  (different nesting events)")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    perl_prefix = sys.argv[1]
    rust_prefix = sys.argv[2]

    # Load all files
    perl_gct_path    = f"{perl_prefix}.GCt.W1000.data"
    perl_kmer_path   = f"{perl_prefix}.kmer.K4.W1000.data"
    perl_rep_path    = f"{perl_prefix}.repeats.W1000.data"
    perl_ins_path    = f"{perl_prefix}.inserts.W1000.data"
    rust_model_path  = f"{rust_prefix}.model"

    print(f"Comparing:")
    print(f"  Perl prefix : {perl_prefix}")
    print(f"  Rust prefix : {rust_prefix}")

    print("\nLoading files...")
    perl_bins, perl_gct = parse_perl_gct(perl_gct_path)
    print(f"  Perl GCt     : {len(perl_bins)} bins, {len(perl_gct)} transitions")
    perl_kmer = parse_perl_kmer(perl_kmer_path)
    print(f"  Perl kmer    : {len(perl_kmer)} bins, "
          f"{sum(len(v) for v in perl_kmer.values())} k-mer entries")
    perl_reps = parse_repeats_file(perl_rep_path)
    print(f"  Perl repeats : {sum(len(v) for v in perl_reps.values())} records")
    perl_ins  = parse_inserts_file(perl_ins_path)
    print(f"  Perl inserts : {sum(len(v) for v in perl_ins.values())} pairs")
    rust_model = load_rust_model(rust_model_path)
    rust_ranges = rust_model['gc_ranges']
    print(f"  Rust model   : {len(rust_ranges)} bins")
    rust_trans = rust_model['gc_trans']
    n_rust_trans = sum(len(row) for row in rust_trans)
    print(f"  Rust GCt     : {len(rust_trans)} bins, {n_rust_trans} transitions")
    rust_kmer_bins = rust_model['gc_binned_kmer_flank_freqs']
    n_rust_kmers = sum(len(BASES) * len(b) for b in rust_kmer_bins)
    print(f"  Rust kmer    : {len(rust_kmer_bins)} bins, {n_rust_kmers} k-mer entries")
    rust_reps = rust_repeats_from_json(rust_model)
    print(f"  Rust repeats : {sum(len(v) for v in rust_reps.values())} records")
    rust_ins  = rust_inserts_from_json(rust_model)
    print(f"  Rust inserts : {sum(len(v) for v in rust_ins.values())} pairs")

    # ── 1. GC Bins ─────────────────────────────────────────────────────────
    section("1. GC BIN BOUNDARIES")
    compare_gc_bins(perl_bins, rust_ranges)

    # ── 2. GC Transition matrix ───────────────────────────────────────────
    section("2. GC TRANSITION MATRIX")
    compare_gc_transitions(perl_gct, perl_bins, rust_model)

    # ── 3. K-mer frequencies ───────────────────────────────────────────────
    section("3. K-MER CONDITIONAL FREQUENCIES")
    compare_kmer_freqs(perl_kmer, rust_model, perl_bins)

    # ── 4. Repeat profiles ─────────────────────────────────────────────────
    section("4. REPEAT PROFILE (TE + SIMPLE)")
    compare_repeat_profiles(perl_reps, rust_reps, perl_bins, rust_ranges)

    # ── 5. Insert table ────────────────────────────────────────────────────
    section("5. NESTED INSERT TABLE")
    compare_insert_tables(perl_ins, rust_ins, perl_bins, rust_ranges)

    print()


if __name__ == '__main__':
    main()
