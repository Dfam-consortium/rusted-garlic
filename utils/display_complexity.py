#!/usr/bin/env python3
"""
display_complexity.py - Display complexity stats from rusted-garlic complexity TSV output.

Port of SeqComplexity/displayStats.pl to read the TSV format produced by
`rusted-garlic complexity` directly.

Usage:
    display_complexity.py -in <complexity.tsv> [-graph] [-bins N] > output

Without --graph: prints a summary table (mean ± stdev per method per file).
With --graph:    outputs a self-contained Google Charts HTML file showing
                 histograms for each complexity metric.

Author:  Ported from displayStats.pl (Robert Hubley / Juan Caballero, ISB 2010-2015)
License: GNU General Public License v3
"""

import sys
import csv
import math
import argparse
from collections import defaultdict

# ---------------------------------------------------------------------------
# Method metadata — mirrors %methodMetaData in displayStats.pl
# Keys with a trailing digit (cm1-cm6 etc.) share the prefix entry.
# ---------------------------------------------------------------------------
METHOD_META = {
    "at":  {"title": "Percent AT",                        "x_units": "AT%"},
    "gc":  {"title": "Percent GC",                        "x_units": "GC%"},
    "ats": {"title": "AT Skew",                           "x_units": "(A-T)/(A+T)"},
    "gcs": {"title": "GC Skew",                           "x_units": "(G-C)/(G+C)"},
    "ce":  {"title": "Shannon Entropy",                   "x_units": "bits"},
    "cm":  {"title": "Markov Model Complexity",           "x_units": "entropy"},
    "cl":  {"title": "Linguistic Complexity",             "x_units": "word fraction"},
    "ct":  {"title": "Trifonov Complexity",               "x_units": "word fraction"},
    "cz":  {"title": "Compression Factor ( gzip )",       "x_units": "compression factor"},
    "cwf": {"title": "Wootton and Federhen Complexity",   "x_units": "complexity"},
    "cpg": {"title": "Percent CpGs",                      "x_units": "CpG%"},
    "ket": {"title": "Keto Skew",                         "x_units": "(G+T)/(A+C)"},
    "pur": {"title": "Purine Skew",                       "x_units": "(A+G)/(C+T)"},
}

# Columns to skip — bookkeeping, not complexity metrics
SKIP_COLS = {"file", "seq", "start", "end", "win_len", "acgt_count"}

# Methods whose x-axis is displayed as a percentage (multiply by 100)
PCT_METHODS = {"at", "gc", "cpg"}


def method_prefix(col):
    """Return the prefix and optional word-size for a column name, e.g. 'cm3' -> ('cm', 3)."""
    i = len(col)
    while i > 0 and col[i - 1].isdigit():
        i -= 1
    prefix = col[:i]
    word_size = int(col[i:]) if i < len(col) else 0
    return prefix, word_size


def read_tsv(paths):
    """
    Read one or more complexity TSV files and return:
      raw_data[method][file] = [float, float, ...]
    When multiple files are given, series from all files are merged into the
    same dict (keyed by the 'file' column value inside each TSV).
    """
    raw_data = defaultdict(lambda: defaultdict(list))
    for path in paths:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                fname = row["file"]
                for col, val in row.items():
                    if col in SKIP_COLS:
                        continue
                    if val == "" or val is None:
                        continue
                    try:
                        raw_data[col][fname].append(float(val))
                    except ValueError:
                        pass
    return raw_data


def stdev(values, mean):
    if len(values) < 2:
        return 0.0
    variance = sum((v - mean) ** 2 for v in values) / (len(values) - 1)
    return math.sqrt(variance)


def print_table(raw_data):
    """Print a summary table: method, file, n, mean, stdev, min, max."""
    header = f"{'method':<8}  {'file':<40}  {'n':>8}  {'mean':>12}  {'stdev':>12}  {'min':>12}  {'max':>12}"
    print(header)
    print("-" * len(header))
    for method in sorted(raw_data):
        for fname in sorted(raw_data[method]):
            vals = raw_data[method][fname]
            n = len(vals)
            mn = sum(vals) / n
            sd = stdev(vals, mn)
            print(
                f"{method:<8}  {fname:<40}  {n:>8}  {mn:>12.6f}  {sd:>12.6f}"
                f"  {min(vals):>12.6f}  {max(vals):>12.6f}"
            )


def bin_data(raw_data, n_bins):
    """
    Bin each method's values across all files.
    Returns:
      bin_counts[method][file][bin_index] = count
      bin_meta[method]  = {"min": float, "bin_size": float}
    """
    # Find global min/max per method across all files
    ranges = {}
    for method, file_data in raw_data.items():
        all_vals = [v for vals in file_data.values() for v in vals]
        if not all_vals:
            continue
        lo, hi = min(all_vals), max(all_vals)
        bin_size = (hi - lo) / n_bins if hi != lo else 0.0
        if bin_size == 0.0:
            # All values identical — mirror Perl logic
            bin_size = lo / n_bins if lo != 0 else 1.0 / n_bins
            lo = lo - (n_bins * bin_size)
        ranges[method] = {"min": lo, "bin_size": bin_size}

    bin_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for method, file_data in raw_data.items():
        if method not in ranges:
            continue
        lo = ranges[method]["min"]
        bs = ranges[method]["bin_size"]
        for fname, vals in file_data.items():
            for v in vals:
                if bs == 0:
                    b = 0
                else:
                    b = int(round((v - lo) / bs))
                b = max(0, min(b, n_bins))  # clamp
                bin_counts[method][fname][b] += 1

    return bin_counts, ranges


def js_safe(name):
    """Make a JavaScript-safe identifier from a column name."""
    return name.replace("-", "_").replace(".", "_")


def print_html(raw_data, n_bins):
    bin_counts, ranges = bin_data(raw_data, n_bins)
    methods = sorted(bin_counts.keys())

    print("<!DOCTYPE html>")
    print("<html>")
    print("<head>")
    print('<meta charset="utf-8">')
    print('<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>')
    print('<script type="text/javascript">')
    print("  google.charts.load('current', {packages: ['corechart', 'line']});")
    print("  google.charts.setOnLoadCallback(drawChart);")
    print("  function drawChart() {")

    for method in methods:
        prefix, word_size = method_prefix(method)
        meta = METHOD_META.get(prefix, {"title": method, "x_units": method})
        jid = js_safe(method)
        files = sorted(bin_counts[method].keys())
        lo = ranges[method]["min"]
        bs = ranges[method]["bin_size"]
        is_pct = prefix in PCT_METHODS

        print(f"    var data_{jid} = new google.visualization.DataTable();")
        x_units = meta["x_units"]
        print(f"    data_{jid}.addColumn('number', '{x_units}');")
        for f in files:
            print(f"    data_{jid}.addColumn('number', '{f}');")

        rows = []
        if bs == 0:
            # All values identical
            x = lo * 100 if is_pct else lo
            counts = ", ".join(str(bin_counts[method][f].get(0, 0)) for f in files)
            rows.append(f"[{x:.6g}, {counts}]")
        else:
            x = lo
            for i in range(n_bins + 1):
                xval = x * 100 if is_pct else x
                counts = ", ".join(str(bin_counts[method][f].get(i, 0)) for f in files)
                rows.append(f"[{xval:.6g}, {counts}]")
                x += bs

        print(f"    data_{jid}.addRows([")
        print("      " + ",\n      ".join(rows))
        print("    ]);")

        title = meta["title"]
        if word_size > 0:
            title += f" ( word size = {word_size} )"

        print(f"    var options_{jid} = {{")
        print(f"      chart: {{ title: '{title}' }},")
        print(f"      width: 900, height: 500,")
        print(f"      series: {{ 0: {{axis: 'bins'}} }},")
        print(f"      axes: {{ y: {{ bins: {{label: 'counts'}} }} }}")
        print(f"    }};")
        print(
            f"    var chart_{jid} = new google.charts.Line("
            f"document.getElementById('chart_{jid}'));"
        )
        print(f"    chart_{jid}.draw(data_{jid}, options_{jid});")

    print("  }")
    print("</script>")
    print("</head>")
    print("<body>")
    for method in methods:
        jid = js_safe(method)
        prefix, word_size = method_prefix(method)
        meta = METHOD_META.get(prefix, {"title": method, "x_units": method})
        title = meta["title"]
        if word_size > 0:
            title += f" (word size {word_size})"
        print(f'  <h3>{title}</h3>')
        print(f'  <div id="chart_{jid}"></div>')
    print("</body>")
    print("</html>")


def main():
    parser = argparse.ArgumentParser(
        description="Display complexity stats from rusted-garlic complexity TSV output."
    )
    parser.add_argument("-in", dest="infiles", required=True, nargs="+",
                        metavar="TSV",
                        help="One or more TSV files from `rusted-garlic complexity`; "
                             "each unique 'file' column value becomes a separate series")
    parser.add_argument("-graph", "-g", dest="graph", action="store_true",
                        help="Output Google Charts HTML instead of a summary table")
    parser.add_argument("-bins", dest="bins", type=int, default=100,
                        help="Number of histogram bins (default: 100, only used with -graph)")
    args = parser.parse_args()

    raw_data = read_tsv(args.infiles)
    if not raw_data:
        print(f"ERROR: No data found in input file(s)", file=sys.stderr)
        sys.exit(1)

    if args.graph:
        print_html(raw_data, args.bins)
    else:
        print_table(raw_data)


if __name__ == "__main__":
    main()
