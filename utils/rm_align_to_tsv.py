#!/usr/bin/env python3
"""
rm_align_to_tsv.py — Convert a RepeatMasker .align (or .out) file to the
rusted-garlic repeat annotation TSV format.

Usage:
    python3 rm_align_to_tsv.py genome.fa.align > genome.repeats.tsv
    python3 rm_align_to_tsv.py genome.fa.align.gz > genome.repeats.tsv
    zcat genome.fa.align.gz | python3 rm_align_to_tsv.py - > genome.repeats.tsv

Output columns (tab-separated, 0-based half-open coordinates):
    chrom  start  end  name  family  strand  div  ins  del  cons_start  cons_end  rep_id

Simple_repeat, Low_complexity, Unknown, and Satellite entries are excluded.
"""

import gzip
import sys


# Classes excluded from the interspersed-repeat model (handled by --trf instead).
SKIP_CLASSES = {"Simple_repeat", "Low_complexity", "Unknown", "Satellite"}

HEADER = (
    "#chrom\tstart\tend\tname\tfamily\tstrand"
    "\tdiv\tins\tdel\tcons_start\tcons_end\trep_id"
)


def open_input(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def is_data_line(s):
    """True if the stripped line looks like a RM alignment header row
    (starts with integer score then decimal divergence)."""
    parts = s.split(None, 2)
    if len(parts) < 2:
        return False
    try:
        int(parts[0])
        float(parts[1])
        return True
    except ValueError:
        return False


def split_name_fam(token):
    """Split 'Name#Class/Family' into (name, family)."""
    if "#" in token:
        name, fam = token.split("#", 1)
    else:
        name, fam = token, "Unknown"
    return name, fam


def convert(path):
    print(HEADER)
    with open_input(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            stripped = line.lstrip()

            if not is_data_line(stripped):
                continue

            # Skip simple / low-complexity / satellite classes.
            if any(cls in stripped for cls in SKIP_CLASSES):
                continue

            f = stripped.split()

            # RepeatMasker .align header:
            #   Positive strand (13 fields):
            #     score div del ins chrom begin end (left) name#fam rbegin rend (rleft) id
            #     [0]   [1] [2] [3]  [4]  [5]  [6]  [7]    [8]     [9]  [10]   [11]  [12]
            #
            #   Negative strand (14 fields) — 'C' inserted at [8]:
            #     score div del ins chrom begin end (left) C name#fam (rleft) rend rbegin id
            #     [0]   [1] [2] [3]  [4]  [5]  [6]  [7]  [8]  [9]    [10]  [11]  [12]  [13]

            if len(f) < 13:
                continue

            div, ins, del_ = f[1], f[2], f[3]
            chrom = f[4]
            # Convert 1-based closed [start, end] → 0-based half-open.
            start = int(f[5]) - 1
            end   = int(f[6])
            rep_id = f[-1]

            if f[8] == "C":
                if len(f) < 14:
                    continue
                name, family = split_name_fam(f[9])
                strand = "-"
                # Negative strand: rbegin is at [12], rend at [11]
                # (consensus coords are always leftmost..rightmost).
                cons_start = f[12]
                cons_end   = f[11]
            else:
                name, family = split_name_fam(f[8])
                strand = "+"
                cons_start = f[9]
                cons_end   = f[10]

            print(
                f"{chrom}\t{start}\t{end}\t{name}\t{family}\t{strand}"
                f"\t{div}\t{ins}\t{del_}\t{cons_start}\t{cons_end}\t{rep_id}"
            )


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    convert(sys.argv[1])


if __name__ == "__main__":
    main()
