#!/bin/sh
# Check `rusted-garlic complexity` against SeqComplex.pm, the Perl module it replaces.
#
# The test windows and SeqComplex.pm's values for them are committed in
# tests/data/complexity/, so the check itself needs no Perl. It also runs as part of
# `cargo test` (tests/complexity.rs); this script prints a per-measure report.
#
# Usage:  sh validation/complexity.sh               compare against the reference
#         sh validation/complexity.sh --regenerate  rerun SeqComplex.pm on windows.fa and
#                                                   rewrite the reference (needs SEQCOMPLEX_HOME)
# Output: out/complexity/compare.txt
# Expect: column 2 (max |Rust - Perl|) is 0 for every measure except cz. cz differs
#         because Perl's gzip header stores a random temp-file name. NA:1 is the 1-bp
#         record, which SeqComplex.pm cannot compute.
#
# windows.fa was made with
#   python3 validation/scripts/make_windows.py chr22.fa simulated.fa > windows.fa
# (150 hg38 chr22 windows, 40 windows of a rusted-garlic chr22 simulation, extremes).
set -e
. "$(dirname "$0")/config.sh"
need_rust
REF=$REPO/tests/data/complexity
C=$OUT/complexity; mkdir -p "$C"

if [ "$1" = "--regenerate" ]; then
    [ -n "$SEQCOMPLEX_HOME" ] || die "set SEQCOMPLEX_HOME to the directory containing SeqComplex.pm"
    # perl_complexity.pl reports both Perl code paths: s_<m> (per-measure subroutines,
    # the reference) and i_<m> (runAllMethodsInline, which disagrees with them).
    SEQCOMPLEX_HOME=$SEQCOMPLEX_HOME perl "$SCRIPTS/perl_complexity.pl" "$REF/windows.fa" \
        > "$C/perl_both_paths.tsv" 2> "$C/perl.log"   # SeqComplex.pm warns on short / N-rich records
    python3 - "$C/perl_both_paths.tsv" "$REF/seqcomplex_reference.tsv" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
cols = [c for c in rows[0] if c.startswith("s_")]
with open(sys.argv[2], "w") as f:
    f.write("\t".join(["seq"] + [c[2:] for c in cols]) + "\n")
    for r in rows:
        f.write("\t".join([r["seq"]] + [r[c] for c in cols]) + "\n")
EOF
    echo "rewrote $REF/seqcomplex_reference.tsv"
fi

"$RG" --log-level error complexity "$REF/windows.fa" > "$C/rust.tsv"   # one window per record
python3 "$SCRIPTS/compare_complexity.py" "$C/rust.tsv" "$REF/seqcomplex_reference.tsv" | tee "$C/compare.txt"
