#!/bin/sh
# Wall-clock time and peak memory of GARLIC 1.4 vs rusted-garlic on human chr22:
# model build (createModel.pl vs build) and emission of one 5 Mb sequence with TE and
# simple-repeat insertion (createFakeSequence.pl vs emit, default and --perl_compat).
# Same inputs and repeat targets as tier2_emit.sh. All steps are single-threaded and
# run one after another. Needs GNU time (/usr/bin/time -v).
#
# Usage:  sh validation/speed.sh
# Time:   ~2.5 min
# Output: out/speed/summary.tsv (step, wall clock, max RSS in kB); out/speed/*.time
# Result (2026-11-05, 0.3.0): build Perl 54.3 s / 217 MB vs Rust 6.0 s / 42 MB;
#         emit 5 Mb Perl 52.6 s / 91 MB vs Rust 5.2 s / 44 MB (--perl_compat 6.1 s).
set -e
. "$(dirname "$0")/config.sh"
need_rust
need_perl_garlic
set_dataset chr22
S=$OUT/speed; rm -rf "$S"; mkdir -p "$S/perl_build" "$S/perl_emit"
T="/usr/bin/time -v"

$T "$RG" build --in_seq "$FA" --tandem_annot "$TRF" --te_annot "$TE_ANNOT" \
    --model "$S/rust.model" > "$S/rust_build.log" 2> "$S/rust_build.time"
(cd "$S/perl_build" && $T perl "$GARLIC_HOME/bin/createModel.pl" -m chr22_perl \
    -f "$FA" -r "$ALIGN" -t "$TRF" -g /dev/null > ../perl_build.log 2> ../perl_build.time)

$T "$RG" --log-level error emit -m "$S/rust.model" -o "$S/rust.fa" --out_bed "$S/rust.bed" \
    -s 5m --te_cons "$LIB" --te_frac $TE --sim_frac $SIM --seed 2 2> "$S/rust_emit.time"
$T "$RG" --log-level error emit -m "$S/rust.model" -o "$S/rust_compat.fa" \
    --out_bed "$S/rust_compat.bed" -s 5m --te_cons "$LIB" --te_frac $TE --sim_frac $SIM \
    --seed 1 --perl_compat 2> "$S/rust_compat_emit.time"
(cd "$S/perl_emit" && $T perl "$GARLIC_HOME/bin/createFakeSequence.pl" -m chr22_perl \
    -d "$S/perl_build/data" -s 5m -r $TE -l $SIM -o perl --tecons_file "$LIB" \
    > ../perl_emit.log 2> ../perl_emit.time)

printf 'step\twall\tmax_rss_kb\n' > "$S/summary.tsv"
for f in rust_build perl_build rust_emit rust_compat_emit perl_emit; do
    w=$(sed -n 's/.*Elapsed (wall clock) time.*: //p' "$S/$f.time")
    m=$(sed -n 's/.*Maximum resident set size (kbytes): //p' "$S/$f.time")
    printf '%s\t%s\t%s\n' $f "$w" "$m" >> "$S/summary.tsv"
done
cat "$S/summary.tsv"
