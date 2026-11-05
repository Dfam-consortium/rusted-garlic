#!/bin/sh
# Tier 2: compare sequences emitted by GARLIC 1.4 and rusted-garlic.
# Emits N replicates (chr22: 5 Mb each; sm-example: 1 Mb each) from three pipelines:
#   perl         createModel.pl + createFakeSequence.pl
#   rust_compat  build --perl_compat + emit --perl_compat (should match Perl)
#   rust         default build + emit (intended differences)
# and compares them with each other and with the real sequence.
#
# Usage:  sh validation/tier2_emit.sh chr22|sm-example
# Time:   chr22 ~2 min with 5 cores free (Perl replicates run in parallel); sm-example ~15 s
# Output: out/<dataset>/results/
#           report.tsv   main report; pass = every perl_vs_rust_compat row in the
#                        "## Equivalence" section is OK
#           inserts.tsv  inserted-element counts and lengths by class (Perl vs rust_compat)
#           repeats.tsv  simulated repeats vs real TEs (GC, CpG, 4-mers)
# Rust emits use fixed seeds and are reproducible; Perl's are not.
set -e
. "$(dirname "$0")/config.sh"
set_dataset "${1:?usage: $0 chr22|sm-example}"

perl_model
rust_models
E=$W/emit

echo "emitting rust_compat and rust ($N x $SIZE)"
"$RG" --log-level error emit -m "$W/models/rust_compat.model" -o "$E/rust_compat.fa" \
    --out_bed "$E/rust_compat.bed" -s $SIZE -N $N --te_cons "$LIB" \
    --te_frac $TE --sim_frac $SIM --perl_compat --seed $SEED_COMPAT
"$RG" --log-level error emit -m "$W/models/rust.model" -o "$E/rust.fa" \
    --out_bed "$E/rust.bed" -s $SIZE -N $N --te_cons "$LIB" \
    --te_frac $TE --sim_frac $SIM --seed $SEED_DEFAULT

# createFakeSequence.pl writes one sequence per run; -d is the directory that
# contains the model directory.
echo "emitting perl ($N x $SIZE, in parallel)"
rm -f "$E"/perl_*.fasta "$E"/perl_*.inserts
i=1
while [ $i -le $N ]; do
    (cd "$E" && perl "$GARLIC_HOME/bin/createFakeSequence.pl" -m "$PERL_MODEL" \
        -d "$PERL_MODEL_DIR" -s $SIZE -r $TE -l $SIM -o perl_$i \
        --tecons_file "$LIB" > perl_$i.log 2>&1) &
    i=$((i + 1))
done
wait
ls "$E"/perl_*.fasta > /dev/null 2>&1 || die "createFakeSequence.pl produced no output, see $E/perl_*.log"

echo "comparing"
REFM=
[ -n "$MASKED" ] && REFM="--ref_masked $MASKED"
python3 "$SCRIPTS/compare_emitted.py" --ref "$FA" $REFM \
    --group perl "$E"/perl_*.fasta --group rust_compat "$E/rust_compat.fa" \
    --group rust "$E/rust.fa" --complexity "$RG" --target_te $TE > "$W/results/report.tsv"
python3 "$SCRIPTS/compare_inserts.py" --perl "$E"/perl_*.inserts \
    --rust "$E/rust_compat.bed" --mb_per_seq $MB_PER_SEQ > "$W/results/inserts.tsv"
python3 "$SCRIPTS/compare_repeats.py" "$FA" "$TE_ANNOT" \
    perl="$E/perl_*.fasta" rust_compat="$E/rust_compat.fa" rust="$E/rust.fa" \
    > "$W/results/repeats.tsv"

echo
sed -n '/^## Equivalence/,/^## /p' "$W/results/report.tsv" | grep -E '^## Equivalence|perl_vs_rust_compat' || true
echo
echo "full results in $W/results/"
