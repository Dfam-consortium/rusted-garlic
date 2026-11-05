#!/bin/sh
# Tier 1: compare the model GARLIC 1.4 builds with the models rusted-garlic builds
# from the same inputs (GC bins, GC transitions, k-mer tables, repeat records,
# nested-insert table).
#
# Usage:  sh validation/tier1_models.sh chr22|sm-example
# Output: out/<dataset>/results/tier1_compat.txt   (build --perl_compat vs Perl)
#         out/<dataset>/results/tier1_default.txt  (default build vs Perl)
# Expect: compat - bins, transitions, repeat records and inserts identical; k-mer
#         max diff ~0.005 (GARLIC 1.4's repeat-end masking bug, not reproduced).
set -e
. "$(dirname "$0")/config.sh"
set_dataset "${1:?usage: $0 chr22|sm-example}"

perl_model
rust_models

P=$PERL_MODEL_DIR/$PERL_MODEL/$PERL_MODEL
for m in compat default; do
    r=$W/models/rust_$m; [ $m = default ] && r=$W/models/rust
    python3 "$SCRIPTS/compare_models.py" "$P" "$r" > "$W/results/tier1_$m.txt"
    echo "wrote $W/results/tier1_$m.txt"
done
