# Shared settings for the validation scripts; sourced, not run.
# Every setting can be overridden from the environment.

VAL=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$VAL/.." && pwd)
SCRIPTS=$VAL/scripts

# rusted-garlic release build
RG=${RG:-$REPO/target/release/rusted-garlic}
# GARLIC 1.4 (Perl) checkout: must contain bin/createModel.pl and bin/createFakeSequence.pl
GARLIC_HOME=${GARLIC_HOME:-}
# Directory containing SeqComplex.pm (complexity check only)
SEQCOMPLEX_HOME=${SEQCOMPLEX_HOME:-}
# Human chr22 inputs you prepare yourself (see README.md)
CHR22_DIR=${CHR22_DIR:-$VAL/data/chr22}
# Dfam 4.0 Homo sapiens TE consensi as FASTA (see README.md)
TE_LIB=${TE_LIB:-$VAL/data/te_lib/dfam-4.0-homo-sapiens.fa}
# All outputs go here, one subdirectory per dataset
OUT=${OUT:-$VAL/out}
# Number of emitted replicates per program
N=${N:-5}

die() { echo "error: $*" >&2; exit 1; }

need_rust() {
    [ -x "$RG" ] || die "no rusted-garlic binary at $RG (run 'cargo build --release' in $REPO, or set RG)"
}

need_perl_garlic() {
    [ -n "$GARLIC_HOME" ] || die "set GARLIC_HOME to your GARLIC 1.4 checkout"
    [ -f "$GARLIC_HOME/bin/createModel.pl" ] || die "no bin/createModel.pl under GARLIC_HOME=$GARLIC_HOME"
}

need_file() {
    [ -s "$1" ] || die "missing $1${2:+ ($2)}"
}

# set_dataset chr22|sm-example: input files, repeat targets and emit settings.
# TE / SIM are % of non-N bases covered by interspersed and simple/low-complexity
# repeats in the real sequence (chr22: RepeatMasker .tbl; sm-example: union of
# annotated intervals).
set_dataset() {
    DS=$1
    case $DS in
    chr22)
        FA=$CHR22_DIR/chr22.fa
        ALIGN=$CHR22_DIR/chr22.fa.align
        MASKED=$CHR22_DIR/chr22.fa.masked
        TRF=$CHR22_DIR/chr22.fa.trf.bed
        TE_ANNOT=$CHR22_DIR/chr22.fa.te_annot.tsv
        LIB=$TE_LIB
        TE=${TE:-44.3}; SIM=${SIM:-1.7}; SIZE=${SIZE:-5m}; MB_PER_SEQ=5
        PERL_MODEL_OPTS=""
        SEED_COMPAT=1; SEED_DEFAULT=2
        need_file "$FA" "hg38 chr22, see README.md"
        need_file "$ALIGN" "RepeatMasker -a output, see README.md"
        need_file "$MASKED" "RepeatMasker .masked output"
        need_file "$TRF" "UCSC simpleRepeat for chr22, see README.md"
        need_file "$LIB" "Dfam 4.0 human library as FASTA, see README.md"
        if [ ! -s "$TE_ANNOT" ]; then
            echo "converting $ALIGN -> $TE_ANNOT"
            python3 "$REPO/utils/rm_align_to_tsv.py" "$ALIGN" > "$TE_ANNOT"
        fi
        ;;
    sm-example)
        D=$REPO/tests/data
        FA=$D/sm-example.fa
        ALIGN=$VAL/data/sm-example.fa.align
        MASKED=
        TRF=$D/sm-example.fa.trf.bed
        TE_ANNOT=$D/sm-example.te_annot.tsv
        # The 38 consensi named in the annotation are enough; override with the full
        # Dfam library to reproduce the published run exactly.
        LIB=${SM_TE_LIB:-$D/sm-example.cons.fa}
        TE=${TE:-38.8}; SIM=${SIM:-3.6}; SIZE=${SIZE:-1m}; MB_PER_SEQ=1
        PERL_MODEL_OPTS="--no_mask_gene"
        SEED_COMPAT=11; SEED_DEFAULT=12
        ;;
    *) die "unknown dataset '$DS' (chr22 or sm-example)" ;;
    esac
    W=$OUT/$DS
    PERL_MODEL_DIR=$W/perl/data          # createModel.pl writes data/<model>/ under its cwd
    PERL_MODEL=${DS}_perl
    mkdir -p "$W/models" "$W/emit" "$W/results" "$W/perl"
}

# Build the GARLIC 1.4 model for the current dataset, unless it already exists.
perl_model() {
    need_perl_garlic
    if [ -s "$PERL_MODEL_DIR/$PERL_MODEL/$PERL_MODEL.GCt.W1000.data" ]; then
        echo "Perl model exists: $PERL_MODEL_DIR/$PERL_MODEL"
        return
    fi
    echo "building Perl model ($PERL_MODEL)"
    (cd "$W/perl" && perl "$GARLIC_HOME/bin/createModel.pl" -m "$PERL_MODEL" \
        -f "$FA" -r "$ALIGN" -t "$TRF" -g /dev/null $PERL_MODEL_OPTS > build.log 2>&1) \
        || die "createModel.pl failed, see $W/perl/build.log"
}

# Build the two rusted-garlic models: default and --perl_compat.
rust_models() {
    need_rust
    echo "building rusted-garlic models"
    "$RG" build --in_seq "$FA" --tandem_annot "$TRF" --te_annot "$TE_ANNOT" \
        --model "$W/models/rust.model" > "$W/models/rust.build.log" 2>&1
    "$RG" build --in_seq "$FA" --tandem_annot "$TRF" --te_annot "$TE_ANNOT" \
        --perl_compat --model "$W/models/rust_compat.model" > "$W/models/rust_compat.build.log" 2>&1
}
