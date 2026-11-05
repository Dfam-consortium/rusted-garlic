# Test fixtures

Used by `tests/determinism.rs` (Tier 3 regression tests).

| File | Contents |
|---|---|
| `sm-example.fa` | 50 kb, 3 sequences (`ex1`–`ex3`); the GARLIC `sm-example` input |
| `sm-example.te_annot.tsv` | RepeatMasker TE annotation of it (`utils/rm_align_to_tsv.py` format) |
| `sm-example.fa.trf.bed` | TRF tandem-repeat annotation |
| `sm-example.cons.fa` | The 38 Dfam 4.0 consensus sequences named in the TE annotation (CC0). `AluJb_short_` and `FLAM_C_short_` are not in Dfam 4.0 and are skipped, as in GARLIC 1.4. |
| `golden/sm-example.model` | Expected `build` output |
| `golden/sm-example.perl_compat.model` | Expected `build --perl_compat` output |
| `golden/emit_seed1.{fa,bed}` | Expected `emit --seed 1` (2 × 20 kb, `--te_frac 40 --sim_frac 2`) |
| `golden/emit_perl_compat_seed1.{fa,bed}` | Same with `--perl_compat` |

The emit goldens are generated from `golden/sm-example.model`, so a build change and an
emit change fail separate tests.

Goldens depend on the exact random streams of `rand` (`StdRng`, pinned in `Cargo.lock`).
Upgrading `rand` will change them; that is expected. After any intended output change:

    GARLIC_UPDATE_GOLDEN=1 cargo test --release --test determinism build_
    GARLIC_UPDATE_GOLDEN=1 cargo test --release --test determinism golden

(the model golden first, because the emit goldens are generated from it). Then check
`git diff --stat tests/data/golden` and commit.

## `complexity/` (used by `tests/complexity.rs`)

| File | Contents |
|---|---|
| `complexity/windows.fa` | 199 one-window records: 150 hg38 chr22 windows (10 with N), 40 windows of a rusted-garlic chr22 simulation, poly-A, di/tri/tetranucleotide repeats, 25% N, and 50, 7 and 1 bp records. Made by `validation/scripts/make_windows.py`. |
| `complexity/seqcomplex_reference.tsv` | SeqComplex.pm's per-measure subroutines on each record (`NA` where Perl cannot compute it). |

To regenerate the reference (needs SeqComplex.pm): `SEQCOMPLEX_HOME=<dir>
sh validation/complexity.sh --regenerate`. `cz` and some 15th-digit values change from
run to run on the Perl side; the test ignores `cz` and compares to 1e-6.
