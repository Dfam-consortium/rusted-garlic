# Validation of rusted-garlic against GARLIC 1.4

`rusted-garlic` was checked against the original Perl GARLIC 1.4 (`createModel.pl`,
`createFakeSequence.pl`, `SeqComplex.pm`) at three levels: the model it builds (Tier 1),
the sequences it emits (Tier 2), and its own reproducibility over time (Tier 3).

The approach had two stages:

1. **Reproduce Perl first.** `build --perl_compat` and `emit --perl_compat` follow the Perl
   code path, including its bugs, and must be indistinguishable from GARLIC 1.4. This
   separates porting errors from intended changes.
2. **Then improve.** Every default that differs from Perl was measured against the real
   genome, and kept only if it moved the simulation closer to it. The list is in
   [Intentional differences](#intentional-differences-from-garlic-14).

Results below are from 2026-10-07 (version 0.3.0) unless noted.

## Datasets

| Dataset | Size | Annotation | Role |
|---|---|---|---|
| Human chr22 (hg38) | 50.8 Mb | RepeatMasker 4.2.4 + Dfam 4.0; UCSC TRF | main validation set |
| `sm-example` (GARLIC's example input) | 47.8 kb, 3 sequences | RepeatMasker; TRF | out-of-sample check; test fixtures |

TE consensi: Dfam 4.0 *Homo sapiens* library (FamDB export converted to FASTA; GARLIC's
EMBL parser keys records by `ID`, which is the accession in Dfam 4.0). 13 of 1,024 chr22
model names (mostly `*_short_` Alus) are not in Dfam 4.0; both tools skip them.
Repeat targets for emission are the real fractions of non-N bases: chr22 `--te_frac 44.3
--sim_frac 1.7`; sm-example `38.8 / 3.6`.

## Tier 1: model equivalence

The Perl model files and the Rust JSON model built from the same chr22 inputs were
compared table by table (GC bin boundaries, GC transition matrix, per-bin conditional
k-mer frequencies, TE / simple-repeat records per bin, nested-insert table).

| Table | `build --perl_compat` vs Perl | default `build` vs Perl |
|---|---|---|
| GC bins (25) | identical | 1 boundary differs (lowest bin) |
| GC transitions (625 cells, 39,395 counts) | identical | max cell diff 0.038 |
| Conditional k-mer frequencies (6,400) | max diff 0.0052 | max diff 0.0095 |
| TE / simple-repeat records (56,790 / 5,692) | identical per bin | same totals; 4 mid-GC bins shift |
| TE divergence / ins / del means | identical | identical |
| Nested-insert pairs (12,018) | identical | 98.8% shared |

The residual k-mer difference in compat mode is a GARLIC 1.4 masking bug (the last 2 bases
of each TE and the last base of each TRF interval are left unmasked). Emulating it brings
the max diff to 0.0019, which confirms the cause; it is deliberately not reproduced.

## Tier 2: emitted sequences

Emission is random, so two runs never produce the same sequence. The test asks whether
Perl and Rust differ **more than two runs of the same program do**. Each program emits
5 replicates (chr22: 5 × 5 Mb; sm-example: 5 × 1 Mb) and three groups are compared:
`perl`, `rust_compat` (`build`/`emit --perl_compat`) and `rust` (defaults).

- **Distributions** — window-GC histogram, 4-mer spectra per GC bin (all bases, and
  background only): ratio = mean Jensen–Shannon divergence between Perl and Rust runs ÷
  mean JSD between runs of the same program. ≈ 1 means indistinguishable. Pass: < 1.5.
- **Scalars** — repeat coverage, insert counts and lengths, GC mean and sd, GC
  autocorrelation (lag 1 and 10 windows), CpG o/e overall and inside repeats:
  z = difference of means / standard error over replicates. Pass: z < 3.
- **Complexity** — distributions of the `SeqComplex` measures per 1 kb window,
  two-sample KS statistic.

**Pass criterion:** every Perl vs `rust_compat` comparison passes. With about 36 tests per
run, a single z of 3–4 can be chance; a real difference repeats across runs.

| | chr22: Perl vs `rust_compat` | chr22: Perl vs `rust` | sm-example: Perl vs `rust_compat` | sm-example: Perl vs `rust` |
|---|---|---|---|---|
| 4-mer spectra, all bases (ratio) | 1.06 | **2.97** | 1.01 | **2.57** |
| 4-mer spectra, background (ratio) | 1.06 | 1.20 | 1.04 | **3.31** |
| window-GC histogram (ratio) | 1.01 | **1.75** | 1.01 | **2.30** |
| CpG o/e overall / in repeats (z) | 1.8 / 1.1 | **40 / 63** | 1.3 / 1.9 | 1.5 / **32** |
| GC mean (z) | 0.2 | **3.2** | **3.0** | 0.2 |
| GC sd / lag-1 autocorrelation (z) | 2.5 / 1.7 | 2.4 / 1.1 | 0.5 / 2.2 | **10.3 / 7.7** |
| repeat coverage (z) | 0.9 | 0.7 | **3.3** | **3.7** |
| complexity distributions (max KS) | 0.026 | — | 0.037 | — |

**Compat mode matches Perl** on chr22 on every metric, and on sm-example apart from two
borderline z values (see open issue below). **The default differs from Perl by design**,
and on chr22 every difference moves toward the real chromosome:

| chr22 | real | Perl | `rust_compat` | `rust` (default) |
|---|---|---|---|---|
| Mean GC % | 47.0 | 49.9 | — | 48.7 |
| Per-window 4-mer JSD vs real | 0 | 0.0089 | — | 0.0082 |
| Repeats: CpG per kb | 14.81 | 24.06 | 23.95 | **14.94** |
| Repeats: CG-containing 4-mers | 0.0445 | 0.0711 | 0.0708 | **0.0445** |
| Repeats: 4-mer JSD vs real TEs | 0 | 0.0488 | 0.0504 | **0.0412** |

The CpG improvement was calibrated on chr22 and checked out of sample on sm-example:
repeat 4-mer JSD vs real TEs 0.020 (Rust) vs 0.027 (Perl).

On sm-example (47 windows) the GC transition table holds ~45 counts against ~600
pseudocounts, so the model has almost no GC memory. Perl's lag-1 autocorrelation there
(0.25) is an artifact of its 998-bp windows overlapping the 1 kb measurement windows;
Rust's exact windows show the model's real value (0.08).

### `complexity` vs `SeqComplex.pm`

199 single-window records (150 chr22 windows including 10 with N, 40 simulated, poly-A,
di/tri/tetranucleotide repeats, 25% N, and records of 50, 7 and 1 bp) were run through
both. Every measure except `cz` is **identical** to the `SeqComplex.pm` per-measure
subroutines on every record (the subroutines fail on the 1-bp record). `cz` differs because Perl gzips a temp file whose random
name is stored in the gzip header (~13 bytes, ~3.6% low on 1 kb windows); Rust
compresses in memory without a name. Perl's inline path (`gatherStats.pl`) disagrees with
its own subroutines in several places; Rust follows the subroutines.

The windows and SeqComplex.pm's values for them are committed
(`tests/data/complexity/`), and `tests/complexity.rs` repeats this check on every
`cargo test`.

## Tier 3: regression tests

`cargo test --release` (run on Linux and macOS by the release workflow):

- 36 unit tests, mostly of the emitter: mutation operators and their acceptance rules,
  retry rounds, insertion placement (default and `--perl_compat`), simple and TE repeat
  evolution, nested inserts, bombardment.
- 7 golden-file tests (`tests/determinism.rs`): `build` and `build --perl_compat` on
  sm-example, and `emit --seed 1` (default and `--perl_compat`), must be byte-identical
  across runs and match `tests/data/golden/`. Fixtures are described in
  `tests/data/README.md`.
- 1 complexity test (`tests/complexity.rs`): `complexity` on the 199 test windows must
  match the committed SeqComplex.pm values to 1e-6 for every measure except `cz`.

The emitter tests were checked by fault injection: 7 behaviour-changing faults were
planted one at a time in `evolve.rs` (mutation operators, acceptance rules, insertion
placement), and each made at least one test fail.

## Intentional differences from GARLIC 1.4

All are reverted by `--perl_compat`. They are explained in the README's
[Differences from the original Perl GARLIC](README.md#differences-from-the-original-perl-garlic).

| Area | GARLIC 1.4 | rusted-garlic default | Evidence |
|---|---|---|---|
| GC bin cut rule | one pass, fixed bucket size; uneven bins on narrow GC distributions | equal-population bins | identical bins on chr22 and typical genomes at 1 kb windows |
| Repeat masking | last 2 bases of each TE / last base of each TRF repeat left unmasked | full interval masked | Tier 1 |
| Deletions in evolved repeats | scores deleting `seq[pos+k-1]`, deletes `seq[pos]` | deletes the scored base | Tier 2 |
| Insertions in evolved repeats | base sampled for the context before `pos` is placed after it | placed where sampled | CpG o/e in repeats 0.357 → 0.333 |
| Truncated full-length TEs | for ~2.6% of records (mostly L1, LTR) a random 3′ piece replaces the full consensus | full consensus used | ~2% fewer, longer TEs |
| Window length | each "1 kb" window is 998 bp; shortfall padded with a copy of a random segment | exact windows | sm-example autocorrelation |
| CpG decay (`--cpg_rate`, default 7) | all repeat sites mutate at the same rate | CpG sites get 7× the transition rate | CpG/kb 14.94 vs real 14.81 |

Not changed: ties in the k-mer acceptance test follow Perl's rule; a stricter rule was
measured and fit real repeat 4-mers worse.

## Known gaps

- Simulated repeats are 4–6 GC points richer than real TEs in both programs (chr22: real
  45.8%, Perl 52.1%, Rust 50.1%), so CpG o/e in repeats undershoots even though CpG
  density matches.
- On sm-example, Perl stops inserting TEs early more often (31.8–37.6% per 1 Mb replicate
  vs Rust 35.6–39.5%, target 38.8%). Model repeat tables are identical and both skip the
  same consensi; the cause has not been found. Not seen on chr22.
- Overall GC mean, whole-sequence 4-mer JSD and lag-10 GC autocorrelation are equally far
  from real chr22 in both programs: they are limits of GARLIC's model, not of the port.

## Performance

Human chr22 (50.8 Mb), same inputs and repeat targets as Tier 2; one run of each, both
single-threaded, on a shared Linux server (2026-11-05, version 0.3.0):

| Step | GARLIC 1.4 (Perl) | rusted-garlic | Speed-up |
|---|---|---|---|
| Build model (`createModel.pl` / `build`) | 54.3 s, 217 MB | 6.0 s, 42 MB | 9.1× |
| Emit one 5 Mb sequence with TE and simple-repeat bombardment (`createFakeSequence.pl` / `emit`) | 52.6 s, 91 MB | 5.2 s, 44 MB | 10.0× |
| Same, `emit --perl_compat` | | 6.1 s, 44 MB | 8.6× |

Times are wall clock; memory is peak resident set size. The server was fully loaded
during the runs, which slows both programs alike.

## Reproducing

The scripts are in [`validation/`](validation/); its README lists the requirements and
how to read each report. sm-example runs from files in this repository; chr22 has to be
downloaded and annotated with RepeatMasker first (instructions there). Perl GARLIC 1.4
is needed for Tiers 1 and 2 and the performance comparison; the complexity check and Tier 3
need nothing beyond this repository.

| Check | Command (from the repository root) |
|---|---|
| Tier 1 | `sh validation/tier1_models.sh chr22` (or `sm-example`) |
| Tier 2 | `sh validation/tier2_emit.sh chr22` (or `sm-example`) |
| Complexity | `sh validation/complexity.sh` |
| Tier 3 | `cargo test --release` |
| Performance | `sh validation/speed.sh` |
