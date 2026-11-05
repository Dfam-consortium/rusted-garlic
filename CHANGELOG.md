# Changelog

All notable changes to rusted-garlic. Validation details: `VALIDATION.md`.

## [Unreleased]

### Added

- `emit --seed <INT>`: all randomness comes from one seeded RNG, so a run can be
  reproduced exactly. Without `--seed` a random seed is drawn and logged.
- Regression tests (`tests/determinism.rs`, `cargo test --release`): build and seeded
  emit must be byte-identical across runs and match golden files in `tests/data/golden/`.
- Complexity test (`tests/complexity.rs`): `complexity` must match SeqComplex.pm on 199
  test windows; the Perl values are committed in `tests/data/complexity/`, so no Perl is
  needed.
- Unit tests for the emitter (`evolve.rs`): mutation operators and their acceptance
  rules, retry rounds, insertion placement (default vs `--perl_compat`), simple and TE
  repeat evolution, nested inserts, and TE/simple-repeat bombardment.

### Added (CpG)

- `emit --cpg_rate` (default **7**): in evolved repeats, a site in a CpG dinucleotide
  receives transitions 7× as often as other sites (CG → TG / CA), within the element's
  annotated divergence (RepeatMasker's `div` is raw mismatch, so it already includes
  CpG decay). Calibrated on human chr22: simulated repeats now match real TEs in CpG
  density (14.9 vs 14.8 CG/kb) and CG-containing 4-mers (0.0445 vs 0.0445); repeat
  4-mer JSD vs real TEs 0.049 → 0.041. `--cpg_rate 1` disables; ignored with
  `--perl_compat`.

### Changed

- Repeat bombardment gives up after 100 attempts in a row without a placement,
  counting as GARLIC 1.4 does (was: 101).
- **Repeat deletions** (`emit`) remove the base whose removal was scored; GARLIC 1.4
  scored the k-mer for deleting `seq[pos + k - 1]` but deleted `seq[pos]`. Deletions
  are now applied one at a time against the current sequence. On chr22 (Tier 2):
  CpG o/e in repeats 0.330 → 0.314 (real TEs 0.285), repeat 4-mer spectrum slightly
  closer to real TEs, background unchanged. `emit --perl_compat` keeps GARLIC 1.4
  deletions.
- `build --perl_compat`: GC binning and TE bin positions follow GARLIC 1.4 exactly
  (`calcBinGC`/`classGC`/`profileRMAlign`). On chr22 the GC bins, all 625 GC
  transition counts, and the TE, simple-repeat and nested-insert tables are then
  identical to GARLIC 1.4. Intended for validation; to be removed eventually.
- Smoothed window GC is now the mean of the neighbours' full-precision GC, rounded
  once (was: mean of values already rounded to 0.1%). As Perl.
- A chromosome's trailing partial window no longer counts as a window in the GC
  histogram, bin counts or transitions; it is still binned by its own GC.
- K-mers are binned by the window of their first base (was: last base), as Perl.
- Models built with earlier versions remain loadable; rebuild to pick up the changes.
- `build` output is now byte-identical across runs. Model maps are written in sorted
  key order and TE records in genomic order instead of hash order. Existing
  `model_version` 2 models stay valid; no rebuild needed.

### Fixed

- `complexity --min_kmer` > 1 panicked (`ct`/`cl` need the counts of every word size
  from 1); `min_kmer` now only selects output columns.
- `complexity` options are `--min_kmer`/`--max_kmer` as documented (`--min-kmer`/
  `--max-kmer` still accepted); duplicated defaults removed from `--help`.
- Validated against `SeqComplex.pm`: all measures identical to its per-measure
  subroutines except `cz`, where Perl's gzip header carries a temp-file name
  (`VALIDATION.md`).

## [0.2.0] — 2026-10-07

### ⚠ Model incompatibility — rebuild all models

Models built with 0.1.x **must be rebuilt** with `rusted-garlic build`.

- Models now carry `model_version` (currently **2**). `emit` refuses any other
  version, including 0.1.x models (no version field), with a "rebuild the model" error.
- 0.1.x models are also *wrong*, not just old: if built with `--te_annot`, their
  per-GC-bin k-mer tables are corrupted by the masking-drift bug below.
- New model fields: `gc_trans_counts` (raw GC transition counts) and `gc_kmer_totals`
  (per-bin k-mer totals).
- Simple-repeat records: the third element is now the TRF copy number
  (`copy_num`, float) instead of a truncated integer misnamed `period`.

### Fixed — model building (`build`)

- **k-mer table / GC-bin misassignment with masking.** FASTA newlines inside masked
  regions were counted as sequence positions, so the mask and GC-bin lookups drifted
  ahead (~0.9%, ≈450 kb by the end of chr22 with TE masking). Per-bin k-mer tables were
  smeared toward average GC; emitted window-GC sd was 4.6 instead of ~7.4.
- **Phantom window per chromosome.** The end-of-chromosome GC smoothing emitted one
  extra window, adding a spurious transition. chr22 transitions now equal GARLIC 1.4.
- **Simple-repeat profile matches GARLIC 1.4:** tandem repeats inside TEs are kept by
  default; overlapping repeats with equal unit length keep the later record.
  Simple-repeat records are now identical to GARLIC 1.4 (5,692 on chr22).

### Fixed — sequence generation (`emit`)

- **Repeat mutation follows GARLIC 1.4.** Rejected k-mer-guided mutations are retried
  for up to 10 rounds against a newly drawn GC bin, then dropped. Previously one guided
  round was followed by unguided forced mutations, leaving too much CpG in repeats.
- **Nested TE insertions** ported from GARLIC 1.4 `getInsert`: up to `nfrg` younger
  inserts, chosen from the age-sorted pool. Previously at most one, in one attempt,
  which made elements too short.
- **Simple repeats use the TRF copy number** (`int(copy_num) + 1` copies) instead of a
  random 3–20; no forced deletion when no indels are requested.
- **GC transitions** get a +1 pseudocount per cell; the GC bin for retry rounds is drawn
  by per-bin k-mer totals; `N` in TE consensi becomes a random base (all as GARLIC 1.4).

### Added

- `emit --perl_compat` (alias `--perl-compat`): reproduces three GARLIC 1.4 bugs for
  side-by-side comparison. By default rusted-garlic fixes them:
  - inserted bases in evolved repeats placed one base past the context they were
    sampled for (raises CpG in repeats);
  - TE records longer than their consensus cut to a random 3′ piece (negative `rand()`
    offset) instead of the full consensus;
  - background "1 kb" windows of 998 bp with a duplicated 3-bp seed, and the length
    shortfall padded with a copy of a random internal segment.
- `build --filter_te_tandems`: opt-in exclusion of tandem repeats contained in TEs
  (the 0.1.x default).
- Tier 2 validation (`../tier2/`): `run_chr22.sh` compares emitted sequences from
  GARLIC 1.4, `--perl_compat` and default. Equivalence is judged against run-to-run noise:
  Perl vs `--perl_compat` passes every check (e.g. 4-mer spectra ratio 0.99, where 1 = no
  more different than two Perl runs).

## [0.1.0]

Initial Rust port of GARLIC 1.4 (`createModel.pl` / `createFakeSequence.pl`).
