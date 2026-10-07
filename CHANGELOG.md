# Changelog

All notable changes to rusted-garlic. Validation details: `../tier2/README.md` and
`../PERL_VS_RUST_DIFFS.md`.

## [Unreleased]

### Added

- `emit --seed <INT>`: all randomness comes from one seeded RNG, so a run can be
  reproduced exactly. Without `--seed` a random seed is drawn and logged.
- Regression tests (`tests/determinism.rs`, `cargo test --release`): build and seeded
  emit must be byte-identical across runs and match golden files in `tests/data/golden/`.

### Changed

- `build` output is now byte-identical across runs. Model maps are written in sorted
  key order and TE records in genomic order instead of hash order. Existing
  `model_version` 2 models stay valid; no rebuild needed.

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
