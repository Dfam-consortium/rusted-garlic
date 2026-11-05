# Validation scripts

Scripts that compare `rusted-garlic` with the original Perl GARLIC 1.4. The method and
the published results are in [`../VALIDATION.md`](../VALIDATION.md); this directory lets
you rerun them.

| Script | What it checks | Needs |
|---|---|---|
| `tier1_models.sh <dataset>` | model tables: GC bins, GC transitions, k-mers, repeat records, nested inserts | GARLIC 1.4 |
| `tier2_emit.sh <dataset>` | statistics of emitted sequences: Perl vs `--perl_compat` vs default vs the real sequence | GARLIC 1.4 |
| `complexity.sh` | `rusted-garlic complexity` vs `SeqComplex.pm` | nothing (reference values are committed); `--regenerate` needs SeqComplex |
| `speed.sh` | wall clock and peak memory, Perl vs Rust | GARLIC 1.4, chr22, GNU `time` |

`<dataset>` is `sm-example` or `chr22`. Tier 3 (regression tests) is `cargo test --release`,
which also runs the complexity check (`tests/complexity.rs`).

**sm-example needs no data downloads:** its sequence and annotations are in
`../tests/data/` and its RepeatMasker alignment is in `data/`. Start there to check your
setup. chr22 is the main validation set, but you need to fetch and annotate it yourself
(below).

## Requirements

- A release build: `cargo build --release` in the repository root (or set `RG`).
- Python 3 with numpy.
- Perl 5 (core modules only).
- GARLIC 1.4 (Perl), for every check except `complexity.sh`. Set `GARLIC_HOME` to the
  directory that contains `bin/createModel.pl` and `bin/createFakeSequence.pl`.
- `SeqComplex.pm`, only to regenerate the complexity reference
  (`complexity.sh --regenerate`). Set `SEQCOMPLEX_HOME` to the directory that contains it.

All settings are environment variables with defaults in `config.sh`:

| Variable | Default | Meaning |
|---|---|---|
| `RG` | `../target/release/rusted-garlic` | rusted-garlic binary |
| `GARLIC_HOME` | — | GARLIC 1.4 checkout |
| `SEQCOMPLEX_HOME` | — | directory with `SeqComplex.pm` |
| `CHR22_DIR` | `data/chr22` | chr22 inputs (below) |
| `TE_LIB` | `data/te_lib/dfam-4.0-homo-sapiens.fa` | TE consensi for chr22 |
| `SM_TE_LIB` | `../tests/data/sm-example.cons.fa` | TE consensi for sm-example |
| `OUT` | `out` | outputs, one subdirectory per dataset |
| `N` | 5 | replicates per program in Tier 2 |

## Preparing chr22

The published results used hg38 chr22, RepeatMasker 4.2.4 with Dfam 4.0 (rmblastn
2.17.1), and UCSC's TRF track. Put these files in `data/chr22/` (or set `CHR22_DIR`):

| File | How to make it |
|---|---|
| `chr22.fa` | UCSC hg38 chr22, upper-cased (see below) |
| `chr22.fa.align`, `chr22.fa.masked` | `RepeatMasker -species human -a chr22.fa` |
| `chr22.fa.trf.bed` | UCSC hg38 TRF track for chr22 (see below) |
| `chr22.fa.te_annot.tsv` | made from `chr22.fa.align` automatically on first use (`../utils/rm_align_to_tsv.py`) |

```sh
mkdir -p data/chr22 && cd data/chr22

# Sequence. UCSC's file is soft-masked; GARLIC expects upper case.
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip -c chr22.fa.gz | tr acgtn ACGTN > chr22.fa

# RepeatMasker. The published run was `RepeatMasker -pa 32 -a chr22.fa`, whose
# .tbl reports the species as homo sapiens; -species human makes that explicit.
# Writes chr22.fa.align and chr22.fa.masked.
RepeatMasker -pa 16 -species human -a chr22.fa

# Tandem repeats: the bigZips TRF file (period <= 12), which is what GARLIC
# downloads itself, not the larger simpleRepeat table.
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz
gunzip -c hg38.trf.bed.gz | awk -F'\t' '$1=="chr22"' | sort -k2,2n -k3,3n > chr22.fa.trf.bed
```

These steps reproduce the published inputs: the sequence is identical, and the TRF file
has the same 6,440 records (line order may differ slightly). A different RepeatMasker or
Dfam version gives slightly different annotation and therefore slightly different
numbers, but the Perl-vs-Rust comparison is still valid because both programs read the
same files.

### TE library

Both programs need the TE consensus sequences as FASTA keyed by family name. The
published runs used the curated Dfam 4.0 *Homo sapiens* families (with ancestors):

```sh
mkdir -p data/te_lib
famdb.py -i <Dfam 4.0 FamDB directory> families -ad --curated --format embl 'homo sapiens' \
    > data/te_lib/dfam-4.0-homo-sapiens.embl
python3 scripts/embl2fa.py data/te_lib/dfam-4.0-homo-sapiens.embl \
    > data/te_lib/dfam-4.0-homo-sapiens.fa
```

The EMBL file can't be used directly: GARLIC 1.4's EMBL parser keys records by the `ID`
field, which is the accession in Dfam 4.0, not the family name. 13 of the 1,024 family
names in the chr22 annotation (mostly `*_short_` Alus) are not in Dfam 4.0; both programs
skip them.

## Running

From the repository root:

```sh
export GARLIC_HOME=/path/to/Garlic

sh validation/tier1_models.sh sm-example     # seconds
sh validation/tier2_emit.sh  sm-example      # ~15 s
sh validation/tier1_models.sh chr22          # ~1 min (Perl model build)
sh validation/tier2_emit.sh  chr22           # ~2 min with 5 free cores
sh validation/complexity.sh                  # no Perl needed
sh validation/speed.sh                       # ~2.5 min; run on a quiet machine
```

The Perl model is built once per dataset and reused; delete `out/<dataset>/perl/` to
rebuild it.

## Reading the results

**Tier 1** — `out/<dataset>/results/tier1_compat.txt` compares `build --perl_compat` with
Perl. On chr22 the GC bins, transition counts, TE and simple-repeat records and nested
inserts are identical; the k-mer tables differ by at most ~0.005 because GARLIC 1.4 leaves
the last bases of each repeat unmasked (a bug `rusted-garlic` does not reproduce).
`tier1_default.txt` shows the default build, which differs by design.

**Tier 2** — `out/<dataset>/results/report.tsv`. Emission is random, so the test asks
whether Perl and Rust differ more than two runs of the same program do. In the
`## Equivalence` section, `ratio` (distributions) should be below 1.5 and `z` (scalars)
below 3, and every `perl_vs_rust_compat` row should read `OK`. With ~12 checks per pair
an occasional z of 3–4 can be chance; rerun and see whether it repeats. `perl_vs_rust`
rows are expected to show `DIFF`: those are the intentional changes. `repeats.tsv`
compares simulated repeats with the real TEs; `inserts.tsv` compares inserted-element
counts and lengths by class.

On sm-example, Perl sometimes stops inserting TEs early (see Known gaps in
`../VALIDATION.md`), so `lower%` and `gc_mean` can show borderline `DIFF`s for
`perl_vs_rust_compat`.

**Complexity** — `out/complexity/compare.txt`: column 2 (max difference from
`SeqComplex.pm`) should be 0 for every measure except `cz`. `NA:1` is the 1-bp record,
which `SeqComplex.pm` cannot compute.

The test windows (`tests/data/complexity/windows.fa`) and SeqComplex.pm's values for them
(`seqcomplex_reference.tsv`) are committed, so neither this script nor `cargo test` needs
Perl. `--regenerate` reruns SeqComplex.pm on the same windows and rewrites the reference.
The result is not byte-identical from run to run: `cz` varies (Perl compresses a randomly
named temp file) and some values change in the 15th digit (Perl's hash order changes
the order of floating-point additions).

## Files

| Path | Contents |
|---|---|
| `config.sh` | shared settings and helper functions |
| `scripts/compare_models.py` | Tier 1 model comparison |
| `scripts/compare_emitted.py` | Tier 2 sequence statistics and equivalence test |
| `scripts/compare_inserts.py` | inserted-element counts and lengths by class |
| `scripts/compare_repeats.py` | simulated repeats vs real TEs |
| `scripts/perl_complexity.pl`, `compare_complexity.py` | complexity check |
| `scripts/make_windows.py` | how `tests/data/complexity/windows.fa` was made (chr22 + a simulated sequence) |
| `scripts/embl2fa.py` | FamDB EMBL export → FASTA keyed by family name |
| `data/sm-example.fa.align` | RepeatMasker alignment of `../tests/data/sm-example.fa` (input to `createModel.pl`) |
