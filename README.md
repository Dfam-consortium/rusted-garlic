# rusted-garlic

A Rust port of [GARLIC](https://doi.org/10.1093/nar/gku356) — the **G**enome **AR**tificia**L** **I**ntergenic **C**ompendium sequence generator, originally written by Juan Caballero and Gwenlyn Glusman.

`rusted-garlic` reads a real genome (or any FASTA), builds a statistical model of its sequence composition, and then generates arbitrarily large synthetic DNA sequences that match the original genome's:

- GC isochore structure (configurable window size, default 1 kb, smoothed over ±5 flanking windows)
- k-mer Markov base frequencies per GC bin (configurable k, default 4)
- GC bin transition probabilities between adjacent windows
- Interspersed repeat (TE) content, evolved from consensus sequences [optional]
- Simple/tandem repeat content [optional]

The model is saved to a compact JSON file so you only need to process the real genome once.

---

## Building

```
cargo build --release
```

The binary is at `target/release/rusted-garlic`.

Run the tests with `cargo test --release`. `tests/determinism.rs` checks that `build`
and seeded `emit` reproduce the golden files in `tests/data/golden/`; see
`tests/data/README.md` for how to update them after an intended output change.

### Dependencies

All dependencies are declared in `Cargo.toml` and fetched automatically by Cargo:

| Crate | Purpose |
|---|---|
| `clap 4` | CLI argument parsing |
| `rand 0.8` | Random number generation |
| `serde` / `serde_json` | Model serialization |
| `flate2` | Transparent gzip decompression; in-process compression for `cz` complexity measure |
| `log` / `env_logger` | Structured logging with `--log-level` and `RUST_LOG` override |
| `itertools` | Iterator utilities |

---

## Quick start

`rusted-garlic` has three subcommands: `build` (create a model from a genome), `emit` (generate sequences from a model), and `complexity` (compute sequence complexity measures over sliding windows).

### 1. Build a model from a genome FASTA

```bash
rusted-garlic build --in_seq genome.fa --model genome.model
```

This makes a single pass over the FASTA, computes GC content per window (default 1 kb), assigns windows to up to 25 equal-population GC bins, and counts 4-mer frequencies within each bin.  Progress is printed to stderr.

Use `-w`, `-k`, and `-n` to override the defaults:

```bash
rusted-garlic build --in_seq genome.fa --model genome.model \
  --win_size 500 --kmer_size 5 --num_gc_bins 20
```

> **Note:** window size, k-mer size, and bin count are embedded in the model file and automatically applied when generating sequences — you do not need to repeat them at emit time.

In this simple case, the sequence model is being trained on the full-genomic sequence including (if present), transposable element sequences, simple/tandem repeats,
and gene regions.  See the full workflow below to generate true-background sequence models.

### 2. Generate an artificial sequence from the saved model

```bash
rusted-garlic emit --model genome.model --out_seq artificial.fa --seq_size 10m
```

`--seq_size` accepts bare integers or the suffixes `k`, `m`, `g` (case-insensitive).  Without `--seq_size` the output is the same length as the input genome.

### 2a. Generate multiple sequences in one run

Use `-N` / `--num_seqs` to write several independent sequences into a single FASTA file:

```bash
rusted-garlic emit --model genome.model --out_seq artificial.fa \
  --seq_size 5m --num_seqs 10
```

The sequences are named `artificial_sequence_1`, `artificial_sequence_2`, … `artificial_sequence_N`.  When `--num_seqs 1` (the default) the name is just `artificial_sequence`.  Each sequence is drawn independently from the same model with its own random seed, so repeat fractions and GC bin walks differ between sequences.

---

## Full workflow with true-background sequence model, and feature simulation

This is the intended high-fidelity workflow for genome benchmark data.

### Step 1 — Obtain and prepare annotation files

You need:

| File | Description |
|---|---|
| `genome.fa[.gz]` | Reference genome FASTA |
| `genome.te_annot.tsv[.gz]` | Interspersed repeat (TE) annotations (see below) |
| `genome.tandem_annot.tsv[.gz]` | Tandem/simple repeat annotations (UCSC TRF BED format) |
| `genome.gene_annot.tsv[.gz]` | Gene body + exon annotations (see format below) |
| `RepeatMaskerLib.fa[.gz]` | TE consensus FASTA (from RepBase / Dfam) |

All input files are auto-detected as gzip-compressed when their name ends in `.gz` and are decompressed transparently — no manual `gunzip` required.

**Preparing the TE annotation TSV from RepeatMasker output:**

`--te_annot` expects a 12-column tab-separated file (see [TE annotation format](#te-annotation-format) below).  Use the provided converter to produce it from a RepeatMasker `.align` file:

```bash
python3 utils/rm_align_to_tsv.py genome.fa.align > genome.te_annot.tsv
# or directly from a gzip-compressed file:
python3 utils/rm_align_to_tsv.py genome.fa.align.gz > genome.te_annot.tsv
```

**Preparing the tandem repeat annotation file:**

`--tandem_annot` expects **UCSC TRF BED format** — this is both the masking source and the simulation profile source.  The format requires at minimum 9 tab-separated columns with the literal string `trf` in column 3:

```
chrom  start  end  trf  score  period  copies  percent_match  percent_indel  ...  consensus
```

All coordinates are 0-based half-open.  The consensus sequence must be in the last column.  UCSC TRF track downloads are already in this format:

```bash
# From a UCSC TRF BED download (already in the right format):
cp hg38.trf.bed genome.tandem_annot.tsv
```

> **Note:** If you supply a file with fewer than 9 columns or without `trf` in column 3, the masking pass will still work (it only needs columns 0–2), but no tandem repeat profile will be built and `--sim_frac` will have no effect.

**Preparing the gene annotation file:**

`--gene_annot` uses a 7-column tab-separated format (0-based half-open):

```
chrom  txStart  txEnd  name  strand  exonStarts  exonEnds
```

Including exon coordinates (columns 5–6) is recommended — they enable `--mask_exon` and make the profiling filter exon-aware.  Columns 3–4 (`name`, `strand`) are informational only.  Use `.` for any optional column when the data are not available; gene-body masking and the profiling filter still work with only columns 0–2.

The gene annotation serves two purposes:
1. **Masking** — gene bodies (default) or exon intervals (`--mask_exon`) are excluded from k-mer counting
2. **Profiling filter** — the same intervals used for masking also act as an exclusion filter during repeat profiling: TEs and tandem repeats that overlap those intervals are excluded from their respective simulation profiles (matching Perl's intron filter behavior)

```bash
# Recommended: include exon coordinates (enables --mask_exon and exon-aware profiling filter)
awk 'NR>1 {
    name=$2; chrom=$3; strand=$4; txStart=$5; txEnd=$6;
    exStarts=$10; exEnds=$11;
    print chrom"\t"txStart"\t"txEnd"\t"name"\t"strand"\t"exStarts"\t"exEnds
}' refGene.txt > genome.gene_annot.tsv

# Alternative: gene bodies only (exon columns omitted — gene-body masking and profiling filter still work)
awk 'NR>1 {print $3"\t"$5"\t"$6"\t"$2"\t"$4"\t.\t."}' refGene.txt > genome.gene_annot.tsv
```

### Step 2 — Build a masked model with repeat profiling

```bash
rusted-garlic build \
  --in_seq       genome.fa \
  --model        genome.model \
  --te_annot     genome.te_annot.tsv \
  --tandem_annot genome.tandem_annot.tsv \
  --gene_annot   genome.gene_annot.tsv
```

What each flag does:

- `--te_annot` — masks interspersed repeats before k-mer counting; profiles TE divergence/insertion/deletion rates per GC bin; TEs overlapping gene regions are excluded from the profile when `--gene_annot` is provided
- `--tandem_annot` — masks tandem/simple repeats before k-mer counting; profiles tandem repeat consensus sequences per GC bin; tandem repeats contained within a TE interval *or* overlapping a gene region are excluded from the profile
- `--gene_annot` — masks gene bodies (default) or exons (`--mask_exon`) before k-mer counting; also acts as a profiling filter for both TE and tandem annotations
- `--mask_exon` — mask only exon intervals instead of whole gene bodies (requires exon coordinates in `--gene_annot`)

Two Perl-compatible sidecar files are written alongside the model:

```
genome.repeats.W1000.data   # TE repeat records per GC bin
genome.inserts.W1000.data   # Nested insertion frequency table
```

### Step 3 — Generate a sequence with TE bombardment

```bash
rusted-garlic emit \
  --model    genome.model \
  --out_seq  artificial.fa \
  --seq_size 10m \
  --te_cons  RepeatMaskerLib.fa \
  --te_frac  30 \
  --sim_frac 1
```

`--te_cons` is required to enable TE insertion.  Without it the output is a base sequence only (upper-case, no repeats).

The bombardment process:

1. A base sequence is generated by the 4-mer Markov model with GC isochore transitions.
2. Evolved TE copies are inserted (lower-cased) until `--te_frac` percent of the sequence is covered.  Each copy is drawn from the per-GC-bin repeat profile, a random fragment of the consensus is taken, and substitutions/indels are applied using the local background k-mer model to bias mutations toward contextually plausible sequences.
3. Simple/tandem repeats are inserted similarly until `--sim_frac` percent coverage.

If `--te_frac` or `--sim_frac` are omitted the fractions are drawn randomly each run (10–60 % and 0–2 % respectively, matching the Perl defaults).

---

## All command-line options

### Global options

| Flag | Default | Description |
|---|---|---|
| `--log-level <LEVEL>` | `info` | Log verbosity: `error`, `warn`, `info`, `debug`, `trace`, or `off`; overridden by `RUST_LOG` environment variable |

### `build` subcommand

```
rusted-garlic build --in_seq <FASTA_FILE> --model <MODEL_FILE> [OPTIONS]
```

**Required**

| Flag | Description |
|---|---|
| `-i, --in_seq <FASTA_FILE>` | Input genome FASTA |
| `-m, --model <MODEL_FILE>` | Model file to save |

**Annotation / masking**

| Flag | Description |
|---|---|
| `--te_annot <TSV_FILE>` | 12-column TE annotation TSV for masking + TE profiling (see [TE annotation format](#te-annotation-format)) |
| `--tandem_annot <TSV_FILE>` | UCSC TRF BED file for masking + tandem repeat profiling (see [Annotation file formats](#annotation-file-formats)) |
| `--gene_annot <TSV_FILE>` | 7-column gene annotation file for masking + profiling filter (see [Annotation file formats](#annotation-file-formats)) |
| `--mask_exon` | Mask exon intervals only instead of whole gene bodies (requires exon columns in `--gene_annot`) |
| `--filter_te_tandems` | Drop tandem repeats fully contained in a TE from the simple-repeat profile (avoids double counting; default keeps them, like GARLIC 1.4) |
| `--perl_compat` | Reproduce GARLIC 1.4 GC binning and TE bin positions exactly, for comparison with the Perl version (see below) |

All annotation files use 0-based half-open coordinates.  TEs and tandem repeats that overlap gene regions are excluded from their simulation profiles.  Tandem repeats fully contained within a TE interval are also excluded from the tandem profile.

**Model parameters (stored in the model file)**

| Flag | Default | Description |
|---|---|---|
| `-w, --win_size <INT>` | `1000` | Window size in bp for GC isochore binning |
| `-k, --kmer_size <INT>` | `4` | k-mer size for the Markov base-frequency model |
| `-n, --num_gc_bins <INT>` | `25` | Maximum number of GC isochore bins |

### `emit` subcommand

```
rusted-garlic emit --model <MODEL_FILE> --out_seq <FASTA_FILE> [OPTIONS]
```

**Required**

| Flag | Description |
|---|---|
| `-m, --model <MODEL_FILE>` | Model file to load |
| `-o, --out_seq <FASTA_FILE>` | Write artificial sequence(s) here |

**Sequence generation**

| Flag | Default | Description |
|---|---|---|
| `-s, --seq_size <SIZE>` | input length | Output length per sequence: integer with optional `k`/`m`/`g` suffix; not valid with `evenly-spaced` mode |
| `-N, --num_seqs <INT>` | `1` | Number of independent sequences to generate |

**Internal repeat bombardment (legacy)**

| Flag | Default | Description |
|---|---|---|
| `--te_cons <FASTA_FILE>` | — | TE consensus FASTA — enables internal TE bombardment |
| `--te_frac <PERCENT>` | random 10–60 | Target TE coverage 0–100; also used by `coverage` insertion mode |
| `--sim_frac <PERCENT>` | random 0–2 | Target simple-repeat fraction 0–100 |

**External feature insertion** *(flag validation implemented; insertion logic not yet complete)*

| Flag | Default | Description |
|---|---|---|
| `--insert_features <FASTA_FILE>` | — | Feature sequences to insert into the background |
| `--insert_mode <MODE>` | `coverage` | `coverage`, `full`, or `evenly-spaced` (see below) |
| `--spacing <BP>` | — | Background bp between features; required for `evenly-spaced`, disallowed otherwise |
| `--allow_nesting` | off | Allow features to overlap prior insertions; `coverage` and `full` only |

> **Note:** `--insert_features` and its associated flags are parsed and validated but the insertion logic is not yet implemented.  The flags will be accepted without error but the feature sequences will not be inserted into the output.

**Output**

| Flag | Default | Description |
|---|---|---|
| `--no_soft_mask` | off | Keep inserted repeats uppercase (default: lower-case) |
| `--perl_compat` | off | Reproduce known GARLIC 1.4 bugs, for side-by-side comparison with the Perl version (see below) |
| `--out_bed <BED_FILE>` | — | Write legacy bombardment repeat insertion records as a 6-column BED file |
| `--seed <INT>` | random | Random seed. The same seed, model, inputs and rusted-garlic version reproduce the output exactly. The seed used is always logged. |

#### Insertion modes *(planned)*

| Mode | Required extras | Disallowed | Behaviour |
|---|---|---|---|
| `coverage` | `--te_frac` (or accepts random default) | `--spacing` | Draw features until `--te_frac` percent coverage is reached |
| `full` | — | `--spacing` | Insert every feature in the FASTA exactly once |
| `evenly-spaced` | `--spacing` | `--seq_size`, `--allow_nesting` | Place features in order with exactly `--spacing` bp of background between each; output length is determined by the feature set |

### `complexity` subcommand

```
rusted-garlic complexity [OPTIONS] <FASTA_FILE> [<FASTA_FILE> ...]
```

Computes a suite of composition and complexity measures for every sliding window across one or more FASTA files.  Output is a tab-separated table written to stdout (or `--out`) with one row per window.

**Inputs**

| Argument | Description |
|---|---|
| `<FASTA_FILE> …` | One or more FASTA files (plain or `.gz`); each sequence is processed independently |

**Windowing**

| Flag | Default | Description |
|---|---|---|
| `-w, --win <INT>` | sequence length | Window size in bp; omit to treat each sequence as one window |
| `-s, --step <INT>` | `--win` | Step between window starts; defaults to `--win` (non-overlapping); set smaller for overlapping windows |
| `-n, --nmax <FLOAT>` | `0.3` | Skip windows where the non-ACGT fraction exceeds this threshold |

**K-mer sizes**

| Flag | Default | Description |
|---|---|---|
| `--min_kmer <INT>` | `1` | Smallest k-mer word size to include in k-mer measures (1–6) |
| `--max_kmer <INT>` | `6` | Largest k-mer word size to include in k-mer measures (1–6) |

**Output**

| Flag | Default | Description |
|---|---|---|
| `-o, --out <FILE>` | stdout | Write TSV here instead of stdout |

#### Output columns

The first six columns identify the window; the remainder are measures:

| Column | Description |
|---|---|
| `file` | Source FASTA filename |
| `seq` | Sequence (record) name |
| `start` | 0-based start position |
| `end` | 0-based end position (exclusive) |
| `win_len` | Window length in bp |
| `acgt_count` | Number of ACGT bases (non-ACGT excluded from denominators) |
| `gc` | GC content: (C+G) / total |
| `at` | AT content: (A+T) / total |
| `gcs` | GC skew: \|G−C\| / (G+C) |
| `ats` | AT skew: \|A−T\| / (A+T) |
| `cpg` | CpG ratio: observed CG dinucleotides / (C_count × G_count) |
| `ket` | Keto skew: \|(G+T)−(A+C)\| / total |
| `pur` | Purine skew: \|(G+A)−(T+C)\| / total |
| `cwf` | Wootton-Federhen complexity |
| `ce` | Shannon entropy (bits) |
| `cz` | Compression ratio (raw length / gzip length); higher = more repetitive |
| `cm<k>` | Markov complexity at word size k (one column per k in min_kmer..=max_kmer) |
| `ct<k>` | Trifonov complexity at word size k |
| `cl<k>` | Linguistic complexity at word size k |

#### Examples

Whole-sequence measures for every record in a FASTA:

```bash
rusted-garlic complexity genome.fa > genome.complexity.tsv
```

Sliding 1 kb windows, 500 bp step, skipping windows with > 10% non-ACGT:

```bash
rusted-garlic complexity --win 1000 --step 500 --nmax 0.1 chr22.fa > chr22.W1000.tsv
```

Non-overlapping 500 bp windows, only 2-mers through 4-mers, output to file:

```bash
rusted-garlic complexity --win 500 --min_kmer 2 --max_kmer 4 --out out.tsv genome.fa
```

Multiple FASTA files in one run:

```bash
rusted-garlic complexity --win 1000 chr1.fa chr2.fa chr3.fa > all_chrs.tsv
```

---

#### `--perl_compat`

`build --perl_compat` and `emit --perl_compat` reproduce GARLIC 1.4 behaviour where
`rusted-garlic` deliberately differs, so the two programs can be compared (see
`tier2/` in the project root).  They are meant for validation and will eventually be
removed.

`build --perl_compat` (`createModel.pl`): GC bins, GC transitions, and the TE,
simple-repeat and nested-insert tables are then identical to GARLIC 1.4 on human
chr22.  K-mer frequencies still differ slightly because of Perl's repeat-masking bug
(item 2 under "Differences from the original Perl GARLIC"), which is not reproduced.

| GARLIC 1.4 behaviour | Default | `build --perl_compat` |
|---|---|---|
| GC bin cuts: one pass with bucket size `total / num_gc_bins`; a GC value holding several buckets' worth of windows merges them into one large bin, and the highest observed GC value never starts a bin | bin count lowered until every bin is non-empty, keeping populations even | Perl cuts |
| Window GC histogram slot | `round(1000 * gc)` | `int(1000 * gc)` |
| Chromosome's trailing partial window | binned by its own GC | last full window's GC repeated |
| Smoothing span | ±5 windows | ±`int(10 kb / win_size / 2)` windows |
| Bin ranges (and the k-mer fallback GC, which uses the range end) | midpoints between adjacent bins' observed GC | Perl cut points |
| TE GC-bin position | midpoint of the span (0-based) | Perl: 1-based start, end of the last fragment read |

`emit --perl_compat` (`createFakeSequence.pl`) reproduces these GARLIC 1.4 bugs:

| GARLIC 1.4 bug | Default | `--perl_compat` |
|---|---|---|
| Background windows: `createSubSeq` re-appends the (k-1)-base seed, so each "1 kb" window is win-k+2 bp; the shortfall is padded with a copy of a random internal segment | exact `win_size` windows | 998-bp windows with duplicated seed, padded end |
| Evolved-repeat insertions: base sampled for the context before `pos` is placed after `seq[pos]` (raises CpG in repeats) | inserted in context | Perl placement |
| TE records longer than their consensus: negative `rand()` offset makes `substr` return a random 3′ piece | full consensus | random 3′ piece |

## Model file format

The model is a single JSON file containing:

- `win_size` and `kmer_size` — the parameters used when the model was built (applied automatically at emit time)
- Per-bin k-mer conditional frequencies (Markov model)
- GC bin transition matrix (normalized, plus raw counts; `emit` adds a +1 pseudocount per cell like GARLIC 1.4)
- Per-bin k-mer totals (weights for drawing a random GC bin during repeat mutation)
- GC bin boundaries and window counts
- Optional `repeat_model` object with per-bin TE and simple repeat records plus the nested-insertion frequency table

The file is human-readable but not designed for hand-editing.

> **Format version:** each model records `model_version`; `emit` refuses models from other versions with a "rebuild the model" error.  Models built by rusted-garlic 0.1.x (no version field) must be rebuilt with 0.2.0 — see `CHANGELOG.md`.

> **Compatibility note:** Models built with different `win_size` or `kmer_size` values are not interchangeable.  If you rebuild a model with new parameters, regenerate any downstream sequences from the new model file.

---

## Tutorial: human chr22

This walks through a realistic small example using chr22.

### Prepare files

```bash
# Download chr22 FASTA (example assumes you have it locally)
# Run RepeatMasker:
RepeatMasker -pa 16 -species human -align chr22.fa

# Convert .align output to the TE annotation TSV:
python3 utils/rm_align_to_tsv.py chr22.fa.align > chr22.te_annot.tsv

# Get tandem repeat intervals from UCSC TRF track (already BED-like, 0-based half-open):
cp hg38_trf.bed chr22.tandem_annot.tsv

# Extract gene annotation from a UCSC refGene table (7-column format with exon coordinates):
awk 'NR>1 {print $3"\t"$5"\t"$6"\t"$2"\t"$4"\t"$10"\t"$11}' refGene.txt > chr22.gene_annot.tsv

# Obtain consensus sequences — many labs use the Dfam release FASTA or
# the RepeatMaskerLib.embl converted to FASTA with a one-liner.
```

### Build the model

```bash
rusted-garlic build \
  --in_seq       chr22.fa \
  --model        chr22.model \
  --te_annot     chr22.te_annot.tsv \
  --tandem_annot chr22.tandem_annot.tsv \
  --gene_annot   chr22.gene_annot.tsv
```

Typical stderr output:

```
Building model from "chr22.fa"...
Loading gene annotations from "chr22.gene_annot.tsv"...
  34091 gene body intervals loaded
Loading TE annotations from "chr22.te_annot.tsv"...
  184302 masked intervals loaded
Loading tandem repeat annotations from "chr22.tandem_annot.tsv"...
  22187 masked intervals loaded
  34091 gene intervals added to mask
Total mask: 240580 intervals across all sequences
Saving model to "chr22.model"...
  Profiling TE repeats from "chr22.te_annot.tsv"...
  184302 TE records profiled
  Profiling tandem repeats from "chr22.tandem_annot.tsv"...
  21954 simple repeat records profiled
Writing repeat profile to "chr22.repeats.W1000.data"...
Writing insert table to "chr22.inserts.W1000.data"...
```

### Generate a 10 Mb test sequence

```bash
rusted-garlic emit \
  --model    chr22.model \
  --out_seq  chr22_artificial.fa \
  --seq_size 10m \
  --te_cons  RepeatMaskerLib.fa \
  --te_frac  32 \
  --sim_frac 1
```

### Generate multiple independent sequences

To produce a benchmark set of 20 independent 5 Mb sequences in one pass:

```bash
rusted-garlic emit \
  --model    chr22.model \
  --out_seq  chr22_bench.fa \
  --seq_size 5m \
  --num_seqs 20 \
  --te_cons  RepeatMaskerLib.fa \
  --te_frac  32 \
  --sim_frac 1
```

The output FASTA will contain 20 records named `artificial_sequence_1` … `artificial_sequence_20`, each generated independently from the same model.

### Build with non-default model parameters

To use 5-mers and 500 bp windows (e.g. for a GC-rich or compact genome):

```bash
rusted-garlic build \
  --in_seq       chr22.fa \
  --model        chr22_k5w500.model \
  --te_annot     chr22.te_annot.tsv \
  --tandem_annot chr22.tandem_annot.tsv \
  --gene_annot   chr22.gene_annot.tsv \
  --kmer_size    5 \
  --win_size     500 \
  --num_gc_bins  20
```

Then generate sequences as usual — the stored parameters are applied automatically:

```bash
rusted-garlic emit \
  --model    chr22_k5w500.model \
  --out_seq  artificial_k5.fa \
  --seq_size 10m \
  --te_cons  RepeatMaskerLib.fa
```

### Verify repeat content

The inserted regions are lower-case; you can count them with:

```bash
python3 -c "
import re, sys
seq = ''.join(l.strip() for l in open('chr22_artificial.fa') if not l.startswith('>'))
total = len(seq)
lc = sum(1 for c in seq if c.islower())
print(f'Total: {total} bp')
print(f'Repeat: {lc} bp  ({100*lc/total:.1f}%)')
"
```

---

## Annotation file formats

All annotation files use 0-based half-open coordinates and may be plain or gzip-compressed.  Lines beginning with `#` are treated as comments and skipped.

### Tandem repeat annotation format (`--tandem_annot`)

`--tandem_annot` accepts **UCSC TRF BED format**.  This file is used for both masking (columns 0–2 only) and building the tandem repeat simulation profile (columns 0–2 + 3–last).  At minimum 9 whitespace-delimited fields are required; the literal string `trf` must appear in column 3:

| Col | Name | Description |
|---|---|---|
| 0 | `chrom` | Sequence / chromosome name |
| 1 | `start` | 0-based half-open start |
| 2 | `end` | 0-based half-open end |
| 3 | `trf` | Literal string `trf` (identifies the format) |
| 4 | (score) | Ignored |
| 5 | `period` | Repeat unit period length |
| 6 | (copies) | Ignored |
| 7 | `pct_match` | Percent match (used to derive divergence: 100 − pct_match) |
| 8 | `pct_indel` | Percent indel rate |
| … | … | Additional columns ignored |
| last | `consensus` | Tandem repeat consensus sequence |

> If the file has fewer than 9 columns or column 3 is not `trf`, every line will be skipped for profiling (a warning is issued).  Masking still works from columns 0–2.

UCSC TRF track downloads are already in this format and require no conversion.

### Gene annotation format (`--gene_annot`)

`--gene_annot` accepts a 7-column tab-separated file.  The full format is recommended — exon coordinates (columns 5–6) enable `--mask_exon` and exon-aware profiling filters.  Columns 0–2 are the only strictly required fields; use `.` for any optional column you don't have:

| Col | Name | Description |
|---|---|---|
| 0 | `chrom` | Sequence / chromosome name |
| 1 | `txStart` | Gene body 0-based start |
| 2 | `txEnd` | Gene body 0-based end (exclusive) |
| 3 | `name` | Gene/transcript name — informational, may be `.` |
| 4 | `strand` | `+` or `-` — informational, may be `.` |
| 5 | `exonStarts` | Comma-separated 0-based exon start positions; `.` if not available |
| 6 | `exonEnds` | Comma-separated 0-based exon end positions (exclusive); `.` if not available |

The gene annotation serves two roles:
1. **Masking** — gene bodies (default) or exon intervals (`--mask_exon`) are excluded from k-mer counting
2. **Profiling filter** — the same intervals used for masking also filter repeat profiling: TE instances and tandem repeats overlapping those regions are excluded from their simulation profiles

### TE annotation format (`--te_annot`) {#te-annotation-format}

The `--te_annot` flag accepts a 12-column tab-separated file with 0-based half-open coordinates (BED convention).  Lines beginning with `#` are treated as comments.

| Col | Name | Description |
|---|---|---|
| 0 | `chrom` | Sequence / chromosome name |
| 1 | `start` | 0-based half-open start |
| 2 | `end` | 0-based half-open end |
| 3 | `name` | Repeat name (e.g. `AluYa5`) |
| 4 | `family` | Class/family string (e.g. `SINE/Alu`) |
| 5 | `strand` | `+` or `-` |
| 6 | `div` | % substitution divergence |
| 7 | `ins` | % insertion rate |
| 8 | `del` | % deletion rate |
| 9 | `cons_start` | Consensus start coordinate |
| 10 | `cons_end` | Consensus end coordinate |
| 11 | `rep_id` | Unique instance ID (use `.` if not available) |

`rep_id` serves two purposes: it ties multiple fragments of the same repeat instance together (so divergence stats are taken as the max across all fragments), and it enables nested-insertion detection (elements appearing between two occurrences of the same `rep_id` in the file are recorded as nested inside that repeat).  Use `.` when the source data does not track instance identity; each row is then treated as an independent unfragmented element and nesting is not recorded.

Simple/tandem repeats must be excluded from this file — they are handled separately by `--tandem_annot`.

**Example rows:**
```
#chrom	start	end	name	family	strand	div	ins	del	cons_start	cons_end	rep_id
chr22	999	1306	AluYa5	SINE/Alu	+	12.40	0.60	1.20	1	308	m_b1s502i1
chr22	5000	5312	L1HS	LINE/L1	-	18.20	1.10	0.80	1	312	m_b2s001i3
```

### Converting from RepeatMasker

`utils/rm_align_to_tsv.py` converts a RepeatMasker `.align` file to this format:

```bash
python3 utils/rm_align_to_tsv.py genome.fa.align    > genome.te_annot.tsv
python3 utils/rm_align_to_tsv.py genome.fa.align.gz > genome.te_annot.tsv
```

The script filters out Simple_repeat, Low_complexity, Unknown, and Satellite entries and handles both positive- and negative-strand `.align` header layouts.

---

## Sidecar files

When `build` is run with `--te_annot` (and/or `--tandem_annot`), two text files are written alongside the model that are compatible with the original Perl `createFakeSequence.pl`.  The filenames always use the `W1000` suffix for Perl compatibility regardless of `--win_size`:

**`<stem>.repeats.W1000.data`** — one record per TE / simple repeat per GC bin:
```
#GC=0.0-42.3
AluYa5:SINE/Alu:+:12.4:0.6:1.2:308:1
...
SIMPLE:AGAT:+:4:3.2:1.1
```

**`<stem>.inserts.W1000.data`** — nested insertion frequency table:
```
#GC=0.0-42.3
AluSx#SINE/Alu    L1HS#LINE/L1    0.142857
```

---

## Architecture overview

```
src/
  main.rs          CLI (build / emit / complexity subcommands), orchestration
  seq_model.rs     GC binning, k-mer model, GC transition model, emit()
  repeat_model.rs  TE and simple repeat profiling from annotation files
  evolve.rs        Mutation operators, evolveRepeat, evolveSimple, bombardment loops
  consensus.rs     FASTA consensus loader with IUPAC resolution
  masking.rs       MaskSet and MaskCursor (streaming interval checker)
  complexity.rs    WindowCounts accumulator and all 12 complexity measures
  io/
    mod.rs         IntervalMap type, sort_and_merge helper
    fasta.rs       Streaming FASTA iterator (plain + gz)
    repeats.rs     12-column repeat TSV parser
    bed.rs         Tandem repeat (TRF BED) and gene annotation parsers
```

The build pipeline is three sequential passes over the input FASTA:

1. **Pass 1** — raw GC content per `win_size` window → GC bin assignments (up to `num_gc_bins` bins), window counts, transition counts
2. **Pass 2** — masked k-mer counting per GC bin (`kmer_size`-mer Markov model)
3. **Pass 3** — TE and simple repeat profiling from annotation files, assigned to GC bins by position (TE midpoint, simple-repeat start)

The three model parameters (`win_size`, `kmer_size`, `num_gc_bins`) are stored in the JSON model file and re-used automatically at emit time.

Emission repeats the following for each of the `num_seqs` requested sequences:

1. Markov walk using per-GC-bin k-mer tables with GC bin transitions every `win_size` bases
2. (Optional) TE bombardment — random insertions until `te_frac` percent coverage
3. (Optional) Simple repeat bombardment — random insertions until `sim_frac` percent coverage

All sequences are written to a single FASTA file.  When `num_seqs > 1` each record is named `artificial_sequence_<N>`.

---

## Differences from the original Perl GARLIC

`rusted-garlic` is a faithful re-implementation of the original `createModel.pl` / `createFakeSequence.pl` workflow, but several deliberate changes were made to fix known bugs, improve display clarity, and modernise the tooling.  The differences below affect output and are worth understanding when comparing results between the two versions.

### Structural / workflow changes

| Aspect | Perl GARLIC | `rusted-garlic` |
|---|---|---|
| Entry points | Two scripts: `createModel.pl` (build) and `createFakeSequence.pl` (emit) | Single binary with `build` and `emit` subcommands |
| Repeat input format | RepeatMasker `.align` / `.out` parsed directly | 12-column TSV via `--te_annot` (convert with `utils/rm_align_to_tsv.py`); UCSC TRF BED via `--tandem_annot`; 7-column gene annotation via `--gene_annot` |
| Model format | Multiple text sidecar files (`.kmer`, `.gc`, `.transition`, …) | Single compact JSON file; two Perl-compatible sidecar files (`.repeats.W1000.data`, `.inserts.W1000.data`) are also written when `--te_annot` is provided |
| Genome data download | `createModel.pl` can auto-fetch annotation files from UCSC | Not supported — provide files explicitly |
| Gzip support | Manual decompression required | All input files auto-detected and decompressed transparently |

### Algorithmic differences that change output

#### 1. GC binning

Both versions compute per-window GC, build a 0.1%-resolution histogram, cut it into up
to `num_gc_bins` equal-population bins, smooth each window's GC over its neighbours
and assign it to a bin.  The cut rule differs.  GARLIC 1.4 makes one pass with a fixed
bucket size, so when a single GC value holds more than one bucket's worth of windows
(very narrow isochores, or large `win_size`) it produces one oversized bin and keeps
many small ones.  `rusted-garlic` lowers the bin count until every bin is non-empty,
which keeps populations even at the cost of fewer bins.  With 1 kb windows on typical
genomes (human, AT-rich, bimodal) both rules give the same bins.  Smaller differences
(histogram rounding, the trailing partial window, smoothing span, bin range
boundaries) are listed under `--perl_compat` above; `build --perl_compat` reproduces
GARLIC 1.4 exactly.

In both modes the smoothed GC is the mean of the neighbours' full-precision GC,
rounded once to 0.1% (as Perl), and k-mers are binned by the window of their first
base (as Perl).

#### 2. Repeat masking covers the whole repeat

`maskRepeat` and `maskTRF` in GARLIC 1.4 mask `end − start − 1` bases, leaving the last
2 bases of each RepeatMasker TE and the last base of each TRF repeat unmasked, so k-mers
overlapping repeat ends are counted.  Rust masks the full annotated interval in both
modes.  On chr22 this gives about 0.5% fewer k-mers per bin; with `build --perl_compat`
the conditional k-mer frequencies are within 0.006 of Perl's.

#### 3. K-mer fallback for unobserved prefixes

When a k-mer prefix has no observations in a bin, both versions use a GC-based prior:
`p` = the bin's upper GC bound, A and T get `(1 − p) / 2`, C and G get `p / 2` (Perl:
`$c / 100` from the bin's range string).  The upper bound itself differs by default
(see bin ranges under `--perl_compat`).

#### 4. Repeat GC bin position

Each TE is assigned to the GC bin at the midpoint of its span; simple repeats use their
start.  GARLIC 1.4 also uses the TE midpoint, but computes it from the 1-based
RepeatMasker start and the end of the last fragment read (`--perl_compat` reproduces
this).  The effect is a few TEs per bin on chr22.

### Features in `rusted-garlic` not present in Perl GARLIC

| Feature | Description |
|---|---|
| `complexity` subcommand | Ports the separate `SeqComplexity.pm` / `gatherStats.pl` workflow into the binary; computes all 12 sequence complexity measures (composition, entropy, k-mer diversity, compression) over sliding windows with a single command |
| Profiling filter for TEs and tandem repeats | TEs overlapping gene regions and tandem repeats that are either contained within a TE interval or overlap a gene region are excluded from simulation profiles (the Perl code only filtered tandem/TE repeats that overlapped gene-masked sequence) |
| `--mask_exon` with unified profiling filter | Masking and profiling exclusion use the same interval set (gene bodies or exons); in Perl the exon-masking and gene-body profiling filter were separate code paths |
| Streaming FASTA reader | Processes arbitrarily large FASTA files without loading them into memory |
| In-process gzip compression | `cz` measure uses `flate2` in-process rather than shelling out to `gzip` with temp files |

### Features not yet fully implemented in `rusted-garlic`

| Feature | Status |
|---|---|
| `--insert_features` (external feature insertion) | CLI flags parsed and validated; insertion logic not yet complete |
| `--gc_post_mask` | Recompute GC bins after masking (Pass 1 uses masked sequence) — Perl only |
| `--revcomp_kmer` | Average k-mer frequencies with their reverse-complement — Perl only |
| `--write_mask_seq` | Write masked genome FASTA alongside the model — Perl only |
| Exclusion patterns | Skip sequences matching a regex (e.g. `_random`, `chrUn`) — Perl only |
| `--binsize` | Override the GC bin step size directly — Perl only |




## AI Assistance Disclosure
This tool was written with the assistance of AI coding agents (Claude, GitHub Copilot).

Correctness is validated by comparing output against Garlic v1.5 on a suite of real assembly datasets - not by manual code review alone. The background sequence modeling/simulation algorithm was hand-coded, however AI generated the full Garlic interface implementation; humans defined the validation criteria and verified the results.

See VALIDATION.md for full validation methodology.

AI rewrite stewardship based on https://rewrites.bio/ standard.
