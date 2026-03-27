# rusted-garlic

A Rust port of [GARLIC](https://doi.org/10.1093/nar/gku356) — the **G**enome **AR**tificia**L** **I**ntergenic **C**ompendium sequence generator, originally written by Juan Caballero and Gwenlyn Glusman.

`rusted-garlic` reads a real genome (or any FASTA), builds a statistical model of its sequence composition, and then generates arbitrarily large synthetic DNA sequences that match the original genome's:

- GC isochore structure (configurable window size, default 1 kb, smoothed over ±5 flanking windows)
- k-mer Markov base frequencies per GC bin (configurable k, default 4)
- GC bin transition probabilities between adjacent windows
- Interspersed repeat (TE) content, evolved from consensus sequences
- Simple/tandem repeat content

The model is saved to a compact JSON file so you only need to process the real genome once.

---

## Building

```
cargo build --release
```

The binary is at `target/release/rusted-garlic`.

### Dependencies

All dependencies are declared in `Cargo.toml` and fetched automatically by Cargo:

| Crate | Purpose |
|---|---|
| `clap 4` | CLI argument parsing |
| `rand 0.8` | Random number generation |
| `serde` / `serde_json` | Model serialization |
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

## Full workflow with repeat masking and TE bombardment

This is the intended high-fidelity workflow for genome benchmark data.

### Step 1 — Obtain and prepare annotation files

You need:

| File | Description |
|---|---|
| `genome.fa[.gz]` | Reference genome FASTA |
| `genome.repeats.tsv[.gz]` | Interspersed repeat annotations (see below) |
| `genome.trf.bed[.gz]` | UCSC TRF BED file (0-based half-open) |
| `refGene.txt[.gz]` | UCSC refGene/ensGene table |
| `RepeatMaskerLib.fa[.gz]` | TE consensus FASTA (from RepBase / Dfam) |

All input files are auto-detected as gzip-compressed when their name ends in `.gz` and are decompressed transparently — no manual `gunzip` required.

**Preparing the repeat TSV from RepeatMasker output:**

`--repeats` expects a 12-column tab-separated file (see [Repeat TSV format](#repeat-tsv-format) below).  Use the provided converter to produce it from a RepeatMasker `.align` file:

```bash
python3 utils/rm_align_to_tsv.py genome.fa.align > genome.repeats.tsv
# or directly from a gzip-compressed file:
python3 utils/rm_align_to_tsv.py genome.fa.align.gz > genome.repeats.tsv
```

### Step 2 — Build a masked model with repeat profiling

```bash
rusted-garlic build \
  --in_seq    genome.fa \
  --model     genome.model \
  --repeats   genome.repeats.tsv \
  --trf       genome.trf.bed \
  --genes     refGene.txt
```

What each masking flag does:

- `--repeats` — masks interspersed repeats before k-mer counting; also profiles TE divergence/insertion/deletion rates per GC bin
- `--trf` — masks simple/tandem repeats before k-mer counting; also profiles simple repeat units and rates per GC bin
- `--genes` — masks gene bodies before k-mer counting (no repeat profiling)

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
| `-r, --repeats <TSV_FILE>` | Repeat annotation TSV for masking + TE profiling (see [Repeat TSV format](#repeat-tsv-format)) |
| `-t, --trf <BED_FILE>` | UCSC TRF BED for masking + simple repeat profiling |
| `-g, --genes <GENE_FILE>` | UCSC refGene/ensGene table for gene-body masking |

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

**External feature insertion**

| Flag | Default | Description |
|---|---|---|
| `--insert_features <FASTA_FILE>` | — | Feature sequences to insert into the background |
| `--insert_mode <MODE>` | `coverage` | `coverage`, `full`, or `evenly-spaced` (see below) |
| `--spacing <BP>` | — | Background bp between features; required for `evenly-spaced`, disallowed otherwise |
| `--allow_nesting` | off | Allow features to overlap prior insertions; `coverage` and `full` only |

**Output**

| Flag | Default | Description |
|---|---|---|
| `--no_soft_mask` | off | Keep inserted repeats uppercase (default: lower-case) |
| `--out_bed <BED_FILE>` | — | Write repeat insertion records as a 6-column BED file |

#### Insertion modes

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

## Model file format

The model is a single JSON file containing:

- `win_size` and `kmer_size` — the parameters used when the model was built (applied automatically at emit time)
- Per-bin k-mer conditional frequencies (Markov model)
- GC bin transition matrix
- GC bin boundaries and window counts
- Optional `repeat_model` object with per-bin TE and simple repeat records plus the nested-insertion frequency table

The file is human-readable but not designed for hand-editing.

> **Compatibility note:** Models built with different `win_size` or `kmer_size` values are not interchangeable.  If you rebuild a model with new parameters, regenerate any downstream sequences from the new model file.

---

## Tutorial: human chr22

This walks through a realistic small example using chr22.

### Prepare files

```bash
# Download chr22 FASTA (example assumes you have it locally)
# Run RepeatMasker:
RepeatMasker -pa 16 -species human -align chr22.fa

# Convert .align output to the repeat TSV format:
python3 utils/rm_align_to_tsv.py chr22.fa.align > chr22.repeats.tsv

# Get TRF BED from UCSC (or run TRF yourself and convert to BED)
# Get refGene from UCSC Table Browser

# Obtain consensus sequences — many labs use the Dfam release FASTA or
# the RepeatMaskerLib.embl converted to FASTA with a one-liner.
```

### Build the model

```bash
rusted-garlic build \
  --in_seq  chr22.fa \
  --model   chr22.model \
  --repeats chr22.repeats.tsv \
  --trf     chr22.trf.bed \
  --genes   refGene.txt
```

Typical stderr output:

```
Building model from "chr22.fa"...
Loading repeat annotations from "chr22.fa.align"...
  184302 masked intervals loaded
Loading TRF annotations from "chr22.trf.bed"...
  22187 masked intervals loaded
Loading gene annotations from "refGene.txt"...
  34091 masked intervals loaded
Total mask: 187044 intervals across all sequences
Saving model to "chr22.model"...
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
  --in_seq      chr22.fa \
  --model       chr22_k5w500.model \
  --repeats     chr22.repeats.tsv \
  --trf         chr22.trf.bed \
  --genes       refGene.txt \
  --kmer_size   5 \
  --win_size    500 \
  --num_gc_bins 20
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

## Repeat TSV format

The `--repeats` flag accepts a 12-column tab-separated file with 0-based half-open coordinates (BED convention).  Lines beginning with `#` are treated as comments.

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

Simple/tandem repeats must be excluded from this file — they are handled separately by `--trf`.

**Example rows:**
```
#chrom	start	end	name	family	strand	div	ins	del	cons_start	cons_end	rep_id
chr22	999	1306	AluYa5	SINE/Alu	+	12.40	0.60	1.20	1	308	m_b1s502i1
chr22	5000	5312	L1HS	LINE/L1	-	18.20	1.10	0.80	1	312	m_b2s001i3
```

### Converting from RepeatMasker

`utils/rm_align_to_tsv.py` converts a RepeatMasker `.align` file to this format:

```bash
python3 utils/rm_align_to_tsv.py genome.fa.align    > genome.repeats.tsv
python3 utils/rm_align_to_tsv.py genome.fa.align.gz > genome.repeats.tsv
```

The script filters out Simple_repeat, Low_complexity, Unknown, and Satellite entries and handles both positive- and negative-strand `.align` header layouts.

---

## Sidecar files

When `build` is run with `--repeats`, two text files are written alongside the model that are compatible with the original Perl `createFakeSequence.pl`.  The filenames always use the `W1000` suffix for Perl compatibility regardless of `--win_size`:

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
    bed.rs         TRF BED and UCSC gene-table parsers
```

The build pipeline is three sequential passes over the input FASTA:

1. **Pass 1** — raw GC content per `win_size` window → GC bin assignments (up to `num_gc_bins` bins), window counts, transition counts
2. **Pass 2** — masked k-mer counting per GC bin (`kmer_size`-mer Markov model)
3. **Pass 3** — TE and simple repeat profiling from annotation files, assigned to GC bins by midpoint position

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
| Repeat input format | RepeatMasker `.align` / `.out` parsed directly | Generic 12-column TSV (convert with `utils/rm_align_to_tsv.py`) |
| Model format | Multiple text sidecar files (`.kmer`, `.gc`, `.transition`, …) | Single compact JSON file; two Perl-compatible sidecar files (`.repeats.W1000.data`, `.inserts.W1000.data`) are also written when `--repeats` is provided |
| Genome data download | `createModel.pl` can auto-fetch annotation files from UCSC | Not supported — provide files explicitly |
| Gzip support | Manual decompression required | All input files auto-detected and decompressed transparently |

### Algorithmic differences that change output

#### 1. Trailing partial window at end of chromosome

When the last window of a chromosome is shorter than `win_size`, the two implementations handle it differently.

- **Perl**: reuses the GC value from the previous complete window for the partial window.
- **Rust**: computes the actual GC fraction of the partial window's bases and uses that value.

The difference is small (at most one window per chromosome) but will cause minor divergence in bin assignments for genomes with many short chromosomes/scaffolds.

#### 2. Bin boundary display (cosmetic — does not affect sequence generation)

Both implementations use the same equal-population binning algorithm and the same smoothed-GC lookup table for window assignment during sequence generation.  Only the displayed GC range differs.

- **Perl**: displays `end[i] = start[i+1] − 1` (raw integer boundary).
- **Rust**: displays the midpoint between the highest observed GC in bin *i* and the original start of bin *i+1* — `midpoint = (bin_actual_max[i] + original_start[i+1]) / 2`.  This gives a more intuitive "halfway between the clusters" boundary.

Because `gc_lookup` (the actual per-window bin assignment table) is identical in both versions, the midpoint display change has no effect on sequences generated.

#### 3. Bug fix — floating-point GC fraction for k-mer fallback

When a k-mer prefix has no observations in a bin, a prior probability is used.  The GC fraction of a bin is `p_val = gc_ranges[bin].end / 1000.0`.

- **Perl**: performed integer division (`/ 1000`), always yielding 0 for bins below 1000, so the fallback was always `A/T = 0.5, C/G = 0.5` regardless of actual bin GC content.
- **Rust**: uses floating-point division (`/ 1000.0`), so the fallback correctly reflects each bin's GC fraction.

This only affects rare k-mer prefixes that were never observed in a given bin; common prefixes are unaffected.

#### 4. Bug fix — AT/GC frequency symmetry in k-mer fallback

In the same fallback path, the Perl code computed frequencies for A, T, C, and G, but contained an error where A and T were both assigned `q_val / 2` (the AT-fraction half) while C and G were also assigned `q_val / 2` instead of `p_val / 2`.  The result was that all four bases received equal probability (each ≈ 0.25) regardless of the bin's GC content.

- **Perl**: `freqA = freqT = freqC = freqG = q_val / 2` (bug — all equal).
- **Rust**: `freqA = freqT = q_val / 2`; `freqC = freqG = p_val / 2` — correctly GC-biased.

#### 5. Repeat midpoint vs. start for GC bin assignment

During repeat profiling (Pass 3), each annotated repeat is assigned to a GC bin based on its position within the genome.

- **Perl**: uses the repeat's start coordinate to look up the GC bin.
- **Rust**: uses the midpoint `(start + end) / 2` of the repeat, which is more representative for long elements that may span a GC bin boundary.

### Features in `rusted-garlic` not present in Perl GARLIC

| Feature | Description |
|---|---|
| `complexity` subcommand | Ports the separate `SeqComplexity.pm` / `gatherStats.pl` workflow into the binary; computes all 12 sequence complexity measures (composition, entropy, k-mer diversity, compression) over sliding windows with a single command |
| Streaming FASTA reader | Processes arbitrarily large FASTA files without loading them into memory |
| In-process gzip compression | `cz` measure uses `flate2` in-process rather than shelling out to `gzip` with temp files |

### Features present in Perl GARLIC not yet in `rusted-garlic`

The following `createModel.pl` / `createFakeSequence.pl` options have no equivalent in the current Rust version:

| Perl flag | Description |
|---|---|
| `--gc_post_mask` | Recompute GC bins after masking (Pass 1 uses masked sequence) |
| `--mask_exon` | Mask only exons rather than whole gene bodies |
| `--revcomp_kmer` | Average k-mer frequencies with their reverse-complement |
| `--write_mask_seq` | Write masked genome FASTA alongside the model |
| Exclusion patterns | Skip sequences matching a regex (e.g. `_random`, `chrUn`) |
| `--binsize` | Override the GC bin step size directly |
