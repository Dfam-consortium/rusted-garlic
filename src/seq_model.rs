use rand::distributions::WeightedIndex;
use rand::prelude::*;
use std::cmp;
use std::collections::{BTreeMap, HashMap};
use std::collections::VecDeque;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io::{BufReader, ErrorKind};
use std::path::Path;
use serde::{Deserialize, Serialize};
use serde_json::json;

use std::error::Error;

use crate::masking::MaskSet;
use crate::repeat_model::{profile_repeats, profile_trf, RepeatModel};
use crate::evolve::{insert_te_repeats, insert_simple_repeats, InsertRecord};

/// Equal-population GC bin cuts (per mille), as GARLIC 1.4 `calcBinGC`.
///
/// `hist[v]` counts windows with `int(1000 * gc) == v`.  A cut is placed at each
/// GC value `v` where `int(cumulative[v] / (total / num_bins))` increases; the
/// first cut is 0 and the last is 1000.  Bin `k` spans `cuts[k]..cuts[k + 1]`.
/// Buckets skipped in one step (a single value holding more than a bucket's
/// worth of windows) are merged, so there can be fewer than `num_bins` bins.
fn gc_bin_cuts(hist: &[u32; 1001], num_bins: usize) -> Vec<u32> {
    let total: u64 = hist.iter().map(|&c| c as u64).sum();
    assert!(total > 0, "no sequence windows with A/C/G/T to define GC bins");
    // Perl iterates up to the highest GC value observed.
    let max_v = hist.iter().rposition(|&c| c > 0).unwrap_or(0);
    let size = total as f64 / num_bins as f64;
    let mut cuts: Vec<u32> = Vec::new();
    let (mut cum, mut prev_bucket) = (0u64, -1i64);
    for (v, &c) in hist.iter().enumerate().take(max_v + 1) {
        cum += c as u64;
        let bucket = (cum as f64 / size) as i64;
        if bucket > prev_bucket {
            cuts.push(if v < max_v { v as u32 } else { 1000 });
            prev_bucket = bucket;
        }
    }
    // Perl leaves the windows above the last cut unbinned if float rounding keeps
    // the final bucket from increasing; close the range at 1000 instead.
    if *cuts.last().unwrap() != 1000 {
        cuts.push(1000);
    }
    // Perl yields no bins at all when every window is below 0.1% GC (the only
    // cut is 1000); use a single 0-1000 bin.
    if cuts[0] != 0 {
        cuts.insert(0, 0);
    }
    cuts
}

/// GC bin of a smoothed GC fraction: rounded to 0.1% (as Perl `classGC`'s
/// `sprintf "%0.3f"`), then the last bin whose start is at or below it.  A value
/// equal to a bin start belongs to that bin, matching `classGC`.
fn gc_bin_of(gc: f64, bin_starts: &[u32]) -> usize {
    let permille = (format!("{:.3}", gc).parse::<f64>().unwrap() * 1000.0).round() as u32;
    bin_starts.partition_point(|&s| s <= permille).saturating_sub(1)
}

/// Equal-population GC bins (the default).
///
/// Tries `num_bins`, then fewer, until every bin is non-empty: with target size
/// `T = ceil(total / ts)`, GC value `v` goes to bin `ceil(cumulative[v] / T) - 1`.
/// A GC value holding more than one bin's worth of windows therefore lowers the
/// bin count rather than producing one oversized bin, keeping populations even.
///
/// Returns each bin's start (lowest observed GC, per mille; used for assignment)
/// and its displayed range, whose boundaries are the midpoints between adjacent
/// bins' observed GC values.
fn equal_population_bins(hist: &[u32; 1001], num_bins: usize) -> (Vec<u32>, Vec<GCRange>) {
    let mut cum = [0u64; 1001];
    let mut acc = 0u64;
    for v in 0..=1000 {
        acc += hist[v] as u64;
        cum[v] = acc;
    }
    let total = cum[1000];
    assert!(total > 0, "no sequence windows with A/C/G/T to define GC bins");
    for ts in (1..=num_bins).rev() {
        let target = (total as f64 / ts as f64).ceil();
        let mut starts: Vec<Option<u32>> = vec![None; ts];
        let mut max_seen = vec![0u32; ts];
        let mut counts = vec![0u64; ts];
        for v in 0..=1000usize {
            if hist[v] == 0 {
                continue;
            }
            let bin = (cum[v] as f64 / target).ceil() as usize - 1;
            starts[bin].get_or_insert(v as u32);
            max_seen[bin] = v as u32;
            counts[bin] += hist[v] as u64;
        }
        if counts.iter().any(|&c| c == 0) {
            continue;
        }
        let starts: Vec<u32> = starts.into_iter().map(|s| s.unwrap()).collect();
        let mut ranges: Vec<GCRange> = Vec::with_capacity(ts);
        for k in 0..ts {
            let start = if k == 0 { starts[0] } else { (max_seen[k - 1] + starts[k]) / 2 };
            let end = if k + 1 == ts { 1000 } else { (max_seen[k] + starts[k + 1]) / 2 };
            ranges.push(GCRange { start, end, count: 0 });
        }
        return (starts, ranges);
    }
    unreachable!("a single bin always has windows")
}

/// Configuration for `SeqModel::emit()` controlling repeat bombardment.
pub struct EmitConfig<'a> {
    /// Map from TE consensus name → DNA bytes (from a RepeatMaskerLib FASTA).
    /// Required to enable TE insertion; pass `None` to skip bombardment.
    pub consensus: Option<&'a std::collections::HashMap<String, Vec<u8>>>,
    /// Target TE fraction (0–100 %).  `None` → random 10–60 % per sequence.
    pub te_frac: Option<f64>,
    /// Target simple-repeat fraction (0–100 %).  `None` → random 0–2 %.
    pub sim_frac: Option<f64>,
    /// Lowercase inserted repeat sequences (Perl default).
    /// Set to `false` to keep insertions uppercase.
    pub soft_mask: bool,
    /// Reproduce GARLIC 1.4 behaviour that is a known Perl bug (see
    /// `evolve.rs`: `add_insertions` placement, over-long TE fragments).
    /// Off by default.
    pub perl_compat: bool,
}

/// Configuration passed to `SeqModel::model()` controlling masking and
/// repeat profiling.  All fields are optional — omit the whole struct (pass
/// `None`) to get the same behaviour as the original unmasked model build.
pub struct ModelConfig<'a> {
    /// Annotation-based mask applied during k-mer counting (Pass 2).
    pub mask: Option<&'a MaskSet>,
    /// 12-column TE annotation TSV for TE repeat profiling.
    pub te_annot_file: Option<&'a std::path::Path>,
    /// UCSC TRF BED file for simple/tandem repeat profiling.
    pub tandem_annot_file: Option<&'a std::path::Path>,
    /// Gene intervals (body or exon) used as a profiling filter:
    /// TEs / tandem repeats overlapping these regions are excluded from profiles.
    pub gene_filter: Option<&'a crate::io::IntervalMap>,
    /// Window size in bp for GC binning (default: 1000).
    pub win_size: usize,
    /// K-mer size for base frequency modelling (default: 4).
    pub kmer_size: usize,
    /// Maximum number of GC bins to create (default: 25).
    pub num_gc_bins: usize,
    /// Drop tandem repeats fully contained in a TE interval when profiling
    /// simple repeats.  Off by default (Perl's `profileTRF` keeps them).
    pub filter_te_tandems: bool,
    /// Reproduce GARLIC 1.4 GC binning (`calcBinGC`/`classGC`) and TE bin
    /// positions exactly, for comparison with the Perl version.  Off by default.
    pub perl_compat: bool,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
struct GCRange {
    start: u32,
    end: u32,
    count: u32,
}

/// Model file format version written by `build` and required by `emit`.
/// Bump whenever a change makes older model files unusable or wrong.
///   2 (rusted-garlic 0.2.0): masking-drift fix, raw GC transition counts,
///     per-bin k-mer totals, simple-repeat copy number.
///   0: models without a version field (rusted-garlic 0.1.x).
pub const MODEL_FORMAT_VERSION: u32 = 2;

#[derive(Serialize, Deserialize)]
pub struct SeqModel {
    #[serde(default)]
    model_version: u32,
    a_count: usize,
    c_count: usize,
    g_count: usize,
    t_count: usize,
    n_count: usize,
    source_sequence_size: usize,
    pub win_size: usize,
    pub kmer_size: usize,
    gc_ranges: Vec<GCRange>,
    gc_trans: Vec<Vec<f64>>,
    /// Raw GC-bin transition counts (before normalization).  Used by `emit` to
    /// apply Perl's +1 pseudocount per cell (`loadGCt`).
    gc_trans_counts: Vec<Vec<f64>>,
    /// Per-bin k-mer totals (+1 per possible k-mer, as Perl's `loadKmers`);
    /// weights for drawing a random GC bin (Perl's `newGC()`).
    gc_kmer_totals: Vec<f64>,
    /// BTreeMap so the serialized model is byte-identical across runs.
    gc_binned_kmer_flank_freqs: Vec<BTreeMap<String, [f64; 4]>>,
    pub repeat_model: Option<RepeatModel>,
}

impl SeqModel {
    pub fn new() -> Self {
        SeqModel {
            model_version: MODEL_FORMAT_VERSION,
            a_count: 0,
            c_count: 0,
            g_count: 0,
            t_count: 0,
            n_count: 0,
            source_sequence_size: 0,
            win_size: 1000,
            kmer_size: 4,
            gc_ranges: Vec::new(),
            gc_trans: Vec::new(),
            gc_trans_counts: Vec::new(),
            gc_kmer_totals: Vec::new(),
            gc_binned_kmer_flank_freqs: Vec::new(),
            repeat_model: None,
        }
    }

    /// Draw a random GC bin weighted by its k-mer total, like Perl's `newGC()`.
    pub fn sample_gc_bin(&self, rng: &mut impl Rng) -> usize {
        WeightedIndex::new(&self.gc_kmer_totals).unwrap().sample(rng)
    }

    /// Look up the conditional frequencies `[P(A|pre), P(C|pre), P(G|pre), P(T|pre)]`
    /// for a k-mer prefix string in the given GC bin.
    /// Returns `None` when the prefix is absent in that bin.
    pub fn gc_kmer_freq(&self, gc_bin: usize, prefix: &str) -> Option<&[f64; 4]> {
        self.gc_binned_kmer_flank_freqs
            .get(gc_bin)
            .and_then(|m| m.get(prefix))
    }

    /// Return GC range strings (e.g. `["0.0-35.6", "35.7-66.6", ...]`) for
    /// the active bins — used when writing Perl-compatible repeat files.
    pub fn gc_range_strings(&self) -> Vec<String> {
        let mut out = Vec::new();
        for r in &self.gc_ranges {
            if r.count == 0 { break; }
            out.push(format!("{:.1}-{:.1}", r.start as f64 / 10.0, r.end as f64 / 10.0));
        }
        out
    }

    pub fn from_file(path: impl AsRef<Path>) -> Result<SeqModel, Box<dyn Error>> {
        // Open the json file
        let mut file_content = match File::open(path) {
            Ok(file) => file,
            Err(_) => panic!("Could not read the json file"),
        };

        // Transform content of the file into a string
        let mut contents = String::new();
        match file_content.read_to_string(&mut contents) {
            Ok(_) => {}
            Err(err) => panic!("Could not deserialize the file, error code: {}", err),
        };

        SeqModel::from_json(&contents)
    }

    /// Parse a model from JSON, refusing files from other model format versions.
    /// The version is checked before full deserialization so that older models
    /// get a clear "rebuild" message rather than a missing-field error.
    pub fn from_json(contents: &str) -> Result<SeqModel, Box<dyn Error>> {
        let value: serde_json::Value = serde_json::from_str(contents)?;
        let version = value.get("model_version").and_then(|v| v.as_u64()).unwrap_or(0);
        if version != MODEL_FORMAT_VERSION as u64 {
            return Err(format!(
                "model format version {} is not supported (expected {}); \
                 rebuild the model with this version of rusted-garlic",
                version, MODEL_FORMAT_VERSION
            )
            .into());
        }
        Ok(serde_json::from_value(value)?)
    }

    pub fn to_file(&self, path: impl AsRef<Path>) {
        // Open the json file
        let mut f = File::create(path).expect("Unable to create file");
        let json = serde_json::to_string(&self).expect("Serialization failed");
        write!(f, "{}", json).expect("Unable to write model file");
    }

    pub fn print_json(&self) {
        let gc_ranges: Vec<serde_json::Value> = self
            .gc_ranges
            .iter()
            .enumerate()
            .map(|(idx, gc_range)| {
                json!({
                    "id": format!("gc_range-{}", idx),
                    "start": gc_range.start,
                    "end": gc_range.end,
                    "count": gc_range.count,
                })
            })
            .collect();
        println!(
            "GC_RANGES: {}\n\n",
            &serde_json::to_string(&gc_ranges).expect("failed to serialize JSON data")
        );

        let json1 = serde_json::to_string(&self.gc_trans).expect("failed to serialize JSON data");
        println!("GC_TRANS: {}\n\n", json1);
        let json2 = serde_json::to_string(&self.gc_binned_kmer_flank_freqs)
            .expect("failed to serialize JSON data");
        println!("GC_BINNED_KMER_FLANK_FREQS: {}", json2);
    }

    /// Emit an artificial sequence.
    ///
    /// When `emit_cfg` is provided and the model contains a repeat profile,
    /// interspersed and simple repeats are inserted after the base sequence is
    /// generated (inserted regions are lower-cased).
    ///
    /// All randomness is drawn from `rng`, so a seeded RNG gives reproducible output.
    pub fn emit(
        &self,
        seq_size: Option<usize>,
        emit_cfg: Option<&EmitConfig>,
        rng: &mut impl Rng,
    ) -> (String, Vec<InsertRecord>) {
        // Build weighted samplers from the stored kmer frequencies.
        let mut gc_kmer_prefix_weights: Vec<HashMap<String, WeightedIndex<f64>>> = Vec::new();
        for _v in self.gc_ranges.iter() {
            gc_kmer_prefix_weights.push(HashMap::new());
        }
        for (gc_idx, gc_bin) in self.gc_binned_kmer_flank_freqs.iter().enumerate() {
            for (kmer_prefix, kmer_suffix_freqs) in gc_bin.iter() {
                gc_kmer_prefix_weights[gc_idx].insert(
                    kmer_prefix.clone(),
                    WeightedIndex::new(kmer_suffix_freqs).unwrap(),
                );
            }
        }

        // GC-bin transitions: like Perl's loadGCt, add one pseudocount to every
        // cell so unobserved transitions remain possible.
        let mut gc_to_gc_weights: Vec<WeightedIndex<f64>> = Vec::new();
        for row in self.gc_trans_counts.iter() {
            let w: Vec<f64> = row.iter().map(|c| c + 1.0).collect();
            gc_to_gc_weights.push(WeightedIndex::new(&w).unwrap());
        }

        // Sample initial GC bin weighted by bin window-count (Perl's newGC()).
        let gc_bin_counts: Vec<u32> = self.gc_ranges.iter().map(|r| r.count).collect();
        let gc_bin_dist = WeightedIndex::new(&gc_bin_counts).unwrap();
        let mut current_gc_bin = gc_bin_dist.sample(rng);

        let bp_to_gen: usize = seq_size.unwrap_or(self.source_sequence_size);
        let mut out_bytes: Vec<u8> = Vec::with_capacity(bp_to_gen + 1);

        // Per-window GC bin assignments: gc_bins[pos / WIN_SIZE] = gc_bin used there.
        let mut gc_bins: Vec<u32> = Vec::new();
        gc_bins.push(current_gc_bin as u32); // window 0

        // Seed with a random (k-1)-mer prefix from the starting bin (Perl's $fseed).
        // Sorted: HashMap order varies between runs and would defeat a fixed seed.
        let mut prefixes: Vec<&String> = gc_kmer_prefix_weights[current_gc_bin].keys().collect();
        prefixes.sort();
        let seed_prefix = prefixes[rng.gen_range(0..prefixes.len())].clone();
        out_bytes.extend_from_slice(seed_prefix.as_bytes());

        let perl_compat = emit_cfg.map_or(false, |c| c.perl_compat);
        if perl_compat {
            // GARLIC 1.4 createSeq/createSubSeq, bug for bug: each window starts
            // from the previous window's last k-1 bases and appends that seed
            // *plus* win-2k+3 new bases, i.e. win-k+2 bp per window (998 for
            // k=4) with the seed duplicated.  Windows therefore drift against
            // pos/win_size, and checkSeqSize pads the shortfall by copying a
            // random internal segment onto the end.
            gc_bins.clear();
            let k1 = self.kmer_size - 1;
            let w = self.win_size - self.kmer_size + 1;
            let mut i = out_bytes.len();
            while i <= bp_to_gen + 100 {
                gc_bins.push(current_gc_bin as u32);
                let mut s: Vec<u8> = out_bytes[out_bytes.len() - k1..].to_vec();
                while s.len() <= w {
                    let km = std::str::from_utf8(&s[s.len() - k1..]).unwrap();
                    let nextbase = gc_kmer_prefix_weights[current_gc_bin]
                        .get(km)
                        .expect("missing kmer prefix in model")
                        .sample(rng);
                    s.push(b"ACGT"[nextbase]);
                }
                out_bytes.extend_from_slice(&s);
                current_gc_bin = gc_to_gc_weights[current_gc_bin].sample(rng);
                i += self.win_size;
            }
            if out_bytes.len() > bp_to_gen {
                out_bytes.truncate(bp_to_gen);
            } else if out_bytes.len() < bp_to_gen {
                let add = bp_to_gen - out_bytes.len();
                let pos = rng.gen_range(0..out_bytes.len().saturating_sub(add).max(1));
                let patch: Vec<u8> = out_bytes[pos..(pos + add).min(out_bytes.len())].to_vec();
                out_bytes.extend(patch);
            }
        } else {
            let mut km = seed_prefix;

            while out_bytes.len() < bp_to_gen {
                let nextbase = gc_kmer_prefix_weights[current_gc_bin]
                    .get(&km)
                    .expect("missing kmer prefix in model")
                    .sample(rng);
                let b: u8 = match nextbase {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => continue,
                };
                out_bytes.push(b);
                km.push(b as char);
                km.remove(0);
                // Transition GC bin every win_size bases (Perl's createSeq outer loop).
                if out_bytes.len() % self.win_size == 0 {
                    current_gc_bin = gc_to_gc_weights[current_gc_bin].sample(rng);
                    gc_bins.push(current_gc_bin as u32);
                }
            }
        }

        // Repeat bombardment — requires a repeat model.
        // TE insertion additionally requires consensus sequences.
        // Simple repeat insertion only needs the model's simple_profile.
        let mut all_records: Vec<InsertRecord> = Vec::new();
        if let (Some(cfg), Some(rm)) = (emit_cfg, &self.repeat_model) {
            let soft = cfg.soft_mask;
            let te_total: usize = rm.te_profile.iter().map(|v| v.len()).sum();
            let sim_total: usize = rm.simple_profile.iter().map(|v| v.len()).sum();
            if let Some(consensus) = cfg.consensus {
                if te_total > 0 {
                    let recs = insert_te_repeats(
                        &mut out_bytes,
                        &gc_bins,
                        rm,
                        consensus,
                        self,
                        cfg.te_frac,
                        soft,
                        rng,
                        cfg.perl_compat,
                    );
                    all_records.extend(recs);
                } else {
                    log::warn!("--te_cons provided but the model has no TE repeat records; \
                                rebuild with --te_annot <tsv_file> to enable TE insertion");
                }
            }
            if sim_total > 0 {
                let recs = insert_simple_repeats(
                    &mut out_bytes,
                    &gc_bins,
                    rm,
                    self,
                    cfg.sim_frac,
                    soft,
                    rng,
                    cfg.perl_compat,
                );
                all_records.extend(recs);
            } else if cfg.sim_frac.is_some() {
                log::warn!("--sim_frac provided but the model has no simple repeat records; \
                            rebuild with --tandem_annot <tsv_file> to enable simple repeat insertion");
            }
        }

        (String::from_utf8(out_bytes).expect("sequence is not valid UTF-8"), all_records)
    }

    /// Build a sequence model from a FASTA file.
    ///
    /// Pass 1 computes GC bins from the *raw* sequence (unmasked), matching the
    /// Perl default where `calcBinGC()` runs before masking.
    ///
    /// Pass 2 counts k-mers only over unmasked [ACGT] positions.  If `mask` is
    /// provided, any position falling inside a mask interval is treated as N
    /// (kmer buffer cleared) so it does not contribute to the frequency model.
    pub fn model(path: impl AsRef<Path>, config: Option<&ModelConfig>) -> SeqModel {
        let mut seq_model = SeqModel::new();
        // Calculate the running smoothed average GC level for in 1kb windows
        //
        // smoothed_avg[0] = avg( win_avg[0] .. win_avg[5] )
        //
        // Convert non ACGTacgt to N
        //
        // Skip any kmer containing an N
        // Add count to GC bin (1000bp non-overlapping windows)
        //            Uses bingc ( 10kb averages of 1kb windows )
        // GC:
        //    |----1000bp----|
        //                   |-------1000bp----|
        //                                     |...
        //
        //      gc = count_gc / count_acgt
        //      count_acgt != length(seq) if there are Ns
        //      gcHash{seq_id} = [bin1, bin2, ...]
        //      count[ 1000*gc ]++  (1000 slots)
        //      cumulative[ 1 ] = cumulative[0] + count[1]; (1000 slots)
        //      total = total number of bins
        //      gcBuckets = 25
        //      gcBins[ [i, cumulative[i]]
        //      gc[] = 25 gc ranges that evenly divide the genome
        //      blk = 10000/1000 / 2 = 5
        //      bingc[
        //      |----1kb---|----1kb---|...|----1kb---|
        //      calculate average GC in -5/+5 around each 1kb window
        //
        let print_garlic = true;
        let win_size = config.map_or(1000, |c| c.win_size);
        let kmer_size = config.map_or(4, |c| c.kmer_size);
        let max_gc_bins = config.map_or(25, |c| c.num_gc_bins);
        seq_model.win_size = win_size;
        seq_model.kmer_size = kmer_size;
        let window_size = win_size;
        let (raw, _) = crate::io::open_raw(path.as_ref());
        let mut buf = BufReader::new(raw);

        let mut bytes = [0; 512];
        let mut in_header = false;
        let mut win_at_count = 0;
        let mut win_cg_count = 0;
        let mut win_bp_count = 0;

        let perl_compat = config.map_or(false, |c| c.perl_compat);

        // Pass 1: raw GC of every window at full precision, a histogram of window
        // GC in 0.1% slots (full windows only), then equal-population GC bins and
        // smoothed per-window bin assignments (below).
        //
        // `perl_compat` follows GARLIC 1.4 calcBinGC/classGC exactly: histogram
        // slot int(1000*gc) instead of rounding, and a trailing partial window is
        // ignored with the chromosome's last full-window GC repeated in its place.
        let mut gc_hist = [0u32; 1001];
        // Per chromosome: name, raw GC per window (None = no A/C/G/T), and the
        // number of full windows (a trailing entry may follow for the tail).
        let mut chrom_raw_gc: Vec<(String, Vec<Option<f64>>, usize)> = Vec::new();
        let mut cur_raw_gc: Vec<Option<f64>> = Vec::new();

        let mut p1_chr_name = String::new();
        let mut p1_in_header_name = false;

        fn window_gc(cg: usize, at: usize) -> Option<f64> {
            if cg + at > 0 { Some(cg as f64 / (cg + at) as f64) } else { None }
        }
        // Record a full window's raw GC.
        fn push_window(cg: usize, at: usize, raw: &mut Vec<Option<f64>>, hist: &mut [u32; 1001],
                       perl_compat: bool) {
            let gc = window_gc(cg, at);
            if let Some(g) = gc {
                let slot = if perl_compat { (1000.0 * g) as usize } else { (1000.0 * g).round() as usize };
                hist[slot] += 1;
            }
            raw.push(gc);
        }
        // Close a chromosome.  The trailing partial window (if any) gets its own
        // GC; Perl instead repeats the last full window's GC, partial or not.
        fn finish_chrom(name: &str, raw: &mut Vec<Option<f64>>, partial: Option<Option<f64>>,
                        perl_compat: bool, out: &mut Vec<(String, Vec<Option<f64>>, usize)>) {
            let n_full = raw.len();
            if perl_compat {
                let last = raw.last().copied().flatten();
                raw.push(last);
            } else if let Some(gc) = partial {
                raw.push(gc);
            }
            out.push((name.to_string(), std::mem::take(raw), n_full));
        }

        loop {
            match buf.read(&mut bytes) {
                // EOF Reached
                Ok(0) => {
                    if win_bp_count == window_size {
                        push_window(win_cg_count, win_at_count, &mut cur_raw_gc, &mut gc_hist, perl_compat);
                        win_bp_count = 0;
                    }
                    if !p1_chr_name.is_empty() || !cur_raw_gc.is_empty() {
                        let partial = (win_bp_count > 0).then(|| window_gc(win_cg_count, win_at_count));
                        finish_chrom(&p1_chr_name, &mut cur_raw_gc, partial, perl_compat, &mut chrom_raw_gc);
                    }
                    break;
                }
                Ok(n) => {
                    for i in 0..n {
                        if win_bp_count == window_size {
                            push_window(win_cg_count, win_at_count, &mut cur_raw_gc, &mut gc_hist, perl_compat);
                            win_at_count = 0;
                            win_cg_count = 0;
                            win_bp_count = 0;
                        }

                        if bytes[i] == '>' as u8 {
                            if !p1_chr_name.is_empty() || !cur_raw_gc.is_empty() {
                                let partial = (win_bp_count > 0).then(|| window_gc(win_cg_count, win_at_count));
                                finish_chrom(&p1_chr_name, &mut cur_raw_gc, partial, perl_compat, &mut chrom_raw_gc);
                            }
                            in_header = true;
                            p1_in_header_name = true;
                            p1_chr_name.clear();
                            win_cg_count = 0;
                            win_at_count = 0;
                            win_bp_count = 0;
                            continue;
                        }
                        if bytes[i] == '\n' as u8 {
                            if in_header { p1_in_header_name = false; }
                            in_header = false;
                            continue;
                        }
                        if in_header {
                            if p1_in_header_name {
                                if bytes[i] == b' ' || bytes[i] == b'\t' {
                                    p1_in_header_name = false;
                                } else {
                                    p1_chr_name.push(bytes[i] as char);
                                }
                            }
                            continue;
                        }
                        // Sequence body (not in header — all header cases continued above)
                        match bytes[i] as char {
                                'A' => {
                                    seq_model.a_count += 1;
                                    win_at_count += 1;
                                    win_bp_count += 1;
                                }
                                'T' => {
                                    seq_model.t_count += 1;
                                    win_at_count += 1;
                                    win_bp_count += 1;
                                }
                                'C' => {
                                    seq_model.c_count += 1;
                                    win_cg_count += 1;
                                    win_bp_count += 1;
                                }
                                'G' => {
                                    seq_model.g_count += 1;
                                    win_cg_count += 1;
                                    win_bp_count += 1;
                                }
                                'a' => {
                                    seq_model.a_count += 1;
                                    win_at_count += 1;
                                    win_bp_count += 1;
                                }
                                't' => {
                                    seq_model.t_count += 1;
                                    win_at_count += 1;
                                    win_bp_count += 1;
                                }
                                'c' => {
                                    seq_model.c_count += 1;
                                    win_cg_count += 1;
                                    win_bp_count += 1;
                                }
                                'g' => {
                                    seq_model.g_count += 1;
                                    win_cg_count += 1;
                                    win_bp_count += 1;
                                }
                                // IUB Codes (all translate to N's)
                                'n' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'N' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'R' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'Y' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'K' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'M' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'S' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'W' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'B' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'D' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'H' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'V' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'r' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'y' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'k' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'm' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                's' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'w' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'b' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'd' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'h' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                'v' => {
                                    seq_model.n_count += 1;
                                    win_bp_count += 1;
                                }
                                _other => (),
                            }
                    }
                }
                Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                Err(e) => panic!("{:?}", e),
            };
        }
        seq_model.source_sequence_size = seq_model.a_count
            + seq_model.c_count
            + seq_model.g_count
            + seq_model.t_count
            + seq_model.n_count;
        if print_garlic {
            println!("Sequence Stats");
            println!("    A : {}", seq_model.a_count);
            println!("    C : {}", seq_model.c_count);
            println!("    G : {}", seq_model.g_count);
            println!("    T : {}", seq_model.t_count);
            println!("    N : {}", seq_model.n_count);
            println!(" Total: {}", seq_model.source_sequence_size);
        }

        // GC bins: `bin_starts[k]` is the lowest GC (per mille) assigned to bin k.
        let (bin_starts, gc_ranges) = if perl_compat {
            let cuts = gc_bin_cuts(&gc_hist, max_gc_bins);
            let ranges = cuts.windows(2).map(|c| GCRange { start: c[0], end: c[1], count: 0 }).collect();
            (cuts[..cuts.len() - 1].to_vec(), ranges)
        } else {
            equal_population_bins(&gc_hist, max_gc_bins)
        };
        let actual_bins = bin_starts.len();
        let sentinel_bin = actual_bins as u32; // "no bin" (window without A/C/G/T)
        seq_model.gc_ranges = gc_ranges;
        seq_model.gc_trans = vec![vec![0.0f64; actual_bins]; actual_bins];

        // Smooth each window's raw GC over +/- flank windows (full precision,
        // rounded once to 0.1%), assign bins, and count transitions between
        // consecutive full windows.  Perl's flank is int((10 kb / win) / 2).
        let flank = if perl_compat { (10_000 / win_size) / 2 } else { 5 };
        // chrom_gc_bins[chrom][w] = smoothed GC bin of window w (sentinel_bin if none),
        // including an entry for a trailing partial window.
        let mut chrom_gc_bins: HashMap<String, Vec<u32>> = HashMap::new();
        for (name, raw, n_full) in &chrom_raw_gc {
            let n = raw.len();
            let mut bins: Vec<u32> = Vec::with_capacity(n);
            for i in 0..n {
                let (lo, hi) = (i.saturating_sub(flank), cmp::min(i + flank, n - 1));
                let (mut sum, mut num) = (0.0f64, 0usize);
                for g in raw[lo..=hi].iter().flatten() {
                    sum += g;
                    num += 1;
                }
                let bin = if num > 0 { gc_bin_of(sum / num as f64, &bin_starts) as u32 } else { sentinel_bin };
                bins.push(bin);
            }
            let full = &bins[..*n_full];
            for &b in full {
                if b < sentinel_bin {
                    seq_model.gc_ranges[b as usize].count += 1;
                }
            }
            for w in full.windows(2) {
                if w[0] < sentinel_bin && w[1] < sentinel_bin {
                    seq_model.gc_trans[w[0] as usize][w[1] as usize] += 1.0;
                }
            }
            chrom_gc_bins.insert(name.clone(), bins);
        }

        // Print a summary of GC bins with their window counts and approximate base-pair
        // coverage.  This appears before Pass 2 so the user can sanity-check the input
        // sequence composition before the heavier k-mer counting begins.
        {
            let nbins = seq_model.gc_ranges.len();
            println!("\nGC Bin Summary ({} bin{}, win_size={}):",
                nbins, if nbins == 1 { "" } else { "s" }, win_size);
            println!("  {:>3}  {:>16}  {:>9}  {:>12}",
                "Bin", "GC Range", "Windows", "~bp");
            println!("  {}  {}  {}  {}",
                "-".repeat(3), "-".repeat(16), "-".repeat(9), "-".repeat(12));
            for v in 0..nbins {
                let r = &seq_model.gc_ranges[v];
                println!("  {:>3}  {:>6.1}% - {:>5.1}%  {:>9}  {:>12}",
                    v,
                    r.start as f64 / 10.0,
                    r.end   as f64 / 10.0,
                    r.count,
                    r.count as u64 * win_size as u64);
            }
            println!();
        }

        // Initialize a hashes to hold per-gc-bin kmer counts
        let mut maps = Vec::new();
        for _v in seq_model.gc_ranges.iter() {
            let hm: HashMap<String, u64> = HashMap::new();
            maps.push(hm);
        }

        // Pass 2: k-mer counting
        //
        // Rewind to beginning of FASTA file and count k-mers per GC bin,
        // skipping any position that is masked (treated as N).
        //
        // Masking is applied only in this pass — Pass 1 (GC binning) uses the
        // raw sequence, matching the Perl default where calcBinGC() runs before
        // masking.
        //
        // To support per-position masking we track:
        //   bp_count     — total bases seen across all sequences (maps into bingc[])
        //   chr_bp_count — 0-based position within the current sequence (for mask lookup)
        //   chr_name     — current sequence id (first word of FASTA header)
        // Reopen the file for Pass 2 — seeking is not possible when the reader
        // wraps a GzDecoder, so we simply open a fresh reader.
        let (raw2, _) = crate::io::open_raw(path.as_ref());
        buf = BufReader::new(raw2);
        in_header = false;
        let kmer_len = kmer_size;
        let mut kmer: VecDeque<char> = VecDeque::with_capacity(kmer_size);
        let mut chr_bp_count: u32 = 0;
        let mut chr_name = String::new();
        let mut in_header_name = false; // true until first whitespace after '>'

        // mask_cursor is refreshed per chromosome as headers are encountered.
        // When config has no mask it stays None and masking is skipped entirely.
        let mask = config.and_then(|c| c.mask);
        let mut mask_cursor: Option<crate::masking::MaskCursor> = None;

        // GC-bin lookup helper: uses the per-chromosome smoothed bins built during
        // the smoothing pass.  This avoids indexing into bingc with a global base
        // counter — that approach accumulated a deficit of one slot per chromosome
        // and panicked on genomes with many sequences.
        let gc_bin_for = |chr: &str, chr_pos: u32| -> usize {
            let win = chr_pos as usize / win_size;
            chrom_gc_bins
                .get(chr)
                .and_then(|v| v.get(win))
                .copied()
                .unwrap_or(sentinel_bin) as usize
        };

        loop {
            match buf.read(&mut bytes) {
                Ok(0) => {
                    // EOF: flush the last kmer if one is pending.
                    if kmer.len() == kmer_len && chr_bp_count >= kmer_len as u32 {
                        let gc_bin = gc_bin_for(&chr_name, chr_bp_count - kmer_len as u32);
                        if gc_bin < actual_bins {
                            let s: String = kmer.iter().collect();
                            *maps[gc_bin].entry(s).or_default() += 1;
                        }
                        kmer.pop_front();
                    }
                    break;
                }
                Ok(n) => {
                    for i in 0..n {
                        // Count the full k-mer that just became available.
                        // chr_bp_count reflects all bases pushed so far in this chromosome;
                        // the most-recently-added base is at position chr_bp_count-1.
                        if kmer.len() == kmer_len && chr_bp_count >= kmer_len as u32 {
                            let gc_bin = gc_bin_for(&chr_name, chr_bp_count - kmer_len as u32);
                            if gc_bin < actual_bins {
                                let s: String = kmer.iter().collect();
                                *maps[gc_bin].entry(s).or_default() += 1;
                            }
                            kmer.pop_front();
                        }

                        let ch = bytes[i];

                        if ch == b'>' {
                            // New sequence: discard any partial kmer, reset per-chr counters.
                            kmer.clear();
                            chr_bp_count = 0;
                            chr_name.clear();
                            in_header = true;
                            in_header_name = true;
                            mask_cursor = None;
                            continue;
                        }

                        if in_header {
                            if ch == b'\n' {
                                in_header = false;
                                in_header_name = false;
                                // Now that we have the full chr_name, get a cursor
                                mask_cursor = mask.map(|m| m.cursor_for(&chr_name));
                            } else if in_header_name {
                                if ch == b' ' || ch == b'\t' {
                                    in_header_name = false; // stop at first whitespace
                                } else {
                                    chr_name.push(ch as char);
                                }
                            }
                            continue;
                        }

                        // Skip line breaks and other non-sequence bytes *before* the mask
                        // check: counting them as positions in masked regions made
                        // chr_bp_count drift ahead (one per FASTA line inside a mask),
                        // shifting both mask lookups and GC-bin assignment downstream.
                        // The accepted set matches the bytes Pass 1 counts as positions.
                        if !b"ACGTNRYKMSWBDHVacgtnrykmswbdhv".contains(&ch) {
                            continue;
                        }

                        // Sequence data — check mask before deciding how to handle the base
                        let masked = mask_cursor
                            .as_mut()
                            .map_or(false, |c| c.is_masked(chr_bp_count));

                        if masked {
                            kmer.clear();
                            chr_bp_count += 1;
                            continue;
                        }

                        match ch as char {
                            'A' | 'a' => { kmer.push_back('A'); chr_bp_count += 1; }
                            'C' | 'c' => { kmer.push_back('C'); chr_bp_count += 1; }
                            'G' | 'g' => { kmer.push_back('G'); chr_bp_count += 1; }
                            'T' | 't' => { kmer.push_back('T'); chr_bp_count += 1; }
                            // N and all IUB ambiguity codes break the k-mer context
                            'N' | 'n' | 'R' | 'r' | 'Y' | 'y' | 'K' | 'k' |
                            'M' | 'm' | 'S' | 's' | 'W' | 'w' | 'B' | 'b' |
                            'D' | 'd' | 'H' | 'h' | 'V' | 'v' => {
                                kmer.clear();
                                chr_bp_count += 1;
                            }
                            _ => (),
                        }
                    }
                }
                Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                Err(e) => panic!("{:?}", e),
            };
        }

        // Pass 3: Repeat profiling
        // Profile TE and simple repeats from annotation files (if provided).
        // Uses chrom_gc_bins built during Pass 1 smoothing.
        if let Some(cfg) = config {
            if cfg.te_annot_file.is_some() || cfg.tandem_annot_file.is_some() {
                let mut rm = RepeatModel::new(actual_bins);
                let te_ivs = if let Some(rpath) = cfg.te_annot_file {
                    log::info!("  Profiling TE repeats from {:?}...", rpath);
                    let te_ivs = profile_repeats(rpath, &chrom_gc_bins, win_size as u32,
                                                 actual_bins, &mut rm, cfg.gene_filter,
                                                 cfg.perl_compat);
                    log::info!("  {} TE records profiled", rm.te_count());
                    te_ivs
                } else {
                    crate::io::IntervalMap::new()
                };
                if let Some(tpath) = cfg.tandem_annot_file {
                    let te_opt = if cfg.filter_te_tandems && !te_ivs.is_empty() {
                        Some(&te_ivs)
                    } else {
                        None
                    };
                    log::info!("  Profiling tandem repeats from {:?}...", tpath);
                    profile_trf(tpath, &chrom_gc_bins, win_size as u32, &mut rm,
                                te_opt, cfg.gene_filter);
                    log::info!("  {} simple repeat records profiled", rm.simple_count());
                }
                seq_model.repeat_model = Some(rm);
            }
        }

        // Need:  tot_sum : total number of words counted in the input sequence (over all GC bins)
        //        gc_sum{GC} : total number of words counted in a given GC bin

        // Initialize a hashes to hold per-gc-bin kmer flanking base frequencies
        for _v in seq_model.gc_ranges.iter() {
            let hm: BTreeMap<String, [f64; 4]> = BTreeMap::new();
            seq_model.gc_binned_kmer_flank_freqs.push(hm);
        }

        let n_kmers = 4f64.powi(kmer_size as i32);
        seq_model.gc_kmer_totals = maps
            .iter()
            .map(|m| m.values().sum::<u64>() as f64 + n_kmers)
            .collect();

        let bases = vec!['A', 'C', 'G', 'T'];
        let mut bin_idx = 0;
        for kmer_bin in maps.iter() {
            if print_garlic {
                println!(
                    "#GC={:.1}-{:.1}",
                    seq_model.gc_ranges[bin_idx].start as f64 / 10.0,
                    seq_model.gc_ranges[bin_idx].end as f64 / 10.0
                );
            }
            // Garlic uses the highest end of the GC range as the GC fraction for this bin.
            // Fix 4: use floating-point division; integer division truncated to 0 for all
            //        bins except 100% GC (e.g. end=400 → 400/1000=0, not 0.4).
            let p_val = seq_model.gc_ranges[bin_idx].end as f64 / 1000.0;
            let q_val = 1.0 - p_val;
            let mut bin_word_count = 0u64;
            for (_kmer, count) in kmer_bin.iter() {
                bin_word_count += count;
            }
            log::debug!("BIN_WORDS bin={} gc={:.1}-{:.1} words={}",
                bin_idx,
                seq_model.gc_ranges[bin_idx].start as f64 / 10.0,
                seq_model.gc_ranges[bin_idx].end as f64 / 10.0,
                bin_word_count);
            let mut kmer = String::with_capacity(kmer_size);
            for x in 0..(bases.len().pow((kmer_len - 1) as u32)) {
                kmer.clear();
                for y in (0..(kmer_len - 1)).rev() {
                    let twobit_base = (x & (3 << (y * 2))) >> (y * 2);
                    match twobit_base {
                        0 => kmer.push('A'),
                        1 => kmer.push('C'),
                        2 => kmer.push('G'),
                        3 => kmer.push('T'),
                        _other => (),
                    }
                }
                let mut group_total = 0;
                for b in bases.iter() {
                    kmer.push(*b);
                    match kmer_bin.get(&kmer) {
                        Some(val) => group_total += val,
                        None => (),
                    }
                    kmer.pop();
                }


                let kmer_entry = seq_model.gc_binned_kmer_flank_freqs[bin_idx]
                    .entry(kmer.clone())
                    .or_default();
                let mut bidx = 0;
                for b in bases.iter() {
                    kmer.push(*b);
                    if group_total > 0 {
                        match kmer_bin.get(&kmer) {
                            // lossy conversion from usize to f64...hmmm
                            Some(val) => {
                                let vf = *val as f64;
                                kmer_entry[bidx] = vf / (group_total as f64);
                                if print_garlic {
                                    if bin_word_count > 0 {
                                        println!(
                                            "{}\t{:.8}\t{:0.8}\t{}",
                                            kmer,
                                            (vf / (group_total as f64)),
                                            (vf / (bin_word_count as f64)),
                                            *val
                                        );
                                    } else {
                                        println!(
                                            "{}\t{:.8}\t-\t{}",
                                            kmer,
                                            (vf / (group_total as f64)),
                                            *val
                                        );
                                    }
                                }
                            }
                            None => {
                                kmer_entry[bidx] = 0.0 as f64;
                                if print_garlic {
                                    println!("{}\t{:.8}\t{:0.8}\t{}", kmer, 0, 0, 0);
                                }
                            }
                        }
                    } else {
                        // Fix 5: A/T share the AT-fraction (q_val); C/G share the GC-fraction
                        //        (p_val). Previously both branches used q_val.
                        let freq = if *b == 'A' || *b == 'T' {
                            q_val / 2.0
                        } else {
                            p_val / 2.0
                        };
                        kmer_entry[bidx] = freq;
                        if print_garlic {
                            println!("{}\t{:.8}\t{:0.8}\t{}", kmer, freq, 0, 0);
                        }
                    }
                    kmer.pop();
                    bidx += 1;
                }
            }
            bin_idx += 1;
        }

        if print_garlic {
            println!("");
            println!("GC to GC counts");
        }
        let def_frq = 1 as f64 / seq_model.gc_ranges.len() as f64;
        seq_model.gc_trans_counts = seq_model.gc_trans.clone();
        for i in 0..seq_model.gc_ranges.len() {
            let mut tot = 0.0 as f64;
            for j in 0..seq_model.gc_ranges.len() {
                tot += seq_model.gc_trans[i][j];
            }

            // Do this before we convert the counts into frequencies
            if print_garlic {
                for j in 0..seq_model.gc_ranges.len() {
                    if tot > 0.0 as f64 {
                        println!(
                            "{:.1}-{:.1}\t{:.1}-{:.1}\t{:.8}\t{}",
                            seq_model.gc_ranges[i].start as f64 / 10.0,
                            seq_model.gc_ranges[i].end as f64 / 10.0,
                            seq_model.gc_ranges[j].start as f64 / 10.0,
                            seq_model.gc_ranges[j].end as f64 / 10.0,
                            seq_model.gc_trans[i][j] as f64 / tot as f64,
                            seq_model.gc_trans[i][j]
                        );
                    } else {
                        println!(
                            "{:.1}-{:.1}\t{:.1}-{:.1}\t{:.8}\t{}",
                            seq_model.gc_ranges[i].start as f64 / 10.0,
                            seq_model.gc_ranges[i].end as f64 / 10.0,
                            seq_model.gc_ranges[j].start as f64 / 10.0,
                            seq_model.gc_ranges[j].end as f64 / 10.0,
                            def_frq,
                            0
                        );
                    }
                }
            }

            for j in 0..seq_model.gc_ranges.len() {
                if tot > 0.0 as f64 {
                    // Should there be pseudocounts for the unobserved transitions?
                    seq_model.gc_trans[i][j] = seq_model.gc_trans[i][j] / tot as f64;
                } else {
                    seq_model.gc_trans[i][j] = def_frq;
                }
            }
        }
        seq_model
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Write a FASTA file with the given (name, sequence) pairs.
    fn write_fasta(path: &std::path::Path, sequences: &[(&str, Vec<u8>)]) {
        let mut file = fs::File::create(path).expect("could not create temp FASTA");
        for (name, seq) in sequences {
            writeln!(file, ">{}", name).unwrap();
            for chunk in seq.chunks(60) {
                file.write_all(chunk).unwrap();
                writeln!(file).unwrap();
            }
        }
    }

    /// Build a sequence of `num_windows` windows, each of `win_size` bp, with
    /// exactly `gc_per_1000 * win_size / 1000` C bases and the rest A bases.
    ///
    /// Because the raw-GC calculation is:
    ///   gc = (cg_count as f64 / acgt_count as f64 * 1000.0) as u32
    ///
    /// choosing values where `gc_per_1000 * win_size` is divisible by 1000 gives
    /// an exact integer GC score — e.g. gc_per_1000=250, win_size=100 → 25 C's
    /// per window → raw GC score = 250.
    fn make_sequence(gc_per_1000: u32, num_windows: usize, win_size: usize) -> Vec<u8> {
        let gc_count = (gc_per_1000 as usize * win_size) / 1000;
        let at_count = win_size - gc_count;
        let window: Vec<u8> = vec![b'C'; gc_count]
            .into_iter()
            .chain(vec![b'A'; at_count])
            .collect();
        window.repeat(num_windows)
    }

    /// Minimal ModelConfig for bin-boundary tests (no masking, no repeat files).
    fn basic_config(win_size: usize, num_gc_bins: usize) -> ModelConfig<'static> {
        ModelConfig {
            mask: None,
            te_annot_file: None,
            tandem_annot_file: None,
            gene_filter: None,
            win_size,
            kmer_size: 2,
            num_gc_bins,
            filter_te_tandems: false,
            perl_compat: false,
        }
    }

    // -----------------------------------------------------------------------
    // GC bin-boundary validation tests
    //
    // Strategy: put each distinct GC value on its own chromosome so the ±5-window
    // smoothing never averages across GC levels.  All outputs for a uniform
    // chromosome therefore stay in the same bin, keeping expected values simple.
    // -----------------------------------------------------------------------

    /// A single chromosome where every window has 50 % GC (raw score = 500).
    /// The equal-population algorithm must fall all the way back to ts=1 → 1 bin.
    ///
    /// Expected bin:
    ///   bin 0 : start=500, end=1000 (last bin always gets end=1000)
    #[test]
    fn test_single_bin_uniform_gc() {
        let tmp = std::env::temp_dir().join("garlic_test_single_bin.fa");
        // 20 windows × 100 bp, each window 50 C + 50 A → raw GC = 500
        write_fasta(&tmp, &[("chr1", make_sequence(500, 20, 100))]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 1,
            "all windows at 50% GC should collapse to a single bin");
        assert_eq!(model.gc_ranges[0].start, 500,
            "the only bin starts at the observed GC value (500 = 50.0%)");
        assert_eq!(model.gc_ranges[0].end, 1000,
            "the last bin always ends at 1000 (100.0%)");
        assert!(model.gc_ranges[0].count > 0, "bin must contain windows");
    }

    /// Two chromosomes: one all-AT (raw GC = 0) and one all-GC (raw GC = 1000).
    ///
    /// cumulative histogram: cumul[0]=10, cumul[1..999]=10, cumul[1000]=20
    /// With ts=2, target_bin_size=10:
    ///   GC 0..999 → bin 0   (ceil(10/10)−1 = 0)
    ///   GC 1000   → bin 1   (ceil(20/10)−1 = 1)
    ///
    /// Midpoint boundary: (actual_max_bin0=0 + actual_min_bin1=1000) / 2 = 500
    ///
    /// Expected bins:
    ///   bin 0 : start=0,   end=500
    ///   bin 1 : start=500, end=1000
    #[test]
    fn test_two_bins_bimodal_gc() {
        let tmp = std::env::temp_dir().join("garlic_test_two_bins.fa");
        write_fasta(&tmp, &[
            ("chrAT", make_sequence(0,    10, 100)),   // raw GC = 0
            ("chrGC", make_sequence(1000, 10, 100)),   // raw GC = 1000
        ]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 2,
            "bimodal GC (0% and 100%) should produce exactly 2 bins");
        // Midpoint between observed GC=0 and GC=1000 is 500.
        assert_eq!(model.gc_ranges[0].start, 0,    "bin 0: start = 0 (0.0% GC)");
        assert_eq!(model.gc_ranges[0].end,   500,  "bin 0: end at midpoint 50.0%");
        assert_eq!(model.gc_ranges[1].start, 500,  "bin 1: start at midpoint 50.0%");
        assert_eq!(model.gc_ranges[1].end,   1000, "last bin end is always 1000");
        assert!(model.gc_ranges[0].count > 0, "bin 0 must contain windows");
        assert!(model.gc_ranges[1].count > 0, "bin 1 must contain windows");
    }

    /// Three chromosomes at GC = 0 %, 50 %, and 100 % (raw scores 0, 500, 1000).
    ///
    /// cumulative histogram: cumul[0]=10, cumul[1..499]=10, cumul[500]=20,
    ///                       cumul[501..999]=20, cumul[1000]=30
    /// With ts=3, target_bin_size=10:
    ///   GC 0..499   → bin 0   (ceil(10/10)−1 = 0)
    ///   GC 500..999 → bin 1   (ceil(20/10)−1 = 1)
    ///   GC 1000     → bin 2   (ceil(30/10)−1 = 2)
    ///
    /// Midpoints: (0+500)/2=250, (500+1000)/2=750
    ///
    /// Expected bins (equal population):
    ///   bin 0 : start=0,   end=250
    ///   bin 1 : start=250, end=750
    ///   bin 2 : start=750, end=1000
    #[test]
    fn test_three_bins_trimodal_gc() {
        let tmp = std::env::temp_dir().join("garlic_test_three_bins.fa");
        write_fasta(&tmp, &[
            ("chrAT",  make_sequence(0,    10, 100)),
            ("chrMid", make_sequence(500,  10, 100)),
            ("chrGC",  make_sequence(1000, 10, 100)),
        ]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 3,
            "three distinct GC levels should produce 3 bins");
        // Midpoint between observed GC=0 and GC=500 is 250; between 500 and 1000 is 750.
        assert_eq!(model.gc_ranges[0].start, 0,    "bin 0: start = 0 (0.0% GC)");
        assert_eq!(model.gc_ranges[0].end,   250,  "bin 0: end at midpoint 25.0%");
        assert_eq!(model.gc_ranges[1].start, 250,  "bin 1: start at midpoint 25.0%");
        assert_eq!(model.gc_ranges[1].end,   750,  "bin 1: end at midpoint 75.0%");
        assert_eq!(model.gc_ranges[2].start, 750,  "bin 2: start at midpoint 75.0%");
        assert_eq!(model.gc_ranges[2].end,   1000, "last bin end is always 1000");
        // All chromosomes have 10 windows → equal-population bins.
        assert_eq!(model.gc_ranges[0].count, model.gc_ranges[1].count,
            "bins 0 and 1 should be equal-population");
        assert_eq!(model.gc_ranges[1].count, model.gc_ranges[2].count,
            "bins 1 and 2 should be equal-population");
    }

    /// Five chromosomes at GC = 0 %, 25 %, 50 %, 75 %, 100 %
    /// (raw scores 0, 250, 500, 750, 1000), 10 windows each → 50 windows total.
    ///
    /// With ts=5, target_bin_size=10:
    ///   GC 0..249   → bin 0   (ceil(10/10)−1 = 0)
    ///   GC 250..499 → bin 1   (ceil(20/10)−1 = 1)
    ///   GC 500..749 → bin 2   (ceil(30/10)−1 = 2)
    ///   GC 750..999 → bin 3   (ceil(40/10)−1 = 3)
    ///   GC 1000     → bin 4   (ceil(50/10)−1 = 4)
    ///
    /// Midpoints: (0+250)/2=125, (250+500)/2=375, (500+750)/2=625, (750+1000)/2=875
    ///
    /// Expected bins (equal population):
    ///   bin 0 : start=0,   end=125
    ///   bin 1 : start=125, end=375
    ///   bin 2 : start=375, end=625
    ///   bin 3 : start=625, end=875
    ///   bin 4 : start=875, end=1000
    #[test]
    fn test_five_bins_equal_population() {
        let tmp = std::env::temp_dir().join("garlic_test_five_bins.fa");
        write_fasta(&tmp, &[
            ("chr0",    make_sequence(0,    10, 100)),
            ("chr250",  make_sequence(250,  10, 100)),
            ("chr500",  make_sequence(500,  10, 100)),
            ("chr750",  make_sequence(750,  10, 100)),
            ("chr1000", make_sequence(1000, 10, 100)),
        ]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 5,
            "five distinct GC levels should produce 5 bins");
        assert_eq!(model.gc_ranges[0].start, 0,    "bin 0: start = 0 (0.0% GC)");
        assert_eq!(model.gc_ranges[0].end,   125,  "bin 0: end at midpoint 12.5%");
        assert_eq!(model.gc_ranges[1].start, 125,  "bin 1: start at midpoint 12.5%");
        assert_eq!(model.gc_ranges[1].end,   375,  "bin 1: end at midpoint 37.5%");
        assert_eq!(model.gc_ranges[2].start, 375,  "bin 2: start at midpoint 37.5%");
        assert_eq!(model.gc_ranges[2].end,   625,  "bin 2: end at midpoint 62.5%");
        assert_eq!(model.gc_ranges[3].start, 625,  "bin 3: start at midpoint 62.5%");
        assert_eq!(model.gc_ranges[3].end,   875,  "bin 3: end at midpoint 87.5%");
        assert_eq!(model.gc_ranges[4].start, 875,  "bin 4: start at midpoint 87.5%");
        assert_eq!(model.gc_ranges[4].end,   1000, "last bin end is always 1000");
        // Equal-population property: all bins get the same window count.
        let expected = model.gc_ranges[0].count;
        for (i, r) in model.gc_ranges.iter().enumerate() {
            assert_eq!(r.count, expected, "bin {} has unexpected window count", i);
        }
    }

    /// Same five-GC data but `num_gc_bins` capped at 3.
    ///
    /// The algorithm now tries ts = 3, 2, 1.  With ts=3, target_bin_size=17:
    ///   GC 0..249   → bin 0   (ceil(10/17)−1 = 0);  actual_max[0]=0
    ///   GC 250..749 → bin 1   (ceil(20/17)−1 = 1 and ceil(30/17)−1 = 1)
    ///                         actual_max[1]=500 (last observed GC in bin 1)
    ///   GC 750..1000→ bin 2   (ceil(40/17)−1 = 2 and ceil(50/17)−1 = 2)
    ///
    /// Midpoints: (0+250)/2=125, (500+750)/2=625
    ///
    /// Expected bins:
    ///   bin 0 : start=0,   end=125
    ///   bin 1 : start=125, end=625
    ///   bin 2 : start=625, end=1000
    #[test]
    fn test_max_bins_cap() {
        let tmp = std::env::temp_dir().join("garlic_test_max_bins.fa");
        write_fasta(&tmp, &[
            ("chr0",    make_sequence(0,    10, 100)),
            ("chr250",  make_sequence(250,  10, 100)),
            ("chr500",  make_sequence(500,  10, 100)),
            ("chr750",  make_sequence(750,  10, 100)),
            ("chr1000", make_sequence(1000, 10, 100)),
        ]);

        let cfg = basic_config(100, 3);   // cap at 3 bins
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 3,
            "num_gc_bins=3 must be respected even when more GC levels are present");
        // bin 1 spans GC=250 and GC=500; its actual_max is 500.
        // Midpoints: (actual_max[0]=0 + original_start[1]=250)/2=125
        //            (actual_max[1]=500 + original_start[2]=750)/2=625
        assert_eq!(model.gc_ranges[0].start, 0,    "bin 0: start = 0");
        assert_eq!(model.gc_ranges[0].end,   125,  "bin 0: end at midpoint 12.5%");
        assert_eq!(model.gc_ranges[1].start, 125,  "bin 1: start at midpoint 12.5%");
        assert_eq!(model.gc_ranges[1].end,   625,  "bin 1: end at midpoint 62.5%");
        assert_eq!(model.gc_ranges[2].start, 625,  "bin 2: start at midpoint 62.5%");
        assert_eq!(model.gc_ranges[2].end,   1000, "last bin end is always 1000");
    }

    // -----------------------------------------------------------------------
    // GC-to-GC transition matrix tests
    //
    // gc_trans[i][j] is the probability of moving from GC bin i to GC bin j
    // between consecutive 1kb windows.  Two structural properties are tested:
    //
    //   1. Diagonal matrix for separate chromosomes
    //      prev_gc_bin is reset to sentinel_bin at every SQ_SENTINEL, so
    //      transitions are never recorded across chromosome boundaries.
    //      Each chromosome is uniform → its windows only see self-transitions.
    //      Result: a diagonal matrix where gc_trans[i][i] = 1.0 and all
    //      off-diagonal entries are 0.0.
    //
    //   2. No reverse transition in a monotonic AT → GC sequence
    //      Within a single chromosome whose raw GC rises monotonically from 0
    //      to 1000, the smoothed (±5 window) average is also non-decreasing.
    //      Once the smoothed average crosses the bin boundary into bin 1 it
    //      never drops back below it, so gc_trans[1][0] = 0.0 exactly.
    // -----------------------------------------------------------------------

    /// Build an AC-alternating sequence (ACACAC…) of `num_windows` windows of
    /// `win_size` bases each.  Every window is exactly 50 % GC (score = 500).
    fn make_alternating_ac(num_windows: usize, win_size: usize) -> Vec<u8> {
        let total = num_windows * win_size;
        (0..total).map(|i| if i % 2 == 0 { b'A' } else { b'C' }).collect()
    }

    /// Helper: assert every row of `gc_trans` sums to 1.0 within floating-point
    /// tolerance.  Embedded in each transition test as a baseline sanity check.
    fn assert_rows_sum_to_one(gc_trans: &Vec<Vec<f64>>) {
        for (i, row) in gc_trans.iter().enumerate() {
            let sum: f64 = row.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "gc_trans row {} sums to {:.15} (expected 1.0)", i, sum
            );
        }
    }

    /// Two separate chromosomes (all-AT and all-GC) never share a transition:
    /// `prev_gc_bin` is reset to sentinel at each SQ_SENTINEL, so no
    /// cross-bin transition is ever recorded.  The result is a 2×2 diagonal
    /// matrix [ [1, 0], [0, 1] ].
    #[test]
    fn test_gc_trans_separate_chromosomes_diagonal() {
        let tmp = std::env::temp_dir().join("garlic_test_gc_trans_diag.fa");
        write_fasta(&tmp, &[
            ("chrAT", make_sequence(0,    10, 100)),
            ("chrGC", make_sequence(1000, 10, 100)),
        ]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 2);
        assert_rows_sum_to_one(&model.gc_trans);

        // Each bin transitions only to itself — strict diagonal.
        assert!((model.gc_trans[0][0] - 1.0).abs() < 1e-10,
            "bin 0 (all-AT chromosome) should be a pure self-loop");
        assert!((model.gc_trans[1][1] - 1.0).abs() < 1e-10,
            "bin 1 (all-GC chromosome) should be a pure self-loop");
        assert_eq!(model.gc_trans[0][1], 0.0,
            "no bin 0 → bin 1 transitions across chromosome boundary");
        assert_eq!(model.gc_trans[1][0], 0.0,
            "no bin 1 → bin 0 transitions across chromosome boundary");

        // Raw counts are kept for emit's Perl-style +1 pseudocount:
        // 10 windows per chromosome → 9 self-transitions each.
        assert_eq!(model.gc_trans_counts, vec![vec![9.0, 0.0], vec![0.0, 9.0]]);
    }

    /// A single chromosome whose first half is all-AT (raw GC = 0) and
    /// second half is all-GC (raw GC = 1000).
    ///
    /// The smoothed (±5 window) average is non-decreasing along this
    /// sequence — each step replaces one AT window with one GC window in
    /// the neighbourhood — so once the average crosses the bin boundary
    /// (midpoint = 500) into bin 1 it never drops back.
    ///
    /// Therefore gc_trans[1][0] must be exactly 0.0, and there must be at
    /// least one 0 → 1 transition at the boundary.
    #[test]
    fn test_gc_trans_unidirectional_crossing() {
        let tmp = std::env::temp_dir().join("garlic_test_gc_trans_unidir.fa");
        // 15 AT windows then 15 GC windows in a single chromosome.
        // 15 is large enough that the boundary region is well clear of the
        // chromosome edges and smoothing artefacts don't affect the interior.
        let seq: Vec<u8> = make_sequence(0, 15, 100)
            .into_iter()
            .chain(make_sequence(1000, 15, 100))
            .collect();
        write_fasta(&tmp, &[("chr1", seq)]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 2,
            "one AT block + one GC block should produce exactly 2 bins");
        assert_rows_sum_to_one(&model.gc_trans);

        // Once the smoothed GC crosses into bin 1 it never returns to bin 0.
        assert_eq!(model.gc_trans[1][0], 0.0,
            "no 1→0 transitions: smoothed GC is non-decreasing in AT→GC sequence");

        // There must be at least one forward crossing at the AT/GC boundary.
        assert!(model.gc_trans[0][1] > 0.0,
            "at least one 0→1 transition must occur at the AT/GC boundary");

        // Both bins must have some self-transitions (interior windows).
        assert!(model.gc_trans[0][0] > 0.0, "bin 0 must have self-transitions");
        assert!(model.gc_trans[1][1] > 0.0, "bin 1 must have self-transitions");
    }

    // -----------------------------------------------------------------------
    // K-mer frequency tests
    //
    // gc_binned_kmer_flank_freqs[bin][prefix] = [P(A|pre), P(C|pre), P(G|pre), P(T|pre)]
    // (accessed via the public gc_kmer_freq(bin, prefix) method).
    //
    // With kmer_size = 2 the prefix is a single base and the table captures
    // first-order Markov dependencies.  Three structural properties are tested:
    //
    //   1. Pure-sequence self-loop
    //      An all-A chromosome produces only "AA" kmers → P(A|"A") = 1.0.
    //
    //   2. Directional dependency
    //      An alternating ACACAC… chromosome produces only "AC" and "CA" kmers
    //      → P(C|"A") = 1.0  and  P(A|"C") = 1.0.
    //
    //   3. Bin isolation
    //      Kmers from separate chromosomes are counted only in the GC bin of
    //      the window containing their first base — they must not contaminate
    //      other bins.  When a prefix has no observed kmers in a bin, the code
    //      applies a GC-fraction fallback that differs from the data-driven 1.0.
    // -----------------------------------------------------------------------

    /// All-A chromosome (GC = 0, one bin): every kmer is "AA", so
    /// P(A | prefix="A") must be exactly 1.0 and all other conditionals 0.0.
    #[test]
    fn test_kmer_freq_pure_at_self_loop() {
        let tmp = std::env::temp_dir().join("garlic_test_kmer_pure_at.fa");
        write_fasta(&tmp, &[("chrAT", make_sequence(0, 10, 100))]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 1);

        // [P(A|A), P(C|A), P(G|A), P(T|A)]
        let f = model.gc_kmer_freq(0, "A")
            .expect("prefix 'A' must have an entry in an all-A sequence");
        assert!((f[0] - 1.0).abs() < 1e-10, "P(A|A) = 1.0 for all-A sequence, got {}", f[0]);
        assert_eq!(f[1], 0.0, "P(C|A) = 0.0 for all-A sequence");
        assert_eq!(f[2], 0.0, "P(G|A) = 0.0 for all-A sequence");
        assert_eq!(f[3], 0.0, "P(T|A) = 0.0 for all-A sequence");
        let sum: f64 = f.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "row sum must be 1.0");
    }

    /// Alternating ACACAC… chromosome (GC = 50 %, one bin).
    /// "A" is always followed by "C" and "C" is always followed by "A", so:
    ///   P(C | "A") = 1.0     P(A | "C") = 1.0
    /// all other conditionals = 0.0.
    #[test]
    fn test_kmer_freq_alternating_ac_dependency() {
        let tmp = std::env::temp_dir().join("garlic_test_kmer_altac.fa");
        write_fasta(&tmp, &[("chrAC", make_alternating_ac(10, 100))]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 1,
            "uniform 50% GC should produce a single bin");

        // A is always followed by C.
        let a = model.gc_kmer_freq(0, "A")
            .expect("prefix 'A' must have an entry in ACAC sequence");
        assert_eq!(a[0], 0.0, "P(A|A) = 0 in ACAC");
        assert!((a[1] - 1.0).abs() < 1e-10, "P(C|A) = 1.0 in ACAC, got {}", a[1]);
        assert_eq!(a[2], 0.0, "P(G|A) = 0 in ACAC");
        assert_eq!(a[3], 0.0, "P(T|A) = 0 in ACAC");

        // C is always followed by A.
        let c = model.gc_kmer_freq(0, "C")
            .expect("prefix 'C' must have an entry in ACAC sequence");
        assert!((c[0] - 1.0).abs() < 1e-10, "P(A|C) = 1.0 in ACAC, got {}", c[0]);
        assert_eq!(c[1], 0.0, "P(C|C) = 0 in ACAC");
        assert_eq!(c[2], 0.0, "P(G|C) = 0 in ACAC");
        assert_eq!(c[3], 0.0, "P(T|C) = 0 in ACAC");

        // Row sums for both observed prefixes.
        for prefix in &["A", "C"] {
            let sum: f64 = model.gc_kmer_freq(0, prefix).unwrap().iter().sum();
            assert!((sum - 1.0).abs() < 1e-10,
                "frequencies for prefix '{}' must sum to 1.0", prefix);
        }
    }

    /// Two separate chromosomes (all-AT → bin 0, all-GC → bin 1).
    ///
    /// Kmers are assigned to the bin of the window containing their first base,
    /// and each chromosome is uniform, so:
    ///   bin 0: only "AA" kmers counted → P(A | "A") = 1.0
    ///   bin 1: only "CC" kmers counted → P(C | "C") = 1.0
    ///
    /// For the opposite bins the prefix is never observed, so the code falls
    /// back to a GC-fraction prior (p_val / q_val) that is clearly ≠ 1.0,
    /// confirming no cross-bin contamination.
    #[test]
    fn test_kmer_freq_bin_isolation() {
        let tmp = std::env::temp_dir().join("garlic_test_kmer_bisolation.fa");
        write_fasta(&tmp, &[
            ("chrAT", make_sequence(0,    10, 100)),
            ("chrGC", make_sequence(1000, 10, 100)),
        ]);

        let cfg = basic_config(100, 25);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 2);

        // --- Bin 0 (all-AT): AA is the only kmer ---
        let bin0_a = model.gc_kmer_freq(0, "A")
            .expect("prefix 'A' must have an entry in bin 0");
        assert!((bin0_a[0] - 1.0).abs() < 1e-10,
            "bin 0: P(A|A) = 1.0 — AA is the only observed kmer");
        assert_eq!(bin0_a[1], 0.0,
            "bin 0: P(C|A) = 0.0 — no C bases in all-A sequence");

        // --- Bin 1 (all-GC): CC is the only kmer ---
        let bin1_c = model.gc_kmer_freq(1, "C")
            .expect("prefix 'C' must have an entry in bin 1");
        assert!((bin1_c[1] - 1.0).abs() < 1e-10,
            "bin 1: P(C|C) = 1.0 — CC is the only observed kmer");
        assert_eq!(bin1_c[0], 0.0,
            "bin 1: P(A|C) = 0.0 — no A bases in all-C sequence");

        // --- Cross-bin isolation ---
        // "A" prefix in bin 1 has no observed kmers → GC-fraction fallback,
        // which gives a value clearly different from the data-driven 1.0.
        let bin1_a = model.gc_kmer_freq(1, "A").unwrap();
        assert!((bin1_a[0] - 1.0).abs() > 0.01,
            "bin 1: P(A|A) must use fallback prior, not the data-driven 1.0");

        // "C" prefix in bin 0 has no observed kmers → same fallback logic.
        let bin0_c = model.gc_kmer_freq(0, "C").unwrap();
        assert!((bin0_c[1] - 1.0).abs() > 0.01,
            "bin 0: P(C|C) must use fallback prior, not the data-driven 1.0");

        // Row sums must hold for all entries (both data-driven and fallback).
        for (bin, prefix) in [(0usize, "A"), (1, "C")] {
            let sum: f64 = model.gc_kmer_freq(bin, prefix).unwrap().iter().sum();
            assert!((sum - 1.0).abs() < 1e-10,
                "bin {}, prefix '{}': frequencies must sum to 1.0", bin, prefix);
        }
    }

    /// Masked regions spanning FASTA line breaks must not advance the position
    /// counter on newline bytes.  Previously each line break inside a mask was
    /// counted as a base, so the counter drifted ahead and the mask (and GC-bin
    /// lookup) ended early.  Here 600 masked A's span ten 60-bp lines; with the
    /// drift the last ten A's escaped the mask and produced AA k-mers.
    #[test]
    fn test_mask_position_ignores_line_breaks() {
        let tmp = std::env::temp_dir().join("garlic_test_mask_linebreaks.fa");
        let mut seq = vec![b'A'; 600];
        seq.extend(b"ACGT".repeat(150));
        write_fasta(&tmp, &[("chr1", seq)]);

        let mut ivs = crate::io::IntervalMap::new();
        ivs.insert("chr1".to_string(), vec![(0, 600)]);
        let mut mask = crate::masking::MaskSet::new();
        mask.add(ivs);

        let mut cfg = basic_config(100, 1);
        cfg.mask = Some(&mask);
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        assert_eq!(model.gc_ranges.len(), 1);
        let a = model.gc_kmer_freq(0, "A").expect("prefix 'A' must be observed");
        assert_eq!(a[0], 0.0, "P(A|A) must be 0: every A-run is masked");
        assert!((a[1] - 1.0).abs() < 1e-10, "P(C|A) = 1.0 in the ACGT repeat");
    }

    /// Models from other format versions (or with no version field, i.e.
    /// rusted-garlic 0.1.x) are refused with a "rebuild" message; current
    /// models round-trip.
    #[test]
    fn test_model_format_version_check() {
        let old = SeqModel::from_json(r#"{"a_count": 0}"#);
        let msg = old.err().expect("unversioned model must be refused").to_string();
        assert!(msg.contains("rebuild"), "unexpected error: {}", msg);

        let tmp = std::env::temp_dir().join("garlic_test_model_version.fa");
        write_fasta(&tmp, &[("chr1", make_sequence(500, 20, 100))]);
        let model = SeqModel::model(&tmp, Some(&basic_config(100, 1)));
        let _ = fs::remove_file(&tmp);
        let json = serde_json::to_string(&model).unwrap();
        let back = SeqModel::from_json(&json).expect("current model must load");
        assert_eq!(back.model_version, MODEL_FORMAT_VERSION);
    }

    // -----------------------------------------------------------------------
    // Default vs --perl_compat GC binning
    // -----------------------------------------------------------------------

    fn hist(levels: &[(usize, u32)]) -> [u32; 1001] {
        let mut h = [0u32; 1001];
        for &(v, n) in levels {
            h[v] = n;
        }
        h
    }

    /// A GC value holding more than one bin's worth of windows: the default
    /// lowers the bin count so populations stay even; Perl keeps more bins but
    /// one absorbs the spike.
    #[test]
    fn test_skewed_histogram_default_vs_perl() {
        // 40 windows at 40.0% GC, then 10 each at 30, 50, 60 and 70%.
        let h = hist(&[(300, 10), (400, 40), (500, 10), (600, 10), (700, 10)]);
        let (starts, ranges) = equal_population_bins(&h, 8);
        // 80 windows: 8..4 bins leave a bin empty (the spike spans several);
        // 3 bins of target 27 work: 10 / 40 / 30 windows.
        assert_eq!(starts, vec![300, 400, 500]);
        let r: Vec<(u32, u32)> = ranges.iter().map(|r| (r.start, r.end)).collect();
        assert_eq!(r, vec![(300, 350), (350, 450), (450, 1000)]);
        // Perl keeps 5 bins: [0,300) is empty, [400,500) holds the spike (40),
        // and [600,1000] takes the top two levels (20).
        assert_eq!(gc_bin_cuts(&h, 8), vec![0, 300, 400, 500, 600, 1000]);
    }

    /// Perl cut rule: a cut wherever int(cumulative / (total / bins)) rises;
    /// the first cut is 0 and the highest observed value closes at 1000.
    #[test]
    fn test_gc_bin_cuts_perl_rule() {
        let h = hist(&[(100, 10), (300, 10), (500, 10), (700, 10)]);
        assert_eq!(gc_bin_cuts(&h, 4), vec![0, 100, 300, 500, 1000]);
        let h = hist(&[(200, 10), (400, 30), (600, 10)]);
        assert_eq!(gc_bin_cuts(&h, 5), vec![0, 200, 400, 1000]);
    }

    /// Degenerate inputs where Perl yields no usable bins give one 0-1000 bin.
    #[test]
    fn test_gc_bin_cuts_degenerate() {
        assert_eq!(gc_bin_cuts(&hist(&[(0, 20)]), 25), vec![0, 1000]);
        assert_eq!(gc_bin_cuts(&hist(&[(500, 20)]), 25), vec![0, 1000]);
        assert_eq!(gc_bin_cuts(&hist(&[(1000, 20)]), 25), vec![0, 1000]);
    }

    /// Rounded to 0.1%, then the last bin starting at or below the value.
    #[test]
    fn test_gc_bin_of() {
        let starts = [0, 342, 365];
        assert_eq!(gc_bin_of(0.0, &starts), 0);
        assert_eq!(gc_bin_of(0.3414, &starts), 0);
        assert_eq!(gc_bin_of(0.3416, &starts), 1, "rounds to 0.342 = start of bin 1");
        assert_eq!(gc_bin_of(0.3649, &starts), 2, "rounds to 0.365");
        assert_eq!(gc_bin_of(1.0, &starts), 2);
        assert_eq!(gc_bin_of(0.1, &[300, 500]), 0, "below the first start -> bin 0");
    }

    /// --perl_compat on five GC levels: Perl's cuts, with the top level (100%)
    /// sharing the 75% bin.
    #[test]
    fn test_perl_compat_five_levels() {
        let tmp = std::env::temp_dir().join("garlic_test_perl_compat_bins.fa");
        write_fasta(&tmp, &[
            ("chr0",    make_sequence(0,    10, 100)),
            ("chr250",  make_sequence(250,  10, 100)),
            ("chr500",  make_sequence(500,  10, 100)),
            ("chr750",  make_sequence(750,  10, 100)),
            ("chr1000", make_sequence(1000, 10, 100)),
        ]);
        let mut cfg = basic_config(100, 25);
        cfg.perl_compat = true;
        let model = SeqModel::model(&tmp, Some(&cfg));
        let _ = fs::remove_file(&tmp);

        let ranges: Vec<(u32, u32, u32)> =
            model.gc_ranges.iter().map(|r| (r.start, r.end, r.count)).collect();
        assert_eq!(ranges, vec![(0, 250, 10), (250, 500, 10), (500, 750, 10), (750, 1000, 20)]);
    }

    /// A trailing partial window never counts as a window for bins or
    /// transitions.  By default it is binned by its own GC (so repeats in the
    /// tail use it); --perl_compat repeats the last full window's GC instead.
    #[test]
    fn test_partial_window() {
        let tmp = std::env::temp_dir().join("garlic_test_partial_window.fa");
        // chrA: 20 x 50% windows + a 60-bp all-C tail; chrB: 20 x 100% windows.
        let mut a = make_sequence(500, 20, 100);
        a.extend(vec![b'C'; 60]);
        write_fasta(&tmp, &[("chrA", a), ("chrB", make_sequence(1000, 20, 100))]);

        for perl_compat in [false, true] {
            let mut cfg = basic_config(100, 2);
            cfg.perl_compat = perl_compat;
            let model = SeqModel::model(&tmp, Some(&cfg));
            let counts: Vec<u32> = model.gc_ranges.iter().map(|r| r.count).collect();
            let total: u32 = counts.iter().sum();
            assert_eq!(total, 40, "only full windows are counted (perl_compat={perl_compat})");
            let trans: f64 = model.gc_trans_counts.iter().flatten().sum();
            assert_eq!(trans, 38.0, "19 transitions per chromosome (perl_compat={perl_compat})");
        }
        let _ = fs::remove_file(&tmp);
    }
}
