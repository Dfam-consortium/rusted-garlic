use rand::distributions::WeightedIndex;
use rand::prelude::*;
use std::cmp;
use std::collections::HashMap;
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
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
struct GCRange {
    start: u32,
    end: u32,
    count: u32,
}

#[derive(Serialize, Deserialize)]
pub struct SeqModel {
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
    gc_binned_kmer_flank_freqs: Vec<HashMap<String, [f64; 4]>>,
    pub repeat_model: Option<RepeatModel>,
}

impl SeqModel {
    pub fn new() -> Self {
        SeqModel {
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
            gc_binned_kmer_flank_freqs: Vec::new(),
            repeat_model: None,
        }
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

        let module: SeqModel = serde_json::from_str(&contents.as_str())?;
        Ok(module)
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
    pub fn emit(&self, seq_size: Option<usize>, emit_cfg: Option<&EmitConfig>) -> (String, Vec<InsertRecord>) {
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

        let mut gc_to_gc_weights: Vec<WeightedIndex<f64>> = Vec::new();
        for i in 0..self.gc_ranges.len() {
            let dist = WeightedIndex::new(&self.gc_trans[i]).unwrap();
            gc_to_gc_weights.push(dist);
        }

        let mut rng = thread_rng();

        // Sample initial GC bin weighted by bin window-count (Perl's newGC()).
        let gc_bin_counts: Vec<u32> = self.gc_ranges.iter().map(|r| r.count).collect();
        let gc_bin_dist = WeightedIndex::new(&gc_bin_counts).unwrap();
        let mut current_gc_bin = gc_bin_dist.sample(&mut rng);

        let bp_to_gen: usize = seq_size.unwrap_or(self.source_sequence_size);
        let mut out_bytes: Vec<u8> = Vec::with_capacity(bp_to_gen + 1);

        // Per-window GC bin assignments: gc_bins[pos / WIN_SIZE] = gc_bin used there.
        let mut gc_bins: Vec<u32> = Vec::new();
        gc_bins.push(current_gc_bin as u32); // window 0

        // Seed with a random (k-1)-mer prefix from the starting bin (Perl's $fseed).
        let prefixes: Vec<&String> = gc_kmer_prefix_weights[current_gc_bin].keys().collect();
        let seed_prefix = prefixes[rng.gen_range(0..prefixes.len())].clone();
        out_bytes.extend_from_slice(seed_prefix.as_bytes());
        let mut km = seed_prefix;

        while out_bytes.len() < bp_to_gen {
            let nextbase = gc_kmer_prefix_weights[current_gc_bin]
                .get(&km)
                .expect("missing kmer prefix in model")
                .sample(&mut rng);
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
                current_gc_bin = gc_to_gc_weights[current_gc_bin].sample(&mut rng);
                gc_bins.push(current_gc_bin as u32);
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
                        &mut rng,
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
                    &mut rng,
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
        let (raw, raw_size) = crate::io::open_raw(path.as_ref());
        let est_windows: usize = raw_size / window_size;
        let mut bingc: Vec<u32> = Vec::with_capacity(est_windows);
        let mut buf = BufReader::new(raw);

        let mut bytes = [0; 512];
        let mut in_header = false;
        let mut win_at_count = 0;
        let mut win_cg_count = 0;
        let mut win_bp_count = 0;
        const NA_SENTINEL: u32 = 1001;
        const SQ_SENTINEL: u32 = 1002;
        let mut cumul_gc_win_counts = [0u32; (NA_SENTINEL + 1) as usize];

        // Per-chromosome tracking for repeat profiling (Pass 3).
        // chrom_names[k] = the sequence id whose SQ_SENTINEL is bingc[k_th_sentinel].
        let mut chrom_names: Vec<String> = Vec::new();
        let mut p1_chr_name = String::new();
        let mut p1_in_header_name = false;

        // Calculate:
        //
        //   bingc[]: GC content of "window_size" windows as a scaled fraction (to limit
        //            the precision: int(cg_count/total_count * 1000).
        //
        //   cumul_gc_win_counts[]: Count of windows for each GC bin to later be converted
        //                          into a cumulative count.
        loop {
            match buf.read(&mut bytes) {
                // EOF Reached
                Ok(0) => {
                    // Process partial window
                    if win_bp_count > 0 {
                        // NOTE: Garlic uses GC of the previous window...as the
                        //       trailing case. I am not sure if this was to avoid
                        //       calling GC if the window was small ( e.g one bp ).
                        let mut gc = NA_SENTINEL;
                        let acgt_count = win_cg_count + win_at_count;
                        // Perhaps only calculate GC iff there is enough data and use previous
                        // GC otherwise
                        if acgt_count > 0 {
                            gc = ((win_cg_count as f64 / acgt_count as f64) * 1000.0) as u32;
                        }

                        // To simulate GARLIC approach
                        //if let Some(val) = bingc.last() {
                        //  gc = val;
                        //}

                        // Emit a window
                        bingc.push(gc);
                        cumul_gc_win_counts[gc as usize] += 1;
                    }
                    // Emit sequence end sentinel
                    bingc.push(SQ_SENTINEL);
                    chrom_names.push(p1_chr_name.clone());
                    break;
                }
                Ok(n) => {
                    for i in 0..n {
                        //num_bytes += 1;

                        if win_bp_count == window_size {
                            let mut gc = NA_SENTINEL;
                            let acgt_count = win_cg_count + win_at_count;
                            if acgt_count > 0 {
                                gc = ((win_cg_count as f64 / acgt_count as f64) * 1000.0) as u32;
                            }
                            // Emit a window
                            bingc.push(gc);
                            cumul_gc_win_counts[gc as usize] += 1;
                            win_at_count = 0;
                            win_cg_count = 0;
                            win_bp_count = 0;
                        }

                        if bytes[i] == '>' as u8 {
                            // Process partial window
                            if win_bp_count > 0 {
                                // See EOF section about trailing cases
                                let mut gc = NA_SENTINEL;
                                let acgt_count = win_cg_count + win_at_count;
                                if acgt_count > 0 {
                                    gc =
                                        ((win_cg_count as f64 / acgt_count as f64) * 1000.0) as u32;
                                }
                                bingc.push(gc);
                                cumul_gc_win_counts[gc as usize] += 1;
                            }
                            // Always emit sequence end sentinel (fixes chromosome
                            // boundary tracking when a chr ends on a window edge).
                            if !p1_chr_name.is_empty() || !bingc.is_empty() {
                                bingc.push(SQ_SENTINEL);
                                chrom_names.push(p1_chr_name.clone());
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

        // Read stats
        //let dur = t1.elapsed().as_secs_f64();
        //let num_bytes = (num_bytes as f64) / 1_048_576f64;
        //println!(
        //    "Read speed, Pass 1: {:.1} MiB / s for {}",
        //    num_bytes / dur,
        //    num_bytes
        //);

        // Currently cumul_gc_win_counts contains absolute counts of windows
        // for each gc value 0-1000,NA_SENTINEL.  Convert this to a cumulative
        // count from gc=0-1000 inplace.
        let mut min_gc = 1000;
        let mut max_gc = 0;
        if cumul_gc_win_counts[0 as usize] > 0 {
            min_gc = 0;
        }
        for v in 1..NA_SENTINEL {
            if cumul_gc_win_counts[v as usize] > 0 {
                if v < min_gc {
                    min_gc = v;
                }
                max_gc = v;
            }
            cumul_gc_win_counts[v as usize] += cumul_gc_win_counts[(v - 1) as usize];
        }

        println!(
            "min_gc = {}, max_gc = {} dynamic_range={}",
            min_gc,
            max_gc,
            max_gc - min_gc + 1
        );
        println!("max_bins = {}", (max_gc - min_gc + 1 / 25) + 1);
        println!(
            "new max bins = {}",
            std::cmp::min((max_gc - min_gc + 1 / 25) + 1, 25)
        );

        let _new_target = std::cmp::min((max_gc - min_gc + 1 / 25) + 1, 25);

        // Identify the GC ranges that divide the genome into MAX_GC_BUCKETS number
        // of roughly equally sized groups of windows.
        //
        // Calculates:
        //   gc_ranges[] : A vector of structures containing start, end and bp count for each
        //                 generated range.
        //   gc_lookup[] : A convenience vector to lookup bins by any scaled GC value.
        //
        let mut gc_lookup = vec![0u32; SQ_SENTINEL as usize];
        let mut local_gc_ranges: Vec<GCRange> =
            vec![GCRange { start: 0, end: 0, count: 0 }; max_gc_bins];
        // Tracks the highest GC value (0–1000 scale) that has at least one window
        // in each bin for the current ts trial.  Used to compute midpoint boundaries.
        let mut bin_actual_max = vec![0u32; max_gc_bins];
        let mut bad = 0;
        let mut num_buckets: u32 = 0;
        for ts in (1..=max_gc_bins as u32).rev() {
            num_buckets = ts;
            for r in local_gc_ranges.iter_mut() {
                *r = GCRange { start: 0, end: 0, count: 0 };
            }
            for m in bin_actual_max.iter_mut() { *m = 0; }
            println!("Trying bin count = {}", ts);
            let target_bin_size = (cumul_gc_win_counts[1000] as f64 / ts as f64).ceil() as u32;
            println!(
                "cumul_gc_win_counts[1000] = {} , target_bin_size = {}",
                cumul_gc_win_counts[1000], target_bin_size
            );
            let mut prev_bin = 10000000;
            for v in 0..NA_SENTINEL {
                let mut bin = 0;
                if cumul_gc_win_counts[v as usize] > 0 {
                    bin = ((cumul_gc_win_counts[v as usize] as f64 / target_bin_size as f64).ceil()
                        - 1.0)
                        .abs() as u32;
                }
                if bin != prev_bin {
                    println!(
                        "  - bin switch @{} counts={}",
                        v, cumul_gc_win_counts[v as usize]
                    );
                }
                gc_lookup[v as usize] = bin;
                // Track the highest GC score that has actual windows in each bin.
                let window_count_at_v = cumul_gc_win_counts[v as usize]
                    - if v > 0 { cumul_gc_win_counts[(v - 1) as usize] } else { 0 };
                if window_count_at_v > 0 {
                    bin_actual_max[bin as usize] = v;
                }
                if local_gc_ranges[bin as usize].count > 0 {
                    local_gc_ranges[bin as usize].count +=
                        cumul_gc_win_counts[v as usize] - cumul_gc_win_counts[(v - 1) as usize];
                } else {
                    local_gc_ranges[bin as usize].start = v;
                    if v == 0 {
                        local_gc_ranges[bin as usize].count = cumul_gc_win_counts[v as usize];
                    } else {
                        local_gc_ranges[bin as usize].count =
                            cumul_gc_win_counts[v as usize] - cumul_gc_win_counts[(v - 1) as usize];
                    }
                }
                prev_bin = bin;
            }
            bad = 0;
            for v in 0..ts as usize {
                println!(
                    "bin[{}]:  {:.3} to {:.3} (inclusive) [len={}]",
                    v,
                    (local_gc_ranges[v].start as f64 / 1000 as f64),
                    (local_gc_ranges[v].end as f64 / 1000 as f64),
                    local_gc_ranges[v].count
                );
                if local_gc_ranges[v].count == 0 {
                    bad = 1;
                }
                if v > 0 && local_gc_ranges[v].start <= local_gc_ranges[v - 1].start {
                    bad = 1;
                }
                if v > 0 && local_gc_ranges[v].start > local_gc_ranges[v - 1].start {
                    local_gc_ranges[v - 1].end = local_gc_ranges[v].start - 1;
                }
            }
            local_gc_ranges[(ts - 1) as usize].end = 1000 as u32;
            if bad == 0 {
                break;
            }
            // TODO: Calculate the max GC range for all bins...and warn if it seems excessive
        }
        if bad == 1 {
            panic!("Could not identify a GC bucket count that satisfies the constraints of zero-empty buckets, and .1% seperation between buckets");
        }
        // Sentinel value for "no bin" — equals actual_bins (one past last valid index).
        let actual_bins = num_buckets as usize;
        gc_lookup[NA_SENTINEL as usize] = actual_bins as u32;

        // Refine bin boundaries to midpoints between observed GC clusters.
        //
        // The raw boundary is  end[i] = start[i+1] - 1, which can look arbitrary when
        // the two neighbouring bins have a large gap between their observed GC values
        // (e.g. bins at 0 % and 100 % get "0–99.9%" and "100–100%" instead of the
        // cleaner "0–50%" and "50–100%").  Using the midpoint makes the ranges easier
        // to interpret without changing which windows fall in which bin (the gc_lookup
        // table governs actual assignment and is unaffected).
        let original_starts: Vec<u32> = local_gc_ranges[..actual_bins].iter().map(|r| r.start).collect();
        for i in 0..actual_bins - 1 {
            let midpoint = (bin_actual_max[i] + original_starts[i + 1]) / 2;
            local_gc_ranges[i].end = midpoint;
            local_gc_ranges[i + 1].start = midpoint;
        }
        local_gc_ranges[actual_bins - 1].end = 1000;

        // Copy valid GC range boundaries into the model (counts recomputed in smoothing pass).
        seq_model.gc_ranges = local_gc_ranges[..actual_bins]
            .iter()
            .map(|r| GCRange { start: r.start, end: r.end, count: 0 })
            .collect();
        seq_model.gc_trans = vec![vec![0.0f64; actual_bins]; actual_bins];

        if print_garlic {
            println!("Initial bin counts");
            for v in 0..actual_bins {
                println!(
                    "bin[{}]:  {:.3} to {:.3} (inclusive) [len={}]",
                    v,
                    (local_gc_ranges[v].start as f64 / 1000 as f64),
                    (local_gc_ranges[v].end as f64 / 1000 as f64),
                    local_gc_ranges[v].count
                );
            }
        }

        // Calculate the smoothed GC content of each window and update the GC range counts
        // This alters gc_lookup using the smoothed GC window values.

        // Counts are recomputed during the smoothing pass below.
        const FLANK_WINDOWS: usize = 5;
        let mut next_dest_idx = 0;
        let mut gc_queue: VecDeque<u32> = VecDeque::with_capacity(14);
        let mut gc_sum = 0;
        let mut gc_cnt = 0;
        let sentinel_bin = actual_bins as u32; // "no bin" sentinel value
        let mut prev_gc_bin = sentinel_bin;

        // chrom_gc_bins[chrom] = smoothed GC bin index for each 1kb window.
        // Used by the repeat profiler to look up GC bin at any genomic position.
        let mut chrom_gc_bins: HashMap<String, Vec<u32>> = HashMap::new();
        let mut sentinel_idx: usize = 0; // index into chrom_names

        for i in 0..bingc.len() {
            if bingc[i] == SQ_SENTINEL {
                let cur_chrom = chrom_names.get(sentinel_idx).cloned().unwrap_or_default();
                // finish up trailing windows for the current chromosome
                for _n in 0..cmp::min(gc_queue.len(), 6) {
                    let mut gc_avg = NA_SENTINEL;
                    if gc_cnt > 0 {
                        gc_avg = (gc_sum as f64 / gc_cnt as f64) as u32;
                    }
                    let gc_bin = gc_lookup[gc_avg as usize];
                    if prev_gc_bin != sentinel_bin && gc_bin != sentinel_bin {
                        seq_model.gc_trans[prev_gc_bin as usize][gc_bin as usize] += 1.0 as f64;
                    }
                    prev_gc_bin = gc_bin;
                    if gc_bin < sentinel_bin {
                        seq_model.gc_ranges[gc_bin as usize].count += 1;
                        bingc[next_dest_idx] = gc_bin;
                    }
                    // Always push so window index stays aligned (same fix as the
                    // left-edge case above).
                    chrom_gc_bins.entry(cur_chrom.clone()).or_default().push(gc_bin);
                    next_dest_idx += 1;

                    if gc_queue.len() >= FLANK_WINDOWS + 1 {
                        let val = gc_queue.pop_front().unwrap();
                        if val < NA_SENTINEL {
                            gc_sum -= val;
                            gc_cnt -= 1;
                        }
                    }
                }
                sentinel_idx += 1;
                gc_sum = 0;
                gc_cnt = 0;
                prev_gc_bin = sentinel_bin;
                gc_queue.clear();
                continue;
            }

            if bingc[i] < NA_SENTINEL {
                gc_sum += bingc[i];
                gc_cnt += 1;
            }
            gc_queue.push_back(bingc[i]);

            // Complete windows case:
            // FLANK_WINDOWS .. |--------| .. FLANK_WINDOWS
            if gc_queue.len() > ((FLANK_WINDOWS * 2) + 1) {
                let val = gc_queue.pop_front().unwrap();
                if val < NA_SENTINEL {
                    gc_sum -= val;
                    gc_cnt -= 1;
                }
            }

            // left-edge case: flank on right reached
            // |--------| .. FLANK_WINDOWS
            if gc_queue.len() > FLANK_WINDOWS {
                let mut gc_avg = NA_SENTINEL;
                if gc_cnt > 0 {
                    gc_avg = (gc_sum as f64 / gc_cnt as f64) as u32;
                }
                let gc_bin = gc_lookup[gc_avg as usize];
                if prev_gc_bin != sentinel_bin && gc_bin != sentinel_bin {
                    seq_model.gc_trans[prev_gc_bin as usize][gc_bin as usize] += 1.0 as f64;
                }
                prev_gc_bin = gc_bin;
                let cur_chrom = chrom_names.get(sentinel_idx).cloned().unwrap_or_default();
                if gc_bin < sentinel_bin {
                    seq_model.gc_ranges[gc_bin as usize].count += 1;
                }
                // Always push (including sentinel_bin) so that chrom_gc_bins[w]
                // corresponds to window index w within the chromosome —
                // previously sentinel windows were skipped, causing an offset.
                chrom_gc_bins.entry(cur_chrom).or_default().push(gc_bin);
                bingc[next_dest_idx] = gc_bin;
                next_dest_idx += 1;
            }
        }
        bingc.resize(next_dest_idx, NA_SENTINEL);

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
        // (due to SQ_SENTINEL entries) and panicked on genomes with many sequences.
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
                        let gc_bin = gc_bin_for(&chr_name, chr_bp_count - 1);
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
                            let gc_bin = gc_bin_for(&chr_name, chr_bp_count - 1);
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
                                                 actual_bins, &mut rm, cfg.gene_filter);
                    log::info!("  {} TE records profiled", rm.te_count());
                    te_ivs
                } else {
                    crate::io::IntervalMap::new()
                };
                if let Some(tpath) = cfg.tandem_annot_file {
                    let te_opt = if te_ivs.is_empty() { None } else { Some(&te_ivs) };
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
            let hm: HashMap<String, [f64; 4]> = HashMap::new();
            seq_model.gc_binned_kmer_flank_freqs.push(hm);
        }

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
    //      the window containing their last base — they must not contaminate
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
    /// Kmers are assigned to the bin of the window containing their last base,
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
}
