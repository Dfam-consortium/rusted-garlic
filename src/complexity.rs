/// Sequence complexity and composition measures.
///
/// Mirrors the SeqComplex.pm Perl module (Caballero & Hubley, ISB 2010-2015).
///
/// All measures are computed from a single pass over the window bytes via
/// [`WindowCounts::new`].  Non-ACGT characters (N, gaps, …) are counted
/// separately and break k-mer runs but do not cause errors.
///
/// # Measures
///
/// **Composition** — `gc`, `at`, `gcs`, `ats`, `cpg`, `ket`, `pur`
///
/// **Entropy** — `cwf` (Wootton-Federhen), `ce` (Shannon)
///
/// **K-mer diversity** — `cm`, `ct`, `cl` (each at word sizes `min_k..=max_k`)
///
/// **Compression** — `cz` (gzip ratio, uses `flate2` in-process)

use flate2::{write::GzEncoder, Compression};
use std::io::Write;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Precomputed 4^k for k = 0..=6.
const POW4: [usize; 7] = [1, 4, 16, 64, 256, 1024, 4096];

// ---------------------------------------------------------------------------
// Base encoding
// ---------------------------------------------------------------------------

/// Map an uppercase ACGT byte to 0/1/2/3; return `None` for anything else.
#[inline]
fn base_val(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// WindowCounts
// ---------------------------------------------------------------------------

/// Per-window accumulators built in a single pass over the sequence.
///
/// Holds everything needed to compute all complexity measures without
/// re-scanning the bytes.
pub struct WindowCounts {
    /// Counts of A, C, G, T (indices 0-3).
    pub bases: [u32; 4],
    /// Total ACGT bases (`bases.iter().sum()`).
    pub total: u32,
    /// Total non-ACGT characters (N, gaps, …).
    pub n_count: u32,
    /// Count of adjacent "CG" dinucleotides (for CpG ratio).
    pub cg_dinuc: u32,
    /// `kmer_counts[ki]` is a flat `Vec<u32>` of length `4^(min_k + ki)`.
    /// Entry `kmer_counts[ki][code]` is the number of times the k-mer
    /// encoded as `code` appears in the window.
    pub kmer_counts: Vec<Vec<u32>>,
    pub min_k: usize,
    pub max_k: usize,
}

impl WindowCounts {
    /// Build counts for `seq` (must be uppercase ACGT/N) for k-mer sizes
    /// `min_k..=max_k` (1 ≤ min_k ≤ max_k ≤ 6).
    pub fn new(seq: &[u8], min_k: usize, max_k: usize) -> Self {
        assert!(min_k >= 1 && max_k <= 6 && min_k <= max_k,
                "k-mer sizes must satisfy 1 ≤ min_k ≤ max_k ≤ 6");

        let n_k = max_k - min_k + 1;
        let mut kmer_counts: Vec<Vec<u32>> = (min_k..=max_k)
            .map(|k| vec![0u32; POW4[k]])
            .collect();

        let mut bases = [0u32; 4];
        let mut n_count = 0u32;
        let mut cg_dinuc = 0u32;

        // Running state for the k-mer sliding window.
        // kmer_state[ki] = encoded value of the last (min_k+ki) valid bases.
        // kmer_valid[ki] = consecutive valid bases seen so far.
        let mut kmer_state = vec![0u32; n_k];
        let mut kmer_valid = vec![0usize; n_k];
        let mut prev_c = false; // was the previous ACGT base a 'C'?

        for &byte in seq {
            match base_val(byte) {
                Some(v) => {
                    bases[v as usize] += 1;

                    // CpG dinucleotide detection.
                    if v == 2 && prev_c {
                        cg_dinuc += 1;
                    }
                    prev_c = v == 1;

                    // Update running k-mer states for each requested word size.
                    for ki in 0..n_k {
                        let k = min_k + ki;
                        let mask = (POW4[k] - 1) as u32;
                        kmer_state[ki] = ((kmer_state[ki] << 2) | v as u32) & mask;
                        kmer_valid[ki] += 1;
                        if kmer_valid[ki] >= k {
                            kmer_counts[ki][kmer_state[ki] as usize] += 1;
                        }
                    }
                }
                None => {
                    n_count += 1;
                    prev_c = false;
                    // Any non-ACGT base breaks the k-mer run.
                    for ki in 0..n_k {
                        kmer_valid[ki] = 0;
                        kmer_state[ki] = 0;
                    }
                }
            }
        }

        let total: u32 = bases.iter().sum();
        WindowCounts { bases, total, n_count, cg_dinuc, kmer_counts, min_k, max_k }
    }

    /// Total window length (ACGT + non-ACGT).
    #[inline]
    pub fn len(&self) -> usize {
        (self.total + self.n_count) as usize
    }

    /// Fraction of the window that is non-ACGT.
    #[inline]
    pub fn n_frac(&self) -> f64 {
        let l = self.len();
        if l == 0 { 0.0 } else { self.n_count as f64 / l as f64 }
    }

    /// K-mer count table for word size `k`.
    pub fn kmers(&self, k: usize) -> &[u32] {
        &self.kmer_counts[k - self.min_k]
    }

    /// Number of distinct k-mers observed for word size `k`.
    pub fn kmer_diversity(&self, k: usize) -> usize {
        self.kmers(k).iter().filter(|&&c| c > 0).count()
    }

    // -----------------------------------------------------------------------
    // Composition measures
    // -----------------------------------------------------------------------

    /// GC content: (C + G) / total.
    pub fn gc(&self) -> f64 {
        if self.total < 2 { return 0.0; }
        (self.bases[1] + self.bases[2]) as f64 / self.total as f64
    }

    /// AT content: (A + T) / total.
    pub fn at(&self) -> f64 {
        if self.total < 2 { return 0.0; }
        (self.bases[0] + self.bases[3]) as f64 / self.total as f64
    }

    /// GC skew: |G - C| / (G + C).
    pub fn gcs(&self) -> f64 {
        let gc = self.bases[1] + self.bases[2];
        if gc < 2 { return 0.0; }
        (self.bases[2] as i64 - self.bases[1] as i64).unsigned_abs() as f64 / gc as f64
    }

    /// AT skew: |A - T| / (A + T).
    pub fn ats(&self) -> f64 {
        let at = self.bases[0] + self.bases[3];
        if at < 2 { return 0.0; }
        (self.bases[0] as i64 - self.bases[3] as i64).unsigned_abs() as f64 / at as f64
    }

    /// CpG ratio: observed_CG / (C_count × G_count).
    pub fn cpg(&self) -> f64 {
        let expected = self.bases[1] as u64 * self.bases[2] as u64;
        if expected < 2 { return 0.0; }
        self.cg_dinuc as f64 / expected as f64
    }

    /// Keto skew: |(G + T) - (A + C)| / total.
    pub fn ket(&self) -> f64 {
        if self.total < 2 { return 0.0; }
        let gt = (self.bases[2] + self.bases[3]) as i64;
        let ac = (self.bases[0] + self.bases[1]) as i64;
        (gt - ac).unsigned_abs() as f64 / self.total as f64
    }

    /// Purine skew: |(G + A) - (T + C)| / total.
    pub fn pur(&self) -> f64 {
        if self.total < 2 { return 0.0; }
        let ga = (self.bases[2] + self.bases[0]) as i64;
        let tc = (self.bases[3] + self.bases[1]) as i64;
        (ga - tc).unsigned_abs() as f64 / self.total as f64
    }

    // -----------------------------------------------------------------------
    // Entropy measures
    // -----------------------------------------------------------------------

    /// Wootton-Federhen complexity.
    ///
    /// `cwf = (win × log₄(win) − Σ log₄(nᵢ)) / total`
    ///
    /// `win` is the raw window size (including non-ACGT characters), matching
    /// the Perl implementation.
    pub fn cwf(&self, win: usize) -> f64 {
        if self.total < 2 { return 0.0; }
        let ln4 = (4.0_f64).ln();
        let log4_win = (win as f64).ln() / ln4;
        let up = log4_win * win as f64;
        let dw: f64 = self.bases.iter()
            .filter(|&&c| c > 0)
            .map(|&c| (c as f64).ln() / ln4)
            .sum();
        (up - dw) / self.total as f64
    }

    /// Shannon entropy (bits): −Σ pᵢ × log₂(pᵢ).
    pub fn ce(&self) -> f64 {
        if self.total < 2 { return 0.0; }
        let tot = self.total as f64;
        self.bases.iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / tot;
                -p * p.log2()
            })
            .sum()
    }

    // -----------------------------------------------------------------------
    // K-mer-based measures
    // -----------------------------------------------------------------------

    /// Markov complexity for word size `k`.
    ///
    /// `cm = −Σ p_w × log_{4^k}(p_w)` where `p_w = count_w / (win − k + 1)`.
    pub fn cm(&self, k: usize, win: usize) -> f64 {
        let dw = (win - k + 1) as f64;
        let ln_m = ((POW4[k]) as f64).ln(); // ln(4^k) = k × ln(4)
        self.kmers(k).iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let r = c as f64 / dw;
                -r * r.ln() / ln_m
            })
            .sum()
    }

    /// Linguistic complexity for word sizes 1..=`k`.
    ///
    /// `cl = Σ vl / Σ vm`  where `vm = min(4^l, win − l + 1)` and `vl` is the
    /// observed distinct l-mer count.
    pub fn cl(&self, k: usize, win: usize) -> f64 {
        let mut sum_vl = 0usize;
        let mut sum_vm = 0usize;
        for l in 1..=k {
            sum_vm += POW4[l].min(win - l + 1);
            sum_vl += self.kmer_diversity(l);
        }
        if sum_vm == 0 { 0.0 } else { sum_vl as f64 / sum_vm as f64 }
    }

    /// Trifonov complexity for word sizes 1..=`k`.
    ///
    /// `ct = Π (vl / vm)` for l = 1..=k.
    pub fn ct(&self, k: usize, win: usize) -> f64 {
        let mut ct = 1.0f64;
        for l in 1..=k {
            let vm = POW4[l].min(win - l + 1);
            if vm == 0 { return 0.0; }
            ct *= self.kmer_diversity(l) as f64 / vm as f64;
        }
        ct
    }

    // -----------------------------------------------------------------------
    // Compression measure
    // -----------------------------------------------------------------------

    /// Compression complexity: raw_length / gzip_compressed_length.
    ///
    /// Higher values = more repetitive (compresses well).
    /// Uses `flate2` in-process — no temporary files.
    pub fn cz(&self, seq: &[u8]) -> f64 {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        if enc.write_all(seq).is_err() { return 0.0; }
        match enc.finish() {
            Ok(compressed) if !compressed.is_empty() => {
                seq.len() as f64 / compressed.len() as f64
            }
            _ => 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// TSV output helpers
// ---------------------------------------------------------------------------

/// Build the TSV header line for the given k-mer range.
pub fn tsv_header(min_k: usize, max_k: usize) -> String {
    let cols = vec![
        "file", "seq", "start", "end", "win_len", "acgt_count",
        "gc", "at", "gcs", "ats", "cpg", "ket", "pur",
        "cwf", "ce", "cz",
    ];
    let mut kmer_cols: Vec<String> = Vec::new();
    for prefix in &["cm", "ct", "cl"] {
        for k in min_k..=max_k {
            kmer_cols.push(format!("{}{}", prefix, k));
        }
    }
    let col_refs: Vec<&str> = cols.iter().map(|s| *s)
        .chain(kmer_cols.iter().map(|s| s.as_str()))
        .collect();
    col_refs.join("\t")
}

/// Compute all measures for one window and return a tab-separated data row.
///
/// `seq` must be the raw window bytes (uppercase ACGT/N).
/// `win` is the nominal window size (= `seq.len()` unless truncated at end of sequence).
/// `file` and `seq_name` are carried through for the TSV columns.
pub fn tsv_row(
    seq: &[u8],
    file: &str,
    seq_name: &str,
    start: usize,
    min_k: usize,
    max_k: usize,
) -> String {
    let win = seq.len();
    // ct and cl at word size k use the counts of every size 1..=k, so count
    // from 1; min_k only selects the output columns.
    let wc = WindowCounts::new(seq, 1, max_k);
    let end = start + win;

    let mut vals: Vec<f64> = vec![
        wc.gc(),
        wc.at(),
        wc.gcs(),
        wc.ats(),
        wc.cpg(),
        wc.ket(),
        wc.pur(),
        wc.cwf(win),
        wc.ce(),
        wc.cz(seq),
    ];

    for k in min_k..=max_k { vals.push(wc.cm(k, win)); }
    for k in min_k..=max_k { vals.push(wc.ct(k, win)); }
    for k in min_k..=max_k { vals.push(wc.cl(k, win)); }

    let float_cols: String = vals.iter()
        .map(|v| format!("{:.6}", v))
        .collect::<Vec<_>>()
        .join("\t");

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        file, seq_name, start, end, win, wc.total, float_cols
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Restricting the output columns (min_k > 1) must not change ct/cl,
    /// which need the counts of every word size from 1.
    #[test]
    fn test_min_k_only_selects_columns() {
        let seq = b"ACGTTGCAAGGCTTACGATCGATCGGATCCATGCA".repeat(10);
        let all = tsv_row(&seq, "f", "s", 0, 1, 4);
        let some = tsv_row(&seq, "f", "s", 0, 3, 4);
        let cols = |row: &str| row.split('\t').map(String::from).collect::<Vec<_>>();
        let hdr_all = cols(&tsv_header(1, 4));
        let hdr_some = cols(&tsv_header(3, 4));
        let (a, b) = (cols(&all), cols(&some));
        for (i, h) in hdr_some.iter().enumerate() {
            let j = hdr_all.iter().position(|x| x == h).unwrap();
            assert_eq!(b[i], a[j], "column {h}");
        }
    }
}
