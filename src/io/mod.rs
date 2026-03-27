pub mod bed;
pub mod fasta;
pub mod repeats;

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// A map of sequence-id → sorted list of 0-based half-open intervals [start, end)
/// that should be masked (treated as non-ACGT) during k-mer profiling.
pub type IntervalMap = std::collections::HashMap<String, Vec<(u32, u32)>>;

// ---------------------------------------------------------------------------
// Transparent gz / plain file openers
// ---------------------------------------------------------------------------

fn is_gzip(path: &Path) -> bool {
    path.extension().and_then(|e| e.to_str()) == Some("gz")
}

/// Open a file for buffered line-by-line reading.
/// Transparently decompresses `.gz` files using gzip streaming.
pub fn open_reader(path: impl AsRef<Path>) -> Box<dyn BufRead> {
    let path = path.as_ref();
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("Cannot open {:?}: {}", path, e));
    if is_gzip(path) {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

/// Open a file for raw chunk reading (e.g. the FASTA pass in `seq_model`).
/// Returns the reader and an estimated *uncompressed* byte count
/// (used only for `Vec::with_capacity` hints — not exact for `.gz`).
pub fn open_raw(path: impl AsRef<Path>) -> (Box<dyn Read>, usize) {
    let path = path.as_ref();
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("Cannot open {:?}: {}", path, e));
    let raw_size = file.metadata().map(|m| m.len()).unwrap_or(0) as usize;
    if is_gzip(path) {
        // FASTA .gz files typically compress 3–5×; use 4× as a capacity hint.
        (Box::new(flate2::read::GzDecoder::new(file)), raw_size * 4)
    } else {
        (Box::new(file), raw_size)
    }
}

/// Sort and merge overlapping intervals within each chromosome of an IntervalMap.
pub fn sort_and_merge(map: &mut IntervalMap) {
    for intervals in map.values_mut() {
        intervals.sort_unstable_by_key(|&(s, _)| s);
        let mut merged: Vec<(u32, u32)> = Vec::new();
        for &(s, e) in intervals.iter() {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        *intervals = merged;
    }
}
