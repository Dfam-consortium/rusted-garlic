/// MaskSet — the union of all annotation-derived mask intervals.
///
/// Intervals are stored per sequence-id as sorted, merged, 0-based half-open
/// ranges [start, end).  During a streaming FASTA pass a MaskCursor is
/// obtained for each chromosome and advanced position-by-position; it answers
/// `is_masked(pos)` in O(1) amortised time.
use std::collections::HashMap;

use crate::io::IntervalMap;

pub struct MaskSet {
    intervals: HashMap<String, Vec<(u32, u32)>>,
}

impl MaskSet {
    pub fn new() -> Self {
        MaskSet {
            intervals: HashMap::new(),
        }
    }

    /// Merge another IntervalMap into this MaskSet (e.g. from a second annotation file).
    pub fn add(&mut self, map: IntervalMap) {
        for (chrom, mut ivs) in map {
            let entry = self.intervals.entry(chrom).or_default();
            entry.append(&mut ivs);
        }
        // Re-sort and merge after every addition so cursor invariants hold.
        crate::io::sort_and_merge(&mut self.intervals);
    }

    /// Return the total number of masked intervals across all chromosomes.
    pub fn interval_count(&self) -> usize {
        self.intervals.values().map(|v| v.len()).sum()
    }

    /// Obtain a forward-only cursor for the given chromosome.
    /// If the chromosome has no intervals the cursor is empty.
    pub fn cursor_for<'a>(&'a self, chrom: &str) -> MaskCursor<'a> {
        let intervals: &'a [(u32, u32)] = self
            .intervals
            .get(chrom)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        MaskCursor { intervals, idx: 0 }
    }
}

/// A forward-only cursor over a chromosome's sorted mask intervals.
///
/// Positions **must** be queried in non-decreasing order (streaming).
/// The cursor advances its internal index past intervals that ended before
/// the current position, so each call is O(1) amortised.
pub struct MaskCursor<'a> {
    intervals: &'a [(u32, u32)],
    idx: usize,
}

impl<'a> MaskCursor<'a> {
    /// Returns true if 0-based position `pos` falls inside any mask interval.
    pub fn is_masked(&mut self, pos: u32) -> bool {
        // Advance past intervals that ended before this position.
        while self.idx < self.intervals.len() && self.intervals[self.idx].1 <= pos {
            self.idx += 1;
        }
        // Check whether pos falls inside the current interval.
        self.idx < self.intervals.len() && self.intervals[self.idx].0 <= pos
    }

    /// Reset the cursor to the beginning (call when starting a new chromosome).
    pub fn reset(&mut self) {
        self.idx = 0;
    }
}
