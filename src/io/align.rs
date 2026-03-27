/// Parser for RepeatMasker `.align` and `.out` annotation files.
///
/// Both formats share the same fixed-column header line layout.  A data line
/// (after stripping leading whitespace) begins with:
///
///   score  pct_div  pct_del  pct_ins  query_id  query_start  query_end  ...
///
/// Coordinates are 1-based, fully-closed [start, end].
///
/// Lines whose class/family field contains Simple_repeat, Low_complexity,
/// Unknown, or Satellite are skipped — those are handled separately via the
/// TRF BED file.  The returned intervals are 0-based half-open [start-1, end).
use std::io::BufRead;
use std::path::Path;

use super::IntervalMap;

pub fn parse(path: impl AsRef<Path>) -> IntervalMap {
    let reader = super::open_reader(path.as_ref());
    let mut map: IntervalMap = IntervalMap::new();

    for line in reader.lines() {
        let line = line.expect("Read error in repeat file");

        // Data lines start with optional whitespace then a score (integer)
        // followed by a decimal divergence value.
        let trimmed = line.trim_start();
        if !looks_like_data_line(trimmed) {
            continue;
        }

        // Skip classes that belong to the simple-repeat / TRF model.
        if trimmed.contains("Simple_repeat")
            || trimmed.contains("Low_complexity")
            || trimmed.contains("Unknown")
            || trimmed.contains("Satellite")
        {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_ascii_whitespace().collect();
        if fields.len() < 7 {
            continue;
        }

        let seq_id = fields[4].to_string();
        let start: u32 = match fields[5].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: u32 = match fields[6].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if end < start {
            continue;
        }

        // Convert 1-based closed [start, end] → 0-based half-open [start-1, end)
        let start0 = start - 1;
        map.entry(seq_id).or_default().push((start0, end));
    }

    super::sort_and_merge(&mut map);
    map
}

/// True if the trimmed line begins with an integer score followed by a
/// decimal divergence value — the signature of a RM data line.
fn looks_like_data_line(s: &str) -> bool {
    // Expect: digits, whitespace, digits.digits
    let mut chars = s.chars().peekable();
    // Must start with a digit
    if !chars.peek().map_or(false, |c| c.is_ascii_digit()) {
        return false;
    }
    while chars.peek().map_or(false, |c| c.is_ascii_digit()) {
        chars.next();
    }
    // whitespace
    if !chars.peek().map_or(false, |c| c.is_ascii_whitespace()) {
        return false;
    }
    while chars.peek().map_or(false, |c| c.is_ascii_whitespace()) {
        chars.next();
    }
    // digits . digits
    if !chars.peek().map_or(false, |c| c.is_ascii_digit()) {
        return false;
    }
    while chars.peek().map_or(false, |c| c.is_ascii_digit()) {
        chars.next();
    }
    chars.peek() == Some(&'.')
}
