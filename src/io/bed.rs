/// Parsers for BED-like annotation files used in masking.
///
/// All three annotation types (TE, tandem repeat, gene) share the same
/// 3-column BED-like format (0-based half-open coordinates):
///
///   chrom  start  end  [additional columns ignored]
///
/// Files may be tab- or whitespace-delimited.  Comment lines (`#`) and
/// lines with fewer than 3 fields are skipped.  No coordinate conversion
/// is performed — coordinates are expected to already be 0-based half-open.
use std::io::BufRead;
use std::path::Path;

use super::IntervalMap;

/// Parse a BED-like annotation file with chrom/start/end in the first three columns.
/// Used for tandem repeat annotations, gene annotations, or any other interval source.
pub fn parse_annot3(path: impl AsRef<Path>) -> IntervalMap {
    parse_bed_columns(path, 0, 1, 2)
}

/// General BED-column reader.  `id_col`, `start_col`, `end_col` specify which
/// whitespace-delimited fields hold the sequence id, start, and end.
/// Comment lines (starting with `#`) and lines with too few fields are skipped.
fn parse_bed_columns(
    path: impl AsRef<Path>,
    id_col: usize,
    start_col: usize,
    end_col: usize,
) -> IntervalMap {
    let reader = super::open_reader(path.as_ref());
    let mut map: IntervalMap = IntervalMap::new();
    let min_cols = id_col.max(start_col).max(end_col) + 1;

    for line in reader.lines() {
        let line = line.expect("Read error in BED file");
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_ascii_whitespace().collect();
        if fields.len() < min_cols {
            continue;
        }

        let seq_id = fields[id_col].to_string();
        let start: u32 = match fields[start_col].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: u32 = match fields[end_col].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if end <= start {
            continue;
        }

        // Already 0-based half-open — no conversion needed.
        map.entry(seq_id).or_default().push((start, end));
    }

    super::sort_and_merge(&mut map);
    map
}
