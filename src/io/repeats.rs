/// Parser for the rusted-garlic repeat annotation TSV format.
///
/// The TSV has 12 tab-separated columns with a 0-based half-open coordinate
/// convention (BED-style).  Lines beginning with `#` are comments and are
/// ignored.
///
/// Columns (0-based index):
///   0  chrom       sequence name
///   1  start       0-based half-open start
///   2  end         0-based half-open end
///   3  name        repeat name (e.g. AluYa5)
///   4  family      class/family (e.g. SINE/Alu)
///   5  strand      + or -
///   6  div         % substitution divergence
///   7  ins         % insertion rate
///   8  del         % deletion rate
///   9  cons_start  consensus start coordinate
///  10  cons_end    consensus end coordinate
///  11  rep_id      unique instance ID; use "." for no fragment grouping
///
/// `rep_id` ties multiple fragments of the same repeat instance together and
/// drives nested-insertion detection.  Use `.` when the source data does not
/// track instance identity; each such row is treated as an independent
/// unfragmented element.
///
/// Simple/tandem repeats should be excluded by the upstream converter script.
use std::io::BufRead;
use std::path::Path;

use super::IntervalMap;

/// Build a masked-interval map from a repeat TSV for use during k-mer
/// counting (Pass 2).  Only the first three columns (chrom, start, end)
/// are read; coordinates are already 0-based half-open.
pub fn parse(path: impl AsRef<Path>) -> IntervalMap {
    let reader = super::open_reader(path.as_ref());
    let mut map: IntervalMap = IntervalMap::new();

    for line in reader.lines() {
        let line = line.expect("Read error in repeat TSV");
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // Only need first three columns — use splitn to avoid scanning the rest.
        let fields: Vec<&str> = trimmed.splitn(4, '\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let chrom = fields[0].to_string();
        let start: u32 = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: u32 = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if end <= start {
            continue;
        }
        // Coordinates are already 0-based half-open — no conversion needed.
        map.entry(chrom).or_default().push((start, end));
    }

    super::sort_and_merge(&mut map);
    map
}
