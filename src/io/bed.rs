/// Parsers for BED-format annotation files used in masking.
///
/// Two formats are supported:
///
/// 1. TRF BED (UCSC format, 0-based half-open):
///       chrom  start  end  ...
///    Columns 0, 1, 2.  All remaining columns are ignored.
///
/// 2. UCSC gene table (ensGene / refGene, 0-based half-open):
///       bin  name  chrom  strand  txStart  txEnd  ...
///    Columns 2, 4, 5 give chrom, txStart, txEnd.
///
/// Both return intervals as 0-based half-open [start, end) ready for the
/// MaskSet, so no coordinate conversion is needed.
use std::io::BufRead;
use std::path::Path;

use super::IntervalMap;

/// Parse a TRF BED file (UCSC format, 0-based half-open).
/// Fields: chrom(0)  start(1)  end(2)  [...]
pub fn parse_trf(path: impl AsRef<Path>) -> IntervalMap {
    parse_bed_columns(path, 0, 1, 2)
}

/// Parse a UCSC gene/transcript table (ensGene / refGene format).
/// Fields: bin(0)  name(1)  chrom(2)  strand(3)  txStart(4)  txEnd(5)  [...]
pub fn parse_gene_table(path: impl AsRef<Path>) -> IntervalMap {
    parse_bed_columns(path, 2, 4, 5)
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
