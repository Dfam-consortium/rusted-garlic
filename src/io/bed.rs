/// Parsers for BED-like annotation files used in masking and profiling filters.
use std::io::BufRead;
use std::path::Path;

use super::IntervalMap;

/// Parse a BED-like annotation file with chrom/start/end in the first three columns.
/// Used for tandem repeat annotations or any interval-only source.
/// Additional columns are ignored.  Files may be tab- or whitespace-delimited.
/// Comment lines (`#`) and lines with fewer than 3 fields are skipped.
/// Coordinates must already be 0-based half-open.
pub fn parse_annot3(path: impl AsRef<Path>) -> IntervalMap {
    parse_bed_columns(path, 0, 1, 2)
}

/// Parse a gene annotation file and return two interval maps:
///   (gene_body_intervals, exon_intervals)
///
/// Format (tab-separated, 0-based half-open):
///   chrom  txStart  txEnd  [name  strand  exonStarts  exonEnds]
///
/// Columns 0–2 are required (gene body).  Columns 5–6 are the exon coordinate
/// lists (comma-separated integers, trailing comma allowed).  Use "." for either
/// exon column when exon coordinates are not available — the exon map will be
/// empty in that case.
///
/// Example line:
///   chr1  11868  14409  DDX11L1  +  11868,12612,13220  12227,12721,14409
pub fn parse_gene_annot(path: impl AsRef<Path>) -> (IntervalMap, IntervalMap) {
    let reader = super::open_reader(path.as_ref());
    let mut bodies: IntervalMap = IntervalMap::new();
    let mut exons: IntervalMap = IntervalMap::new();

    for line in reader.lines() {
        let line = line.expect("Read error in gene annotation file");
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let f: Vec<&str> = trimmed.split('\t').collect();
        if f.len() < 3 {
            continue;
        }

        let chrom = f[0];
        let tx_start: u32 = match f[1].parse() { Ok(v) => v, Err(_) => continue };
        let tx_end: u32   = match f[2].parse() { Ok(v) => v, Err(_) => continue };
        if tx_end <= tx_start {
            continue;
        }

        bodies.entry(chrom.to_string()).or_default().push((tx_start, tx_end));

        // Parse exon columns (cols 5 and 6) when present and not "."
        if f.len() >= 7 && f[5] != "." && f[6] != "." {
            let starts: Vec<u32> = f[5].split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse().ok())
                .collect();
            let ends: Vec<u32> = f[6].split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse().ok())
                .collect();
            let n = starts.len().min(ends.len());
            let entry = exons.entry(chrom.to_string()).or_default();
            for i in 0..n {
                if ends[i] > starts[i] {
                    entry.push((starts[i], ends[i]));
                }
            }
        }
    }

    super::sort_and_merge(&mut bodies);
    super::sort_and_merge(&mut exons);
    (bodies, exons)
}

/// General BED-column reader.
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

        map.entry(seq_id).or_default().push((start, end));
    }

    super::sort_and_merge(&mut map);
    map
}
