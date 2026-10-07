/// Repeat profile model — built from RepeatMasker .align and TRF BED files.
///
/// Mirrors the Perl `profileRMAlign` / `profileTRF` / `parseRepeats` logic in
/// `createModel.pl`.  The profile is stored per GC bin (indexed 0..num_bins-1)
/// and serialised alongside the sequence model in the JSON model file.
use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::io::IntervalMap;

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

/// One interspersed (TE) repeat record — the fields needed to simulate an
/// evolved copy.
///
/// Serialized as a compact JSON array `[name, family, dir, div, ins, del, length, nfrg]`
/// rather than a named-field object, to avoid repeating field names across the
/// ~2 M records typically present in a genome model.
#[derive(Debug, Clone)]
pub struct TeRepeat {
    pub name: String,   // repeat name   e.g. "AluYa5"
    pub family: String, // class/family  e.g. "SINE/Alu"
    pub dir: char,      // strand        '+' or '-'
    pub div: f64,       // max % divergence across fragments
    pub ins: f64,       // max % insertion rate
    pub del: f64,       // max % deletion rate
    pub length: u32,    // consensus length spanned (end - start of repeat coords)
    pub nfrg: u32,      // number of fragments (1 = unfragmented)
}

impl Serialize for TeRepeat {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let mut seq = s.serialize_seq(Some(8))?;
        seq.serialize_element(&self.name)?;
        seq.serialize_element(&self.family)?;
        seq.serialize_element(&self.dir.to_string())?;
        seq.serialize_element(&self.div)?;
        seq.serialize_element(&self.ins)?;
        seq.serialize_element(&self.del)?;
        seq.serialize_element(&self.length)?;
        seq.serialize_element(&self.nfrg)?;
        seq.end()
    }
}

impl<'de> Deserialize<'de> for TeRepeat {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> serde::de::Visitor<'de> for V {
            type Value = TeRepeat;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("[name, family, dir, div, ins, del, length, nfrg]")
            }
            fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut s: A) -> Result<TeRepeat, A::Error> {
                use serde::de::Error;
                let name:   String = s.next_element()?.ok_or_else(|| A::Error::invalid_length(0, &self))?;
                let family: String = s.next_element()?.ok_or_else(|| A::Error::invalid_length(1, &self))?;
                let dirs:   String = s.next_element()?.ok_or_else(|| A::Error::invalid_length(2, &self))?;
                let dir = dirs.chars().next().ok_or_else(|| A::Error::custom("empty dir"))?;
                let div:    f64    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(3, &self))?;
                let ins:    f64    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(4, &self))?;
                let del:    f64    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(5, &self))?;
                let length: u32    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(6, &self))?;
                let nfrg:   u32    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(7, &self))?;
                Ok(TeRepeat { name, family, dir, div, ins, del, length, nfrg })
            }
        }
        d.deserialize_seq(V)
    }
}

/// One simple/tandem repeat record.
///
/// Serialized as a compact JSON array `[consensus, dir, copy_num, div, indel]`.
#[derive(Debug, Clone)]
pub struct SimpleRepeat {
    pub consensus: String, // repeat unit e.g. "AGAT"
    pub dir: char,         // canonical strand after revcomp check
    pub copy_num: f64,     // TRF copy number (BED column 6; Perl $period / $exp)
    pub div: f64,          // divergence (100 - percent_match from TRF)
    pub indel: f64,        // indel rate from TRF
}

impl Serialize for SimpleRepeat {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let mut seq = s.serialize_seq(Some(5))?;
        seq.serialize_element(&self.consensus)?;
        seq.serialize_element(&self.dir.to_string())?;
        seq.serialize_element(&self.copy_num)?;
        seq.serialize_element(&self.div)?;
        seq.serialize_element(&self.indel)?;
        seq.end()
    }
}

impl<'de> Deserialize<'de> for SimpleRepeat {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> serde::de::Visitor<'de> for V {
            type Value = SimpleRepeat;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("[consensus, dir, copy_num, div, indel]")
            }
            fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut s: A) -> Result<SimpleRepeat, A::Error> {
                use serde::de::Error;
                let consensus: String = s.next_element()?.ok_or_else(|| A::Error::invalid_length(0, &self))?;
                let dirs:      String = s.next_element()?.ok_or_else(|| A::Error::invalid_length(1, &self))?;
                let dir = dirs.chars().next().ok_or_else(|| A::Error::custom("empty dir"))?;
                let copy_num:  f64    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(2, &self))?;
                let div:       f64    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(3, &self))?;
                let indel:     f64    = s.next_element()?.ok_or_else(|| A::Error::invalid_length(4, &self))?;
                Ok(SimpleRepeat { consensus, dir, copy_num, div, indel })
            }
        }
        d.deserialize_seq(V)
    }
}

/// The full repeat profile for one genome, indexed by GC bin.
#[derive(Debug, Serialize, Deserialize)]
pub struct RepeatModel {
    /// `te_profile[gc_bin]` = all TE records assigned to that bin
    pub te_profile: Vec<Vec<TeRepeat>>,
    /// `simple_profile[gc_bin]` = all simple repeat records assigned to that bin
    pub simple_profile: Vec<Vec<SimpleRepeat>>,
    /// `insert_table[gc_bin][outer_type#family][inner_type#family]` = frequency
    /// (BTreeMap: deterministic serialization and sampling order)
    pub insert_table: Vec<BTreeMap<String, BTreeMap<String, f64>>>,
}

impl RepeatModel {
    pub fn new(num_bins: usize) -> Self {
        RepeatModel {
            te_profile: vec![Vec::new(); num_bins],
            simple_profile: vec![Vec::new(); num_bins],
            insert_table: vec![BTreeMap::new(); num_bins],
        }
    }

    pub fn te_count(&self) -> usize {
        self.te_profile.iter().map(|v| v.len()).sum()
    }

    pub fn simple_count(&self) -> usize {
        self.simple_profile.iter().map(|v| v.len()).sum()
    }
}

// ---------------------------------------------------------------------------
// GC bin lookup helper
// ---------------------------------------------------------------------------

/// Returns the smoothed GC bin index for a given (seq_id, 0-based position).
/// Returns None when the chromosome is unknown or the window is out of range.
pub fn get_gc_bin(
    chrom_gc_bins: &HashMap<String, Vec<u32>>,
    seq_id: &str,
    pos: u32,
    win_size: u32,
) -> Option<u32> {
    let bins = chrom_gc_bins.get(seq_id)?;
    let win = (pos / win_size) as usize;
    bins.get(win).copied()
}

// ---------------------------------------------------------------------------
// Interval filter helpers
// ---------------------------------------------------------------------------

/// Returns true if [start, end) overlaps any interval in the sorted list.
/// Requires the list to be sorted and non-overlapping (output of sort_and_merge).
fn overlaps_any(intervals: &[(u32, u32)], start: u32, end: u32) -> bool {
    if intervals.is_empty() || end <= start {
        return false;
    }
    // First interval whose end is > start (the first candidate that could overlap).
    let pos = intervals.partition_point(|&(_, e)| e <= start);
    pos < intervals.len() && intervals[pos].0 < end
}

/// Returns true if [start, end) is fully contained within any interval in the sorted list.
/// Requires the list to be sorted and non-overlapping (output of sort_and_merge).
fn contained_in_any(intervals: &[(u32, u32)], start: u32, end: u32) -> bool {
    if intervals.is_empty() || end <= start {
        return false;
    }
    // Last interval whose start is <= start.
    let pos = intervals.partition_point(|&(s, _)| s <= start);
    if pos == 0 {
        return false;
    }
    intervals[pos - 1].1 >= end
}

// ---------------------------------------------------------------------------
// Repeat TSV profiler
// ---------------------------------------------------------------------------

/// Parse a rusted-garlic repeat annotation TSV and populate `model` with
/// TE repeat records and the nested-insert table.
///
/// Expected columns (tab-separated, 0-based half-open coordinates):
///   chrom  start  end  name  family  strand  div  ins  del  cons_start  cons_end  rep_id
///
/// When `rep_id` is `.` each row is treated as an independent unfragmented
/// element (no fragment assembly, no nesting detection).
///
/// `gene_filter`: when provided, TE instances that overlap any interval in this
/// map are excluded from the TE profile (intron filter).
///
/// Returns an `IntervalMap` of all TE genomic intervals (sorted, merged per
/// chromosome) for use by `profile_trf` to exclude tandem repeats that fall
/// entirely within a TE.
pub fn profile_repeats(
    path: impl AsRef<Path>,
    chrom_gc_bins: &HashMap<String, Vec<u32>>,
    win_size: u32,
    num_bins: usize,
    model: &mut RepeatModel,
    gene_filter: Option<&IntervalMap>,
) -> IntervalMap {
    let reader = crate::io::open_reader(path.as_ref());

    // Accumulates per-instance data while streaming.
    // Key: rep_id; value: accumulated record built from one or more fragments.
    struct RepData {
        label: String,    // semicolon-joined fragment params
        ini: u32,         // earliest genomic start (0-based)
        end: u32,         // latest genomic end (0-based half-open)
        seq_id: String,
        type_fam: String, // "name#family"
    }
    let mut repdata: HashMap<String, RepData> = HashMap::new();
    // Collect all TE genomic intervals for use by profile_trf.
    let mut te_intervals: IntervalMap = IntervalMap::new();

    // For nested-insert detection: track the type#family at each data-row
    // index, and the last row index where each rep_id appeared.
    let mut row_rep_name: Vec<Option<String>> = Vec::new();
    let mut id_last_row: HashMap<String, usize> = HashMap::new();

    // Temporary insert counts: [gc_bin][outer][inner] = count
    let mut insert_counts: Vec<HashMap<String, HashMap<String, u32>>> =
        vec![HashMap::new(); num_bins];

    let mut row_idx: usize = 0;
    let mut skipped_lines: usize = 0;

    for raw in reader.lines() {
        let line = raw.expect("Read error in repeat TSV");
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        row_rep_name.push(None);
        let cur_row = row_idx;
        row_idx += 1;

        let f: Vec<&str> = trimmed.split('\t').collect();
        if f.len() < 12 {
            skipped_lines += 1;
            continue;
        }

        let seq_id    = f[0];
        let ini: u32  = f[1].parse().unwrap_or(0);  // 0-based half-open start
        let end: u32  = f[2].parse().unwrap_or(0);  // 0-based half-open end
        let name      = f[3];
        let family    = f[4];
        let dir       = f[5].chars().next().unwrap_or('+');
        let div: f64  = f[6].parse().unwrap_or(0.0);
        let ins: f64  = f[7].parse().unwrap_or(0.0);
        let del: f64  = f[8].parse().unwrap_or(0.0);
        let rini: u32 = f[9].parse().unwrap_or(0);
        let rend: u32 = f[10].parse().unwrap_or(0);
        let rep_id    = f[11];

        let type_fam = format!("{}#{}", name, family);
        let fragment_label = format!(
            "{}:{}:{}:{:.2}:{:.2}:{:.2}:{}:{}",
            name, family, dir, div, ins, del, rini, rend
        );

        row_rep_name[cur_row] = Some(type_fam.clone());

        // Collect TE interval for use by profile_trf regardless of rep_id.
        te_intervals.entry(seq_id.to_string()).or_default().push((ini, end));

        // When rep_id is "." treat as independent — no fragment assembly or nesting.
        if rep_id == "." {
            // Apply gene filter: skip if this TE overlaps a gene region.
            if let Some(gf) = gene_filter {
                if let Some(gene_ivs) = gf.get(seq_id) {
                    if overlaps_any(gene_ivs, ini, end) {
                        continue;
                    }
                }
            }
            let pos = ini + (end.saturating_sub(ini)) / 2;
            let gc = match get_gc_bin(chrom_gc_bins, seq_id, pos, win_size) {
                Some(g) => g as usize,
                None => continue,
            };
            if let Some(te) = parse_te_label(&fragment_label) {
                model.te_profile[gc].push(te);
            }
            continue;
        }

        // Nested insert detection: elements appearing between two occurrences
        // of the same rep_id are nested inside that repeat.
        if let Some(&prev_row) = id_last_row.get(rep_id) {
            if let Some(rd) = repdata.get(rep_id) {
                let pos = rd.ini + (end.saturating_sub(rd.ini)) / 2;
                if let Some(gc_bin) = get_gc_bin(chrom_gc_bins, rd.seq_id.as_str(), pos, win_size) {
                    let outer = rd.type_fam.clone();
                    for ln in (prev_row + 1)..cur_row {
                        if let Some(Some(inner)) = row_rep_name.get(ln) {
                            if !inner.is_empty() {
                                *insert_counts[gc_bin as usize]
                                    .entry(outer.clone())
                                    .or_default()
                                    .entry(inner.clone())
                                    .or_default() += 1;
                            }
                        }
                    }
                }
            }
        }
        id_last_row.insert(rep_id.to_string(), cur_row);

        // Accumulate or create the repeat record for this rep_id.
        repdata
            .entry(rep_id.to_string())
            .and_modify(|rd| {
                rd.label.push(';');
                rd.label.push_str(&format!(
                    "{:.2}:{:.2}:{:.2}:{}:{}",
                    div, ins, del, rini, rend
                ));
                if end > rd.end {
                    rd.end = end;
                }
            })
            .or_insert_with(|| RepData {
                label: fragment_label,
                ini,
                end,
                seq_id: seq_id.to_string(),
                type_fam,
            });
    }

    // Sort and merge TE intervals so the returned map is ready for binary-search queries.
    crate::io::sort_and_merge(&mut te_intervals);

    // Assign each accumulated repeat instance to a GC bin, applying gene filter.
    // Genomic order (rep_id breaks ties), not HashMap order, so te_profile is
    // the same on every run.
    let mut repdata: Vec<(String, RepData)> = repdata.into_iter().collect();
    repdata.sort_by(|(ka, a), (kb, b)| {
        (&a.seq_id, a.ini, a.end, ka).cmp(&(&b.seq_id, b.ini, b.end, kb))
    });
    for (_, rd) in repdata.iter() {
        // Apply gene filter: skip if this TE overlaps a gene region.
        if let Some(gf) = gene_filter {
            if let Some(gene_ivs) = gf.get(&rd.seq_id) {
                if overlaps_any(gene_ivs, rd.ini, rd.end) {
                    continue;
                }
            }
        }
        let pos = rd.ini + (rd.end.saturating_sub(rd.ini)) / 2;
        let gc = match get_gc_bin(chrom_gc_bins, &rd.seq_id, pos, win_size) {
            Some(g) => g as usize,
            None => continue,
        };
        if let Some(te) = parse_te_label(&rd.label) {
            model.te_profile[gc].push(te);
        }
    }

    // Convert insert counts to frequencies.
    for (gc_bin, outer_map) in insert_counts.into_iter().enumerate() {
        for (outer, inner_map) in outer_map {
            let total: u32 = inner_map.values().sum();
            if total == 0 {
                continue;
            }
            let freq_map: BTreeMap<String, f64> = inner_map
                .into_iter()
                .map(|(k, v)| (k, v as f64 / total as f64))
                .collect();
            model.insert_table[gc_bin].insert(outer, freq_map);
        }
    }

    if skipped_lines > 0 && model.te_count() == 0 {
        log::error!(
            "0 TE records loaded from {:?} — all {} non-comment lines were skipped \
             (fewer than 12 tab-separated fields). \
             Did you pass a raw RepeatMasker .align file? \
             Convert it first with: python3 utils/rm_align_to_tsv.py <file> > <file>.tsv",
            path.as_ref(),
            skipped_lines
        );
        std::process::exit(1);
    }

    te_intervals
}

/// Parse a semicolon-joined multi-fragment TE label into a `TeRepeat`.
/// Format: `name:fam:dir:div:ins:del:rini:rend[;div:ins:del:rini:rend...]`
fn parse_te_label(label: &str) -> Option<TeRepeat> {
    let frags: Vec<&str> = label.split(';').collect();
    let first: Vec<&str> = frags[0].split(':').collect();
    if first.len() < 8 {
        return None;
    }
    let name = first[0].to_string();
    let family = first[1].to_string();
    let dir = first[2].chars().next().unwrap_or('+');
    let mut div: f64 = first[3].parse().unwrap_or(0.0);
    let mut ins: f64 = first[4].parse().unwrap_or(0.0);
    let mut del: f64 = first[5].parse().unwrap_or(0.0);
    let mut rini: u32 = first[6].parse().unwrap_or(0);
    let mut rend: u32 = first[7].parse().unwrap_or(0);

    for frag in &frags[1..] {
        let f: Vec<&str> = frag.split(':').collect();
        if f.len() < 5 {
            continue;
        }
        let fd: f64 = f[0].parse().unwrap_or(0.0);
        let fi: f64 = f[1].parse().unwrap_or(0.0);
        let fd2: f64 = f[2].parse().unwrap_or(0.0);
        let fri: u32 = f[3].parse().unwrap_or(0);
        let fre: u32 = f[4].parse().unwrap_or(0);
        if fd > div { div = fd; }
        if fi > ins { ins = fi; }
        if fd2 > del { del = fd2; }
        if fri < rini { rini = fri; }
        if fre > rend { rend = fre; }
    }

    let length = rend.saturating_sub(rini);
    let nfrg = frags.len() as u32;
    Some(TeRepeat { name, family, dir, div, ins, del, length, nfrg })
}

// ---------------------------------------------------------------------------
// TRF profiler
// ---------------------------------------------------------------------------

/// Parse a UCSC TRF BED file and populate `model` with simple/tandem repeat records.
///
/// Expected columns (whitespace-delimited):
///   chrom(0) start(1) end(2) "trf"(3) period(4) copies(5) ?(6)
///   percent_match(7) percent_indel(8) ... consensus(last)
///
/// `te_intervals`: when provided, tandem repeats that are *fully contained
/// within* a TE interval are excluded (they are already accounted for by the
/// TE profile and would distort the simple repeat model).
///
/// `gene_filter`: when provided, tandem repeats that *overlap* any gene interval
/// are excluded (intron filter — mirrors Perl's `checkGene` logic).
///
/// Mirrors Perl's `profileTRF`.
pub fn profile_trf(
    path: impl AsRef<Path>,
    chrom_gc_bins: &HashMap<String, Vec<u32>>,
    win_size: u32,
    model: &mut RepeatModel,
    te_intervals: Option<&IntervalMap>,
    gene_filter: Option<&IntervalMap>,
) {
    let reader = crate::io::open_reader(path.as_ref());

    let mut last_ini: i64 = -1;
    let mut last_end: i64 = -1;
    let mut last_consensus = String::new();
    let mut last_gc: Option<usize> = None;

    for raw in reader.lines() {
        let line = raw.expect("Read error");
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let f: Vec<&str> = line.split_ascii_whitespace().collect();
        if f.len() < 9 {
            continue;
        }
        if f[3] != "trf" {
            log::warn!("TRF file line not in expected format, skipping: {}", line);
            continue;
        }

        let seq_id = f[0];
        let ini: u32 = match f[1].parse() { Ok(v) => v, Err(_) => continue };
        let end: u32 = match f[2].parse() { Ok(v) => v, Err(_) => continue };
        // Column 6 is TRF's copy number (column 5 is the period).
        let copy_num: f64 = match f[5].parse() { Ok(v) => v, Err(_) => continue };
        let pct_match: f64 = match f[7].parse() { Ok(v) => v, Err(_) => continue };
        let indel: f64 = match f[8].parse() { Ok(v) => v, Err(_) => continue };
        let div = 100.0 - pct_match;
        let consensus = f[f.len() - 1];

        // Skip tandem repeats that are fully contained within a TE interval.
        if let Some(te_ivs) = te_intervals {
            if let Some(te_chrom) = te_ivs.get(seq_id) {
                if contained_in_any(te_chrom, ini, end) {
                    continue;
                }
            }
        }

        // Skip tandem repeats that overlap a gene region (intron filter).
        if let Some(gf) = gene_filter {
            if let Some(gene_ivs) = gf.get(seq_id) {
                if overlaps_any(gene_ivs, ini, end) {
                    continue;
                }
            }
        }

        let gc = match get_gc_bin(chrom_gc_bins, seq_id, ini, win_size) {
            Some(g) => g as usize,
            None => continue,
        };

        // Canonical direction: prefer the lexicographically smaller of
        // consensus and its reverse complement.
        let (canon_con, dir) = canonical_consensus(consensus);

        // Handle overlapping repeats: if this repeat overlaps the last one,
        // keep whichever has the longer consensus; on a tie the newer record
        // replaces the previous one (Perl's profileTRF).
        let ini_i = ini as i64;
        let end_i = end as i64;
        if ini_i >= last_ini && ini_i <= last_end {
            if canon_con.len() < last_consensus.len() {
                continue; // discard this one
            }
            // Replace the previous entry
            if let Some(prev_gc) = last_gc {
                model.simple_profile[prev_gc].pop();
            }
        }

        model.simple_profile[gc].push(SimpleRepeat {
            consensus: canon_con.clone(),
            dir,
            copy_num,
            div,
            indel,
        });

        last_ini = ini_i;
        last_end = end_i;
        last_consensus = canon_con;
        last_gc = Some(gc);
    }
}

/// Returns (canonical_consensus, strand) where canonical is the
/// lexicographically smaller of the consensus and its reverse complement.
fn canonical_consensus(seq: &str) -> (String, char) {
    let rc = reverse_complement(seq);
    if rc.as_str() < seq {
        (rc, '-')
    } else {
        (seq.to_string(), '+')
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            other => other,
        })
        .collect()
}
