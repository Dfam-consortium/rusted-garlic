/// Repeat profile model — built from RepeatMasker .align and TRF BED files.
///
/// Mirrors the Perl `profileRMAlign` / `profileTRF` / `parseRepeats` logic in
/// `createModel.pl`.  The profile is stored per GC bin (indexed 0..num_bins-1)
/// and serialised alongside the sequence model in the JSON model file.
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::Path;

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

/// One interspersed (TE) repeat record — the fields that `createFakeSequence`
/// needs to simulate an evolved copy.
///
/// Corresponds to a line in the `.repeats.W1000.data` file:
///   `name:family:dir:div:ins:del:length:nfrg`
#[derive(Debug, Clone, Serialize, Deserialize)]
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

/// One simple/tandem repeat record.
///
/// Corresponds to a `SIMPLE:...` line in the `.repeats.W1000.data` file:
///   `SIMPLE:consensus:dir:period:div:indel`
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimpleRepeat {
    pub consensus: String, // repeat unit e.g. "AGAT"
    pub dir: char,         // canonical strand after revcomp check
    pub period: u32,       // repeat unit length
    pub div: f64,          // divergence (100 - percent_match from TRF)
    pub indel: f64,        // indel rate from TRF
}

/// The full repeat profile for one genome, indexed by GC bin.
#[derive(Debug, Serialize, Deserialize)]
pub struct RepeatModel {
    /// `te_profile[gc_bin]` = all TE records assigned to that bin
    pub te_profile: Vec<Vec<TeRepeat>>,
    /// `simple_profile[gc_bin]` = all simple repeat records assigned to that bin
    pub simple_profile: Vec<Vec<SimpleRepeat>>,
    /// `insert_table[gc_bin][outer_type#family][inner_type#family]` = frequency
    pub insert_table: Vec<HashMap<String, HashMap<String, f64>>>,
}

impl RepeatModel {
    pub fn new(num_bins: usize) -> Self {
        RepeatModel {
            te_profile: vec![Vec::new(); num_bins],
            simple_profile: vec![Vec::new(); num_bins],
            insert_table: vec![HashMap::new(); num_bins],
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
// RM .align profiler
// ---------------------------------------------------------------------------

/// Parse a RepeatMasker `.align` (or `.out`) file and populate `model` with
/// TE repeat records and the nested-insert table.
///
/// Mirrors Perl's `profileRMAlign`.
pub fn profile_rm_align(
    path: impl AsRef<Path>,
    chrom_gc_bins: &HashMap<String, Vec<u32>>,
    win_size: u32,
    num_bins: usize,
    model: &mut RepeatModel,
) {
    let reader = crate::io::open_reader(path.as_ref());

    // `repdata`: accumulates per-repeat-ID data while streaming.
    // Key: repeat ID string (last field of .align line, e.g. "m_b1s502i1")
    // Value: accumulated record
    struct RepData {
        label: String,  // semicolon-joined fragment params
        ini: u32,       // earliest genomic start
        end: u32,       // latest genomic end
        seq_id: String,
        type_fam: String, // "AluYa5#SINE/Alu"
    }
    let mut repdata: HashMap<String, RepData> = HashMap::new();

    // For nested-insert detection we need to track the type#family seen at
    // each line number, and the last line where each repeat ID appeared.
    let mut line_rep_name: Vec<Option<String>> = Vec::new(); // index = line number
    let mut id_hash: HashMap<String, usize> = HashMap::new(); // rep_id -> line number

    // Temporary insert counts: [gc_bin][outer][inner] = count
    let mut insert_counts: Vec<HashMap<String, HashMap<String, u32>>> =
        vec![HashMap::new(); num_bins];

    let mut line_num: usize = 0;

    for raw in reader.lines() {
        let line = raw.expect("Read error");
        line_num += 1;
        line_rep_name.push(None);

        let trimmed = line.trim_start();
        if !looks_like_data_line(trimmed) {
            continue;
        }
        if trimmed.contains("Simple_repeat")
            || trimmed.contains("Low_complexity")
            || trimmed.contains("Unknown")
            || trimmed.contains("Satellite")
        {
            continue;
        }

        let f: Vec<&str> = trimmed.split_ascii_whitespace().collect();
        // Two variants of the .align header line:
        //   Positive strand (13 fields):  ... (left)  name#fam  rbegin  rend  (rleft)  id
        //   Negative strand (14 fields):  ... (left)  C  name#fam  (rleft)  rend  rbegin  id
        // Both need at least 13 fields.
        if f.len() < 13 {
            continue;
        }

        let div: f64 = f[1].parse().unwrap_or(0.0);
        let ins: f64 = f[2].parse().unwrap_or(0.0);
        let del: f64 = f[3].parse().unwrap_or(0.0);
        let seq_id = f[4];
        let ini: u32 = f[5].parse().unwrap_or(0);
        let end: u32 = f[6].parse().unwrap_or(0);

        // Parse strand + name + family + consensus coords.
        //
        // RepeatMasker .align header format:
        //   Positive strand (13 fields):
        //     score div del ins seqid begin end (left)  name#fam  rbegin rend (rleft) id
        //     [0]   [1] [2] [3]  [4]  [5]  [6]   [7]     [8]      [9]  [10]   [11]  [12]
        //
        //   Negative strand (14 fields) — 'C' is inserted at [8], coords reversed:
        //     score div del ins seqid begin end (left) C name#fam (rleft) rend rbegin id
        //     [0]   [1] [2] [3]  [4]  [5]  [6]   [7] [8]  [9]    [10]  [11]  [12]  [13]
        //
        // Note: positive strand has no explicit strand character; 'C' marks negative.
        let (dir, type_str, fam_str, rini, rend) = if f[8] == "C" {
            // Negative strand — need at least 14 fields
            if f.len() < 14 { continue; }
            let (t, fa) = split_name_fam(f[9]);
            // f[10] = (rleft) [parenthesized, parse fails → ignored]
            let ri: u32 = f[12].parse().unwrap_or(0); // rbegin (leftmost in consensus)
            let re: u32 = f[11].parse().unwrap_or(0); // rend   (rightmost in consensus)
            ('-', t, fa, ri, re)
        } else {
            // Positive strand — name#fam is at [8], no explicit strand character
            let (t, fa) = split_name_fam(f[8]);
            let ri: u32 = f[9].parse().unwrap_or(0);  // rbegin
            let re: u32 = f[10].parse().unwrap_or(0); // rend
            ('+', t, fa, ri, re)
        };

        let type_fam = format!("{}#{}", type_str, fam_str);
        let rep_id = f[f.len() - 1].to_string();
        let fragment_label = format!(
            "{}:{}:{}:{:.2}:{:.2}:{:.2}:{}:{}",
            type_str, fam_str, dir, div, ins, del, rini, rend
        );

        // Update the name at this line for nested-insert tracking
        line_rep_name[line_num - 1] = Some(type_fam.clone());

        // Nested insert detection: if we've seen this rep_id before, the lines
        // between the two occurrences are elements inserted inside this repeat.
        if let Some(&prev_line) = id_hash.get(&rep_id) {
            if let Some(rd) = repdata.get(&rep_id) {
                let pos = rd.ini + (end.saturating_sub(rd.ini)) / 2;
                let gc = get_gc_bin(chrom_gc_bins, rd.seq_id.as_str(), pos, win_size);
                if let Some(gc_bin) = gc {
                    let outer = rd.type_fam.clone();
                    for ln in (prev_line + 1)..line_num.saturating_sub(1) {
                        if let Some(Some(inner)) = line_rep_name.get(ln) {
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
        id_hash.insert(rep_id.clone(), line_num - 1);

        // Accumulate or create repeat record
        repdata
            .entry(rep_id)
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

    // Assign each accumulated repeat to a GC bin and parse into TeRepeat
    for rd in repdata.values() {
        let pos = rd.ini + (rd.end.saturating_sub(rd.ini)) / 2;
        let gc = match get_gc_bin(chrom_gc_bins, &rd.seq_id, pos, win_size) {
            Some(g) => g as usize,
            None => continue,
        };
        if let Some(te) = parse_te_label(&rd.label) {
            model.te_profile[gc].push(te);
        }
    }

    // Convert insert counts to frequencies
    for (gc_bin, outer_map) in insert_counts.into_iter().enumerate() {
        for (outer, inner_map) in outer_map {
            let total: u32 = inner_map.values().sum();
            if total == 0 {
                continue;
            }
            let freq_map: HashMap<String, f64> = inner_map
                .into_iter()
                .map(|(k, v)| (k, v as f64 / total as f64))
                .collect();
            model.insert_table[gc_bin].insert(outer, freq_map);
        }
    }
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

/// Split "Name#Class/Family" into (name, family).
fn split_name_fam(s: &str) -> (&str, &str) {
    if let Some(idx) = s.find('#') {
        (&s[..idx], &s[idx + 1..])
    } else {
        (s, "Unknown")
    }
}

fn looks_like_data_line(s: &str) -> bool {
    let mut chars = s.chars().peekable();
    if !chars.peek().map_or(false, |c| c.is_ascii_digit()) {
        return false;
    }
    while chars.peek().map_or(false, |c| c.is_ascii_digit()) { chars.next(); }
    if !chars.peek().map_or(false, |c| c.is_ascii_whitespace()) { return false; }
    while chars.peek().map_or(false, |c| c.is_ascii_whitespace()) { chars.next(); }
    if !chars.peek().map_or(false, |c| c.is_ascii_digit()) { return false; }
    while chars.peek().map_or(false, |c| c.is_ascii_digit()) { chars.next(); }
    chars.peek() == Some(&'.')
}

// ---------------------------------------------------------------------------
// TRF profiler
// ---------------------------------------------------------------------------

/// Parse a GARLIC-format TRF BED file and populate `model` with simple repeat records.
///
/// Expected columns (whitespace-delimited):
///   chrom(0) start(1) end(2) "trf"(3) period(4) copies(5) ?(6)
///   percent_match(7) percent_indel(8) ... consensus(last)
///
/// Mirrors Perl's `profileTRF`.
pub fn profile_trf(
    path: impl AsRef<Path>,
    chrom_gc_bins: &HashMap<String, Vec<u32>>,
    win_size: u32,
    model: &mut RepeatModel,
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
            eprintln!("warning: TRF file line not in expected format, skipping: {}", line);
            continue;
        }

        let seq_id = f[0];
        let ini: u32 = match f[1].parse() { Ok(v) => v, Err(_) => continue };
        let end: u32 = match f[2].parse() { Ok(v) => v, Err(_) => continue };
        let period: u32 = match f[5].parse::<f64>() { Ok(v) => v as u32, Err(_) => continue };
        let pct_match: f64 = match f[7].parse() { Ok(v) => v, Err(_) => continue };
        let indel: f64 = match f[8].parse() { Ok(v) => v, Err(_) => continue };
        let div = 100.0 - pct_match;
        let consensus = f[f.len() - 1];

        let gc = match get_gc_bin(chrom_gc_bins, seq_id, ini, win_size) {
            Some(g) => g as usize,
            None => continue,
        };

        // Canonical direction: prefer the lexicographically smaller of
        // consensus and its reverse complement.
        let (canon_con, dir) = canonical_consensus(consensus);

        // Handle overlapping repeats: if this repeat overlaps the last one,
        // keep whichever has the longer consensus (prefer more informative).
        let ini_i = ini as i64;
        let end_i = end as i64;
        if ini_i >= last_ini && ini_i <= last_end {
            if canon_con.len() <= last_consensus.len() {
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
            period,
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

// ---------------------------------------------------------------------------
// Text file output (Perl-compatible)
// ---------------------------------------------------------------------------

/// Write a `.repeats.W1000.data` compatible text file for interoperability
/// with the Perl `createFakeSequence.pl`.
pub fn write_repeats_file(
    path: impl AsRef<Path>,
    model: &RepeatModel,
    gc_range_strings: &[String],
) {
    use std::io::Write;
    let mut f = File::create(path.as_ref())
        .unwrap_or_else(|e| panic!("Cannot create repeats file {:?}: {}", path.as_ref(), e));

    for (gc_bin, gc_str) in gc_range_strings.iter().enumerate() {
        writeln!(f, "#GC={}", gc_str).expect("Write error");
        for te in &model.te_profile[gc_bin] {
            writeln!(
                f,
                "{}:{}:{}:{:.2}:{:.2}:{:.2}:{}:{}",
                te.name, te.family, te.dir, te.div, te.ins, te.del, te.length, te.nfrg
            )
            .expect("Write error");
        }
        for sr in &model.simple_profile[gc_bin] {
            writeln!(
                f,
                "SIMPLE:{}:{}:{}:{:.2}:{:.2}",
                sr.consensus, sr.dir, sr.period, sr.div, sr.indel
            )
            .expect("Write error");
        }
        writeln!(f).expect("Write error");
    }
}

/// Write a `.inserts.W1000.data` compatible text file.
pub fn write_inserts_file(
    path: impl AsRef<Path>,
    model: &RepeatModel,
    gc_range_strings: &[String],
) {
    use std::io::Write;
    let mut f = File::create(path.as_ref())
        .unwrap_or_else(|e| panic!("Cannot create inserts file {:?}: {}", path.as_ref(), e));

    for (gc_bin, gc_str) in gc_range_strings.iter().enumerate() {
        writeln!(f, "#GC={}", gc_str).expect("Write error");
        let mut outer_keys: Vec<&String> = model.insert_table[gc_bin].keys().collect();
        outer_keys.sort();
        for outer in outer_keys {
            let inner_map = &model.insert_table[gc_bin][outer];
            let mut inner_keys: Vec<&String> = inner_map.keys().collect();
            inner_keys.sort();
            for inner in inner_keys {
                let frq = inner_map[inner];
                writeln!(f, "{}\t{}\t{:.8}", outer, inner, frq).expect("Write error");
            }
        }
    }
}
