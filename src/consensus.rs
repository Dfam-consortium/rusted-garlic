/// Load TE consensus sequences from a FASTA file.
///
/// Resolves IUPAC ambiguity codes to random concrete bases (matching Perl's
/// `checkBases`).  Unknown codes default to 'A'.
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use rand::prelude::*;

/// Load a FASTA file into a name → uppercase-DNA map.
/// Only the first word of each FASTA header is used as the key.
/// Transparently decompresses `.gz` files.
pub fn load_fasta(path: impl AsRef<Path>, rng: &mut impl Rng) -> HashMap<String, Vec<u8>> {
    let reader = crate::io::open_reader(path.as_ref());

    let mut map: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Read error in consensus FASTA");
        let trimmed = line.trim();
        if trimmed.starts_with('>') {
            if !current_name.is_empty() {
                map.insert(current_name.clone(), std::mem::take(&mut current_seq));
            }
            // Use the first whitespace-delimited word of the header, then
            // strip any trailing "#CLASS/FAMILY" suffix (Dfam-format FASTAs
            // use ">AluYa5#SINE/Alu" while the model stores just "AluYa5").
            current_name = trimmed[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .split('#')
                .next()
                .unwrap_or("")
                .to_string();
            current_seq = Vec::new();
        } else if !trimmed.is_empty() && !trimmed.starts_with('#') {
            for b in trimmed.bytes() {
                let resolved = resolve_iupac(b, rng);
                current_seq.push(resolved);
            }
        }
    }
    if !current_name.is_empty() {
        map.insert(current_name, current_seq);
    }
    map
}

/// Resolve one IUPAC base byte to a concrete ACGT byte (uppercase).
/// Matches Perl's `checkBases` sub.
fn resolve_iupac(b: u8, rng: &mut impl Rng) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'A',
        b'C' => b'C',
        b'G' => b'G',
        b'T' | b'U' => b'T',
        b'R' => *[b'A', b'G'].choose(rng).unwrap(),
        b'Y' => *[b'C', b'T'].choose(rng).unwrap(),
        b'M' => *[b'A', b'C'].choose(rng).unwrap(),
        b'K' => *[b'G', b'T'].choose(rng).unwrap(),
        b'W' => *[b'A', b'T'].choose(rng).unwrap(),
        b'S' => *[b'C', b'G'].choose(rng).unwrap(),
        b'B' => *[b'C', b'T', b'G'].choose(rng).unwrap(),
        b'D' => *[b'A', b'T', b'G'].choose(rng).unwrap(),
        b'H' => *[b'A', b'T', b'C'].choose(rng).unwrap(),
        b'V' => *[b'A', b'C', b'G'].choose(rng).unwrap(),
        _ => b'A',
    }
}
