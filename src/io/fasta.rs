/// Streaming FASTA reader.
///
/// Yields one `(name, sequence)` pair per record.  The sequence is returned
/// as an uppercase `Vec<u8>`.  Supports both plain and gzip-compressed files
/// via the shared `open_reader` helper.
use std::io::BufRead;
use std::path::Path;

pub struct FastaIter {
    reader: Box<dyn BufRead>,
    /// Name of the next record, read ahead when we hit a `>` line while
    /// assembling the previous record.
    pending_name: Option<String>,
}

impl FastaIter {
    pub fn new(path: &Path) -> Self {
        FastaIter {
            reader: super::open_reader(path),
            pending_name: None,
        }
    }
}

impl Iterator for FastaIter {
    type Item = (String, Vec<u8>);

    fn next(&mut self) -> Option<Self::Item> {
        let mut seq: Vec<u8> = Vec::new();
        let mut name: Option<String> = self.pending_name.take();

        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(0) | Err(_) => break,
                Ok(_) => {}
            }
            let trimmed = line.trim_end_matches(['\n', '\r', ' ']);
            if let Some(rest) = trimmed.strip_prefix('>') {
                let new_name = rest.trim().to_string();
                if name.is_some() {
                    self.pending_name = Some(new_name);
                    break;
                } else {
                    name = Some(new_name);
                }
            } else if name.is_some() {
                seq.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
            }
        }

        name.map(|n| (n, seq))
    }
}
