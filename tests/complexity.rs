//! `complexity` must reproduce SeqComplex.pm, the Perl module it replaces.
//!
//! `tests/data/complexity/windows.fa` holds 199 one-window records (real chr22 windows,
//! some with N; windows of a simulated sequence; low-complexity extremes; 50, 7 and 1 bp
//! records). `seqcomplex_reference.tsv` is SeqComplex.pm's output for each record (its
//! per-measure subroutines, as used by profileComplexSeq.pl). See
//! `tests/data/README.md` and `validation/complexity.sh` for how they were made.
//!
//! Every measure must match to Rust's 6 printed decimals, except `cz`: Perl gzips a temp
//! file whose random name ends up in the gzip header, so its compression ratio is a few
//! percent low; Rust compresses in memory. `NA` in the reference (records SeqComplex.pm
//! cannot compute, e.g. the 1-bp record) is skipped.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_rusted-garlic");
const TOL: f64 = 1e-6;
const SKIP: &[&str] = &["cz"];

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/complexity").join(name)
}

/// Parse a TSV with a header into {seq: {column: value}}.
fn parse_tsv(text: &str) -> HashMap<String, HashMap<String, String>> {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("empty TSV").split('\t').collect();
    let seq_col = header.iter().position(|&c| c == "seq").expect("no 'seq' column");
    lines
        .filter(|l| !l.is_empty())
        .map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            let row = header
                .iter()
                .zip(&fields)
                .map(|(h, v)| (h.to_string(), v.to_string()))
                .collect();
            (fields[seq_col].to_string(), row)
        })
        .collect()
}

#[test]
fn complexity_matches_seqcomplex() {
    let out = Command::new(BIN)
        .args(["--log-level", "error", "complexity"])
        .arg(data("windows.fa"))
        .output()
        .expect("failed to run rusted-garlic");
    assert!(out.status.success(), "complexity failed:\n{}", String::from_utf8_lossy(&out.stderr));
    let rust = parse_tsv(&String::from_utf8(out.stdout).unwrap());

    let reference_text = std::fs::read_to_string(data("seqcomplex_reference.tsv")).unwrap();
    let measures: Vec<String> = reference_text
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .skip(1)
        .map(String::from)
        .collect();
    let reference = parse_tsv(&reference_text);

    assert_eq!(rust.len(), reference.len(), "record count differs from the reference");

    let mut failures = Vec::new();
    let mut compared = 0;
    for (seq, perl_row) in &reference {
        let rust_row = rust.get(seq).unwrap_or_else(|| panic!("record {seq} missing from Rust output"));
        for m in &measures {
            if SKIP.contains(&m.as_str()) || perl_row[m] == "NA" {
                continue;
            }
            let p: f64 = perl_row[m].parse().unwrap();
            let r: f64 = rust_row
                .get(m)
                .unwrap_or_else(|| panic!("measure {m} missing from Rust output"))
                .parse()
                .unwrap();
            compared += 1;
            if (r - p).abs() > TOL {
                failures.push(format!("{seq} {m}: rust {r} perl {p}"));
            }
        }
    }
    assert!(compared > 5000, "only {compared} values compared; reference or output malformed?");
    assert!(
        failures.is_empty(),
        "{} of {compared} values differ from SeqComplex.pm:\n{}",
        failures.len(),
        failures.iter().take(20).cloned().collect::<Vec<_>>().join("\n")
    );
}
