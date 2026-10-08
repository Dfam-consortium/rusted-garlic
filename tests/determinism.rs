//! Tier 3 regression tests: determinism and golden outputs.
//!
//! - `build` must produce a byte-identical model on every run.
//! - `emit --seed N` must produce byte-identical FASTA and BED on every run.
//! - Both must match the golden files in `tests/data/golden/`, so any change to
//!   model building or sequence generation shows up as a test failure.
//!
//! When a change in output is intended, regenerate the goldens and commit them:
//!
//!     GARLIC_UPDATE_GOLDEN=1 cargo test --release --test determinism
//!
//! Fixtures are described in `tests/data/README.md`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_rusted-garlic");

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn golden(name: &str) -> PathBuf {
    data("golden").join(name)
}

/// A fresh scratch directory per test (tests run in parallel).
fn scratch(test: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("determinism").join(test);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(args: &[&str]) {
    let out = Command::new(BIN)
        .arg("--log-level")
        .arg("error")
        .args(args)
        .output()
        .expect("failed to run rusted-garlic");
    assert!(
        out.status.success(),
        "rusted-garlic {:?} failed:\n{}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
}

fn build(model: &Path) {
    build_with(model, &[]);
}

fn build_with(model: &Path, extra: &[&str]) {
    let (fa, te, trf) = (data("sm-example.fa"), data("sm-example.te_annot.tsv"), data("sm-example.fa.trf.bed"));
    run(&[
        "build",
        "--in_seq", fa.to_str().unwrap(),
        "--te_annot", te.to_str().unwrap(),
        "--tandem_annot", trf.to_str().unwrap(),
        "--model", model.to_str().unwrap(),
    ].into_iter().chain(extra.iter().copied()).collect::<Vec<_>>());
}

/// Emit 2 x 20 kb with TE and simple-repeat bombardment; returns (fasta, bed).
fn emit(model: &Path, dir: &Path, tag: &str, seed: u64, perl_compat: bool) -> (Vec<u8>, Vec<u8>) {
    let fa = dir.join(format!("{tag}.fa"));
    let bed = dir.join(format!("{tag}.bed"));
    let seed = seed.to_string();
    let cons = data("sm-example.cons.fa");
    let mut args = vec![
        "emit",
        "--model", model.to_str().unwrap(),
        "--out_seq", fa.to_str().unwrap(),
        "--out_bed", bed.to_str().unwrap(),
        "--seq_size", "20k",
        "--num_seqs", "2",
        "--te_cons", cons.to_str().unwrap(),
        "--te_frac", "40",
        "--sim_frac", "2",
        "--seed", &seed,
    ];
    if perl_compat {
        args.push("--perl_compat");
    }
    run(&args);
    (fs::read(fa).unwrap(), fs::read(bed).unwrap())
}

/// Compare `actual` to a golden file, or rewrite it when GARLIC_UPDATE_GOLDEN is set.
fn check_golden(name: &str, actual: &[u8]) {
    let path = golden(name);
    if std::env::var_os("GARLIC_UPDATE_GOLDEN").is_some() {
        fs::write(&path, actual).unwrap();
        return;
    }
    let expected = fs::read(&path).unwrap_or_else(|e| {
        panic!("missing golden {path:?} ({e}); run with GARLIC_UPDATE_GOLDEN=1 to create it")
    });
    assert!(
        expected == actual,
        "{name} differs from golden {path:?}. If the change is intended, rerun with \
         GARLIC_UPDATE_GOLDEN=1 and commit the new golden files."
    );
}

#[test]
fn build_is_deterministic() {
    let dir = scratch("build_is_deterministic");
    let (m1, m2) = (dir.join("a.model"), dir.join("b.model"));
    build(&m1);
    build(&m2);
    assert_eq!(fs::read(&m1).unwrap(), fs::read(&m2).unwrap(), "two builds differ");
}

#[test]
fn build_matches_golden() {
    let dir = scratch("build_matches_golden");
    let m = dir.join("sm-example.model");
    build(&m);
    check_golden("sm-example.model", &fs::read(&m).unwrap());
}

#[test]
fn build_perl_compat_matches_golden() {
    let dir = scratch("build_perl_compat_matches_golden");
    let m = dir.join("sm-example.perl_compat.model");
    build_with(&m, &["--perl_compat"]);
    check_golden("sm-example.perl_compat.model", &fs::read(&m).unwrap());
}

#[test]
fn emit_same_seed_is_identical() {
    let dir = scratch("emit_same_seed_is_identical");
    let model = golden("sm-example.model");
    for compat in [false, true] {
        let a = emit(&model, &dir, &format!("a{compat}"), 7, compat);
        let b = emit(&model, &dir, &format!("b{compat}"), 7, compat);
        assert!(a.0 == b.0, "FASTA differs for same seed (perl_compat={compat})");
        assert!(a.1 == b.1, "BED differs for same seed (perl_compat={compat})");
    }
}

#[test]
fn emit_different_seed_differs() {
    let dir = scratch("emit_different_seed_differs");
    let model = golden("sm-example.model");
    let a = emit(&model, &dir, "a", 7, false);
    let b = emit(&model, &dir, "b", 8, false);
    assert!(a.0 != b.0, "different seeds gave the same FASTA");
}

#[test]
fn emit_matches_golden() {
    let dir = scratch("emit_matches_golden");
    let model = golden("sm-example.model");
    let (fa, bed) = emit(&model, &dir, "default", 1, false);
    check_golden("emit_seed1.fa", &fa);
    check_golden("emit_seed1.bed", &bed);
}

#[test]
fn emit_perl_compat_matches_golden() {
    let dir = scratch("emit_perl_compat_matches_golden");
    let model = golden("sm-example.model");
    let (fa, bed) = emit(&model, &dir, "compat", 1, true);
    check_golden("emit_perl_compat_seed1.fa", &fa);
    check_golden("emit_perl_compat_seed1.bed", &bed);
}
