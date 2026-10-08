/// Sequence evolution and repeat bombardment.
///
/// Mirrors the Perl `createFakeSequence.pl` functions:
///   `addTransitions`, `addTransversions`, `addInsertions`, `addDeletions`,
///   `evolveRepeat`, `evolveSimple`, `insertRepeat`, `insertLowComplex`.
///
/// Key design choices:
///   - Mutations are accepted only when they make the local k-mer context
///     *more* compatible with the background model (matching Perl's eval loop).
///   - Nested insertions are supported via the insert_table from the model.
///   - Inserted sequences are lower-cased; the base sequence is upper-case.
use std::collections::HashMap;

use rand::prelude::*;

use crate::repeat_model::{RepeatModel, SimpleRepeat, TeRepeat};
use crate::seq_model::SeqModel;

// ---------------------------------------------------------------------------
// Insertion record (for BED output)
// ---------------------------------------------------------------------------

/// One inserted repeat element — enough to write a BED line.
#[derive(Debug, Clone)]
pub struct InsertRecord {
    /// 0-based start in the final sequence.
    pub start: usize,
    /// 0-based exclusive end in the final sequence.
    pub end: usize,
    /// Repeat name (e.g. "AluYa5").
    pub name: String,
    /// Repeat family (e.g. "SINE/Alu").  Empty string for simple repeats.
    pub family: String,
    /// Strand: '+' or '-'.
    pub strand: char,
}

/// Options for repeat evolution.
#[derive(Debug, Clone, Copy)]
pub struct EvolveOpts {
    /// Reproduce GARLIC 1.4 behaviour (see `--perl_compat`); also disables
    /// the CpG rate.
    pub perl_compat: bool,
    /// Transition rate at CpG sites relative to other sites (1 = no CpG effect).
    pub cpg_rate: f64,
}

impl EvolveOpts {
    /// CpG rate in effect: 1 under `perl_compat` (GARLIC 1.4 has no CpG effect).
    fn effective_cpg_rate(&self) -> f64 {
        if self.perl_compat { 1.0 } else { self.cpg_rate }
    }
}

// Maximum number of k-mer-guided mutation rounds (Perl's $mut_cyc).
const MUT_CYCLES: usize = 10;
// Maximum insertion-attempt cycles (Perl's $ins_cyc).
const INS_CYCLES: usize = 100;

// ---------------------------------------------------------------------------
// Helper: kmer frequency lookup
// ---------------------------------------------------------------------------

/// Look up P(suffix | prefix) in the given GC bin from the model's kmer table.
/// `kmer` must be exactly `model.kmer_size` bytes of ASCII DNA (uppercase or lower).
/// Returns 0.0 when the prefix is not in the model.
pub fn kmer_freq(model: &SeqModel, gc_bin: usize, kmer: &[u8]) -> f64 {
    let kmer_size = model.kmer_size;
    if kmer.len() != kmer_size {
        return 0.0;
    }
    let suffix_idx = match kmer[kmer_size - 1].to_ascii_uppercase() {
        b'A' => 0usize,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => return 0.0,
    };
    // Build an uppercase prefix string for the lookup.
    let prefix: String = kmer[..kmer_size - 1]
        .iter()
        .map(|b| b.to_ascii_uppercase() as char)
        .collect();
    model
        .gc_kmer_freq(gc_bin, &prefix)
        .map(|freqs| freqs[suffix_idx])
        .unwrap_or(0.0)
}

// ---------------------------------------------------------------------------
// Helper: reverse complement
// ---------------------------------------------------------------------------

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Position selection (Perl's selPosition / randSel)
// ---------------------------------------------------------------------------

/// Return a randomly-ordered list of approximately len/2 unique positions from
/// 0..len-1.  Mirrors Perl's `selPosition` → `randSel`.
fn sel_positions(len: usize, rng: &mut impl Rng) -> Vec<usize> {
    let want = (len / 2).saturating_add(2).min(len);
    // Partial Fisher-Yates: shuffle the first `want` elements in place.
    let mut positions: Vec<usize> = (0..len).collect();
    for i in 0..want {
        let j = rng.gen_range(i..len);
        positions.swap(i, j);
    }
    positions.truncate(want);
    positions
}

// ---------------------------------------------------------------------------
// Mutation operators
// ---------------------------------------------------------------------------

/// Whether a mutation that changes the frequencies of the k-mers it touches
/// from `old` to `new` is acceptable: none may get less likely (GARLIC 1.4,
/// `next if old >= new + 1e-15`).  Neutral changes are accepted; rejecting
/// them too was tried and made repeat 4-mer spectra less like real TEs
/// (PERL_VS_RUST_DIFFS.md, 2026-10-07).
fn improves(old: &[f64], new: &[f64]) -> bool {
    old.iter().zip(new).all(|(o, n)| *o < n + 1e-15)
}

/// Apply up to `n` A↔G / T↔C transitions to `seq`.
///
/// Each candidate position is accepted only when neither k-mer the transition
/// touches becomes less likely under the background model (see `improves`).
/// With `eval >= MUT_CYCLES` the check is skipped (unused by `mutate_guided`).
///
/// Returns the number of transitions *not* applied (skipped).
fn add_transitions(
    seq: &mut Vec<u8>,
    mut n: usize,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    eval: usize,
    _perl_compat: bool,
) -> usize {
    let kmer_size = model.kmer_size;
    let len = seq.len();
    for pos in sel_positions(len, rng) {
        if n == 0 {
            break;
        }
        // Perl guard: skip the first kmer positions and the last position.
        if pos < kmer_size || pos >= len {
            continue;
        }
        let old = seq[pos - 1];
        let new = match old.to_ascii_uppercase() {
            b'A' => b'G',
            b'G' => b'A',
            b'C' => b'T',
            b'T' => b'C',
            _ => continue,
        };

        if eval < MUT_CYCLES {
            // pre context: kmer ending at pos-1
            let pre_old = &seq[pos - kmer_size..pos];
            let mut pre_new = pre_old.to_vec();
            pre_new[kmer_size - 1] = new;
            // post context: kmer starting at pos-(kmer_size-1), including seq[pos]
            let post_old = &seq[pos - kmer_size + 1..pos + 1];
            let mut post_new = post_old.to_vec();
            post_new[kmer_size - 2] = new; // mutation is at index kmer_size-2 in post slice
            let old_f = [kmer_freq(model, gc_bin, pre_old), kmer_freq(model, gc_bin, post_old)];
            let new_f = [kmer_freq(model, gc_bin, &pre_new), kmer_freq(model, gc_bin, &post_new)];
            if !improves(&old_f, &new_f) {
                continue;
            }
        }
        seq[pos - 1] = new;
        n -= 1;
    }
    n // skipped count
}

/// Apply up to `n` A↔T / C↔G transversions to `seq`.
fn add_transversions(
    seq: &mut Vec<u8>,
    mut n: usize,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    eval: usize,
    _perl_compat: bool,
) -> usize {
    let kmer_size = model.kmer_size;
    let len = seq.len();
    for pos in sel_positions(len, rng) {
        if n == 0 {
            break;
        }
        if pos < kmer_size || pos >= len {
            continue;
        }
        let old = seq[pos - 1];
        let new = match old.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => continue,
        };

        if eval < MUT_CYCLES {
            let pre_old = &seq[pos - kmer_size..pos];
            let mut pre_new = pre_old.to_vec();
            pre_new[kmer_size - 1] = new;
            let post_old = &seq[pos - kmer_size + 1..pos + 1];
            let mut post_new = post_old.to_vec();
            post_new[kmer_size - 2] = new;
            let old_f = [kmer_freq(model, gc_bin, pre_old), kmer_freq(model, gc_bin, post_old)];
            let new_f = [kmer_freq(model, gc_bin, &pre_new), kmer_freq(model, gc_bin, &post_new)];
            if !improves(&old_f, &new_f) {
                continue;
            }
        }
        seq[pos - 1] = new;
        n -= 1;
    }
    n
}

/// Place CpG transitions in `seq` and return how many of the `n_sit`
/// requested transitions they used.
///
/// The element's transitions are spread so that a site in a CpG dinucleotide
/// mutates `cpg_rate` times as often as any other site: with `L` sites, `c`
/// CpGs and per-site rate `r = n_sit / (L + (cpg_rate - 1) * 2c)`, the CpGs
/// receive `cpg_rate * r * 2c` transitions (at most one per CpG, since the
/// first destroys it).  Each turns CG into TG or CA (C->T or G->A, equally
/// likely).  CpG decay is a mutation bias, not selection, so these are not
/// checked against the k-mer model; the remaining transitions are.
fn add_cpg_transitions(seq: &mut [u8], n_sit: usize, cpg_rate: f64, rng: &mut impl Rng) -> usize {
    if cpg_rate <= 1.0 || n_sit == 0 {
        return 0;
    }
    let mut cpgs: Vec<usize> = (0..seq.len().saturating_sub(1))
        .filter(|&i| seq[i].to_ascii_uppercase() == b'C' && seq[i + 1].to_ascii_uppercase() == b'G')
        .collect();
    if cpgs.is_empty() {
        return 0;
    }
    let (l, c) = (seq.len() as f64, cpgs.len() as f64);
    let r = n_sit as f64 / (l + (cpg_rate - 1.0) * 2.0 * c);
    let expected = (cpg_rate * r * 2.0 * c).min(c).min(n_sit as f64);
    // Stochastic rounding keeps the expected count exact.
    let mut k = expected.floor() as usize;
    if rng.gen::<f64>() < expected.fract() {
        k += 1;
    }
    for i in 0..k {
        let j = rng.gen_range(i..cpgs.len());
        cpgs.swap(i, j);
        let pos = cpgs[i];
        if rng.gen::<bool>() {
            seq[pos] = b'T';
        } else {
            seq[pos + 1] = b'A';
        }
    }
    k
}

/// Signature shared by the guided point-mutation operators.
type MutationOp<R> = fn(&mut Vec<u8>, usize, usize, &SeqModel, &mut R, usize, bool) -> usize;

/// Apply `n` mutations with Perl's retry semantics (`addTransitions` et al.).
///
/// Each round tries the remaining count with k-mer-guided acceptance.  Rejected
/// mutations are retried in the next round against a newly drawn GC bin
/// (Perl's `newGC()`), for at most `MUT_CYCLES` rounds; whatever is still
/// rejected after that is dropped, so the realised divergence can fall short of
/// the target.  Returns the number of mutations dropped.
fn mutate_guided<R: Rng>(
    op: MutationOp<R>,
    seq: &mut Vec<u8>,
    mut n: usize,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut R,
    perl_compat: bool,
) -> usize {
    let mut bin = gc_bin;
    for eval in 0..MUT_CYCLES {
        if n == 0 {
            break;
        }
        n = op(seq, n, bin, model, rng, eval, perl_compat);
        if n > 0 {
            bin = model.sample_gc_bin(rng);
        }
    }
    n
}

/// Apply up to `n` single-base insertions to `seq`, growing the vector.
///
/// Mirrors Perl's `addInsertions`: one pass over ~len/2 random positions.  At
/// each position a base is sampled from the k-mer model conditioned on the
/// (k-1)-mer ending just before `pos`, and kept when the dice roll is at least
/// the frequency of the downstream k-mer `seed[1..] + base + seq[pos]`.
///
/// The base is inserted at `pos`, between the seed and `seq[pos]` where it was
/// evaluated.  With `perl_compat` it is placed *after* `seq[pos]` instead, as
/// GARLIC 1.4 does (`substr($seq, $pos, 1) = "$old$n"`) — a Perl bug that leaves
/// inserted bases out of context and raises CpG content in repeats.
fn add_insertions(
    seq: &mut Vec<u8>,
    mut n: usize,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    perl_compat: bool,
) {
    let kmer_size = model.kmer_size;
    let dna: [u8; 4] = [b'A', b'C', b'G', b'T'];
    for pos in sel_positions(seq.len(), rng) {
        if n == 0 {
            break;
        }
        if pos < kmer_size || pos >= seq.len() {
            continue;
        }
        let seed = &seq[pos - kmer_size + 1..pos];
        let post = seq[pos];

        // Sample new base from kmer distribution conditioned on seed
        let dice: f64 = rng.gen();
        let mut chosen = dna[3];
        let mut p = 0.0f64;
        for &b in &dna {
            let kmer: Vec<u8> = seed.iter().copied().chain(std::iter::once(b)).collect();
            let q = p + kmer_freq(model, gc_bin, &kmer);
            if dice >= p {
                chosen = b;
            }
            if dice >= p && dice <= q + 1e-15 {
                break;
            }
            p = q;
        }

        // Downstream acceptance: seed[1..] + chosen + post
        let mut downstream: Vec<u8> = seed[1..].to_vec();
        downstream.push(chosen);
        downstream.push(post);
        if dice < kmer_freq(model, gc_bin, &downstream) + 1e-15 {
            continue;
        }

        seq.insert(if perl_compat { pos + 1 } else { pos }, chosen);
        n -= 1;
    }
}

/// Apply up to `n` single-base deletions to `seq`, shrinking the vector.
///
/// For a candidate base `d = seq[i]` (i = pos + k - 1) the k-mer ending at it,
/// `seq[i-k+1..=i]`, is compared with the k-mer that replaces it once `d` is
/// gone, `seq[i-k+1..i] + seq[i+1]`; the deletion is accepted unless the latter
/// is less likely (see `improves`).  Deletions are applied one at a time, so each
/// is judged on the current sequence.
///
/// With `perl_compat` this follows GARLIC 1.4 `addDeletions`, which scores the
/// same k-mers but deletes `seq[pos]` (k - 1 bases earlier than the base it
/// evaluated), with all candidates judged on the unmodified sequence.
fn add_deletions(
    seq: &mut Vec<u8>,
    mut n: usize,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    eval: usize,
    perl_compat: bool,
) -> usize {
    let kmer_size = model.kmer_size;
    let mut to_delete: Vec<usize> = Vec::new();
    for pos in sel_positions(seq.len(), rng) {
        if n == 0 {
            break;
        }
        if pos < kmer_size || pos + kmer_size >= seq.len() {
            continue;
        }
        if eval < MUT_CYCLES {
            let old_kmer = &seq[pos..pos + kmer_size];
            let mut new_kmer = seq[pos..pos + kmer_size - 1].to_vec();
            new_kmer.push(seq[pos + kmer_size]);
            let f_old = kmer_freq(model, gc_bin, old_kmer);
            let f_new = kmer_freq(model, gc_bin, &new_kmer);
            if !improves(&[f_old], &[f_new]) {
                continue;
            }
        }
        if perl_compat {
            to_delete.push(pos);
        } else {
            seq.remove(pos + kmer_size - 1);
        }
        n -= 1;
    }
    // Perl marks deletions and removes them afterwards: highest index first.
    to_delete.sort_unstable_by(|a, b| b.cmp(a));
    for i in to_delete {
        seq.remove(i);
    }
    n // skipped count
}

// ---------------------------------------------------------------------------
// Repeat evolution
// ---------------------------------------------------------------------------

/// TE instances per (GC bin, repeat name), sorted by age ascending — Perl's
/// `%byGCByRepeat`.  Used to pick nested inserts younger than their host.
pub struct AgeIndex<'a> {
    by_bin: Vec<HashMap<&'a str, Vec<(f64, &'a TeRepeat)>>>,
}

impl<'a> AgeIndex<'a> {
    pub fn new(rm: &'a RepeatModel) -> Self {
        let by_bin = rm
            .te_profile
            .iter()
            .map(|reps| {
                let mut m: HashMap<&str, Vec<(f64, &TeRepeat)>> = HashMap::new();
                for r in reps {
                    m.entry(r.name.as_str()).or_default().push((te_age(r, r.nfrg), r));
                }
                for v in m.values_mut() {
                    v.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
                }
                m
            })
            .collect();
        AgeIndex { by_bin }
    }
}

/// Perl's "age" of a TE record: div + ins + del + 10 * fragments.
fn te_age(rep: &TeRepeat, nfrg: u32) -> f64 {
    rep.div + rep.ins + rep.del + nfrg as f64 * 10.0
}

/// Evolve one interspersed TE repeat and return its mutated sequence.
///
/// Mirrors Perl's `evolveRepeat`.  Returns `None` when the consensus is not
/// found or the fragment is too short.
///
/// Fragmented records (`nfrg > 1`) receive up to `nfrg` younger nested
/// inserts drawn from the model's `insert_table` (see `get_insert`).
pub fn evolve_repeat(
    rep: &TeRepeat,
    consensus: &HashMap<String, Vec<u8>>,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    rm: &RepeatModel,
    ages: &AgeIndex,
    opts: EvolveOpts,
) -> Option<(Vec<u8>, char)> {
    evolve_repeat_with(rep, rep.nfrg, f64::INFINITY, consensus, gc_bin, model, rng, rm, ages, opts)
}

/// `evolve_repeat` with Perl's extra arguments: `nfrg` overrides the record's
/// fragment count (nested inserts use 1, so they are not nested further) and
/// records older than `max_age` are rejected.
fn evolve_repeat_with(
    rep: &TeRepeat,
    nfrg: u32,
    max_age: f64,
    consensus: &HashMap<String, Vec<u8>>,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    rm: &RepeatModel,
    ages: &AgeIndex,
    opts: EvolveOpts,
) -> Option<(Vec<u8>, char)> {
    let perl_compat = opts.perl_compat;
    let cons = consensus.get(&rep.name)?;
    if cons.len() < 10 {
        return None;
    }

    // Fragment length: cap at consensus length.
    let frag_len = (rep.length as usize).min(cons.len());
    if frag_len < 10 {
        return None;
    }

    let age = te_age(rep, nfrg);
    if age > max_age {
        return None;
    }

    // Pick a random start position in the consensus.
    let (ini, end) = if perl_compat && rep.length as usize > cons.len() {
        // GARLIC 1.4 bug: Perl only caps the fragment length for range values,
        // so `int(rand(len - frag))` gets a negative argument and `substr` takes a
        // window anchored m bases before the end (m uniform in 0..frag-len),
        // clipped to the consensus — usually a short 3' piece, not the full copy.
        let len = cons.len() as i64;
        let m = rng.gen_range(0..rep.length as i64 - len);
        let start = if m == 0 { 0 } else { len - m };
        let stop = start + rep.length as i64;
        (start.clamp(0, len) as usize, stop.clamp(0, len) as usize)
    } else {
        let max_start = cons.len().saturating_sub(frag_len);
        let ini = if max_start > 0 { rng.gen_range(0..max_start) } else { 0 };
        (ini, ini + frag_len)
    };
    let mut seq: Vec<u8> = cons[ini..end].to_vec();

    // Random orientation.
    let orient_minus = rng.gen::<f64>() > 0.5;
    let strand = if orient_minus { '-' } else { '+' };
    if orient_minus {
        seq = revcomp(&seq);
    }

    // Apply mutations.
    let mut_count = (rep.div * seq.len() as f64 / 100.0) as usize;
    let n_sit = mut_count / 2;
    let n_ver = mut_count - n_sit;
    let n_ins = (rep.ins * seq.len() as f64 / 100.0) as usize;
    let n_del = (rep.del * seq.len() as f64 / 100.0) as usize;

    let n_sit = n_sit - add_cpg_transitions(&mut seq, n_sit, opts.effective_cpg_rate(), rng);
    mutate_guided(add_transitions, &mut seq, n_sit, gc_bin, model, rng, perl_compat);
    mutate_guided(add_transversions, &mut seq, n_ver, gc_bin, model, rng, perl_compat);
    if n_ins > 0 {
        add_insertions(&mut seq, n_ins, gc_bin, model, rng, perl_compat);
    }
    mutate_guided(add_deletions, &mut seq, n_del, gc_bin, model, rng, perl_compat);

    // Nested insertions (Perl: up to `break` inserts, `ins_cyc + 1` tries).
    // Targets are drawn on the un-nested sequence and applied from the highest
    // position down so earlier targets stay valid.
    if nfrg > 1 && age > 10.0 && !seq.is_empty() {
        let outer_key = format!("{}#{}", rep.name, rep.family);
        let mut pending: Vec<(usize, Vec<u8>)> = Vec::new();
        let mut num = 1u32;
        for _ in 0..=INS_CYCLES {
            num += 1;
            match get_insert(&outer_key, age, consensus, gc_bin, model, rng, rm, ages, opts) {
                Some(ins) if !ins.is_empty() => {
                    let target = rng.gen_range(0..seq.len());
                    pending.push((target, ins));
                    if num > nfrg {
                        break;
                    }
                }
                _ => num -= 1,
            }
        }
        // Stable sort: inserts sharing a target keep Perl's stacking order.
        pending.sort_by(|x, y| y.0.cmp(&x.0));
        for (target, ins) in pending {
            let lc: Vec<u8> = ins.iter().map(|b| b.to_ascii_lowercase()).collect();
            seq.splice(target..target, lc);
        }
    }

    Some((seq, strand))
}

/// Draw one nested insert for a host repeat — Perl's `getInsert`.
///
/// The inner type is sampled from `insert_table[gc_bin][host]`; an instance of
/// it no older than max(host age - 20, 10) is chosen from the age-sorted list
/// and evolved without further nesting.  Up to `MUT_CYCLES` tries.
fn get_insert(
    outer_key: &str,
    age: f64,
    consensus: &HashMap<String, Vec<u8>>,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    rm: &RepeatModel,
    ages: &AgeIndex,
    opts: EvolveOpts,
) -> Option<Vec<u8>> {
    let new_age = (age - 20.0).max(10.0);
    let inner_map = rm.insert_table.get(gc_bin)?.get(outer_key)?;
    // Perl requires at least two candidate inner types (`last unless $#rep > 0`).
    if inner_map.len() < 2 {
        return None;
    }
    for _ in 0..MUT_CYCLES {
        let dice: f64 = rng.gen();
        let mut p = 0.0f64;
        let mut chosen: Option<&str> = None;
        for (inner_key, &freq) in inner_map.iter() {
            chosen = Some(inner_key.as_str());
            p += freq;
            if dice <= p {
                break;
            }
        }
        let inner_name = chosen?.split('#').next().unwrap_or("");
        let empty = Vec::new();
        let cands = ages.by_bin.get(gc_bin).and_then(|m| m.get(inner_name)).unwrap_or(&empty);

        // Perl's binary search for the youngest-enough prefix of the age-sorted
        // list, kept verbatim: an index of 0 counts as "none young enough".
        let (mut min_idx, mut max_idx) = (0isize, cands.len() as isize - 1);
        while max_idx - min_idx > 1 {
            let thresh = min_idx + (max_idx - min_idx + 1) / 2;
            if cands[thresh as usize].0 > new_age {
                max_idx = thresh - 1;
            } else {
                min_idx = thresh;
            }
        }
        if min_idx == 0 {
            return None;
        }
        let inner = cands[rng.gen_range(0..=min_idx as usize)].1;
        if let Some((s, _)) =
            evolve_repeat_with(inner, 1, new_age, consensus, gc_bin, model, rng, rm, ages, opts)
        {
            return Some(s);
        }
    }
    None
}

/// Evolve one simple/tandem repeat and return its mutated sequence.
///
/// Mirrors Perl's `evolveSimple`.
pub fn evolve_simple(
    sim: &SimpleRepeat,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    opts: EvolveOpts,
) -> Option<(Vec<u8>, char)> {
    let perl_compat = opts.perl_compat;
    // Perl: `$seed x (int($exp) + 1)` with $exp the TRF copy number.
    let exp: usize = sim.copy_num as usize + 1;
    let unit = sim.consensus.as_bytes();
    let mut seq: Vec<u8> = unit.iter().copied().cycle().take(unit.len() * exp).collect();
    if seq.len() < 10 {
        return None;
    }

    // Random orientation.
    let orient_minus = rng.gen::<f64>() > 0.5;
    let strand = if orient_minus { '-' } else { '+' };
    if orient_minus {
        seq = revcomp(&seq);
    }

    let mut_count = (sim.div * seq.len() as f64 / 100.0) as usize;
    let n_sit = mut_count / 2;
    let n_ver = mut_count - n_sit;
    let n_id = (sim.indel * seq.len() as f64 / 100.0) as usize;
    // Perl: $ndel = int(rand $nid); the rest are insertions.
    let n_del = if n_id > 0 { rng.gen_range(0..n_id) } else { 0 };
    let n_ins = n_id - n_del;

    let n_sit = n_sit - add_cpg_transitions(&mut seq, n_sit, opts.effective_cpg_rate(), rng);
    mutate_guided(add_transitions, &mut seq, n_sit, gc_bin, model, rng, perl_compat);
    mutate_guided(add_transversions, &mut seq, n_ver, gc_bin, model, rng, perl_compat);
    if n_ins > 0 {
        add_insertions(&mut seq, n_ins, gc_bin, model, rng, perl_compat);
    }
    mutate_guided(add_deletions, &mut seq, n_del, gc_bin, model, rng, perl_compat);

    Some((seq, strand))
}

// ---------------------------------------------------------------------------
// Bombardment loops
// ---------------------------------------------------------------------------

/// Insert interspersed TE repeats into `seq` until `te_frac` percent is covered.
///
/// `gc_bins[window_idx]` gives the GC bin for each WIN_SIZE window.
/// When `soft_mask` is true, inserted sequences are lower-cased (Perl default).
/// Returns a record for each successful insertion (for BED output).
/// Mirrors Perl's `insertRepeat`.
pub fn insert_te_repeats(
    seq: &mut Vec<u8>,
    gc_bins: &[u32],
    repeat_model: &RepeatModel,
    consensus: &HashMap<String, Vec<u8>>,
    seq_model: &SeqModel,
    te_frac: Option<f64>,
    soft_mask: bool,
    rng: &mut impl Rng,
    opts: EvolveOpts,
) -> Vec<InsertRecord> {
    let target_frac = te_frac.unwrap_or_else(|| rng.gen_range(10.0f64..60.0));
    log::info!("  => inserting {:.1}% TE elements", target_frac);

    let mut records: Vec<InsertRecord> = Vec::new();
    let mut repeat_bases: usize = 0;
    let mut tot_try: usize = 0;
    let ages = AgeIndex::new(repeat_model);

    loop {
        let cur_frac = 100.0 * repeat_bases as f64 / seq.len() as f64;
        if cur_frac >= target_frac {
            break;
        }
        tot_try += 1;
        if tot_try >= INS_CYCLES {  // Perl: last if ++$tot_try >= $ins_cyc
            break;
        }

        // Pick a random position; need at least 100 unmasked bases.
        let pos = rng.gen_range(0..seq.len().saturating_sub(100));
        // Check for 100 upper-case bases at this position.
        if seq[pos..pos.saturating_add(100).min(seq.len())]
            .iter()
            .any(|b| b.is_ascii_lowercase())
        {
            continue;
        }

        let win_idx = pos / seq_model.win_size;
        let gc_bin = match gc_bins.get(win_idx) {
            Some(&b) => b as usize,
            None => continue,
        };

        // Pick a random TE from this GC bin's profile.
        let bin_reps = &repeat_model.te_profile[gc_bin];
        if bin_reps.is_empty() {
            continue;
        }
        let rep = &bin_reps[rng.gen_range(0..bin_reps.len())];

        // Evolve the repeat.
        let (evolved, strand) = match evolve_repeat(rep, consensus, gc_bin, seq_model, rng, repeat_model, &ages, opts) {
            Some(pair) if pair.0.len() >= 10 => pair,
            _ => continue,
        };

        let ins_len = evolved.len();
        let end = pos + ins_len;
        if end > seq.len() {
            continue;
        }

        // Reject if insertion window is already occupied.
        if seq[pos..end].iter().any(|b| b.is_ascii_lowercase()) {
            continue;
        }

        // Write sequence (soft-masked or uppercase).
        if soft_mask {
            let lc: Vec<u8> = evolved.iter().map(|b| b.to_ascii_lowercase()).collect();
            seq[pos..end].copy_from_slice(&lc);
        } else {
            let uc: Vec<u8> = evolved.iter().map(|b| b.to_ascii_uppercase()).collect();
            seq[pos..end].copy_from_slice(&uc);
        }

        records.push(InsertRecord {
            start: pos,
            end,
            name: rep.name.clone(),
            family: rep.family.clone(),
            strand,
        });
        repeat_bases += ins_len;
        tot_try = 0;
    }
    records
}

/// Insert simple/tandem repeats into `seq` until `sim_frac` percent is covered.
///
/// When `soft_mask` is true, inserted sequences are lower-cased (Perl default).
/// Returns a record for each successful insertion (for BED output).
/// Mirrors Perl's `insertLowComplex`.
pub fn insert_simple_repeats(
    seq: &mut Vec<u8>,
    gc_bins: &[u32],
    repeat_model: &RepeatModel,
    seq_model: &SeqModel,
    sim_frac: Option<f64>,
    soft_mask: bool,
    rng: &mut impl Rng,
    opts: EvolveOpts,
) -> Vec<InsertRecord> {
    let target_frac = sim_frac.unwrap_or_else(|| rng.gen_range(0.0f64..2.0));
    log::info!("  => inserting {:.1}% low complexity sequences", target_frac);

    let mut records: Vec<InsertRecord> = Vec::new();
    let mut repeat_bases: usize = 0;
    let mut tot_try: usize = 0;

    loop {
        let cur_frac = 100.0 * repeat_bases as f64 / seq.len() as f64;
        if cur_frac >= target_frac {
            break;
        }
        tot_try += 1;
        if tot_try >= INS_CYCLES {  // Perl: last if ++$tot_try >= $ins_cyc
            break;
        }

        let pos = rng.gen_range(0..seq.len().saturating_sub(100));
        if seq[pos..pos.saturating_add(100).min(seq.len())]
            .iter()
            .any(|b| b.is_ascii_lowercase())
        {
            continue;
        }

        let win_idx = pos / seq_model.win_size;
        let gc_bin = match gc_bins.get(win_idx) {
            Some(&b) => b as usize,
            None => continue,
        };

        let bin_sims = &repeat_model.simple_profile[gc_bin];
        if bin_sims.is_empty() {
            continue;
        }
        let sim = &bin_sims[rng.gen_range(0..bin_sims.len())];

        let (evolved, strand) = match evolve_simple(sim, gc_bin, seq_model, rng, opts) {
            Some(pair) if pair.0.len() >= 10 => pair,
            _ => continue,
        };

        let ins_len = evolved.len();
        let end = pos + ins_len;
        if end > seq.len() {
            continue;
        }

        if seq[pos..end].iter().any(|b| b.is_ascii_lowercase()) {
            continue;
        }

        if soft_mask {
            let lc: Vec<u8> = evolved.iter().map(|b| b.to_ascii_lowercase()).collect();
            seq[pos..end].copy_from_slice(&lc);
        } else {
            let uc: Vec<u8> = evolved.iter().map(|b| b.to_ascii_uppercase()).collect();
            seq[pos..end].copy_from_slice(&uc);
        }

        records.push(InsertRecord {
            start: pos,
            end,
            name: sim.consensus.clone(),
            family: "Simple_repeat".to_string(),
            strand,
        });
        repeat_bases += ins_len;
        tot_try = 0;
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq_model::MODEL_FORMAT_VERSION;
    use rand::rngs::StdRng;
    use std::collections::BTreeMap;

    /// A one-bin model with k = 2: `freqs` maps each 1-base prefix to
    /// [P(A|pre), P(C|pre), P(G|pre), P(T|pre)].  Unlisted prefixes are uniform.
    fn model_k2(freqs: &[(&str, [f64; 4])]) -> SeqModel {
        let mut table: BTreeMap<String, [f64; 4]> =
            ["A", "C", "G", "T"].iter().map(|p| (p.to_string(), [0.25; 4])).collect();
        for (p, f) in freqs {
            table.insert(p.to_string(), *f);
        }
        let json = serde_json::json!({
            "model_version": MODEL_FORMAT_VERSION,
            "a_count": 0, "c_count": 0, "g_count": 0, "t_count": 0, "n_count": 0,
            "source_sequence_size": 1000, "win_size": 100, "kmer_size": 2,
            "gc_ranges": [{"start": 0, "end": 1000, "count": 10}],
            "gc_trans": [[1.0]], "gc_trans_counts": [[9.0]], "gc_kmer_totals": [100.0],
            "gc_binned_kmer_flank_freqs": [table],
            "repeat_model": null,
        });
        SeqModel::from_json(&json.to_string()).unwrap()
    }

    /// Default options without the CpG effect, so tests see the plain operators.
    const PLAIN: EvolveOpts = EvolveOpts { perl_compat: false, cpg_rate: 1.0 };
    const COMPAT: EvolveOpts = EvolveOpts { perl_compat: true, cpg_rate: 1.0 };

    const ALWAYS_A: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
    const ALWAYS_G: [f64; 4] = [0.0, 0.0, 1.0, 0.0];
    const ALWAYS_T: [f64; 4] = [0.0, 0.0, 0.0, 1.0];

    fn all(f: [f64; 4]) -> Vec<(&'static str, [f64; 4])> {
        ["A", "C", "G", "T"].iter().map(|&p| (p, f)).collect()
    }

    fn rng(seed: u64) -> StdRng {
        StdRng::seed_from_u64(seed)
    }

    fn count(seq: &[u8], b: u8) -> usize {
        seq.iter().filter(|&&x| x == b).count()
    }

    #[test]
    fn test_revcomp_and_kmer_freq() {
        assert_eq!(revcomp(b"AACGTn"), b"NACGTT".to_vec());
        let m = model_k2(&[("C", [0.1, 0.2, 0.3, 0.4])]);
        assert_eq!(kmer_freq(&m, 0, b"CG"), 0.3);
        assert_eq!(kmer_freq(&m, 0, b"cg"), 0.3, "case-insensitive");
        assert_eq!(kmer_freq(&m, 0, b"CN"), 0.0, "non-ACGT suffix");
        assert_eq!(kmer_freq(&m, 0, b"CGT"), 0.0, "wrong length");
    }

    /// About len/2 + 2 distinct positions, all in range (never more than len).
    #[test]
    fn test_sel_positions() {
        let mut r = rng(1);
        let p = sel_positions(1000, &mut r);
        assert_eq!(p.len(), 502);
        let uniq: std::collections::HashSet<_> = p.iter().collect();
        assert_eq!(uniq.len(), p.len(), "positions are distinct");
        assert!(p.iter().all(|&x| x < 1000));
        assert_eq!(sel_positions(3, &mut r).len(), 3);
        assert!(sel_positions(0, &mut r).is_empty());
    }

    /// The model only ever emits A, so G->A transitions always raise the
    /// context k-mer frequencies: all n are applied, nothing else changes.
    #[test]
    fn test_transitions_applied_when_model_favours_them() {
        let m = model_k2(&all(ALWAYS_A));
        let mut seq = vec![b'G'; 200];
        let left = add_transitions(&mut seq, 20, 0, &m, &mut rng(2), 0, false);
        assert_eq!(left, 0);
        assert_eq!(count(&seq, b'A'), 20);
        assert_eq!(count(&seq, b'G'), 180);
        assert_eq!(seq[0], b'G', "the first k-1 bases are never mutated");
    }

    /// The base that is evaluated is the one that changes.  On GTGT... with
    /// P(A|.) = P(T|.) = 0.5, G->A raises the incoming k-mer and T->C lowers
    /// it, so only G's mutate and every T survives.
    #[test]
    fn test_transition_mutates_the_evaluated_base() {
        let m = model_k2(&all([0.5, 0.0, 0.0, 0.5]));
        let mut seq: Vec<u8> = b"GT".repeat(100);
        let dropped = mutate_guided(add_transitions, &mut seq, 30, 0, &m, &mut rng(9), false);
        assert_eq!(dropped, 0);
        assert_eq!(count(&seq, b'A'), 30);
        assert_eq!(count(&seq, b'T'), 100);
        assert_eq!(count(&seq, b'C'), 0);
    }

    /// One round visits only ~len/2 positions; further rounds (up to
    /// MUT_CYCLES) are needed to place more mutations than that.
    #[test]
    fn test_mutations_retried_over_rounds() {
        let m = model_k2(&all(ALWAYS_A));
        let mut seq = vec![b'G'; 200];
        let dropped = mutate_guided(add_transitions, &mut seq, 150, 0, &m, &mut rng(10), false);
        assert_eq!(dropped, 0);
        assert_eq!(count(&seq, b'A'), 150);
    }

    /// The model only emits G, so every G->A is rejected; after MUT_CYCLES
    /// rounds the mutations are dropped and the sequence is unchanged.
    #[test]
    fn test_rejected_mutations_are_dropped() {
        let m = model_k2(&all(ALWAYS_G));
        let mut seq = vec![b'G'; 200];
        let dropped = mutate_guided(add_transitions, &mut seq, 20, 0, &m, &mut rng(3), false);
        assert_eq!(dropped, 20);
        assert_eq!(seq, vec![b'G'; 200]);
    }

    #[test]
    fn test_transversions_applied_when_model_favours_them() {
        let m = model_k2(&all(ALWAYS_T));
        let mut seq = vec![b'A'; 200];
        let dropped = mutate_guided(add_transversions, &mut seq, 30, 0, &m, &mut rng(4), false);
        assert_eq!(dropped, 0);
        assert_eq!(count(&seq, b'T'), 30);
        assert_eq!(seq.len(), 200);
    }

    /// Inserted bases are sampled from the context before `pos`: after C the
    /// model emits A, after G it emits T.  By default the base goes right
    /// after that context, so A follows C and T follows G.  With perl_compat
    /// it lands one base later (after seq[pos]), so A follows G and T follows C.
    #[test]
    fn test_insertions_placed_in_context() {
        let m = model_k2(&[("C", ALWAYS_A), ("G", ALWAYS_T)]);
        let base: Vec<u8> = b"CG".repeat(100);
        for perl_compat in [false, true] {
            let mut seq = base.clone();
            add_insertions(&mut seq, 10, 0, &m, &mut rng(5), perl_compat);
            assert_eq!(seq.len(), 210, "all 10 insertions applied (perl_compat={perl_compat})");
            // Pairs (preceding base, inserted base) for every inserted A or T.
            let pairs: Vec<(u8, u8)> = seq.windows(2)
                .filter(|w| w[1] == b'A' || w[1] == b'T')
                .map(|w| (w[0], w[1]))
                .collect();
            assert_eq!(pairs.len(), 10);
            let in_context = pairs.iter()
                .filter(|&&p| p == (b'C', b'A') || p == (b'G', b'T'))
                .count();
            if perl_compat {
                assert!(in_context <= 1, "Perl places bases after seq[pos]: {pairs:?}");
            } else {
                assert!(in_context >= 9, "default places bases in context: {pairs:?}");
            }
        }
    }

    /// A deletion is accepted when removing a base makes the k-mer ending
    /// there more likely.  On ACAC... with P(A|A) = P(A|C) = 1, "AC" -> "AA"
    /// passes, so the C's are deleted (by default deletions apply one at a
    /// time, and an occasional A next to an earlier deletion is a neutral,
    /// accepted change).  GARLIC 1.4 (`perl_compat`) scores the same k-mers
    /// but removes the base k-1 earlier, an A.
    #[test]
    fn test_deletions() {
        let m = model_k2(&[("A", ALWAYS_A), ("C", ALWAYS_A)]);
        for perl_compat in [false, true] {
            let mut seq: Vec<u8> = b"AC".repeat(100);
            let dropped = mutate_guided(add_deletions, &mut seq, 15, 0, &m, &mut rng(6), perl_compat);
            assert_eq!(dropped, 0, "perl_compat={perl_compat}");
            assert_eq!(seq.len(), 185);
            let (a_del, c_del) = (100 - count(&seq, b'A'), 100 - count(&seq, b'C'));
            if perl_compat {
                assert_eq!((a_del, c_del), (15, 0));
            } else {
                assert!(c_del >= 13, "mostly the scored C's: {a_del} A / {c_del} C deleted");
            }
        }

        // A less likely k-mer is always rejected.
        let worse = model_k2(&[("A", [0.0, 1.0, 0.0, 0.0]), ("C", ALWAYS_A)]);
        for perl_compat in [false, true] {
            let mut seq: Vec<u8> = b"AC".repeat(100);
            assert_eq!(mutate_guided(add_deletions, &mut seq, 5, 0, &worse, &mut rng(6), perl_compat), 5);
            assert_eq!(seq.len(), 200);
        }
    }

    /// Neutral changes (no k-mer more or less likely) are accepted; a change
    /// that makes any touched k-mer less likely is not.
    #[test]
    fn test_improves() {
        assert!(improves(&[0.2, 0.5], &[0.3, 0.5]));
        assert!(improves(&[0.2, 0.5], &[0.2, 0.5]), "ties accepted");
        assert!(!improves(&[0.2, 0.5], &[0.9, 0.4]));
        let m = model_k2(&[("A", ALWAYS_A)]);
        let mut seq = vec![b'A'; 200];
        assert_eq!(mutate_guided(add_deletions, &mut seq, 5, 0, &m, &mut rng(11), false), 0);
        assert_eq!(seq.len(), 195);
    }

    fn simple(consensus: &str, copy_num: f64) -> SimpleRepeat {
        SimpleRepeat { consensus: consensus.into(), dir: '+', copy_num, div: 0.0, indel: 0.0 }
    }

    /// int(copy_num) + 1 copies of the unit, either strand; < 10 bp -> None.
    #[test]
    fn test_evolve_simple() {
        let m = model_k2(&[]);
        let mut strands = std::collections::HashSet::new();
        for seed in 0..20 {
            let (s, strand) = evolve_simple(&simple("CAG", 4.7), 0, &m, &mut rng(seed), PLAIN).unwrap();
            let fwd = b"CAG".repeat(5);
            let expect = if strand == '+' { fwd.clone() } else { revcomp(&fwd) };
            assert_eq!(s, expect);
            strands.insert(strand);
        }
        assert_eq!(strands.len(), 2, "both orientations occur");
        assert!(evolve_simple(&simple("AC", 1.0), 0, &m, &mut rng(0), PLAIN).is_none());
    }

    fn te(name: &str, length: u32, div: f64, nfrg: u32) -> TeRepeat {
        TeRepeat { name: name.into(), family: "SINE/Alu".into(), dir: '+', div, ins: 0.0, del: 0.0,
                   length, nfrg }
    }

    fn consensus(names: &[&str], len: usize) -> HashMap<String, Vec<u8>> {
        let mut r = rng(99);
        names.iter()
            .map(|n| (n.to_string(), (0..len).map(|_| b"ACGT"[r.gen_range(0..4)]).collect()))
            .collect()
    }

    fn oriented(s: &[u8], strand: char) -> Vec<u8> {
        if strand == '-' { revcomp(s) } else { s.to_vec() }
    }

    fn contains(hay: &[u8], needle: &[u8]) -> bool {
        hay.windows(needle.len()).any(|w| w == needle)
    }

    /// Without divergence the fragment is a piece of the consensus of the
    /// record's length; a record longer than its consensus gives the full
    /// consensus, or with perl_compat a 3' piece (GARLIC 1.4 bug).
    #[test]
    fn test_evolve_repeat_fragment() {
        let m = model_k2(&[]);
        let cons = consensus(&["X"], 50);
        let x = &cons["X"];
        let rm = RepeatModel::new(1);
        let ages = AgeIndex::new(&rm);

        for seed in 0..10 {
            let (s, strand) = evolve_repeat(&te("X", 20, 0.0, 1), &cons, 0, &m, &mut rng(seed), &rm, &ages, PLAIN).unwrap();
            assert_eq!(s.len(), 20);
            assert!(contains(x, &oriented(&s, strand)), "fragment comes from the consensus");

            let (s, strand) = evolve_repeat(&te("X", 80, 0.0, 1), &cons, 0, &m, &mut rng(seed), &rm, &ages, PLAIN).unwrap();
            assert_eq!(oriented(&s, strand), *x, "over-long record -> full consensus");

            if let Some((s, strand)) = evolve_repeat(&te("X", 80, 0.0, 1), &cons, 0, &m, &mut rng(seed), &rm, &ages, COMPAT) {
                assert!(x.ends_with(&oriented(&s, strand)), "perl_compat -> 3' piece");
            }
        }
        assert!(evolve_repeat(&te("Missing", 20, 0.0, 1), &cons, 0, &m, &mut rng(0), &rm, &ages, PLAIN).is_none());
    }

    /// A fragmented, old host receives younger nested inserts (lower-case)
    /// chosen through the insert table; its own bases stay upper-case and,
    /// with no indels, keep the host fragment length.
    #[test]
    fn test_nested_inserts() {
        let m = model_k2(&[]);
        let cons = consensus(&["Host", "X", "Y"], 300);
        let mut rm = RepeatModel::new(1);
        for name in ["X", "Y"] {
            for div in [5.0, 6.0, 7.0] {
                rm.te_profile[0].push(te(name, 40, div, 1));
            }
        }
        let host = TeRepeat { family: "LINE/L1".into(), ..te("Host", 200, 40.0, 3) };
        rm.insert_table[0].insert(
            "Host#LINE/L1".into(),
            [("X#SINE/Alu".to_string(), 0.5), ("Y#SINE/Alu".to_string(), 0.5)].into_iter().collect(),
        );
        let ages = AgeIndex::new(&rm);

        let mut nested = 0;
        for seed in 0..10 {
            let (s, _) = evolve_repeat(&host, &cons, 0, &m, &mut rng(seed), &rm, &ages, PLAIN).unwrap();
            let upper = s.iter().filter(|b| b.is_ascii_uppercase()).count();
            let lower = s.len() - upper;
            assert_eq!(upper, 200, "host bases keep the fragment length");
            assert!(lower % 40 == 0 && lower <= 3 * 40, "whole 40-bp inserts, at most nfrg: {lower}");
            nested += (lower > 0) as usize;
        }
        assert!(nested >= 8, "an old fragmented host usually gets inserts ({nested}/10)");

        // Unfragmented hosts are never nested.
        let single = TeRepeat { nfrg: 1, ..host.clone() };
        let (s, _) = evolve_repeat(&single, &cons, 0, &m, &mut rng(0), &rm, &ages, PLAIN).unwrap();
        assert!(s.iter().all(|b| b.is_ascii_uppercase()));
    }

    /// Bombardment overwrites background in place: length unchanged, target
    /// coverage reached, records non-overlapping and exactly the lower-case
    /// spans; without soft-masking nothing is lower-case.
    #[test]
    fn test_insert_te_repeats() {
        let m = model_k2(&[]);
        let cons = consensus(&["X", "Y"], 300);
        let mut rm = RepeatModel::new(1);
        rm.te_profile[0].push(te("X", 200, 10.0, 1));
        rm.te_profile[0].push(te("Y", 300, 20.0, 1));
        let background: Vec<u8> = b"ACGT".repeat(5000);
        let gc_bins = vec![0u32; 200];

        let mut seq = background.clone();
        let mut recs = insert_te_repeats(&mut seq, &gc_bins, &rm, &cons, &m, Some(30.0), true, &mut rng(7), PLAIN);
        assert_eq!(seq.len(), background.len());
        let lower = seq.iter().filter(|b| b.is_ascii_lowercase()).count();
        let covered: usize = recs.iter().map(|r| r.end - r.start).sum();
        assert_eq!(lower, covered, "lower-case bases are exactly the inserted spans");
        assert!(covered as f64 >= 0.30 * seq.len() as f64, "target coverage reached");
        recs.sort_by_key(|r| r.start);
        for w in recs.windows(2) {
            assert!(w[0].end <= w[1].start, "inserts do not overlap");
        }
        for r in &recs {
            assert!(seq[r.start..r.end].iter().all(|b| b.is_ascii_lowercase()));
            assert!(r.name == "X" || r.name == "Y");
        }

        let mut seq = background.clone();
        let recs = insert_te_repeats(&mut seq, &gc_bins, &rm, &cons, &m, Some(30.0), false, &mut rng(7), PLAIN);
        assert!(!recs.is_empty());
        assert!(seq.iter().all(|b| b.is_ascii_uppercase()), "no soft-masking");
        assert_ne!(seq, background);
    }

    #[test]
    fn test_insert_simple_repeats() {
        let m = model_k2(&[]);
        let mut rm = RepeatModel::new(1);
        rm.simple_profile[0].push(simple("CAG", 9.0));
        let mut seq: Vec<u8> = b"ACGT".repeat(2500);
        let recs = insert_simple_repeats(&mut seq, &vec![0u32; 100], &rm, &m, Some(5.0), true, &mut rng(8), PLAIN);
        let covered: usize = recs.iter().map(|r| r.end - r.start).sum();
        assert!(covered >= 500);
        assert!(recs.iter().all(|r| r.end - r.start == 30 && r.family == "Simple_repeat"));
        assert_eq!(seq.iter().filter(|b| b.is_ascii_lowercase()).count(), covered);
    }

    /// CpG transitions: none without CpGs or with rate 1; otherwise each hit
    /// turns one CG into TG or CA, at most one per CpG and never more than the
    /// transitions requested, and nothing else changes.
    #[test]
    fn test_cpg_transitions() {
        let mut no_cpg = b"ACAGTT".repeat(50);
        assert_eq!(add_cpg_transitions(&mut no_cpg, 20, 10.0, &mut rng(1)), 0);
        let mut seq = b"ACGTTACGGA".repeat(30); // 300 bp, 60 CpGs
        assert_eq!(add_cpg_transitions(&mut seq.clone(), 20, 1.0, &mut rng(1)), 0);

        let before = seq.clone();
        let k = add_cpg_transitions(&mut seq, 20, 10.0, &mut rng(2));
        // r = 20 / (300 + 9 * 120) = 0.0145; CpG share = 10 * r * 120 = 17.4
        assert!((17..=18).contains(&k), "k = {k}");
        let cg = |s: &[u8]| s.windows(2).filter(|w| w == b"CG").count();
        assert_eq!(cg(&before) - cg(&seq), k);
        for (i, (a, b)) in before.iter().zip(&seq).enumerate() {
            if a != b {
                let ok = (*a == b'C' && *b == b'T' && before[i + 1] == b'G')
                      || (*a == b'G' && *b == b'A' && before[i - 1] == b'C');
                assert!(ok, "only CG -> TG / CA at {i}");
            }
        }
        // A high rate hits every CpG but never exceeds the requested count.
        let mut seq2 = before.clone();
        assert_eq!(add_cpg_transitions(&mut seq2, 30, 1000.0, &mut rng(3)), 30);
        let mut seq3 = before.clone();
        assert_eq!(add_cpg_transitions(&mut seq3, 500, 1000.0, &mut rng(3)), 60);
        assert_eq!(cg(&seq3), 0);
    }

    /// --perl_compat disables the CpG rate.
    #[test]
    fn test_perl_compat_disables_cpg_rate() {
        let o = EvolveOpts { perl_compat: true, cpg_rate: 10.0 };
        assert_eq!(o.effective_cpg_rate(), 1.0);
        assert_eq!(EvolveOpts { perl_compat: false, ..o }.effective_cpg_rate(), 10.0);
    }
}
