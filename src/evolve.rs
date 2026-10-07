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

/// Apply up to `n` A↔G / T↔C transitions to `seq`.
///
/// Each candidate position is accepted only when the transition increases the
/// k-mer frequency in the background model (Perl's acceptance criterion).
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
            let f_pre_old = kmer_freq(model, gc_bin, pre_old);
            let f_pre_new = kmer_freq(model, gc_bin, &pre_new);
            if f_pre_old >= f_pre_new + 1e-15 {
                continue;
            }
            // post context: kmer starting at pos-(kmer_size-1), including seq[pos]
            let post_old = &seq[pos - kmer_size + 1..pos + 1];
            let mut post_new = post_old.to_vec();
            post_new[kmer_size - 2] = new; // mutation is at index kmer_size-2 in post slice
            let f_post_old = kmer_freq(model, gc_bin, post_old);
            let f_post_new = kmer_freq(model, gc_bin, &post_new);
            if f_post_old >= f_post_new + 1e-15 {
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
            if kmer_freq(model, gc_bin, pre_old) >= kmer_freq(model, gc_bin, &pre_new) + 1e-15 {
                continue;
            }
            let post_old = &seq[pos - kmer_size + 1..pos + 1];
            let mut post_new = post_old.to_vec();
            post_new[kmer_size - 2] = new;
            if kmer_freq(model, gc_bin, post_old) >= kmer_freq(model, gc_bin, &post_new) + 1e-15 {
                continue;
            }
        }
        seq[pos - 1] = new;
        n -= 1;
    }
    n
}

/// Signature shared by the guided point-mutation operators.
type MutationOp<R> = fn(&mut Vec<u8>, usize, usize, &SeqModel, &mut R, usize) -> usize;

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
) -> usize {
    let mut bin = gc_bin;
    for eval in 0..MUT_CYCLES {
        if n == 0 {
            break;
        }
        n = op(seq, n, bin, model, rng, eval);
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
/// A position is accepted when the resulting kmer has higher background
/// frequency than the current kmer (Perl's deletion acceptance criterion).
fn add_deletions(
    seq: &mut Vec<u8>,
    mut n: usize,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    eval: usize,
) -> usize {
    let kmer_size = model.kmer_size;
    let len = seq.len();
    // Collect positions to delete, then remove them in reverse order.
    let mut to_delete: Vec<usize> = Vec::new();
    for pos in sel_positions(len, rng) {
        if n == 0 {
            break;
        }
        if pos < kmer_size || pos + kmer_size >= len {
            continue;
        }
        if eval < MUT_CYCLES {
            // old kmer = seq[pos..pos+kmer_size]
            let old_kmer = &seq[pos..pos + kmer_size];
            // new kmer = seq[pos..pos+kmer_size-1] + seq[pos+kmer_size]
            let mut new_kmer = seq[pos..pos + kmer_size - 1].to_vec();
            new_kmer.push(seq[pos + kmer_size]);
            let f_old = kmer_freq(model, gc_bin, old_kmer);
            let f_new = kmer_freq(model, gc_bin, &new_kmer);
            if f_old >= f_new + 1e-15 {
                continue;
            }
        }
        to_delete.push(pos);
        n -= 1;
    }
    // Remove from highest index to lowest so indices stay valid.
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
    perl_compat: bool,
) -> Option<(Vec<u8>, char)> {
    evolve_repeat_with(rep, rep.nfrg, f64::INFINITY, consensus, gc_bin, model, rng, rm, ages,
                       perl_compat)
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
    perl_compat: bool,
) -> Option<(Vec<u8>, char)> {
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

    mutate_guided(add_transitions, &mut seq, n_sit, gc_bin, model, rng);
    mutate_guided(add_transversions, &mut seq, n_ver, gc_bin, model, rng);
    if n_ins > 0 {
        add_insertions(&mut seq, n_ins, gc_bin, model, rng, perl_compat);
    }
    mutate_guided(add_deletions, &mut seq, n_del, gc_bin, model, rng);

    // Nested insertions (Perl: up to `break` inserts, `ins_cyc + 1` tries).
    // Targets are drawn on the un-nested sequence and applied from the highest
    // position down so earlier targets stay valid.
    if nfrg > 1 && age > 10.0 && !seq.is_empty() {
        let outer_key = format!("{}#{}", rep.name, rep.family);
        let mut pending: Vec<(usize, Vec<u8>)> = Vec::new();
        let mut num = 1u32;
        for _ in 0..=INS_CYCLES {
            num += 1;
            match get_insert(&outer_key, age, consensus, gc_bin, model, rng, rm, ages, perl_compat) {
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
    perl_compat: bool,
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
            evolve_repeat_with(inner, 1, new_age, consensus, gc_bin, model, rng, rm, ages, perl_compat)
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
    perl_compat: bool,
) -> Option<(Vec<u8>, char)> {
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

    mutate_guided(add_transitions, &mut seq, n_sit, gc_bin, model, rng);
    mutate_guided(add_transversions, &mut seq, n_ver, gc_bin, model, rng);
    if n_ins > 0 {
        add_insertions(&mut seq, n_ins, gc_bin, model, rng, perl_compat);
    }
    mutate_guided(add_deletions, &mut seq, n_del, gc_bin, model, rng);

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
    perl_compat: bool,
) -> Vec<InsertRecord> {
    let target_frac = te_frac.unwrap_or_else(|| rng.gen_range(10.0f64..60.0));
    log::info!("  => inserting {:.1}% TE elements", target_frac);

    let mut records: Vec<InsertRecord> = Vec::new();
    let mut repeat_bases: usize = 0;
    let mut tot_try: usize = 0;
    let ins_cycles = INS_CYCLES;
    let ages = AgeIndex::new(repeat_model);

    loop {
        let cur_frac = 100.0 * repeat_bases as f64 / seq.len() as f64;
        if cur_frac >= target_frac {
            break;
        }
        tot_try += 1;
        if tot_try > ins_cycles {
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
        let (evolved, strand) = match evolve_repeat(rep, consensus, gc_bin, seq_model, rng, repeat_model, &ages, perl_compat) {
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
    perl_compat: bool,
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
        if tot_try > INS_CYCLES {
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

        let (evolved, strand) = match evolve_simple(sim, gc_bin, seq_model, rng, perl_compat) {
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
