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

// Maximum mutation retry cycles before relaxing acceptance (Perl's $mut_cyc).
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
/// After `MUT_CYCLES` passes the check is skipped.
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

/// Apply up to `n` single-base insertions to `seq`, growing the vector.
///
/// A base is sampled from the kmer model conditioned on the 3-mer prefix.
/// The insertion is accepted when the downstream kmer frequency is below the
/// sampled dice roll (biased toward locally-unusual positions, matching Perl).
fn add_insertions(
    seq: &mut Vec<u8>,
    mut n: usize,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
) {
    let kmer_size = model.kmer_size;
    let dna: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut pos_idx = 0usize;
    let mut tried = 0usize;
    while n > 0 && tried < seq.len().max(1) * 2 {
        tried += 1;
        let len = seq.len();
        if len < kmer_size + 1 {
            break;
        }
        let pos = rng.gen_range(kmer_size..len);
        // seed = seq[pos-(kmer_size-1)..pos]
        let seed = &seq[pos - kmer_size + 1..pos];
        // post = seq[pos]
        let post = seq[pos];

        // Sample new base from kmer distribution conditioned on seed
        let dice: f64 = rng.gen();
        let mut chosen = dna[3];
        let mut p = 0.0f64;
        let mut found = false;
        for &b in &dna {
            let kmer: Vec<u8> = seed.iter().copied().chain(std::iter::once(b)).collect();
            let q = p + kmer_freq(model, gc_bin, &kmer);
            if dice >= p {
                chosen = b;
            }
            if dice >= p && dice <= q + 1e-15 {
                found = true;
                break;
            }
            p = q;
        }
        if !found {
            pos_idx += 1;
            continue;
        }

        // Downstream acceptance: seq[pos-(kmer_size-2)..pos] + chosen + post
        let eed_start = pos - kmer_size + 2;
        let mut downstream: Vec<u8> = seq[eed_start..pos].to_vec();
        downstream.push(chosen);
        downstream.push(post);
        let f_down = kmer_freq(model, gc_bin, &downstream);
        if dice < f_down + 1e-15 {
            pos_idx += 1;
            continue;
        }

        seq.insert(pos, chosen);
        n -= 1;
        pos_idx = pos_idx.saturating_add(1);
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

/// Evolve one interspersed TE repeat and return its mutated sequence.
///
/// Mirrors Perl's `evolveRepeat`.  Returns `None` when the consensus is not
/// found or the fragment is too short.
///
/// Nested insertions (from `insert_table`) are attempted when `rep.nfrg > 1`.
pub fn evolve_repeat(
    rep: &TeRepeat,
    consensus: &HashMap<String, Vec<u8>>,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
    rm: &RepeatModel,
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

    // Pick a random start position in the consensus.
    let max_start = cons.len().saturating_sub(frag_len);
    let ini = if max_start > 0 { rng.gen_range(0..max_start) } else { 0 };
    let mut seq: Vec<u8> = cons[ini..ini + frag_len].to_vec();

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

    let skip_sit = add_transitions(&mut seq, n_sit, gc_bin, model, rng, 0);
    if skip_sit > 0 {
        add_transitions(&mut seq, skip_sit, gc_bin, model, rng, MUT_CYCLES);
    }
    let skip_ver = add_transversions(&mut seq, n_ver, gc_bin, model, rng, 0);
    if skip_ver > 0 {
        add_transversions(&mut seq, skip_ver, gc_bin, model, rng, MUT_CYCLES);
    }
    if n_ins > 0 {
        add_insertions(&mut seq, n_ins, gc_bin, model, rng);
    }
    let skip_del = add_deletions(&mut seq, n_del, gc_bin, model, rng, 0);
    if skip_del > 0 {
        add_deletions(&mut seq, skip_del, gc_bin, model, rng, MUT_CYCLES);
    }

    // Nested insertions: if nfrg > 1, try to insert a younger element inside.
    let age = rep.div + rep.ins + rep.del + (rep.nfrg as f64 * 10.0);
    if rep.nfrg > 1 && age > 10.0 {
        let outer_key = format!("{}#{}", rep.name, rep.family);
        if let Some(inner_map) = rm.insert_table.get(gc_bin).and_then(|t| t.get(&outer_key)) {
            if !inner_map.is_empty() {
                // Sample one inner repeat type proportionally.
                let dice: f64 = rng.gen();
                let mut p = 0.0f64;
                let mut chosen_inner: Option<&str> = None;
                for (inner_key, &freq) in inner_map.iter() {
                    p += freq;
                    if dice <= p {
                        chosen_inner = Some(inner_key.as_str());
                        break;
                    }
                }
                if let Some(inner_key) = chosen_inner {
                    // Find a repeat instance of this inner type in the GC bin.
                    let inner_name = inner_key.split('#').next().unwrap_or(inner_key);
                    let candidates: Vec<&TeRepeat> = rm.te_profile[gc_bin]
                        .iter()
                        .filter(|r| r.name == inner_name)
                        .collect();
                    if !candidates.is_empty() {
                        let inner_rep = *candidates.choose(rng).unwrap();
                        // Only insert if younger than parent.
                        let inner_age =
                            inner_rep.div + inner_rep.ins + inner_rep.del + (inner_rep.nfrg as f64 * 10.0);
                        if inner_age < age - 20.0 {
                            if let Some((inner_seq, _)) =
                                evolve_repeat(inner_rep, consensus, gc_bin, model, rng, rm)
                            {
                                if inner_seq.len() >= 10 && seq.len() > 1 {
                                    let target = rng.gen_range(0..seq.len());
                                    // Nested inserts are always lowercase (they are older elements).
                                    let lc_inner: Vec<u8> =
                                        inner_seq.iter().map(|b| b.to_ascii_lowercase()).collect();
                                    seq.splice(target..target, lc_inner);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Some((seq, strand))
}

/// Evolve one simple/tandem repeat and return its mutated sequence.
///
/// Mirrors Perl's `evolveSimple`.
pub fn evolve_simple(
    sim: &SimpleRepeat,
    gc_bin: usize,
    model: &SeqModel,
    rng: &mut impl Rng,
) -> Option<(Vec<u8>, char)> {
    // Pick a random expansion count (simulate stored ranges).
    let exp: usize = rng.gen_range(3..=20);
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
    let n_del = rng.gen_range(0..=n_id.max(1));
    let n_ins = n_id.saturating_sub(n_del);

    if n_sit > 0 {
        add_transitions(&mut seq, n_sit, gc_bin, model, rng, 0);
    }
    if n_ver > 0 {
        add_transversions(&mut seq, n_ver, gc_bin, model, rng, 0);
    }
    if n_ins > 0 {
        add_insertions(&mut seq, n_ins, gc_bin, model, rng);
    }
    if n_del > 0 {
        add_deletions(&mut seq, n_del, gc_bin, model, rng, 0);
    }

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
) -> Vec<InsertRecord> {
    let target_frac = te_frac.unwrap_or_else(|| rng.gen_range(10.0f64..60.0));
    eprintln!("  => inserting {:.1}% TE elements", target_frac);

    let mut records: Vec<InsertRecord> = Vec::new();
    let mut repeat_bases: usize = 0;
    let mut tot_try: usize = 0;
    let ins_cycles = INS_CYCLES;

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
        let (evolved, strand) = match evolve_repeat(rep, consensus, gc_bin, seq_model, rng, repeat_model) {
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
) -> Vec<InsertRecord> {
    let target_frac = sim_frac.unwrap_or_else(|| rng.gen_range(0.0f64..2.0));
    eprintln!("  => inserting {:.1}% low complexity sequences", target_frac);

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

        let (evolved, strand) = match evolve_simple(sim, gc_bin, seq_model, rng) {
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
