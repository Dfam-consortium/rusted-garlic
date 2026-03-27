pub mod consensus;
pub mod evolve;
pub mod io;
pub mod masking;
pub mod repeat_model;
pub mod seq_model;

use clap::{Args, Parser, Subcommand, ValueEnum};
use log::{error, info, warn};
use masking::MaskSet;
use seq_model::{EmitConfig, ModelConfig, SeqModel};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// Top-level CLI
// ---------------------------------------------------------------------------

#[derive(ValueEnum, Debug, Clone, Default)]
enum LogLevel {
    Error,
    Warn,
    #[default]
    Info,
    Debug,
    Trace,
    Off,
}

impl LogLevel {
    fn as_str(&self) -> &'static str {
        match self {
            LogLevel::Error => "error",
            LogLevel::Warn  => "warn",
            LogLevel::Info  => "info",
            LogLevel::Debug => "debug",
            LogLevel::Trace => "trace",
            LogLevel::Off   => "off",
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "rusted_garlic",
    version,
    author,
    about = "Rusted_garlic is a Rust port of the GARLIC artificial sequence generator\n\
             developed by Juan Caballero and G. Glusman (https://doi.org/10.1093/nar/gku356)."
)]
struct Cli {
    /// Log verbosity level (overridden by RUST_LOG if set)
    #[arg(long = "log-level", global = true, default_value = "info", value_name = "LEVEL")]
    log_level: LogLevel,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Build a sequence model from a genome FASTA
    Build(BuildArgs),
    /// Emit artificial sequence(s) from a saved model
    Emit(EmitArgs),
}

// ---------------------------------------------------------------------------
// Build subcommand
// ---------------------------------------------------------------------------

#[derive(Args, Debug)]
#[command(after_help = "\
Examples:
  # Build a basic model:
  rusted_garlic build --in_seq hg38.fa --model hg38.model

  # Build a masked model (repeats + simple repeats + genes):
  rusted_garlic build --in_seq hg38.fa --model hg38.model \\
      --repeats hg38.fa.align --trf hg38.trf.bed --genes refGene.txt

  # Build with non-default parameters:
  rusted_garlic build --in_seq hg38.fa --model hg38.model \\
      --kmer_size 5 --win_size 500 --num_gc_bins 20
")]
struct BuildArgs {
    /// Input genome sequence in FASTA format
    #[arg(short, long = "in_seq", value_name = "FASTA_FILE")]
    in_seq: PathBuf,

    /// Model file to save
    #[arg(short, long, value_name = "MODEL_FILE")]
    model: PathBuf,

    /// RepeatMasker .align or .out file for masking interspersed repeats
    #[arg(short, long, value_name = "ALIGN_FILE")]
    repeats: Option<PathBuf>,

    /// TRF BED file (UCSC format) for masking simple/tandem repeats
    #[arg(short, long, value_name = "BED_FILE")]
    trf: Option<PathBuf>,

    /// UCSC gene table (ensGene/refGene format) for masking gene bodies
    #[arg(short, long, value_name = "GENE_FILE")]
    genes: Option<PathBuf>,

    /// Window size in bp for GC isochore binning
    #[arg(short = 'w', long = "win_size", value_name = "INT", default_value = "1000")]
    win_size: usize,

    /// K-mer size for base frequency modelling
    #[arg(short = 'k', long = "kmer_size", value_name = "INT", default_value = "4")]
    kmer_size: usize,

    /// Maximum number of GC isochore bins
    #[arg(short = 'n', long = "num_gc_bins", value_name = "INT", default_value = "25")]
    num_gc_bins: usize,
}

// ---------------------------------------------------------------------------
// Emit subcommand
// ---------------------------------------------------------------------------

#[derive(ValueEnum, Debug, Clone, PartialEq)]
enum InsertMode {
    /// Insert features drawn from the feature FASTA until --te_frac coverage is reached
    Coverage,
    /// Insert every feature in the FASTA exactly once
    Full,
    /// Insert features with exactly --spacing bp of background between each one;
    /// --seq_size and --allow_nesting are not permitted in this mode
    #[value(name = "evenly-spaced")]
    EvenlySpaced,
}

#[derive(Args, Debug)]
#[command(after_help = "\
Examples:
  # Emit a 10 Mb sequence:
  rusted_garlic emit --model hg38.model --out_seq artificial.fa --seq_size 10m

  # Emit 20 independent 5 Mb sequences:
  rusted_garlic emit --model hg38.model --out_seq artificial.fa \\
      --seq_size 5m --num_seqs 20

  # Emit with internal TE bombardment (legacy):
  rusted_garlic emit --model hg38.model --out_seq artificial.fa \\
      --te_cons RepeatMaskerLib.fa --te_frac 30 --sim_frac 1

  # Insert external features, coverage-driven:
  rusted_garlic emit --model hg38.model --out_seq artificial.fa --seq_size 10m \\
      --insert_features tes.fa --insert_mode coverage --te_frac 35 --allow_nesting

  # Insert external features, full set:
  rusted_garlic emit --model hg38.model --out_seq artificial.fa --seq_size 10m \\
      --insert_features tes.fa --insert_mode full

  # Insert external features, evenly spaced (seq_size determined by feature set):
  rusted_garlic emit --model hg38.model --out_seq artificial.fa \\
      --insert_features tes.fa --insert_mode evenly-spaced --spacing 5000
")]
struct EmitArgs {
    /// Model file to load
    #[arg(short, long, value_name = "MODEL_FILE")]
    model: PathBuf,

    /// Output artificial sequence FASTA
    #[arg(short, long = "out_seq", value_name = "FASTA_FILE")]
    out_seq: PathBuf,

    /// Size of each sequence to generate (e.g. 1000, 500k, 2m, 1g)
    ///   [default: same size as the modeled sequence; not valid with evenly-spaced]
    #[arg(short, long = "seq_size", value_name = "SIZE")]
    seq_size: Option<String>,

    /// Number of artificial sequences to generate
    #[arg(short = 'N', long = "num_seqs", value_name = "INT", default_value = "1")]
    num_seqs: usize,

    // ---- Legacy internal bombardment ----------------------------------------

    /// TE consensus sequences in FASTA format (enables internal repeat bombardment)
    #[arg(long = "te_cons", value_name = "FASTA_FILE")]
    te_cons: Option<PathBuf>,

    /// Target simple-repeat fraction 0-100 [default: random 0-2 per sequence]
    #[arg(long = "sim_frac", value_name = "PERCENT")]
    sim_frac: Option<f64>,

    // ---- External feature insertion -----------------------------------------

    /// Feature sequences in FASTA format to insert into the background sequence
    #[arg(long = "insert_features", value_name = "FASTA_FILE")]
    insert_features: Option<PathBuf>,

    /// Insertion mode for --insert_features [default: coverage]
    ///   coverage:      draw features until --te_frac coverage is reached
    ///   full:          insert every feature in the FASTA exactly once
    ///   evenly-spaced: place features with --spacing bp of background between each
    #[arg(long = "insert_mode", value_name = "MODE", default_value = "coverage")]
    insert_mode: InsertMode,

    /// Target feature coverage 0-100; used by coverage mode and legacy --te_cons
    ///   [default: random 10-60 per sequence]
    #[arg(long = "te_frac", value_name = "PERCENT")]
    te_frac: Option<f64>,

    /// Background bp to place between each feature (evenly-spaced mode only)
    #[arg(long, value_name = "BP")]
    spacing: Option<usize>,

    /// Allow features to be inserted inside previously inserted regions
    ///   (coverage and full modes only)
    #[arg(long = "allow_nesting")]
    allow_nesting: bool,

    // ---- Output options ------------------------------------------------------

    /// Keep inserted repeats uppercase (default: lowercase / soft-masked)
    #[arg(long = "no_soft_mask")]
    no_soft_mask: bool,

    /// Write a BED file of inserted repeat locations alongside the sequence
    #[arg(long = "out_bed", value_name = "BED_FILE")]
    out_bed: Option<PathBuf>,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Parse a size string with optional suffix k/m/g (case-insensitive) into bytes.
fn parse_size(s: &str) -> Result<usize, String> {
    let s = s.trim();
    let (digits, suffix) = match s.chars().last() {
        Some(c) if c.is_alphabetic() => s.split_at(s.len() - c.len_utf8()),
        _ => (s, ""),
    };
    let base: usize = digits
        .parse()
        .map_err(|_| format!("invalid size value '{}'", digits))?;
    match suffix.to_lowercase().as_str() {
        "" => Ok(base),
        "k" => Ok(base * 1_000),
        "m" => Ok(base * 1_000_000),
        "g" => Ok(base * 1_000_000_000),
        other => Err(format!("unknown size suffix '{}' (use k, m, or g)", other)),
    }
}

fn die(msg: &str) -> ! {
    error!("{}", msg);
    std::process::exit(1);
}

// ---------------------------------------------------------------------------
// Subcommand handlers
// ---------------------------------------------------------------------------

fn run_build(args: BuildArgs) {
    let mask = build_mask(
        args.repeats.as_deref(),
        args.trf.as_deref(),
        args.genes.as_deref(),
    );
    let config = ModelConfig {
        mask: mask.as_ref(),
        repeat_file: args.repeats.as_deref(),
        trf_file: args.trf.as_deref(),
        win_size: args.win_size,
        kmer_size: args.kmer_size,
        num_gc_bins: args.num_gc_bins,
    };

    info!("Building model from {:?}...", args.in_seq);
    let m = SeqModel::model(&args.in_seq, Some(&config));

    info!("Saving model to {:?}...", args.model);
    m.to_file(&args.model);

    if m.repeat_model.is_some() {
        let gc_strs = m.gc_range_strings();
        let rm = m.repeat_model.as_ref().unwrap();
        let stem = args.model.with_extension("");
        let rep_path = stem.with_extension("repeats.W1000.data");
        let ins_path = stem.with_extension("inserts.W1000.data");
        info!("Writing repeat profile to {:?}...", rep_path);
        repeat_model::write_repeats_file(&rep_path, rm, &gc_strs);
        info!("Writing insert table to {:?}...", ins_path);
        repeat_model::write_inserts_file(&ins_path, rm, &gc_strs);
    }
}

fn run_emit(args: EmitArgs) {
    // --- Validate mode/flag combinations ------------------------------------
    if args.insert_features.is_some() {
        match args.insert_mode {
            InsertMode::EvenlySpaced => {
                if args.seq_size.is_some() {
                    die("--seq_size is not permitted with --insert_mode evenly-spaced \
                         (output length is determined by the feature set and --spacing)");
                }
                if args.allow_nesting {
                    die("--allow_nesting is not permitted with --insert_mode evenly-spaced");
                }
                if args.spacing.is_none() {
                    die("--spacing is required with --insert_mode evenly-spaced");
                }
            }
            InsertMode::Coverage | InsertMode::Full => {
                if args.spacing.is_some() {
                    die("--spacing is only valid with --insert_mode evenly-spaced");
                }
                if args.insert_mode == InsertMode::Coverage && args.te_frac.is_none() {
                    // te_frac has a random default — only warn, don't die
                    warn!(
                        "--te_frac not set; coverage mode will use a random target (10–60%)"
                    );
                }
            }
        }
    } else {
        // No insert_features — flags that only make sense with it are errors
        if args.spacing.is_some() {
            die("--spacing requires --insert_features");
        }
        if args.allow_nesting {
            die("--allow_nesting requires --insert_features");
        }
    }

    // --- Load model ---------------------------------------------------------
    info!("Loading model from {:?}...", args.model);
    let model = SeqModel::from_file(&args.model).unwrap_or_else(|e| {
        die(&format!("could not load model: {}", e));
    });

    // --- Parse seq_size -----------------------------------------------------
    let seq_size: Option<usize> = match &args.seq_size {
        Some(s) => Some(parse_size(s).unwrap_or_else(|e| die(&e))),
        None => None,
    };

    // --- Load legacy TE consensus if provided --------------------------------
    let consensus_map = args.te_cons.as_deref().map(|path| {
        info!("Loading TE consensus sequences from {:?}...", path);
        let mut rng = rand::thread_rng();
        let map = consensus::load_fasta(path, &mut rng);
        info!("  {} consensus sequences loaded", map.len());
        map
    });

    let emit_cfg = EmitConfig {
        consensus: consensus_map.as_ref(),
        te_frac: args.te_frac,
        sim_frac: args.sim_frac,
        soft_mask: !args.no_soft_mask,
    };

    // --- Open output files --------------------------------------------------
    let fasta_file =
        File::create(&args.out_seq).expect("Unable to create output FASTA file");
    let mut fasta_writer = BufWriter::new(fasta_file);

    let mut bed_writer: Option<BufWriter<File>> = args.out_bed.as_deref().map(|p| {
        BufWriter::new(File::create(p).expect("Unable to create output BED file"))
    });

    // --- Emit sequences -----------------------------------------------------
    let num_seqs = args.num_seqs.max(1);
    for seq_idx in 0..num_seqs {
        let seq_name = if num_seqs == 1 {
            "artificial_sequence".to_string()
        } else {
            format!("artificial_sequence_{}", seq_idx + 1)
        };

        if num_seqs > 1 {
            info!("Generating sequence {}/{}...", seq_idx + 1, num_seqs);
        } else {
            info!("Generating sequence...");
        }

        let (seq, records) = model.emit(seq_size, Some(&emit_cfg));
        info!("Writing {} bp ({})...", seq.len(), seq_name);

        writeln!(fasta_writer, ">{}", seq_name).expect("Write error");
        for chunk in seq.as_bytes().chunks(60) {
            fasta_writer.write_all(chunk).expect("Write error");
            writeln!(fasta_writer).expect("Write error");
        }

        if let Some(ref mut bw) = bed_writer {
            for rec in &records {
                writeln!(
                    bw,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    seq_name, rec.start, rec.end, rec.name, rec.family, rec.strand
                )
                .expect("Write error");
            }
        }
    }
    info!("Done. Wrote {} sequence(s) to {:?}", num_seqs, args.out_seq);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

fn main() {
    let cli = Cli::parse();

    env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or(cli.log_level.as_str()),
    )
    .format(|buf, record| {
        use std::io::Write;
        match record.level() {
            log::Level::Info => writeln!(buf, "{}", record.args()),
            _ => writeln!(buf, "{}: {}", record.level().as_str().to_lowercase(), record.args()),
        }
    })
    .init();

    match cli.command {
        Commands::Build(args) => run_build(args),
        Commands::Emit(args) => run_emit(args),
    }
}

// ---------------------------------------------------------------------------
// build_mask helper
// ---------------------------------------------------------------------------

/// Load annotation files into a MaskSet.  Returns None when no files were provided.
fn build_mask(
    repeats: Option<&Path>,
    trf: Option<&Path>,
    genes: Option<&Path>,
) -> Option<MaskSet> {
    if repeats.is_none() && trf.is_none() && genes.is_none() {
        return None;
    }
    let mut mask = MaskSet::new();

    if let Some(path) = repeats {
        info!("Loading repeat annotations from {:?}...", path);
        let intervals = io::align::parse(path);
        info!(
            "  {} masked intervals loaded",
            intervals.values().map(|v| v.len()).sum::<usize>()
        );
        mask.add(intervals);
    }
    if let Some(path) = trf {
        info!("Loading TRF annotations from {:?}...", path);
        let intervals = io::bed::parse_trf(path);
        info!(
            "  {} masked intervals loaded",
            intervals.values().map(|v| v.len()).sum::<usize>()
        );
        mask.add(intervals);
    }
    if let Some(path) = genes {
        info!("Loading gene annotations from {:?}...", path);
        let intervals = io::bed::parse_gene_table(path);
        info!(
            "  {} masked intervals loaded",
            intervals.values().map(|v| v.len()).sum::<usize>()
        );
        mask.add(intervals);
    }

    info!(
        "Total mask: {} intervals across all sequences",
        mask.interval_count()
    );
    Some(mask)
}
