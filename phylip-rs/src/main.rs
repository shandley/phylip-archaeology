//! Command-line interface for phylip-rs.
//!
//! Provides access to phylogenetic analysis methods from PHYLIP via subcommands:
//! ml, parsimony, distance, bootstrap, and evaluate.

use std::env;
use std::fs;
use std::process;

use phylip_rs::distance::ml_distances::ml_distance_matrix;
use phylip_rs::distance::{fitch_margoliash, neighbor_joining};
use phylip_rs::io::fasta::read_fasta;
use phylip_rs::io::phylip_output::{format_distance_report, format_ml_report, format_parsimony_report};
use phylip_rs::likelihood::models::{F84Model, Jc69Model, SubstitutionModel};
use phylip_rs::likelihood::pruning;
use phylip_rs::likelihood::search;
use phylip_rs::likelihood::ttratio::{estimate_base_frequencies, estimate_ttratio_counting};
use phylip_rs::tree::newick::{parse_newick, write_newick};

// ============================================================
// Usage strings
// ============================================================

const MAIN_USAGE: &str = "\
phylip-rs - Modern Rust reimplementations of PHYLIP algorithms

Usage: phylip-rs <command> [options] <input_file>

Commands:
  ml          Maximum likelihood tree search (DNAML-like)
  parsimony   Wagner parsimony tree search (DNAPARS-like)
  distance    Distance-based tree building (NEIGHBOR-like)
  bootstrap   ML bootstrap analysis
  evaluate    Evaluate likelihood of a tree

Common options:
  -o, --output FILE     Output file (default: stdout)
  -s, --seed N          Random seed (default: 42)
  -h, --help            Print help information

Run 'phylip-rs <command> --help' for command-specific options.";

const ML_USAGE: &str = "\
phylip-rs ml - Maximum likelihood tree search

Usage: phylip-rs ml [options] <input.fasta>

Options:
  -o, --output FILE     Output file (default: stdout)
  -s, --seed N          Random seed (default: 42)
  --model MODEL         Substitution model: jc69, f84 (default: f84)
  --ttratio N           Transition/transversion ratio (default: auto-estimate)
  --gamma-cats N        Number of gamma rate categories (default: 0 = no gamma)
  -h, --help            Print this help message";

const PARSIMONY_USAGE: &str = "\
phylip-rs parsimony - Wagner parsimony tree search

Usage: phylip-rs parsimony [options] <input.fasta>

Options:
  -o, --output FILE     Output file (default: stdout)
  -s, --seed N          Random seed (default: 42)
  -h, --help            Print this help message";

const DISTANCE_USAGE: &str = "\
phylip-rs distance - Distance-based tree building

Usage: phylip-rs distance [options] <input.fasta>

Options:
  -o, --output FILE     Output file (default: stdout)
  --method METHOD       Tree method: nj, fm (default: nj)
  --model MODEL         Distance model: jc69, f84 (default: jc69)
  -h, --help            Print this help message";

const BOOTSTRAP_USAGE: &str = "\
phylip-rs bootstrap - ML bootstrap analysis

Usage: phylip-rs bootstrap [options] <input.fasta>

Options:
  -o, --output FILE     Output file (default: stdout)
  -s, --seed N          Random seed (default: 42)
  --model MODEL         Substitution model: jc69, f84 (default: f84)
  --ttratio N           Transition/transversion ratio (default: auto-estimate)
  --replicates N        Number of bootstrap replicates (default: 100)
  -h, --help            Print this help message";

const EVALUATE_USAGE: &str = "\
phylip-rs evaluate - Evaluate likelihood of a tree

Usage: phylip-rs evaluate [options] --tree <tree.nwk> <input.fasta>

Options:
  -o, --output FILE     Output file (default: stdout)
  --tree FILE           Newick tree file to evaluate (required)
  --model MODEL         Substitution model: jc69, f84 (default: f84)
  --ttratio N           Transition/transversion ratio (default: auto-estimate)
  -h, --help            Print this help message";

// ============================================================
// Argument parsing
// ============================================================

#[derive(Debug)]
struct Args {
    command: String,
    input_file: Option<String>,
    output_file: Option<String>,
    seed: u64,
    model: String,
    ttratio: Option<f64>,
    gamma_cats: usize,
    replicates: usize,
    distance_method: String,
    tree_file: Option<String>,
}

impl Args {
    fn default() -> Self {
        Self {
            command: String::new(),
            input_file: None,
            output_file: None,
            seed: 42,
            model: "f84".to_string(),
            ttratio: None,
            gamma_cats: 0,
            replicates: 100,
            distance_method: "nj".to_string(),
            tree_file: None,
        }
    }
}

fn parse_args() -> Result<Args, String> {
    let argv: Vec<String> = env::args().collect();

    if argv.len() < 2 {
        return Err(MAIN_USAGE.to_string());
    }

    // Check for top-level --help
    if argv[1] == "--help" || argv[1] == "-h" {
        return Err(MAIN_USAGE.to_string());
    }

    let mut args = Args::default();
    args.command = argv[1].clone();

    // Validate command
    match args.command.as_str() {
        "ml" | "parsimony" | "distance" | "bootstrap" | "evaluate" => {}
        other => {
            return Err(format!("Unknown command: '{}'\n\n{}", other, MAIN_USAGE));
        }
    }

    // Parse remaining arguments
    let rest = &argv[2..];
    let mut i = 0;
    while i < rest.len() {
        match rest[i].as_str() {
            "-h" | "--help" => {
                let usage = match args.command.as_str() {
                    "ml" => ML_USAGE,
                    "parsimony" => PARSIMONY_USAGE,
                    "distance" => DISTANCE_USAGE,
                    "bootstrap" => BOOTSTRAP_USAGE,
                    "evaluate" => EVALUATE_USAGE,
                    _ => MAIN_USAGE,
                };
                return Err(usage.to_string());
            }
            "-o" | "--output" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--output requires a file argument".to_string());
                }
                args.output_file = Some(rest[i].clone());
            }
            "-s" | "--seed" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--seed requires a numeric argument".to_string());
                }
                args.seed = rest[i]
                    .parse()
                    .map_err(|_| format!("Invalid seed: '{}'", rest[i]))?;
            }
            "--model" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--model requires an argument (jc69 or f84)".to_string());
                }
                let m = rest[i].to_lowercase();
                match m.as_str() {
                    "jc69" | "f84" => args.model = m,
                    _ => return Err(format!("Unknown model: '{}'. Use jc69 or f84.", rest[i])),
                }
            }
            "--ttratio" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--ttratio requires a numeric argument".to_string());
                }
                let val: f64 = rest[i]
                    .parse()
                    .map_err(|_| format!("Invalid ttratio: '{}'", rest[i]))?;
                if val <= 0.0 {
                    return Err("--ttratio must be positive".to_string());
                }
                args.ttratio = Some(val);
            }
            "--gamma-cats" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--gamma-cats requires a numeric argument".to_string());
                }
                args.gamma_cats = rest[i]
                    .parse()
                    .map_err(|_| format!("Invalid gamma-cats: '{}'", rest[i]))?;
            }
            "--replicates" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--replicates requires a numeric argument".to_string());
                }
                args.replicates = rest[i]
                    .parse()
                    .map_err(|_| format!("Invalid replicates: '{}'", rest[i]))?;
                if args.replicates == 0 {
                    return Err("--replicates must be at least 1".to_string());
                }
            }
            "--method" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--method requires an argument (nj or fm)".to_string());
                }
                let m = rest[i].to_lowercase();
                match m.as_str() {
                    "nj" | "fm" => args.distance_method = m,
                    _ => return Err(format!("Unknown method: '{}'. Use nj or fm.", rest[i])),
                }
            }
            "--tree" => {
                i += 1;
                if i >= rest.len() {
                    return Err("--tree requires a file argument".to_string());
                }
                args.tree_file = Some(rest[i].clone());
            }
            arg if arg.starts_with('-') => {
                return Err(format!("Unknown option: '{}'", arg));
            }
            _ => {
                // Positional argument: input file
                if args.input_file.is_some() {
                    return Err(format!("Unexpected argument: '{}'", rest[i]));
                }
                args.input_file = Some(rest[i].clone());
            }
        }
        i += 1;
    }

    // Validate required arguments
    if args.input_file.is_none() {
        let usage = match args.command.as_str() {
            "ml" => ML_USAGE,
            "parsimony" => PARSIMONY_USAGE,
            "distance" => DISTANCE_USAGE,
            "bootstrap" => BOOTSTRAP_USAGE,
            "evaluate" => EVALUATE_USAGE,
            _ => MAIN_USAGE,
        };
        return Err(format!("Missing input file\n\n{}", usage));
    }

    if args.command == "evaluate" && args.tree_file.is_none() {
        return Err(format!(
            "The evaluate command requires --tree <file>\n\n{}",
            EVALUATE_USAGE
        ));
    }

    Ok(args)
}

// ============================================================
// Output helper
// ============================================================

fn write_output(content: &str, output_file: &Option<String>) -> Result<(), String> {
    match output_file {
        Some(path) => {
            fs::write(path, content).map_err(|e| format!("Error writing to '{}': {}", path, e))
        }
        None => {
            print!("{}", content);
            Ok(())
        }
    }
}

// ============================================================
// Build model from args
// ============================================================

fn build_model(
    args: &Args,
    alignment: &phylip_rs::tree::types::Alignment,
) -> Box<dyn SubstitutionModel> {
    match args.model.as_str() {
        "jc69" => Box::new(Jc69Model),
        "f84" | _ => {
            let freqs = estimate_base_frequencies(alignment);
            let ttratio = match args.ttratio {
                Some(t) => t,
                None => {
                    let counting_est = estimate_ttratio_counting(alignment);
                    eprintln!("Auto-estimated ttratio: {:.4}", counting_est);
                    counting_est
                }
            };
            Box::new(F84Model::new(freqs, ttratio))
        }
    }
}

// ============================================================
// Command implementations
// ============================================================

fn cmd_ml(args: &Args) -> Result<(), String> {
    let input = fs::read_to_string(args.input_file.as_ref().unwrap())
        .map_err(|e| format!("Error reading input file: {}", e))?;
    let alignment = read_fasta(&input).map_err(|e| format!("Error parsing FASTA: {}", e))?;

    eprintln!(
        "ML search: {} taxa, {} sites",
        alignment.ntaxa(),
        alignment.nsites()
    );

    let model = build_model(args, &alignment);
    eprintln!("Model: {}", model.name());

    let result = search::search(&alignment, model.as_ref(), Some(args.seed))
        .map_err(|e| format!("ML search error: {}", e))?;

    let newick = write_newick(&result.tree);
    let report = format_ml_report(&result, model.name());

    let mut output = String::new();
    output.push_str(&report);
    output.push_str("\nBest tree (Newick):\n");
    output.push_str(&newick);
    output.push('\n');

    write_output(&output, &args.output_file)
}

fn cmd_parsimony(args: &Args) -> Result<(), String> {
    let input = fs::read_to_string(args.input_file.as_ref().unwrap())
        .map_err(|e| format!("Error reading input file: {}", e))?;
    let alignment = read_fasta(&input).map_err(|e| format!("Error parsing FASTA: {}", e))?;

    eprintln!(
        "Parsimony search: {} taxa, {} sites",
        alignment.ntaxa(),
        alignment.nsites()
    );

    let result =
        phylip_rs::parsimony::wagner::search(&alignment, Some(args.seed));

    let newick = write_newick(&result.tree);
    let report = format_parsimony_report(&result);

    let mut output = String::new();
    output.push_str(&report);
    output.push_str("\nBest tree (Newick):\n");
    output.push_str(&newick);
    output.push('\n');

    write_output(&output, &args.output_file)
}

fn cmd_distance(args: &Args) -> Result<(), String> {
    let input = fs::read_to_string(args.input_file.as_ref().unwrap())
        .map_err(|e| format!("Error reading input file: {}", e))?;
    let alignment = read_fasta(&input).map_err(|e| format!("Error parsing FASTA: {}", e))?;

    eprintln!(
        "Distance analysis: {} taxa, {} sites",
        alignment.ntaxa(),
        alignment.nsites()
    );

    // Build the distance model (used for computing the distance matrix)
    let dist_model: Box<dyn SubstitutionModel> = match args.model.as_str() {
        "f84" => {
            let freqs = estimate_base_frequencies(&alignment);
            let ttratio = args.ttratio.unwrap_or(2.0);
            Box::new(F84Model::new(freqs, ttratio))
        }
        _ => Box::new(Jc69Model),
    };

    let matrix = ml_distance_matrix(&alignment, dist_model.as_ref());

    let method_name = match args.distance_method.as_str() {
        "fm" => "Fitch-Margoliash",
        _ => "Neighbor-Joining",
    };

    eprintln!("Method: {}", method_name);

    let tree = match args.distance_method.as_str() {
        "fm" => fitch_margoliash(&matrix),
        _ => neighbor_joining(&matrix),
    };

    let newick = write_newick(&tree);
    let report = format_distance_report(&tree, &matrix, method_name);

    let mut output = String::new();
    output.push_str(&report);
    output.push_str("\nTree (Newick):\n");
    output.push_str(&newick);
    output.push('\n');

    write_output(&output, &args.output_file)
}

fn cmd_bootstrap(args: &Args) -> Result<(), String> {
    let input = fs::read_to_string(args.input_file.as_ref().unwrap())
        .map_err(|e| format!("Error reading input file: {}", e))?;
    let alignment = read_fasta(&input).map_err(|e| format!("Error parsing FASTA: {}", e))?;

    eprintln!(
        "Bootstrap analysis: {} taxa, {} sites, {} replicates",
        alignment.ntaxa(),
        alignment.nsites(),
        args.replicates
    );

    let model = build_model(args, &alignment);
    eprintln!("Model: {}", model.name());

    let result = phylip_rs::bootstrap::ml::ml_bootstrap(
        &alignment,
        model.as_ref(),
        args.replicates,
        args.seed,
    )
    .map_err(|e| format!("Bootstrap error: {}", e))?;

    let newick = write_newick(&result.tree);

    let mut output = String::new();
    output.push_str(&format!(
        "Bootstrap analysis complete ({} replicates)\n",
        args.replicates
    ));
    output.push_str(&format!("Log-likelihood of best tree: {:.6}\n", result.lnl));
    output.push_str(&format!(
        "Number of taxa: {}\n",
        result.tree.num_leaves()
    ));
    output.push('\n');

    // Show support values
    let supports = result.support_values.all_supports();
    if !supports.is_empty() {
        output.push_str("Bootstrap support values:\n");
        for (bipartition, support) in &supports {
            output.push_str(&format!("  {}: {:.1}%\n", bipartition, support * 100.0));
        }
        output.push('\n');
    }

    output.push_str("Best tree with bootstrap support (Newick):\n");
    output.push_str(&newick);
    output.push('\n');

    write_output(&output, &args.output_file)
}

fn cmd_evaluate(args: &Args) -> Result<(), String> {
    let input = fs::read_to_string(args.input_file.as_ref().unwrap())
        .map_err(|e| format!("Error reading input file: {}", e))?;
    let alignment = read_fasta(&input).map_err(|e| format!("Error parsing FASTA: {}", e))?;

    let tree_str = fs::read_to_string(args.tree_file.as_ref().unwrap())
        .map_err(|e| format!("Error reading tree file: {}", e))?;
    let mut tree =
        parse_newick(&tree_str).map_err(|e| format!("Error parsing Newick tree: {}", e))?;

    eprintln!(
        "Evaluating tree: {} taxa, {} sites",
        alignment.ntaxa(),
        alignment.nsites()
    );

    let model = build_model(args, &alignment);
    eprintln!("Model: {}", model.name());

    // Optimize branch lengths
    let lnl = pruning::optimize_branch_lengths(&mut tree, &alignment, model.as_ref())
        .map_err(|e| format!("Likelihood error: {}", e))?;

    let site_lnl = pruning::site_log_likelihoods(&tree, &alignment, model.as_ref())
        .map_err(|e| format!("Likelihood error: {}", e))?;

    let newick = write_newick(&tree);

    let mut output = String::new();
    output.push_str("Tree evaluation\n");
    output.push_str(&format!("Model: {}\n", model.name()));
    output.push_str(&format!("Number of taxa: {}\n", tree.num_leaves()));
    output.push_str(&format!("Number of sites: {}\n", alignment.nsites()));
    output.push_str(&format!("Log-likelihood: {:.6}\n", lnl));

    // Per-site summary
    if !site_lnl.is_empty() {
        let mean = site_lnl.iter().sum::<f64>() / site_lnl.len() as f64;
        let min = site_lnl.iter().copied().fold(f64::INFINITY, f64::min);
        let max = site_lnl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        output.push_str(&format!("Per-site lnL mean: {:.6}\n", mean));
        output.push_str(&format!("Per-site lnL min:  {:.6}\n", min));
        output.push_str(&format!("Per-site lnL max:  {:.6}\n", max));
    }

    output.push_str(&format!("\nOptimized tree (Newick):\n{}\n", newick));

    write_output(&output, &args.output_file)
}

// ============================================================
// Main
// ============================================================

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(msg) => {
            // If the message looks like usage/help text, print to stdout and exit 0.
            // Otherwise it is an error; print to stderr and exit 1.
            if msg.starts_with("phylip-rs") {
                println!("{}", msg);
                process::exit(0);
            } else {
                eprintln!("{}", msg);
                process::exit(1);
            }
        }
    };

    let result = match args.command.as_str() {
        "ml" => cmd_ml(&args),
        "parsimony" => cmd_parsimony(&args),
        "distance" => cmd_distance(&args),
        "bootstrap" => cmd_bootstrap(&args),
        "evaluate" => cmd_evaluate(&args),
        _ => unreachable!(),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
