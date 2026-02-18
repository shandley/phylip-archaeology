// supplement_bound.rs
//
// Demonstrates the power of branch-and-bound pruning for exact parsimony search.
//
// In phylogenetics, the number of possible unrooted tree topologies grows
// super-exponentially with the number of taxa:
//
//   n taxa  ->  (2n-5)!! = (2n-5) x (2n-3) x ... x 3 x 1  unrooted topologies
//
// For 5 taxa there are 15 topologies; for 9 taxa there are 135,135. Exhaustive
// enumeration quickly becomes infeasible. The branch-and-bound algorithm
// (Hendy & Penny 1982) avoids evaluating most topologies by computing a
// lower bound on the parsimony score of any partial tree: if the partial tree
// already has a score >= the best complete tree found so far, the entire branch
// of the search space can be pruned.
//
// This example shows that branch-and-bound typically examines only a tiny
// fraction of all possible topologies, yet is guaranteed to find the globally
// optimal tree(s) -- unlike heuristic search which may get stuck in local optima.
//
// Usage:
//   cargo run --example supplement_bound

use phylip_rs::parsimony::branch_and_bound::{branch_and_bound_fitch, BranchAndBoundResult};
use phylip_rs::parsimony::wagner::search_with;
use phylip_rs::parsimony::traits::FitchScorer;
use phylip_rs::tree::newick::write_newick;
use phylip_rs::tree::{Alignment, Base, Sequence};

// ---------------------------------------------------------------------------
// Manually designed DNA alignments with clear phylogenetic signal
// ---------------------------------------------------------------------------

/// Build an alignment from (name, sequence_string) pairs.
fn make_alignment(data: &[(&str, &str)]) -> Alignment {
    let seqs: Vec<Sequence> = data
        .iter()
        .map(|(name, seq)| {
            let bases: Vec<Base> = seq
                .chars()
                .filter(|c| !c.is_whitespace())
                .map(|c| Base::from_char(c).unwrap())
                .collect();
            Sequence::new(*name, bases)
        })
        .collect();
    Alignment::new(seqs).unwrap()
}

/// Compute the double factorial: n!! = n x (n-2) x (n-4) x ... x 1.
/// For n taxa, the number of unrooted binary tree topologies is (2n-5)!!.
fn double_factorial(n: u64) -> u64 {
    if n <= 1 {
        return 1;
    }
    let mut result: u64 = 1;
    let mut k = n;
    while k > 1 {
        result = result.saturating_mul(k);
        k -= 2;
    }
    result
}

/// 5 taxa: clear grouping (AB)(CDE) with some homoplasy.
/// Taxa A,B share many sites; C,D,E share others, with a few convergent sites.
fn alignment_5_taxa() -> Alignment {
    make_alignment(&[
        //       signal: AB | CDE          homoplasy
        ("A", "AACCAACC GGTT AACCGGTT ACAC"),
        ("B", "AACCAACC GGTT AACCGGTT CACA"),
        ("C", "CCAACCAA TTGG CCAATTGG ACAC"),
        ("D", "CCAACCAA TTGG CCAATTGG GAGA"),
        ("E", "CCAACCAA TTGG GGTTAACC TGTG"),
    ])
}

/// 7 taxa: (AB)(CD)(EFG) structure with informative sites.
fn alignment_7_taxa() -> Alignment {
    make_alignment(&[
        //       AB signal   CD signal   EFG signal  noise
        ("A", "AACC AACC GGTT GGTT CCCC AACCGG ACGT"),
        ("B", "AACC AACC GGTT GGTT CCCC AACCGG CGTA"),
        ("C", "CCAA CCAA AACC AACC GGGG CCAAGG ACGT"),
        ("D", "CCAA CCAA AACC AACC GGGG CCAAGG GTAC"),
        ("E", "GGTT GGTT CCAA CCAA AACC GGTTAA ACGT"),
        ("F", "GGTT GGTT CCAA CCAA AACC GGTTAA TGCA"),
        ("G", "GGTT GGTT CCAA CCAA TTGG GGTTCC ACGT"),
    ])
}

/// 8 taxa: ((AB)(CD))((EF)(GH)) with informative variation.
fn alignment_8_taxa() -> Alignment {
    make_alignment(&[
        //       AB|CD deep  AB    CD    EF|GH deep  EF    GH    noise
        ("A", "AACC AACC AACC GGTT CCAA GGTT GGTT CCCC ACGT ATAT"),
        ("B", "AACC AACC AACC GGTT CCAA GGTT GGTT CCCC CGTA TATA"),
        ("C", "AACC AACC CCAA AACC CCAA GGTT GGTT CCCC ACGT GAGA"),
        ("D", "AACC AACC CCAA AACC CCAA GGTT GGTT CCCC GTAC AGAG"),
        ("E", "CCAA GGTT GGTT CCAA AACC AACC AACC GGGG ACGT ACAC"),
        ("F", "CCAA GGTT GGTT CCAA AACC AACC AACC GGGG TGCA CACA"),
        ("G", "CCAA GGTT GGTT CCAA GGTT CCAA CCAA TTTT ACGT GTGT"),
        ("H", "CCAA GGTT GGTT CCAA GGTT CCAA CCAA TTTT CGTA TGTG"),
    ])
}

/// 9 taxa: ((AB)(CD))((EF)(GHI)) with phylogenetic signal.
fn alignment_9_taxa() -> Alignment {
    make_alignment(&[
        //       ABCD deep  AB     CD     EFGHI deep EF     GHI    noise
        ("A", "AACC AACC AACC GGTT CCAA GGTT GGTT CCCC ACGT ATAT AC"),
        ("B", "AACC AACC AACC GGTT CCAA GGTT GGTT CCCC CGTA TATA CA"),
        ("C", "AACC AACC CCAA AACC CCAA GGTT GGTT CCCC ACGT GAGA AC"),
        ("D", "AACC AACC CCAA AACC CCAA GGTT GGTT CCCC GTAC AGAG GT"),
        ("E", "CCAA GGTT GGTT CCAA AACC AACC AACC GGGG ACGT ACAC AC"),
        ("F", "CCAA GGTT GGTT CCAA AACC AACC AACC GGGG TGCA CACA GT"),
        ("G", "CCAA GGTT GGTT CCAA GGTT CCAA CCAA TTTT ACGT GTGT AC"),
        ("H", "CCAA GGTT GGTT CCAA GGTT CCAA CCAA TTTT CGTA TGTG CA"),
        ("I", "CCAA GGTT GGTT CCAA GGTT TTGG TTGG AAAA ACGT GAGA GT"),
    ])
}

// ---------------------------------------------------------------------------
// Main demonstration
// ---------------------------------------------------------------------------

fn main() {
    println!();
    println!("BRANCH-AND-BOUND: Exact Search Through Exponential Space");
    println!("=========================================================");
    println!();
    println!("The number of possible unrooted tree topologies for n taxa is:");
    println!();
    println!("  T(n) = (2n-5)!! = (2n-5) x (2n-3) x ... x 3 x 1");
    println!();
    println!("This grows super-exponentially:");
    for n in [4u64, 5, 6, 7, 8, 9, 10, 12, 15, 20] {
        if 2 * n < 5 {
            continue;
        }
        let t = double_factorial(2 * n - 5);
        if t < u64::MAX / 2 {
            println!("  {:>2} taxa -> {:>15} topologies", n, format_number(t));
        } else {
            println!("  {:>2} taxa -> overflow (astronomically large)", n);
        }
    }
    println!();
    println!("Exhaustive enumeration is infeasible beyond ~10 taxa. The");
    println!("branch-and-bound algorithm (Hendy & Penny 1982) avoids this by");
    println!("pruning branches of the search tree whose lower bound on the");
    println!("parsimony score already exceeds the best complete tree found so far.");
    println!();
    println!("=========================================================");
    println!("EXPERIMENT: Branch-and-bound vs. total search space");
    println!("=========================================================");
    println!();

    // Define the test cases
    let cases: Vec<(&str, Alignment)> = vec![
        ("5 taxa", alignment_5_taxa()),
        ("7 taxa", alignment_7_taxa()),
        ("8 taxa", alignment_8_taxa()),
        ("9 taxa", alignment_9_taxa()),
    ];

    // Table header
    println!(
        "{:<10} {:>8} {:>12} {:>10} {:>10} {:>12} {:>10}",
        "Dataset", "Sites", "Total trees", "Examined", "Pruned", "Pruned %", "Optimal"
    );
    println!("{}", "-".repeat(78));

    let mut detailed_results: Vec<(String, usize, BranchAndBoundResult, u64)> = Vec::new();

    for (label, alignment) in &cases {
        let ntaxa = alignment.ntaxa() as u64;
        let nsites = alignment.nsites();
        let total_topologies = double_factorial(2 * ntaxa - 5);

        // Run branch-and-bound
        let bb_result = branch_and_bound_fitch(alignment, None);

        let examined = bb_result.trees_examined;
        let pruned = bb_result.trees_pruned;
        let total_visited = examined + pruned;
        let pruned_pct = if total_visited > 0 {
            100.0 * pruned as f64 / total_visited as f64
        } else {
            0.0
        };

        println!(
            "{:<10} {:>8} {:>12} {:>10} {:>10} {:>11.1}% {:>10}",
            label,
            nsites,
            format_number(total_topologies),
            format_number(examined as u64),
            format_number(pruned as u64),
            pruned_pct,
            bb_result.score
        );

        detailed_results.push((label.to_string(), nsites, bb_result, total_topologies));
    }

    println!();
    println!("=========================================================");
    println!("DETAILED RESULTS");
    println!("=========================================================");

    for (label, _nsites, bb_result, total_topologies) in &detailed_results {
        let ntrees_found = bb_result.trees.len();
        let examined = bb_result.trees_examined;

        println!();
        println!("--- {} ---", label);
        println!("  Total possible topologies:  {}", format_number(*total_topologies));
        println!("  Topologies examined by B&B: {}", format_number(examined as u64));
        println!("  Topologies pruned by B&B:   {}", format_number(bb_result.trees_pruned as u64));
        println!("  Optimal parsimony score:    {}", bb_result.score);
        println!("  Number of optimal trees:    {}", ntrees_found);

        if examined as u64 <= *total_topologies {
            let fraction = examined as f64 / *total_topologies as f64;
            if fraction < 0.01 {
                println!(
                    "  Search efficiency:          examined only {:.2}% of all topologies",
                    fraction * 100.0
                );
            } else {
                println!(
                    "  Search efficiency:          examined {:.1}% of all topologies",
                    fraction * 100.0
                );
            }
        }

        // Show the first optimal tree
        if let Some(tree) = bb_result.trees.first() {
            println!("  Best tree (Newick):         {}", write_newick(tree));
        }
    }

    println!();
    println!("=========================================================");
    println!("HEURISTIC vs. EXACT: Does heuristic find the optimum?");
    println!("=========================================================");
    println!();
    println!("Heuristic search (stepwise addition + SPR) is fast but not");
    println!("guaranteed to find the global optimum. Branch-and-bound is");
    println!("exact. Let's compare them:");
    println!();

    let cases2: Vec<(&str, Alignment)> = vec![
        ("5 taxa", alignment_5_taxa()),
        ("7 taxa", alignment_7_taxa()),
        ("8 taxa", alignment_8_taxa()),
        ("9 taxa", alignment_9_taxa()),
    ];

    println!(
        "{:<10} {:>14} {:>14} {:>10}",
        "Dataset", "Heuristic", "B&B (exact)", "Match?"
    );
    println!("{}", "-".repeat(52));

    for (label, alignment) in &cases2 {
        // Heuristic search with a fixed seed
        let heuristic = search_with(alignment, &FitchScorer, Some(42));
        let bb = branch_and_bound_fitch(alignment, None);

        let matches = if heuristic.score == bb.score {
            "YES"
        } else {
            "NO"
        };

        println!(
            "{:<10} {:>14} {:>14} {:>10}",
            label, heuristic.score, bb.score, matches
        );
    }

    println!();
    println!("For small datasets with clear signal, the heuristic typically");
    println!("finds the optimum. But with more taxa or weaker signal, the");
    println!("heuristic can get trapped in local optima. Branch-and-bound");
    println!("guarantees the global optimum -- at exponential cost.");

    println!();
    println!("=========================================================");
    println!("THE COMBINATORIAL EXPLOSION");
    println!("=========================================================");
    println!();
    println!("To appreciate why branch-and-bound pruning matters, consider");
    println!("the raw numbers:");
    println!();

    for (label, _nsites, bb_result, total_topologies) in &detailed_results {
        let examined = bb_result.trees_examined as u64;
        if examined > 0 && *total_topologies > examined {
            let speedup = *total_topologies as f64 / examined as f64;
            println!(
                "  {}: B&B examined {} of {} topologies ({:.0}x speedup)",
                label,
                format_number(examined),
                format_number(*total_topologies),
                speedup
            );
        }
    }

    println!();
    println!("The pruning is especially dramatic for larger datasets. With");
    println!("9 taxa, there are 135,135 possible unrooted topologies, but");
    println!("branch-and-bound can prune away the vast majority of them.");
    println!("This is what makes exact search practical for moderate numbers");
    println!("of taxa (up to ~12), while heuristic methods are needed for");
    println!("larger datasets.");

    println!();
    println!("=========================================================");
    println!("References");
    println!("=========================================================");
    println!();
    println!("  Hendy, M.D. & Penny, D. (1982). Branch and bound algorithms to");
    println!("  determine minimal evolutionary trees. Mathematical Biosciences,");
    println!("  59, 277-290.");
    println!();
    println!("  Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.");
    println!("  Chapter 7: Counting evolutionary changes on trees.");
    println!();
    println!("  Fitch, W.M. (1971). Toward defining the course of evolution:");
    println!("  minimum change for a specific tree topology. Systematic Zoology,");
    println!("  20, 406-416.");
    println!();
}

/// Format a number with comma separators for readability.
fn format_number(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}
