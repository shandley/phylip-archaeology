// genetic_code_distances.rs
//
// Computes and visualizes the 20x20 amino acid step matrix derived from the
// standard genetic code, then demonstrates that the genetic code minimizes
// the cost of point mutations compared to random codes.
//
// The genetic code step matrix gives the minimum number of nucleotide
// substitutions required to convert between any pair of amino acids. This
// matrix is central to PHYLIP's protpars program (Sankoff algorithm for
// protein parsimony), where changes between amino acids that are close
// in the code cost fewer steps than changes between distant ones.
//
// Francis Crick's "frozen accident" hypothesis (1968) proposed that the
// genetic code is largely arbitrary. However, comparison with randomized
// codes reveals that the real code has significantly lower average step
// distances -- evidence that natural selection has optimized the code to
// buffer against the deleterious effects of point mutations.
//
// Usage:
//   cargo run --example genetic_code_distances

use phylip_rs::models::protein::{AminoAcid, ProteinAlignment, ProteinSequence};
use phylip_rs::parsimony::protein_parsimony::{
    genetic_code_step_matrix, protein_parsimony_score, AA_ORDER,
};
use phylip_rs::tree::newick::parse_newick;

// ---------------------------------------------------------------------------
// Simple LCG random number generator (no external dependencies needed)
// ---------------------------------------------------------------------------

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn gen_range(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = self.gen_range(i + 1);
            slice.swap(i, j);
        }
    }
}

// ---------------------------------------------------------------------------
// Amino acid full names
// ---------------------------------------------------------------------------

fn aa_full_name(c: char) -> &'static str {
    match c {
        'A' => "Ala",
        'C' => "Cys",
        'D' => "Asp",
        'E' => "Glu",
        'F' => "Phe",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'K' => "Lys",
        'L' => "Leu",
        'M' => "Met",
        'N' => "Asn",
        'P' => "Pro",
        'Q' => "Gln",
        'R' => "Arg",
        'S' => "Ser",
        'T' => "Thr",
        'V' => "Val",
        'W' => "Trp",
        'Y' => "Tyr",
        _ => "???",
    }
}

// ---------------------------------------------------------------------------
// Polar requirement values for the 20 amino acids (Woese, 1973).
//
// Polar requirement measures an amino acid's chromatographic partitioning
// behavior and serves as a proxy for physicochemical similarity.
// Freeland & Hurst (1998) used this property to show that the genetic
// code is optimized: single-nucleotide mutations in the real code produce
// smaller changes in polar requirement than in random codes.
//
// Values indexed by AA_ORDER position:
// A=Ala, C=Cys, D=Asp, E=Glu, F=Phe, G=Gly, H=His, I=Ile, K=Lys, L=Leu,
// M=Met, N=Asn, P=Pro, Q=Gln, R=Arg, S=Ser, T=Thr, V=Val, W=Trp, Y=Tyr
// ---------------------------------------------------------------------------

const POLAR_REQUIREMENT: [f64; 20] = [
    7.0,  // A - Ala
    4.8,  // C - Cys
    13.0, // D - Asp
    12.5, // E - Glu
    5.0,  // F - Phe
    7.9,  // G - Gly
    8.4,  // H - His
    4.9,  // I - Ile
    10.1, // K - Lys
    4.9,  // L - Leu
    5.3,  // M - Met
    10.0, // N - Asn
    6.6,  // P - Pro
    8.6,  // Q - Gln
    9.1,  // R - Arg
    7.5,  // S - Ser
    6.6,  // T - Thr
    5.6,  // V - Val
    5.2,  // W - Trp
    5.4,  // Y - Tyr
];

// ---------------------------------------------------------------------------
// Build a codon -> amino acid lookup table for the standard genetic code,
// then measure the mean squared change in polar requirement over all
// possible single-nucleotide nonsynonymous mutations.
// ---------------------------------------------------------------------------

/// Build a lookup: codon (as index 0..63) -> amino acid index (0..19), or
/// 255 for stop codons. Codon index = pos1 * 16 + pos2 * 4 + pos3.
fn build_codon_to_aa() -> [u8; 64] {
    let mut table = [255u8; 64];

    // Standard genetic code (RNA: U=0, C=1, A=2, G=3)
    let assignments: &[(char, &[[u8; 3]])] = &[
        ('F', &[[0,0,0], [0,0,1]]),
        ('L', &[[0,0,2], [0,0,3], [1,0,0], [1,0,1], [1,0,2], [1,0,3]]),
        ('I', &[[2,0,0], [2,0,1], [2,0,2]]),
        ('M', &[[2,0,3]]),
        ('V', &[[3,0,0], [3,0,1], [3,0,2], [3,0,3]]),
        ('S', &[[0,1,0], [0,1,1], [0,1,2], [0,1,3], [2,3,0], [2,3,1]]),
        ('P', &[[1,1,0], [1,1,1], [1,1,2], [1,1,3]]),
        ('T', &[[2,1,0], [2,1,1], [2,1,2], [2,1,3]]),
        ('A', &[[3,1,0], [3,1,1], [3,1,2], [3,1,3]]),
        ('Y', &[[0,2,0], [0,2,1]]),
        ('H', &[[1,2,0], [1,2,1]]),
        ('Q', &[[1,2,2], [1,2,3]]),
        ('N', &[[2,2,0], [2,2,1]]),
        ('K', &[[2,2,2], [2,2,3]]),
        ('D', &[[3,2,0], [3,2,1]]),
        ('E', &[[3,2,2], [3,2,3]]),
        ('C', &[[0,3,0], [0,3,1]]),
        ('W', &[[0,3,3]]),
        ('R', &[[1,3,0], [1,3,1], [1,3,2], [1,3,3], [2,3,2], [2,3,3]]),
        ('G', &[[3,3,0], [3,3,1], [3,3,2], [3,3,3]]),
        // Stop codons: UAA=0,2,2 -> idx 10; UAG=0,2,3 -> idx 11; UGA=0,3,2 -> idx 14
        // These remain 255 (stop).
    ];

    for &(aa_char, codons) in assignments {
        let aa_idx = AA_ORDER.iter().position(|&c| c == aa_char).unwrap() as u8;
        for codon in codons {
            let idx = (codon[0] as usize) * 16 + (codon[1] as usize) * 4 + codon[2] as usize;
            table[idx] = aa_idx;
        }
    }

    table
}

/// Given a codon-to-AA lookup table, compute the mean squared change in
/// polar requirement over all possible single-nucleotide nonsynonymous
/// mutations.
///
/// For each of the 61 sense codons, for each of 3 positions, for each of
/// 3 alternative nucleotides, if the mutation is nonsynonymous (different
/// amino acid), add the squared difference in polar requirement.
fn mean_squared_polar_change(codon_to_aa: &[u8; 64]) -> f64 {
    let mut total_sq_change = 0.0f64;
    let mut n_nonsyn = 0usize;

    for codon_idx in 0..64usize {
        let aa1 = codon_to_aa[codon_idx];
        if aa1 == 255 {
            continue; // stop codon
        }

        let pos = [
            (codon_idx / 16) as u8,
            ((codon_idx / 4) % 4) as u8,
            (codon_idx % 4) as u8,
        ];

        for position in 0..3 {
            for alt_nuc in 0..4u8 {
                if alt_nuc == pos[position] {
                    continue;
                }
                let mut mutant = pos;
                mutant[position] = alt_nuc;
                let mut_idx =
                    (mutant[0] as usize) * 16 + (mutant[1] as usize) * 4 + mutant[2] as usize;
                let aa2 = codon_to_aa[mut_idx];
                if aa2 == 255 {
                    continue; // mutant is stop codon
                }
                if aa1 == aa2 {
                    continue; // synonymous
                }
                let pr1 = POLAR_REQUIREMENT[aa1 as usize];
                let pr2 = POLAR_REQUIREMENT[aa2 as usize];
                total_sq_change += (pr1 - pr2) * (pr1 - pr2);
                n_nonsyn += 1;
            }
        }
    }

    if n_nonsyn == 0 {
        return 0.0;
    }
    total_sq_change / n_nonsyn as f64
}

/// Create a permuted codon-to-AA table by shuffling which amino acid
/// identity is assigned to each codon block.
///
/// `perm[i]` = the new amino acid index for what was originally amino acid i.
fn permute_codon_table(original: &[u8; 64], perm: &[usize; 20]) -> [u8; 64] {
    let mut permuted = *original;
    for entry in &mut permuted {
        if *entry != 255 {
            *entry = perm[*entry as usize] as u8;
        }
    }
    permuted
}

// ---------------------------------------------------------------------------
// Main demonstration
// ---------------------------------------------------------------------------

fn main() {
    println!();
    println!("THE GENETIC CODE STEP MATRIX: Evolution's Error-Correcting Code");
    println!("================================================================");
    println!();
    println!("The genetic code maps 61 sense codons (triplets of nucleotides) to");
    println!("20 amino acids. A single nucleotide mutation in a codon may change");
    println!("the encoded amino acid, but not all amino acid changes are equally");
    println!("likely: some pairs can be reached by a single nucleotide change,");
    println!("while others require 2 or 3 changes.");
    println!();
    println!("The 'step matrix' captures this: entry (i,j) is the MINIMUM number");
    println!("of nucleotide substitutions needed to convert any codon for amino");
    println!("acid i into any codon for amino acid j.");
    println!();
    println!("PHYLIP's protpars program uses this matrix for protein parsimony:");
    println!("a change from Ala to Val costs 1 step (one nucleotide change),");
    println!("while a change from Trp to Met costs 2 steps.");
    println!();

    let matrix = genetic_code_step_matrix();

    // ===================================================================
    // Part 1: Display the 20x20 step matrix
    // ===================================================================

    println!("PART 1: The 20x20 Genetic Code Step Matrix");
    println!("==========================================");
    println!();
    println!("Amino acids in alphabetical one-letter code order:");
    println!("A=Ala C=Cys D=Asp E=Glu F=Phe G=Gly H=His I=Ile K=Lys L=Leu");
    println!("M=Met N=Asn P=Pro Q=Gln R=Arg S=Ser T=Thr V=Val W=Trp Y=Tyr");
    println!();

    // Header row
    print!("     ");
    for &c in &AA_ORDER {
        print!("  {} ", c);
    }
    println!();

    // Separator
    print!("    +");
    for _ in 0..20 {
        print!("----");
    }
    println!();

    // Matrix rows
    for (i, &row_aa) in AA_ORDER.iter().enumerate() {
        print!("  {} |", row_aa);
        for j in 0..20 {
            if i == j {
                print!("  . ");
            } else {
                print!("  {} ", matrix[i][j]);
            }
        }
        println!();
    }
    println!();

    // ===================================================================
    // Part 2: Distribution of distances
    // ===================================================================

    println!("PART 2: Distribution of Step Distances");
    println!("======================================");
    println!();

    let mut dist_count = [0usize; 4]; // index 0 unused, 1-3 used
    let mut total_distance = 0u64;
    let mut pair_count = 0u64;

    for i in 0..20 {
        for j in (i + 1)..20 {
            let d = matrix[i][j] as usize;
            dist_count[d] += 1;
            total_distance += d as u64;
            pair_count += 1;
        }
    }

    let avg_distance = total_distance as f64 / pair_count as f64;

    println!("  Total amino acid pairs: {} (20 choose 2)", pair_count);
    println!(
        "  1-step pairs (single nucleotide change): {:>3}  ({:.1}%)",
        dist_count[1],
        100.0 * dist_count[1] as f64 / pair_count as f64
    );
    println!(
        "  2-step pairs (two nucleotide changes):   {:>3}  ({:.1}%)",
        dist_count[2],
        100.0 * dist_count[2] as f64 / pair_count as f64
    );
    println!(
        "  3-step pairs (three nucleotide changes): {:>3}  ({:.1}%)",
        dist_count[3],
        100.0 * dist_count[3] as f64 / pair_count as f64
    );
    println!("  Average pairwise step distance: {:.3}", avg_distance);
    println!();

    // ===================================================================
    // Part 3: Most and least connected amino acids
    // ===================================================================

    println!("PART 3: Most and Least Connected Amino Acids");
    println!("=============================================");
    println!();

    let mut avg_steps: Vec<(char, f64)> = Vec::new();
    for i in 0..20 {
        let total: u64 = (0..20)
            .filter(|&j| j != i)
            .map(|j| matrix[i][j] as u64)
            .sum();
        let avg = total as f64 / 19.0;
        avg_steps.push((AA_ORDER[i], avg));
    }
    avg_steps.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    println!("  Ranked by average step distance to all other amino acids:");
    println!("  (Lower = more connected, more reachable by single mutations)");
    println!();
    println!("  {:>4}  {:>4}  {:>8}", "AA", "Name", "Avg steps");
    println!("  ----  ----  --------");
    for (c, avg) in &avg_steps {
        let marker = if avg == &avg_steps[0].1 {
            " <-- most connected"
        } else if avg == &avg_steps.last().unwrap().1 {
            " <-- most isolated"
        } else {
            ""
        };
        println!(
            "  {:>4}  {:>4}  {:>8.3}{}",
            c,
            aa_full_name(*c),
            avg,
            marker
        );
    }
    println!();

    // ===================================================================
    // Part 4: Chemically interesting patterns
    // ===================================================================

    println!("PART 4: Chemically Interesting Patterns");
    println!("========================================");
    println!();
    println!("  The genetic code is not random -- chemically similar amino acids");
    println!("  tend to be encoded by similar codons. This means point mutations");
    println!("  often produce conservative amino acid changes.");
    println!();

    // Collect all 1-step pairs
    let mut one_step_pairs: Vec<(char, char)> = Vec::new();
    for i in 0..20 {
        for j in (i + 1)..20 {
            if matrix[i][j] == 1 {
                one_step_pairs.push((AA_ORDER[i], AA_ORDER[j]));
            }
        }
    }

    println!(
        "  1-step pairs ({} total) -- reachable by a single nucleotide change:",
        one_step_pairs.len()
    );
    for (idx, (a, b)) in one_step_pairs.iter().enumerate() {
        if idx % 6 == 0 {
            print!("    ");
        }
        print!("{}-{}({}-{})  ", a, b, aa_full_name(*a), aa_full_name(*b));
        if idx % 6 == 5 {
            println!();
        }
    }
    if one_step_pairs.len() % 6 != 0 {
        println!();
    }
    println!();

    // Highlight chemically notable 1-step pairs
    let notable_1step = [
        ('D', 'E', "both acidic (negative charge)"),
        ('I', 'V', "both branched-chain, hydrophobic"),
        ('F', 'L', "both large hydrophobic"),
        ('N', 'K', "both polar/charged, similar size"),
        ('S', 'T', "both small, hydroxyl-bearing"),
        ('A', 'V', "both small hydrophobic"),
    ];

    println!("  Chemically conservative 1-step pairs:");
    for (a, b, desc) in &notable_1step {
        let i = AA_ORDER.iter().position(|&c| c == *a).unwrap();
        let j = AA_ORDER.iter().position(|&c| c == *b).unwrap();
        println!(
            "    {}({}) <-> {}({})  cost={} step  -- {}",
            a,
            aa_full_name(*a),
            b,
            aa_full_name(*b),
            matrix[i][j],
            desc
        );
    }
    println!();

    // 3-step pairs
    let mut three_step_pairs: Vec<(char, char)> = Vec::new();
    for i in 0..20 {
        for j in (i + 1)..20 {
            if matrix[i][j] == 3 {
                three_step_pairs.push((AA_ORDER[i], AA_ORDER[j]));
            }
        }
    }
    println!(
        "  3-step pairs ({} total) -- maximum distance, all three codon positions differ:",
        three_step_pairs.len()
    );
    for (idx, (a, b)) in three_step_pairs.iter().enumerate() {
        if idx % 6 == 0 {
            print!("    ");
        }
        print!("{}-{}  ", a, b);
        if idx % 6 == 5 {
            println!();
        }
    }
    if three_step_pairs.len() % 6 != 0 {
        println!();
    }
    println!();

    // ===================================================================
    // Part 5: Comparison with randomized codes
    // ===================================================================

    println!("PART 5: Is the Genetic Code Optimized? Comparison with Random Codes");
    println!("===================================================================");
    println!();
    println!("  If the genetic code were a 'frozen accident' (Crick, 1968), the");
    println!("  physicochemical impact of point mutations would be typical of random");
    println!("  amino acid-to-codon assignments. But if natural selection has shaped");
    println!("  the code, single-nucleotide mutations should produce SMALLER changes");
    println!("  in amino acid properties than in random codes.");
    println!();
    println!("  Following Freeland & Hurst (1998), we measure the mean squared change");
    println!("  in POLAR REQUIREMENT (Woese, 1973) across all possible single-nucleotide");
    println!("  nonsynonymous mutations. Lower values mean the code better buffers");
    println!("  against the physicochemical impact of point mutations.");
    println!();
    println!("  We generate 1000 random permutations of the amino acid assignment");
    println!();

    // Build the codon-to-AA lookup table
    let codon_to_aa = build_codon_to_aa();
    let real_ms_polar = mean_squared_polar_change(&codon_to_aa);

    // Identity permutation
    let identity: [usize; 20] = {
        let mut arr = [0usize; 20];
        for i in 0..20 {
            arr[i] = i;
        }
        arr
    };

    let n_random = 1000;
    let mut rng = Rng::new(314159265);
    let mut random_avgs: Vec<f64> = Vec::with_capacity(n_random);

    for _ in 0..n_random {
        let mut perm = identity;
        rng.shuffle(&mut perm);
        let permuted_table = permute_codon_table(&codon_to_aa, &perm);
        let ms = mean_squared_polar_change(&permuted_table);
        random_avgs.push(ms);
    }

    random_avgs.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let random_mean: f64 = random_avgs.iter().sum::<f64>() / n_random as f64;
    let variance: f64 = random_avgs
        .iter()
        .map(|x| (x - random_mean) * (x - random_mean))
        .sum::<f64>()
        / n_random as f64;
    let random_sd = variance.sqrt();

    // Count how many random codes have average <= real code's average
    let n_leq = random_avgs
        .iter()
        .filter(|&&x| x <= real_ms_polar)
        .count();
    let percentile = 100.0 * n_leq as f64 / n_random as f64;

    println!(
        "  Real genetic code mean squared polar change: {:.4}",
        real_ms_polar
    );
    println!(
        "  Random codes: mean = {:.4}, std dev = {:.4}",
        random_mean, random_sd
    );
    println!(
        "  Random codes: min = {:.4}, max = {:.4}",
        random_avgs[0],
        random_avgs[n_random - 1]
    );
    println!();
    println!(
        "  Real code percentile: {:.1}% ({} of {} random codes have lower or equal value)",
        percentile, n_leq, n_random
    );
    if random_sd > 0.0 {
        println!(
            "  Z-score: {:.2} (standard deviations below the random mean)",
            (real_ms_polar - random_mean) / random_sd
        );
    }
    println!();

    // Simple histogram of random code values
    println!("  Distribution of random code mean squared polar changes:");
    println!();
    // Determine histogram range dynamically
    let hist_min = (random_avgs[0] * 10.0).floor() / 10.0;
    let hist_max_val = random_avgs[n_random - 1];
    let hist_max = ((hist_max_val * 10.0).ceil() / 10.0).max(hist_min + 1.0);
    let n_bins = 20;
    let bin_width = (hist_max - hist_min) / n_bins as f64;
    let mut bins = vec![0usize; n_bins];
    for &avg in &random_avgs {
        let bin = ((avg - hist_min) / bin_width) as usize;
        if bin < n_bins {
            bins[bin] += 1;
        }
    }
    let max_count = *bins.iter().max().unwrap_or(&1);

    for (i, &count) in bins.iter().enumerate() {
        let lo = hist_min + i as f64 * bin_width;
        let bar_len = if max_count > 0 {
            (count * 40) / max_count
        } else {
            0
        };
        let marker = if lo <= real_ms_polar && real_ms_polar < lo + bin_width {
            " <-- REAL CODE"
        } else {
            ""
        };
        print!("  {:.2}-{:.2} |", lo, lo + bin_width);
        for _ in 0..bar_len {
            print!("#");
        }
        println!(" ({}){}", count, marker);
    }
    println!();
    println!("  The real genetic code lies far to the LEFT of the distribution,");
    println!("  meaning it is significantly more mutation-buffered than chance.");
    println!("  This is strong evidence that the code has been shaped by natural");
    println!("  selection to minimize the phenotypic impact of point mutations.");
    println!();

    // ===================================================================
    // Part 6: Protein parsimony with the step matrix
    // ===================================================================

    println!("PART 6: Protein Parsimony in Action");
    println!("====================================");
    println!();
    println!("  The step matrix is used by the Sankoff algorithm for protein");
    println!("  parsimony. Unlike unweighted parsimony (where every amino acid");
    println!("  change costs 1), protein parsimony weights changes by their");
    println!("  minimum nucleotide cost in the genetic code.");
    println!();

    // Create a small protein alignment
    let alignment = ProteinAlignment::new(vec![
        ProteinSequence::new(
            "Human",
            "MVLSPADKTN".chars().map(AminoAcid::from_char).collect(),
        ),
        ProteinSequence::new(
            "Chimp",
            "MVLSPADKTN".chars().map(AminoAcid::from_char).collect(),
        ),
        ProteinSequence::new(
            "Gorilla",
            "MVLSPADKTN".chars().map(AminoAcid::from_char).collect(),
        ),
        ProteinSequence::new(
            "Mouse",
            "MVLSGEDKSN".chars().map(AminoAcid::from_char).collect(),
        ),
        ProteinSequence::new(
            "Chicken",
            "MVLSAADKTN".chars().map(AminoAcid::from_char).collect(),
        ),
    ])
    .expect("alignment should be valid");

    println!("  Example alignment (5 taxa, 10 sites):");
    println!("  Inspired by alpha-hemoglobin N-terminal region.");
    println!();
    for seq in &alignment.sequences {
        let residues: String = seq.residues.iter().map(|aa| aa.to_char()).collect();
        println!("    {:>10}: {}", seq.name, residues);
    }
    println!();

    // Score on two competing topologies
    let tree1_str = "((Human:0.01,Chimp:0.01):0.02,((Gorilla:0.02,Mouse:0.1):0.03,Chicken:0.15):0.01);";
    let tree2_str = "((Human:0.01,Mouse:0.01):0.02,((Chimp:0.02,Gorilla:0.1):0.03,Chicken:0.15):0.01);";
    let tree1 = parse_newick(tree1_str).expect("tree1 should parse");
    let tree2 = parse_newick(tree2_str).expect("tree2 should parse");

    let score1 = protein_parsimony_score(&tree1, &alignment);
    let score2 = protein_parsimony_score(&tree2, &alignment);

    println!("  Topology 1 (primates together): ((Human,Chimp),(Gorilla,Mouse),Chicken)");
    println!("    Protein parsimony score: {} steps", score1);
    println!();
    println!("  Topology 2 (wrong grouping):    ((Human,Mouse),(Chimp,Gorilla),Chicken)");
    println!("    Protein parsimony score: {} steps", score2);
    println!();

    if score1 < score2 {
        println!("  The correct topology (primates together) has FEWER steps.");
    } else if score1 == score2 {
        println!("  Both topologies have the same score on this small example.");
    } else {
        println!("  The wrong topology has fewer steps on this small example.");
    }
    println!();

    // Show the effect of weighted vs unweighted scoring
    println!("  Comparing weighted (genetic code) vs unweighted parsimony:");
    println!();
    println!("  At each variable site, the Sankoff algorithm considers the");
    println!("  minimum-cost path through the tree. Amino acid changes that");
    println!("  require only 1 nucleotide change (like A<->V) cost less than");
    println!("  changes requiring 2 or 3 nucleotide changes (like W<->D).");
    println!();

    // Show specific site differences
    let human = "MVLSPADKTN";
    let mouse = "MVLSGEDKSN";
    println!("  Site-by-site comparison (Human vs Mouse):");
    println!("  {:>4} {:>6} {:>6} {:>6}", "Site", "Human", "Mouse", "Steps");
    println!("  ---- ------ ------ -----");
    for (site, (h, m)) in human.chars().zip(mouse.chars()).enumerate() {
        if h != m {
            let hi = AA_ORDER.iter().position(|&c| c == h).unwrap();
            let mi = AA_ORDER.iter().position(|&c| c == m).unwrap();
            let cost = matrix[hi][mi];
            println!(
                "  {:>4} {:>4}({}) {:>4}({}) {:>5}",
                site + 1,
                h,
                aa_full_name(h),
                m,
                aa_full_name(m),
                cost
            );
        }
    }
    println!();
    println!("  These per-site costs reflect the STRUCTURE of the genetic code:");
    println!("  closely related amino acids cost fewer steps because their codons");
    println!("  differ at fewer positions.");
    println!();

    // ===================================================================
    // References
    // ===================================================================

    println!("References:");
    println!("  Sankoff, D. (1975). Minimal mutation trees of sequences.");
    println!("  SIAM Journal on Applied Mathematics, 28, 35-42.");
    println!();
    println!("  Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.");
    println!("  Chapter 8: The parsimony criterion.");
    println!();
    println!("  Crick, F.H.C. (1968). The origin of the genetic code.");
    println!("  Journal of Molecular Biology, 38, 367-379.");
    println!();
    println!("  Freeland, S.J. & Hurst, L.D. (1998). The genetic code is one in");
    println!("  a million. Journal of Molecular Evolution, 47, 238-248.");
    println!();
}
