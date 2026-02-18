// kirchhoff_contrasts.rs
//
// Demonstrates that Felsenstein's (1985) independent contrasts algorithm
// computes the exact same quantities as solving a resistor network using
// Kirchhoff's laws. This is not a metaphor -- the mathematics are identical.
//
// A phylogenetic tree with branch lengths is isomorphic to an electrical
// circuit where each branch is a resistor. Branch lengths are resistances.
// Trait values at the tips are voltages. The independent contrasts algorithm
// computes "currents" (contrasts) and "effective resistances" (variances)
// using the exact same parallel-resistance formulas used in circuit theory.
//
// Usage:
//   cargo run --example kirchhoff_contrasts

use phylip_rs::comparative::contrasts::{independent_contrasts, ContrastResult};
use phylip_rs::comparative::ContinuousData;
use phylip_rs::tree::newick::parse_newick;

// ---------------------------------------------------------------------------
// Kirchhoff's circuit computation (manual, from first principles)
// ---------------------------------------------------------------------------

/// Represents a node in our manual circuit computation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct CircuitNode {
    name: String,
    voltage: f64,          // Trait value ("voltage") at this node
    resistance: f64,       // Effective resistance from this node upward
    contrast: f64,         // "Current" flowing at this junction
    std_contrast: f64,     // Standardized contrast (current / sqrt(resistance))
}

/// Compute the circuit equivalent for two children being joined.
///
/// When two resistors R_L and R_R are connected in series (child branches)
/// and their junction is connected to a parent, the effective circuit
/// values at the junction are:
///
///   Effective voltage (weighted mean):
///     V_node = (V_L * R_R + V_R * R_L) / (R_L + R_R)
///
///   This is the "voltage divider" -- the junction voltage is the
///   resistance-weighted average of the child voltages, just as a
///   voltage divider circuit produces a weighted mean of its inputs.
///
///   Contrast ("current" across the junction):
///     I = (V_L - V_R) / sqrt(R_L + R_R)
///
///   This is Ohm's law: current = voltage difference / resistance.
///   The total resistance between L and R is R_L + R_R (series).
///
///   Effective resistance added to parent branch:
///     R_additional = (R_L * R_R) / (R_L + R_R)
///
///   This is the parallel resistance formula! When two resistors
///   are "used up" to form a contrast, the remaining variance
///   added to the parent branch equals their parallel combination.
fn compute_junction(
    name: &str,
    left: &CircuitNode,
    right: &CircuitNode,
    branch_to_parent: f64,
) -> (CircuitNode, f64) {
    let r_l = left.resistance;
    let r_r = right.resistance;
    let v_l = left.voltage;
    let v_r = right.voltage;

    let r_total = r_l + r_r;

    // Voltage divider: weighted mean
    let v_node = (v_l * r_r + v_r * r_l) / r_total;

    // Current (raw contrast)
    let contrast = v_l - v_r;

    // Standardized contrast (Ohm's law)
    let std_contrast = contrast / r_total.sqrt();

    // Parallel resistance formula
    let r_parallel = (r_l * r_r) / r_total;

    // Effective resistance of this node = branch to parent + parallel component
    let r_effective = branch_to_parent + r_parallel;

    let node = CircuitNode {
        name: name.to_string(),
        voltage: v_node,
        resistance: r_effective,
        contrast,
        std_contrast,
    };

    (node, r_total)
}

fn main() {
    println!();
    println!("KIRCHHOFF'S CONTRASTS: Independent Contrasts = Resistor Networks");
    println!("=================================================================");
    println!();
    println!("Felsenstein's (1985) independent contrasts algorithm is mathematically");
    println!("identical to solving a resistor network using Kirchhoff's laws. This");
    println!("demonstration proves this equivalence with concrete numerical results.");
    println!();
    println!("The mapping between phylogenetics and circuit theory:");
    println!();
    println!("  Phylogenetics             Circuit Theory");
    println!("  ----------------------    ----------------------");
    println!("  Branch length             Resistance (Ohms)");
    println!("  Trait value at tip        Voltage source (Volts)");
    println!("  Weighted average          Voltage divider");
    println!("  Contrast                  Current (Amperes)");
    println!("  Standardized contrast     Current / sqrt(R)");
    println!("  Variance propagation      Parallel resistance");
    println!();

    // =========================================================================
    // Define the tree and trait data
    // =========================================================================
    println!("=================================================================");
    println!("  THE TREE AND DATA");
    println!("=================================================================");
    println!();

    let tree_newick = "((A:0.5,B:0.3):0.2,(C:0.4,(D:0.1,E:0.6):0.3):0.1);";
    let tree = parse_newick(tree_newick).expect("newick should parse");

    println!("  Tree: ((A:0.5, B:0.3):0.2, (C:0.4, (D:0.1, E:0.6):0.3):0.1)");
    println!();
    println!("  As a circuit diagram:");
    println!();
    println!("         0.5 Ohms      0.3 Ohms");
    println!("    A ---/\\/\\/\\---+---/\\/\\/\\--- B");
    println!("                  |");
    println!("             0.2 Ohms");
    println!("                  |");
    println!("               [ROOT]");
    println!("                  |");
    println!("             0.1 Ohms");
    println!("                  |");
    println!("    C ---/\\/\\/\\---+---/\\/\\/\\---+---/\\/\\/\\--- D");
    println!("         0.4 Ohms     0.3 Ohms |   0.1 Ohms");
    println!("                               |");
    println!("                          0.6 Ohms");
    println!("                               |");
    println!("                               E");
    println!();

    let tip_values: [(&str, f64); 5] = [
        ("A", 2.0),
        ("B", 4.0),
        ("C", 1.0),
        ("D", 5.0),
        ("E", 3.0),
    ];

    println!("  Tip voltages (trait values):");
    for (name, val) in &tip_values {
        println!("    {} = {:.1} V", name, val);
    }
    println!();

    // =========================================================================
    // Run Felsenstein's independent contrasts
    // =========================================================================
    println!("=================================================================");
    println!("  FELSENSTEIN'S INDEPENDENT CONTRASTS (phylip-rs computation)");
    println!("=================================================================");
    println!();

    let data = ContinuousData::new(
        tip_values.iter().map(|(n, _)| n.to_string()).collect(),
        tip_values.iter().map(|(_, v)| vec![*v]).collect(),
    )
    .expect("data should be valid");

    let result: ContrastResult =
        independent_contrasts(&tree, &data).expect("contrasts should succeed");

    println!("  Felsenstein's algorithm produces {} contrasts for {} taxa:",
             result.contrasts[0].len(), tip_values.len());
    println!();
    for (i, (&contrast, &variance)) in result
        .contrasts[0]
        .iter()
        .zip(result.variances.iter())
        .enumerate()
    {
        println!(
            "    Contrast {}: standardized = {:>8.6}, variance = {:.6}",
            i + 1,
            contrast,
            variance
        );
    }
    println!("    Root value: {:.6}", result.root_values[0]);
    println!();

    // =========================================================================
    // Manual Kirchhoff circuit computation
    // =========================================================================
    println!("=================================================================");
    println!("  KIRCHHOFF CIRCUIT COMPUTATION (manual, step by step)");
    println!("=================================================================");
    println!();

    // Leaf nodes (voltage sources)
    let leaf_a = CircuitNode {
        name: "A".to_string(),
        voltage: 2.0,
        resistance: 0.5, // branch length to parent
        contrast: 0.0,
        std_contrast: 0.0,
    };
    let leaf_b = CircuitNode {
        name: "B".to_string(),
        voltage: 4.0,
        resistance: 0.3,
        contrast: 0.0,
        std_contrast: 0.0,
    };
    let leaf_c = CircuitNode {
        name: "C".to_string(),
        voltage: 1.0,
        resistance: 0.4,
        contrast: 0.0,
        std_contrast: 0.0,
    };
    let leaf_d = CircuitNode {
        name: "D".to_string(),
        voltage: 5.0,
        resistance: 0.1,
        contrast: 0.0,
        std_contrast: 0.0,
    };
    let leaf_e = CircuitNode {
        name: "E".to_string(),
        voltage: 3.0,
        resistance: 0.6,
        contrast: 0.0,
        std_contrast: 0.0,
    };

    // We need to process in postorder matching the tree structure.
    // The tree is: ((A:0.5,B:0.3):0.2,(C:0.4,(D:0.1,E:0.6):0.3):0.1)
    //
    // Postorder traversal visits:
    //   1. Join D and E -> node_DE (branch to parent = 0.3)
    //   2. Join A and B -> node_AB (branch to parent = 0.2)
    //   3. Join C and node_DE -> node_CDE (branch to parent = 0.1)
    //   4. Join node_AB and node_CDE -> ROOT

    // The order of contrasts depends on postorder traversal. Let's figure out
    // the tree structure by examining what parse_newick produces.
    // From the newick: ((A:0.5,B:0.3):0.2,(C:0.4,(D:0.1,E:0.6):0.3):0.1)
    // The postorder should be: A, B, [AB], D, E, [DE], C, [CDE], [ROOT]
    // Or possibly: A, B, [AB], C, D, E, [DE], [CDE], [ROOT]

    // Let's compute all possible junction orderings and match to Felsenstein.
    // Actually, the postorder is determined by the tree structure, and contrasts
    // are computed at internal nodes in postorder. Let's compute both orderings
    // and figure out which one matches.

    println!("  Processing nodes in postorder (bottom-up), computing circuit values:");
    println!();

    // Junction 1: A + B
    println!("  --- Junction 1: Join A and B ---");
    println!("    R_A = {:.1} Ohms (branch length A)      V_A = {:.1} V", leaf_a.resistance, leaf_a.voltage);
    println!("    R_B = {:.1} Ohms (branch length B)      V_B = {:.1} V", leaf_b.resistance, leaf_b.voltage);
    println!();
    let (node_ab, var_ab) = compute_junction("AB", &leaf_a, &leaf_b, 0.2);
    println!("    Voltage divider (weighted average):");
    println!("      V_AB = (V_A * R_B + V_B * R_A) / (R_A + R_B)");
    println!("           = ({:.1} * {:.1} + {:.1} * {:.1}) / ({:.1} + {:.1})",
             leaf_a.voltage, leaf_b.resistance, leaf_b.voltage, leaf_a.resistance,
             leaf_a.resistance, leaf_b.resistance);
    println!("           = {:.6}", node_ab.voltage);
    println!();
    println!("    Current (contrast):");
    println!("      I = V_A - V_B = {:.1} - {:.1} = {:.6}", leaf_a.voltage, leaf_b.voltage, node_ab.contrast);
    println!("      I_std = I / sqrt(R_A + R_B) = {:.1} / sqrt({:.1}) = {:.6}",
             node_ab.contrast, var_ab, node_ab.std_contrast);
    println!();
    let r_parallel_ab = (leaf_a.resistance * leaf_b.resistance) / (leaf_a.resistance + leaf_b.resistance);
    println!("    Parallel resistance (variance propagation):");
    println!("      R_parallel = (R_A * R_B) / (R_A + R_B) = ({:.1} * {:.1}) / {:.1} = {:.6}",
             leaf_a.resistance, leaf_b.resistance, var_ab, r_parallel_ab);
    println!("      R_effective_AB = branch_to_parent + R_parallel = {:.1} + {:.6} = {:.6}",
             0.2, r_parallel_ab, node_ab.resistance);
    println!();

    // Junction 2: D + E
    println!("  --- Junction 2: Join D and E ---");
    println!("    R_D = {:.1} Ohms      V_D = {:.1} V", leaf_d.resistance, leaf_d.voltage);
    println!("    R_E = {:.1} Ohms      V_E = {:.1} V", leaf_e.resistance, leaf_e.voltage);
    println!();
    let (node_de, var_de) = compute_junction("DE", &leaf_d, &leaf_e, 0.3);
    println!("    Voltage divider:");
    println!("      V_DE = ({:.1} * {:.1} + {:.1} * {:.1}) / ({:.1} + {:.1}) = {:.6}",
             leaf_d.voltage, leaf_e.resistance, leaf_e.voltage, leaf_d.resistance,
             leaf_d.resistance, leaf_e.resistance, node_de.voltage);
    println!("    Current:");
    println!("      I = {:.1} - {:.1} = {:.6}", leaf_d.voltage, leaf_e.voltage, node_de.contrast);
    println!("      I_std = {:.6} / sqrt({:.6}) = {:.6}", node_de.contrast, var_de, node_de.std_contrast);
    let r_parallel_de = (leaf_d.resistance * leaf_e.resistance) / (leaf_d.resistance + leaf_e.resistance);
    println!("    Parallel resistance:");
    println!("      R_parallel = ({:.1} * {:.1}) / {:.1} = {:.6}",
             leaf_d.resistance, leaf_e.resistance, var_de, r_parallel_de);
    println!("      R_effective_DE = {:.1} + {:.6} = {:.6}", 0.3, r_parallel_de, node_de.resistance);
    println!();

    // Junction 3: C + DE
    println!("  --- Junction 3: Join C and DE ---");
    println!("    R_C = {:.1} Ohms            V_C = {:.1} V", leaf_c.resistance, leaf_c.voltage);
    println!("    R_DE = {:.6} Ohms      V_DE = {:.6} V", node_de.resistance, node_de.voltage);
    println!();
    let (node_cde, var_cde) = compute_junction("CDE", &leaf_c, &node_de, 0.1);
    println!("    Voltage divider:");
    println!("      V_CDE = ({:.1} * {:.6} + {:.6} * {:.1}) / ({:.1} + {:.6}) = {:.6}",
             leaf_c.voltage, node_de.resistance, node_de.voltage, leaf_c.resistance,
             leaf_c.resistance, node_de.resistance, node_cde.voltage);
    println!("    Current:");
    println!("      I = {:.1} - {:.6} = {:.6}", leaf_c.voltage, node_de.voltage, node_cde.contrast);
    println!("      I_std = {:.6} / sqrt({:.6}) = {:.6}", node_cde.contrast, var_cde, node_cde.std_contrast);
    let r_parallel_cde = (leaf_c.resistance * node_de.resistance) / (leaf_c.resistance + node_de.resistance);
    println!("    Parallel resistance:");
    println!("      R_parallel = ({:.1} * {:.6}) / {:.6} = {:.6}",
             leaf_c.resistance, node_de.resistance, var_cde, r_parallel_cde);
    println!("      R_effective_CDE = {:.1} + {:.6} = {:.6}", 0.1, r_parallel_cde, node_cde.resistance);
    println!();

    // Junction 4: AB + CDE (root)
    println!("  --- Junction 4: Join AB and CDE (ROOT) ---");
    println!("    R_AB = {:.6} Ohms      V_AB = {:.6} V", node_ab.resistance, node_ab.voltage);
    println!("    R_CDE = {:.6} Ohms     V_CDE = {:.6} V", node_cde.resistance, node_cde.voltage);
    println!();
    let (root, var_root) = compute_junction("ROOT", &node_ab, &node_cde, 0.0);
    println!("    Voltage divider:");
    println!("      V_ROOT = ({:.6} * {:.6} + {:.6} * {:.6}) / ({:.6} + {:.6})",
             node_ab.voltage, node_cde.resistance, node_cde.voltage, node_ab.resistance,
             node_ab.resistance, node_cde.resistance);
    println!("             = {:.6}", root.voltage);
    println!("    Current:");
    println!("      I = {:.6} - {:.6} = {:.6}", node_ab.voltage, node_cde.voltage, root.contrast);
    println!("      I_std = {:.6} / sqrt({:.6}) = {:.6}", root.contrast, var_root, root.std_contrast);
    println!();

    // =========================================================================
    // Side-by-side comparison
    // =========================================================================
    println!("=================================================================");
    println!("  SIDE-BY-SIDE COMPARISON: Felsenstein vs Kirchhoff");
    println!("=================================================================");
    println!();

    // The contrasts from Felsenstein are computed in postorder. We need to match
    // the order. The tree postorder visits internal nodes as:
    // 1. AB junction, 2. DE junction, 3. CDE junction, 4. ROOT
    // So contrasts[0] = AB, contrasts[1] = DE, contrasts[2] = CDE, contrasts[3] = ROOT

    let felsenstein_contrasts = &result.contrasts[0];
    let felsenstein_variances = &result.variances;

    // The circuit contrasts in our computation order
    let circuit_contrasts = [node_ab.std_contrast, node_de.std_contrast, node_cde.std_contrast, root.std_contrast];
    let circuit_variances = [var_ab, var_de, var_cde, var_root];
    let circuit_voltages = [node_ab.voltage, node_de.voltage, node_cde.voltage, root.voltage];
    let junction_names = ["A vs B", "D vs E", "C vs DE", "AB vs CDE"];

    println!(
        "  {:>12} | {:>16} {:>16} | {:>7}",
        "Junction", "Felsenstein", "Kirchhoff", "Match?"
    );
    println!(
        "  -------------|{:-<34}|--------",
        ""
    );

    println!();
    println!("  Standardized contrasts (= currents):");

    let mut all_match = true;
    for i in 0..4 {
        let f_val = felsenstein_contrasts[i];
        let k_val = circuit_contrasts[i];
        let matches = (f_val - k_val).abs() < 1e-8;
        if !matches {
            all_match = false;
        }
        println!(
            "  {:>12} | {:>16.10} {:>16.10} | {}",
            junction_names[i],
            f_val,
            k_val,
            if matches { "  YES" } else { "  NO" }
        );
    }

    println!();
    println!("  Variances (= total resistance at junction):");

    for i in 0..4 {
        let f_val = felsenstein_variances[i];
        let k_val = circuit_variances[i];
        let matches = (f_val - k_val).abs() < 1e-8;
        if !matches {
            all_match = false;
        }
        println!(
            "  {:>12} | {:>16.10} {:>16.10} | {}",
            junction_names[i],
            f_val,
            k_val,
            if matches { "  YES" } else { "  NO" }
        );
    }

    println!();
    println!("  Root value (= final voltage at root):");
    let root_match = (result.root_values[0] - root.voltage).abs() < 1e-8;
    if !root_match {
        all_match = false;
    }
    println!(
        "  {:>12} | {:>16.10} {:>16.10} | {}",
        "Root",
        result.root_values[0],
        root.voltage,
        if root_match { "  YES" } else { "  NO" }
    );

    println!();
    if all_match {
        println!("  ALL VALUES MATCH TO 8+ DECIMAL PLACES.");
        println!("  Felsenstein's algorithm and Kirchhoff's laws produce IDENTICAL results.");
    } else {
        println!("  Some values differ (check postorder traversal ordering).");
    }
    println!();

    // =========================================================================
    // The physics of why this works
    // =========================================================================
    println!("=================================================================");
    println!("  WHY THIS WORKS: The Mathematics Are Identical");
    println!("=================================================================");
    println!();
    println!("  At each internal node, Felsenstein's algorithm computes:");
    println!();
    println!("    X_parent = (X_L * v_R + X_R * v_L) / (v_L + v_R)");
    println!("    contrast = (X_L - X_R) / sqrt(v_L + v_R)");
    println!("    v_additional = (v_L * v_R) / (v_L + v_R)");
    println!();
    println!("  In a resistor network, Kirchhoff's laws give:");
    println!();
    println!("    V_junction = (V_1 * R_2 + V_2 * R_1) / (R_1 + R_2)   [voltage divider]");
    println!("    I = (V_1 - V_2) / (R_1 + R_2)                        [Ohm's law, series]");
    println!("    R_parallel = (R_1 * R_2) / (R_1 + R_2)               [parallel formula]");
    println!();
    println!("  These are the SAME equations with different variable names:");
    println!();
    println!("    v_L <-> R_1      (branch length = resistance)");
    println!("    v_R <-> R_2");
    println!("    X   <-> V        (trait value = voltage)");
    println!();
    println!("  The parallel resistance formula R_parallel = (R1*R2)/(R1+R2) appears");
    println!("  because when two branches of variance v_L and v_R are 'resolved' by");
    println!("  computing a contrast, the remaining variance to propagate upward is");
    println!("  their harmonic mean -- exactly the parallel combination of two resistors.");
    println!();
    println!("  This is not a coincidence. Both problems involve propagating information");
    println!("  through a tree-structured network where the 'signal' at each junction");
    println!("  is a weighted average, with weights inversely proportional to the");
    println!("  'path length' (resistance / branch length) from each source.");
    println!();

    // =========================================================================
    // Numerical verification summary
    // =========================================================================
    println!("=================================================================");
    println!("  NUMERICAL SUMMARY");
    println!("=================================================================");
    println!();
    println!("  Tree: ((A:0.5,B:0.3):0.2,(C:0.4,(D:0.1,E:0.6):0.3):0.1)");
    println!("  Tips: A=2.0, B=4.0, C=1.0, D=5.0, E=3.0");
    println!();
    println!("  Junction  | Contrast (I_std) | Variance (R_total) | Voltage (V)");
    println!("  ----------|------------------|--------------------|-----------");
    for i in 0..4 {
        println!(
            "  {:>9} | {:>16.10} | {:>18.10} | {:>10.6}",
            junction_names[i],
            circuit_contrasts[i],
            circuit_variances[i],
            circuit_voltages[i],
        );
    }
    println!();
    println!("  Root voltage = {:.10}", root.voltage);
    println!("  (This is the BLUE estimate of the ancestral trait value,", );
    println!("   Best Linear Unbiased Estimate under Brownian motion.)");
    println!();

    // =========================================================================
    // References
    // =========================================================================
    println!("=================================================================");
    println!("  REFERENCES");
    println!("=================================================================");
    println!();
    println!("  Felsenstein, J. (1985). Phylogenies and the comparative method.");
    println!("  American Naturalist, 125, 1-15.");
    println!();
    println!("  Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.");
    println!("  Chapter 25: Quantitative characters, phylogenies, and morphometrics.");
    println!();
    println!("  Kirchhoff, G. (1847). Ueber die Auflosung der Gleichungen, auf");
    println!("  welche man bei der Untersuchung der linearen Vertheilung");
    println!("  galvanischer Strome gefuehrt wird. Annalen der Physik, 148, 497-508.");
    println!();
}
