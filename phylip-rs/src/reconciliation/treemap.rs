//! TREEMAP: host-parasite tree reconciliation via dynamic programming.
//!
//! Implements the reconciliation algorithm of Page (1994) for mapping a
//! parasite (symbiont) phylogeny onto a host phylogeny, given a known
//! leaf-level association.  The algorithm assigns each internal parasite
//! node to a host node and classifies the corresponding evolutionary event
//! (cospeciation, duplication, host-switch, or sorting/loss), minimising
//! total weighted cost.
//!
//! # Algorithm overview
//!
//! 1. **Leaf mapping**: each parasite leaf is mapped to its known host leaf
//!    via the user-supplied association `phi`.
//!
//! 2. **LCA reconciliation** (postorder on P): for each internal parasite
//!    node `p` with children `p1`, `p2`, the image `sigma(p)` is set to
//!    `LCA_H(sigma(p1), sigma(p2))`.  The event type is determined by
//!    whether `sigma(p1)` and `sigma(p2)` lie in different subtrees of
//!    `sigma(p)` (cospeciation) or coincide with `sigma(p)` (duplication).
//!
//! 3. **Sorting costs**: between `sigma(p)` and `sigma(child)` the number
//!    of sorting (lineage-loss) events equals `depth(sigma(child)) -
//!    depth(sigma(p)) - 1`.
//!
//! # References
//!
//! - Page, R.D.M. (1994). Maps between trees and cladistic analysis of
//!   historical associations among genes, organisms, and areas.
//!   *Systematic Biology*, 43:58-77.
//! - Page, R.D.M. (1994). Parallel phylogenies: reconstructing the history
//!   of host-parasite assemblages. *Cladistics*, 10:155-173.

use std::collections::HashMap;

use crate::tree::types::Tree;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Costs assigned to the four cophylogenetic event types.
///
/// Default costs follow the standard TREEMAP weighting: cospeciation is free,
/// duplication and sorting cost 1, and host-switching costs 2.
#[derive(Debug, Clone, Copy)]
pub struct EventCosts {
    /// Cost of a cospeciation event (default 0).
    pub cospeciation: f64,
    /// Cost of a duplication event (default 1).
    pub duplication: f64,
    /// Cost of a host-switching event (default 2).
    pub host_switch: f64,
    /// Cost of a sorting (lineage loss) event (default 1).
    pub sorting: f64,
}

impl Default for EventCosts {
    fn default() -> Self {
        Self {
            cospeciation: 0.0,
            duplication: 1.0,
            host_switch: 2.0,
            sorting: 1.0,
        }
    }
}

/// Classification of the evolutionary event at an internal parasite node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Event {
    /// Leaf node (no event).
    Leaf,
    /// Host and parasite speciated together.
    Cospeciation,
    /// Parasite speciated within the same host lineage.
    Duplication,
    /// Parasite lineage was lost from a host lineage.
    Sorting,
    /// Parasite jumped to a different host lineage.
    HostSwitch,
}

/// The complete result of a host-parasite reconciliation.
#[derive(Debug, Clone)]
pub struct ReconciliationResult {
    /// `mapping[p]` = host node index that parasite node `p` maps to.
    pub mapping: Vec<usize>,
    /// `events[p]` = the event classified at parasite node `p`.
    pub events: Vec<Event>,
    /// Total reconciliation cost (sum of all event costs).
    pub total_cost: f64,
    /// Number of cospeciation events.
    pub n_cospeciations: usize,
    /// Number of duplication events.
    pub n_duplications: usize,
    /// Number of sorting (lineage-loss) events.
    pub n_sortings: usize,
    /// Number of host-switch events.
    pub n_host_switches: usize,
}

// ---------------------------------------------------------------------------
// LCA computation
// ---------------------------------------------------------------------------

/// Compute the depth of every node in `tree` (root has depth 0).
fn compute_depths(tree: &Tree) -> Vec<usize> {
    let mut depths = vec![0usize; tree.num_nodes()];
    for &id in &tree.preorder() {
        if let Some(parent) = tree.nodes[id].parent {
            depths[id] = depths[parent] + 1;
        }
    }
    depths
}

/// Compute the Lowest Common Ancestor of two nodes in a rooted tree.
///
/// The algorithm walks the deeper node upward until both nodes are at the
/// same depth, then walks both upward in tandem until they meet.
///
/// # Panics
///
/// Panics if `node_a` or `node_b` is out of range for `tree`.
pub fn lca(tree: &Tree, node_a: usize, node_b: usize) -> usize {
    let depths = compute_depths(tree);
    lca_with_depths(tree, node_a, node_b, &depths)
}

/// LCA using pre-computed depth information (avoids redundant work during
/// reconciliation).
fn lca_with_depths(tree: &Tree, mut a: usize, mut b: usize, depths: &[usize]) -> usize {
    // Equalise depths.
    while depths[a] > depths[b] {
        a = tree.nodes[a].parent.expect("node has no parent during LCA walk");
    }
    while depths[b] > depths[a] {
        b = tree.nodes[b].parent.expect("node has no parent during LCA walk");
    }
    // Walk both up until they meet.
    while a != b {
        a = tree.nodes[a].parent.expect("node has no parent during LCA walk");
        b = tree.nodes[b].parent.expect("node has no parent during LCA walk");
    }
    a
}

/// Return `true` if `ancestor` is an ancestor of `descendant` (or they are
/// the same node) in the given tree, using pre-computed depths.
fn is_ancestor_of(
    tree: &Tree,
    ancestor: usize,
    descendant: usize,
    depths: &[usize],
) -> bool {
    if depths[ancestor] > depths[descendant] {
        return false;
    }
    let mut node = descendant;
    while depths[node] > depths[ancestor] {
        node = tree.nodes[node].parent.expect("node has no parent during ancestor check");
    }
    node == ancestor
}

// ---------------------------------------------------------------------------
// Reconciliation
// ---------------------------------------------------------------------------

/// Reconcile a parasite tree onto a host tree.
///
/// Given a host tree `host_tree`, a parasite tree `parasite_tree`, and a
/// leaf-level mapping `leaf_mapping` (parasite leaf name -> host leaf name),
/// compute the LCA reconciliation that maps every parasite node to a host
/// node and classifies each internal parasite node's event.
///
/// # Arguments
///
/// * `host_tree` -- the rooted host phylogeny.
/// * `parasite_tree` -- the rooted parasite phylogeny.
/// * `leaf_mapping` -- maps each parasite leaf name to its host leaf name.
///   Every parasite leaf must appear as a key; host leaves may be reused
///   (multiple parasites on one host).
/// * `costs` -- the event cost scheme.
///
/// # Errors
///
/// Returns `Err` if:
/// - A parasite leaf has no entry in `leaf_mapping`.
/// - A host leaf name in `leaf_mapping` does not exist in `host_tree`.
/// - Either tree is empty.
///
/// # References
///
/// Page, R.D.M. (1994). Parallel phylogenies: reconstructing the history
/// of host-parasite assemblages. *Cladistics*, 10:155-173.
pub fn reconcile(
    host_tree: &Tree,
    parasite_tree: &Tree,
    leaf_mapping: &HashMap<String, String>,
    costs: &EventCosts,
) -> Result<ReconciliationResult, String> {
    // --- Validate inputs ---------------------------------------------------
    if host_tree.num_nodes() == 0 {
        return Err("host tree is empty".into());
    }
    if parasite_tree.num_nodes() == 0 {
        return Err("parasite tree is empty".into());
    }

    // Build name -> node-id index for host leaves.
    let host_name_to_id: HashMap<&str, usize> = host_tree
        .leaves()
        .iter()
        .filter_map(|n| n.name.as_deref().map(|name| (name, n.id)))
        .collect();

    // Pre-compute host depths.
    let host_depths = compute_depths(host_tree);

    // --- Initialise mapping arrays -----------------------------------------
    let n_parasite = parasite_tree.num_nodes();
    let mut mapping: Vec<usize> = vec![0; n_parasite];
    let mut events: Vec<Event> = vec![Event::Leaf; n_parasite];

    // Map parasite leaves to host leaves.
    for leaf in parasite_tree.leaves() {
        let pname = leaf
            .name
            .as_deref()
            .ok_or_else(|| format!("parasite leaf node {} has no name", leaf.id))?;
        let hname = leaf_mapping
            .get(pname)
            .ok_or_else(|| format!("no host mapping for parasite leaf '{}'", pname))?;
        let &hid = host_name_to_id
            .get(hname.as_str())
            .ok_or_else(|| format!("host leaf '{}' not found in host tree", hname))?;
        mapping[leaf.id] = hid;
    }

    // --- Postorder reconciliation ------------------------------------------
    let mut n_cospeciations = 0usize;
    let mut n_duplications = 0usize;
    let mut n_sortings = 0usize;
    let n_host_switches = 0usize;
    let mut total_cost = 0.0f64;

    for &pid in &parasite_tree.postorder() {
        let pnode = &parasite_tree.nodes[pid];
        if pnode.is_leaf() {
            continue;
        }

        // Internal node: must have exactly 2 children for a fully resolved
        // (binary) tree.  We handle the general case by taking the LCA over
        // all children.
        let child_ids = &pnode.children;
        if child_ids.is_empty() {
            continue;
        }

        // Compute LCA of all children's host mappings.
        let mut ancestor = mapping[child_ids[0]];
        for &cid in &child_ids[1..] {
            ancestor = lca_with_depths(host_tree, ancestor, mapping[cid], &host_depths);
        }
        mapping[pid] = ancestor;

        // Classify the event (binary case: exactly two children).
        if child_ids.len() == 2 {
            let h1 = mapping[child_ids[0]];
            let h2 = mapping[child_ids[1]];
            let h_parent = ancestor;

            if h1 == h2 {
                // Both children map to the same host node: DUPLICATION.
                events[pid] = Event::Duplication;
                n_duplications += 1;
                total_cost += costs.duplication;
            } else if h1 == h_parent || h2 == h_parent {
                // One child maps to the LCA itself: check whether the other
                // is in a different subtree (host-switch) or same subtree
                // (duplication).
                //
                // If one child maps directly to the ancestor and the other
                // descends from the ancestor, this is a duplication: the
                // parasite speciated within the host lineage leading to
                // the ancestor, and one copy remained there while the other
                // descended.
                events[pid] = Event::Duplication;
                n_duplications += 1;
                total_cost += costs.duplication;
            } else {
                // Children map to different descendants of the LCA:
                // check whether they are in different subtrees
                // (cospeciation) or the same subtree (duplication).
                let lca_node = &host_tree.nodes[h_parent];
                if lca_node.children.len() >= 2 {
                    // Determine which subtree each child falls into.
                    let left_subtree = lca_node.children[0];
                    let right_subtree = lca_node.children[1];

                    let h1_in_left =
                        is_ancestor_of(host_tree, left_subtree, h1, &host_depths);
                    let h2_in_right =
                        is_ancestor_of(host_tree, right_subtree, h2, &host_depths);
                    let h1_in_right =
                        is_ancestor_of(host_tree, right_subtree, h1, &host_depths);
                    let h2_in_left =
                        is_ancestor_of(host_tree, left_subtree, h2, &host_depths);

                    if (h1_in_left && h2_in_right) || (h1_in_right && h2_in_left) {
                        // Children in different subtrees: COSPECIATION.
                        events[pid] = Event::Cospeciation;
                        n_cospeciations += 1;
                        total_cost += costs.cospeciation;
                    } else {
                        // Same subtree: DUPLICATION.
                        events[pid] = Event::Duplication;
                        n_duplications += 1;
                        total_cost += costs.duplication;
                    }
                } else {
                    // Degenerate host node (fewer than 2 children).
                    events[pid] = Event::Duplication;
                    n_duplications += 1;
                    total_cost += costs.duplication;
                }
            }

            // Count sorting events along each child lineage.
            // Sorting events = number of host internal nodes traversed
            // between sigma(p) and sigma(child), i.e.,
            //   depth(sigma(child)) - depth(sigma(p)) - 1
            // (clamped to zero for safety).
            for &cid in child_ids {
                let h_child = mapping[cid];
                let d_parent = host_depths[h_parent];
                let d_child = host_depths[h_child];
                if d_child > d_parent + 1 {
                    let n_sort = d_child - d_parent - 1;
                    n_sortings += n_sort;
                    total_cost += n_sort as f64 * costs.sorting;
                }
            }
        } else {
            // Multifurcating case: classify as duplication (conservative).
            events[pid] = Event::Duplication;
            n_duplications += 1;
            total_cost += costs.duplication;

            for &cid in child_ids {
                let h_child = mapping[cid];
                let d_parent = host_depths[ancestor];
                let d_child = host_depths[h_child];
                if d_child > d_parent + 1 {
                    let n_sort = d_child - d_parent - 1;
                    n_sortings += n_sort;
                    total_cost += n_sort as f64 * costs.sorting;
                }
            }
        }
    }

    Ok(ReconciliationResult {
        mapping,
        events,
        total_cost,
        n_cospeciations,
        n_duplications,
        n_sortings,
        n_host_switches,
    })
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::types::Tree;

    // -----------------------------------------------------------------------
    // Helper: build a balanced 3-taxon tree
    //
    //        root (0)
    //       /        \
    //    int (1)     C (4)
    //   /     \
    //  A (2)  B (3)
    // -----------------------------------------------------------------------
    fn make_three_taxon_tree(a: &str, b: &str, c: &str) -> Tree {
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let internal = tree.add_node(None, None, Some(root));
        tree.add_node(Some(a.into()), None, Some(internal)); // 2
        tree.add_node(Some(b.into()), None, Some(internal)); // 3
        tree.add_node(Some(c.into()), None, Some(root)); // 4
        tree.root = root;
        tree
    }

    // -----------------------------------------------------------------------
    // Helper: build a 4-taxon symmetric tree
    //
    //            root (0)
    //           /        \
    //       int1 (1)    int2 (2)
    //       /   \       /    \
    //     A (3) B (4) C (5)  D (6)
    // -----------------------------------------------------------------------
    fn make_four_taxon_symmetric(a: &str, b: &str, c: &str, d: &str) -> Tree {
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let int1 = tree.add_node(None, None, Some(root));
        let int2 = tree.add_node(None, None, Some(root));
        tree.add_node(Some(a.into()), None, Some(int1)); // 3
        tree.add_node(Some(b.into()), None, Some(int1)); // 4
        tree.add_node(Some(c.into()), None, Some(int2)); // 5
        tree.add_node(Some(d.into()), None, Some(int2)); // 6
        tree.root = root;
        tree
    }

    // -----------------------------------------------------------------------
    // Test 1: Perfect cospeciation -- identical topologies, mirrored mapping
    // -----------------------------------------------------------------------
    #[test]
    fn test_perfect_cospeciation_3_taxa() {
        let host = make_three_taxon_tree("H1", "H2", "H3");
        let parasite = make_three_taxon_tree("P1", "P2", "P3");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H2".into());
        phi.insert("P3".into(), "H3".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        assert_eq!(result.n_cospeciations, 2, "expected 2 cospeciation events");
        assert_eq!(result.n_duplications, 0, "expected 0 duplications");
        assert_eq!(result.n_sortings, 0, "expected 0 sorting events");
        assert_eq!(result.n_host_switches, 0, "expected 0 host-switches");
        assert!(
            (result.total_cost - 0.0).abs() < 1e-10,
            "perfect cospeciation should have zero cost"
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: Perfect cospeciation -- 4 taxa, symmetric tree
    // -----------------------------------------------------------------------
    #[test]
    fn test_perfect_cospeciation_4_taxa() {
        let host = make_four_taxon_symmetric("H1", "H2", "H3", "H4");
        let parasite = make_four_taxon_symmetric("P1", "P2", "P3", "P4");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H2".into());
        phi.insert("P3".into(), "H3".into());
        phi.insert("P4".into(), "H4".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        assert_eq!(result.n_cospeciations, 3, "expected 3 cospeciation events");
        assert_eq!(result.n_duplications, 0);
        assert_eq!(result.n_sortings, 0);
        assert!(result.total_cost.abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // Test 3: Single duplication -- two parasites on same host
    //
    // Host:       root - H1
    //                  \ H2
    //
    // Parasite:   root - P1
    //                  \ P2
    //
    // Both P1 and P2 map to H1 => duplication at parasite root
    // -----------------------------------------------------------------------
    #[test]
    fn test_single_duplication() {
        // Host: simple 2-leaf tree
        let mut host = Tree::new();
        let hroot = host.add_node(None, None, None);
        host.add_node(Some("H1".into()), None, Some(hroot));
        host.add_node(Some("H2".into()), None, Some(hroot));
        host.root = hroot;

        // Parasite: simple 2-leaf tree
        let mut parasite = Tree::new();
        let proot = parasite.add_node(None, None, None);
        parasite.add_node(Some("P1".into()), None, Some(proot));
        parasite.add_node(Some("P2".into()), None, Some(proot));
        parasite.root = proot;

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H1".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        assert_eq!(result.n_duplications, 1, "expected 1 duplication");
        assert_eq!(result.n_cospeciations, 0);
        assert!((result.total_cost - 1.0).abs() < 1e-10, "duplication costs 1");
    }

    // -----------------------------------------------------------------------
    // Test 4: Sorting events -- parasite missing from some host lineages
    //
    // Host (4-taxon symmetric):
    //            root
    //           /    \
    //        int1    int2
    //       /   \   /    \
    //      H1  H2  H3   H4
    //
    // Parasite (2-leaf):
    //        root
    //       /    \
    //      P1    P2
    //
    // P1 -> H1, P2 -> H3
    //
    // LCA of H1 and H3 = root.
    // Cospeciation at root; sorting H1->int1 (depth 2 - depth 0 - 1 = 1)
    // and H3->int2 (depth 2 - depth 0 - 1 = 1) = 2 sorting events.
    // -----------------------------------------------------------------------
    #[test]
    fn test_sorting_events() {
        let host = make_four_taxon_symmetric("H1", "H2", "H3", "H4");

        let mut parasite = Tree::new();
        let proot = parasite.add_node(None, None, None);
        parasite.add_node(Some("P1".into()), None, Some(proot));
        parasite.add_node(Some("P2".into()), None, Some(proot));
        parasite.root = proot;

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H3".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        assert_eq!(result.n_cospeciations, 1, "expected 1 cospeciation at root");
        assert_eq!(result.n_sortings, 2, "expected 2 sorting events");
        // cost = 0 (cospeciation) + 2 * 1 (sorting) = 2
        assert!(
            (result.total_cost - 2.0).abs() < 1e-10,
            "total cost should be 2.0, got {}",
            result.total_cost
        );
    }

    // -----------------------------------------------------------------------
    // Test 5: 3-taxon host + 3-taxon parasite with known reconciliation
    //
    // Host:        root(0)
    //             /       \
    //          int(1)    H3(4)
    //         /     \
    //       H1(2)  H2(3)
    //
    // Parasite:    root(0)
    //             /       \
    //          int(1)    P3(4)
    //         /     \
    //       P1(2)  P2(3)
    //
    // Mapping: P1->H1, P2->H3, P3->H2
    //
    // Parasite internal node int(1): LCA(H1, H3) = root(0)
    //   - H1 is under int(1) which is left child of root
    //   - H3 is right child of root
    //   => cospeciation at root(0)
    //   - sorting for P1: depth(H1)=2, depth(root)=0 => 2-0-1 = 1 sorting
    //   - sorting for P2: depth(H3)=1, depth(root)=0 => 1-0-1 = 0 sorting
    //
    // Parasite root(0): LCA(sigma(int)=root(0), H2) = root(0)
    //   - sigma(int)=root(0) which equals LCA => duplication
    //   - sorting for int: depth(root(0))=0, depth(root(0))=0 => 0
    //   - sorting for P3: depth(H2)=2, depth(root)=0 => 2-0-1 = 1 sorting
    // -----------------------------------------------------------------------
    #[test]
    fn test_3_taxon_known_reconciliation() {
        let host = make_three_taxon_tree("H1", "H2", "H3");
        let parasite = make_three_taxon_tree("P1", "P2", "P3");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H3".into());
        phi.insert("P3".into(), "H2".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        // Parasite int(1): maps to host root(0) via LCA(H1,H3)=root
        assert_eq!(result.mapping[1], 0, "parasite int should map to host root");
        // Parasite root(0): maps to host root(0) via LCA(root,H2)=root
        assert_eq!(result.mapping[0], 0, "parasite root should map to host root");

        // Events: 1 cospeciation (at parasite int) + 1 duplication (at parasite root)
        assert_eq!(result.n_cospeciations, 1);
        assert_eq!(result.n_duplications, 1);
        // Sorting: 1 (P1 path) + 0 (P2 path) + 0 (int path) + 1 (P3 path) = 2
        assert_eq!(result.n_sortings, 2);

        // cost = 0 + 1 + 2*1 = 3
        assert!(
            (result.total_cost - 3.0).abs() < 1e-10,
            "total cost should be 3.0, got {}",
            result.total_cost
        );
    }

    // -----------------------------------------------------------------------
    // Test 6: Mismatched tree sizes (more parasites than hosts)
    //
    // Host:  root - H1
    //             \ H2
    //
    // Parasite (3 taxa):
    //          root
    //         /    \
    //       int    P3
    //      /   \
    //     P1   P2
    //
    // P1->H1, P2->H1, P3->H2
    //
    // Parasite int: LCA(H1, H1) = H1 => duplication
    // Parasite root: LCA(H1, H2) = root => cospeciation
    // -----------------------------------------------------------------------
    #[test]
    fn test_mismatched_tree_sizes() {
        let mut host = Tree::new();
        let hroot = host.add_node(None, None, None);
        host.add_node(Some("H1".into()), None, Some(hroot));
        host.add_node(Some("H2".into()), None, Some(hroot));
        host.root = hroot;

        let parasite = make_three_taxon_tree("P1", "P2", "P3");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H1".into());
        phi.insert("P3".into(), "H2".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        assert_eq!(result.n_cospeciations, 1, "expected 1 cospeciation");
        assert_eq!(result.n_duplications, 1, "expected 1 duplication");
        assert_eq!(result.n_sortings, 0, "expected 0 sorting events");
        // cost = 0 (cospeciation) + 1 (duplication) = 1
        assert!(
            (result.total_cost - 1.0).abs() < 1e-10,
            "total cost should be 1.0, got {}",
            result.total_cost
        );
    }

    // -----------------------------------------------------------------------
    // Test 7: Custom event costs change total cost
    // -----------------------------------------------------------------------
    #[test]
    fn test_custom_event_costs() {
        // Same scenario as test_sorting_events: 1 cospeciation + 2 sortings
        let host = make_four_taxon_symmetric("H1", "H2", "H3", "H4");

        let mut parasite = Tree::new();
        let proot = parasite.add_node(None, None, None);
        parasite.add_node(Some("P1".into()), None, Some(proot));
        parasite.add_node(Some("P2".into()), None, Some(proot));
        parasite.root = proot;

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H3".into());

        // With default costs: 0 + 2*1 = 2
        let result_default =
            reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();
        assert!((result_default.total_cost - 2.0).abs() < 1e-10);

        // Increase sorting cost to 5
        let expensive_sorting = EventCosts {
            cospeciation: 0.0,
            duplication: 1.0,
            host_switch: 2.0,
            sorting: 5.0,
        };
        let result_expensive =
            reconcile(&host, &parasite, &phi, &expensive_sorting).unwrap();
        // 0 + 2*5 = 10
        assert!(
            (result_expensive.total_cost - 10.0).abs() < 1e-10,
            "expected cost 10.0 with sorting=5, got {}",
            result_expensive.total_cost
        );

        // Make cospeciation expensive too
        let all_expensive = EventCosts {
            cospeciation: 3.0,
            duplication: 1.0,
            host_switch: 2.0,
            sorting: 2.0,
        };
        let result_all =
            reconcile(&host, &parasite, &phi, &all_expensive).unwrap();
        // 3 (cospeciation) + 2*2 (sorting) = 7
        assert!(
            (result_all.total_cost - 7.0).abs() < 1e-10,
            "expected cost 7.0, got {}",
            result_all.total_cost
        );
    }

    // -----------------------------------------------------------------------
    // Test 8: Error -- missing leaf mapping
    // -----------------------------------------------------------------------
    #[test]
    fn test_error_missing_leaf_mapping() {
        let host = make_three_taxon_tree("H1", "H2", "H3");
        let parasite = make_three_taxon_tree("P1", "P2", "P3");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        // P2 and P3 missing!

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default());
        assert!(result.is_err(), "should error on missing leaf mapping");
    }

    // -----------------------------------------------------------------------
    // Test 9: Error -- host name not in tree
    // -----------------------------------------------------------------------
    #[test]
    fn test_error_host_not_found() {
        let host = make_three_taxon_tree("H1", "H2", "H3");
        let parasite = make_three_taxon_tree("P1", "P2", "P3");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H2".into());
        phi.insert("P3".into(), "NONEXISTENT".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default());
        assert!(result.is_err(), "should error when host leaf not found");
    }

    // -----------------------------------------------------------------------
    // Test 10: Error -- empty tree
    // -----------------------------------------------------------------------
    #[test]
    fn test_error_empty_tree() {
        let host = Tree::new();
        let parasite = make_three_taxon_tree("P1", "P2", "P3");
        let phi = HashMap::new();

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default());
        assert!(result.is_err(), "should error on empty host tree");
    }

    // -----------------------------------------------------------------------
    // Test 11: LCA correctness -- standalone tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_lca_same_node() {
        let tree = make_four_taxon_symmetric("A", "B", "C", "D");
        assert_eq!(lca(&tree, 3, 3), 3, "LCA of a node with itself is itself");
    }

    #[test]
    fn test_lca_siblings() {
        let tree = make_four_taxon_symmetric("A", "B", "C", "D");
        // A(3) and B(4) are siblings under int1(1)
        assert_eq!(lca(&tree, 3, 4), 1);
    }

    #[test]
    fn test_lca_across_subtrees() {
        let tree = make_four_taxon_symmetric("A", "B", "C", "D");
        // A(3) and C(5) are in different subtrees => LCA is root(0)
        assert_eq!(lca(&tree, 3, 5), 0);
    }

    #[test]
    fn test_lca_leaf_and_ancestor() {
        let tree = make_four_taxon_symmetric("A", "B", "C", "D");
        // A(3) and root(0)
        assert_eq!(lca(&tree, 3, 0), 0, "LCA of leaf and root is root");
    }

    // -----------------------------------------------------------------------
    // Test 12: Event classification at each node
    // -----------------------------------------------------------------------
    #[test]
    fn test_event_vector_perfect_cospeciation() {
        let host = make_four_taxon_symmetric("H1", "H2", "H3", "H4");
        let parasite = make_four_taxon_symmetric("P1", "P2", "P3", "P4");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H2".into());
        phi.insert("P3".into(), "H3".into());
        phi.insert("P4".into(), "H4".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        // Leaves should be Event::Leaf
        assert_eq!(result.events[3], Event::Leaf);
        assert_eq!(result.events[4], Event::Leaf);
        assert_eq!(result.events[5], Event::Leaf);
        assert_eq!(result.events[6], Event::Leaf);

        // Internal nodes should be Event::Cospeciation
        assert_eq!(result.events[0], Event::Cospeciation); // root
        assert_eq!(result.events[1], Event::Cospeciation); // int1
        assert_eq!(result.events[2], Event::Cospeciation); // int2
    }

    // -----------------------------------------------------------------------
    // Test 13: Mapping vector correctness
    // -----------------------------------------------------------------------
    #[test]
    fn test_mapping_vector_perfect_cospeciation() {
        let host = make_four_taxon_symmetric("H1", "H2", "H3", "H4");
        let parasite = make_four_taxon_symmetric("P1", "P2", "P3", "P4");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H2".into());
        phi.insert("P3".into(), "H3".into());
        phi.insert("P4".into(), "H4".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        // Leaf mappings
        assert_eq!(result.mapping[3], 3, "P1 -> H1");
        assert_eq!(result.mapping[4], 4, "P2 -> H2");
        assert_eq!(result.mapping[5], 5, "P3 -> H3");
        assert_eq!(result.mapping[6], 6, "P4 -> H4");

        // Internal mappings: parasite int1 -> host int1, etc.
        assert_eq!(result.mapping[1], 1, "parasite int1 -> host int1");
        assert_eq!(result.mapping[2], 2, "parasite int2 -> host int2");
        assert_eq!(result.mapping[0], 0, "parasite root -> host root");
    }

    // -----------------------------------------------------------------------
    // Test 14: Multiple parasites on one host (widespread duplication)
    //
    // Host:  root - H1
    //             \ H2
    //
    // Parasite (3 taxa, all on H1):
    //          root
    //         /    \
    //       int    P3
    //      /   \
    //     P1   P2
    //
    // All parasites on H1 => all internal nodes are duplications.
    // -----------------------------------------------------------------------
    #[test]
    fn test_all_duplications() {
        let mut host = Tree::new();
        let hroot = host.add_node(None, None, None);
        host.add_node(Some("H1".into()), None, Some(hroot));
        host.add_node(Some("H2".into()), None, Some(hroot));
        host.root = hroot;

        let parasite = make_three_taxon_tree("P1", "P2", "P3");

        let mut phi = HashMap::new();
        phi.insert("P1".into(), "H1".into());
        phi.insert("P2".into(), "H1".into());
        phi.insert("P3".into(), "H1".into());

        let result = reconcile(&host, &parasite, &phi, &EventCosts::default()).unwrap();

        assert_eq!(result.n_duplications, 2, "all internal nodes are duplications");
        assert_eq!(result.n_cospeciations, 0);
        assert_eq!(result.n_sortings, 0);
        // cost = 2 * 1 = 2
        assert!((result.total_cost - 2.0).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // Test 15: Default costs are correct
    // -----------------------------------------------------------------------
    #[test]
    fn test_default_event_costs() {
        let costs = EventCosts::default();
        assert!((costs.cospeciation - 0.0).abs() < 1e-10);
        assert!((costs.duplication - 1.0).abs() < 1e-10);
        assert!((costs.host_switch - 2.0).abs() < 1e-10);
        assert!((costs.sorting - 1.0).abs() < 1e-10);
    }
}
