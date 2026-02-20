#!/usr/bin/env python3
"""
JC69 sequence simulator for phylogenetics benchmarking.

Generates random tree topologies via iterative leaf addition, simulates
DNA sequences under the Jukes-Cantor (1969) model, and writes FASTA
alignments + true Newick trees.

No external dependencies — uses only Python's standard library.

JC69 transition probabilities for branch length t:
    P(same)      = 1/4 + 3/4 * exp(-4t/3)
    P(diff_each) = 1/4 - 1/4 * exp(-4t/3)
"""

from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------------
# Tree data structures
# ---------------------------------------------------------------------------

@dataclass
class Node:
    """A node in an unrooted phylogenetic tree (stored as rooted for simulation)."""
    label: str | None = None
    children: list[tuple[Node, float]] = field(default_factory=list)  # (child, branch_length)

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def newick(self) -> str:
        """Return Newick string for the subtree rooted at this node."""
        if self.is_leaf:
            return self.label or ""
        parts = [f"{child.newick()}:{bl:.8f}" for child, bl in self.children]
        inner = ",".join(parts)
        label = self.label or ""
        return f"({inner}){label}"


# ---------------------------------------------------------------------------
# Tree generation — birth-death-like iterative leaf addition
# ---------------------------------------------------------------------------

def _collect_edges(node: Node) -> list[tuple[Node, int, Node, float]]:
    """Collect all (parent, child_index, child, branch_length) edges."""
    edges: list[tuple[Node, int, Node, float]] = []
    for i, (child, bl) in enumerate(node.children):
        edges.append((node, i, child, bl))
        edges.extend(_collect_edges(child))
    return edges


def generate_tree(ntaxa: int, rng: random.Random, branch_scale: float = 0.1) -> Node:
    """
    Build a random binary tree by iterative taxon addition.

    Start with a root connecting two taxa. For each subsequent taxon,
    pick a random existing edge, break it with a new internal node, and
    attach the new leaf there. Random branch lengths are drawn from an
    exponential distribution and later rescaled so the mean root-to-tip
    distance is approximately `branch_scale`.
    """
    labels = [f"t{i+1:03d}" for i in range(ntaxa)]

    # Start: root with two children
    root = Node()
    bl1 = rng.expovariate(1.0)
    bl2 = rng.expovariate(1.0)
    root.children = [
        (Node(label=labels[0]), bl1),
        (Node(label=labels[1]), bl2),
    ]

    for idx in range(2, ntaxa):
        edges = _collect_edges(root)
        # Also include edges from root itself (already covered by _collect_edges)
        parent, ci, child, old_bl = rng.choice(edges)

        # Split point along the chosen edge
        frac = rng.random()
        bl_above = old_bl * frac
        bl_below = old_bl * (1.0 - frac)

        # New internal node
        new_internal = Node()
        new_internal.children = [
            (child, bl_below),
            (Node(label=labels[idx]), rng.expovariate(1.0)),
        ]

        # Replace the old child in the parent
        parent.children[ci] = (new_internal, bl_above)

    # Rescale branch lengths so mean root-to-tip ~ branch_scale
    _rescale_tree(root, branch_scale)
    return root


def _root_to_tip_distances(node: Node, depth: float = 0.0) -> list[float]:
    """Collect root-to-tip distances."""
    if node.is_leaf:
        return [depth]
    dists: list[float] = []
    for child, bl in node.children:
        dists.extend(_root_to_tip_distances(child, depth + bl))
    return dists


def _scale_branches(node: Node, factor: float) -> None:
    """Multiply all branch lengths by factor."""
    node.children = [
        (child, bl * factor) for child, bl in node.children
    ]
    for child, _ in node.children:
        _scale_branches(child, factor)


def _rescale_tree(root: Node, target_mean: float) -> None:
    """Rescale tree so mean root-to-tip distance equals target_mean."""
    dists = _root_to_tip_distances(root)
    current_mean = sum(dists) / len(dists) if dists else 1.0
    if current_mean > 0:
        _scale_branches(root, target_mean / current_mean)


# ---------------------------------------------------------------------------
# JC69 sequence simulation
# ---------------------------------------------------------------------------

BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}


def jc69_transition_probs(t: float) -> list[float]:
    """
    Return [P(same), P(diff_each)] for branch length t under JC69.

    P(same)      = 1/4 + 3/4 * exp(-4t/3)
    P(diff_each) = 1/4 - 1/4 * exp(-4t/3)
    """
    if t <= 0.0:
        return [1.0, 0.0]
    exp_term = math.exp(-4.0 * t / 3.0)
    p_same = 0.25 + 0.75 * exp_term
    p_diff = 0.25 - 0.25 * exp_term
    return [p_same, p_diff]


def evolve_sequence(parent_seq: list[int], branch_length: float, rng: random.Random) -> list[int]:
    """Evolve a sequence along a branch under JC69."""
    p_same, p_diff = jc69_transition_probs(branch_length)
    # CDF: [p_same, p_same+p_diff, p_same+2*p_diff, 1.0]
    # For each site, draw uniform and assign new base
    child_seq: list[int] = []
    for base in parent_seq:
        r = rng.random()
        if r < p_same:
            child_seq.append(base)
        else:
            # Pick one of the 3 different bases
            offset = r - p_same
            diff_idx = int(offset / p_diff) if p_diff > 0 else 0
            diff_idx = min(diff_idx, 2)  # clamp to 0, 1, 2
            # The 3 "other" bases in order
            others = [b for b in range(4) if b != base]
            child_seq.append(others[diff_idx])
    return child_seq


def simulate_sequences(
    root: Node,
    nsites: int,
    rng: random.Random,
) -> dict[str, str]:
    """
    Simulate DNA alignment on the tree under JC69.

    Generates a root sequence with equal base frequencies, then evolves
    down every branch. Returns {taxon_label: sequence_string}.
    """
    # Root sequence: uniform base frequencies
    root_seq = [rng.randint(0, 3) for _ in range(nsites)]

    sequences: dict[str, str] = {}

    def _recurse(node: Node, seq: list[int]) -> None:
        if node.is_leaf:
            sequences[node.label] = "".join(BASES[b] for b in seq)
            return
        for child, bl in node.children:
            child_seq = evolve_sequence(seq, bl, rng)
            _recurse(child, child_seq)

    _recurse(root, root_seq)
    return sequences


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Write sequences in FASTA format."""
    with open(path, "w") as f:
        for name, seq in sorted(sequences.items()):
            f.write(f">{name}\n")
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def write_newick(root: Node, path: Path) -> None:
    """Write rooted tree in Newick format."""
    with open(path, "w") as f:
        f.write(root.newick() + ";\n")


# ---------------------------------------------------------------------------
# Dataset matrix
# ---------------------------------------------------------------------------

DATASET_MATRIX: list[tuple[int, int]] = [
    (10, 500),
    (10, 1000),
    (20, 500),
    (20, 1000),
    (20, 5000),
    (50, 1000),
    (50, 5000),
    (100, 1000),
    (100, 5000),
    (200, 1000),
    (200, 5000),
    (500, 1000),
]

NUM_REPLICATES = 3
BRANCH_SCALE = 0.1
BASE_SEED = 42


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    script_dir = Path(__file__).resolve().parent
    out_dir = script_dir / "datasets"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest: list[dict] = []
    total = len(DATASET_MATRIX) * NUM_REPLICATES
    count = 0
    t0 = time.time()

    print(f"Generating {total} datasets in {out_dir}/")
    print(f"  Replicates per condition: {NUM_REPLICATES}")
    print(f"  Branch scale (mean root-to-tip): {BRANCH_SCALE}")
    print(f"  Base seed: {BASE_SEED}")
    print()

    for ntaxa, nsites in DATASET_MATRIX:
        for rep in range(1, NUM_REPLICATES + 1):
            count += 1
            seed = BASE_SEED + ntaxa * 10000 + nsites * 10 + rep
            rng = random.Random(seed)

            prefix = f"sim_{ntaxa}_{nsites}_rep{rep}"
            fasta_path = out_dir / f"{prefix}.fasta"
            tree_path = out_dir / f"{prefix}.true.nwk"

            t_start = time.time()
            sys.stdout.write(
                f"  [{count:3d}/{total}] {ntaxa:>4d} taxa, {nsites:>5d} sites, rep {rep} ... "
            )
            sys.stdout.flush()

            # Generate tree
            tree = generate_tree(ntaxa, rng, branch_scale=BRANCH_SCALE)

            # Simulate sequences
            sequences = simulate_sequences(tree, nsites, rng)

            # Write outputs
            write_fasta(sequences, fasta_path)
            write_newick(tree, tree_path)

            elapsed = time.time() - t_start

            # Collect stats
            dists = _root_to_tip_distances(tree)
            mean_depth = sum(dists) / len(dists) if dists else 0.0

            entry = {
                "prefix": prefix,
                "ntaxa": ntaxa,
                "nsites": nsites,
                "replicate": rep,
                "seed": seed,
                "branch_scale": BRANCH_SCALE,
                "mean_root_to_tip": round(mean_depth, 6),
                "fasta": f"{prefix}.fasta",
                "tree": f"{prefix}.true.nwk",
                "fasta_bytes": fasta_path.stat().st_size,
                "tree_bytes": tree_path.stat().st_size,
                "generation_time_sec": round(elapsed, 3),
            }
            manifest.append(entry)

            print(
                f"done ({elapsed:.2f}s, "
                f"fasta {entry['fasta_bytes']:,}B, "
                f"tree {entry['tree_bytes']:,}B)"
            )

    # Write manifest
    manifest_path = out_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(
            {
                "description": "Simulated JC69 datasets for phylogenetics benchmarking",
                "generator": "generate_data.py",
                "base_seed": BASE_SEED,
                "branch_scale": BRANCH_SCALE,
                "num_replicates": NUM_REPLICATES,
                "datasets": manifest,
            },
            f,
            indent=2,
        )
        f.write("\n")

    total_time = time.time() - t0
    total_fasta_bytes = sum(e["fasta_bytes"] for e in manifest)
    total_tree_bytes = sum(e["tree_bytes"] for e in manifest)

    print()
    print("=" * 60)
    print(f"  Datasets generated:  {total}")
    print(f"  Total FASTA size:    {total_fasta_bytes:,} bytes ({total_fasta_bytes / 1024 / 1024:.1f} MB)")
    print(f"  Total tree size:     {total_tree_bytes:,} bytes ({total_tree_bytes / 1024:.1f} KB)")
    print(f"  Manifest:            {manifest_path}")
    print(f"  Total time:          {total_time:.1f}s")
    print("=" * 60)


if __name__ == "__main__":
    main()
