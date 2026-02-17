# CONSENSE Analysis - Consensus Tree Methods

**Source file**: `phylip-3.698/src/consense.c` (443 lines)
**Core logic**: `phylip-3.698/src/cons.c` (the actual consensus algorithms)
**Dependencies**: `phylip.h`, `cons.h`

## Overview

CONSENSE computes consensus trees from a set of input trees. A consensus tree
summarizes the agreement among multiple phylogenetic trees, typically bootstrap
replicates or trees from different methods. The program supports four consensus
methods:

1. **Strict consensus** (`strict`): Include only splits present in ALL input trees
2. **Majority rule (MR)** (`mr`): Include splits present in >50% of input trees
3. **Extended majority rule (MRe)** (`mre`, default): Majority rule plus
   additional compatible splits in decreasing frequency order
4. **Ml consensus** (`ml`): Include splits present in >= l fraction of trees
   (user-specified threshold)

The program is architecturally split: `consense.c` handles I/O and options,
while `cons.c` contains the actual consensus algorithms. `cons.c` is also
used by other programs (treedist, etc.), making it PHYLIP's shared consensus
computation library.

## Key Data Structures

### Bipartition (Split) Representation
```c
group_type **grouping;   // grouping[group_index][word] - bit-encoded splits
group_type **grping2;    // secondary grouping storage
group_type **group2;     // another grouping array
group_type *fullset;     // the full taxon set (all bits set)
```

Each bipartition (split) is represented as a bit vector. With n taxa, the
bit vector has `ceil(n / (8*sizeof(long)))` words. Taxon i is in the split
if bit i is set. By convention, the side containing taxon 1 is always
represented (the other side is the complement).

### Pattern Hash Table (from cons.h)
```c
typedef struct pattern_elm {
    group_type *apattern;    // the bipartition bit vector
    long *patternsize;       // number of taxa in this side
    double *length;          // branch length information
} pattern_elm;
```

The `pattern_array` is a hash table that maps bipartitions to their
frequency counts. This allows O(1) lookup of whether a given split has
been seen before.

### Tree Node Array (lines 395-403)
```c
pointarray nodep;   // nodep[0..2*(1+spp)-1] - all nodes
```

Allocated with 2*(1+spp) entries to accommodate the consensus tree, which
may have up to spp-1 internal nodes plus spp tips. Each tip node stores
the taxon name in `node->nayme[]`.

### Key Global Variables (from cons.c, declared extern)
- `ntrees`: total number of input trees (as double)
- `maxgrp`: maximum number of groups/splits (initial hash table size: 32767)
- `lasti`: index tracking for tree construction
- `tree_pairing`: set to `NO_PAIRING` for consensus (vs. tree comparison)
- `noroot`: whether input trees are unrooted
- `strict`, `mr`, `mre`, `ml`: consensus type flags (mutually exclusive)
- `mlfrac`: threshold fraction for Ml consensus (default 0.5)

## Core Algorithms

### 1. read_groups() (in cons.c)

Called from `main()` (line 392):
```c
read_groups(&pattern_array, trees_in, tip_count, intree);
```

This function:
1. Reads each input tree from the tree file
2. For each tree, extracts all bipartitions (splits) by traversing the tree
3. Stores each bipartition in the `pattern_array` hash table
4. Increments the frequency count for each observed split
5. Tracks branch lengths associated with each split

The bipartition extraction works by post-order traversal: at each internal
node, the set of descendant taxa defines one side of the split.

### 2. consensus() (in cons.c)

Called from `main()` (line 404):
```c
consensus(pattern_array, trees_in);
```

This function implements the core consensus algorithm:

1. **Sort splits by frequency**: All observed splits are sorted by how
   often they appear across the input trees.

2. **Filter splits by criterion**:
   - Strict: keep only if frequency == ntrees (100%)
   - Majority rule: keep if frequency > ntrees/2 (>50%)
   - Extended majority rule: keep all majority rule splits, then greedily
     add compatible splits in decreasing frequency order
   - Ml: keep if frequency >= mlfrac * ntrees

3. **Build consensus tree**: Starting with the trivial splits (each taxon
   vs. rest), add accepted splits one by one. Each split defines an internal
   edge in the consensus tree. The tree is built by:
   - For each accepted split, find where it fits in the growing tree
   - Create a new internal node that separates the taxa on each side
   - Store the frequency as the branch "length" (for annotation)

4. **Compatibility check** (for MRe): A new split is compatible with the
   existing tree if it is either a subset of or disjoint from every
   already-accepted split. This is checked via bitwise operations on the
   split representations.

### 3. treeout() - Newick Output (lines 276-350)

Writes the consensus tree in Newick format with split frequencies as
branch annotations:

```c
void treeout(node *p) {
    if (p->tip) {
        // print taxon name (spaces -> underscores)
    } else {
        putc('(', outtree);
        // recursively print children, comma-separated
        putc(')', outtree);
    }
    // print branch annotation (frequency)
    if (!strict) {
        fprintf(outtree, ":%5.1f", x);  // x = frequency count
    }
}
```

For strict consensus, no branch lengths are printed (they would all be
`ntrees`). For other methods, the frequency (number of trees containing
the split) is printed as a branch length.

### 4. count_siblings() (lines 250-273)

A traversal helper that counts the number of children at a node (walking
the ring structure). Uses a safety limit of 1000 iterations to prevent
infinite loops.

## Consensus Methods in Detail

### Strict Consensus
The most conservative method. A split appears in the strict consensus if
and only if it appears in every single input tree. The strict consensus
can be very poorly resolved (few internal edges) if the input trees
disagree on most splits.

### Majority Rule (MR)
A split is included if it appears in more than 50% of input trees. By the
"splits equivalence theorem," any set of splits each with >50% frequency
is guaranteed to be pairwise compatible, so the MR consensus is always
well-defined.

### Extended Majority Rule (MRe)
Starts with the MR consensus, then greedily adds compatible splits in
decreasing frequency order. This produces a fully resolved tree when
possible. It is the default consensus method in CONSENSE.

### Ml Consensus
Like MR but with a user-specified threshold `mlfrac` (between 0.5 and 1.0).
Setting `mlfrac = 0.5` gives the standard MR consensus. Setting
`mlfrac = 0.95` gives a very conservative consensus requiring near-unanimous
agreement.

## I/O

**Input**: A file of Newick trees (one per line, semicolon-terminated).
The number of trees is determined by counting semicolons (`countsemic`).
The number of taxa is determined by counting commas plus one (`countcomma`).
Trees can be rooted or unrooted (controlled by `noroot` flag).

**Output**:
- List of species (line 383: "Species in order:")
- Set membership for each split (via `printset()` in cons.c)
- ASCII tree diagram (via `printree()` in cons.c)
- Newick tree with frequency annotations (via `treeout()`)

## Program Architecture

The `main()` function (lines 353-443) orchestrates the workflow:

```c
int main() {
    openfiles();
    getoptions();           // interactive option setting
    trees_in = countsemic(&intree);
    tip_count = countcomma(&intree) + 1;

    // Phase 1: Read all trees, extract and hash bipartitions
    read_groups(&pattern_array, trees_in, tip_count, intree);

    // Phase 2: Compute consensus
    consensus(pattern_array, trees_in);

    // Phase 3: Output
    treeout(root);          // write Newick tree

    // Cleanup
    free all nodes and arrays
}
```

The split between `consense.c` and `cons.c` means the consensus algorithms
are reusable. Other PHYLIP programs can call `read_groups()` and
`consensus()` directly.

## Complexity

- **Time**: O(T * n * s) for reading T trees with n taxa, where s is the
  number of unique splits per tree (up to n-3 for binary unrooted trees).
  Hash table operations are O(1) amortized. Sorting splits is
  O(S * log(S)) where S is total unique splits. Compatibility checking
  for MRe is O(S^2 * w) where w is the number of words in the bit vector.
- **Space**: O(S * w) for the hash table of splits, plus O(T) for reading
  trees. With many trees and many taxa, the hash table can be large.
  Initial size is 32767 (`maxgrp`, line 377).

## Numerical Notes

1. **Split frequency as branch length**: The consensus tree's "branch
   lengths" are actually split frequencies (number of input trees
   containing that split). Values range from 1 to ntrees.

2. **Hash table sizing**: The initial hash table size of 32767 (line 377)
   is rehashed dynamically (via `rehash()` in cons.c) when load factor
   exceeds a threshold.

3. **Bit vector operations**: Split compatibility is checked via bitwise
   AND: two splits are compatible if their intersection is empty, or one
   is a subset of the other.

## Modernization Notes for Rust Reimplementation

1. **Proper hash map**: Replace the custom hash table (`pattern_array`)
   with Rust's `HashMap<BitVec, SplitInfo>` using a proper bitvector type
   from the `bitvec` crate.

2. **Split as a first-class type**: Create a `Split` struct with methods
   for intersection, union, complement, compatibility testing, and
   comparison. Implement `Hash` and `Eq` for use as HashMap keys.

3. **Consensus method as enum/trait**:
   ```rust
   enum ConsensusMethod {
       Strict,
       MajorityRule,
       ExtendedMajorityRule,
       Threshold(f64),
   }
   ```

4. **Streaming tree input**: Instead of reading all trees into memory,
   process trees one at a time, updating split counts incrementally.
   This reduces memory usage from O(T*n) to O(n^2) for the split table.

5. **Separate cons.c logic**: In Rust, the consensus algorithms should be
   a standalone library crate, not just a shared C file. This enables
   reuse across programs without the global state coupling.

6. **Tree builder**: Replace the procedural tree construction with a
   builder pattern:
   ```rust
   let tree = ConsensusTreeBuilder::new(taxa)
       .add_splits(accepted_splits)
       .build()?;
   ```

7. **Parallel tree reading**: When processing hundreds of bootstrap trees,
   tree parsing and split extraction can be parallelized across trees
   using `rayon`.

8. **Support for branch length averaging**: In addition to split
   frequencies, compute average branch lengths across trees that contain
   each split. This is standard in modern consensus tools but not in
   PHYLIP's CONSENSE.
