# NEIGHBOR Analysis - Neighbor-Joining and UPGMA Clustering

**Source file**: `phylip-3.698/src/neighbor.c` (629 lines)
**Dependencies**: `phylip.h`, `dist.h`

## Overview

NEIGHBOR implements two distance-based tree construction algorithms:
1. **Neighbor-Joining (NJ)** -- Saitou & Nei (1987), an agglomerative clustering
   algorithm that produces an unrooted tree from a distance matrix. NJ is the
   most widely used distance-based method in phylogenetics due to its speed and
   statistical consistency (it converges to the true tree as sequence length
   grows, given a correct distance measure).
2. **UPGMA** -- Unweighted Pair Group Method with Arithmetic Mean (Sokal &
   Michener, 1958), producing a rooted ultrametric tree. UPGMA assumes a
   molecular clock (equal rates across all lineages).

Both algorithms are greedy: they iteratively join the pair of nodes that
minimizes a criterion, never revisiting previous decisions. This makes them
O(n^3) time and O(n^2) space, where n = number of taxa.

The program takes a precomputed distance matrix as input (e.g., from DNADIST)
and outputs a tree in Newick format.

## Key Data Structures

### Distance Matrix (lines 56, 213-215)
```c
vector *x;       // x[spp][spp] -- the distance matrix (modified in place)
intvector *reps;  // reps[spp][spp] -- replicate counts for each distance
```
The distance matrix `x` is a full `spp x spp` matrix of doubles. It is
symmetrized at the start of `jointree()` (lines 353-359) by averaging
`x[i][j]` and `x[j][i]`. During the algorithm, distances to newly formed
clusters replace the row/column of the "absorbed" node.

### Cluster Array (line 66)
```c
node **cluster;   // cluster[spp] -- active cluster representatives
```
Each entry points to the current representative node for that cluster. When
two clusters are joined, one entry is updated to the new internal node while
the other is set to NULL, marking it as inactive.

### Tree Structure (line 60)
```c
tree curtree;     // the tree being constructed
```
Uses the standard PHYLIP tree structure with `nodep` array, `start` pointer,
and ring-structure internal nodes (each internal node has 3 nodes linked via
`next` pointers forming a ring).

### Other Key Variables
- `njoin` (line 59): boolean flag, true for NJ, false for UPGMA
- `enterorder` (line 62): permutation array for input order randomization
- `R[]` (line 350): NJ sum-of-distances array (added by Y. Ina revision)
- `av[]` (line 363): accumulated branch lengths for each active cluster
- `oc[]` (line 364): cluster sizes (number of OTUs in each cluster)
- `fotu2` (line 361): remaining OTUs minus 2, used in NJ criterion

## Core Algorithms

### 1. Neighbor-Joining Algorithm (jointree, lines 341-536)

The NJ algorithm works by iteratively finding the pair (i,j) that minimizes
the modified distance criterion Q(i,j), then joining them into a new node.

**Initialization** (lines 353-368):
1. Symmetrize the distance matrix: `x[i][j] = (x[i][j] + x[j][i]) / 2`
2. Set `fotu2 = spp - 2` (number of active OTUs minus 2)
3. Initialize `av[i] = 0`, `oc[i] = 1` for all taxa
4. Set iteration count: `iter = spp - 3` for NJ, `spp - 1` for UPGMA

**Main cycle** (lines 374-487, repeated `iter` times):

Step 1 -- Compute R[i] sums (lines 381-396, NJ only):
```
R[i] = SUM over all active j != i of x[i][j]
```
This is the sum of distances from taxon i to all other active taxa.
(Revision by Y. Ina for computational efficiency.)

Step 2 -- Find minimum pair (lines 397-416):
For NJ, the criterion is:
```
Q(i,j) = fotu2 * x[i][j] - R[i] - R[j]
```
For UPGMA, the criterion is simply `x[i][j]`.
The pair (mini, minj) minimizing this criterion is selected.

Step 3 -- Compute branch lengths (lines 418-436):
For NJ (lines 418-431):
```
dio = (SUM_k x[i][k] - x[i][j]) / fotu2
djo = (SUM_k x[j][k] - x[i][j]) / fotu2
bi = (x[i][j] + dio - djo) / 2 - av[i]
bj = x[i][j] - bi - av[j]    (derived: bj = x[i][j] - (bi + av[i]) ... simplified)
```
For UPGMA (lines 432-436):
```
bi = x[i][j] / 2 - av[i]
bj = x[i][j] / 2 - av[j]
```
This assigns half the distance to each branch, maintaining the ultrametric
property.

Step 4 -- Join clusters (lines 453-461):
```c
hookup(curtree.nodep[nextnode-1]->next, cluster[mini-1]);
hookup(curtree.nodep[nextnode-1]->next->next, cluster[minj-1]);
cluster[mini-1]->v = bi;
cluster[minj-1]->v = bj;
cluster[mini-1] = curtree.nodep[nextnode-1];  // new cluster replaces mini
cluster[minj-1] = NULL;                        // minj is absorbed
```

Step 5 -- Update distance matrix (lines 466-486):
For NJ (line 469):
```
x[new][k] = (x[mini][k] + x[minj][k]) / 2
```
For UPGMA (lines 475-476):
```
x[new][k] = (x[mini][k] * oc[mini] + x[minj][k] * oc[minj]) / (oc[mini] + oc[minj])
```
UPGMA uses a weighted average proportional to cluster sizes, ensuring the
distance from the new cluster to any other cluster is the average of all
pairwise distances between their constituent OTUs.

Step 6 -- Update bookkeeping:
```
fotu2 -= 1
oc[mini] += oc[minj]
av[mini] = dmin / 2  (NJ only)
```

**Last cycle** (lines 488-536):
When 3 clusters remain (NJ) or 1 (UPGMA):

For NJ: The final 3 clusters are connected to a single internal node
(the trifurcation root). Branch lengths are computed from the triangle
of distances:
```
bi = (x[0][1] + x[0][2] - x[1][2]) / 2 - av[0]
bj = x[0][1] - bi - av[1]    (simplified from the code)
bk = x[0][2] - bi - av[2]    (simplified)
```

For UPGMA: The last remaining cluster becomes the root with `back = NULL`.

### 2. maketree() - Tree Construction Entry Point (lines 539-578)

1. Read input distance matrix via `inputdata()` (from dist.c)
2. Optionally randomize species order (`jumble`)
3. Initialize clusters: `cluster[i] = curtree.nodep[i]` (each tip is its own cluster)
4. Call `jointree()` to construct the tree
5. For NJ: reroot at outgroup species
6. Print tree and branch lengths
7. Write Newick tree to file

### 3. describe() - Output Branch Information (lines 274-304)

Recursively prints branch length information. For NJ, prints "Between/And/Length"
format. For UPGMA, additionally prints cumulative height from root.

## Distance Updating: NJ vs UPGMA

The key algorithmic difference between NJ and UPGMA lies in the distance
update formula and the selection criterion:

| Aspect | NJ | UPGMA |
|--------|----|----|
| Criterion | Q(i,j) = (n-2)*d(i,j) - R(i) - R(j) | d(i,j) |
| Branch lengths | Asymmetric (different for each child) | Symmetric (d/2 each) |
| Distance update | Simple average | Weighted average by cluster size |
| Result | Unrooted tree | Rooted ultrametric tree |
| Assumption | None (consistent estimator) | Molecular clock |

## I/O

**Input**: PHYLIP distance matrix format (square or triangular). Read via
`inputdata()` from `dist.c`. Supports:
- Full square matrix
- Lower-triangular matrix (`lower` flag)
- Upper-triangular matrix (`upper` flag)
- Replicate counts (`replicates` flag)

**Output**: Text report of branch lengths (`summarize()`, lines 307-329),
plus Newick tree via `treeout()` (line 568). The Newick output uses the
constant 0.43429448222 as a scaling factor (the reciprocal of ln(10), for
converting natural log branch lengths to log10).

## Randomization

Input order randomization (`jumble`, lines 555-556) uses `randumize()` to
create a random permutation of `enterorder[]`. This affects which pair is
selected when ties occur in the NJ/UPGMA criterion. The NJ algorithm is
deterministic given a distance matrix, but ties can lead to different trees.

## Complexity

- **Time**: O(n^3) where n = number of taxa.
  - The main loop runs n-3 (NJ) or n-1 (UPGMA) iterations.
  - Each iteration scans O(n^2) pairs and updates O(n) distances.
  - The R[i] computation adds O(n^2) per iteration.
  - Total: O(n * n^2) = O(n^3).
- **Space**: O(n^2) for the distance matrix.
- Note: Studier & Keppler (1988) showed NJ can be implemented in O(n^2) time,
  but PHYLIP uses the straightforward O(n^3) implementation.

## Numerical Notes

1. **Negative branch lengths**: NJ can produce negative branch lengths
   (line 270: "Negative branch lengths allowed"). These are biologically
   meaningless but statistically valid -- they indicate the data is not
   perfectly tree-like.

2. **Distance symmetrization**: The matrix is explicitly symmetrized at the
   start (lines 353-359), handling any asymmetry in the input.

3. **DBL_MAX sentinel**: `tmin` is initialized to `DBL_MAX` (line 379) rather
   than an arbitrary large number, using the proper floating-point maximum.

## Modernization Notes for Rust Reimplementation

1. **Separate NJ and UPGMA**: The single `jointree()` function mixes NJ and
   UPGMA logic via boolean branches. These should be separate implementations
   sharing a common trait (e.g., `DistanceTreeBuilder`) in Rust.

2. **Distance matrix as owned type**: The in-place modification of the
   distance matrix is error-prone. Use a dedicated `DistanceMatrix` struct
   that manages its own storage and provides safe update operations.

3. **Eliminate global state**: The ~15 global variables (x, reps, cluster,
   curtree, enterorder, etc.) should be encapsulated in a `NeighborState`
   struct.

4. **Active cluster tracking**: Instead of using NULL entries in `cluster[]`
   and scanning for active entries, use a `Vec` or `HashSet` of active
   cluster indices for O(1) membership testing.

5. **O(n^2) NJ**: Implement the Studier-Keppler O(n^2) variant, which
   maintains R[i] incrementally rather than recomputing from scratch each
   iteration.

6. **Result type for tree**: Return a proper tree structure rather than
   modifying a global `curtree`. The tree should own its nodes and provide
   iteration over branches.

7. **Error handling**: Replace `exxit(-1)` calls with Rust's `Result` type
   for proper error propagation (e.g., "must have at least 3 species" on
   line 546).

8. **Branch length type**: Use a newtype wrapper `BranchLength(f64)` that
   can optionally enforce non-negativity for UPGMA while allowing negative
   values for NJ.
