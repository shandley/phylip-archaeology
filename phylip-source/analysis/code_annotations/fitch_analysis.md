# FITCH Analysis - Fitch-Margoliash Distance Method

**Source file**: `phylip-3.698/src/fitch.c` (1203 lines)
**Dependencies**: `phylip.h`, `dist.h`

## Overview

FITCH implements the Fitch-Margoliash method (Fitch & Margoliash 1967) for
phylogenetic tree construction from a distance matrix. It finds the tree
topology and branch lengths that minimize a weighted least-squares criterion:

```
S = SUM_{i<j} w(i,j) * (d_obs(i,j) - d_tree(i,j))^2
```

where:
- `d_obs(i,j)` is the observed pairwise distance
- `d_tree(i,j)` is the tree-implied (patristic) distance
- `w(i,j) = 1 / d_obs(i,j)^power` is the weight

The program also supports **Minimum Evolution** (`minev`), where instead of
minimizing the weighted sum of squares, it minimizes the total sum of branch
lengths (while still fitting distances).

Key differences from NEIGHBOR:
- FITCH performs tree search (not just one-pass clustering)
- FITCH optimizes branch lengths iteratively
- FITCH can use local and global rearrangements
- FITCH is statistically more rigorous but much slower

## Key Data Structures

### Per-Node Distance and Weight Arrays
```c
node->d[]   // double[nonodes2] - distances from this node to all others
node->w[]   // double[nonodes2] - weights for distance fitting
node->v     // branch length of edge leading to parent
node->dist  // computed distance used in branch length estimation
node->iter  // boolean - whether this branch length is optimizable
```

Each node maintains its own view of distances (`d[]`) and weights (`w[]`)
to all other nodes in the tree. These views are updated as the tree topology
changes.

### Distance Matrix (lines 87-88)
```c
vector *x;       // x[spp][nonodes2] - input distances (tip-indexed)
intvector *reps;  // reps[spp][spp] - replicate counts
```

### Tree Copies (line 100)
```c
tree curtree, priortree, bestree, bestree2;
```
Four tree copies are maintained for the search:
- `curtree`: current working tree
- `priortree`: copy before last modification (for backtracking)
- `bestree`: best tree found so far
- `bestree2`: best tree across multiple jumble replicates

### Key Constants (lines 35-38)
```c
#define zsmoothings     10      // iterations for zero-branch correction
#define epsilonf        0.000001 // convergence threshold
#define delta           0.0001   // convergence for likelihood
```

### Key Global Variables
- `power` (line 92): the exponent in the weight function (default 2.0)
- `minev` (line 89): use Minimum Evolution criterion
- `global` (line 89): enable global rearrangements
- `negallowed` (line 89): allow negative branch lengths
- `addwhere` (line 96): best insertion point found during traversal

## Core Algorithms

### 1. evaluate() - Weighted Least-Squares Score (lines 429-440)

Computes the fit criterion by traversing the tree twice:

```c
double evaluate(tree *t) {
    double sum = 0.0;
    long nx = 0;
    firsttraverse(t->start->back, &nx, &sum);
    firsttraverse(t->start, &nx, &sum);
    t->likelihood = -sum;   // negate because search maximizes
    return (-sum);
}
```

**firsttraverse()** (lines 412-426): Visits all tips. For each tip in
non-minev mode, calls `secondtraverse()` to compute distances.

**secondtraverse()** (lines 394-409): From each tip, traverses the tree
accumulating the path length `z = y + q->v`, then at each other tip:
```c
TEMP = q->d[nx-1] - z;           // observed minus tree distance
sum += q->w[nx-1] * (TEMP * TEMP); // weighted squared difference
```

For Minimum Evolution (`minev`), `firsttraverse()` simply sums branch
lengths (line 416: `sum += p->v`).

The criterion is negated to `-sum` so the tree search framework
(which maximizes) works correctly.

### 2. nudists() - Update Node Distances (lines 443-477)

When the tree topology changes, distances between interior nodes and tips
must be recomputed. `nudists()` updates the distance and weight views for
a pair of nodes `x` and `y`:

```c
// Distance from x to y, combining info from x's two children
x->d[ny-1] = ((dil - vi) * wil + (djl - vj) * wjl) / (wil + wjl)
x->w[ny-1] = wil + wjl

// Distance from y to x (using y's view of x's children)
y->d[nx-1] = ((dil - vi) * wil + (djl - vj) * wjl) / (wil + wjl)
y->w[nx-1] = wil + wjl
```

This is a weighted average of the two child-to-y distances, adjusted by
subtracting the child branch lengths. The weight is the sum of the
children's weights.

### 3. makedists() - Three-Way Distance Computation (lines 480-501)

Computes `p->dist` for each of the three edges at an internal node, using
weighted averages of the observed distances:

```c
p->dist = (s->w[nr-1] * s->d[nr-1] + r->w[ns-1] * r->d[ns-1])
          / (s->w[nr-1] + r->w[ns-1])
```

where `r` and `s` are the nodes across from edge `p`. This gives the
"best estimate" of the total path length through this edge.

### 4. makebigv() - Initial Branch Length Estimation (lines 504-522)

Computes initial branch lengths from the three-way distances using the
standard formula:

```c
p->v = (p->dist + r->dist - q->dist) / 2.0
```

This is the same formula used in the three-point method: the branch
leading to `p->back` equals half the sum of the two distances that include
it minus the distance that doesn't.

### 5. correctv() - Iterative Branch Length Refinement (lines 525-560)

Refines branch lengths iteratively (10 rounds, `zsmoothings`), handling
the constraint that branches should be non-negative (unless `negallowed`):

```c
for (i = 1; i <= zsmoothings; i++) {
    for each of the 3 edges:
        wr = weights from one subtree
        wq = weights from other subtree
        p->v = ((p->dist - q->v) * wq + (r->dist - r->v) * wr) / (wr + wq)
        if (p->v < 0 && !negallowed)
            p->v = 0.0
}
```

This iteratively adjusts each branch length given the current values of
the other two, converging to a local optimum. The non-negativity constraint
is enforced by clamping.

### 6. alter() / nuview() - Subtree Distance Propagation (lines 563-590)

**alter()** recursively updates distances between a node and all nodes in
a subtree after a branch length change.

**nuview()** calls `alter()` from each of the three edges at an internal
node, propagating distance information through the entire tree.

### 7. update() - Full Node Update (lines 593-605)

Combines distance computation, branch length estimation, and view
propagation for a single internal node:

```c
void update(node *p) {
    makedists(p);
    if (any edge is optimizable)
        makebigv(p);
        correctv(p);
    nuview(p);
}
```

### 8. smooth() - Tree-Wide Optimization (lines 608-616)

Recursively visits all nodes, calling `update()` at each:

```c
void smooth(node *p) {
    if (p->tip) return;
    update(p);
    smooth(p->next->back);
    smooth(p->next->next->back);
}
```

### 9. insert_() - Insert Node and Optimize (lines 646-667)

Inserts a subtree and iteratively optimizes branch lengths:

```c
hookup(p->next->next, q->back);
hookup(p->next, q);
x = q->v / 2.0;               // split existing branch
// set initial branch lengths
fillin(p->back, p, contin_);   // compute distance views
evaluate(&curtree);
do {
    oldlike = curtree.likelihood;
    smooth(p);
    smooth(p->back);
    evaluate(&curtree);
} while (fabs(curtree.likelihood - oldlike) > delta);
```

### 10. addtraverse() - Tree Search (lines 783-800)

Tries inserting a subtree at position `q` and all positions below it:

```c
void addtraverse(node *p, node *q, boolean contin, ...) {
    insert_(p, q, true);
    if (evaluate(&curtree) > bestree.likelihood + epsilon)
        copy_(&curtree, &bestree);    // save better tree
    copy_(&priortree, &curtree);      // restore original
    if (!q->tip && contin) {
        addtraverse(p, q->next->back, contin, ...);
        addtraverse(p, q->next->next->back, contin, ...);
    }
}
```

### 11. rearrange() - Local SPR Rearrangement (lines 885-904)

For each internal node, removes one child and tries reinserting it at
neighboring edges:

```c
void rearrange(node *p, ...) {
    if (!p->tip && !p->back->tip) {
        r = p->next->next;
        re_move(&r, &q);
        addtraverse(r, q->next->back, false, ...);
        addtraverse(r, q->next->next->back, false, ...);
        copy_(&bestree, &curtree);
    }
    rearrange(p->next->back, ...);
    rearrange(p->next->next->back, ...);
}
```

### 12. globrearrange() - Global Rearrangement (lines 815-882)

More thorough search: for each internal node, removes each child subtree
and tries reinserting it at every edge in the tree:

```c
for each internal node i:
    for each sibling j of i:
        remove j
        for each possible insertion point k:
            addtraverse(j, k, true, ...)
        if found better tree: save to globtree
    restore original tree
copy globtree to curtree
```

This allocates temporary tree copies (`globtree`, `oldtree`) for the search,
then frees them afterward.

### 13. maketree() - Complete Search (lines 1020-1152)

**For heuristic search** (lines 1062-1145):
1. Build a 3-taxon starting tree (`buildsimpletree`)
2. Stepwise addition: add each remaining taxon at the best position
3. After each addition, apply local rearrangements
4. When all taxa added and `global`, apply `globrearrange()`
5. Repeat for `njumble` random addition orders
6. Keep best tree across jumbles in `bestree2`

**For user trees** (lines 1027-1061):
1. Read each user tree from input file
2. Set up tip distances (`setuptipf`)
3. Iterate branch lengths (`treevaluate`)
4. Report fit statistics

## Weight Function and Power Parameter

The weight function `w(i,j) = n(i,j) / d(i,j)^power` (line 744) controls
how much each distance contributes to the criterion:

- `power = 0`: Ordinary least squares (all distances weighted equally)
- `power = 1`: Cavalli-Sforza & Edwards criterion
- `power = 2`: Fitch-Margoliash criterion (default) -- gives more weight
  to smaller, more accurately estimated distances

The weights are computed in `setuptipf()` (lines 730-757) and stored in
the `w[]` array for each tip.

## Minimum Evolution Option

When `minev` is true (line 89), the criterion changes from weighted least
squares to minimizing the total tree length (sum of all branch lengths).
In `firsttraverse()` (line 416), only `sum += p->v` is accumulated instead
of the full WLS calculation. This is Rzhetsky & Nei's (1992) minimum
evolution principle.

## I/O

**Input**: PHYLIP distance matrix format (square or triangular). Same
format as NEIGHBOR. Supports replicate counts.

**Output**:
- Sum of squares or sum of branch lengths (line 938/940)
- Average percent standard deviation (for power=2, line 950-951)
- Branch length table (lines 953-957)
- Newick tree with branch lengths scaled by log10(e) = 0.43429 (line 961)

## Numerical Methods

1. **Iterative branch length optimization**: The `smooth()` function
   iterates until the likelihood change is less than `delta = 0.0001`
   (line 666). This is a coordinate descent approach.

2. **Zero-branch correction**: `correctv()` runs 10 iterations
   (`zsmoothings`) to handle the non-negativity constraint, using a
   form of projected gradient descent.

3. **Distance propagation**: The `nudists()`/`alter()` system propagates
   distance information through the tree after topological changes,
   maintaining O(n) distance views per node.

## Complexity

- **Time**: O(n^3) per tree evaluation (each of n tips must compute
  distances to all other tips via O(n) traversal). Stepwise addition
  is O(n^2) insertions, each with O(n^3) evaluation, giving O(n^5)
  for construction. Global rearrangement adds O(n^2) SPR moves per
  round, each O(n^3), so O(n^5) per round.
- **Space**: O(n^2) for the distance arrays at all nodes (each node
  stores distances to all n nodes).

## Modernization Notes for Rust Reimplementation

1. **Eliminate per-node distance arrays**: Each node storing O(n) distances
   is wasteful. Consider computing distances on-the-fly from branch lengths
   during evaluation, trading time for space.

2. **Trait-based criterion**: Abstract the WLS and MinEv criteria via a
   trait:
   ```rust
   trait DistanceCriterion {
       fn evaluate(&self, tree: &Tree, distances: &DistMatrix) -> f64;
   }
   ```

3. **Tree copy optimization**: The repeated `copy_()` calls (lines 791-795,
   843-864) are expensive. Use an undo-stack approach instead of full tree
   copies.

4. **Power parameter**: Move from a global `power` variable to a
   configuration struct. Consider supporting arbitrary weight functions
   via closures.

5. **Branch length optimizer**: Replace the ad-hoc `makebigv`/`correctv`
   cycle with a proper constrained optimization library call (e.g.,
   L-BFGS-B for non-negative constraints).

6. **Parallel evaluation**: The `firsttraverse`/`secondtraverse` pair
   computes independent pairwise contributions. These can be parallelized
   across tip pairs.

7. **Global rearrangement memory**: The current implementation allocates
   and frees full tree copies within `globrearrange()` (lines 825-881).
   Use a persistent tree pool instead.

8. **Remove redundant tree structure**: The `curtree`/`priortree`/`bestree`/
   `bestree2` pattern appears in FITCH, DNAML, and other programs. Abstract
   this into a `TreeSearchState` struct.
