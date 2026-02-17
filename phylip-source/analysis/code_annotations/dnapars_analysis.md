# DNAPARS Analysis - DNA Parsimony (Fitch Algorithm)

**Source file**: `phylip-3.698/src/dnapars.c` (1663 lines)
**Dependencies**: `phylip.h`, `seq.h`

## Overview

DNAPARS implements maximum parsimony phylogenetic inference for DNA sequences.
It finds the tree(s) requiring the minimum total number of nucleotide
substitutions (steps) to explain the observed sequences. The core algorithm is
the Fitch parsimony algorithm (Fitch 1971), which computes the minimum number
of changes at each site independently using a post-order tree traversal.

The program performs:
1. Tree search via stepwise addition + local/global rearrangements
2. Step counting via the Fitch algorithm (set intersection/union)
3. Storage and enumeration of multiple equally parsimonious trees
4. Support for threshold parsimony and transversion-only parsimony
5. Ancestral state reconstruction via the Fitch method

Unlike DNAML which optimizes a continuous likelihood function, DNAPARS
optimizes a discrete step count, making the algorithm fundamentally different
in its optimization landscape.

## Key Data Structures

### Per-Node Arrays (from phylip.h, allocated in allocnode)
```c
node->base[]       // long[endsite] - bitwise-encoded base sets
node->numsteps[]   // long[endsite] - step counts per site
node->numnuc[]     // nucarray[endsite] - per-base counts at each site
node->numdesc      // number of descendants (for n-ary trees)
node->sumsteps     // accumulated weighted step count
```

The `base[]` array uses bit encoding for nucleotide sets: each site stores
the set of possible ancestral states as a bitmask. The Fitch algorithm works
by intersecting sets going up the tree and taking unions when intersections
are empty.

### bestrees Array (line 97)
```c
bestelm *bestrees;  // bestrees[maxtrees] - stores topology encodings
```
Each `bestelm` contains a `btree[]` array encoding a tree topology as a
sequence of attachment points, plus flags for whether global rearrangement
has been performed. This allows storing thousands of equally parsimonious
trees compactly.

### Temporary Node Arrays (lines 101-103)
```c
node *temp, *temp1, *temp2, *tempsum, *temprm, *tempadd, *tempf;
node *tmp, *tmp1, *tmp2, *tmp3, *tmprm, *tmpadd;
```
Thirteen temporary node structures used during tree rearrangement to avoid
modifying the actual tree until a better arrangement is confirmed. This
large number reflects the complexity of the n-ary parsimony SPR operations.

### Key Global Variables (lines 74-105)
- `root`: pointer to the current root node
- `treenode[]`: array of all node pointers (tips + internal)
- `enterorder[]`: randomized species addition order
- `threshwt[]`: threshold weights per site (for threshold parsimony)
- `like`, `bestyet`, `bestlike`, `bstlike2`: step count tracking
- `there`: best insertion/rearrangement point found
- `mulf`: flag indicating multifurcation insertion
- `lastrearr`: flag for final rearrangement pass
- `grbg`: garbage list for recycling freed nodes
- `zeros[]`: zero-filled array for initializations

## Core Algorithms

### 1. evaluate() - Parsimony Score (lines 506-537)

Computes the total weighted parsimony score (number of steps) for the tree.

```c
for each site i:
    steps = root->numsteps[i]
    if steps <= threshwt[i]:
        term = steps
    else:
        term = threshwt[i]    // cap at threshold
    sum += term
like = -sum    // negative because search maximizes "likelihood"
```

The score is negated because the search framework maximizes `like`, but
parsimony minimizes steps. The threshold capping (when `thresh` is true)
implements threshold parsimony, where additional steps beyond a threshold
do not count.

For user trees, per-site step counts are stored in `fsteps[][]` and compared
across trees using `standev()` to compute Templeton test statistics.

### 2. Fitch Algorithm (via sumnsteps/multifillin in seq.c)

The core parsimony computation uses the Fitch algorithm, implemented through
helper functions in `seq.c`:

- `multifillin()`: Computes the Fitch intersection/union at an internal
  node given its children's base sets.
- `sumnsteps()`: Combines two child nodes' base/step data.
- `sumnsteps2()`: Same but with threshold-weighted step accumulation.

For each site, the algorithm works as follows:
1. At each tip: `base[site]` = the observed nucleotide (as a bit)
2. At each internal node:
   - Intersection = base[child1] AND base[child2]
   - If intersection is non-empty: base[node] = intersection, steps unchanged
   - If intersection is empty: base[node] = base[child1] OR base[child2],
     steps += 1

### 3. tryadd() - Evaluate Insertion Point (lines 540-639)

This function evaluates the parsimony score when inserting a new species
(`item`) at a given position (`p`) in the tree using fork node `nufork`.

The function handles two cases:
1. **Below insertion** (standard): Insert between `p` and `p->back`, creating
   a new internal node. Computed via `sumnsteps()` and `sumnsteps2()`.
2. **Multifurcation insertion** (`multf`): Add as an additional child of
   an existing internal node (when `belowsum >= p->sumsteps`).

During the last rearrangement pass (`lastrearr`), it also manages the
`bestrees[]` array, storing equally parsimonious trees via `addbestever()`
and `addtiedtree()`. Collapsible branches (zero-length) are detected and
handled to avoid storing equivalent trees.

### 4. addpreorder() - Traverse for Best Insertion (lines 642-658)

Pre-order traversal of the tree, calling `tryadd()` at each node to find
the best insertion point for a new species:

```c
void addpreorder(node *p, node *item, node *nufork) {
    tryadd(p, item, nufork);
    if (!p->tip) {
        q = p->next;
        while (q != p) {
            addpreorder(q->back, item, nufork);
            q = q->next;
        }
    }
}
```

### 5. trylocal() / trylocal2() - Local Rearrangement (lines 836-1050)

These functions implement SPR (Subtree Pruning and Regrafting) moves for
local rearrangement:

**trylocal()** (lines 836-938): Used when the fork node has more than 2
descendants. It:
1. Removes the item from its fork node (reducing its degree)
2. Tries reinserting at the fork node itself, at `forknode->back`, and at
   all descendants of the fork node's siblings
3. Uses `trydescendants()` for recursive exploration

**trylocal2()** (lines 941-1050): Used when the fork node has exactly 2
descendants (binary node). It:
1. Removes the item, collapsing the binary fork
2. Tries reinserting along the "other" child's subtree

Both functions maintain careful bookkeeping of `temprm` (removed state)
and `tempadd` (added state) to compute parsimony scores without actually
modifying the tree.

### 6. trydescendants() - Explore Rearrangement Neighborhood (lines 661-833)

Recursively explores potential rearrangement positions below a given node.
For each position, it computes the parsimony score and compares against
`bstlike2` (current best). It handles both "above" insertions (as a new
sibling) and "below" insertions (splitting an edge).

This is the most complex function in the file, reflecting the difficulty of
efficient SPR moves in the context of potentially multifurcating trees.

### 7. rearrange() - Iterative Local Search (lines 1128-1141)

Repeatedly applies local rearrangements until no improvement is found:

```c
while (success) {
    success = false;
    clearvisited(treenode);
    repreorder(root, &success);
}
```

Uses a `visited` flag on each node to avoid redundant rearrangement
attempts within a single pass.

### 8. globrearrange() - Global Rearrangement (lines 1231-1287)

More thorough SPR search: for every node in the tree, removes it and tries
reinserting it at every possible position:

```c
for each node j in tree:
    remove j
    addpreorder(root, j, nufork)   // try all positions
    reinsert at best position
```

Repeats until no further improvement is found (`bestlike > gotlike`).

### 9. grandrearr() - Grand Rearrangement on Stored Trees (lines 1317-1337)

Applies global rearrangement to each stored best tree:

```c
do {
    find next unrearranged tree in bestrees
    load_tree(treei)
    globrearrange()
} while (!done);
```

This ensures all equally parsimonious trees have been globally rearranged.

### 10. maketree() - Main Tree Search (lines 1340-1533)

The complete search procedure:

1. **Stepwise addition** (lines 1358-1418): Add species one at a time,
   trying all positions via `addpreorder()`, with local rearrangement
   after each addition.

2. **Global rearrangement** (lines 1416-1423): After all species are added,
   apply `globrearrange()` and `grandrearr()`.

3. **Output** (lines 1429-1469): Load each stored tree, reroot at outgroup,
   compute parsimony score, print tree and descriptions.

4. **User trees** (lines 1471-1521): If user trees are provided, evaluate
   each, compute step counts, and perform Templeton's test for comparing
   trees.

## Parsimony Variants

### Threshold Parsimony (lines 412-418)
When `thresh` is true, step contributions per site are capped at `threshold`:
```c
threshwt[i] = (long)(threshold * weight[i] + 0.5);
```
Steps beyond the threshold are not counted, giving less weight to highly
variable (potentially saturated) sites.

### Transversion Parsimony (line 121, `transvp`)
When `transvp` is true, only transversions (purine <-> pyrimidine changes)
are counted. Transitions (A<->G, C<->T) are ignored. This is implemented
by modifying the base encoding: purines and pyrimidines are grouped together.

## N-ary Trees

Unlike many parsimony programs that work only with binary trees, DNAPARS
supports n-ary (multifurcating) trees. This is important because:
1. Multiple equally parsimonious binary resolutions of a node may exist
2. Collapsing zero-length branches produces multifurcations
3. The `numdesc` field tracks the number of descendants at each node

The `multifillin()`, `collapsetree()`, and `collapsebestrees()` functions
handle the bookkeeping for multifurcating trees.

## Tree Storage

Equally parsimonious trees are stored compactly in the `bestrees[]` array.
Each tree is encoded as a sequence of attachment points:
```c
bestrees[i].btree[j]  // where species j was attached
```
Positive values indicate attachment to a node; negative values indicate
attachment as a sibling (multifurcation). This encoding allows storing
up to `maxtrees` (default 10000) trees without full tree copies.

## I/O

**Input**: PHYLIP format DNA sequences (interleaved or sequential).
Site patterns are compressed via `makeweights()` (lines 392-422):
identical site patterns are combined with accumulated weights.

**Output**:
- Parsimony score: "requires a total of X.XXX" (line 1150)
- Branch lengths (measured in steps): `printbranchlengths()` (line 1153)
- Per-site step counts (optional): `writesteps()` (line 1156)
- Ancestral state reconstruction (optional): `hypstates()` (line 1158)
- Newick trees: `treeout3()` (line 1164)

## Complexity

- **Time**: O(n^2 * s) per tree evaluation (n taxa, s unique sites).
  Stepwise addition is O(n^2) insertions, each requiring O(n*s) for
  the Fitch algorithm, giving O(n^3 * s) for construction.
  Global rearrangement adds O(n^2) SPR moves per round, each O(n*s),
  so O(n^3 * s) per round. Multiple rounds may be needed.
- **Space**: O(n * s) for the base/step arrays at all nodes, plus
  O(maxtrees * n) for stored tree topologies.

## Modernization Notes for Rust Reimplementation

1. **Eliminate temporary nodes**: The 13 temporary node variables (temp,
   temp1, ..., tmpadd) should be replaced with stack-allocated buffers
   or a small pool. In Rust, these could be local variables or a
   `TempNodePool` struct.

2. **Bit-packed base sets**: The `base[]` array already uses bitwise
   operations. In Rust, use explicit bitflags for clarity:
   ```rust
   bitflags! { struct BaseSet: u8 { const A=1; const C=2; const G=4; const T=8; } }
   ```

3. **Tree topology encoding**: Replace the implicit `btree[]` encoding
   with an explicit enum-based topology representation. The positive/negative
   convention is brittle.

4. **Separate concerns**: The monolithic `tryadd()` function (100 lines)
   mixes score computation, tree comparison, and tree storage. These should
   be separate functions.

5. **Iterator-based traversal**: Replace the recursive `addpreorder()` with
   an iterator that yields insertion candidates, separating traversal from
   evaluation.

6. **Parallel site evaluation**: The Fitch algorithm evaluates each site
   independently. Use `rayon` for parallel iteration over sites in
   `evaluate()` and `sumnsteps()`.

7. **SIMD for Fitch operations**: The bitwise AND/OR operations on base
   sets across sites are naturally parallelizable with SIMD. Process
   multiple sites simultaneously using wide registers.

8. **Proper error types**: Replace the implicit `-10 * spp * chars`
   sentinel for "worst possible score" (lines 1069, 1244, 1371) with
   `Option<Score>` or a dedicated `Score` type.

9. **Tree collapse as post-processing**: The inline collapse detection
   during rearrangement (lines 610-614) complicates the code. Consider
   separating tree search (on binary trees) from post-search collapse.

10. **Generic parsimony framework**: The Fitch algorithm is not DNA-specific.
    Abstract it to work with any character type via a trait:
    ```rust
    trait ParsimonyCharacter {
        fn intersection(&self, other: &Self) -> Option<Self>;
        fn union(&self, other: &Self) -> Self;
    }
    ```
