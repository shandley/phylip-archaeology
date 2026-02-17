# DNAML Analysis - DNA Maximum Likelihood (Felsenstein's Pruning Algorithm)

**Source file**: `phylip-3.698/src/dnaml.c` (2619 lines)
**Dependencies**: `phylip.h`, `seq.h`

## Overview

DNAML implements DNA sequence maximum likelihood phylogenetic inference using
Felsenstein's pruning algorithm (1981). This is the single most important
algorithm in computational phylogenetics. It computes the likelihood of
observing a set of DNA sequences given a tree topology and a substitution
model (F84, which subsumes HKY85, K2P, and JC69).

The program performs:
1. Tree search via stepwise addition + local/global rearrangements (SPR)
2. Branch length optimization via Newton-Raphson
3. Likelihood evaluation via the post-order pruning algorithm
4. Support for rate heterogeneity (gamma-distributed rates, invariant sites, HMM)

## Key Data Structures

### valrec (line 35-38)
```c
typedef struct valrec {
  double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1, xiy1xv;
  double *ww, *zz, *wwzz, *vvzz;
} valrec;
```
Precomputed transition probability components for each combination of rate
category (rcategs) and site category (categs). The `ww`/`zz`/`wwzz`/`vvzz`
arrays are indexed by sibling, storing per-child transition probability
components during the nuview computation.

- `rat`: combined rate for this rate-category x site-category
- `ratxi`: rat * xi (transition rate parameter)
- `ratxv`: rat * xv (transversion rate parameter)
- `orig_zz`, `z1`, `z1zz`, `z1yy`: precomputed exponentials of branch length

### node (from phylip.h, lines 467-534)
```c
typedef struct node {
  struct node *next, *back;    // ring structure (next) and connection (back)
  long index;                   // 1-based node index
  boolean tip;                  // true if leaf node
  boolean iter;                 // whether branch length is optimizable
  boolean initialized;          // whether conditional likelihoods are current
  phenotype x;                  // x[site][ratecateg][base] - conditional likelihoods
  double v;                     // branch length
  double *underflows;           // log-underflow corrections per site
} node;
```

### phenotype / sitelike / ratelike (from phylip.h)
```c
typedef double sitelike[(long)T - (long)A + 1];  // 4-element array [A,C,G,T]
typedef sitelike *ratelike;                        // array over rate categories
typedef ratelike *phenotype;                       // array over sites
```
So `node->x[site][ratecateg][base]` gives the conditional likelihood of
observing the data at the subtree below this node, for a given site, rate
category, and ancestral base.

### tree (from phylip.h, lines 541-570)
```c
typedef struct tree {
  pointarray nodep;  // array of node pointers, index 0..nonodes2-1
  node *start;       // starting node for traversal
  node *root;        // root node (for rooted trees)
  double likelihood;
} tree;
```

### Global State (lines 99-134)
Massive amount of global state including:
- `curtree`, `bestree`, `bestree2`, `priortree` - tree copies for search
- `xi`, `xv` - F84 model parameters (transition/transversion eigenvalue components)
- `freqa`, `freqc`, `freqg`, `freqt` - base frequencies
- `ttratio` - transition/transversion ratio
- `rcategs`, `categs` - number of rate/site categories
- `probcat`, `rrate` - rate category probabilities and rates
- `alpha`, `cv`, `invarfrac` - gamma/invariant sites parameters
- `tbl[rcategs][categs]` - precomputed valrec lookup table
- `term`, `slopeterm`, `curveterm` - per-site likelihood derivative terms
- `contribution` - per-site per-category likelihood contributions
- `lambda` - autocorrelation parameter for HMM of rates
- `fracchange` - expected fraction of changed sites per unit branch length

## Core Algorithm: Felsenstein's Pruning Algorithm

### 1. nuview() - The Pruning/Peeling Step (lines 957-1096)

This is THE core function. It computes conditional likelihood vectors at an
interior node by combining information from its children using the pruning
algorithm.

**Algorithm:**
For each site i, rate category j:
1. For each child c with branch length v_c, compute transition probabilities:
   - `ww = exp(ratxi * (-v_c))` -- transition component
   - `zz = exp(ratxv * (-v_c))` -- transversion component
   - `wwzz = ww * zz`
   - `vvzz = (1 - ww) * zz`

2. For each child, compute intermediate sums:
   - `sum[c] = yy[c] * (freqa*x_c[A] + freqc*x_c[C] + freqg*x_c[G] + freqt*x_c[T])`
   - `sumr[c] = freqar*x_c[A] + freqgr*x_c[G]` (purine sum)
   - `sumy[c] = freqcy*x_c[C] + freqty*x_c[T]` (pyrimidine sum)
   - `vzsumr[c] = vvzz[c] * sumr[c]`
   - `vzsumy[c] = vvzz[c] * sumy[c]`

3. Multiply across children for each base b:
   ```
   p_xx[A] = PRODUCT over children of (sum[c] + wwzz[c]*x_c[A] + vzsumr[c])
   p_xx[C] = PRODUCT over children of (sum[c] + wwzz[c]*x_c[C] + vzsumy[c])
   p_xx[G] = PRODUCT over children of (sum[c] + wwzz[c]*x_c[G] + vzsumr[c])
   p_xx[T] = PRODUCT over children of (sum[c] + wwzz[c]*x_c[T] + vzsumy[c])
   ```

4. Store result in `p->x[i][j]`.

**Underflow handling** (lines 1087-1091): If the maximum conditional likelihood
value drops below `MIN_DOUBLE` (10e-100), `fix_x()` rescales and accumulates
the correction in `p->underflows[i]`.

This is the F84 substitution model's transition probability matrix applied
in the Felsenstein factored form. The three components correspond to:
- Same base (wwzz term)
- Same purine/pyrimidine class transition (vvzz term)
- Any base (yy term = 1-zz, the transversion component)

### 2. evaluate() - Likelihood Computation (lines 825-923)

Computes the log-likelihood of the entire tree by evaluating across the
edge connecting nodes p and p->back.

**Algorithm:**
1. Compute transition probabilities for the connecting edge (lines 840-846)
2. For each site i (lines 847-884):
   - For each rate category j, compute:
     ```
     tterm[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2
     ```
     where:
     - `prod12 = SUM_b(freq_b * x1[b] * x2[b])` -- same base
     - `prod3 = (purines1)(purines2) + (pyrimidines1)(pyrimidines2)` -- same class
     - `prod1 = SUM_b(freq_b * x1[b])`, `prod2 = SUM_b(freq_b * x2[b])` -- marginals
   - Sum over rate categories weighted by `probcat[j]`
   - Take log and add underflow corrections
   - Accumulate weighted by `aliasweight[i]`

3. HMM autocorrelation (lines 886-907): If `auto_`, apply Hidden Markov Model
   transition between rate categories across sites using parameter `lambda`.

### 3. slopecurv() - Likelihood, Slope, and Curvature (lines 1099-1236)

Computes the log-likelihood AND its first and second derivatives with respect
to branch length, using the chain rule through the same site-by-site
computation as evaluate(). The derivatives of the exponentials are:
- `zzs = -rat * zz` (first derivative of zz)
- `z1s = -ratxv * z1` (first derivative of z1)
- `zzc = rat^2 * zz` (second derivative of zz)
- `z1c = ratxv^2 * z1` (second derivative of z1)

These are needed for the Newton-Raphson branch length optimization.

### 4. makenewv() - Newton-Raphson Branch Optimization (lines 1239-1286)

Optimizes a single branch length using modified Newton-Raphson:
```
y_new = y_old + slope / |curve|
```
The division by `|curve|` (not `curve`) forces uphill movement. If the new
point is worse, it retracts 95% of the way back (line 1277). Convergence
when `|y - yold| < 0.1 * epsilon`.

### 5. smooth() / update() - Full Tree Optimization (lines 1289-1337)

`smooth()` recursively visits all nodes, calling `update()` which calls
`makenewv()` on each branch. Multiple passes through `smooth()` are
performed until convergence (controlled by `smoothings = 4`).

### 6. Tree Search: addtraverse() / rearrange() / globrearrange()

**Stepwise addition** (maketree, lines 2401-2468): Species are added one at a
time. For each new species, `addtraverse()` tries inserting it at every edge
in the tree, keeping the best insertion point.

**Local rearrangements** (rearrange, lines 1587-1663): SPR moves where a
subtree is pruned and regrafted nearby. For each internal node, remove one
child subtree and try reinserting it at neighboring edges.

**Global rearrangements** (globrearrange, lines 1504-1584): More thorough SPR
where every subtree is tried at every possible position in the tree.

## Substitution Model: F84

The F84 model (Felsenstein 1984) is parameterized by:
- Base frequencies: `freqa`, `freqc`, `freqg`, `freqt`
- Transition/transversion ratio: `ttratio`

Derived parameters (computed in `getbasefreqs()` in seq.c):
- `freqr = freqa + freqg` (purine frequency)
- `freqy = freqc + freqt` (pyrimidine frequency)
- `freqar = freqa / freqr`, `freqgr = freqg / freqr`
- `freqcy = freqc / freqy`, `freqty = freqt / freqy`
- `xi` and `xv` - eigenvalue-related parameters

The transition probability matrix P(t) is factored as:
```
P(i,j;t) = freq_j * (1 - exp(-xv*t))             [any base]
          + freq_j/freq_class * (exp(-xv*t) - exp(-(ratxi+ratxv)*t))  [same class]
          + delta(i,j) * exp(-(ratxi+ratxv)*t)      [same base]
```

## Rate Heterogeneity

Three rate variation models:
1. **Gamma-distributed rates** (`gama`): Discretized gamma with `rcategs`
   categories (lines 490-491). Shape parameter `alpha = 1/cv^2`.
2. **Gamma + invariant sites** (`invar`): Gamma for variable sites plus a
   proportion `invarfrac` of invariant sites (lines 493-506).
3. **User-defined HMM**: Custom rate categories with probabilities and optional
   autocorrelation (lambda parameter).

The autocorrelation HMM (when `auto_` is true) models spatial correlation of
rates along the sequence using transition probability `lambda` between
adjacent sites (lines 886-903).

## I/O

**Input**: PHYLIP interleaved or sequential format. Read via `inputdata()`
in seq.c. Site patterns are compressed via `makeweights()` (lines 616-649) --
identical site patterns are combined with weights for efficiency.

**Output**: Newick tree format via `dnaml_treeout()` (lines 2144-2207), plus
detailed text output including branch lengths, confidence intervals, and
optionally reconstructed ancestral sequences.

## Numerical Methods

1. **Log-likelihood in log-space**: Underflow prevention via `underflows[]`
   array and `fix_x()` rescaling (lines 1087-1091).
2. **Newton-Raphson with safeguards**: Modified NR that forces uphill movement
   and retracts on failure (lines 1239-1286).
3. **Confidence intervals via curvature**: Uses the second derivative of the
   log-likelihood at the MLE to compute approximate confidence intervals for
   branch lengths using chi-squared critical value 3.841 (line 1794-1795).
4. **Site pattern compression**: Identical patterns combined with weights to
   avoid redundant computation (lines 616-649).
5. **Precomputed lookup tables**: `tbl[rcategs][categs]` stores precomputed
   rate*xi, rate*xv products (lines 745-822).

## Complexity

- **Time**: O(n * s * r * c) per likelihood evaluation, where n = taxa,
  s = unique site patterns, r = rate categories, c = site categories.
  Tree search is O(n^2) insertion attempts, each requiring O(n) traversal
  and likelihood evaluation, giving O(n^3 * s * r) overall.
- **Space**: O(n * s * r) for conditional likelihood vectors at all nodes.

## Modernization Notes for Rust Reimplementation

1. **Global state to structs**: The ~50 global variables should become fields
   of a `DnamlConfig` struct and a `DnamlState` struct. The tree should own
   its likelihood computation state.

2. **Fixed arrays to Vec/ndarray**: `sitelike[4]` should become a fixed-size
   array `[f64; 4]`, but the dynamic dimensions (sites, rate categories)
   should use `Vec<Vec<[f64; 4]>>` or an ndarray.

3. **Node ring structure**: The circular linked list (next pointers forming
   a ring at each internal node) should become either a proper tree structure
   with `children: Vec<NodeId>` or an arena-based approach with explicit
   parent/child relationships.

4. **Separation of concerns**: The monolithic file mixes I/O, options parsing,
   tree manipulation, likelihood computation, and optimization. These should
   be separate modules.

5. **Tree copying**: The repeated `dnamlcopy()` calls during search are
   expensive. Consider a more efficient tree representation that supports
   undo operations instead of full copies.

6. **Parallelism**: Site-level parallelism is natural -- each site's
   likelihood contribution is independent. Use rayon for parallel iteration
   over sites.

7. **SIMD**: The inner loop in nuview() operates on 4-element sitelike arrays,
   which map perfectly to SIMD instructions (SSE/AVX).

8. **Numerical precision**: Use Kahan summation or pairwise summation for
   the log-likelihood accumulation. Consider using log-space throughout
   instead of the underflow correction approach.

9. **Generic substitution model**: The F84-specific transition probability
   factoring should be abstracted to support pluggable models (GTR, etc.)
   via a trait.

10. **Memory layout**: The current `x[site][ratecateg][base]` layout causes
    poor cache behavior when iterating over sites. Consider
    `x[ratecateg][site][base]` or structure-of-arrays layout for better
    vectorization.
