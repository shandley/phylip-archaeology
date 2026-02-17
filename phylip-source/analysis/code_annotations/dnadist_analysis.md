# DNADIST Analysis - DNA Distance Matrix Computation

**Source file**: `phylip-3.698/src/dnadist.c` (1335 lines)
**Dependencies**: `phylip.h`, `seq.h`

## Overview

DNADIST computes pairwise evolutionary distances between DNA sequences under
various substitution models. The output is a distance matrix suitable for
input to distance-based tree-building programs (NEIGHBOR, FITCH, KITSCH).

The program supports five distance measures:
1. **F84** (default): Felsenstein 1984 model with unequal base frequencies
   and different transition/transversion rates
2. **Kimura 2-parameter (K2P)**: Kimura 1980, distinguishes transitions
   and transversions but assumes equal base frequencies
3. **Jukes-Cantor (JC69)**: Jukes & Cantor 1969, equal rates for all
   substitutions and equal base frequencies
4. **LogDet**: Lake 1994 / Lockhart et al. 1994, robust to base frequency
   variation across lineages
5. **Similarity**: Simple fraction of identical sites (not a true distance)

Each distance measure can optionally incorporate:
- **Gamma-distributed rate variation**: Models among-site rate heterogeneity
- **Invariant sites**: A fraction of sites are assumed to be invariable
- **Site-specific rate categories**: User-defined rate classes

## Key Data Structures

### valrec (lines 35-37)
```c
typedef struct valrec {
    double rat, ratxv, z1, y1, z1zz, z1yy, z1xv;
} valrec;
```
Precomputed rate-category-dependent values for the F84 distance estimation.
Unlike the `valrec` in DNAML (which has per-child arrays), this version is
simpler because distances are computed pairwise without a tree structure.

- `rat`: rate for this category
- `ratxv`: rate * xv (transversion rate eigenvalue)
- `z1`, `y1`, `z1zz`, `z1yy`, `z1xv`: precomputed exponentials used in
  the EM-like distance iteration

### Per-Taxon Conditional Likelihood (lines 668-787)
```c
nodep[i]->x[site][0][base]  // conditional "likelihood" at each site
```
For unambiguous bases, this is a 0/1 indicator. For ambiguous bases (IUPAC
codes), multiple positions are set to 1.0. This allows proper handling of
ambiguity in distance computation.

### Distance Matrix (lines 388-390)
```c
double **d;  // d[spp][spp] - the computed distance matrix
```

### Precomputed Tables (lines 49-53)
```c
double rate[maxcategs];           // rate for each category
valrec tbl[maxcategs];            // precomputed rate values
double *weightrat;                // weight[i] * rate[category[i]]
double sumweightrat;              // sum of weighted rates
```

### Key Global Variables (lines 41-48)
- `xi`, `xv`: F84 model eigenvalue parameters (transition/transversion)
- `ttratio`: transition/transversion ratio (default 2.0)
- `freqa`, `freqc`, `freqg`, `freqt`: base frequencies
- `freqr`, `freqy`: purine/pyrimidine frequencies
- `freqar`, `freqcy`, `freqgr`, `freqty`: conditional base frequencies
- `cvi`: gamma distribution shape parameter (1/alpha)
- `invarfrac`: fraction of invariant sites
- `fracchange`: expected fraction of changed sites per unit time

## Core Algorithms

### 1. makev() - Compute One Pairwise Distance (lines 922-1213)

This is the central function, computing the distance between species `m`
and species `n`. It dispatches to different algorithms based on the model.

**Overlap check** (lines 942-954): First verifies that the two sequences
share at least one informative site (both are non-gap, non-ambiguous).

**Quick path detection** (lines 956-981): Determines if all sites are
unambiguous (sum of indicators = 1.0) or completely ambiguous (sum = 4.0).
If so, a faster "quick" computation is used that works with integer counts
rather than floating-point likelihoods.

#### Jukes-Cantor Distance (lines 1002-1011)

**Without rate variation** (line 1004):
```
d = -0.75 * ln((4p - 1) / 3)
```
where `p = numerator/denominator` is the fraction of identical sites.

**With gamma** (lines 1006-1007):
```
d = 0.75 * alpha * (((4p - 1) / 3)^(-1/alpha) - 1)
```

**With gamma + invariant sites** (lines 1009-1011):
```
d = 0.75 * alpha * (((4 * (p - invarfrac)/(1 - invarfrac) - 1) / 3)^(-1/alpha) - 1)
```

#### Kimura 2-Parameter Distance (lines 1013-1088)

Uses an EM-like iterative procedure to estimate the distance:

1. Count transitions (`num2`) and identical sites (`num1`) (lines 1013-1037)
2. Initialize `tt = 0.1`, `delta = 0.1` (lines 1040-1041)
3. Iterate up to 100 times (`iterationsd`):
   a. Compute expected probabilities under current distance:
      ```
      p1 = exp(-tt)                    // same base
      p2 = exp(-xv*tt) - exp(-tt)      // same class (transition)
      p3 = 1 - exp(-xv*tt)             // transversion
      ```
      Or with gamma:
      ```
      p1 = (1 + tt/alpha)^(-alpha)
      p2 = (1 + xv*tt/alpha)^(-alpha) - (1 + tt/alpha)^(-alpha)
      p3 = 1 - (1 + xv*tt/alpha)^(-alpha)
      ```
   b. Compute expected proportions q1, q2, q3 with invariant site correction
   c. Compute slope of log-likelihood with respect to tt
   d. Update: `if slope < 0: delta = -|delta|/2; else: delta = |delta|`
   e. `tt += delta`

This is a bisection-like method guided by the sign of the slope.

#### F84 Distance (lines 1089-1172)

The general case for F84 (and non-quick K2P) uses a similar iterative
approach but works with the full probability arrays:

1. Precompute per-site products (lines 1093-1110):
   - `prod[i]`: product of marginal base frequency sums
   - `prod2[i]`: purine-purine + pyrimidine-pyrimidine product
   - `prod3[i]`: same-base product (weighted by frequencies)

2. Iterate distance estimation (lines 1114-1158):
   ```c
   for each category:
       z1 = exp(ratxv * lz)         // transversion component
       z1zz = exp(rat * lz)         // total rate component
       y1 = 1 - z1
       z1yy = z1 - z1zz
   for each site:
       slope += weightrat[i] * (z1zz*(bb-aa) + z1xv*(cc-bb))
                / (aa*z1zz + bb*z1yy + cc*y1)
   ```
   With gamma, the exponentials are replaced by gamma-transformed versions:
   ```
   z1 = (1 - ratxv*lz/alpha)^(-alpha)
   ```

3. Update: same bisection-like approach as K2P.

#### LogDet Distance (lines 1173-1202)

LogDet is a non-parametric distance based on the determinant of the
base substitution matrix:

1. Build a 4x4 joint base frequency table (lines 1180-1189):
   ```c
   basetable[k][l] += weight[i]  // where k=base in seq1, l=base in seq2
   ```

2. Compute log-determinant via Gauss-Jordan elimination (line 1191):
   ```c
   vv = lndet(basetable)
   ```

3. Correct for marginal frequencies (lines 1198-1201):
   ```
   d = -0.25 * (lndet - 0.5 * (sum of log marginal freqs))
   ```

#### lndet() - Log-Determinant Computation (lines 891-919)

Computes the log of the determinant of a 4x4 matrix using Gauss-Jordan
elimination (in-place inversion):

```c
for each row i:
    ld *= a[i][i]           // accumulate determinant
    normalize row i by a[i][i]
    for each other row j:
        subtract scaled row i from row j
return log(ld)
```

Returns 99.0 as a sentinel if the determinant is non-positive (indicating
the sequences are too divergent or have anomalous composition).

#### Similarity (lines 1203-1211)

Simple fraction of identical sites: `vv = numerator / denominator`.
Not a true evolutionary distance -- included for comparison purposes.

### 2. makedists() - Compute Full Distance Matrix (lines 1216-1276)

Iterates over all pairs of species, calling `makev()` for each:

```c
for i = 0 to spp-1:
    for j = i+1 to spp:
        makev(i, j, &v)
        d[i][j] = d[j][i] = |v|
```

Bad distances (from convergence failure, negative determinant, etc.) are
set to -1.0 and a warning is printed.

### 3. dnadist_empiricalfreqs() - EM Base Frequency Estimation (lines 790-824)

Estimates base frequencies from the data using 8 iterations of EM:

```c
for k = 1 to 8:
    for each species i, each site j:
        sum = freqa*x[A] + freqc*x[C] + freqg*x[G] + freqt*x[T]
        suma += w * freqa * x[A] / sum
        // similarly for C, G, T
    freqa = suma / total_sum
    // similarly for C, G, T
```

This handles ambiguous bases correctly: when `x[A] = x[G] = 1.0` (R = purine),
the EM algorithm apportions the observation between A and G in proportion
to their current estimated frequencies.

### 4. makeweights() - Site Pattern Compression (lines 633-665)

Compresses identical site patterns for computational efficiency:

1. `dnadist_sitesort()`: Shell sort of sites by pattern (lines 532-566)
2. `dnadist_sitecombine()`: Merge identical adjacent patterns (lines 569-594)
3. `dnadist_sitescrunch()`: Move representatives to front (lines 597-630)
4. Compute `weight[]`, `location[]`, and normalize rates (lines 648-664)

After compression, `endsite` unique patterns remain. Each pattern's weight
equals the sum of the original sites' weights that map to it.

### 5. inittable() - Precompute Rate Tables (lines 879-888)

```c
for each category i:
    tbl[i].rat = rate[i]
    tbl[i].ratxv = rate[i] * xv
```

### 6. getbasefreqs() - F84 Model Parameters (from seq.c)

Computes derived model parameters from base frequencies and the
transition/transversion ratio:

- `freqr = freqa + freqg` (purine frequency)
- `freqy = freqc + freqt` (pyrimidine frequency)
- `xi` and `xv`: eigenvalue-related parameters of the F84 rate matrix
- `fracchange`: expected fraction of changed sites per unit branch length

## Rate Heterogeneity

### Gamma-distributed Rates
When `gama` is true, the gamma distribution (shape = alpha = 1/cvi^2) is
used to model among-site rate variation. The distance formula uses the
gamma-transformed probability:

```
P(invariant | rate r, time t) = (1 + r*t/alpha)^(-alpha)
```

instead of the exponential `exp(-r*t)`. This accounts for the fact that
some sites evolve faster than others.

### Invariant Sites
When `invar` is true, a fraction `invarfrac` of sites is assumed invariable.
The probability of observing identity becomes:

```
P(identity) = invarfrac + (1-invarfrac) * P(identity | variable)
```

### Combined Gamma + Invariant
Both can be active simultaneously (`invar` implies rate variation for the
variable sites follows a gamma distribution).

## Ambiguity Handling

DNADIST handles all IUPAC ambiguity codes (lines 686-783):
- Standard bases: A, C, G, T/U -> single 1.0 in the appropriate position
- Two-fold degenerate: M(AC), R(AG), W(AT), S(CG), Y(CT), K(GT)
- Three-fold degenerate: B(CGT), D(AGT), H(ACT), V(ACG)
- Fully ambiguous: N, X, ?, O, - -> all positions 1.0

For JC69, K2P, and LogDet, ambiguous sites are excluded from the "quick"
computation and handled by the general F84-like iteration. For LogDet
with ambiguous sites, the program issues a warning and writes -1.0 (line
982-987), as LogDet requires unambiguous data.

## I/O

**Input**: PHYLIP format DNA sequences (interleaved or sequential).
Optional weights and categories files.

**Output**: Distance matrix in one of three formats (controlled by
`matrix_flags`):
- `MAT_MACHINE`: Full square matrix (default, machine-readable)
- `MAT_LOWERTRI`: Lower-triangular matrix
- `MAT_HUMAN`: Human-readable format with formatting

Output is via `output_matrix_d()` from `phylip.c`, using the
`stringnames_new()`/`stringnames_delete()` utilities (lines 1283-1286).

## Numerical Methods

1. **EM for base frequencies**: 8 iterations (line 800) of EM are used to
   estimate base frequencies from ambiguous data. This converges quickly
   because the problem is convex.

2. **Bisection-like distance iteration**: The iterative distance estimation
   (up to 100 iterations, `iterationsd`) uses a method where:
   - If the slope is positive, step forward
   - If the slope is negative, halve the step size and reverse direction
   This is more robust than Newton-Raphson for this problem but slower.

3. **Convergence criterion**: `|delta| > 0.0000002` (line 1043/1114).
   Divergence detection: if `delta >= 0.1` after all iterations, the
   distance is flagged as bad.

4. **Log-determinant stability**: The `lndet()` function uses Gauss-Jordan
   elimination without pivoting (lines 900-914). This is acceptable for
   4x4 matrices but would be numerically unstable for larger matrices.

5. **Zero-frequency protection**: Base frequencies below 1e-8 are bumped
   to 1e-6 (lines 841-864) to prevent division by zero in the F84
   formula.

## Complexity

- **Time**: O(spp^2 * endsite * categs * iterationsd) for F84/K2P distances.
  JC69 and LogDet are O(spp^2 * endsite) (no iteration needed).
  With n species and s sites: O(n^2 * s) for simple models,
  O(n^2 * s * 100) worst case for iterative models.
- **Space**: O(n * s) for the sequence data, O(n^2) for the distance matrix.

## Modernization Notes for Rust Reimplementation

1. **Distance model as trait**:
   ```rust
   trait DistanceModel {
       fn compute(&self, seq1: &[Base], seq2: &[Base], weights: &[f64]) -> Result<f64, DistanceError>;
   }
   ```
   Implement for JC69, K2P, F84, LogDet, and Similarity.

2. **Proper Newton-Raphson**: Replace the bisection-like iteration with
   actual Newton-Raphson (using both slope and curvature) for faster
   convergence. Fall back to bisection when NR diverges.

3. **Parallel distance computation**: Each pair (i,j) is independent.
   Use `rayon` to parallelize the O(n^2) pair loop. This is embarrassingly
   parallel and gives near-linear speedup.

4. **SIMD for site summation**: The inner loops over sites (computing
   numerator, denominator, prod, prod2, prod3) are reduction operations
   over arrays, ideal for SIMD vectorization.

5. **Ambiguity as enum**: Replace the switch statement (lines 686-783) with
   an enum:
   ```rust
   enum Base { A, C, G, T, Ambiguous(BaseSet) }
   ```

6. **Error handling for bad distances**: Replace printf + sentinel (-1.0)
   with `Result<f64, DistanceError>`:
   ```rust
   enum DistanceError {
       NoOverlap(usize, usize),
       InfiniteDistance(usize, usize),
       NegativeDeterminant(usize, usize),
       ConvergenceFailure(usize, usize),
   }
   ```

7. **Matrix output format**: Use a proper matrix formatting library or
   trait rather than the C-style `output_matrix_d()` function. Support
   CSV and TSV in addition to PHYLIP format.

8. **Site pattern compression**: Factor this out as a reusable utility,
   shared with DNAML and other programs that benefit from pattern
   compression.

9. **Rate heterogeneity as composition**: Model rate heterogeneity as
   a wrapper around a base model:
   ```rust
   struct GammaRates<M: DistanceModel> {
       base_model: M,
       alpha: f64,
       invariant_fraction: Option<f64>,
   }
   ```

10. **LogDet with pivoting**: Use LU decomposition with partial pivoting
    for the 4x4 determinant computation, or compute it analytically
    (the 4x4 determinant has a closed-form expression).
