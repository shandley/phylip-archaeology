# Hidden Insights in PHYLIP: What Felsenstein Understood That We're Losing

*A technical analysis of the algorithms in Joe Felsenstein's PHYLIP software, reimplemented in Rust as part of the PHYLIP Archaeology project.*

## 1. The Pruning Algorithm Is a General-Purpose Computational Pattern

The most important algorithm in PHYLIP is Felsenstein's pruning algorithm (1981). On the surface, it computes the likelihood of DNA sequences on a phylogenetic tree. But its true nature is far more general.

The core recurrence is:

```
L_node(b) = PRODUCT over children c of
              SUM over states j of P_c(b, j) * L_c(j)
```

This is a **message-passing algorithm on a tree-structured graphical model** — a special case of the sum-product algorithm (belief propagation). Felsenstein published this in 1981, **seven years** before Judea Pearl formalized belief propagation (1988) and over a decade before the machine learning community embraced variational inference on graphical models.

The algorithm is not limited to 4 nucleotide states. It works on any discrete state space: 20 amino acids, 61 codons, or arbitrary categorical variables. It is not limited to molecular evolution. It works on **any tree-structured process** where:

- Observed data sits at the leaves
- Hidden states exist at internal nodes
- Transitions between states follow a known probabilistic model
- The tree topology and branch lengths are given

### Modern applications beyond phylogenetics

| Domain | Leaves (observed) | States | Tree structure |
|--------|-------------------|--------|----------------|
| Tumor phylogenetics | Single-cell genotypes | Mutation states | Clonal lineage tree |
| Cell lineage tracing | CRISPR barcode states | Barcode characters | Cell division tree |
| Language evolution | Modern word forms | Proto-language states | Language family tree |
| Gene expression evolution | Expression levels (binned) | Expression categories | Species tree |
| Viral quasispecies | Sampled sequences | Nucleotide states | Transmission tree |
| Cultural evolution | Present/absent traits | Cultural trait states | Population tree |

Most modern tools in these domains reinvent the pruning algorithm (often poorly). Felsenstein had the clean, numerically stable version — with underflow prevention via log-scaling — in 1981.

### Implementation

See [`phylip-rs/src/likelihood/pruning.rs`](phylip-rs/src/likelihood/pruning.rs). Key design features:

- **Underflow prevention**: When conditional likelihoods become extremely small (deep trees, many taxa), values are rescaled and a log-scaling factor is accumulated per site. This avoids floating-point underflow without expensive log-space arithmetic.
- **Model-agnostic**: The algorithm takes a `SubstitutionModel` trait object, decoupling the tree traversal from the specific probability model.
- **Postorder traversal**: The algorithm processes nodes from tips to root, ensuring children are computed before parents — a natural fit for any bottom-up tree computation.

---

## 2. The F84 Model: Exploiting Symmetry to Avoid Matrix Exponentiation

The F84 substitution model ([`phylip-rs/src/likelihood/models.rs`](phylip-rs/src/likelihood/models.rs)) reveals a deep insight about parameterization.

Most modern phylogenetics software computes transition probabilities via general matrix exponentiation: P(t) = exp(Qt), which requires eigendecomposition of the rate matrix Q. This is numerically fragile for ill-conditioned matrices and computationally expensive.

Felsenstein took a different approach. He observed that the purine/pyrimidine biochemical distinction creates a **block symmetry** in the rate matrix. By decomposing transitions into three components — identity, within-class (transition), and between-class (transversion) — each with its own exponential decay:

```
P(i->j; t) = delta(i,j) * exp1
            + same_class(i,j) * (pi_j / pi_class) * (exp2 - exp1)
            + pi_j * (1 - exp2)
```

he obtained **closed-form transition probabilities** that require only two `exp()` calls instead of a full eigendecomposition. This is mathematically equivalent to what physicists do with Lie algebra decompositions, but Felsenstein derived it from biological reasoning.

### Why this matters

For models with biological symmetry — and many biologically meaningful models *do* have symmetry — the Felsenstein approach of exploiting structure analytically is:

1. **Faster**: Two `exp()` calls vs. eigendecomposition + matrix multiply
2. **More numerically stable**: No ill-conditioned eigenvector matrices
3. **More interpretable**: Parameters map directly to biological quantities (transition/transversion ratio, base frequencies)

This principle applies to codon models with synonymous/nonsynonymous structure, amino acid models with chemical groupings, and epigenetic state models with methylation/demethylation asymmetry.

---

## 3. Site-Pattern Compression: A Data Structure Insight Worth Stealing

[`phylip-rs/src/likelihood/optimized.rs`](phylip-rs/src/likelihood/optimized.rs) implements `CompressedAlignment`: group identical alignment columns and weight them.

For a 1000-taxon, 10,000-site alignment, there might be only 500-2000 unique column patterns. Since the pruning algorithm is O(sites × nodes × states²), replacing "sites" with "unique patterns" gives a **5-10x speedup with zero approximation error**.

### The broader principle

This is the same insight as:
- **Suffix arrays** in genome assembly
- **Feature hashing** in machine learning
- **Run-length encoding** in data compression
- **Column-oriented storage** in databases (Parquet, Arrow)

The general rule: **before doing expensive per-column computation on categorical data, compress identical columns**.

### Where modern bioinformatics should use this but doesn't

| Data type | Column = | Typical redundancy |
|-----------|----------|--------------------|
| Variant call matrices (VCF) | Genotype pattern across samples | High (many monomorphic sites) |
| Microbiome OTU tables | Abundance pattern across samples | Moderate |
| Single-cell count matrices | Expression pattern across cells | Very high (many zeros) |
| Methylation arrays | Methylation state across CpG sites | High |
| GWAS genotype matrices | SNP genotype pattern | Moderate to high |

Tools like IQ-TREE and RAxML adopted pattern compression for phylogenetics. But most per-site statistical tests (GWAS, selection scans, differential methylation analysis) still process every column independently.

---

## 4. The Bootstrap: Weight Vectors, Not Physical Resampling

Felsenstein's bootstrap (1985) is widely used but often misunderstood. The key insight in the PHYLIP implementation is that **you never need to physically resample the data**.

A bootstrap replicate is fully described by a **weight vector** — an integer for each site saying how many times it was sampled. A site with weight 3 contributes 3× to the likelihood. Combined with site-pattern compression, this means:

1. Generate a weight vector (microseconds)
2. Re-compress with weights (fast)
3. Recompute likelihood (same cost as one evaluation)

No data copying, no new alignment objects.

### Block bootstrap and spatial correlation

The implementation also includes a **circular block bootstrap**: sample contiguous blocks of sites, wrapping around at sequence boundaries. Felsenstein recognized that adjacent sites in a gene may not be independent due to linkage and secondary structure.

This is the same principle as the **moving block bootstrap** in time series econometrics (Künsch 1989), but Felsenstein applied it to biological sequences in 1985.

**Where this matters today**: Most single-cell and spatial omics methods don't account for spatial or genomic correlation in their confidence estimates. The block bootstrap is the simplest correct approach, and it's been available since 1985.

---

## 5. Discrete Gamma Rates: A General Mixture Model Framework

Yang's (1994) discrete gamma model, implemented in [`phylip-rs/src/likelihood/gamma.rs`](phylip-rs/src/likelihood/gamma.rs), accounts for among-site rate variation. The underlying insight is Felsenstein's: **a single substitution rate across all sites is biologically absurd**.

But look at what the computation actually is:

```
L_site = (1/k) * sum_{c=1}^{k} L(site | rate_c)
```

This is a **mixture model** where the site likelihood is a weighted average over hidden rate categories. The log-sum-exp computation used for numerical stability:

```rust
let max_lnl = site_lnls.iter().copied()
    .fold(f64::NEG_INFINITY, f64::max);
let sum_exp: f64 = site_lnls.iter()
    .map(|&lnl| (lnl - max_lnl).exp()).sum();
```

is **identical to the computation in**:
- Gaussian mixture models (EM algorithm, E-step)
- Hidden Markov models (forward algorithm)
- Variational autoencoders (ELBO computation)

Felsenstein and Yang were doing mixture modeling with log-sum-exp stabilization decades before deep learning made it standard practice.

### The framework generalizes

The gamma rate heterogeneity framework accommodates **any source of hidden heterogeneity** in site-wise analysis:

| Source of heterogeneity | What varies | Alpha interpretation |
|-------------------------|-------------|---------------------|
| Functional constraint | Substitution rate | Low alpha = many constrained sites |
| Selection pressure | dN/dS ratio | Gamma-distributed omega |
| Methylation dynamics | CpG mutation rate | Rate variation across CpG sites |
| Chromatin accessibility | Mutation rate | Open vs. closed chromatin |
| RNA secondary structure | Evolutionary rate | Paired vs. unpaired sites |

### Mathematics from first principles

The implementation builds the entire gamma distribution machinery from scratch — Lanczos approximation for ln(Gamma), series and continued fraction representations for the incomplete gamma function, Halley's method for inversion, and Abramowitz & Stegun's rational approximation for the normal quantile function. This is a complete mathematical toolkit implemented in ~300 lines of Rust with no dependencies.

---

## 6. Fitch Parsimony: Bitwise Computation on Trees

The Fitch algorithm in [`phylip-rs/src/parsimony/wagner.rs`](phylip-rs/src/parsimony/wagner.rs) represents nucleotide state sets as bit masks:

```rust
pub struct StateSet(u8);  // Bit 0=A, 1=C, 2=G, 3=T, 4=gap

pub fn intersection(self, other: StateSet) -> StateSet {
    StateSet(self.0 & other.0)
}
pub fn union(self, other: StateSet) -> StateSet {
    StateSet(self.0 | other.0)
}
```

The entire parsimony score computation reduces to: if AND is non-zero, take AND (free); else take OR (cost +1). This computes minimum character state changes on the whole tree in O(sites × nodes) time — **linear in both dimensions**.

### Where parsimony beats likelihood

Modern genomics has largely abandoned parsimony for probabilistic methods. But for several applications, parsimony gives exact answers to the right question:

- **Somatic mutation counting in tumor phylogenies**: You want minimum changes, not maximum likelihood
- **CRISPR barcode lineage tracing**: Character-based data, not sequence substitution
- **Pangenomics**: Gene presence/absence on a tree — minimum gain/loss events
- **Antimicrobial resistance tracking**: Minimum acquisition/loss events for resistance genes
- **Ancestral state reconstruction**: When the tree is known and you want the most parsimonious explanation

The bitwise representation also scales naturally:
- `u64` → 64-state alphabets
- SIMD (`__m256i`) → 256 sites processed per instruction
- GPU → millions of sites in parallel

---

## 7. Model Selection: The Discipline of Counting Parameters

[`phylip-rs/src/likelihood/model_selection.rs`](phylip-rs/src/likelihood/model_selection.rs) implements AIC, BIC, and AICc — information criteria that penalize model complexity.

The critical function is parameter counting:

| Model | Substitution params | Branch length params | Total |
|-------|-------------------|---------------------|-------|
| JC69 | 0 | 2m - 3 | 2m - 3 |
| F84 | 4 (3 freq + 1 ts/tv) | 2m - 3 | 2m + 1 |
| + Gamma | +1 (alpha) | — | +1 |
| + Invariant sites | +1 (p_inv) | — | +1 |

where m = number of taxa.

The insight: **every free parameter costs you roughly 1-2 log-likelihood units of penalty**. This is Occam's razor made quantitative. A more complex model must improve the likelihood by at least this much per parameter to be justified.

Modern genomics routinely fits models with thousands of parameters (regularized regression, neural networks) but rarely asks: "How many of these parameters are actually justified by the data?" The PHYLIP tradition of explicit parameter counting and model comparison is a discipline the field could use more of.

---

## 8. What Felsenstein Understood That We're Losing

The deepest insight in this codebase isn't any single algorithm. It's the **relationship between the algorithms**.

Felsenstein understood that parsimony, distance methods, and maximum likelihood are all views of the same underlying problem — inferring evolutionary history from observed data. He proved when each one works and when it fails:

- **The Felsenstein zone** (1978): When two long branches are separated by a short internal branch, parsimony converges on the **wrong** tree with probability 1 as you add data. Likelihood gets it right.
- **Consistency**: ML is statistically consistent (converges to the truth with enough data). Parsimony is not, in general.
- **Efficiency**: NJ on good distances gives a reasonable topology in seconds. Use it as a starting point for ML, which refines it in minutes.

### Three things being lost to time

**1. The discipline of derivation from first principles.**
Our `gamma.rs` implements the Lanczos approximation, continued fractions, and Halley's method from scratch. Modern practitioners install scipy and call `gammainc()`. They can't debug it, modify it, or adapt it to new problems.

**2. The understanding that statistical phylogenetics is really about continuous-time Markov chains on trees.**
This abstraction applies to any tree-structured stochastic process — gene regulation cascades, language evolution, cultural transmission, epidemic spread. But it's taught as "a phylogenetics thing" rather than "a fundamental computational framework."

**3. The habit of testing against known analytical results.**
Our tests compare against hand-calculated values: JC69 transition probabilities from the formula, Yang's published gamma rates for alpha=0.5, the exponential CDF as a special case of the incomplete gamma function. This grounds code in theory. Most modern bioinformatics tests are regression tests — they verify the output hasn't changed, not that it's correct.

---

## 9. Independent Contrasts Are Kirchhoff's Circuit Laws

Felsenstein's independent contrasts algorithm (1985) solves phylogenetic regression in O(n) time instead of the naive O(n^3) matrix inversion. The algorithm's secret is that it is secretly an electrical circuit computation.

The core of the algorithm ([`phylip-rs/src/comparative/contrasts.rs`](phylip-rs/src/comparative/contrasts.rs), lines 136-188) performs a postorder traversal where each internal node produces one contrast:

```rust
let v_l = effective_length[left_id];
let v_r = effective_length[right_id];
let v_sum = v_l + v_r;

// The contrast: difference between left and right subtree values
let standardized = (x_l - x_r) / v_sum.sqrt();

// Propagate weighted average upward
node_values[node_id][c] = (x_l * v_r + x_r * v_l) / v_sum;

// The parallel resistor formula
let v_additional = (v_l * v_r) / v_sum;
effective_length[node_id] += v_additional;
```

The variance propagation formula `v_additional = (v_l * v_r) / (v_l + v_r)` is the **parallel resistor formula** from electrical engineering. This is not a metaphor — it is a formal mathematical equivalence.

### The circuit analogy

Under Brownian motion on a tree, the variance of the difference between two tips equals the total branch length along the path connecting them. This is exactly how resistance works in a circuit: the resistance between two nodes in a resistive network equals the sum of resistances along any path connecting them (for series connections) or the parallel combination for branching paths.

| Phylogenetics | Electrical circuit |
|---|---|
| Branch length | Resistance |
| Tip trait value | Boundary voltage |
| Contrast at internal node | Current through a branch |
| Variance of contrast | Equivalent resistance |
| Weighted average propagation | Kirchhoff's voltage law |
| Sum of squared standardized contrasts | Power dissipation |

When an internal node combines its two children, the effective branch length to its parent gains the parallel combination of the two children's effective lengths — exactly as two resistors in parallel produce an effective resistance of `R1*R2/(R1+R2)`. Kirchhoff's current laws state that current is conserved at each node; the independent contrasts equations state that the weighted average is the optimal estimate at each node. These are the same equation.

### Why O(n) instead of O(n^3)

The phylogenetic variance-covariance matrix V (where V_ij = sum of branch lengths from root to MRCA of taxa i and j) is dense and n×n. Inverting it directly costs O(n^3). But V has tree structure — it can be factored as a product of sparse matrices, one per internal node. The independent contrasts algorithm implicitly computes V^{-1}x by propagating messages along the tree, visiting each of the n-1 internal nodes exactly once with O(1) work per node.

This same principle appears in:
- **Gaussian process regression on trees**: O(n) exact inference via message passing
- **Kalman filtering**: O(n) inference on chain-structured models (a special case of trees)
- **Sparse Cholesky factorization**: exploiting the sparsity pattern of the precision matrix (V^{-1} is sparse when V has tree structure)
- **Factor analysis**: tree-structured latent variable models admit O(n) exact inference

Felsenstein published this in 1985. The connection to circuit theory was implicit — he derived it from statistical considerations. But the algorithm he found is the same one physicists use to solve resistor networks, and the same one probabilistic graphical model researchers would later formalize as message passing on trees.

---

## 10. The Contml Stereographic Projection: Differential Geometry on Probability Spaces

The contml program computes maximum likelihood trees for continuous characters and gene frequency data. Before running the likelihood computation, the PHYLIP source (`contml.c`) applies a remarkable coordinate transformation to gene frequency data: a **stereographic projection** from the allele frequency simplex to Euclidean space via the unit hypersphere.

### The problem

Allele frequencies at a locus are constrained: they must be non-negative and sum to 1. This means they live on a simplex, not in unconstrained Euclidean space. Brownian motion on a simplex is not isotropic — movement near the simplex boundary is constrained differently from movement in the interior. Applying a standard Gaussian likelihood model to raw frequencies would give biased results.

### Felsenstein's solution

The PHYLIP code applies a square-root transformation followed by a stereographic projection. For allele frequencies `x_1, x_2, ..., x_k` (summing to 1), define:

```
u_a = sqrt(x_a)
```

These transformed values satisfy `sum(u_a^2) = 1`, so the vector `(u_1, ..., u_k)` lives on the unit hypersphere in k-dimensional space. Under genetic drift (the natural stochastic process for allele frequencies), Brownian motion on the simplex maps approximately to Brownian motion on this sphere. The projection to Euclidean coordinates then justifies a standard Gaussian likelihood.

Our implementation of the Cavalli-Sforza chord distance ([`phylip-rs/src/models/gene_freq.rs`](phylip-rs/src/models/gene_freq.rs)) uses the same transformation:

```rust
let mut cos_angle = 0.0;
for a in 0..n_alleles {
    let product = x[a] * y[a];
    if product > 0.0 {
        cos_angle += product.sqrt();
    }
}
sum += 1.0 - cos_angle;
```

The operation `sum(sqrt(x_a * y_a))` is the inner product of the square-root-transformed vectors — the cosine of the angle between two points on the hypersphere. The "chord" in "chord distance" is literal: it is the Euclidean straight-line distance between two points on the surface of a sphere.

### Why this is differential geometry

The square-root transformation of probability distributions is the **Hellinger embedding**. The resulting metric on the sphere is the **Fisher information metric** for multinomial distributions — the natural Riemannian metric on the statistical manifold of discrete probability distributions.

Felsenstein was computing on a Riemannian manifold in 1973. The modern field of **information geometry** (Amari, 1985) formalized these ideas, but the computational technique — project a constrained probability space to an unconstrained Euclidean space via the Hellinger embedding — was already in PHYLIP.

### Where this appears outside phylogenetics

| Application | Data on simplex | Square-root transformation |
|---|---|---|
| Microbiome composition | Relative abundances | Hellinger distance between samples |
| Topic modeling | Document-topic proportions | Spherical topic embeddings |
| Spectral clustering | Normalized features | Spherical K-means |
| Text analysis | Term frequency distributions | BM25 and related measures |
| Optimal transport | Probability measures | Wasserstein distance via Hellinger |
| Compositional data analysis | Chemical/geological compositions | Centered log-ratio transform (related) |

The general principle: **when your data lives on a probability simplex, map it to a sphere before applying Euclidean methods**. The mapping is one line of code (`sqrt`), but the mathematical justification is deep.

---

## 11. The Hendy-Penny Supplement Bound: Dual Decomposition in Combinatorial Optimization

The branch-and-bound algorithm in PHYLIP's `dnapenny` finds *all* most parsimonious trees — guaranteed globally optimal, not heuristic. The algorithm is implemented in [`phylip-rs/src/parsimony/branch_and_bound.rs`](phylip-rs/src/parsimony/branch_and_bound.rs).

The standard part of branch-and-bound is straightforward: build a partial tree, compute its parsimony score as a lower bound, prune if the bound exceeds the best known complete-tree score. But the clever part — buried in the original C code's `supplement()` function — is a second, independent lower bound computed from the taxa that haven't been placed yet.

### How supplement() works

The algorithm adds taxa to a growing partial tree one at a time. At any point, some taxa are placed (in the tree) and some are not. The `supplement()` function asks: for each character (alignment column), is there any state present in the unplaced taxa that is *absent* from all placed taxa? If so, at least one extra step will be required when that taxon is eventually added, regardless of where it is attached.

```
For each site:
    placed_states    = union of states in placed taxa
    unplaced_states  = union of states in unplaced taxa
    if (unplaced_states AND (NOT placed_states)) != 0:
        guaranteed_extra_steps += 1
```

This produces a lower bound from the unplaced taxa alone: the number of character states that must appear in the final tree but cannot be explained by any taxon already placed.

### Why this is dual decomposition

The total lower bound is:

```
lower_bound = score(partial_tree) + supplement(unplaced_taxa)
```

This is a form of **dual decomposition** (also called Lagrangian relaxation), a technique from combinatorial optimization where a hard problem is decomposed into independent subproblems whose costs are additive:

- **Subproblem 1** (placed taxa): What is the minimum cost of the partial tree so far?
- **Subproblem 2** (unplaced taxa): What is the minimum additional cost that must be incurred, regardless of placement decisions?

Together, these bounds are tighter than either alone. The placed-taxa bound says "we've already spent X steps." The supplement bound says "we'll spend at least Y more steps." If X + Y exceeds the best known solution, the entire subtree of the search space is pruned.

### The ordering heuristic

The original PHYLIP code also sorts candidate insertion positions by their partial scores before trying them:

```
For each edge in the tree:
    try inserting the next taxon on this edge
    record the resulting score
Sort edges from best (lowest) to worst (highest) score
Try insertions in this order
```

This ensures the best placements are tried first, updating the bound early and causing more aggressive pruning of subsequent (worse) placements. It is a **greedy ordering heuristic inside an exact algorithm** — a common technique in integer programming (variable ordering in branch-and-bound) and constraint satisfaction (value ordering in backtracking search).

### Generalizable to

| Problem | Placed items | Supplement bound |
|---|---|---|
| Phylogenetic B&B | Taxa in partial tree | States unique to unplaced taxa |
| Traveling salesman | Cities visited so far | Minimum spanning tree on unvisited cities |
| Job scheduling | Jobs already assigned | Lower bound from unassigned job processing times |
| Set cover | Sets already chosen | Elements not yet covered |
| Constraint satisfaction | Variables assigned | Minimum cost of satisfying remaining constraints |

The general principle: **in any placement problem, compute separate lower bounds from the placed and unplaced components**. The combined bound is always at least as tight as either alone, and often much tighter.

---

## 12. Dollo Parsimony as a Max-Flow Problem

Dollo parsimony (Le Quesne 1974, Farris 1977) encodes a biological assumption as an asymmetric set operation. The assumption: a complex derived trait (like a vertebral column) can originate only once in evolutionary history, but can be lost any number of times independently. The implementation ([`phylip-rs/src/parsimony/dollo.rs`](phylip-rs/src/parsimony/dollo.rs)) reveals that this is not merely a different parsimony criterion — it transforms the problem from symmetric set intersection into an asymmetric coverage problem with a graph-theoretic interpretation.

### The inverted logic

Standard Fitch parsimony combines children by *intersection* (if non-empty) or *union* (if empty): the parent state is the set of states shared by both children, costing nothing if there is overlap.

Dollo inverts this: children combine by **OR** for the derived state. If *either* child has state 1, the parent must have state 1:

```rust
fn combine(&self, left: StateSet, right: StateSet) -> (StateSet, usize) {
    let left_has_one = !left.intersection(StateSet::T).is_empty();
    let right_has_one = !right.intersection(StateSet::T).is_empty();

    if left_has_one || right_has_one {
        // Parent must be 1 (derived state cannot re-arise)
        let mut cost = 0;
        if !left_has_one { cost += 1; }  // loss event on left branch
        if !right_has_one { cost += 1; } // loss event on right branch
        (StateSet::T, cost)
    } else {
        (StateSet::A, 0) // Both 0, parent is 0
    }
}
```

The OR-rule forces the derived state upward: if any descendant has the trait, every ancestor back to the point of origin must also have had it. Each descendant that *lacks* the trait despite having an ancestor with it must have undergone a loss — and each loss is a cost.

### The graph theory

Under Dollo's law, the set of nodes with state 1 forms a connected subtree rooted at the single gain event. The taxa with state 0 are separated from this subtree by at least one "loss" edge. The Dollo parsimony score equals the number of edges that must be "cut" (marked as loss events) to separate all the 0-taxa from the connected 1-subtree.

This is a **minimum edge cut** problem on the tree. Equivalently, by the max-flow/min-cut theorem, the Dollo score equals the maximum flow from the gain node to the set of 0-taxa, where each edge has capacity 1.

For a tree (which has no cycles), the min-cut is trivially computed by the upward OR-pass: each 0-child of a 1-parent requires exactly one cut. But recognizing the graph-theoretic structure reveals why the algorithm works and how it generalizes:

| Application | Derived state (1) | Ancestral state (0) | Loss = |
|---|---|---|---|
| Gene family evolution | Gene present | Gene absent | Gene deletion |
| Antibiotic resistance | Resistance allele present | Susceptible | Loss of resistance element |
| Island biogeography | Species present on island | Species absent | Local extinction |
| Horizontal gene transfer | Gene in genome | Gene absent | Gene loss (after single HGT) |
| Cultural trait evolution | Cultural practice exists | Practice absent | Cultural loss |

### The two-pass necessity

Unlike Fitch parsimony, Dollo requires a **downward correction pass** after the upward pass. The upward pass determines the minimum set of nodes that must carry the derived state (everything on the path from the gain point to the tips with state 1). The downward pass resolves ambiguous cases: if an ancestor is known to be 0, a descendant marked as potentially 1 by the upward pass can be corrected to 0 without introducing a re-derivation. This second pass is structurally necessary — not optional for disambiguation — because the Dollo constraint flows in both directions: upward (derived state propagates to ancestors) and downward (if the ancestor never had the trait, descendants can't have lost it).

---

## 13. LogDet: The Only Compositionally Robust Distance

The LogDet distance ([`phylip-rs/src/models/logdet.rs`](phylip-rs/src/models/logdet.rs)) has a property that no other standard phylogenetic distance measure possesses: it gives correct distances even when different lineages have different base compositions. The mathematical reason is beautiful.

### The factorization

For two sequences, construct the 4×4 divergence matrix F, where F[a][b] = proportion of sites where sequence 1 has base a and sequence 2 has base b. Under a Markov model of sequence evolution:

```
F = diag(π_i) × P(t)
```

where π_i are the base frequencies of sequence i and P(t) is the transition probability matrix over evolutionary time t. Taking the determinant:

```
det(F) = det(diag(π_i)) × det(P(t))
       = (π_A × π_C × π_G × π_T)_i × det(P(t))
```

The LogDet distance isolates the P(t) term by subtracting the frequency contributions:

```rust
let d = (-det_f.ln() + 0.5 * (prod_fi.ln() + prod_fj.ln())) / n;
```

This is:

```
d = -(1/4) × ln(det(F)) + (1/8) × (ln(∏ π_i) + ln(∏ π_j))
  = -(1/4) × ln(det(P(t)))
```

The base frequency terms cancel completely, leaving only the pure evolutionary signal encoded in det(P(t)).

### Why other distances fail

Under the JC69 or K2P models, the correction formulas assume *all* sequences have the same base composition. When GC-content varies across lineages (common in bacteria, mitochondrial DNA, and endosymbiont genomes), these formulas give biased distances — sequences with similar GC-content look artificially close, regardless of their true evolutionary relationship.

LogDet doesn't make this assumption. The divergence matrix F automatically captures whatever compositional differences exist, and the correction terms remove them exactly.

### The determinant from first principles

Our implementation computes the 4×4 determinant using Laplace cofactor expansion — no linear algebra library needed:

```rust
pub fn det4x4(m: &[[f64; 4]; 4]) -> f64 {
    let mut det = 0.0;
    for j in 0..4 {
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        det += sign * m[0][j] * minor3x3(m, 0, j);
    }
    det
}
```

For a fixed-size 4×4 matrix, cofactor expansion is exact and fast. The PHYLIP C code uses Gauss-Jordan elimination instead (tracking the determinant as the product of pivot elements), which generalizes better to larger matrices but is equivalent for the 4-state DNA case.

### Where compositional robustness matters

| Data type | Why composition varies | Impact |
|---|---|---|
| Bacterial genomes | GC-content ranges from 25% to 75% | Standard distances group by GC, not by phylogeny |
| Mitochondrial DNA | Strand-specific mutational bias | MT distances are biased for distant comparisons |
| Endosymbiont genomes | AT-enrichment from Muller's ratchet | Endosymbionts cluster artifactually |
| Ancient DNA | Deamination creates C→T bias | aDNA distances are systematically inflated |
| Metagenomic assemblies | Chimeric contigs with mixed composition | Distance artifacts from assembly errors |

---

## 14. Kitsch and the Pool-Adjacent-Violators Algorithm

The Kitsch algorithm ([`phylip-rs/src/distance/kitsch.rs`](phylip-rs/src/distance/kitsch.rs)) fits a distance matrix to an ultrametric (clock-constrained) tree. The central challenge is enforcing the ultrametric constraint: every parent node must have a height greater than its children. When the optimal unconstrained height for a parent would be *lower* than its child's height, PHYLIP resolves this by an algorithm that is, in disguise, **isotonic regression on a tree**.

### The constraint problem

In an ultrametric tree, branch lengths are determined by node heights: the branch from child to parent has length `h_parent - h_child`. All leaf heights are zero (leaves are at the present), and parent heights must exceed children's heights. The optimization target is weighted least squares:

```
WLS = Σ_{i<j} w_ij × (d_ij - 2×h_LCA(i,j))^2
```

where d_ij is the observed distance and h_LCA is the height of the lowest common ancestor of taxa i and j. The weights are w_ij = 1/d_ij² (Fitch-Margoliash weighting).

### The scrunch algorithm

The original PHYLIP C code handles constraint violations through a function called `scrunch()`. When the optimal height for a parent node P is less than the height of its child C (which would create a negative branch length), `scrunch()` merges P and C into a single "supernode" with a pooled height:

```
h_merged = (h_P × w_P + h_C × w_C) / (w_P + w_C)
```

It repeats until no violations remain:

```
do {
    find the tallest child of current node
    if child_height > node_height:
        merge them (weighted average of heights)
    else:
        stop
} while (violation found)
```

This is the **pool-adjacent-violators algorithm (PAV)**, the standard algorithm for isotonic regression, generalized from a linear sequence to a tree partial order. Standard PAV enforces a monotonicity constraint y_1 ≤ y_2 ≤ ... ≤ y_n by merging adjacent violating pairs into their weighted average. Kitsch's `scrunch()` does the same thing on a tree: if a parent-child pair violates h_parent > h_child, merge them.

### Why isotonic regression on trees matters

The tree-structured PAV algorithm applies whenever you need to enforce monotonicity constraints on a hierarchical structure:

| Application | Hierarchy | Monotonicity constraint |
|---|---|---|
| Molecular dating | Phylogenetic tree | Ancestor dates > descendant dates |
| Probability calibration | Decision tree | Parent class probability ≥ children |
| Hierarchical classification | Taxonomy tree | Superclass confidence ≥ subclass |
| UPGMA/single-linkage | Dendrogram | Merge heights increase monotonically |
| Revenue attribution | Organization hierarchy | Department revenue ≥ sum of sub-departments |
| Dose-response curves | Dose levels (linear tree) | Response increases with dose |

The PHYLIP implementation demonstrates that isotonic regression on trees is a natural operation that arises whenever you fit a hierarchical model subject to ordering constraints. The PAV algorithm guarantees the solution is optimal (minimizes weighted least squares subject to the monotonicity constraints) and runs in O(n) time.

---

## 15. Bron-Kerbosch and the Clique-Tree Equivalence

PHYLIP's `clique` program finds the maximum set of binary characters that can all evolve on a single tree without homoplasy. This sounds like a phylogenetic problem, but the algorithm used — Bron-Kerbosch for maximum clique finding — is one of the most important algorithms in combinatorial optimization. The deep insight is *why* maximum clique equals phylogenetic tree.

### Character compatibility

Two binary characters are **compatible** if at most three of the four possible combinations (0,0), (0,1), (1,0), (1,1) appear among the taxa. Our implementation ([`phylip-rs/src/compatibility/mod.rs`](phylip-rs/src/compatibility/mod.rs)):

```rust
pub fn are_compatible(matrix: &BinaryMatrix, char_i: usize, char_j: usize) -> bool {
    let mut seen = [false; 4]; // (0,0), (0,1), (1,0), (1,1)
    for taxon in 0..matrix.n_taxa {
        let si = matrix.get(taxon, char_i);
        let sj = matrix.get(taxon, char_j);
        let idx = si * 2 + sj;
        seen[idx] = true;
    }
    // Compatible if at most 3 of the 4 types are present
    seen.iter().filter(|&&s| s).count() <= 3
}
```

If all four combinations are present, no single tree can explain both characters without requiring at least one character to change state twice on some lineage (homoplasy).

### The Splits Equivalence Theorem

Each binary character defines a **bipartition** (split) of the taxa: those with state 0 and those with state 1. Buneman (1971) proved that a set of pairwise compatible splits can be represented as the edge bipartitions of a single tree. Conversely, any set containing two incompatible splits cannot be realized by any tree.

Therefore: the **maximum clique** in the compatibility graph (where vertices are characters and edges connect compatible pairs) identifies the largest set of characters that can all evolve on a single tree. And that set of characters directly implies the tree — each character in the clique corresponds to one edge (internal branch) of the inferred tree.

### The Bron-Kerbosch algorithm

Finding a maximum clique is NP-hard in general. The Bron-Kerbosch algorithm with pivoting ([`phylip-rs/src/compatibility/clique.rs`](phylip-rs/src/compatibility/clique.rs)) is an exact backtracking algorithm with aggressive pruning:

```rust
fn bron_kerbosch(adj, r, p, x, best, total_cliques) {
    if p.is_empty() && x.is_empty() {
        // r is a maximal clique
        if r.len() > best.len() { *best = r.clone(); }
        return;
    }
    // Choose pivot vertex with most neighbors in P
    let pivot = choose_pivot(adj, p, x);
    // Only try vertices NOT adjacent to the pivot
    let candidates: Vec<usize> = p.iter()
        .filter(|&&v| !adj[pivot][v])
        .copied().collect();
    for v in candidates {
        r.push(v);
        let new_p = p intersect neighbors(v);
        let new_x = x intersect neighbors(v);
        bron_kerbosch(adj, r, &new_p, &new_x, best, total_cliques);
        r.pop();
        p.remove(v);
        x.push(v);
    }
}
```

The pivot selection is the critical optimization: by choosing the vertex in P ∪ X with the most neighbors in P, the algorithm minimizes the number of recursive calls. For each vertex u in P that *is* adjacent to the pivot, we know that u can be added to any clique containing the pivot — so we don't need to explore them as independent starting points.

### The broader pattern

The **clique = compatible structure** equivalence generalizes:

| Domain | Vertices | Compatibility | Maximum clique = |
|---|---|---|---|
| Phylogenetics | Binary characters | No 4-gamete violation | Largest perfect-phylogeny character set |
| Interval scheduling | Time intervals | Non-overlapping | Maximum independent set of tasks |
| Graph coloring | Vertices | Not adjacent | Maximum independent set |
| Register allocation | Variables | Non-interfering | Variables sharing a register |
| Protein threading | Residue-structure pairs | Consistent | Best structural alignment |

---

## 16. Lake's Invariants and the Birth of Algebraic Phylogenetics

Lake's invariants ([`phylip-rs/src/invariants/lake.rs`](phylip-rs/src/invariants/lake.rs)) are polynomial functions of site-pattern frequencies that equal zero under one tree topology and are non-zero under the alternatives. They work only for 4 taxa, and they are statistically inconsistent. Yet they represent one of the most intellectually deep ideas in phylogenetics — one that connects PHYLIP to modern algebraic geometry.

### Why transversions are the key

For 4 taxa, there are 4^4 = 256 possible site patterns. Lake's method classifies them by the **transversion** pattern: which pairs of taxa differ by a transversion (purine ↔ pyrimidine) versus a transition (purine ↔ purine or pyrimidine ↔ pyrimidine).

Transitions (A↔G, C↔T) occur frequently and saturate rapidly. Two distantly related sequences will show ~50% transitions simply from multiple substitutions. Transversions (A↔C, A↔T, G↔C, G↔T) are rarer and retain phylogenetic signal longer.

The invariant for topology T1: (1,2)|(3,4) asks:

```
Are transversion patterns xxyy (where taxa 1,2 share one purine/pyrimidine
class and taxa 3,4 share the other) equally frequent as xyxy and xyyx?
```

Under topology T1, patterns xxyy are expected to be more frequent than xyxy or xyyx (because the internal branch groups taxa 1,2 together). The linear invariant measures this excess:

```
L1 = count(xxyy-type patterns supporting T1)
   - count(xxyy-type patterns supporting T2 or T3)
```

If T1 is the true topology, L1 should be the largest of the three invariant values.

### The algebraic geometry connection

A phylogenetic model (topology + substitution model) defines a parametric family of probability distributions over site patterns. This family is a **semi-algebraic variety** in the 256-dimensional probability simplex — a set defined by polynomial equalities and inequalities.

An **invariant** is a polynomial that vanishes on one variety but not on others. Finding invariants is equivalent to computing the **ideal** of the variety in polynomial algebra. This is exactly what computational algebraic geometry does.

Cavender's (1978) quadratic invariants make this explicit. Under the symmetric (CFN) two-state model, the invariant for topology T1 is:

```
f(RRRY) × f(RYRR) - f(RRRR) × f(RYYY) = 0
```

This is a 2×2 determinant — the condition for two rows of a contingency table to be independent. Independence of the marginals is exactly the conditional independence implied by the tree: taxa 1,2 are independent of taxa 3,4 given the state at the internal node.

### Why this matters beyond phylogenetics

The connection between conditional independence and polynomial invariants is the foundation of **algebraic statistics** (Sturmfels, Pachter, and collaborators, 2000s onward):

- **Bayesian networks**: Conditional independence relations define polynomial constraints on joint distributions
- **Latent variable models**: Hidden variables create polynomial relations among observed marginals
- **Tensor decomposition**: The rank of a probability tensor encodes model complexity
- **Model identifiability**: Algebraic methods determine whether model parameters can be recovered from data

Felsenstein's PHYLIP implemented the first practical use of phylogenetic invariants. Modern algebraic phylogenetics extends the idea to 5+ taxa, codon models, and network (non-tree) models. The key insight — that evolutionary models are algebraic varieties and trees are their defining ideals — was present in embryonic form in PHYLIP's `dnainvar`.

---

## 17. The O(n) Brownian Motion Likelihood: Edge-by-Edge Gaussians

The contml log-likelihood computation ([`phylip-rs/src/comparative/contml.rs`](phylip-rs/src/comparative/contml.rs)) avoids constructing the n×n variance-covariance matrix by decomposing the multivariate Gaussian into a product of independent univariate Gaussians, one per internal edge of the tree.

### The decomposition

For a Brownian motion process on a tree with n tips and n-1 internal edges, the joint likelihood of all tip values factors as:

```
ln L = Σ_{edges e} [-(df/2) × ln(v_e) - (1/(2v_e)) × Σ_c (view_left[c] - view_right[c])²]
```

where v_e is the total variance along edge e (sum of the two branch lengths on either side), and the "view" vectors are variance-weighted averages of the data in each subtree. Each edge contributes a term that looks like the log-density of a one-dimensional Gaussian:

```
ln N(Δ; 0, v) = -(1/2) × ln(2πv) - Δ²/(2v)
```

where Δ = view_left - view_right is the difference between the two subtree averages.

### The view propagation

The `nuview()` function (in the original C code) and our `independent_contrasts()` function both compute the "view" at each internal node as a variance-weighted average:

```rust
// Weight for left child is proportional to RIGHT child's variance
// (and vice versa) — the more certain subtree gets more weight
node_values[node_id][c] = (x_l * v_r + x_r * v_l) / v_sum;
```

The effective variance propagation uses the parallel combination:

```rust
effective_length[node_id] += (v_l * v_r) / v_sum;
```

This means the view at each node is the **best linear unbiased predictor (BLUP)** of the ancestral value, and the effective variance is the **prediction variance**. The product of all edge-wise Gaussians equals the joint multivariate Gaussian — but computed in O(n) instead of O(n³).

### Why this is profound

The contml algorithm anticipates a major result in spatial statistics and machine learning: **Gaussian process regression on trees is linear-time**. For a Gaussian process with a tree-structured covariance kernel (which is exactly what Brownian motion on a phylogenetic tree defines), the marginal likelihood can be computed by message passing in O(n) time, avoiding the usual O(n³) matrix operations.

This principle extends to:
- **Phylogenetic mixed models** (PHYLOLM, etc.): Fast regression with phylogenetic correlation
- **Spatial statistics on dendrogram-structured domains**: E.g., hierarchical administrative regions
- **Hierarchical Bayesian models**: When the random effects structure is tree-shaped

Felsenstein published the O(n) Brownian ML algorithm in 1973. The general theory of linear-time Gaussian inference on tree-structured precision matrices wasn't formalized until the graphical models literature of the 1990s-2000s.

---

## 18. The Felsenstein-Churchill HMM for Rate Variation

PHYLIP's `dnaml` program implements not one but two models of among-site rate variation. Yang's (1994) discrete gamma model (covered in Insight #5) treats each site's rate as an independent draw from a gamma distribution. The alternative — the Felsenstein-Churchill (1996) hidden Markov model — assumes rates are **autocorrelated along the sequence**.

### The biological motivation

Adjacent nucleotide positions in a gene do not evolve independently. Sites in the same structural element (helix, loop, binding pocket) tend to have similar evolutionary rates. A site buried in the protein core is constrained, and its neighbors are likely also in the core. Yang's model ignores this correlation; the HMM model captures it.

### The algorithm

The HMM model assigns each site to one of k rate categories, but now the rate category at site s+1 depends on the category at site s via a Markov transition matrix. The likelihood computation becomes:

1. **Forward pass**: For each site s, compute the probability of the data at sites 1..s given each possible rate category at site s.
2. **Backward pass**: Compute posterior probabilities of each rate category at each site.
3. **Parameter estimation**: Update the rate category transition probabilities.

The forward pass at each site is:

```
α_s(c) = L(site_s | rate_c) × Σ_{c'} T(c'→c) × α_{s-1}(c')
```

where L(site_s | rate_c) is the site likelihood under rate category c (computed by the pruning algorithm), T is the rate transition matrix, and α is the forward variable.

### Computational cost

The remarkable thing is that the HMM model adds essentially no cost to the per-site pruning computation. The forward-backward algorithm is O(n_sites × k²) where k is the number of rate categories (typically 4-8). Since k is small and fixed, this is O(n_sites) — the same asymptotic cost as the independent-rates model. But the biological realism is substantially better.

### Where rate autocorrelation matters

| Data type | Why rates are autocorrelated | Consequence of ignoring |
|---|---|---|
| Protein-coding genes | Structural domains have uniform constraint | Overestimates rate variation across sites |
| rRNA genes | Stem/loop structure creates rate blocks | Stem sites appear over-conserved |
| Viral genomes | Functional modules evolve as units | Incorrect site-specific rate estimates |
| Whole genomes | Isochore structure, chromatin domains | Phylogenetic signal varies by region |

The Felsenstein-Churchill HMM is the phylogenetic analog of the Baum-Welch algorithm for speech recognition — and it predates the widespread use of HMMs in bioinformatics for gene finding, protein structure prediction, and sequence alignment.

---

## 19. Protein Parsimony and the Genetic Code Step Matrix

The protein parsimony implementation ([`phylip-rs/src/parsimony/protein_parsimony.rs`](phylip-rs/src/parsimony/protein_parsimony.rs)) reveals an underappreciated fact about amino acid evolution: the **minimum number of nucleotide changes** required to convert one amino acid to another is not 1 — it ranges from 1 to 3 depending on the codon assignments. This creates a natural 20×20 weighted cost matrix that is far more biologically meaningful than the simple "equal costs" assumption.

### Building the step matrix from the genetic code

The algorithm constructs the genetic code from first principles, then computes the minimum nucleotide substitution cost between every pair of amino acids:

```
For amino acids i and j:
    cost(i, j) = min over all codons c_i encoding i
                     and all codons c_j encoding j
                 of hamming_distance(c_i, c_j)
```

For example, Phenylalanine (UUU, UUC) → Leucine (UUA, UUG, CUU, CUC, CUA, CUG): the minimum Hamming distance is 1 (UUU→UUA, changing position 3). But Phenylalanine → Tryptophan (UGG): the minimum is 2 (UUU→UGU→UGG or UUU→UUG→UGG). And some pairs like Phenylalanine → Methionine (AUG) require 3 changes.

### The Sankoff algorithm with costs

Standard Fitch parsimony uses binary costs (0 for matching states, 1 for any mismatch). The Sankoff algorithm generalizes this to arbitrary cost matrices. At each internal node, instead of maintaining a state set, it maintains a **cost vector** — the minimum cost of explaining all data below this node, for each possible ancestral state:

```rust
// For each possible parent state p (0..20 amino acids):
// cost[p] = min over children c of
//           (cost_child_left[p_left] + step_matrix[p][p_left])
//         + (cost_child_right[p_right] + step_matrix[p][p_right])
```

The time complexity is O(sites × nodes × states²) — for proteins, states = 20, so each node requires 20² = 400 operations per site. This is 100× more expensive than Fitch's bitwise operations but captures the biological reality that not all amino acid changes are equally likely.

### What the step matrix reveals

The genetic code step matrix has striking structure:

- **Chemically similar amino acids tend to be 1 step apart**: Leu↔Ile (1 step), Asp↔Glu (1 step), Ser↔Thr (1 step)
- **The code minimizes the effect of point mutations**: Most single-nucleotide changes either produce the same amino acid (synonymous) or a chemically similar one
- **Some amino acids are more "connected" than others**: Serine has 6 codons in 2 disconnected codon blocks, making it 1 step from many different amino acids
- **The stop codons create dead ends**: Any amino acid that is only 1 nucleotide change from a stop codon is under stronger purifying selection

This structure in the genetic code has been interpreted as evidence for both **optimization by natural selection** (the code minimizes the fitness cost of translation errors) and **frozen accident** (the code reflects the historical order in which amino acids were recruited into the code).

### Beyond phylogenetics

The Sankoff algorithm with a weighted cost matrix is a general dynamic programming framework for tree-structured optimization with arbitrary transition costs. It applies to:

- **Ancestral genome reconstruction**: Cost = number of rearrangement operations (inversions, translocations)
- **Morphological character evolution**: Cost = morphological distance between states
- **Natural language processing**: Cost = edit distance between word forms in cognate detection
- **Hierarchical classification**: Cost = mis-classification penalty in a taxonomy

---

## 20. Score-Ordered Search and the Greedy Heuristic Inside the Exact Algorithm

A subtle but important optimization in PHYLIP's branch-and-bound implementation is the **ordering** of candidate positions. After all possible insertion points for the next taxon are enumerated, they are sorted by their partial parsimony scores before being explored:

```
positions = all edges in current tree
for each position: compute score if taxon inserted here
sort positions: best (lowest) score first
explore positions in this order
```

### Why ordering matters

Branch-and-bound's efficiency depends entirely on how quickly the bound is tightened. If the first position explored happens to produce a poor partial score, the bound stays loose and few subsequent positions are pruned. If the first position explored produces the (nearly) optimal score, the bound tightens immediately and most subsequent positions are pruned.

By sorting positions best-first, the algorithm ensures that:
1. The first complete tree found is likely close to optimal
2. The bound is tightened early in the search
3. Subsequent (worse) positions are pruned aggressively

This is a **greedy heuristic embedded inside an exact algorithm** — the ordering doesn't affect correctness (all positions will eventually be tried unless pruned), but it dramatically affects efficiency.

### The general principle

The technique of using a greedy heuristic to guide the search order in an exact algorithm appears throughout combinatorial optimization:

| Algorithm | Exact guarantee | Greedy ordering |
|---|---|---|
| PHYLIP B&B | All optimal trees found | Insert at best position first |
| Simplex method | Optimal LP solution | Steepest-edge pivot rule |
| A* search | Optimal shortest path | f(n) = g(n) + h(n) heuristic |
| Alpha-beta pruning | Minimax optimal | Try best moves first |
| Davis-Putnam (SAT) | Satisfying assignment | VSIDS variable ordering |
| IDA* | Optimal path | Iterative deepening with heuristic |

The common pattern: **exact algorithms benefit enormously from good heuristic guidance**, not because the heuristic changes the answer, but because it changes the order in which answers are discovered, enabling more pruning. PHYLIP's implementation demonstrates this principle in a bioinformatic context that predates many of the computer science formalizations.

---

## Summary

| Algorithm | PHYLIP year | "Rediscovered" as | Modern application |
|-----------|-------------|-------------------|-------------------|
| Pruning algorithm | 1981 | Belief propagation (1988) | Any tree-structured inference |
| F84 closed-form P(t) | 1984 | Lie algebra decomposition | Structured substitution models |
| Site-pattern compression | ~1980s | Column-oriented databases | Any per-column computation |
| Bootstrap weights | 1985 | Weighted resampling | Spatial/genomic correlation |
| Block bootstrap | 1985 | Moving block bootstrap (1989) | Correlated site data |
| Discrete gamma rates | 1994 | Mixture models / log-sum-exp | Hidden heterogeneity in omics |
| Fitch parsimony | 1971 | Bitwise set operations | Character-based lineage tracing |
| Model selection (AIC/BIC) | ~1980s | Regularization theory | Parameter counting discipline |
| Independent contrasts | 1985 | Kirchhoff's circuit laws | Gaussian processes on trees |
| Contml projection | 1973 | Hellinger embedding / information geometry | Compositional data analysis |
| Hendy-Penny supplement | 1982 | Dual decomposition / Lagrangian relaxation | Combinatorial optimization bounds |
| Dollo parsimony | 1977 | Min-cut on trees | Gene gain/loss, biogeography |
| LogDet distance | 1994 | Determinantal factorization | Compositionally robust distances |
| Kitsch scrunch | ~1980s | Pool-adjacent-violators on trees | Isotonic regression, calibration |
| Clique analysis | 1986 | Bron-Kerbosch / Splits Equivalence | Constraint satisfaction, scheduling |
| Lake's invariants | 1987 | Algebraic statistics / variety ideals | Model identifiability testing |
| O(n) Brownian ML | 1973 | Tree-structured Gaussian processes | Spatial statistics, hierarchical models |
| Felsenstein-Churchill HMM | 1996 | Baum-Welch for rate correlation | Genome annotation, speech recognition |
| Protein Sankoff | ~1990s | Weighted dynamic programming on trees | Edit distance, ancestral reconstruction |
| Score-ordered B&B | 1982 | Greedy-guided exact search (A*, alpha-beta) | SAT solving, game tree search |

Felsenstein built an intellectual framework, not just software. The algorithms are timeless. The software is archaeology. The *thinking* is what's worth preserving.

---

*This analysis is part of the [PHYLIP Archaeology](https://github.com/shandley/phylip-archaeology) project, which preserves and modernizes Joe Felsenstein's PHYLIP phylogenetics software as a pure-Rust implementation with zero external dependencies.*

## References

- Amari, S. (1985). *Differential-Geometrical Methods in Statistics*. Springer.
- Bron, C. & Kerbosch, J. (1973). Algorithm 457: finding all cliques of an undirected graph. *Communications of the ACM*, 16, 575-577.
- Buneman, P. (1971). The recovery of trees from measures of dissimilarity. In *Mathematics in the Archaeological and Historical Sciences*, Edinburgh University Press.
- Cavalli-Sforza, L.L. & Edwards, A.W.F. (1967). Phylogenetic analysis: models and estimation procedures. *American Journal of Human Genetics*, 19, 233-257.
- Cavender, J.A. (1978). Taxonomy with confidence. *Mathematical Biosciences*, 40, 271-280.
- Farris, J.S. (1977). Phylogenetic analysis under Dollo's law. *Systematic Zoology*, 26, 77-88.
- Felsenstein, J. (1973). Maximum-likelihood estimation of evolutionary trees from continuous characters. *American Journal of Human Genetics*, 25, 471-492.
- Felsenstein, J. (1978). Cases in which parsimony or compatibility methods will be positively misleading. *Systematic Zoology*, 27, 401-410.
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution*, 17, 368-376.
- Felsenstein, J. (1984). A likelihood approach to character weighting and what it tells us about parsimony and compatibility. *Biological Journal of the Linnean Society*, 16, 183-196.
- Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the bootstrap. *Evolution*, 39, 783-791.
- Felsenstein, J. (1985). Phylogenies and the comparative method. *American Naturalist*, 125, 1-15.
- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
- Felsenstein, J. & Churchill, G.A. (1996). A hidden Markov model approach to variation among sites in rate of evolution. *Molecular Biology and Evolution*, 13, 93-104.
- Fitch, W.M. (1971). Toward defining the course of evolution. *Systematic Zoology*, 20, 406-416.
- Fitch, W.M. & Margoliash, E. (1967). Construction of phylogenetic trees. *Science*, 155, 279-284.
- Hendy, M.D. & Penny, D. (1982). Branch and bound algorithms to determine minimal evolutionary trees. *Mathematical Biosciences*, 59, 277-290.
- Lake, J.A. (1987). A rate-independent technique for analysis of nucleic acid sequences: evolutionary parsimony. *Molecular Biology and Evolution*, 4, 167-191.
- Lake, J.A. (1994). Reconstructing evolutionary trees from DNA and protein sequences: paralinear distances. *PNAS*, 91, 1455-1459.
- Le Quesne, W.J. (1974). The uniquely evolved character concept and its cladistic application. *Systematic Zoology*, 23, 513-517.
- Lockhart, P.J., Steel, M.A., Hendy, M.D. & Penny, D. (1994). Recovering evolutionary trees under a more realistic model of sequence evolution. *Molecular Biology and Evolution*, 11, 605-612.
- Nei, M. (1972). Genetic distance between populations. *American Naturalist*, 106, 283-292.
- Pearl, J. (1988). *Probabilistic Reasoning in Intelligent Systems*. Morgan Kaufmann.
- Reynolds, J., Weir, B.S. & Cockerham, C.C. (1983). Estimation of the coancestry coefficient: basis for a short-term genetic distance. *Genetics*, 105, 767-779.
- Saitou, N. & Nei, M. (1987). The neighbor-joining method. *Molecular Biology and Evolution*, 4, 406-425.
- Sankoff, D. (1975). Minimal mutation trees of sequences. *SIAM Journal on Applied Mathematics*, 28, 35-42.
- Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites. *Journal of Molecular Evolution*, 39, 306-314.
