# RA0 · G4 — ASYMPTOTIC METASTABLE ALGEBRA (result)

**Central question.** When dynamics supplies a canonical slow spectral projector through a diverging gap, does its range
converge to a proper unital commutative algebra of observables? The convergence must happen without a threshold, supplied
basin labels or an observable channel.

**Evidence:**
- propositions below, each labelled KNOWN / REDERIVED, PROVED HERE, ASSEMBLED or CONJECTURE;
- `g4/g4_algebra.py` + `g4/g4_algebra.log` (independent code path, not independent reviewer; NUMERICAL ILLUSTRATION).

**Scope.** No quantum extension, collapse model, Born rule or consciousness interpretation.

## G4.1 Mathematical object

**Setting.** (X_N, P_N) is a reversible, irreducible Markov family, in continuous time (generator G_N) or discrete time
(I − P_N). Let E_N be the spectral projector below a cut whose adjacent ratio diverges; it includes the constants. Let
V_N = Ran(E_N), with dimension k = the rank of the cut. V_N is semigroup-invariant.

**Norm: L²(π_N). This is PRICED.**
- π_N is the unique stationary law. It is derived from P_N using **irreducibility**.
- E_N is an orthogonal projector in this inner product because P_N is π-self-adjoint, which uses **reversibility**.
- Nothing in G4 is measure-free.

**Defects (basis-independent).**
- Δ_alg(N) = sup { ‖(I − E_N)(fg)‖_π : f, g ∈ V_N, ‖f‖_π = ‖g‖_π = 1 }.
  - The log's `Delta_sup` is a numerical **lower** estimate of this supremum, obtained by alternating maximisation over
    30 starts.
- Δ_HS(N) = (Σ_{i,j} ‖(I − E_N)(u_i u_j)‖²_π)^{1/2} for any π-orthonormal basis (u_i) of V_N.
  - This is the Hilbert–Schmidt norm of the product-residual map V ⊗ V → V^⊥, so it is basis-free.
  - It is a rigorous **upper** bound: Δ_alg ≤ Δ_HS.
- No tolerance is used. Only Δ_alg(N) → 0 counts as asymptotic closure.

## G4.2 Exact finite theorem

**PROP G4-F (KNOWN / REDERIVED; the finite case of Gelfand duality).** Let X be finite and A ⊂ ℝ^X a unital subalgebra
(a linear subspace containing 1 and closed under pointwise products). Then A is exactly the set of functions constant on
the blocks of a **unique** partition of X.

*Proof.*
1. Define x ~ y iff f(x) = f(y) for all f ∈ A. The blocks of ~ form a partition, and A ⊂ {functions constant on blocks}.
2. Take distinct blocks B ≠ B′. Some f ∈ A separates them. Because A is unital,
   h_{B,B′} = (f − f(B′)) / (f(B) − f(B′)) ∈ A, with h_{B,B′} = 1 on B and 0 on B′.
3. The product over all B′ ≠ B of h_{B,B′} lies in A. It equals 1 on B and 0 on every other block, so 1_B ∈ A.
4. Therefore A = span{1_B}.
5. Uniqueness: the blocks are recovered from A as its minimal nonzero idempotents. ∎

**Corollary.** An exact slow algebra gives a canonical partition at finite N.

**Limits of the theorem.**
- It does **not** apply to a mere vector subspace: SSEP's lowest eigenspace plus constants is a subspace with product
  defect 1.0 (G3).
- At finite N the metastable V_N are **not** exact algebras. Δ_alg > 0 in every non-lumpable case logged.

## G4.3 Approximate-algebra problem

### Literature search

Metadata only, via web search; full text was not re-read in this environment.

| source | what it gives | applies to "Δ_alg → 0 ⇒ nearby partition algebra"? |
|---|---|---|
| Davis & Kahan, *SIAM J. Numer. Anal.* 7 (1970) 1–46, sin Θ theorem | eigenspace rotation ≤ ‖perturbation‖ / separation | **yes, as a tool**: gives PROP G4-T below for the nearly-uncoupled class |
| Deuflhard, Huisinga, Fischer & Schütte, *Linear Algebra Appl.* 315 (2000) 39–59 | in reversible nearly uncoupled chains the dominant eigenvectors are perturbations of block-constant vectors (basis of PCCA) | **same regime as G4-T**; their aggregation step needs a cluster count, which G4 replaces by the certified cut rank |
| Bovier, Eckhoff, Gayrard & Klein (math/0007160; CMP 2001/2002) | low-lying eigenvalues ↔ metastable sets; eigenfunctions ↔ equilibrium potentials | potential-theoretic hypotheses; not a closure theorem |
| Peng, Sun & Zanetti, COLT 2015 / *SIAM J. Comput.* 46 (2017) 710–743, structure theorem | if Υ = λ_{k+1}/ρ(k) is large, the bottom-k eigenspace is close to the indicators of an optimal k-partition | **closest general result**, but its hypothesis uses the k-way expansion ρ(k), not λ_{k+1}/λ_k (see below) |
| Lee, Oveis Gharan & Trevisan, STOC 2012 / *JACM* 61 (2014) | higher-order Cheeger: λ_k/2 ≤ ρ(k) ≤ poly(k)·√λ_k | links ρ(k) to the spectrum |
| Macgregor & Sun, ICML 2022 | tighter structure theorem | same hypothesis type |
| Johnson, *J. London Math. Soc.* 37 (1988) 294–316 (AMNM) | approximately multiplicative **maps** are near multiplicative ones | concerns maps, not almost-closed subspaces: **not directly** |
| Kadison–Kastler / Christensen near-inclusion theory (standard; not re-located this session) | closeness of two **algebras** | both objects must be algebras: **no** |
| Robeva, *SIAM J. Matrix Anal. Appl.* 37 (2016) 86–102; Anandkumar, Ge, Hsu, Kakade & Telgarsky, *JMLR* 15 (2014) 2773–2832 | odeco tensors: their robust eigenvectors are the decomposition vectors; a perturbation theorem for the robust tensor power method | **supports the blind idempotent step** in G4.5; constants not re-read |

**Why Peng–Sun–Zanetti does not settle the question.**
- The condition λ_{k+1}/λ_k → ∞ does not imply Υ → ∞.
- Lee–Oveis Gharan–Trevisan only give ρ(k) ≲ √λ_k.
- Example: in the ladder with β = 3, λ_{k+1}/√λ_k ~ L^{−1/2} → 0, even though the λ-ratio diverges.

**Outcome.** No located theorem states that "Δ_alg → 0 implies a nearby partition algebra" for general rank. What follows is
proved for rank 2, and for the nearly-uncoupled class, and is conjectured in general.

### PROP G4-P (rank 2; PROVED HERE; dimension-free)

**Statement.** Let V = span{1, f} ⊂ L²(π), where f has π-mean 0 and π-variance 1. f is unique up to sign, so it is
basis-free. Let s = π(f³), κ = π(f⁴) and δ_f = (κ − 1 − s²)^{1/2}. Then:

1. Δ_alg(V) ≥ δ_f.
2. Exactness: δ_f = 0 iff f is two-valued, the Pearson equality case.
3. If δ_f < 1, the two-block partition {f > s/2}, {f ≤ s/2} has partition algebra A_f with ‖E_V − E_{A_f}‖_π ≤ δ_f.

There is no dependence on |X|, N or π_min. The partition is computed from V and π alone, and the sign flip f → −f gives
the same partition (up to the tie set {f = s/2}).

**Proof.**
- **Part 1.** E(f²) = 1 + s·f, so (I − E)(f²) = f² − 1 − s f. Its squared norm is κ − 1 − s².
- **Part 3, rounding.** Let t_± = (s ± √(s² + 4))/2 be the roots of t² − s t − 1, and let ρ(t) be the root nearer to t
  (the midpoint is s/2). Then
  |t² − s t − 1| = |t − t₊||t − t₋| ≥ |t − ρ(t)|·(t₊ − t₋)/2 ≥ |t − ρ(t)|, since t₊ − t₋ = √(s² + 4) ≥ 2.
  Hence ‖f − ρ(f)‖_π ≤ δ_f, and ρ(f) ∈ A_f.
- **Both blocks are nonempty.** Otherwise ρ(f) would be a constant c, and ‖f − c‖²_π = 1 + c² ≥ 1 > δ_f².
- **Distance bound.** For a unit h = α + βf ∈ V (with α² + β² = 1), dist(h, A_f) = |β|·dist(f, A_f) ≤ δ_f.
- **Gap metric.** For equal finite dimensions, ‖E_V − E_A‖ = ‖(I − E_A)E_V‖. ∎

**Numerics.** In every rank-2 case logged (F2 coarse level, ladders) the bound holds. Δ_sup coincides with δ_f, and the
partition agrees with the blind argmax partition of G4.5.

### PROP G4-T (nearly-uncoupled class; ASSEMBLED from Weyl and Davis–Kahan; not externally reviewed)

**Hypotheses.** The hypotheses mention a partition 𝔅 = {B_1, …, B_K}. The **procedure never sees it**: 𝔅 enters only the
proof.

- **H1.** The chain is reversible and irreducible, so π is derived.
- **H2.** Write G = G₀ + G₁. G₀ keeps only the intra-block edges and is reversible with respect to the same π, because
  detailed balance holds edge by edge. G₁ is the inter-block part. Let g = the minimum intra-block spectral gap of −G₀
  (each block irreducible), and η = ‖G₁‖_π. Assume η/g → 0.
- **H3** (needed only for part iii). The minimum block mass p_min is bounded below, and C_V = sup_{f ∈ V, ‖f‖_π = 1} ‖f‖_∞
  is bounded. C_V is the basis-free Christoffel constant √(max_x Σ u_i(x)²).

**Conclusions.**
1. **Diverging cut at rank K.** By Weyl, λ_{K−1}(−G) ≤ η and λ_K(−G) ≥ g − η, so the ratio at rank K is
   ≥ (g − η)/η → ∞.
2. **Subspace convergence.** The kernel of −G₀ is A_𝔅 = span{1_{B_i}}. Davis–Kahan sin Θ gives
   s := ‖E_V − E_{A_𝔅}‖_π ≤ η/(g − η) → 0.
3. **Asymptotic closure.**
   - Write f = a + f′ and g = b + g′ with a, b ∈ A_𝔅. For unit f, g ∈ V:
     ‖fg − ab‖ ≤ ‖a g′‖ + ‖f′ b‖ + ‖f′ g′‖ and dist(ab, V) ≤ s‖ab‖.
   - Use ‖a‖_∞ ≤ 1/√p_min on A_𝔅, ‖f′‖_∞ ≤ C_V + 1/√p_min and ‖f′‖, ‖g′‖ ≤ s.
   - Result: **Δ_alg ≤ s·(4/√p_min + C_V) → 0.** ∎

**Analytic certificates for the logged families** (they use the construction, as a proof about the family should):

| family | g | η | η/g |
|---|---|---|---|
| F1 (G4.4) | intra conductance ≥ 0.2/n_s and ν ≤ 1.5; Dirichlet-form comparison gives g ≥ 0.2/1.5 = 2/15 | Schur bound ‖G₁‖ ≤ 2·max exit rate ≤ 2·(1.8·2·3.5M/M²)/0.5 = 50.4/M | ≤ 378/M → 0 (logged η/g ≈ 18/M) |
| ladder β = 3, site-modulated | ring gap with conductance ≥ 0.7 gives g ≳ 0.7·4π²/L² | ≤ 3L^{−3} | = O(1/L) → 0 |
| F2, fine level (4 basins) | as F1 | as F1 | → 0 |
| F2, coarse level (2 super-basins) | numerical only | numerical only | logged η/g ∝ 1/M (0.33, 0.17, 0.084, 0.042) |

For the F2 coarse level, an analytic lower bound on the super-block gap would follow from a standard decomposition
argument, but it is **not written out**, so that certificate is numerical.

### CONJECTURE G4-C (general rank)

For any finite X, π and unital V ⊂ L²(π) of dimension k, the blind partition A_blind of G4.5 satisfies
‖E_V − E_{A_blind}‖_π ≤ c_k·Δ_alg(V), plausibly with c_k = 1.

- **Status:** proved for k = 2 (PROP G4-P, with c = 1); **open for k ≥ 3**.
- **Evidence:** in all 40 logged blind runs, including the stress runs on false positives, gap/Δ_sup ≤ 0.93. In the
  positive controls the ratio is ≈ 0.4–0.5.
- Δ_sup is a lower estimate of Δ_alg, so the measured ratio is conservative.

## G4.4 Deterministic non-symmetric metastable family (repairs D3)

**F1.** Fixed rule across M, no random draws.
- **Sectors:** K = 3, with sizes M, ⌊1.5M⌋, 2M.
- **Stationary weights:** ν(x) = 1 + 0.5 cos(1.7x), so π ∝ ν is non-uniform.
- **Conductances:** with h(x, y) = 1 + 0.5 cos(0.37(x + y)) + 0.3 cos(1.1|x − y|):
  - intra-sector conductance h(x, y)/n_s;
  - inter-sector conductance h(x, y)·κ_st·M^{−2}, with κ = (1, 0.5, 2) on the three sector pairs.
- **Rates:** x → y at C(x, y)/ν(x). The chain is reversible with no permutation symmetry; sizes, weights and couplings all
  differ.

**F2** (nested). 4 basins with sizes M, 1.25M, 1.5M, 1.75M, grouped in super-basins {0,1} and {2,3}. Inter-basin
conductance is M^{−2} within a super-basin and M^{−3} across. The same h and ν are used.

**Ladder (site-modulated).** Hop rates are 1.0 and 1.6 per lane, times (1 + 0.3 cos 1.1x). Rungs are
L^{−β}(1 + 0.5 cos 0.37x).
- The G3 ladder is translation-invariant, which makes the lanes **exactly** lumpable for every β. That is a G1-type exact
  factor, not an emergent one.
- The modulated ladder is not lumpable.

**F1 results.** Blind procedure at the certified rank 3; the audit uses labels afterwards.

| M (n) | ratio at rank 3 | η/g | Δ_HS | Δ_sup | C_V | robust idempotents | gap(V, A_blind) | misassigned π-mass (audit) |
|---|---|---|---|---|---|---|---|---|
| 5 (22) | 1.20 | 3.46 | 0.61 | 0.45 | 3.89 | 3 | 0.177 | 0.31 |
| 10 (45) | 1.03 | 2.33 | 1.97 | 1.87 | 5.40 | 3 | 0.348 | 0.33 |
| 20 (90) | 1.86 | 0.93 | 0.338 | 0.311 | 2.19 | 3 | 0.145 | **0** (exact) |
| 40 (180) | 3.38 | 0.46 | 0.163 | 0.147 | 2.18 | 3 | 0.068 | **0** (exact) |
| 80 (360) | 6.44 | 0.23 | 0.080 | 0.072 | 2.15 | 3 | 0.033 | **0** (exact) |

- For M ≥ 20, Δ_alg ∝ 1/M and the ratio ∝ M, matching G4-T.
- For M = 5 and 10, η > g: the rank-3 cut has not formed yet (the ratio is about 1), and recovery fails. This is a
  **pre-asymptotic** regime, reported as such.

## G4.5 Blind recovery test

**Procedure** (`blind(lam, U, pi, rank)`). The function's signature carries no labels.

1. **Spectral data.** The procedure reads the spectrum of P_N, π (derived) and the certified cut rank k. It sets
   V = span of the k lowest π-orthonormal eigenfunctions.
2. **Defects.** Compute Δ_HS and Δ_sup.
3. **Canonical idempotents.**
   - Define the compressed product f∗g = E(fg) on V and the symmetric tensor τ(a, b, c) = π(f_a f_b f_c).
   - Find the **attracting** fixed points v of the tensor power map v ↦ τ(v, v)/‖τ(v, v)‖, with λ = τ(v, v, v). Each
     gives an idempotent e = (Vv)/λ, satisfying e∗e = e.
   - For an exact partition algebra, τ is orthogonally decomposable. Its attracting fixed points are exactly the
     primitive idempotents 1_{B_i} (Robeva; Anandkumar et al.). The unit and other non-primitive idempotents are repelling.
   - The set of attracting fixed points is fixed by V and π. Random starting points only *find* that set; they do not
     define it.
4. **Canonical rounding.** Assign x ↦ argmax_i e_i(x). The resulting partition algebra is A_blind.
5. **Distance.** Report gap(V, A_blind) = ‖E_V − E_{A_blind}‖_π.

**What is not supplied.** No K, label, coordinate, ε, crispness objective or requested cluster number.
- The **number** of robust idempotents is an output.
- It equals the cut rank in every positive-control run: 3 for F1; 2 and 4 for F2; 2 for the ladders.

**Results.**

| family | blind outcome | audit, after the selector |
|---|---|---|
| F1 | 3 blocks; gap 0.145 → 0.033 | exact match to the hidden sectors for M ≥ 20 |
| F2 coarse (rank 2) | 2 blocks; gap 0.087 → 2.8e-4 (∝ M^{−2}) | exact match to the hidden super-basins at all M |
| F2 fine (rank 4) | 4 blocks; gap 0.65 → 0.014 (∝ M^{−1}) | exact match to the hidden basins for M ≥ 10 |
| ladder β = 3, modulated | 2 blocks; Δ 7.6e-4 → 3.8e-6; gap 3.8e-4 → 1.9e-6 | exact lanes at all L |
| ladder β = 3, translation-invariant | Δ at round-off; exact lanes | an exact G1 factor, not emergent |

**Terminal: BLIND PARTITION RECOVERY.**
- Numerical for F1, F2 and the ladder.
- **Proved for rank 2** by PROP G4-P, because the rank-2 rounding is computable from V alone.
- The subspace convergence that underlies it is proved in class G4-T.
- For k ≥ 3, correctness of the rounding step is CONJECTURE G4-C.

## G4.6 False-positive controls

**Refusals.** G4 acts only on a **certified** diverging cut. Every control below is refused by a proof of
boundedness, not by numerics.

| control | certified cut? | G4 output |
|---|---|---|
| A — SSEP (L = 8, 10, 12) | no: PROP G3-B gives ratios ≤ 4 in the slow window | **REFUSE**: no partition |
| B — independent particles (N = 2 … 8) | no: PROP G3-A gives ratios ≤ 1 + θ | **REFUSE** |
| C — ladder β = 0, 1 (modulated) | no: bounded ratio ≈ 4 | **REFUSE** |
| C — ladder β = 2 (modulated) | no: bounded ratio ≈ 19.4 | **REFUSE** |
| D — ladder β = 3 (modulated) | yes: G4-T | 2-block partition (lanes) |
| E — nested F2 | yes, at ranks 2 and 4 | nested partitions; see G4.7 |

**Stress diagnostic.** Not a G4 output: the largest finite-N cut is forced and the procedure is run anyway.

| control | Δ_sup | other outcome | reading |
|---|---|---|---|
| SSEP | 1.28 / 1.30 / 1.32, flat in L | 138–147 attracting idempotents (≫ rank 3); dimension mismatch, gap = 1 | the slow space is far from any algebra |
| independent particles | ≈ 1.0 – 1.3 | V contains N *exact* idempotents (the single-particle indicators, i.e. G1's independent factors), but they do not sum to 1; gap = 1 for N ≥ 6 | no partition algebra |
| ladder β = 0, 1 | ≈ 0.71, flat | gap ≈ 0.57 – 1 | diffusive modes do not algebraize |
| **ladder β = 2** | **1.9e-2 → 7.7e-4, decreasing** | the 2-block lane partition is recovered | see below |

**The ladder β = 2 result matters.** The procedure refuses there because the cut is bounded (the no-ε rule). But the
forced slow space **does** algebraize. So asymptotic algebraization and a diverging cut are **distinct properties**:
- diverging cut ⇒ algebra in class G4-T;
- algebra without a diverging cut occurs at β = 2.

The cut rule is therefore conservative. It refuses some real partition structure rather than admitting it with a chosen
constant.

**Finite-size heuristics are unreliable.**
- A rule like "the ratio grows monotonically over the tested sizes" would **wrongly accept** SSEP, whose rank-3 ratio
  rises 1.619 → 1.716 → 1.775 toward a bounded limit.
- The same rule would **wrongly reject** F1, whose rank-3 ratio goes 1.20 → 1.03 before diverging.
- Divergence must therefore be **certified analytically**. In this gate, every certificate is a proof about the supplied
  family.

## G4.7 Hierarchy of algebras

On F2 the blind partitions at ranks 2 and 4 are compared with each other. No labels are used.

- **Nesting:** every fine block lies inside a coarse block, with **zero violating π-mass at every M**, including the
  pre-asymptotic M = 5.
- **Order:** the coarse algebra A^(2) is contained in the fine algebra A^(4), in the order of the spectral filtration.
- **Audit:** the nesting matches the hidden super-basins and basins.
- The spectral subspaces are nested automatically (V^(2) ⊂ V^(4)). The partition nesting is a nontrivial output, because
  the rounding at each level is done independently.

**Terminal: CANONICAL DYNAMICAL PARTITION FILTRATION** (numerical, in class G4-T). No level is selected.

## G4.8 Non-circularity audit

| forbidden input | used? |
|---|---|
| basin labels | **no**: labels are passed only to `audit` and `thm_T`, which are called after `blind` returns |
| spatial coordinates | no |
| cluster number K | no: the number of blocks is an output (the count of attracting idempotents) |
| spectral threshold ε | no modelling threshold. Computational equality tolerances only: zero-mode 1e-10 relative, convergence 1e-13, de-duplication 1e-6 |
| PCCA crispness objective | no: nothing is optimized; argmax of canonical idempotents |
| requested number of clusters | no |
| known partition | no |
| manually chosen rank | **no, with one priced caveat**: the rank is the rank of an *analytically certified* diverging cut, and the certificate is a proof about the supplied family (G4-T) |

**Two further priced inputs.** The L²(π) norm requires reversibility and irreducibility. The family itself is supplied.

## G4.9 Terminal

| criterion | status |
|---|---|
| 1. diverging spectral cut from D alone | **yes**, certified in class G4-T (Weyl) |
| 2. basis-free projector | **yes** |
| 3. no supplied K | **yes** |
| 4. no ε | **yes** |
| 5. slow space closes asymptotically | **proved** in class G4-T (under H3) and for rank 2; numerical rates F1 ∝ M^{−1}, F2 ∝ M^{−2} (coarse) / M^{−1} (fine), modulated ladder 7.6e-4 → 3.8e-6 over L = 25 … 200 |
| 6. genuine partition algebra recovered blindly | **numerically yes** (exact for F1 at M ≥ 20, F2 coarse at all M, F2 fine at M ≥ 10, ladder at all L); **proved for rank 2**; rank ≥ 3 rounding is CONJECTURE G4-C |
| 7. nested separations give nested partition algebras | **yes, numerically** (F2) |
| 8. false positives return no partition | **yes**, by the certified-cut rule. Caveat: at ladder β = 2 a real algebra exists and is refused |

**Terminal: CONDITIONAL DYNAMICAL PARTITION DERIVATION**, for the nearly-uncoupled reversible class G4-T. It is
conditional on the supplied family. Its rigorous core is G4-T plus G4-P; its blind-rounding step for rank ≥ 3 is
numerical.

**Precision on what the limit object is.**
- At finite N the blind partition is canonical: an argmax of canonical idempotents.
- The **asymptotic** object is a partition modulo sets of vanishing π-mass. Two canonical roundings may differ on such
  sets, and so may two "true" partitions.

**Not TRUE COMPRESSION.** No frozen residual witness is distinguished.

## G4.10 Consciousness firewall

Even this positive result derives only **dynamics → metastable macro-partition** (a coarse-graining of the state space).
It establishes no awareness, experience, subjectivity, readout, self-model, memory or observer selection.

The partition is a coarse-graining of states. It is **not** a subsystem decomposition (tensor factor or net).
