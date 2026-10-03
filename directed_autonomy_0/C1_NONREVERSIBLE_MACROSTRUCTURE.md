# DA0 · C1 — NON-REVERSIBLE ENDOGENOUS MACROSTRUCTURE (result)

**Central question.** Can the RA0 endogenous macro-algebra / partition construction survive the removal of detailed
balance, and produce a canonically directed coarse structure without supplying:
- a time orientation (beyond the priced forward-kernel convention);
- a subsystem boundary;
- an observable channel?

**Status of results:**
- **Proofs** are INTERNALLY PROVED / NOT EXTERNALLY REVIEWED. "ASSEMBLED" means built from standard tools (Kato's
  contour-integral perturbation theory; the pseudospectral / Bauer–Fike inclusion), with constants derived here.
- **Numerics:** `c1/c1_nonreversible.py` + `c1/c1_nonreversible.log` (independent code path, not independent reviewer;
  NUMERICAL ILLUSTRATION).
- **Preregistration:** DA0 CHARTER REPAIR 01 (DC1-01 … 04) was applied before any computation.

## C1.0 Canonical objects used (DC1-01)

**From (G, π) only:**
- the stationary law π (irreducibility priced);
- the π-adjoint G*(x, y) = π_y G(y, x)/π_x;
- the current J_xy = π_x G_xy − π_y G_yx;
- the entropy production EP = ½ Σ J log[(π_x G_xy)/(π_y G_yx)];
- the cycle-space **dimension** of the support graph;
- Riesz projectors of conjugation-closed spectral sets.

**Mutual support.** Every family has mutually supported edges, so the finite EP formula applies. The code asserts this.

**What is not used:**
- **No cycle basis or fundamental-cycle decomposition** is used anywhere.
- A cycle is called canonical **only** when the cycle space is one-dimensional. Then every basis equals ± the unique
  cycle, and the sign of the current fixes its orientation, so nothing is chosen (PROP C1-C).

## C1.1 Theory

### LEMMA C1-0 (KNOWN / REDERIVED)

- J is divergence-free.
- J(G*) = −J(G), and σ(G*) is the complex conjugate of σ(G), so the decay rates are equal.
- EP ≥ 0, with EP = 0 iff detailed balance holds.
- Coarse-graining: for any partition Π, the lumped π-weighted generator
  Q_Π(a, b) = Σ_{x∈a, y∈b} π_x G_xy / π(a) has π(Π) as its stationary law.
- Its current is the coarse-grained current, J_Q(a, b) = Σ_{x∈a, y∈b} J_xy.
- EP_Q ≤ EP, by the log-sum inequality.

### LEMMA C1-1 (PROVED HERE, elementary)

**Setting.** For a hidden proof partition 𝔅 = {B_1, …, B_K}, split G = G₀ + G₁:
- G₀ keeps the intra-block rates (with its diagonal adjusted), and each block is irreducible under G₀;
- G₁ is the inter-block part.

**Statement.** Then:
- 0 is a semisimple eigenvalue of G₀ with multiplicity K.
- Its Riesz projector is E₀ = Σ_b 1_{B_b} ⟨π_b^{(0)}, ·⟩, where π_b^{(0)} is the stationary law of block b under G₀.
- **Ran E₀ = A_𝔅**, the block-indicator algebra, **with or without reversibility**.
- E₀ is π-orthogonal iff π_b^{(0)} ∝ π|_{B_b}.

### THEOREM C1-T (non-reversible analogue of G4-T; ASSEMBLED, internally proved)

**Hypotheses** (fixed K; priced). All norms are in L²(π).

| | hypothesis | meaning |
|---|---|---|
| N1 | e_N := ‖E₀‖ ≤ e* | block-stationary laws are not wildly non-π |
| N2 | there is g_N > 0 such that a_N := sup_{Re z ≥ −g_N} ‖(z − G₀)^{−1}(I − E₀)‖ satisfies a_N g_N ≤ C_a | **non-normal amplification of the fast intra-block dynamics** is bounded. For reversible G₀ with g = half the intra gap, a·g = 1 |
| N3 | η_N := ‖G₁‖ with η_N/g_N → 0 | blocks are nearly uncoupled |

**Definitions.** ℓ = ηe/(1 − ηa), ρ = √(ℓg) and δ = η(e/ρ + a).

**Conclusions.** For all N with ηa < 1, ℓ < g and δ < 1:
1. **Spectral split.** σ(G) ⊂ {|z| ≤ ℓ} ∪ {Re z ≤ −g}.
2. **Exactly K slow eigenvalues.** They lie in |z| ≤ ℓ and form a conjugation-closed set. Their Riesz projector satisfies
   ‖E − E₀‖ ≤ 2eaη + (e + ρa)·δ²/(1 − δ) = O(η/g).
3. **Certified diverging decay-rate cut at rank K:** r_K/r_{K−1} ≥ g/ℓ → ∞.
4. **No non-normal instability:** ‖E‖ ≤ e + o(1).
5. **Subspace convergence:** gap(V_N, A_𝔅) ≤ ‖E − E₀‖ → 0, where V_N = Ran E_N.

**Proof.**
- **(1)** Since G₀E₀ = 0, (z − G₀)^{−1} = E₀/z + R_⊥(z), with R_⊥ analytic and ‖R_⊥‖ ≤ a on Re z ≥ −g. If λ ∈ σ(G),
  λ ≠ 0 and Re λ > −g, the pseudospectral inclusion (Bauer–Fike type) gives 1/η ≤ ‖(λ − G₀)^{−1}‖ ≤ e/|λ| + a. Hence
  |λ| ≤ ℓ.
- **(2) Contour.** The circle |z| = ρ lies in Re z > −g (ρ < g), and separates the two parts of the spectrum (ℓ < ρ).
  On it, ‖(z − G₀)^{−1}‖ ≤ R₀ := e/ρ + a, and ‖(z − G)^{−1}‖ ≤ R₀/(1 − δ).
- **(2) Expansion.**
  E − E₀ = (2πi)^{−1}∮ (z−G₀)^{−1}G₁(z−G₀)^{−1} dz + (2πi)^{−1}∮ (z−G)^{−1}G₁(z−G₀)^{−1}G₁(z−G₀)^{−1} dz.
- **(2) First term.** Substitute the decomposition. The E₀G₁E₀/z² term and the R_⊥G₁R_⊥ term integrate to 0. The residue
  term is E₀G₁R_⊥(0) + R_⊥(0)G₁E₀, of norm ≤ 2eaη.
- **(2) Second term.** Norm ≤ ρ·R₀/(1 − δ)·η²R₀² = (e + ρa)δ²/(1 − δ).
- **(2) Rank.** ‖E − E₀‖ < 1 forces rank E = rank E₀ = K (Kato).
- **Asymptotics.** δ ≈ √(ηe/g) + C_a·η/g → 0, so ‖E − E₀‖ = O(η/g).
- **(3)** Fast rates are ≥ g; slow rates are ≤ |λ| ≤ ℓ.
- **(5)** For unit x in one range, the distance to the other range is ≤ ‖(E − E₀)x‖. The dimensions are equal. ∎

### COROLLARY C1-R (blind partition recovery transfers; internally proved)

Under C1-T, add p_min ≥ p_* and s_N²·C_{V_N} → 0, where s_N = gap(V_N, A_𝔅). Then **Theorem G5 applies verbatim**:
- R(V_N) has exactly 2^K idempotents and K primitive ones once the rigorous condition holds;
- misclassified π-mass m_N ≤ (256/9)Kε_N² + 8s_N² = O((η/g)²);
- Corollary G5-U (asymptotic uniqueness).

**Why it transfers.** G5's proof uses only s → 0, the algebra structure of A_𝔅 and the probability measure π. It never
uses reversibility, which entered RA0 only through Davis–Kahan in G4-T. C1-T replaces that step.

**Blind / certified split (RA3-03, inherited).** Recovery is blind. **Certification is not:** checking N1 – N3 is a
family-level proof that uses the hidden 𝔅.

### PROP C1-D (macro current recovery; internally proved)

**Setting.** Let 𝔅 be the hidden partition, J_𝔅 its coarse-grained current and Π_N = R(V_N) the blind partition.

**Bound.** After block matching, every entry satisfies |J_Π − J_𝔅| ≤ 4·q_max·m_N, where q_max is the maximum exit rate.
- Every changed pair involves a misassigned state.
- Σ_{x∈S} Σ_y |J_xy| ≤ 2Σ_{x∈S} π_x q_x (stationarity), and this is counted twice.

**Recovery.** Add **N4 (non-degenerate macro drive):** each nonzero entry of J_𝔅 is ≥ c_J·η_N, or J_𝔅 ≡ 0. Then
4q_max·m_N/(c_J η_N) = O(η/g²) → 0. So the **support and sign pattern of the macro current**, and the macro EP, are
recovered blindly.

**Covariance.** G* satisfies N1 – N4 with the same constants:
- ‖G₁*‖ = ‖G₁‖;
- E₀(G*) = E₀(G)*;
- the adjoint resolvent has the same norm at conjugate points;
- the blocks are the same.

So Π(G*) and Π(G) converge in mass to the same 𝔅, and J_Q(G*) = −J_Q(G). **The macro current and the affinities reverse;
the undirected partition is invariant.** ∎

### PROP C1-A (A-pricing; PROVED HERE)

Take two generators with the same π and the same symmetric part S, so they differ only in A.
- J = 2F_a, where F_a is the antisymmetric part of the flux π_x G_xy. So **J_𝔅 depends on A alone**.
- If A carries no inter-block flux, J_𝔅 ≡ 0. The macro-generator then satisfies **detailed balance at the macro level**,
  even though micro EP > 0.
- So macro directedness is a property of the supplied A: **A-PRICED, not TRUE COMPRESSION.** ∎

### PROP C1-C (canonical macro cycle; KNOWN graph theory, applied)

**Canonical cycle.** If the support graph of Q_Π has a one-dimensional cycle space, the unique cycle, oriented by the sign
of J_Q, is canonical, and J_Q is a multiple of it.

**No canonical cycle.** If the dimension is ≥ 2 (e.g. a complete macro graph with K ≥ 4), J_Q is canonical but no cycle
is. The terminal is then **FLUX STRUCTURE ONLY — NO CANONICAL CYCLE**.

**Caveat for this run.** In every family run here, the macro cycle-space dimension is 1 **by construction**: K = 3
complete, or a K = 4 ring. That is a family property, not a derived one.

## C1.2 Results (numerical illustration)

### K1 controls

| control | result | reading |
|---|---|---|
| two-state bit | EP = 0 for any rates | every two-state chain satisfies detailed balance, so **no directedness is possible** |
| driven 3-state loop (thermostat / clock) | EP = 2.07; cycle space dim 1; J = 0.300 on every edge; eigenvalues −1.65 ± 0.78i, 0 | **a canonical directed cycle with no metastable cut.** Directed cycles per se are K1-level |

### R1 matched pair (DC1-03)

R1a and R1b have:
- the same π (difference ≤ 2e-16);
- the same S (difference ≤ 4e-15);
- non-collinear A: the π-weighted cosine between them is ~1e-15, i.e. they are orthogonal.

### Partition families

The blind side reads (G, π) and the certified rank. The audit uses the hidden labels, afterwards.

| family (M = 80 unless noted) | micro EP | rank-K cut | ‖E‖_π / K(Γ) / cond | R(V): idempotents / primitive | misclassified (audit) | macro EP_Q | macro current (blind) | audit of orientation | G → G* |
|---|---|---|---|---|---|---|---|---|---|
| **R0** F1 reversible (K4) | 1e-29 | 6.44 @ 3 | 1.000 / 1.65 / 1.7 | 8 / 3 | 0 | 1e-30 | **zero** | — | same Π; J_Q* = −J_Q (both 0) |
| **R1a** inter-block circulation | 8.9e-4 | 6.45 @ 3 | 1.000 / 1.65 / 1.7 | 8 / 3 | 0 | 6.7e-4 | canonical macro 3-cycle, circulation 1.0e-3 | = hidden 0→1→2 (all M) | same Π; max\|J_Q* + J_Q\| 1.5e-16 |
| **R1b** intra-block circulation | 1.3e-2 | 6.44 @ 3 | 1.000 / 1.65 / 12 | 8 / 3 | 0 | 2e-32 | **zero**: macro detailed balance despite micro EP > 0 | — | same Π |
| **R2** 4 basins on a driven ring | 2.4e-2 | 15.0 @ 4 (2.0, 4.1, 7.7, 15.0 for M = 10 … 80) | 1.000 / 1.35 / 5.6 | 16 / 4 | 0 | 2.1e-2 | canonical macro 4-cycle, circulation 3.2e-3 | = hidden 0→1→2→3 (all M) | same Π; 1.3e-16 |
| **R4** non-normal stress (3 drifting rings) | 0.49 | 126 @ 3 | 1.000 / 1.02 / 20 | 8 / 3 | 0 | 2e-17 | incidental tiny current (2.5e-12), orientation varies with M | matches hidden J_𝔅 | same Π |

**Proof-side quantities for C1-T** (hidden labels; after the blind side):

| family | η/g (M = 20 → 80) | a·g | ‖E − E₀‖ vs bound | reading |
|---|---|---|---|---|
| R0, R1a, R1b | 1.86 → 0.46 | 1.00 | 0.033 vs bound 41.9 at M = 80; bound infinite at M ≤ 40 | **vacuous at tested sizes**; η/g ∝ 1/M, so the hypotheses hold asymptotically |
| R2 | 1.39 → 0.14 (M = 10 → 80) | 1.00 | 0.0070 ≤ 0.96 at M = 80 | bound non-vacuous from M = 40 |
| R4 | 0.14 → 0.016 (M = 10 → 80) | 0.55 → 0.03 | 3.0e-6 ≤ 0.019 at M = 80 | bound non-vacuous at all M |

**The bounds are valid wherever finite, and loose.** The gap(V, A) measured on the blind side equals the audit value.

### R3 — biased site-modulated ring (no metastability)

- **No cut.** The top decay-rate ratio stays bounded at → 4 (3.94, 3.98, 3.995, 3.999 for L = 25 … 200), like G3's
  diffusive case.
- **Canonical cycle.** The cycle space has dimension 1, so the ring is the canonical cycle.
- **Uniform current.** J is identical on every edge, as it must be for a divergence-free current on a cycle graph
  (2.75e-3 at L = 200).
- **Rotating slow mode.** The slowest pair is −0.0007 ± 0.0173i at L = 200.

### DC1-02 (non-normal instability)

**Not triggered in any family.**
- ‖E‖_π ≤ 1.008, K(Γ) ≤ 3.8 and cond_π(eigenvectors) ≤ 37 across all partition families.
- R4, built to stress non-normality, has a·g **decreasing** (0.55 → 0.03).
- A diagnostic probe of driven rings with a slow defect bond gave cond up to 1.2e2 and ‖E(rank 3)‖ up to 7.6, both
  **decreasing** with L.
- No family producing SPECTRAL SEPARATION PRESENT — NONNORMAL RECOVERY UNSTABLE was constructed. That terminal stays
  untested, and is not claimed to be impossible.

## C1.3 Per-family terminals

| family | C1-P (partition) | C1-D (directed macro-generator) | terminal |
|---|---|---|---|
| R0 | yes (= RA0 G5) | **no** (zero macro current) | K4 baseline reproduced; directed content: none |
| R1a | yes | **yes**: canonical macro 3-cycle, covariant | **NONREVERSIBLE PARTITION DERIVATION**, graded CONDITIONAL DERIVATION OF DIRECTED MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D |
| R1b | yes | **no**: macro detailed balance, micro EP > 0 | NONREVERSIBLE PARTITION DERIVATION with **C1-D negative**. The directedness stays microscopic |
| R2 | yes | **yes**: canonical macro 4-cycle; complex slow eigenvalues (macro rotation), i.e. a metastable-set / cycle **mixture** in the Conrad–Weber–Schütte sense | **NONREVERSIBLE PARTITION DERIVATION**, graded CONDITIONAL DERIVATION OF DIRECTED MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D |
| R4 | yes | technically yes, but negligible (current ~1e-12; orientation set by incidental asymmetry) | partition derived; directed content **not interpreted** |
| R3 | **no** (no diverging cut) | — | **CYCLE / FLUX STRUCTURE ONLY — NO PARTITION**, with a canonical cycle (cycle space dim 1). **K1-equivalent** |

**Overall C1 terminal: NONREVERSIBLE PARTITION DERIVATION** (theorem-grade under N1 – N3 plus G5 regularity; macro-current
recovery under N4), graded **CONDITIONAL DERIVATION OF DIRECTED MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D**. **Directed
macrostructure appears iff the supplied A carries net inter-block circulation** (PROP C1-A; R1a vs R1b).

## C1.4 Comparator table (C0 §2.3)

| control | relation to the C1 positive (R1a, R2) |
|---|---|
| **K1** metastable bit / thermostat | **Coincides in kind at the macro level.** The derived macro-dynamics is itself a small driven cycle, the same type of object as the driven 3-state loop. What differs is that the loop's states are **derived, not supplied**. A directed cycle alone is not novel |
| **K2** ε-machine / causal states | **Not exceeded.** The macro process is asymptotically a K-state Markov chain whose rows differ generically, so its causal states are, asymptotically, the macrostates themselves. C1 **supplies the endogenous process and channel** that K2 takes as given (permitted by C0-P3), but adds nothing beyond predictive sufficiency at the macro level |
| **K3** IIT 4.0 | **Structurally different.** No units, candidate subsets, grain or maximisation are used; there is no integration measure |
| **K4** RA0 reversible partition | **Extends it.** Partition recovery now holds without detailed balance (Kato in place of Davis–Kahan), and a directed macro current is recovered, with covariance under G → G* |

## C1.5 Information accounting (DC1-04)

**Derived:**
- the macro-partition, blindly, in the certified nearly-uncoupled non-reversible class;
- the macro-current topology and sign pattern (under N4);
- the canonical macro cycle where the macro cycle space is one-dimensional;
- covariance under time reversal.

**Supplied:**
- the microscopic non-reversible generator D, **including A**, which fixes whether any macro current exists;
- its forward-time convention (one bit);
- irreducibility;
- the family / scaling, and the nearly-uncoupled certificate (N3);
- the non-normal regularity N1 and N2;
- p_min and s²C_V;
- macro-drive non-degeneracy N4;
- the one-dimensional macro cycle space (by construction).

**Not earned:**
- orientation from orientation-free laws;
- inside / outside;
- subsystem factorisation;
- consciousness;
- **TRUE COMPRESSION** relative to frozen GRUT: 0.

## C1.6 Hard-stop answers

See `DA0_ZOOM_OUT_01.md`.
