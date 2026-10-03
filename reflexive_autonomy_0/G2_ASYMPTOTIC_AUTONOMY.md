# RA0 · G2 — ASYMPTOTIC AUTONOMY (result)

**Central question.** Can a nontrivial autonomous macrostructure emerge **uniquely** in an N → ∞ limit, even though no
nontrivial exact autonomous partition exists at finite N?

**Evidence:**
- `g2/g2_asymptotic.py` + `g2/g2_asymptotic.log` (independent code path, not independent reviewer; NUMERICAL ILLUSTRATION);
- known results:
  - hydrodynamic limit of SSEP → heat equation (standard; Kipnis–Landim, *Scaling Limits of Interacting Particle Systems*,
    1999; cited from the standard literature, not re-read);
  - interchange / exclusion spectral gap = random-walk gap (Caputo–Liggett–Richthammer, JAMS 23 (2010) 831; abstract
    located);
  - spectral metastable coarse-graining (PCCA+, Deuflhard–Weber, LAA 398 (2005) 161; located).

**Firewalls in force:**
- RF-2: no ε. Every defect is a sequence in L, judged by its limit only.
- RF-3: macro criterion priced.
- RF-4: scaling priced.

## G2.1 Defects (no ε)

Three defect notions were computed, all as functions of system size:

| defect | definition | measure-free? |
|---|---|---|
| **exact-block sup** δ_N | max over blocks i, j and x, x′ ∈ B_i of \|P^τ(x, B_j) − P^τ(x′, B_j)\| (owner's G2.1 form, at time τ) | yes |
| **mean-density sup** | sup over configurations with equal class counts of \|E[count in class j at τ \| x] − E[… \| x′]\| / \|class\| (first-moment / LLN level) | yes |
| **spectral slowness** | relaxation rate of a mode divided by the intrinsic gap λ₁(N) of the same kernel; a mode is *slow* iff the ratio stays **bounded** as N → ∞ | yes (the reference scale is the gap of D itself, not a chosen ε) |

## G2.4 Controls

### Control A — independent particles (N two-state particles, product kernel)

**NUM.**
- The occupation-number partition and every single-particle partition are **exactly** lumpable (one-step defect ≈ 1e-15)
  for N = 2 … 10. A random 3-block partition has defect 0.16 – 0.43.
- The slowest nontrivial mode is N-independent (|eig| = 0.650), with **multiplicity N**: one single-particle mode per
  particle.

**Reading.**
- Exact autonomy here is **symmetry / independence**: S_N orbits give occupation numbers; product structure gives each
  particle.
- There is **no gap closing**, so the spectral-slowness criterion picks the N **individual particle** observables, not a
  macro variable.

**Grade: ASYMPTOTIC AUTONOMY = SYMMETRY / CONSERVATION RELOCATION; NONUNIQUE.** No macro emergence without spatial
structure in D.

### Control B — SSEP on a ring (L sites, N = L/2, nearest-neighbour bond swaps; reversible, uniform π derived from P)

**B1 — exact-block (TV) autonomy of the two-cell count.**

| L | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|
| sup defect, τ = 0.02 L³ | 0.269 | 0.286 | 0.311 | 0.307 | 0.330 |
| π-average defect, τ = 0.02 L³ | 0.041 | 0.039 | 0.035 | 0.032 | 0.029 |
| sup defect, τ = 0.1 L³ | 0.007 | 0.012 | 0.009 | 0.011 | 0.009 |

- **No sign that the sup defect vanishes** at fixed diffusive time (small L only, so this is not an asymptotic proof).
- Small values at τ = 0.1 L³ reflect near-global equilibration (trivial autonomy), not macro autonomy.
- The π-average decreases slowly, but it is **measure-weighted**. Using it would price the selection in the stationary
  measure.

**B2 — mean-density (LLN-level) autonomy, exact via SSEP duality, τ = 0.05 L³.**

| L | 64 | 256 | 1024 | 4096 |
|---|---|---|---|---|
| 4 contiguous cells (fixed m) | 0.0560 | 0.0559 | 0.0559 | 0.0559 |
| √L-size contiguous cells | 0.0330 | 0.0171 | 0.0087 | 0.0043 |
| 4 interleaved classes (site mod 4) | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

- **Fixed cells:** the defect does not vanish. Within-cell arrangement relaxes on the same diffusive scale as the macro
  dynamics.
- **Cells of size √L:** defect → 0, roughly ∝ L^(−1/2). The local-density field becomes asymptotically autonomous with
  **no ε**, but only under a **supplied cell-scaling** ℓ/L → 0, ℓ → ∞.
- **Interleaved classes:** zero defect, because their densities relax trivially to equality. The zero-defect criterion
  alone **accepts them too**.

**B3 — slow spectrum of the full many-particle kernel.**

| L | states | lowest many-body gaps | one-particle gaps | residual of the two slowest eigenvectors outside span{1, η_x} |
|---|---|---|---|---|
| 8 | 70 | 0.0732, 0.0732, 0.1186, … | 0.0732, 0.0732, 0.25, … | 0, 0 |
| 10 | 252 | 0.0382, 0.0382, 0.0655, … | 0.0382, 0.0382, 0.138, … | 0, 0 |
| 12 | 924 | 0.0223, 0.0223, 0.0396, … | 0.0223, 0.0223, 0.0833, … | 0, 0 |

- The **slowest** modes of the full kernel are **exactly the linear density modes** (Fourier q = ±1). The gap equals the
  one-particle gap, consistent with Caputo–Liggett–Richthammer.
- The next slow modes (residuals 0.93 – 1.0) are **nonlinear correlation modes**: slower than the second density mode, but
  not linear in occupations.
- **One-particle ratios.** λ_q / λ₁ → q² (L = 10 000: 1, 4, 9, 16, 25), while the fast-mode ratio λ_{L/2}/λ₁ ~ L²
  (1.0e7 at L = 10⁴).
- **Scale separation is intrinsic to D.** The long-wavelength density modes and their correlation modes form a slow sector
  with bounded ratios. Everything else has diverging ratios. This uses no ε and **no named empirical field**.

### Control C — symmetry / scaling-rich: two-lane ladder with inter-lane rate γ_L = L^(−β)

**NUM.** Ratio (lane-imbalance rate) / (slowest along-lane density rate):

| β | L = 50 | 200 | 800 | 3200 | trend |
|---|---|---|---|---|---|
| 0 | 127 | 2.0e3 | 3.2e4 | 5.2e5 | → ∞ (lane imbalance fast; macro = density only) |
| 1 | 2.54 | 10.1 | 40.5 | 162 | → ∞ (lane imbalance fast) |
| 2 | 0.0507 | 0.0507 | 0.0507 | 0.0507 | bounded (lane imbalance **joins** the slow sector) |
| 3 | 1.0e-3 | 2.5e-4 | 6.3e-5 | 1.6e-5 | → 0 (relative to the true gap 2γ_L, the **lane imbalance alone** is slow; density modes become fast) |

**Reading.** The emergent macro variables **flip** with the supplied family exponent β: density only (β < 2), density +
lane imbalance (β = 2), lane imbalance only (β > 2).

**Grade: SCALING RELOCATION.** The family (X_N, P_N) is supplied, and its scaling decides the answer.

## G2.2 Macro criterion (priced)

- The exact micro-partition is always autonomous (G1 PROP 4b), and trivially relaxing partitions have zero defect (B2
  interleaved). So "defect → 0" alone **does not** single out a macro structure. **ASYMPTOTIC AUTONOMY NONUNIQUE.**
- The criterion that does isolate the hydrodynamic sector is **spectral slowness with bounded ratio to the intrinsic gap**.
  It is ε-free, but it is a **supplied definition of macroscopicity**, not a theorem-forced consequence of reflexive
  closure.
- Candidate structural alternatives (blocks / |X_N| → 0; block size → ∞) are equally supplied, and B2 shows cell scaling
  matters.

## G2.3 Canonical object

- **Collection of all zero-defect asymptotically autonomous observables:**
  - includes conserved quantities, the full micro-description, and trivially relaxing class densities;
  - collapses toward "everything that is either microscopic or eventually constant";
  - **not** a canonical macro algebra.
- **Spectral slow sector** (modes with bounded rate / gap ratio):
  - **canonical given the supplied family**: defined from P_N alone, relabeling-invariant, and symmetry-covariant (the
    ±q pairs are degenerate, so the sector is canonical as a subspace, not as a preferred variable);
  - **nontrivial** (density modes);
  - **smaller than the full algebra**;
  - B3 shows it contains nonlinear correlation modes. **Closure under products, i.e. whether it is an algebra
    asymptotically generated by density modes, is OPEN**: not shown here.
- **Grade:** **CANONICAL SLOW SUBSPACE — CONDITIONAL ON THE SUPPLIED FAMILY.** The algebra status is OPEN. This is
  **not** UNIQUE NONTRIVIAL ASYMPTOTIC STRUCTURE in the selector sense: control C shows the family scaling changes it.

## G2.5 Hydrodynamic distinction

| | ordinary hydrodynamic derivation | what G2 obtained |
|---|---|---|
| supplied | space, the conserved density, the scaling, the empirical field | the system family (ring, nearest-neighbour bonds, L → ∞), and the slow-sector definition |
| derived | the limiting PDE | **the identity of the macro variables** (density Fourier modes) read off the generator's slow spectrum, **without naming the empirical field** |

This is the **one genuinely positive feature** of G2. The macro variables were not told to the construction. But they
come from **locality already encoded in D** (nearest-neighbour coupling makes the gap close as 1/L² per site-time). That
is the frozen Bridge R1 pattern: the net is re-encoded in the dynamics.

**Grade: CONDITIONAL MACROVARIABLE DERIVATION (from a supplied local family).**

## G2.6 Scaling firewall

| scaling item | status |
|---|---|
| size sequence / family | **supplied** (ring vs ladder; β changes the answer) |
| time rescaling | **derivable** from the intrinsic gap (τ ∝ 1/λ₁(N)), not chosen |
| spatial / cell rescaling | supplied in B2 (ℓ = √L); **not needed** in the spectral route |
| field normalization | not needed in the spectral route |
| boundary conditions | supplied (periodic ring) |

**Grade: SCALING-PRICED** through the family.

## G2.7 Twin-world test (identical microscopic dynamics, two macro decompositions)

**NUM (SSEP ring, contiguous 4 cells vs interleaved classes site mod 4):**

| criterion | contiguous cells | interleaved classes | distinguishes? |
|---|---|---|---|
| zero mean-defect (B2) | not autonomous with fixed m; autonomous with √L cells | autonomous (trivially, 0.0000) | **no**: accepts the interleaved twin |
| spectral slowness (D) | class differences contain the q = 1 mode: ratio **1.0** (slow) | class differences live at q = L/4: ratio 2.1e2 → 8.5e5 (fast) | **yes** |

**Reading.**
- The spectral criterion **does** distinguish these two macro decompositions *of the same dynamics*, through the locality
  of D.
- These are **not frozen selector witnesses.** The frozen Σ witnesses concern subsystem nets of the same generator. The
  slow sector is an operator-level function of the generator, so it is the same object for any two descriptions of one
  generator, and it outputs collective modes, not a subsystem net.
- **NO Σ SELECTOR.**

## G2 terminal

| outcome | applies? | basis |
|---|---|---|
| ASYMPTOTIC AUTONOMY EMPTY | no | slow density sector exists (B3, B2 √L cells) |
| **ASYMPTOTIC AUTONOMY NONUNIQUE** | **yes** | zero-defect accepts micro, conserved and trivially relaxing partitions (A, B2) |
| **ASYMPTOTIC AUTONOMY = SYMMETRY / CONSERVATION RELOCATION** | **yes (control A)** | no locality ⇒ only symmetry / independence structure; micro variables are slow |
| **CONDITIONAL MACROVARIABLE DERIVATION** | **yes (control B)** | density modes from the intrinsic slow spectrum of a supplied local family, ε-free, no named field |
| CANONICAL ASYMPTOTIC AUTONOMY ALGEBRA | **partial** | canonical slow **subspace** given the family; algebra closure OPEN |
| UNIQUE NONTRIVIAL ASYMPTOTIC STRUCTURE | **no** | control C flips the answer with β |
| A_resolution RELOCATION | avoided in the spectral route (no ε); **present** if the π-average defect or a chosen tolerance were used |
| **SCALING RELOCATION** | **yes (control C)** | the family exponent decides the macro variables |

**Positive-result checklist (owner):**

| item | met? |
|---|---|
| no finite ε | **yes** |
| no supplied observable channel | **yes** in the spectral route |
| no supplied target partition | **yes** |
| no arbitrary optimization objective | **yes** |
| nontrivial structure in the limit | **yes** |
| invariant under relabelings | **yes** |
| unique / canonical maximal algebra | **partial** (canonical subspace given the family; not unique across families) |
| discrete micro-partition excluded by a theorem-backed criterion | **no**: excluded only by the **supplied** bounded-ratio slowness definition |

**Net:** **CONDITIONAL DERIVATION OF MACROSTRUCTURE FROM A SUPPLIED SYSTEM FAMILY**, as the owner anticipated. It
distinguishes **no frozen selector witness**.

## Procedural notes

- An empty-block indexing bug in control A (N = 2 random partition) was fixed before the logged run.
- G1 log note: the G1.5 "covariance check" line is tautological ((P*)* = P) and carries no evidential weight. The G1
  orientation conclusions rest on PROP 9 and the explicit forward / backward example.
- B1 is limited to L ≤ 14 (exact state spaces up to 3432). Its "no decay" reading is an illustration, not an asymptotic
  theorem.
