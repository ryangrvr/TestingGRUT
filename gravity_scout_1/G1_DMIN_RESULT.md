# GRAVITY-SCOUT-1 · G1 — DOES GRAVITY FORCE A NONZERO MINIMUM SPLITTING DISTANCE d_min? (result)

**Target.** The split-buffer scale d that QFT-SCOUT-1 left supplied ("H_cross = 0 admissible only with a split buffer
(d, 𝒩)").

**Question (owner).** "Does gravitational consistency force d_min > 0 or otherwise constrain the admissible collar scale?"

**Evidence:**
- the theorem / source chain in §1 – §3;
- `g1/g1_dmin.py` + `g1/g1_dmin.log`, an exact product-state energy SDP. This is an independent code path, not an
  independent reviewer. Lattice runs are illustrations (GF-11).

## Verdict (two levels, per GF-7)

| level | terminal | one line |
|---|---|---|
| **d_min^incl** (split inclusion; Fewster's splitting distance) | **NO CONSTRAINT** | No tested gravitational premise produces a positive infimum for the split inclusion. Perturbative quantum gravity (QG-7) reprices the *notion* of split at every d; it does not produce a scale |
| **d_min^prep** (horizon-free preparation of a normal product state across the collar) | **CONSTRAINED-NONUNIQUE — CONDITIONAL, HEURISTIC GRADE** | Under QG-1 – QG-6 plus a localization premise, a horizon-free product preparation across a spherical collar at radius R needs **d ≥ d\*(R) = (8π α₃ N ℓ_P² R)^{1/3}**. Here α₃ ≈ 0.0068 per free massless scalar, computed below. The bound is region- and species-dependent, is **not** ℓ_P, and selects no d |

- **Planck-scale insertion** (GF-1) and the **Jacobson** cutoff reading (GF-5) are **RELOCATION**.
- The **Bousso** route (GF-3) is **CONJECTURE-CONDITIONED** and yields only a sub-Planckian bound (§3.4).
- **CONDITIONALLY SELECTED: no.** **TRUE COMPRESSION: 0.** **Payoff: none.**

## 1. Flat / AQFT control (GF-9)

### 1.1 Inclusion level

**d_min^QFT = 0 at the level of the split-distance infimum, under the stated split-property assumptions** (owner-stated
control; Fewster arXiv:1501.02682, owner-scoped):
- The sharp boundary d = 0 admits no normal product state (QFT-SCOUT-1 G1-P1).
- Every d > 0 can admit a split inclusion.
- The infimum 0 is **not attained**.

### 1.2 Preparation level

This is new here: **GS-L1, an exact lattice statement.**

**Statement.** For a free scalar on a lattice, the minimum energy above the vacuum of **any** state (mixed, non-Gaussian,
displaced) whose restriction to 𝒜(A) ∨ 𝒜(B) is a product equals the value of the SDP:

    min ½ tr P + ½ tr K X   s.t.   [[X, ½I], [½I, P]] ⪰ 0,   X_AB = P_AB = 0.

**Proof sketch** (GRAVITY-SCOUT PROPOSITION, not externally reviewed):
- A product state has vanishing centered A–B covariance.
- Its Gaussianization has the same second moments and energy, and it is a product state on A ∪ B.
- First moments only add energy.
- Averaging with the time-reversed state (p → −p) removes the x–p blocks. This preserves feasibility, energy and the
  product constraint.
- [[X, ½I], [½I, P]] ⪰ 0 is unitarily equivalent to the uncertainty condition Γ + iΩ/2 ⪰ 0 for Γ = X ⊕ P. ∎

**Run (Part A).** 1+1 dimensions, collar of physical width D = 1, mass M = 0.5, lattice spacing a reduced:

| a | E_prod(D = 1) | E_prod(sharp, D = 0) |
|---|---|---|
| 1/2 | 0.05277 | 0.17973 |
| 1/3 | 0.05404 | 0.30848 |
| 1/4 | 0.05482 | 0.44058 (≈ 0.11 / a, divergent) |

- **Buffered:** the energy converges and stays finite.
- **Sharp:** the energy diverges as a → 0, consistent with G1-P1.
- The minimizer's excess energy lies **≥ 99%** within 3D of the collar.

**Flat preparation-level control: d_min^prep,flat = 0.** Every D > 0 costs finite energy. Without gravity, nothing
converts that energy into an inadmissibility.

## 2. Inclusion level under gravity: NO CONSTRAINT

| premise layer | what it does to d(S) | grade | result |
|---|---|---|---|
| **QFT on a fixed curved background** (QG-4 only; G not dynamical) | Fewster's locally covariant analysis transports split structure across globally hyperbolic spacetimes. Under timeslice plus local quasi-equivalence, sufficiently strong distal splitting gives d = 0 for balls. Curvature enters the collar geometry but supplies **no universal minimum** (GF-6) | owner-scoped source | **NO CONSTRAINT** |
| **Semiclassical backreaction** (QG-5) | Acts on **states** through ⟨T_μν⟩_ω. For each fixed solution the net is again a QFT on a curved background, so the inclusion-level infimum is unchanged | proposition (sketch) | **NO CONSTRAINT** at the inclusion level (see §3 for states) |
| **Perturbative QG, gravitational dressing** (QG-7) | Gauge invariance already conflicts with tensor factorization and with a net of commuting local subalgebras **at leading order** in G. Donnelly–Giddings construct a leading-order **"gravitational splitting"**. The obstruction and its repair are **scale-free**: they concern gravitational charges and dressing, **not a minimum collar width** | SOURCE LOCATED, TEXT NOT RE-READ (Donnelly & Giddings, PRD 98, 086006 (2018), arXiv:1805.11095; abstract only) | **no d_min.** It is a **class modification at every d**: the split notion itself is repriced |
| **Crossed product / observer dressing** | The type III₁ algebras of horizon regions become type II∞ (Witten, arXiv:2112.12828, JHEP 10 (2022) 008) or type II₁ (de Sitter static patch, Chandrasekaran–Longo–Penington–Witten, arXiv:2206.10780). Entropy becomes defined up to a state-independent constant. **Type II is still non-type-I, so G1-P1's sharp obstruction persists.** No collar scale appears, and the finite entropy comes without a literal cutoff (consistent with GF-2) | SOURCE LOCATED, TEXT NOT RE-READ (abstracts) | **NO CONSTRAINT** on d; G1-P1 survives |

## 3. Preparation level under gravity

### 3.1 The chain (GS-P1, conditional)

**Premises:**
- **QG-1 – QG-5**, plus **QG-6** in its spherical (Misner–Sharp) form.
- **F:** N free massless scalars, with the flat local approximation d ≪ R and d well above the semiclassical cutoff.
- **L (localization; HEURISTIC):** the excess energy of a product preparation lies within r ≤ R + O(d).
  - This is **verified for the energy minimizer**: 97–99.7% of the excess lies within 3d in every 3+1 run, and ≥ 99% in
    1+1.
  - It is **not proved for all product states.** Negative local energy densities are not excluded pointwise.

**Steps:**
1. **Flat energy cost** (GS-L1, plus the transverse-mode decomposition; Part B). For a planar collar of width d, the minimum
   product-state energy per unit area is

       E/A = α₃ N ħc / d³,     α₃ = 0.00709, 0.00689, 0.00685, 0.00679   (d = 2, 3, 4, 6 lattice sites).

   - The function f(u) = d·e(d, u/d) is d-independent to a few percent where the integrand is concentrated (u ≲ 2).
     This is the continuum scaling check. Deviations appear at large u for d = 2 (lattice) and at the smallest u for
     d = 6 (finite box).
   - Conservative read: **α₃ = 0.0068 ± 0.0003** per scalar. The spread covers lattice and finite-box effects.
   - The bound is additive over independent fields.
   - Rotation averaging preserves feasibility and energy, so a spherically symmetric minimizer exists.
2. **Spherical collar at R ≫ d:** E ≥ 4πR² α₃ N ħc / d³ · (1 + O(d/R)).
3. **No trapped surface** (QG-6, spherical): 2G m(r) / (c² r) < 1. With L, m(R + O(d)) ≳ E/c².
4. **Result:**

       d³ > 8π α₃ N ℓ_P² R      ⇔      d ≥ d*(R) = (8π α₃ N ℓ_P² R)^{1/3} ≈ (0.171 N ℓ_P² R)^{1/3}.

**Sizes** (N = 1):

| R | d* | d* / ℓ_P |
|---|---|---|
| 10⁻¹⁵ m | 3.6 × 10⁻²⁹ m | 2 × 10⁶ |
| 1 m | 3.6 × 10⁻²⁴ m | 2 × 10¹¹ |
| 6.4 × 10⁶ m | 6.6 × 10⁻²² m | 4 × 10¹³ |
| 1.4 × 10²⁶ m | 1.8 × 10⁻¹⁵ m | 1 × 10²⁰ |

**Self-consistency.** d* ≫ ℓ_P and d* ≪ R both hold whenever R ≫ ℓ_P, so the semiclassical regime is respected.

### 3.2 Firewall audit of GS-P1

| check | result |
|---|---|
| GF-1 (Planck insertion) | **passes.** d* is not set by dimensional analysis: it follows from a consistency condition (no trapped surface), and its form ℓ_P^{2/3} R^{1/3} is not ℓ_P |
| Supplied scales | ħ (QP-5), G (QG-1), c (QG-4), R (the localization choice), N (field content). Every scale in d* is supplied, so this is **not compression** (GF-8) |
| G → 0 control (GF-9) | d* → 0: the flat result is recovered. **The constraint is entirely gravitational** |
| QNEC (GF-4) | not used |
| Susskind–Uglum (GF-2) | not used. The bound uses **energy**, not the entropy divergence |
| Level (GF-7) | a statement about **preparing H_cross = 0 states**, not about the split inclusion. Normal product states still **exist** for every d > 0. Below d* they cannot be prepared horizon-free (conditional on L) |

### 3.3 Grade and terminal

**Terminal: CONSTRAINED-NONUNIQUE (preparation level).** Part of the collar family, d < d*(R), is ruled out for horizon-free
product preparation. Every d ≥ d* remains admissible, so nothing is selected.

**Grade:** CONDITIONAL / HEURISTIC.
- Step 1 is exact on the lattice for free scalars. Its continuum form is a numerical illustration.
- Steps 2 – 3 use the planar approximation and premise L.
- The general-shape version needs the hoop **conjecture**.

**Likely not novel.** It is the gravitational form of the familiar "disentangling a boundary costs divergent energy"
argument; compare the "energetic curtain" of Braunstein, Pirandola & Życzkowski, PRL 110, 101301 (2013) (SOURCE LOCATED,
TEXT NOT RE-READ). The new elements here are the exact lattice minimum, α₃ and the level separation. No literature search
for prior statements of d* was exhaustive.

### 3.4 The other gravity–entropy routes

| route | output | terminal |
|---|---|---|
| **Dimensional analysis** d ~ √(ħG/c³) | d = ℓ_P | **RELOCATION / SCALE INSERTION** (GF-1) |
| **Jacobson entanglement equilibrium** | Identifying the area-law coefficient with 1/(4ħG) and reading it as a cutoff gives ε ~ √N ℓ_P. But ε is a regulator, not a collar. Jacobson assumes the small-ball and equilibrium premises (GF-5), and Susskind–Uglum renormalization blocks the literal-cutoff reading (GF-2) | **RELOCATION** |
| **Bousso-type entropy bound on the collar** | The vacuum mutual information across the collar is I/A = κ_I N/d², with κ_I ≈ 0.005 (Part B: 0.0054 → 0.0048; IR / box-sensitive). Requiring I/2 ≤ A/(4ℓ_P²) gives d ≥ √(2κ_I N) ℓ_P ≈ **0.1 √N ℓ_P**. For N ≲ 100 this is **sub-Planckian**, i.e. outside the semiclassical domain where it could be stated | **CONJECTURE-CONDITIONED** (GF-3); at best Planck-scale RELOCATION; **no admissible constraint** in the semiclassical regime |
| **QNEC** | no gravity | control only (GF-4) |

**Read-off.** The entropy routes return Planck-scale numbers (RELOCATION). The only route that produced a non-trivial
constraint uses **energy plus collapse**, and it constrains states, not inclusions.

## 4. Payoff (seven-criterion bar)

d*(R):
- is not forced at theorem grade (premise L; hoop / spherical);
- is not residual-input invariant (it depends on R and N);
- is ordinary semiclassical GR + QFT physics, so it is **not distinctive** versus the comparison class (criterion 5);
- is not observable at any accessible scale (≤ 10⁻¹⁵ m even for Hubble-scale R).

**NO DISTINCTIVE PAYOFF. ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS (preserved).**

## 5. Residual update

| component | before | after G1 |
|---|---|---|
| d (split collar), inclusion level | supplied; d_min^QFT = 0 | **unchanged: NO CONSTRAINT** |
| (H_cross = 0, d), preparation level | admissible for every d > 0 | **CONSTRAINED-NONUNIQUE (conditional):** horizon-free product preparation needs d ≥ (8πα₃Nℓ_P²R)^{1/3} |
| split notion itself | type-I interpolation 𝒩 | under QG-7, repriced at every d (gravitational splitting modulo dressing / charges). Source located; **feeds G2 (𝒩)** |
| G1-P1 sharp obstruction | non-type-I | survives the III₁ → II crossed-product change |

Every result is **auxiliary to canonical GRUT** (QG-1 – QG-7 and QP-5 are supplied).

## 6. Procedural notes

- The SDP solver (CLARABEL) reports `optimal_inaccurate` for a few modes. The tail modes sit at the solver floor of about
  10⁻⁶ and are dropped below 10⁻⁵, replaced by an exponential tail (contributing ≤ 2 × 10⁻⁵ to α₃).
- A finer a = 1/8 Part-A run was dropped for cost.
- The Part-B box (L = 20 sites per side) truncates the infrared at the smallest u for d = 6. This affects κ_I by about 10%
  and α₃ by about 1%.
