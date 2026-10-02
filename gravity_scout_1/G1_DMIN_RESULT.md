# GRAVITY-SCOUT-1 · G1 — DOES GRAVITY FORCE A NONZERO MINIMUM SPLITTING DISTANCE d_min? (result)

> **Repaired by GRAVITY REPAIR 01 (GR1-01 … GR1-10; `GRAVITY_CORRECTION_LEDGER.md`).** Accepted provisionally by owner
> ruling after that repair. Script and log outputs are kept as emitted; where wording differs, the correction ledger
> takes precedence.

**Target.** The split-buffer scale d that QFT-SCOUT-1 left supplied ("H_cross = 0 admissible only with a split buffer
(d, 𝒩)").

**Question (owner).** "Does gravitational consistency force d_min > 0 or otherwise constrain the admissible collar scale?"

**Evidence:**
- the theorem / source chain in §1 – §3;
- `g1/g1_dmin.py` + `g1/g1_dmin.log`, an exact finite-lattice product-state energy SDP. This is an independent code path,
  not an independent reviewer. Lattice runs are illustrations (GF-11).

## Verdict

### Inclusion level, by regime [GR1-01, GR1-07]

| regime | terminal |
|---|---|
| **fixed-background locally covariant QFT** | **d_min^incl = 0** at the level of the split-distance infimum, under the stated Fewster / split-property assumptions |
| **semiclassical backreaction** | **no tested semiclassical-backreaction argument derives a positive universal collar scale at the inclusion level.** State-dependent backreaction is *not* simply identical to a fixed-background locally covariant QFT theorem, so no stronger statement is made |
| **perturbative dynamical gravity** (Donnelly–Giddings, O(κ)) | **NO POSITIVE COLLAR SCALE DERIVED; SUBSYSTEM NOTION REPLACED / REPRICED.** Ordinary LQFT localization is altered by gravitational gauge invariance. A first-order gravitational splitting is constructed from an arbitrary extended neighbourhood U_ε, with no positive minimum ε at O(κ) |
| **fine-grained / full gravity** (Raju, scoped examples) | **ORDINARY d_min^incl MAY BE NOT APPLICABLE / BLOCKED IN THESE GRAVITATIONAL SETTINGS.** The ordinary split property fails in the specific asymptotically flat / AdS settings Raju treats, via holography of information. This is **not** generalized to every possible theory of quantum gravity: the paper scopes its examples |

**Overall inclusion conclusion.** Gravity did **not** produce a positive universal minimum splitting length. In dynamical
gravity, the stronger issue is that the **ordinary split structure may be replaced or may fail**, rather than acquire a
finite d_min.

### Preparation level [GR1-02, GR1-03]

| object | terminal |
|---|---|
| **(H_cross = 0, d, R, N, state class)**, horizon-free product preparation | **PREPARATION-CLASS CONSTRAINT: CONSTRAINED-NONUNIQUE — CONDITIONAL ON PREMISE L_all — HEURISTIC** |

- Under QG-1 – QG-6 and L_all, a horizon-free product preparation across a spherical collar at radius R needs
  **d ≳ d\*(R)**.
- Robust candidate scaling: **d\* ~ (N ℓ_P² R)^{1/3}**. The model-dependent prefactor is (8π α₃)^{1/3}, with
  α₃ ≈ 0.0068 per free massless scalar [GR1-05].
- **d ≥ d\* is not excluded by this particular collapse test** [GR1-03]. It is a necessary condition under its premises.
  It does not prove that a physically preparable state exists for d ≥ d\*, that it is stable, or that no other
  gravitational obstruction applies.
- The constraint sits on the **H_cross / H_epoch (preparation) side** of the residual. It does **not** narrow Fewster's
  split-distance primitive, and it selects no d [GR1-02].

### Other routes

- **Planck-scale insertion** (GF-1) and the **Jacobson** cutoff reading (GF-5) are **RELOCATION**.
- The **Bousso** route is a **HEURISTIC APPLICATION OF A CONJECTURE** and supplies no G1 evidence [GR1-08].
- **CONDITIONALLY SELECTED: no.** **TRUE COMPRESSION: 0.** **Payoff: none.**

## 1. Flat / AQFT control (GF-9)

### 1.1 Inclusion level

**d_min^QFT = 0 at the level of the split-distance infimum, under the stated split-property assumptions** (owner-stated
control; Fewster arXiv:1501.02682, owner-scoped):
- The sharp boundary d = 0 admits no normal product state (QFT-SCOUT-1 G1-P1).
- Every d > 0 can admit a split inclusion.
- The infimum 0 is not attained.

### 1.2 Preparation level

This is new here: **GS-L1, an exact finite-lattice statement at its stated scope** [GR1-04].

**Statement.** For a free scalar on a finite lattice with K > 0, the minimum energy above the vacuum, over all states whose
restriction to A ∪ B is a product, equals the value of the SDP:

    min ½ tr P + ½ tr K X   s.t.   [[X, ½I], [½I, P]] ⪰ 0,   X_AB = P_AB = 0.

**Proof** (covariance argument; GRAVITY-SCOUT PROPOSITION, not externally reviewed):
1. Start with any state whose AB restriction is a product.
2. Its centered A–B covariance blocks vanish.
3. Gaussianization with the same first and second moments preserves the quadratic energy. It gives a Gaussian state whose
   reduced AB covariance is block diagonal, so its AB Gaussian restriction is a product.
4. Remove the first moments by displacement. Since K > 0, this cannot increase the energy.
5. At the **covariance** level, average Γ with its time-reversed covariance TΓT.
   - The uncertainty cone {Γ : Γ + iΩ/2 ⪰ 0} is convex and T-invariant.
   - The energy is unchanged.
   - The AB cross blocks stay zero.
   - The x–p block vanishes.

   (No mixture of *states* is formed. A convex mixture of product states can create classical A–B correlations; a
   covariance average with zero cross blocks does not.)
6. Realize the averaged covariance as a Gaussian state.
7. With Γ = X ⊕ P, the condition Γ + iΩ/2 ⪰ 0 is unitarily equivalent to [[X, ½I], [½I, P]] ⪰ 0. ∎

**Run (Part A).** 1+1 dimensions, collar of physical width D = 1, mass M = 0.5, lattice spacing a reduced:

| a | E_prod(D = 1) | E_prod(sharp, D = 0) |
|---|---|---|
| 1/2 | 0.05277 | 0.17973 |
| 1/3 | 0.05404 | 0.30848 |
| 1/4 | 0.05482 | 0.44058 (≈ 0.11/a, divergent) |

- **Buffered:** the energy converges and stays finite.
- **Sharp:** the energy diverges as a → 0, consistent with G1-P1.
- The minimizer's excess energy lies ≥ 99% within 3D of the collar.

**Flat preparation-level control: d_min^prep,flat = 0.**

## 2. Inclusion level under gravity

| regime | content | grade | result |
|---|---|---|---|
| **Fixed-background locally covariant QFT** (QG-4 only) | Fewster's analysis transports split structure. Under timeslice plus local quasi-equivalence, sufficiently strong distal splitting gives d = 0 for balls. Curvature supplies no universal minimum (GF-6) | owner-scoped | **d_min^incl = 0** (under those assumptions) |
| **Semiclassical backreaction** (QG-5) [GR1-07] | Acts on states through ⟨T_μν⟩_ω. **No tested argument derives a positive universal inclusion-level collar scale.** The earlier claim that the infimum is unchanged for each fixed semiclassical solution **is withdrawn**: it does not follow from Fewster | — | **no positive scale derived** |
| **Perturbative QG** (QG-7) [GR1-01] | Donnelly–Giddings start from the LQFT split vacuum with an **arbitrary collar U_ε**. Gauge-invariant dressings are nonlocal, so ordinary localization fails at O(κ). They construct a **first-order gravitational splitting** (subspaces labelled by total Poincaré charges). **No positive minimum ε is introduced** | **PRIMARY / SOURCE-TEXT VERIFIED** (owner) [GR1-10] | **NO POSITIVE COLLAR SCALE DERIVED; SUBSYSTEM NOTION REPLACED / REPRICED** |
| **Fine-grained gravity** [GR1-01] | Raju: in the gravity settings treated, observables near the boundary of a Cauchy slice fix the state on the whole slice (holography of information), so **the ordinary split property fails** there. The examples are explicitly scoped | **PRIMARY / SOURCE-TEXT VERIFIED** (owner) [GR1-10] | **ORDINARY d_min^incl MAY BE NOT APPLICABLE / BLOCKED in these settings** |
| **Crossed-product controls** [GR1-06] | **Witten:** in a specific emergent large-N black-hole setting, type III₁ → type II∞ crossed product. **CLPW:** a gravitationally dressed de Sitter static-patch algebra of type II₁. These are strong examples of gravity changing algebraic type, **not a universal theorem for arbitrary local regions** | PRIMARY ARXIV ABSTRACT VERIFIED (each) | **no collar scale.** Whether G1-P1 applies: see below |

**G1-P1 scope under type change** [GR1-06]:
- **If** the gravitational setup supplies the relevant factor / commutant sharp pair (M, M′), non-type-I remains sufficient
  for the normal-product obstruction.
- The crossed-product papers **by themselves** do not establish that geometric pairing for arbitrary complementary
  regions.
- "Type II" is therefore **not** taken to imply G1-P1 for geometric regions.

## 3. Preparation level under gravity

### 3.1 The chain (GS-P1, heuristic)

**Premises:**
- QG-1 – QG-5, plus QG-6 in its spherical (Misner–Sharp / Hayward) form.
- **F:** N free massless scalars, in the flat local approximation d ≪ R.
- **L_all** [GR1-02]: every admissible product-state preparation carries at least the required gravitational mass within
  R + O(d).
  - The numerics establish localization **only for the minimum-energy Gaussian product state**: 97–99.7% of its excess lies
    within 3d in 3+1, and ≥ 99% in 1+1.
  - **They do not establish L_all.** A global bound E_total ≥ E_prod(d) does **not** by itself imply
    m(R + O(d)) ≥ E_prod(d)/c². A higher-energy product state could move positive energy outward or exploit negative local
    energy densities.
  - Owner decision 2: **L_all stays a heuristic premise in this campaign.**

**Steps:**
1. **Flat energy cost** (GS-L1 + transverse-mode decomposition; Part B). Minimum product-state energy per unit area of a
   planar collar:

       E/A = α₃ N ħc / d³,     α₃ = 0.00709, 0.00689, 0.00685, 0.00679   (d = 2, 3, 4, 6 sites)

   - f(u) = d·e(d, u/d) is d-independent to a few percent where the integrand is concentrated (u ≲ 2).
   - Conservative: **α₃ = 0.0068 ± 0.0003** per scalar.
   - **Additivity over fields and the spherically symmetric minimizer use covariance-level averaging only** [GR1-04]: the
     product constraint and the energy are linear in Γ, and the uncertainty cone is convex and invariant under the
     relevant symmetries.
2. **Spherical collar at R ≫ d:** E ≥ 4πR² α₃ N ħc / d³ · (1 + O(d/R)).
3. **Trapped-sphere criterion** (spherical GR: a sphere is trapped when E > r/2 in geometric units). With L_all,
   m(R + O(d)) ≳ E/c².
4. **Result:** d³ ≳ 8π α₃ N ℓ_P² R.

**Scope of the coefficient** [GR1-05]:
- **Robust candidate scaling:** d\* ~ (N ℓ_P² R)^{1/3}.
- **Model-dependent prefactor:** (8π α₃)^{1/3} ≈ 0.56. It is **not theorem-grade**. Substituting the flat collar energy
  for the semiclassical mass near threshold uses all of the following:
  - semiclassical backreaction;
  - the flat local collar estimate;
  - L_all;
  - neglect of gravitational binding and of stress-tensor fluctuations.
- Curvature ≫ ℓ_P⁻² is **not by itself sufficient** for the validity of the semiclassical Einstein equation if
  stress-tensor fluctuations are anomalously large. A product state across a thin collar is a candidate for exactly that.

**Sizes** (N = 1; prefactor-dependent):

| R | d\* | d\* / ℓ_P |
|---|---|---|
| 10⁻¹⁵ m | 3.6 × 10⁻²⁹ m | 2 × 10⁶ |
| 1 m | 3.6 × 10⁻²⁴ m | 2 × 10¹¹ |
| 6.4 × 10⁶ m | 6.6 × 10⁻²² m | 4 × 10¹³ |
| 1.4 × 10²⁶ m | 1.8 × 10⁻¹⁵ m | 1 × 10²⁰ |

### 3.2 Firewall audit of GS-P1

| check | result |
|---|---|
| GF-1 (Planck insertion) | passes. The bound comes from a consistency condition, and the form ℓ_P^{2/3} R^{1/3} is not ℓ_P |
| Supplied scales | ħ, G, c, R, N are all supplied. **Not compression** (GF-8) |
| G → 0 | d\* → 0: the flat result is recovered |
| GF-2, GF-4 | the entropy divergence and QNEC are not used |
| GF-7 (level) | a **preparation-class** statement. Normal product states still **exist** for every d > 0 in the flat net |

### 3.3 Terminal, grade, novelty [GR1-02, GR1-03, GR1-09]

**PREPARATION-CLASS CONSTRAINT: CONSTRAINED-NONUNIQUE — CONDITIONAL ON PREMISE L_all — HEURISTIC.**
- The collar family d < d\*(R) is excluded for horizon-free product preparation, **conditional on L_all**.
- d ≥ d\* is **not excluded by this test**. Nothing is selected.

**Novelty** [GR1-09]:
- **No novelty claim for the scaling.** The cubic-root form (ℓ_P² R)^{1/3} has older analogues in Károlyházy / Ng–van Dam
  distance-uncertainty arguments (Ng & van Dam, gr-qc/9906003, Found. Phys. 30, 795 (2000)).
- The mechanism here differs: disentanglement energy plus collapse. **α₃ and the product-state-energy route have not been
  exhaustively novelty-searched.**
- **Braunstein–Pirandola–Życzkowski** (arXiv:0907.1190; PRL 110, 101301 (2013)) is cited **only as conceptual
  "energetic curtain" precedent**. It concerns black-hole entanglement / information and is **not** a prior derivation of
  α₃ or d\*.

### 3.4 The other gravity–entropy routes

| route | output | terminal |
|---|---|---|
| **Dimensional analysis** | d = ℓ_P | **RELOCATION / SCALE INSERTION** (GF-1) |
| **Jacobson cutoff reading** | ε ~ √N ℓ_P. But ε is a regulator, not a collar; the small-ball and equilibrium premises apply (GF-5), and Susskind–Uglum blocks the literal-cutoff reading (GF-2) | **RELOCATION** |
| **Bousso-type** [GR1-08] | The covariant entropy bound constrains entropy on suitable light-sheets. The step I(A:B)/2 ≤ A/(4ℓ_P²) for the collar's vacuum mutual information (κ_I ≈ 0.005) is **not itself the Bousso conjecture**. It gives d ≥ ≈ 0.1 √N ℓ_P, which is sub-Planckian for N ≲ 100 | **HEURISTIC APPLICATION OF A CONJECTURE**; supplies **no G1 evidence** |
| **QNEC** | no gravity | control only (GF-4) |

## 4. Payoff

d\*(R):
- is not theorem-grade (L_all; prefactor);
- is not input-invariant (it depends on R and N);
- is ordinary semiclassical GR + QFT reasoning;
- is unobservable.

**NO DISTINCTIVE PAYOFF. ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS (preserved).**

## 5. Residual update

| component | after G1 (repaired) |
|---|---|
| d, inclusion level, fixed background | supplied; d_min^incl = 0 under the split assumptions |
| d, inclusion level, perturbative gravity | no positive scale; **subsystem notion replaced / repriced** (gravitational splitting) |
| d, inclusion level, fine-grained gravity (scoped examples) | **ordinary split may fail; d_min^incl not the right variable** |
| **(H_cross = 0, d)\|_horizon-free** | **candidate constrained preparation class** (H_cross / H_epoch side), conditional on L_all, heuristic |
| G1-P1 | applies wherever an exact non-type-I factor / commutant pair is supplied; crossed-product papers do not by themselves supply the geometric pairing |

All results are **auxiliary to canonical GRUT**.

## 6. Procedural notes

- CLARABEL reports `optimal_inaccurate` for a few modes. Modes below 10⁻⁵ are dropped and replaced by an exponential tail
  (≤ 2 × 10⁻⁵ in α₃).
- The Part-A a = 1/8 run was dropped for cost.
- The Part-B box (L = 20) affects κ_I by about 10% and α₃ by about 1%.
- The docstring of `g1/g1_dmin.py` predates GR1-04. Its proof sketch is superseded by §1.2 [GR1-04].
