# SCOUT-1 WAVE-1 PROBE CHARTERS (preregistered before any Wave-1 computation)

**Common rules:**
- Outcome labels are fixed here.
- Every apparent selector gets the BASELINE ten-point hostile test.
- "Supplied value renamed as a principle" ⇒ REJECT.
- Imported theorems are labelled KNOWN RESULT IMPORT, with source grade in `LITERATURE_LEDGER.md`.

---
## W1-C — Scaling/units and deformation-group non-selection (C06 + C07) — RUN FIRST
- **Question.** Is every earned predicate (E-1…E-22, as recorded) invariant under admissible deformations — time/rate rescaling `K → λK`, `t → t/λ`; pin shift within the gapped class; unitary relabeling of hidden sectors — so that no earned principle can fix the quotient components those deformations move?
- **Hypothesis.** All earned predicates are scale-free, so no dimensionful quotient value is earned-fixable. More generally, a selector invariant under G cannot fix a G-moved component, and any selector that does fix one must break G, i.e. supply a G-non-invariant datum.
- **Null.** Some earned predicate carries an intrinsic scale and fixes a dimensionful value.
- **Outcomes:** NON-SELECTION THEOREM (scope stated) / COUNTEREXAMPLE FOUND (an earned intrinsic scale) / UNFORMULABLE.
- **Fence.** No new predicates are invented. The E-list is the record's.

## W1-A — LSM/Oshikawa: symmetry + sector fix the soft-point momenta (C01; C02/C03 hostile siblings)
- **Question.** For U(1) charge conservation × lattice translation at filling ν, is the low-energy soft-point momentum `q* = 2πν` (mod 2π) fixed independent of statistics and interaction, while velocities stay free?
- **Hypothesis.** Yes. This is the LSM/Oshikawa/Yamanaka–Oshikawa–Affleck theorem: a dimensionless Q4 component is fixed by two supplied layers acting jointly.
- **Null.** The soft momentum moves with U or with statistics.
- **Tests.** Exact diagonalization (fermions, hard-core bosons, Bose–Hubbard U = 1, 4) at ν = 1/2 and ν = 1/4. Hostile:
  - break translation (staggered potential; a pinned site);
  - break U(1);
  - vary ν.
- **Outcomes:** CROSS-LAYER SELECTION (dimensionless component fixed; layers supplied) / NOT FIXED / UNRESOLVED.

## W1-I — Access + topology fix the one-particle contraction (C08, C09, C33)
- **Question.** Does the retained-site spectral measure `μ_r` together with path (chain) topology fix `K` uniquely? Without the topology, do cospectral-vertex counterexamples show non-fixing?
- **Hypothesis.** On the chain: yes, by the Jacobi inverse spectral theorem (Lanczos/Stieltjes). Without topology: no.
- **Outcomes:** CROSS-LAYER SELECTION (topology + access) / NON-SELECTION (counterexample) — both may hold at different scopes.

## W1-P — Complete passivity fixes the forcing-law temperature profile (C12)
- **Question.** For a bath of independent modes, does complete passivity — no work extractable from any number of copies under cyclic unitaries — force a single temperature (KMS), collapsing Q2's mode-temperature freedom to `(J(ω), T)`?
- **Hypothesis.** Yes (Pusz–Woronowicz). Single-copy passivity allows `T(ω)` freedom (monotone occupations).
- **Tests.** Explicit work extraction from copies of a non-KMS passive state; single-copy passivity check.
- **Outcomes:** SELECTION PRINCIPLE (premise-priced) / NOT FIXED.

## W1-G — Quantum-foundations selectors for the outcome weight law (C10, C11)
- **Question.** Do Gleason (d ≥ 3), Busch (POVMs, d ≥ 2) or composition plus no-signalling force `h(p) = p`?
- **Hypothesis.** Yes, each given a supplied Hilbert-space/composition structure. Projective-only d = 2 does not (explicit non-Born frame function).
- **Outcomes:** SELECTION PRINCIPLE (premise-priced) / NOT FIXED; plus the d = 2 counterexample.

## W1-R — Unitarity rigidity quantizes the IR central charge (C04)
- **Question.** Does unitarity plus conformal (Virasoro) symmetry with `c < 1` force `c` and the exponents into a discrete set (FQS/Kac)? Does a lattice parent land on a discrete value, independent of its microscopic couplings, at criticality?
- **Tests.**
  - Kac-table enumeration.
  - Exact diagonalization of the critical transverse-field Ising chain: entanglement fit of `c ≈ 1/2` across couplings tuned to criticality.
  - Hostile: an XX chain (c = 1, outside `c < 1`); off-critical gap.
- **Outcomes:** RIGIDITY SELECTION (premise-priced) / NOT FIXED.

## W1-L — Reconstruction axioms select the lift's number field (C32)
- **Question.** Does local tomography (with composition) exclude real and quaternionic quantum theory, fixing complex QM, and so remove the "complex structure" price SCOUT-0 found for the lift?
- **Tests.** Parameter counting `K_AB = K_A·K_B` for real, complex and quaternionic state spaces (explicit dimensions).
- **Outcomes:** SELECTION PRINCIPLE (premise-priced) / NOT FIXED.

## W1-S — Spectral dimension fixes the edge exponent (C20)
- **Question.** For translation-invariant short-range couplings, is the retained-site edge exponent a function only of spectral dimension and boundary type? That would reduce Q3's `γ` from a continuum to a discrete set.
- **Tests.** Chain (d = 1), square lattice (d = 2), cubic lattice (d = 3), bulk and boundary sites, several coupling strengths and anisotropies.
- **Outcomes:** CROSS-LAYER RESTRICTION (`γ ∈` discrete set; dimension supplied) / NOT FIXED.

## W1-F — Prediction-first: fluctuation–dissipation relations (C13, C14)
- **Question.** Does GRUT's earned structure force the Einstein relation or the KMS noise/dissipation ratio `2coth(πω/H)`?
- **Expected.** Both follow from supplied/borrowed FDT/KMS and are standard, so RETIRE as NON-DISTINCTIVE unless a GRUT-earned premise is found.

---

# Wave 2 (preregistered after ZOOM_OUT_02)

Common question: can a **dynamical attractor** replace a supplied *price* (preparation or tuning) found in
Wave 1? If it can, a dynamical structure could earn a weight-0 value. If the attractor needs its own supplied
structure, the wall is confirmed at the dynamical level.

## W2-QE — Born as a dynamical attractor (de Broglie–Bohm relaxation)
- **Question.** Starting from non-equilibrium `ρ₀ ≠ |ψ₀|²` in a 2D box, does Bohmian dynamics drive the
  coarse-grained `ρ̄ → |ψ|²‾` (H̄ → 0) without a supplied equilibrium preparation?
- **Hypothesis.** Yes for many-mode ψ (Valentini–Westman), with three prices:
  - the guidance law;
  - coarse-graining (fine-grained H̄ is conserved);
  - mode complexity — few-mode or stationary ψ do not relax.
- **Controls:**
  - equilibrium start (stays);
  - single mode (no motion);
  - two modes (incomplete);
  - dt halving (numerical convergence).
- **Outcomes:**
  - ATTRACTOR SELECTION (prices stated);
  - NO ATTRACTOR;
  - PARTIAL (relaxation only in a class).

## W2-ETH — KMS without supplied passivity (closed-system thermalization)
- **Question.** Does a closed non-integrable spin chain, started in a non-thermal product state, drive a small
  subsystem to the KMS state at one temperature? Is the temperature fixed by anything but the supplied energy?
- **Hypothesis.** Yes (ETH) for non-integrable chains. No for integrable chains, which go to a GGE. T is set by
  the initial energy density (supplied).
- **Outcomes:** ATTRACTOR SELECTION (prices) / NO.

## W2-FP — criticality without tuning (self-organized criticality)
- **Question.** Does a slowly driven, locally relaxing system (BTW/Manna sandpile) reach a scale-free state with
  no tuned parameter? Is a time-scale separation (drive rate → 0) the hidden tuning?
- **Hypothesis.** Scale-free avalanches appear only in the limit drive/dissipation → 0. The tuning moves into a
  rate ratio.
- **Outcomes:** UNTUNED FIXED POINT / TUNING RELOCATED / NO.

## W2-QE — FINE-GRAINED FIREWALL (pre-registered per external audit, BEFORE the running simulation's output was read)

- **Frozen fact.** Under pilot-wave evolution, `f = ρ/|ψ|²` is advected by the trajectory flow (`∂_t f + v·∇f = 0`).
  The flow preserves the `|ψ|²` measure. Consequences:
  - the fine-grained relative entropy `H_fine = ∫ρ ln f` is **exactly conserved**;
  - there is no fine-grained dissipative attraction to f = 1;
  - Valentini's H-theorem is a **coarse-grained** statement, and it does not prove equilibrium is always reached.
- **Required reports:**
  1. a fine-grained quantity (`H_fine` conservation, via the backward-trajectory evaluation of f);
  2. the coarse-grained H̄;
  3. dependence on coarse-cell size;
  4. dependence on mode count;
  5. dependence on initial microstructure;
  6. recurrence / failure-to-relax controls;
  7. a low-mode control;
  8. whether apparent relaxation survives refinement of the coarse graining.
- **Adjudications:**
  - **A. COARSE-GRAINED RELAXATION — PREMISE-PRICED.** The ledger keeps the pilot-wave law, ψ (mixing), the
    initial ensemble (no microstructure) and the coarse-graining.
  - **B. ROBUST DYNAMICAL RELAXATION.** This would require proof that it is not a coarse-graining artifact.
  - **C. NO RELAXATION.**
- **Fence.** Outcome A must not be described as "Born derived from dynamics without preparation".

## W2-DT — dimensional transmutation / anomalous scale generation vs TC-1 (required before TC-1 promotion)
- **Question.** Can a theory with no explicit classical dimensionful selector dynamically generate a finite
  non-zero physical scale, and so evade TC-1?
- **Cases:**
  - asymptotically free RG (Λ_QCD-type; also the exactly solvable large-N Gross–Neveu / O(N) gap equations);
  - Coleman–Weinberg.
- **Track** where the absolute scale enters, distinguishing:
  1. breaking of classical scale symmetry;
  2. an RG-invariant scale;
  3. an absolute numerical scale in physical units.

  For each, check whether it needs a renormalization boundary condition, a measured coupling at a reference
  scale, a vacuum choice, or another supplied dimensionful datum.
- **Outcomes:**
  - **TC1-SURVIVES-ANOMALY** — a relational scale is generated, but its absolute value needs G-breaking boundary/reference data;
  - **TC1-COUNTEREXAMPLE**;
  - **TC1-NEEDS-REFORMULATION** — the deformation group is anomalous.

---

# D1-SCOUT (owner ruling after ZOOM_OUT_04) — SCOUT-ONLY premise experiment, NOT a GRUT premise adoption

**Fence lift (SCOUT-only, new hypothesis class):** `D1-SCOUT = pin = 0 + spatially extended nonlinear (conserved)
dynamics`. GRUT-RAI is untouched and SCOUT-0 stays frozen.

**Primary question (selector-of-the-selector).** Does structure inherited from GRUT uniquely determine the relevant
nonlinear universality class? Or does D1 merely exchange a supplied IR label for a supplied nonlinear operator/class?

**Inherited objects:**
- C1-a/L0-1 chain `ẋ = −Kx + noise`, with K reciprocal (E-1);
- L0-1d ring with cycle affinity (E-4, asymmetric K, a declared class);
- the S2 C-B drift `dx = [−Kx − 4βx^{∘3}]dt + B dW`, with `Q = 2 diag(T_i)` (E-15: gradient, odd, on-site, additive
  site noise, supplied `T_i` profile).

**Steps:**
- **D1-0 formulability map:** field type, dimension, symmetries, locality, conservation, parity, time reversal,
  noise, access. The most general low-order local equation allowed. No silent KPZ import.
- **D1-1 operator uniqueness:**
  - A. UNIVERSALITY-CLASS-SELECTED;
  - B. OPERATOR-CLASS-SUPPLIED;
  - C. NONLINEARITY-FORBIDDEN.
- **D1-2 RG flow** (only if an operator survives): coupling / crossover / exponents / amplitudes / symmetry labels.
- **D1-3 symmetry hostiles:** parity, number of conserved fields, range, momentum conservation, detailed balance,
  dimension, boundary. Each is mapped `assumption → operator → IR class` and tagged earned / supplied / D1-only.
- **D1-4 pin role:** pin > 0 vs pin = 0 under the same nonlinear dynamics. Does the pin select, or only reveal?
- **D1-5 information accounting:** inputs before and after, separating UV values forgotten, operator content,
  symmetry class, dimension, crossover scales and universal IR numbers. TRUE NET SELECTION only if the
  specification of the IR law becomes smaller.
- **D1-6 GRUT bridge:** can the field nonlinearity be obtained by local/conservative extension of an admitted GRUT
  term (S2's drift)? If it must be newly supplied, say so.

**Theorem target:** if B repeats, a **SELECTOR RELOCATION NO-GO** candidate. If A, the first GRUT-adjacent selector.

---

# W4-OG — operator generation / closure hostile against TC-4 (owner ruling after ZOOM_OUT_05)

**Question:** is the irreducible supplied object the OPERATOR LIST, or only the FIELD / SYMMETRY / REPRESENTATION
content that determines the RG-closed operator space `O_allowed(Φ, d, S)`?

**Steps:**
- **W4-1 closure principle.** Distinguish four kinds of operator:
  1. forbidden by symmetry;
  2. allowed but absent (fine tuning);
  3. protected from generation (non-renormalization / integrability / topology);
  4. redundant (field redefinitions / equations of motion).
- **W4-2 background shift.** The D1 cubic current `J = λ₃u³` expanded about `u = u₀ + φ` has `λ₂,eff = 3λ₃u₀`. Same
  microscopic λ₃, with u₀ = 0 and several u₀ ≠ 0. Does the class cross from EW/marginal to KPZ? What is u₀ (filling,
  preparation, spontaneous, external)?
- **W4-3 generated operator.** A model with `λ₂,bare = 0` that is Z₂-breaking (λ₂ allowed) only through a
  non-current ingredient. Show `λ₂,eff ≠ 0` after coarse-graining (controlled calculation + numerics). Run the Z₂
  control.
- **W4-4 SSB.** Action symmetry vs state symmetry vs effective algebra around the state. What selects the broken
  vacuum (mass sign, boundary condition, thermodynamic limit, infinitesimal source, initial condition)? Price each.
- **W4-5 anomaly.** Anomaly-fixed coefficients (WZ/WZW, π⁰→γγ, Chern–Simons levels, LSM): coefficient fixed vs
  representation supplied.
- **W4-6 `CONTENT_SELECTION_LEDGER.md`.** Levels:
  - C1 couplings;
  - C2 operator basis / Wilson coefficients;
  - C3 symmetry class;
  - C4 field / representation;
  - C5 dimension / locality / composition.

**Outcomes:**
- A. TC4-REFUTED-AT-OPERATOR-LEVEL;
- B. TC4-SURVIVES-SCOPED;
- C. TC4-REFINED-HIERARCHY.

**D3 gate:** decided only after W4 and its zoom-out.
