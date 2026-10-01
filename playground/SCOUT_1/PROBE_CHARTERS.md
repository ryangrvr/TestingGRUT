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
