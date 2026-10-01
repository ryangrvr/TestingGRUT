# S2-3b RESULT — horocycle flow: unique measure vs time-average universality vs preparation forgetting

**Charter:** `probes/PROBE_CHARTERS.md` §S2-3b. Pre-registered in REPAIR 01, before the run. The three notions were
measured separately.

**Files:** `probes/S2-3b/s2_3b_horocycle.py`, with log `probes/S2-3b/s2_3b_horocycle.log`.

**System:**
- the horocycle flow h_t = [[1,t],[0,1]] on X₂ = SL(2,ℝ)/SL(2,ℤ) (unimodular lattices), simulated with Gauss reduction;
- the observable f = 1[|v₁|² < 0.5], whose Haar value is 3s/π = 0.4775 (analytic);
- sampler check: 0.4769, and 0.4776 after the flow.

**Wall attacked:** H.

**Labels:**
- **H-MEASURE SELECTED (compact case, KNOWN RESULT IMPORT) / A-COARSE FORGETTING**;
- **D-PRICED** (rigid homogeneous geometry);
- KNOWN RESULT IMPORT — SECONDARY, with numerics ✓ on X₂:
  - Furstenberg (unique ergodicity, compact quotient);
  - Dani (invariant measures on X₂);
  - Dani–Smillie (equidistribution of non-periodic orbits);
  - Marcus / Ratner (mixing and rates).

## 0. Verdict, notion by notion

**A. UNIQUE MEASURE — FAILS on the non-compact X₂ (the hostile variant).**
- Lattices with a horizontal vector, diag(a, 1/a)ℤ², lie on periodic orbits of period a².
- Each periodic orbit carries its own invariant probability measure. The orbit averages of f are 0.000 (a = 1) and
  1.000 (a = 0.5, 0.2), against 0.4775 for Haar.
- On a **compact** quotient Γ\SL(2,ℝ), A holds (Furstenberg). This is imported, not simulated.
- So A is earned only when D contains **compactness plus a lattice in a Lie group**.

**B. TIME-AVERAGE UNIVERSALITY — HOLDS for every non-periodic start tested; FAILS on the periodic set.**
Birkhoff averages of f (target 0.4775):

| start | T = 10² | 10³ | 10⁴ | 10⁵ |
|---|---|---|---|---|
| Haar-random | 0.4724 | 0.4771 | 0.4778 | 0.4774 |
| ℤ² rotated by √2 − 1 | 0.4756 | 0.4776 | 0.4774 | 0.4775 |
| ℤ² rotated by 10⁻² | 0.2468 | 0.4516 | 0.4780 | 0.4770 |
| ℤ² rotated by 10⁻³ | 0.0000 | 0.2464 | 0.4546 | 0.4751 |
| ℤ² (periodic) | 0 | 0 | 0 | 0 |

- **B holds without A here**, which confirms that the firewall is real: B is not A.
- Convergence is **non-uniform**: its time grows like 1/angle near the periodic set.
- "Every point" is replaced by "every point off an explicitly characterized exceptional set". The set is dense and
  Haar-null. That is a characterization, not a measure-priced "generic".

**C. PREPARATION FORGETTING — three separate measurements.**

*C-i. Decay of correlations (Haar, N = 4·10⁵).*
- C_ff: 0.2495 → 0.0240 (t = 2) → 0.0003 (t = 10) → noise (|·| ≤ 5·10⁻⁴).
- The K-dependent observable cos 2φ decays more slowly: 0.4999 → 0.0459 (t = 10) → 0.0019 (t = 100) → noise.
- So mixing is present, at a polynomial, observable-dependent rate.
- This is a statement **about Haar**. Haar enters here as the reference, which is H. On the compact quotient Haar is
  the *selected* measure.

*C-ii. Weak convergence of the accessible observable ⟨f⟩ under pushforward.*

| t | 0 | 5 | 10 | 20 | 50 | 100 | 500 |
|---|---|---|---|---|---|---|---|
| blob 1 (a.c.) | 0.000 | 0.208 | 0.519 | 0.462 | 0.477 | 0.477 | 0.478 |
| blob 2 (a.c.) | 1.000 | 0.827 | 0.459 | 0.478 | 0.483 | 0.477 | 0.479 |
| Dirac (irrational) | 0 | 0 | 0 | 1 | 0 | 1 | 1 (flickers forever) |
| periodic diag(½,2) | 1 | 1 | 1 | 1 | 1 | 1 | 1 |

Over a 10×10×10 coarse partition, the TV distance between the pushed blob ensembles goes 1.000 → 0.834 (t = 10) →
0.155 (t = 50) → 0.044 (t = 500). The sampling floor is 0.039.
- **Coarse forgetting holds, but only for a.c. preparations.**
- Dirac and periodic preparations are never forgotten.

*C-iii. Fine-grained information under pullback.*
- The pulled-back observable 1_{B₁} ∘ h_{−t}, applied by an agent who holds only the time-t lattice, separates the
  two preparations with **contrast 1.0000 at every t up to 500**.
- The round-trip error is ≤ 7·10⁻¹¹.
- **Fine-grained information is exactly conserved.** The flow is invertible and measure-preserving, as the firewall
  requires.
- The cost of accessing it grows **polynomially**. The number of 20³ cells occupied by h_t(B₁) is:

  | t | 0 | 1 | 2 | 5 | 10 | 20 | 50 | 100 |
  |---|---|---|---|---|---|---|---|---|
  | cells | 12 | 41 | 72 | 350 | 1282 | 3940 | 7651 | saturated (8000) |

  The effective exponent is ≈ 1.5–1.8 before saturation. Compare the doubling map's 2ⁿ.
- Zero entropy therefore means **the preparation is hidden from a fixed coarse readout only polynomially fast**.

## 1. Structure inserted vs forced

| Item | Status |
|---|---|
| unique invariant measure (compact quotient) | **forced by D**: compactness + homogeneous Lie-group geometry (imported) |
| unique invariant measure (X₂) | **not forced**: the periodic family ⇒ A fails |
| time-average value for every non-periodic point | **forced** (Dani–Smillie; numerics ✓) |
| ensemble forgetting | **only coarse (A)**, and **only for a.c. preparations** |
| the class "a.c. preparations" | **inserted**. On the compact quotient its reference measure is *selected* by D; what remains supplied is the *resolution limit* of preparation (no point preparations) → **A** |
| fine-grained distinguishability | **conserved** (forced by invertibility) |
| "horocycle-type" dynamics itself | **inserted, rigid**: Ratner rigidity; not a generic smooth flow → **D-PRICED** |

## 2. Accounting against T2-1′ (C5 → D ⊕ H ⊕ A)

**What moved.** In this toy, H splits in two:
- **the reference measure → D.** It is selected by the dynamics on a compact homogeneous space. No independent
  measure choice remains.
- **the preparation class → A.** Forgetting holds only for preparations an agent can make (spread, a.c.) and only for
  readouts an agent can perform (a fixed coarse observable). Point preparations and fine-grained pullbacks defeat it.

**What did not move.** Forgetting is never fine-grained, so "probability as forgotten preparation" is **A-COARSE**.
The price of selecting the measure is a highly special D: compactness, a lattice in SL(2,ℝ), parabolic flow. That is
the opposite of generic.

**Hostile record:**
- the non-compact case kills A;
- convergence near the periodic set is non-uniform;
- Dirac preparations never relax.

**H is not eliminated. It is partially re-expressed as (special D) + A.** This is the first probe where a piece of H
(the measure) moves to D without supplying a measure. It is scoped to rigid homogeneous dynamics.

**Status: S2-3b COMPLETE — H-MEASURE SELECTED (compact quotient, imported; fails on non-compact X₂) / B universal off
the periodic set / A-COARSE FORGETTING (a.c. preparations, coarse observables; fine-grained contrast stays 1.0000) /
D-PRICED (rigid geometry).**
