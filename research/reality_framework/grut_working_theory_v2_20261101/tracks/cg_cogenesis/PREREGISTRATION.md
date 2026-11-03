# Co-genesis: can one CP-odd source fix both phases? Pre-registration

Committed **before** any candidate below is evaluated. Claude, for the owner (issue #13). It follows the owner's rules
(relayed from their discussion): state the source before evaluating it; it must not use Ω_DM/Ω_b, S/Q, the M2
windows or either observed asymmetry to define itself; it must return both phases or a relation between them; report
every candidate; fix the pass condition in advance.

## What M2 has now (the inputs being replaced)

- **The holder phase.** The holder term is −(g_h e^{iα} s Σ*Y² + h.c.), with s a real spectator. The CP-odd quantity
  is α together with sign(s_*). The held state needs μ_eff at the slice near 14–16 H_i at arg μ_eff ≈ 100–109°. That
  is far from the CP-conserving direction (μ_eff real). This phase must be present **during inflation**, when the
  state is prepared.
- **The lepton phase.** It is δ_L in R(δ_L) = R0 exp(δ_L A) (e2e_sm). η_B needs δ_L ≈ 8–9×10⁻⁴ (window 10⁻⁴–2×10⁻³;
  e2e_card's corrected value). This phase must be present **at leptogenesis**, at T ≈ 2–3×10¹² GeV, t ≈ 10³–10⁴ H_i⁻¹.
- **What stays as it is.** c0 and s_* stay supplied in every candidate; co-genesis does not touch them. H_i stays
  calibrated to Ω_DM. So the test is whether η_B (equivalently Ω_DM/Ω_b, given Ω_DM) is **predicted** once δ_L is
  tied to the holder phase.

## The candidates (all reported)

- **K1. A constant spurion.** One CP-odd spurion with phase θ breaks an otherwise exact CP. The holder phase and the
  lepton phase carry integer charges q_h and q_L ∈ {1, 2, 3}, so the CP-odd holder angle is q_h θ and δ_L = q_L θ.
  - In the CP limit (θ = 0) the holder is CP-conserving, so arg μ_eff − arg(real direction) = q_h θ.
  - θ is fixed by the held-state requirement; δ_L then follows.
  - Charges above 3, or any hierarchical ratio, count as tuning.
- **K2. The χ vacuum phase.** The common source is v_χ/f₂ = 1.31 rad (e2e_sm's axial completion), which already sits
  inside δ_L. The holder phase would be n_h v_χ/f₂ through a holder coupling to e^{i n_h χ/f₂}.
- **K3. The clock angle.** The common source is θ₁ = σ/f₁, which is nonzero during inflation. It enters the holder and
  the N1 phase.
- **K4. A frozen axion-like phase field ϑ.** It is light through leptogenesis and enters both phases with charges
  (q_h, q_L) ∈ {1, 2, 3}. This is K1 made dynamical: θ = ϑ_* is an initial condition. It adds a check of correlated
  isocurvature in both asymmetries.

## The pass condition (fixed now)

The candidate passes only if, with no phase fitted separately:
- the predicted η_B lies within a **factor 2** of the measured 6.1×10⁻¹⁰ (so Ω_DM/Ω_b lies within ×2 of 5.36);
- the sign is right;
- every structural requirement holds (each phase is present when it is needed).

**The sign condition is vacuous in M2** [A]: X carries no SM charge, so the sign of the dark charge is a convention.
The relative sign can always be chosen to give baryons. This is stated in advance, and the sign is not counted as
evidence.

## Failure conditions (from the owner's list)

- The source fixes one phase and leaves the other free.
- It needs separate tuning: a charge ratio outside {1, 2, 3}/{1, 2, 3}, or a cancellation tuned below 10⁻².
- A phase is absent when it is needed.
- It works only after the desired ratio is chosen.

**Chance baseline.** If a natural prediction of η_B spreads over about ±3 decades, the ×2 window covers about 10%.

## Method

- η_B(δ_L) on the M2 histories, using the lepton code already used: `tracks/e2e_sm/code` and
  `tracks/e2e_card/code/lepto_fix.py`, unchanged. Where δ_L is not small, use the full R(δ_L), not the linear
  coefficient.
- The holder angle from the M2-Q / M2-X definitions in `tracks/lf_chain`.
