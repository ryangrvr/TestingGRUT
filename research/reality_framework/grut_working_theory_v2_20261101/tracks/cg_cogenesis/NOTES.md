# Co-genesis: can one CP-odd source fix both phases? Result

Claude, for the owner (issue #13). The candidates, pass condition and failure conditions were committed first, in
[PREREGISTRATION.md](PREREGISTRATION.md) (dfc0a4a). This is my own analysis with no separate skeptic. Every choice is a
model choice. Nothing here validates or excludes anything.

## Answer

**All four candidates fail.** No single CP-odd source of the kinds pre-registered fixes both the holder phase and the
lepton phase within the pass condition. That condition is η_B within ×2 of 6.12×10⁻¹⁰, with no phase fitted.

**The decisive number** [C, `code/eta_vs_delta.py`, which writes `results/eta_vs_delta.json`]. On the M2 histories,
with e2e_card's corrected leptogenesis code run unchanged, η_B is about 10⁻⁷ for any order-one lepton CP parameter.
That is ×90–290 above the measured value.

| lepton CP parameter δ_L | 10⁻³ | 0.1 | 0.3 | 1.0 |
|---|---|---|---|---|
| η_B (range over the five M2 histories) | 5.7–8.9×10⁻¹⁰ | 5.4–8.4×10⁻⁸ | 1.2–1.8×10⁻⁷ | 1.1–1.8×10⁻⁷ |
| ratio to measured | ×0.9–1.5 | ×89–138 | ×190–291 | ×176–288 |

The rates and m̃_j at each δ come from the seesaw row at that δ, so the stronger washout at large δ is included. η_B
saturates near 10⁻⁷. The measured value needs an effective δ_L ≈ 10⁻³.

The holder, by contrast, needs a **large** CP-odd angle: α = arg(μ_eff − μ) = arg(10i − μ) + δ ≈ 2.19 rad ± 0.07
(M1r's α_I = 2.192). Any co-genesis must therefore transmit an O(1) holder phase into a lepton phase of about 10⁻³.

## Candidate by candidate

- **K1, a constant spurion with charges (q_h, q_L) ∈ {1, 2, 3}.** δ_L = (q_L/q_h) × 2.19 rad. That gives effective CP
  parameters |sin(q_L θ)| = 0.67, 0.89, 0.99, 0.81, 0.15, 0.94, 0.29 for the ratios 1/3, 1/2, 2/3, 1, 3/2, 2, 3.
  - The smallest, 0.15, still gives η_B ≈ ×100–200 too large.
  - Fitting would need q_L/q_h ≈ 4×10⁻⁴, or a cancellation to about 10⁻³. Both count as separate tuning under the
    pre-registered rule.
  - **Fails.**
  - Caveat [A]: e2e_sm's δ is a CP-odd direction in the Casas–Ibarra matrix, not a literal Yukawa phase, so the map
    from charge × θ to δ depends on the texture. The saturation at about 10⁻⁷ for any O(0.1–1) CP violation does not.
- **K2, the χ vacuum phase v_χ/f₂ = 1.31 rad.** χ ≈ 0 during inflation, displaced only by the small tadpole bias. So a
  holder phase built from e^{i n χ/f₂} is absent when the held state is prepared. **Fails structurally** [D]: the phase
  is absent when needed. Its leptonic side alone would also overproduce, as in K1.
- **K3, the clock angle θ₁ = σ/f₁.**
  - During inflation θ₁ ≈ σ_c/f₁ ~ 10² rad. It rotates by about |η| θ₁ ≈ 1 rad per e-fold near the slice [E], so it is
    not a held phase.
  - At leptogenesis σ sits at its vacuum σ = 0, so θ₁ = 0 and there is no leptonic phase.
  - **Fails structurally** [D/E]: on the lepton side the phase is absent when needed.
- **K4, a frozen axion-like phase ϑ.** Numerically identical to K1, with θ = ϑ_* an initial condition. **Fails**, as
  K1. Its correlated-isocurvature check is moot.

**Sign.** Declared vacuous in advance: X carries no SM charge, so the sign of the dark charge is a convention. It is
not counted.

**Chance baseline.** About 10% per candidate (the ×2 window against a ±3-decade spread). Observed passes: 0 of 4.

## What the failure teaches

The two phases live at opposite ends of the allowed range:
- **the holder phase is O(1)**, near the maximal-torque region, and is needed during inflation;
- **the lepton phase must be about 10⁻³**, and is needed at T ≈ 10¹² GeV.

A common source would have to bridge that ratio. The one structure that could do so without explicit tuning is
**radiative transmission**: the lepton phase induced at loop level from the holder phase, δ_L ~ (g²/16π²) α, which
gives about 10⁻³ for g ≈ 0.3–0.4 [A]. It is not evaluated here. It would need its own model of how the hidden holder
sector talks to the leptons, and its own pre-registration. It is the only route this test leaves open.
