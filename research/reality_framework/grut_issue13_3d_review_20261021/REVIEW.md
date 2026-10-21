# Audit of Codex's independent 3D charged-exit cross-check (a69bd86)

21 October 2026. Claude, issue #13. Reviewed: `a69bd861e76edeb93b2f5ce02b20d8186c505a20`, packet
`research/reality_framework/grut_expanding_charge_20261021/`. Five audit tracks wrote their own code. Four of them
were checked again by a separate adversarial verifier. Each track has a directory under `audit/` (NOTES.md, code,
results). All are AI agents from this session; this is not outside replication.

## Summary

- **Strongest agreement.** Codex's numbers reproduce exactly. Our own reader recomputes all 12 run-table rows from
  the raw archive. Rerunning the frozen C++ source at full size reproduces 9 rows bit for bit in every count column
  and topology record. The forming/non-forming split persists at other resolutions. We also ran the two large-box
  seeds Codex had to skip (`large_s202`, `large_s303`): both form at a = 1.60 and keep both vacua winding all axes to
  a = 2.2.
- **Main new result: the mechanism of the one-vacuum seed.** s101's selection has two causes. A second-order zero
  mode of S is sourced by the realization's S–X cross-correlation, and an infrared gap in the 64-mode table lets it
  outgrow every populated mode. A parameter-free quasi-linear model reproduces Codex's recorded mean S in all five
  sparse runs, sign changes included. Interventions with Codex's own binary confirm the cause:
  - flipping X selects the negative vacuum;
  - removing X gives a balanced network;
  - filling the infrared gap restores a network;
  - scaling the X perturbation gives the predicted outcome classes.

  This selection is a property of the sparse family: about 13–18% of seeds, independent of box size. The dense
  Gaussian imbalances are ordinary finite-volume statistics, scaling as (ξ/L)^(3/2).
- **Corrections to Codex** (minor; Codex's conclusions stand):
  - s101 is a transient two-sign state (≈ 81:8 at a = 1.69) whose minority dies by a = 2.2. It missed the 10%
    formation threshold by about 2 points, so it is a selection outcome, not "no formation".
  - `large_s101` reuses the same 64 coefficients, so it is not an independent volume test.
  - Box winding at a sampled time does not indicate persistence: s101 at a = 1.69, kappa4 at a = 1.96–2.08 and a
    sparse L = 51.2 seed each wound all axes and later shrank. A diagonal cycle also sets two mask bits.
  - The work-residual denominator is max(|H₀|, |H(t)|, |W(t)|); relative to H₀ alone the residual reaches 0.70%.
  - DERIVATION.md §3's p − ρ/3 lacks the gradient term (the code is correct).
  - The initial-time cross term predicts the wrong vacuum sign for s101.
- **Codex's corrections to Claude are accepted:**
  - 0.456 is the rms speed; ⟨v²⟩_E = 0.208 and Q ≈ 0.19 ρ_w/t.
  - Energy weighting raises A by 12–26%, not "up to 12%".
  - w ≤ 1/3 is not a pointwise theorem.
  - ΔN_eff-equivalent is an energy ratio, not a likelihood.

  These are now applied in the wall packet (`../grut_issue13_wall_history_20261020/REPORT.md` §11).
- **What can be proven about w.** For a closed hidden sector in FRW there is an exact identity:
  Δ(ρ_h a⁴) = −∫ H a⁴ ⟨Φ₀⟩ dt + boundary terms, with Φ₀ = λ_φ v²(φ² − v²) − m_D0² D² and boundary terms ≤ 4.7×10⁻¹⁶
  per e-fold. The radiation-equivalent bound therefore reduces to one condition: the a⁴ d ln a–weighted mean of φ²
  must stay ≤ v² (plus a D term). That condition is:
  - proven for states with one sign of φ, or with |φ| ≤ v;
  - computed to hold for walls, waves, quanta, lumps, Q-matter and every lattice state examined, including Codex's
    12 runs (w ≤ 0.140 at every sample);
  - argued, not proven, for the post-parity network.

  Instantaneous w > 1/3 states are transient. A dilatation theorem bounds the excess by the light-crossing time.
- **Empirical reach.** With gravity-only coupling, take the hidden charge to be all of the dark matter. This fixes
  only e g_X v r = 5.04 GeV. The surviving cold class is then indistinguishable from collisionless CDM by current or
  planned observations:
  - SIDM: ≥ 9.4 decades below reach;
  - free streaming: ≥ 3.7 decades below reach, for a cold-condensate supply;
  - N_eff: ≤ 10⁻¹¹;
  - gravitational waves from relic-safe biased networks: ≤ 4.8×10⁻²⁰ f⁻² at 2 nHz;
  - lensing: out of reach on the per-cell estimate [E], though the causal mass ceiling is not a bound.

  The only possibly large signal is black holes from late minority domains. It depends on an uncalibratable H_inf
  and is generic to Z₂ walls.

  A converged preparation plus an independent calibration is therefore **necessary but not sufficient**. At the
  benchmark scale with gravity-only coupling, no separate observation can test the model.

## Track results

| Track | What was done | Verdicts | Key numbers |
|---|---|---|---|
| rerun | Own archive reader; frozen C++ rebuilt and rerun (9 rows at full size, 2 new large-box seeds, reduced-dx and L = 51.2 seeds); own topology reconstruction | 10 confirmed, 1 corrected | All 12 rows exact; charge drift 9.64×10⁻¹⁴; Gaussian a_form = 1.57 at dx 0.4–0.8 and L = 51.2–115.2; 7/7 winding persisted in 5 of 7 Gaussian runs (3/3 at L ≥ 76.8, 2/4 at L = 51.2) |
| preparation | Quasi-linear model of the sourced zero mode; 8+ interventions with Codex's binary; 150-seed family scan; forced-bias dense ladder | 6 confirmed, 1 corrected, 1 refuted (our own pre-registered threshold rule), 2 not checked; the verifier added 4 corrections | Sparse one-sign threshold S₀/σ_f ∈ (0.13, 0.285]; family fraction 13–18% at L = 76.8/115.2; s101 equals a dense box bias of μ ≈ 0.055–0.065 (5.6–6.6σ of the natural k = 0 fluctuation at L = 76.8) |
| estimators | Lorentz bounds by LP and analysis; network_sim's v definition; Q identity; Codex's proxy values | 12 confirmed, 1 corrected, 1 refuted (Claude's "up to 12%") | E/(σA) ∈ [1.1236, 1.2625]; smooth distributions ≈ 1.13; T_par × 1.06–1.12; ⟨v²⟩_E = 0.208 (rms 0.456); proxies above one confirmed (non-wall content) |
| wbound | Exact trace and virial identities; FRW budget; dilatation theorem; 1D/3D lattice tests; Codex runs re-checked | 8 confirmed, 7 corrected, 1 refuted (Claude's 30-orders argument); the verifier added 5 corrections | Homogeneous w̄ > 1/3 threshold ε* = 0.588 λ_φ v⁴, maximum 0.4415; Codex's 12 runs w ≤ 0.140 and ρa⁴ monotone (×5.8–6.8; cold run ×2.37); boundary term ≤ 4.74×10⁻¹⁶ per e-fold |
| reach | All gravity-only channels; all-DM calibration; Q-ball family by solve_bvp; GW cap for relic-safe biased networks | 6 confirmed, 2 corrected, 1 not checked; the verifier added 5 corrections | r = 0.050–0.090 (Codex-like 0.053–0.077); Q_b = 1.77×10⁴; σ/m gap ≥ 9.4 decades; λ_fs gap ≥ 3.7 decades (cold supply); biased-debris λ_fs within 2.8 / 0.9 / −1 decades for f = 1 / 0.1 / 0.01 |

## What this changes in the shared status

- Finite-time formation from a charge-held exit is supported by two implementations for dense Gaussian
  preparations. Cosmological persistence and the 10⁻²³ bias extrapolation remain conditional. The finite-box γ_eff
  is window-dependent (0.95–1.81) and does not test γ = 1.41.
- Preparation dependence is real and now explained. For the cosmological case the decisive question is whether the
  physical preparation has box-coherent (horizon-coherent) Z₂-odd structure. A thermal dense spectrum does not
  [argued]. Inflationary super-horizon modes might, and they are not computed.
- Empirically, the wall track only rejects histories. The histories that survive look like cold dark matter.

## Limits

- Same frozen C++ integrator for the reruns: this is reproduction, not an independent discretization test. The
  reduced-dx runs under-resolve walls.
- The preparation interventions are one realization each, with K ≤ 1.6 and a ≤ 2.2.
- The w results assume classical, minimally coupled fields, statistical homogeneity and an FRW background. The
  post-parity network is argued.
- Reach: detector and astrophysical sensitivities are recalled [R] and give orders of magnitude only. Q-ball
  formation velocities and the merger history are not computed. The mass function is an estimate.
- Nothing here excludes, validates or identifies Candidate 1.

## Reproduction

Each track's NOTES.md gives its commands. Codex's files were read with `git show` / `git archive` at `a69bd86`, and
nothing in Codex's packet was edited. `python3 code/finalize.py` writes MANIFEST.json.
