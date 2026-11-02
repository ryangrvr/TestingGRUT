# Brief: an explicit inflaton sector for M1r (3 November 2026)

Read [WORKING_THEORY.md](WORKING_THEORY.md) §4, §7–§9 and [BRIEF_DECISIVE.md](BRIEF_DECISIVE.md) first; the earlier
rules apply.

The decisive pass left the inflaton sector as the one place where M1r is unspecified:
- **(1) History before the slice.** The vector number is set at horizon exit, 3–4 e-folds before the end-of-inflation
  slice. It decides the branch at e = 10⁻⁶: static slice S/Q = 4.955, free X; rolling embeddings S/Q = 3.1–3.4,
  Q-matter.
- **(2) Decay history.** It sets K = a²H; R_vec ∝ K².
- **(3) No realisable SM reheating channel.**
  - The reflection-even inflaton forbids φ_I\|H\|²; M1r has zero SM portals; a singlet inflaton has no renormalizable
    Yukawa to chiral SM fermions.
  - Γ_I/m_I = 1/6 needs N_f y² ≈ 4 and is kinematically blocked at φ_I ≈ 0.8 M_P.
- **(4) Higgs survival.** Without a Higgs-coupled plasma within τ ≈ 8–10 H_i⁻¹, it fails at M_t + 1σ.
- **(5) The engine's inflaton is phenomenological.** It is a post-inflation quadratic oscillator with m_I = 3 H_i,
  started at a = 1 at rest with amplitude √6 H M_P/m_I. A potential that inflates at H ≈ 1.3–1.5×10¹¹ GeV with
  r ~ 10⁻⁷ must be a plateau or hilltop shape.

**Goal.** Construct explicit, labelled benchmark inflaton sectors and push them through the whole M1r chain. Every
choice is a model choice, not a consequence of closure.

## Requirements for a benchmark

- **CMB.** n_s = 0.9649 ± 0.0042, A_s = 2.10×10⁻⁹ (Planck 2018), r < 0.036 (BICEP/Keck 2021). H at CMB exit and at
  the end of inflation, so that the dark-matter calibration H_i ≈ 1.1–1.5×10¹¹ GeV can be met. Note H_i is the Hubble
  rate at the M1r slice.
- **Post-inflation oscillation.** Mass m_I/H_end and amplitude. Avoid the Y resonance 2m_I ≈ ω_YL = 22.25 H_i, i.e.
  m_I ≲ 5 or ≳ 30 H_i (REVIEW_REHEATING). Note any fragmentation or oscillon formation for steep minima (small-α
  attractors).
- **Reheating channel to the SM** that preserves the CP-sign bookkeeping. If φ_I → −φ_I reflection is kept, use
  singlet fermions, e.g. right-handed neutrinos N with y φ_I N̄N, N → LH (seesaw consistent with neutrino masses), and
  give the reflection's action on N. Otherwise state how α_I and the sign of φ_I set the asymmetry. Requirements:
  - perturbative width Γ_I and the kinematic condition y\|φ_I\| < m_I/2 during the oscillation;
  - the reheating temperature, K, the matter-era length;
  - the time a Higgs-coupled plasma appears (needed within τ ≈ 8–10 H_i⁻¹ for M_t above about 172.8 GeV);
  - gravitino or other bounds are not needed (no SUSY).
- **The holder.** Its coupling g_I is re-fixed so that μ_eff = 10i H_i at the slice, with the slice inflaton value of
  each model.

## Codes

- **Release engine.** `tracks/dq_entropy/code/rpdq.c`: the M1r engine with the fixed forced-Y conj and the portal
  counterterm option; reproduce its benchmark first. Parent: `tracks/fx_holder/code/rptri.c`; parent of that:
  `research/reality_framework/grut_formation_review_20261025/tracks/rp_reheating/code/rp.py` + `rpreheat.c`. The
  inflaton potential and decay are inside the engine (mode 0 oscillator); new potentials need a modified copy. The
  preparation must be re-solved at the new slice.
- **Vectors.** `tracks/dq_heat/code` (Proca modes on a recorded a(η) and F(η), horizon-exit dependence) and its
  skeptic's frozen-angle formula N_A = (F_v/F_i)² K^{3/2} H² m^{−1/2} I/8π².
- **Entropy.** `tracks/dq_entropy/code` (S/Q(W) closed form) and `tracks/dq_critic/code` (branch, calibration,
  verdict).
- **Higgs.** `tracks/higgs_reheat/code` (with the skeptic's normalisation fix in `tracks/higgs_reheat/skeptic`).
- **SM running.** `research/reality_framework/grut_closure_forcing_20261027/tracks/cf_vacuum`.

## Rules

As in BRIEF.md:
- write only in `tracks/<key>/` (nothing over 2 MB); scratch in `$SP/wt2/<key>/`;
- git read-only; no commits;
- at most 2 threads; foreground only;
- claim labels;
- cite primary literature with arXiv numbers; WebSearch via ToolSearch "select:WebSearch";
- never claim validation or exclusion;
- return a digest of at most about 700 words.
