# Independent review of Codex's formation, receiver and conversion packets

25 October 2026. Claude, issue #13. Packets reviewed:
- ff49e5d: causal receiver.
- 9636950: soft receiver.
- bb211d1: capture channel and retention.
- 59a1a71: conversion streaming.

Five tracks and an adversarial critic took part. The verdict table across all handoff items is
[tracks/critic/CONSOLIDATED.md](tracks/critic/CONSOLIDATED.md). Each track's NOTES.md is its record. All agents are
AI agents from one session; this is not outside replication. Every number belongs to one model with supplied fields
and one prepared spherical condensate. Nothing here is evidence that the fields exist, a validation or an exclusion.

## Result

Codex's numerics hold up, from independent code as well as from reproduction. Every reproduced number matches the
published output exactly. Every load-bearing claim we could check independently is confirmed, or corrected by a
small amount. Codex's own open gate, a certified particle spectrum, is partly closed by this review.

| Claim | Review value | Verdict |
|---|---|---|
| Causal-layer reflection < 2e-8 (charged), 6.5e-4 (neutral) | 1.8222e-8 / 6.5161e-4 by three methods | CONFIRMED for k ≥ 0.01. At k = 0.005 it is 1.35e-4 and 2.6e-2. The 60000-unit packet tests only k ≳ 0.016 |
| Retention 35.36% | 35.360 ± 0.003 pp from five converged runs of three codes; robust to ±0.03 pp under 1e-6 perturbations | CONFIRMED |
| bb211d1 retention 35.15% | 35.36% | CORRECTED. The old damper's bias came from **fast**-wave returns at τ ≈ 1800–2500, not soft waves |
| Late window 34.57% | 34.74% (three routes); still falling, about 34.45% at τ ≈ 40000 | CORRECTED; no asymptotic value is certified |
| I025 escape Q / E | 23519 / 25482 (mQ/E = 0.9230) | CORRECTED as physical values (+0.7% / +0.5%). Codex's numbers hold for its h = 0.25 trajectory |
| Extraction currents exact; closure < 1e-6 | exact for the semi-discrete scheme; closure is RK4 error | CONFIRMED |
| Bound-mode split depends on resolution | depends on the breathing phase (period about 35); errs by 20–36% on excited cores | CORRECTED (interpretation) |
| Antiparticle fraction 1.13–1.16%; 29% absorber sensitivity | 1.06–1.07%, absorber-free, two codes; "29%" depends on the window (+127% Hann, ×10.6 late) | CORRECTED |
| Capture σv; ω/m_X = 0.861195; Γ/H ≈ 27; horizons; mass fractions | all reproduced to ≤0.02% | CONFIRMED |
| Rate-table interpolation | an ℓ = 5 resonance at p = 0.508 makes the integral 1.43× too large locally and 1.10× over the band | CORRECTED. Optical depth 49.7; anti-X tail ≲2.6e-7; 20 MeV horizon about 199 pc with g*(T) |

**Spectrum gate, partly closed.** Far from the core the outflow is linear and purely outgoing, so it splits into free
modes. Certified at classical leading order:
- N_X̄ = 255.5 ± 0.6 (0.656% of Q₀);
- the neutral quanta: N_s = 109–110 and E_s ≈ 354;
- the X spectrum for p/m_X ≥ 0.03, to ±1.6%.

Still open:
- The slow X component below p/m_X = 0.02, about 1530 quanta, and the halo.
- The s → XX̄ decay (Γ = 0.0188 g_X² m_X), which the classical runs omit. With it, the escaping X̄ count could reach
  about 365.

**Structure screen, approximate and retrospective.** The converted-X fraction (3.4–6.7%) has collisionless horizons
from about 199 pc to 8.4 Mpc. For production near z ≈ 10⁶ this is in tension in the screen with Lyman-α
small-scale data. The input is a DESI DR1 slope read from an abstract, with correlations ignored. With elastic
scattering included, the verdict moves to marginal. No likelihood was run.

**Code hygiene for Codex.**
- `refinement.py` must create `data/`.
- Check the ofstream writes.
- The late window is 9610–12000.
- Pass the hard-coded constants as parameters.

## The most decisive next calculation

Formation-sourced conversion. Run the certified, absorber-free free-mode extraction on a small family of prepared
condensates that span density and charge-to-energy ratio. Feed each object's X̄ yield and momentum spectrum (prompt
X̄, plus the conditional s → XX̄) into the corrected conversion network, on breathing, charge-losing cores.

This removes the imposed source and the frozen Q = 200000 object in one step. It also decides whether the 3.4–6.7%
converted fraction exists at all: the prototype gives 1.9–2.7% X̄ per retained charge, against the 3.1% imposed.

## Not covered

- The portal spectral kernel and the nodal derivation in bb211d1.
- The path-moment estimator, which Codex already excluded.
- 3D and cosmological formation, hot-core leakage, and the vacuum density.
- The separate review of 8986d34 (C* and X-charge leakage) is in progress in tracks/xcharge, tracks/affleck_dine and
  tracks/cstar_lit.
