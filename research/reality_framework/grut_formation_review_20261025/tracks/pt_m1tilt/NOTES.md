# Track pt_m1tilt: the φ-condensation transition under Codex's GRUT-M1 potential

1 November 2026. Claude review of Codex 5f5513b (issue #13) and of the GRUT-M1 variant d4747cb
(`grut_matched_candidate_20261027/CLAUDE_HANDOFF.md`). Follows tracks pt_matching and its skeptic.

**Everything here is [C]: computed in one homogeneous, classical X + φ model with a prescribed mass drift. It is not a
cosmological result, not a validation and not an exclusion.** No Codex code was imported or run. The code is a copy of
the pt_matching skeptic's own adaptive DOP853 integrator (`sk_dop.c`) with one potential mode added.

## 1. Question

Codex's packet 5f5513b: slowly matching the coherent X orbit to a massive free field preserves the radial action J_r,
hence about 1.76 X̄ per net charge (τ-20 orbit), "no extra counting postulate". pt_matching and its skeptic modelled the
actual transition (φ condensing while the homogeneous orbit runs) and found, for the untilted Z₂ potential, an
ε-independent **random split** n₋/q = 0.70 ± 0.02, N_φ/q = 0.47 (τ-20) and 0.197 / 0.124 (late orbit): the loss comes
from φ growing out of the Floquet-unstable symmetric state φ = 0. The skeptic also found that a linear tilt β = 0.3
conserves J_r.

GRUT-M1 replaces the φ potential by V_φ = λ_φ/4 (φ² − v²)² − ε_t v (φ − v)³ with ε_t = 0.002, λ_φ = 0.03. Then +v is the
only stationary point (derivative (u−1)[λu² + (λ−3ε_t)u + 3ε_t], discriminant −1.44×10⁻⁴ < 0), V(v) = V′(v) = 0 and
V″(v) = 2λ_φ v² are unchanged, and φ = 0 is no longer stationary (V′(0) = −3ε_t v³). **Does the slow transition under
Codex's own M1 potential return J_r (adiabatic, deterministic), or a reduced random split as in the untilted model?**

## 2. Model and method

Conformal variables χ = aφ, M = av (radiation era, a″ = 0), so V_φ,conf = λ_P/4 (χ² − M²)² − ε_t M (χ − M)³ exactly.
The full potential (`code/m1_dop.c`, mode 3):

  V = λ_X/4 r⁴ + g²/2 χ² r² + λ_P/4 (χ² − M²)² − ε_t M (χ − M)³,  λ_X = 10⁻³, g² = 10⁻², λ_P = 0.03, ε_t = 0.002.

- Force on χ: −g² r² χ − λ_P (χ² − M²) χ + 3ε_t M (χ − M)². Explicit s-derivative for the work integral:
  dV/ds = M′ [−λ_P M (χ² − M²) − ε_t (χ − M)³ + 3ε_t M (χ − M)²]. Both checked against central differences of the
  energy routine (6×10⁻¹⁰ and 4×10⁻⁷ relative, the finite-difference floor).
- Expanded, the tilt is the skeptic's linear tilt with β = 3ε_t/λ_P = **0.2**, plus +3ε_t M² χ² (the tachyonic quadratic
  near χ = 0 is 0.6 of the untilted one) and −ε_t M χ³.
- Drift M = M₀ e^{κs}, κ = ε·ω_r(initial orbit), M₀ = 0.5 M_c,peri with M_c = g r/√λ_P, as in pt_matching. End at
  M_end = 40 M_c,apo; checkpoints at 5, 10, 20 M_c,apo. Same adaptive DOP853 (rtol 10⁻¹¹; the work integral is an extra
  ODE component), same 60-X-period window, same measures with the vacuum at χ = +M exactly:
  n₋/q = ⟨[(p_X²/2 + g²M²r²/2)/(gM) − q]/(2q)⟩, N_φ/q = ⟨(δχ′²/2 + m_φ² δχ²/2)/m_φ⟩/q, δχ = χ − M, m_φ² = 2λ_P M².
- Orbits (own quadrature; identical to the skeptic's): τ-20 from `MATCH_L12.npz light_mean_q/p`: q = 8117.8816,
  J_r/q = 1.763725, T_r = 0.814466, ω_r = 7.71449, M_c,peri = 10.217, M_c,apo = 82.378; "late" (r₀ = 31.4946, same q):
  J_r/q = 0.480492, T_r = 1.06160, ω_r = 5.91858, M_c,peri = 18.183, M_c,apo = 60.668.
- Start of φ: χ₀ = instantaneous minimum of V near χ = 0 at the orbit's starting point (Newton from
  3ε_t M₀³/(g² r² − λ_P M₀² + 6ε_t M₀²)), χ′₀ = 0. χ₀/M_c,peri ranges from 4×10⁻⁴ (apocentre) to 3×10⁻² (pericentre).
  Also tested: χ₀ = 0 and χ₀ = 10⁻⁶ M_c,peri.
- Design: uniform random phases (fraction of T_r; fixed rng, the same phases for every ε and orbit), 32 runs per
  (orbit, ε) for ε = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵, plus 16 phase-pair runs (k/16 and k/16 + 10⁻³) at 10⁻³ and 10⁻⁴, 8 pairs
  at 10⁻⁵, threshold points 3×10⁻⁴ (late) and 3×10⁻⁵ (both), and ε = 10⁻⁶ (2–3 runs to 40 M_c,apo, 6 more to
  10 M_c,apo; n₋ at 10 and 40 M_c,apo differ by ≤ 7×10⁻⁴). Runtime 0.1 s (10⁻²) to 290 s (10⁻⁶) per run.
- Escape diagnostics per run: first grid M with |χ| > 0.1 M ("escape", as in pt_matching), first M with χ > 0.5 M, and
  the "settle" M after which χ > 0.5 M for ever. The grid has 500 log-spaced points from M₀ to M_end.
- Controls: (a) ε_t = 0 in the same code (with the Z₂ vacuum ±M in the measure) against the skeptic's untilted numbers;
  (b) the skeptic's linear-tilt model (mode 1, β = 0.2, vacuum x = 1.08803 M, m_φ/m_X = 2.543) to separate the linear
  tilt from the extra quadratic/cubic terms; (c) J_r/q = 1.7637 / 0.4805 as the adiabatic reference.
- Integration checks: E(s) − E(0) − W(s) relative to |E| ("W" below) and max |q/q₀ − 1| ("q"), and rtol 10⁻¹²/10⁻¹⁰
  reruns.

All numbers: `results/m1_transition.json` (every group, check and control) and `results/m1_frequency_lock.json`
(diagnostics along M). Per-run files (767, with traces) are in `$SP/review/pt_m1tilt/runs/`.

## 3. Results

### 3.1 Main table (M1 potential, χ₀ = instantaneous minimum, rtol 10⁻¹¹, M₀ = 0.5 M_c,peri, window at 40 M_c,apo)

**τ-20 orbit, J_r/q = 1.76373** (random and pair designs pooled; the outcome is deterministic, see 3.2)

| ε | runs | n₋/q mean | sd | min | max | 1 − n₋/J_r | N_φ/q mean | sd | min | max | n₋ + N_φ | W | q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10⁻² | 32 | 1.76333 | 0.00015 | 1.76316 | 1.76361 | 2.3×10⁻⁴ | 3.7×10⁻⁴ | 0.8×10⁻⁴ | 2.1×10⁻⁴ | 4.5×10⁻⁴ | 1.7637 | 1×10⁻¹¹ | 2×10⁻¹¹ |
| 10⁻³ | 64 | 1.76369 | 0.00002 | 1.76366 | 1.76371 | 2.2×10⁻⁵ | 3.8×10⁻⁵ | 1.3×10⁻⁵ | 1.9×10⁻⁵ | 5.9×10⁻⁵ | 1.7637 | 2×10⁻¹⁰ | 2×10⁻¹⁰ |
| 10⁻⁴ | 64 | 1.76350 | 0.00015 | 1.76331 | 1.76372 | 1.3×10⁻⁴ | 1.6×10⁻⁴ | 1.0×10⁻⁴ | 1.5×10⁻⁵ | 2.9×10⁻⁴ | 1.7637 | 2×10⁻⁹ | 3×10⁻⁹ |
| 3×10⁻⁵ | 32 | 1.76011 | 0.00527 | 1.74417 | 1.76373 | 2.1×10⁻³ | 2.4×10⁻³ | 3.5×10⁻³ | 1.3×10⁻⁵ | 1.3×10⁻² | 1.7625 | 7×10⁻⁹ | 8×10⁻⁹ |
| **10⁻⁵** | 48 | **1.5202** | 0.0122 | 1.4982 | 1.5639 | **0.138** | **0.162** | 0.008 | 0.13 | 0.18 | 1.683 | 3×10⁻⁷ | 1×10⁻⁸ |
| **10⁻⁶** | 2 (+6 at 10 M_c,apo) | **1.4888** (pooled 8: 1.4904) | 0.0012 (0.0036) | 1.4880 (1.4875) | 1.4896 (1.4985) | **0.156** | **0.183** | 0.001 (0.002) | | | 1.672 | 3×10⁻⁶ | 1×10⁻⁷ |

**Late orbit, J_r/q = 0.48049**

| ε | runs | n₋/q mean | sd | min | max | 1 − n₋/J_r | N_φ/q mean | sd | min | max | n₋ + N_φ | W | q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10⁻² | 32 | 0.47994 | 0.00070 | 0.47882 | 0.48076 | 1.2×10⁻³ | 4.5×10⁻⁴ | 3.7×10⁻⁴ | 2×10⁻⁵ | 1.0×10⁻³ | 0.4804 | 1×10⁻¹¹ | 3×10⁻¹¹ |
| 10⁻³ | 64 | 0.48017 | 0.00028 | 0.47974 | 0.48052 | 6.7×10⁻⁴ | 2.5×10⁻⁴ | 1.8×10⁻⁴ | 2×10⁻⁵ | 5.3×10⁻⁴ | 0.4804 | 2×10⁻¹⁰ | 2×10⁻¹⁰ |
| 3×10⁻⁴ | 32 | 0.47350 | 0.00502 | 0.46551 | 0.48037 | 1.5×10⁻² | 4.7×10⁻³ | 3.3×10⁻³ | 1×10⁻⁴ | 1.0×10⁻² | 0.4782 | 1×10⁻⁹ | 5×10⁻¹⁰ |
| **10⁻⁴** | 64 | **0.2231** | 0.0449 | 0.1630 | 0.4108 | **0.536** | **0.172** | 0.030 | 0.046 | 0.21 | 0.395 | 8×10⁻⁸ | 2×10⁻⁹ |
| 3×10⁻⁵ | 16 | 0.1785 | 0.0115 | 0.1597 | 0.2059 | 0.629 | 0.200 | 0.007 | 0.18 | 0.21 | 0.379 | 3×10⁻⁷ | 6×10⁻⁹ |
| **10⁻⁵** | 48 | **0.1464** | 0.0276 | 0.1165 | 0.2527 | **0.695** | **0.221** | 0.016 | 0.16 | 0.24 | 0.368 | 1×10⁻⁶ | 2×10⁻⁸ |
| **10⁻⁶** | 3 (+6 at 10 M_c,apo) | **0.1088** (pooled 9: 0.1243) | 0.0094 (0.0292) | 0.1015 (0.0839) | 0.1194 (0.1729) | **0.774** | **0.245** | 0.006 (0.017) | | | 0.353 | 1×10⁻⁵ | 2×10⁻⁷ |

Reading: **the M1 transition is adiabatic for ε ≥ 10⁻⁴ (τ-20) and ε ≥ 10⁻³ (late), to ≤ 2×10⁻⁴ in J_r, but below a
sharp threshold (between 3×10⁻⁵ and 10⁻⁵ for τ-20, between 3×10⁻⁴ and 10⁻⁴ for late) the slower the transition, the
more radial action it loses**, deterministically, and the loss is still creeping up at 10⁻⁶. The ε → 0 limit is
therefore not J_r.

- Work identity: ≤ 2×10⁻⁹ relative to |E| for ε ≥ 10⁻⁴ (2×10⁻¹⁰ at 10⁻³; 2×10⁻¹¹ with rtol 10⁻¹²); it grows with the
  step count to 10⁻⁵ at ε = 10⁻⁶ (4×10⁸ steps). Charge likewise. The rtol reruns (3.2) show the results do not depend
  on this.
- n₋ with the actual mass g|χ| instead of gM differs by < 10⁻⁵ at 40 M_c,apo. n₋ at 10 vs 40 M_c,apo differ by
  ≤ 7×10⁻⁴ (τ-20) and 2×10⁻⁴ (late).

### 3.2 Deterministic or chaotic? Phase pairs, seed perturbation, rtol, χ₀, M₀

| Test | τ-20 | late |
|---|---|---|
| neighbouring phases (Δphase = 10⁻³), max |Δn₋/q| | 3.5×10⁻⁷ (10⁻³), 1.4×10⁻⁶ (10⁻⁴), 1.5×10⁻² (10⁻⁵; n₋ varies smoothly from 1.498 to 1.564 over the phase) | 5.2×10⁻⁶ (10⁻³), 1.1×10⁻² (10⁻⁴; n₋ = 0.163–0.386 over the phase), 3.7×10⁻² (10⁻⁵), 1.4×10⁻² (10⁻⁶, one pair) |
| χ₀ × (1 + 10⁻⁹) | 4×10⁻¹² (10⁻³), 1.7×10⁻¹¹ (10⁻⁴), 2.1×10⁻⁷ (10⁻⁵) | 6×10⁻¹³ (10⁻³), 4.3×10⁻¹⁰ (10⁻⁴), **1.3×10⁻²** (10⁻⁵) |
| rtol 10⁻¹² vs 10⁻¹¹ | 1.6×10⁻⁹ (10⁻³), 5.2×10⁻⁴ (10⁻⁵) | 4.8×10⁻¹⁰ (10⁻³), 4.5×10⁻⁶ (10⁻⁴); rtol 10⁻¹⁰: 5.2×10⁻⁵ (10⁻⁴) |
| χ₀ = 0 or 10⁻⁶ M_c,peri instead of the minimum (ε = 10⁻³, 32 phases each) | ≤ 1.8×10⁻⁵ | – |
| M₀ = 0.25 M_c,peri instead of 0.5 | 10⁻⁵: 1.519 ± 0.014 (8) vs 1.520 ± 0.012 | 10⁻³: 0.4805 vs 0.4802; **10⁻⁴: 0.388 ± 0.048 (16) vs 0.223 ± 0.045**; 10⁻⁵: 0.142 ± 0.010 (8) vs 0.146 ± 0.028 |

So:
- At every ε the outcome is a **smooth, deterministic function of the phase**: nothing like the untilted broad split
  (τ-20: 0.016–1.31 at the same ε = 10⁻³, control (a)). The late orbit at 10⁻⁵ shows a weak sensitivity (a 10⁻⁹ seed
  change moves n₋ by 10⁻²), but the distribution stays narrow (0.146 ± 0.028).
- The lossy result does not come from the integrator (rtol) or the φ start (χ₀). It is insensitive to the start of the
  ramp deep in the lossy regime (M₀, 10⁻⁵), but **at the threshold ε** (late, 10⁻⁴) the start matters (0.39 vs 0.22).
- The escape diagnostics are not informative for a tilted potential: |χ| > 0.1 M at the first or second grid point
  (0.51–0.64 M_c,peri) because the tilt pushes φ out at once; χ > 0.5 M first at 0.81–0.93 M_c,peri; φ settles on the
  +M side at **7.1–7.5 M_c,peri (τ-20) and 2.9–3.1 M_c,peri (late)**, independent of ε. (Untilted control: escape at
  0.98 / 0.89 M_c,peri, as the skeptic's 0.999 / 0.897.)

### 3.3 Controls

**(a) ε_t = 0 in the same code** (ε = 10⁻³, 16 random phases, seed 10⁻⁶ M_c,peri; vacuum ±M in the measure):
τ-20 n₋/q = 0.618 ± 0.104 (sd 0.42, range 0.016–1.31), N_φ/q = 0.54 ± 0.06; late 0.169 ± 0.029 (sd 0.11, range
0.004–0.40), N_φ/q = 0.143 ± 0.015. The skeptic: 0.704 ± 0.017 / 0.475 (560 runs) and 0.197 ± 0.004 / 0.124; the
track: 0.70 ± 0.02 / 0.47. Consistent within the 16-run error. Run-by-run, mode 3 with ε_t = 0 reproduces the skeptic's
deterministic ε = 10⁻² runs to 2×10⁻¹³ in n₋ (same phase and seed).

**(b) Linear tilt only, β = 0.2** (skeptic's model; 16 runs at 10⁻³ and 10⁻⁴, 8 at 10⁻⁵):

| ε | τ-20 n₋/q | N_φ/q | late n₋/q | N_φ/q |
|---|---|---|---|---|
| 10⁻³ | 1.76368 ± 0.00001 | 4×10⁻⁵ | 0.48031 ± 0.00003 | 2×10⁻⁴ |
| 10⁻⁴ | 1.76365 ± 0.00005 | 6×10⁻⁵ | 0.48025 ± 0.00011 | 2×10⁻⁴ |
| 10⁻⁵ | 1.7619 ± 0.0012 [1.7609, 1.7637] | 1.2×10⁻³ | 0.4774 ± 0.0031 [0.4719, 0.4805] | 2×10⁻³ |

The linear piece of the M1 tilt alone keeps J_r to ≤ 10⁻³ (≤ 6×10⁻³ for late) down to ε = 10⁻⁵, where the full M1
potential loses 14% / 70%. **The small-ε loss is produced by the extra terms of the M1 shape** (+3ε_t M² χ², i.e. the
0.6× tachyonic quadratic, and −ε_t M χ³), not by the linear tilt. Which of the two is decisive was not separated.

**(c) Adiabatic reference:** J_r/q = 1.76373 (τ-20), 0.48049 (late). Reached to 2×10⁻⁴ at ε ≥ 10⁻⁴ / 10⁻³.

### 3.4 Where and how the action is lost (diagnostics; `results/m1_frequency_lock.json`)

- **The φ swing is large and ε-independent.** From the start the tilt drives φ into a large forced oscillation: the
  rms of |χ/M − 1| over the grid is 0.8–0.9 from 0.5 to 5 M_c,peri (τ-20; to 2 M_c,peri for late) with excursions
  χ/M ∈ [−0.7, 2.7]. Binned in M, the profile is the same at ε = 10⁻², 10⁻³, 10⁻⁴ and 10⁻⁵ to a few per cent. The
  β = 0.2 model has a similar swing (rms 0.9–1.0).
- **It is slaved.** In 24-radial-period windows the dominant Fourier frequency of χ equals that of r² to 10⁻⁴ and
  ≥ 99.9% of the χ power sits at harmonics of the radial frequency, for both orbits and both regimes, up to
  ≈ 5.6 M_c,peri (τ-20) / 2.0 M_c,peri (late). The measured radial frequency in this regime is far below the
  frozen-orbit estimate (e.g. 1.16 vs 1.9 per unit s at 5 M_c,peri), i.e. the X orbit is strongly modified by φ.
- **The loss happens where the slaved regime ends.** In the lossy τ-20 run (ε = 10⁻⁵) a free φ component at the linear
  m_φ appears at 5.6–7 M_c,peri (harmonic power fraction 1.00 → 0.96 → 0.87 → 0.78 at 5.6, 6.8, 8.2 M_c,peri; the
  dominant χ frequency switches to m_φ by 13 M_c,peri and the ratio to the radial frequency tends to
  m_φ/2m_X = 1.2247). In the adiabatic run (10⁻⁴) the fraction stays ≥ 0.998 through the same M and the free part
  decays to N_φ/q ~ 10⁻⁴. For the late orbit the release is at 2.2–3 M_c,peri (fraction 1.00 → 0.25 at 2.19 at
  ε = 10⁻⁴; ≥ 0.99 at 10⁻³). Both sites coincide with the collapse of the swing (rms δχ/M from 0.7 to 0.3) and with
  the "settle" M of 3.2. The window-averaged n₋ and N_φ are flat after ≈ 10 M_c,peri (τ-20) / 4 M_c,peri (late).
- For τ-20 the linear ratio m_φ/ω_r (frozen orbit) rises to ≈ 2.0 at 6.2 M_c,peri and falls back to 1.22; for late it
  peaks at ≈ 1.8 at 2.2 M_c,peri. A grazing 1:2 commensurability is therefore a candidate for τ-20 but not a uniform
  explanation; the sharp ε threshold, the deterministic phase-smooth outcome and the slow saturation below threshold are
  the signatures of capture into a φ–orbit resonance (autoresonance-like) rather than of a single fast crossing (which
  would give a loss that vanishes with ε) or of chaotic diffusion. **The mechanism is characterised, not derived.**

## 4. Verdict

Under Codex's own M1 potential, in the homogeneous classical X + φ model with the slow drift M ∝ e^{εω_r s} [C]:

1. **The chaotic random split of the untilted model is gone.** φ is pushed to +M from the start (no Floquet-unstable
   symmetric state to escape from), the outcome is a deterministic, phase-smooth function, and for ε ≥ 10⁻⁴ (τ-20)
   or ≥ 10⁻³ (late) the transition returns J_r to ≤ 2×10⁻⁴: n₋/q = 1.7637 and 0.4805 with N_φ/q ≤ 4×10⁻⁴.
2. **But the slow limit is not adiabatic.** Below a sharp threshold (3×10⁻⁵–10⁻⁵ for τ-20, 3×10⁻⁴–10⁻⁴ for late) the
   transition deterministically converts radial action into free φ quanta at the end of the slaved φ regime, and the
   loss grows as ε falls: τ-20 n₋/q = 1.520 ± 0.012 (10⁻⁵) → 1.489–1.490 (10⁻⁶), N_φ/q = 0.162 → 0.183; late
   0.146 ± 0.028 (10⁻⁵) → 0.109–0.124 (10⁻⁶), N_φ/q = 0.221 → 0.245. Codex's premise "if the actual continuation is
   adiabatic" fails in this model also for M1, by a different mechanism from the untilted case.
3. **Resulting coherent-orbit numbers per orbit charge q, taken at ε = 10⁻⁵–10⁻⁶** (the GRUT value is 6–9 decades
   smaller, so these are the slow-limit values of this scan, not converged limits [E]):
   - τ-20 orbit: **n₋/q ≈ 1.49–1.52 X̄ (0.85 J_r) plus 0.16–0.18 φ quanta** (sum 1.67–1.68; if φ → X X̄ each quantum
     adds one X̄, giving ≈ 1.67–1.70), against J_r = 1.76 and the untilted 0.70 + 0.47.
   - late orbit: **n₋/q ≈ 0.11–0.15 (0.23–0.30 J_r) plus 0.22–0.25 φ quanta** (sum 0.35–0.37), against J_r = 0.48 and
     the untilted 0.20 + 0.12.
4. The linear tilt alone (β = 0.2) would conserve J_r to 10⁻³ at all ε tested; the loss is specific to the M1 shape.

## 5. Limits

- One homogeneous coherent mode; φ condenses as a zero mode; classical; no gauge, Y or vector fields; the late orbit is
  a pure-quartic stand-in; the same drift law and start as pt_matching (the start matters only at the threshold ε).
- ε = 10⁻⁶ rests on 8–9 runs per orbit, and the loss is still rising there; the extrapolation to ε ~ 10⁻¹²–10⁻¹⁵ is an
  extrapolation.
- The resonance responsible is characterised (site, slaved-to-free transition, threshold, determinism) but not
  identified analytically; the quadratic and cubic extra terms were not separated.
- Work identity 10⁻⁵ at ε = 10⁻⁶ (step-count limited); rtol reruns bound the effect on n₋ at ≤ 5×10⁻⁴.

## Reproduce

```bash
cd research/reality_framework/grut_formation_review_20261025/tracks/pt_m1tilt/code
for o in tau20 late; do for e in 1e-2 1e-3 1e-4 1e-5; do python3 -B m1_scan.py $o $e rand32; done; done
python3 -B m1_scan.py late 3e-4 rand32; python3 -B m1_scan.py tau20 3e-5 rand32; python3 -B m1_scan.py late 3e-5 rand16
for o in tau20 late; do for e in 1e-3 1e-4; do python3 -B m1_scan.py $o $e pair16; python3 -B m1_scan.py $o $e rand4 --pert 1e-9; done
  python3 -B m1_scan.py $o 1e-5 pair8; python3 -B m1_scan.py $o 1e-5 rand2 --pert 1e-9; python3 -B m1_scan.py $o 1e-3 rand4 --rtol 1e-12; done
python3 -B m1_scan.py late 1e-4 rand4 --rtol 1e-12; python3 -B m1_scan.py late 1e-4 rand4 --rtol 1e-10; python3 -B m1_scan.py tau20 1e-5 rand2 --rtol 1e-12
python3 -B m1_scan.py tau20 1e-3 rand32 --phi0 0; python3 -B m1_scan.py tau20 1e-3 rand32 --phi0 1e-6
for o in tau20 late; do python3 -B m1_scan.py $o 1e-5 rand8 --mstart 0.25; done; python3 -B m1_scan.py late 1e-4 rand16 --mstart 0.25; python3 -B m1_scan.py late 1e-3 rand8 --mstart 0.25
python3 -B m1_scan.py tau20 1e-6 rand2; python3 -B m1_scan.py late 1e-6 rand1; python3 -B m1_scan.py late 1e-6 pair1
for o in tau20 late; do python3 -B m1_scan.py $o 1e-6 rand6 --mend 10; done
for o in tau20 late; do python3 -B m1_scan.py $o 1e-3 rand16 --epst 0 --phi0 1e-6
  for e in 1e-3 1e-4; do python3 -B m1_scan.py $o $e rand16 --mode 1 --beta 0.2 --epst 0 --phi0 1e-6; done
  python3 -B m1_scan.py $o 1e-5 rand8 --mode 1 --beta 0.2 --epst 0 --phi0 1e-6; done
python3 -B m1_diag.py tau20 1e-5 0.3; python3 -B m1_diag.py late 1e-4 0.3
for a in "tau20 1e-5" "tau20 1e-4" "late 1e-4" "late 1e-3"; do python3 -B m1_freq.py $a 0.3; done
python3 -B m1_aggregate.py
```

The C library compiles automatically (gcc -O2) into `$SP/review/pt_m1tilt/`. Two worker processes; about 2.5 h in all,
dominated by ε ≤ 10⁻⁵.
