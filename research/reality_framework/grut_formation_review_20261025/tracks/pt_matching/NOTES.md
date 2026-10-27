# Track pt_matching: coherent matching and whether the binding transition keeps the radial action

27 October 2026. Claude review of Codex 5f5513b, issue #13. Briefs: BRIEF_TRANSPORT.md and BRIEF.md.

All numbers here come from **independent code**: no Codex module was imported or run.
- The integrator is our own: a position-Verlet map in extended phase space, composed to order 8 (Kahan–Li) and, for
  checks, order 6 (Yoshida). Each step depends only on s, so every step is a symplectic map.
- Quadrature uses mpmath tanh–sinh in a variable where the endpoint factors cancel analytically. FFTs use numpy.

The only inputs are:
- the τ = 20 light mean field, `MATCH_L12.npz` `light_mean_q/p`;
- the coherent share quoted from Codex's `WINDOW_base48_24_34.json`.

This is one supplied, reduced, classical model. Nothing here is a validation, an exclusion or new theory.

## 1. The starting orbit and the matching claims [C, D]

**Orbit** (`results/orbit_spectral.json`):

| Quantity | Ours | Codex |
|---|---|---|
| q, E | 8117.881574752555, 105233.42336942146 | same |
| r_min, r_max | 17.697074113681, 142.682453672177 | same |
| T_r, Δθ | 0.81446568177456, 2.99558926699110 | same |
| J_r/q (quadrature) | 1.7637252793970 | 1.763725279396997 |
| A₋/q, the lab negative-frequency statistic (exact quasi-periodic Fourier series, no window) | 1.7377883531122 | 1.7377883531120 |
| J_r and A₋ per original net charge (8113.477) | 1.7646827365883, 1.7387317301648 | 1.7646827, 1.7387317 |

Other checks on the orbit:
- Ω/ω_r = 0.476763, where ω_r = 7.71449.
- Orbit closure is 6×10⁻¹⁴.
- Step convergence: N = 1024 to 8192 steps per period agree to 10⁻¹³.

**Why the two differ: an exact identity [D].** Write ψ = χ₁ + iχ₂ = Σₙ cₙ e^{i(Ω + nω_r)s}, with ωₙ = Ω + nω_r. Here Ω is
the mean angular velocity, so ψe^{−iΩs} is periodic with T_r. Evaluate the action integrals on the invariant torus,
using (1/2π)∮ p·dq = ⟨Re(ψ̇* ∂ψ/∂θᵢ)⟩. This gives

  q = Σₙ ωₙ|cₙ|²,  J_r = Σₙ n ωₙ|cₙ|²   (Percival's spectral form of the actions [R]).

When 0 < Ω < ω_r, ωₙ < 0 exactly when n ≤ −1. Then

  J_r − A₋ = Σ_{n≥1} n ωₙ|cₙ|² + Σ_{n≤−2} (|n| − 1)|ωₙ||cₙ|² ≥ 0.

This vanishes only for the two-harmonic orbit c₀e^{iΩs} + c₋₁e^{i(Ω−ω_r)s}. That is the isotropic harmonic ellipse,
where ω_r = 2Ω = 2m and E/m = 2J_r + q, so n₋ = (E/m − q)/2 = m|B|² = J_r = A₋.
- **Numerically:** the identities hold to 2×10⁻¹⁴.
- **Quartic orbit:** J_r − A₋ = 0.025937 q. Of this, 0.016534 comes from n ≥ 1 (mostly n = 1 at ω = 11.39) and 0.009403
  from n ≤ −2 (mostly n = −2 at ω = −11.75).
- **Interpolating family** (same initial state, potential (1−f)·m²r²/2 + f·λ_X r⁴/4):

  | Quartic fraction f | A₋/q | J_r/q |
  |---|---|---|
  | 0 (harmonic) | 1.98238 | 1.98238 |
  | 0.5 | 1.84509 | 1.85467 |
  | 1 (quartic) | 1.73779 | 1.76373 |

- So A₋ is not an adiabatic invariant. A slow quartic-to-harmonic deformation keeps J_r and turns it into A₋ = n₋.

**Ramps** (`results/matching_ramps.json`). Codex's radial-phase convention is states at p·T_r after the τ = 20 state.

| Duration T | Ours, 8 phases | Codex |
|---|---|---|
| 0 (sudden) | 0.0769391–18.175696 | 0.0769391–18.175696 |
| 0.2 | 0.110720–8.73811 | same |
| 2 | 1.24548–2.33108 | same |
| 20 | spread 3.5612×10⁻⁴ | 3.5612×10⁻⁴ |
| 80 | mean 1.76372527936, spread 1.99889×10⁻⁵ | 1.76372527845, 1.9989×10⁻⁵ |

- **Mass-only control** (λ fixed, T = 80): the final J_r/q − J₀ is −4.80×10⁻⁷, +1.19×10⁻⁷, +4.79×10⁻⁷ and −1.18×10⁻⁷.
  This is Codex's 5×10⁻⁷.
- **Sudden limit over a continuous phase** (512 phases): 0.00535 (phase 0.773) to 18.312 (phase 0.273), with mean 8.42.
  Codex's "0.077–18.18" is the 8-phase sample.
- **Dependence on T** (16 phases):
  - Codex's cosine ramp: max|n₋ − J_r|/q = 9.4×10⁻⁴, 1.9×10⁻⁴, 4.7×10⁻⁵, 1.1×10⁻⁵, 2.8×10⁻⁶ and 7.1×10⁻⁷ at
    T = 10, 20, 40, 80, 160 and 320. That is a power law **T^−2.02**, set by the jump in f″ at the ends of the ramp.
  - The phase-mean error is ≤ 6×10⁻⁹ for T ≥ 20.
  - A C^∞ ramp gives 3.0×10⁻³, 9.7×10⁻⁷, 2.7×10⁻¹⁰ and 4×10⁻¹³ (round-off) at T = 10, 20, 40 and 80, faster than
    any low power.
  - So the residual comes from the ramp's endpoints, not from the matching.
- Work and charge identities hold to ≤ 10⁻¹¹. Halving the step changes results by ≤ 1.2×10⁻¹².

## 2. The physical transition: φ condensation [C/E]

**Model.** The conformal radiation-era potential (a″ = 0) is

  U = λ_φφ⁴/4 − λ_φM(s)²φ²/2 + g²φ²r²/2 + λ_X r⁴/4,  with M = v·a(s),

and λ_X = 10⁻³, g² = 10⁻², λ_φ = 0.03. Late on, φ → M, m_X = gM and m_φ = √(2λ_φ)M = 2.449 m_X.
- The drift is M = M_start·e^{κs}, with κ = ε·ω_r(initial orbit). We also ran the radiation-era form M ∝ 1 + κs.
- Start: M_start = 0.5·M_c,peri, where M_c,peri = g·r_min/√λ_φ. The X state is the orbit at phase p; φ = σ·M_c,peri,
  with φ′ = 0.
- End: M_end = 40·M_c,apo. There, the X field is a nearly free oscillator with an anharmonic correction of about 10⁻⁴.
- **Measures:** n₋/q is the window average of [(p_X²/2 + g²M²r²/2)/(gM) − q]/(2q) over 60 X periods.
  N_φ/q = ⟨(δφ′²/2 + m_φ²δφ²/2)/m_φ⟩/q, the classical number of φ quanta. n₋ is converged by M ≈ 5·M_c,apo (to
  0.01 per run).
- **Scaling.** The system is scale-covariant, so only the orbit shape (J_r/q), ε, phase and ln σ matter.

**Orbits.**
- **τ = 20:** J_r/q = 1.7637. M_c,peri = 10.217 and M_c,apo = 82.378.
- **Late orbit:** pure-quartic with the same q, matched to Codex's coherent share in s 24–34. That share is
  A₋ = 0.47604 per net charge, divided by mean Q₀/Q = 1.00907, which gives 0.47177 per coherent charge. The orbit then
  has J_r/q = 0.48049, r_max/r_min = 3.34, M_c,peri = 18.18 and M_c,apo = 60.67.
- **Comparison:** matching instead the sampled r_max/r_min = 3.96 of the base48 coherent trace gives J_r/q = 0.646. That
  trace feels the Hartree mass of about 20, so the two matchings differ.

**Floquet stability of φ = 0** (`results/floquet.json`; cross-checked with scipy DOP853 to 10⁻¹²).
- φ = 0 first becomes unstable at **M = 0.732·M_c,peri for both orbits**. This is not at the static pericentre
  pitchfork. It is a parametric band driven by the orbit (multiplier +, phase ≈ 2π per radial period), and it merges
  into the pericentre tachyonic region.
- For τ = 20, φ = 0 stays unstable up to 3.79·M_peri. ν/ω_r = 0.05 at 0.87 M_peri and 0.15 at 1.75 M_peri.
- So the onset is orbit-pumped, and in particle language a coherent X X̄ → φφ channel.

**Scan** (`results/transition_scan.json`): τ = 20 orbit, 16 phases × seeds σ = 10⁻², 10⁻⁶, 10⁻¹⁰.

| ε | n₋/q mean ± sem [min, max] | N_φ/q | n₋ + N_φ | M_escape/M_peri |
|---|---|---|---|---|
| 10⁻¹ | 6.50 ± 0.80 [0.19, 20.8] | 6.93 | 13.4 | 12.3 |
| 3×10⁻² | 1.15 ± 0.14 [0.005, 4.51] | 1.16 | 2.30 | 5.6 |
| 10⁻² | 0.735 ± 0.077 [0.043, 1.81] | 0.544 | 1.28 | 2.39 |
| 3×10⁻³ | 0.716 ± 0.064 [0.008, 1.44] | 0.463 | 1.18 | 1.27 |
| 10⁻³ | 0.697 ± 0.061 [0.002, 1.37] | 0.473 | 1.17 | 0.95 |
| 3×10⁻⁴ | 0.756 ± 0.057 [0.048, 1.39] | 0.443 | 1.20 | 0.82 |
| 10⁻⁴ | 0.728 ± 0.056 [0.047, 1.38] | 0.459 | 1.19 | 0.77 |
| 3×10⁻⁵ | 0.726 ± 0.059 [0.005, 1.41] | 0.459 | 1.18 | 0.75 |

Pooled over ε ≤ 3×10⁻³ (240 runs):
- **n₋/q = 0.724 ± 0.026** (std 0.41; 10–90% range 0.16–1.28; never above J₀ = 1.764). That is 0.41 J₀.
- **N_φ/q = 0.459 ± 0.017**.
- n₋ + N_φ = 1.184 ± 0.009.
- Trend: slope −0.011 ± 0.037 per decade of ε.
- Seeds: 0.66 / 0.73 / 0.78 ± 0.045 for σ = 10⁻² / 10⁻⁶ / 10⁻¹⁰.
- The φ oscillation holds about 32% of the late coherent energy.
- The escape point tends to the Floquet edge (0.750 at 3×10⁻⁵, edge 0.732).

**Late orbit** (ε = 10⁻³ and 10⁻⁴, 96 runs):
- n₋/q = 0.166 ± 0.011 (range 0.0015–0.363), against J₀ = 0.4805.
- N_φ/q = 0.141 ± 0.006; the sum is 0.306 ± 0.005.
- ε = 10⁻²: 0.200.

**Controls** (`transition_controls.json`, `transition_chaos.json`):
- **Imposed mass g·M(s) with no φ:** J_r is kept to 4×10⁻⁶ at ε = 10⁻⁴ (error ∝ ε from the abrupt start of the
  exponential).
- **Born–Oppenheimer φ** (φ at its instantaneous minimum, which gives the attractive effective quartic
  λ_X − g⁴/λ_φ < 0 inside r < √3·M): J_r is kept to ≤ 5.3×10⁻⁴ for both orbits.
- **Step halving** at ε = 10⁻³: 0.741 ± 0.059.
- **Radiation-era growth** (ε = 10⁻³ at M_peri): 0.618 ± 0.065.
- **Chaos test:** changing the φ seed by a relative 10⁻¹² changes n₋/q by up to 0.96 at ε ≤ 10⁻³, but by only ≤ 2×10⁻¹²
  at ε = 10⁻².
- Charge holds to ≤ 7×10⁻¹² and the work identity to ≤ 3×10⁻¹¹ in every run.

**Reading.** In this model the radial action is **not conserved** through the transition.
- The loss comes neither from the mass switch-on nor from the binding effective potential: both controls keep J_r. It
  comes from φ's finite inertia. m_φ is comparable to the orbit frequency, so the X–φ system is a strongly coupled,
  chaotic two-degree-of-freedom problem.
- The change does not vanish as ε → 0, unlike a single separatrix crossing (Neishtadt; Cary–Escande–Tennyson; the
  jump there is O(ε ln ε) [R]). Instead it tends to an ε-independent distribution, as in the ergodic-adiabatic picture
  (Ott 1979 [R]).
- **Extrapolating to ε ~ 10⁻¹⁵ [E]:** the distribution is flat from 3×10⁻³ to 3×10⁻⁵ and escape sits at the Floquet
  edge. This assumes no new slow mechanism below 3×10⁻⁵.
- **Fate of the φ quanta [E].** m_φ = 2.449 m_X > 2m_X, so φ → X X̄ is open:
  Γ ≈ g⁴v²√(1 − 4m_X²/m_φ²)/(4πm_φ) ≈ 1.9×10⁻⁴ m_X, so Γ/H ~ 2×10¹¹.
  - If φ has no faster channel, each φ quantum adds one X̄ (and one X), as waves of p ≈ 0.71 m_X.
  - The coherent channel would then give about **1.18 X̄ per coherent charge** (τ = 20 orbit) or **0.31** (late orbit),
    against J_r = 1.76 or 0.48.

**Waves** (`waves_switch_on.json`; linear test modes, m_w² = 20 Hartree [Q], k = 0.25–8, 8 phases):
- **Smooth envelope g·M(s):** |β_k|² ≤ 1.1×10⁻⁵, which is adiabatic and consistent with H/ω ~ ε.
- **Actual φ(s) of the homogeneous runs**, at M = 2·M_peri:
  - |β|² ≈ 10⁻³–3×10⁻² for most k;
  - the k = 6 band (ω ≈ ω_r) reaches 1.9 at ε = 10⁻³ and 300 at 10⁻⁴, growing as ε falls;
  - by 4–8·M_peri, |β|² ≫ 1 for every k (overflow in linear modes with no backreaction).
- **So the wave mass switch-on is not adiabatic** when φ condenses through the orbit-driven resonance. The wave number
  then needs a 3-D lattice treatment through the transition (cf. preheating: Kofman–Linde–Starobinsky; Greene et al.
  [R]).

## 3. Verdict on "slow matching to massive particles gives about 1.76; no extra counting postulate is needed"

- **Counting rule: CONFIRMED [D].** In a free massive basis, coherent counter-rotation is n₋ = m|B|².
- **Matching numbers: CONFIRMED [C].** 1.737788 and 1.763725 reproduce, and the slow external ramps return J_r.
- **The number as a physical statement: CORRECTED [C].** In the homogeneous X + φ model, the GRUT mass appears through φ
  condensation, and that transition does not continue adiabatically.
  - It gives 0.72 ± 0.03 coherent X̄ per charge (0.41 J_r, broad and chaotic), plus 0.46 φ quanta. For the late orbit it
    gives 0.17 + 0.14.
  - That is about 1.18 or 0.31 if φ → X X̄, not 1.76 or 0.48.
- Codex's own caveat ("if the actual continuation is adiabatic") is the premise that fails here.

**Limits.**
- One homogeneous coherent mode; φ condenses as a zero mode. Real condensation is inhomogeneous: spinodal growth, and
  Z₂ domains unless φ is biased.
- Classical; no gauge, Y or vector fields; m_w² is held fixed in the wave test.
- The late orbit is a pure-quartic stand-in.
- The ε → 0 statement is an extrapolation over 10 decades.

## Reproduce

```bash
cd research/reality_framework/grut_formation_review_20261025/tracks/pt_matching/code
python3 -B part1_orbit.py; python3 -B part1_ramps.py; python3 -B part2_floquet.py
for e in 1e-1 3e-2 1e-2 3e-3 1e-3 3e-4 1e-4 3e-5; do python3 -B part2_scan.py tau20 $e 16; done   # 3e-5: ~16 min, 2 procs
for e in 1e-1 1e-2 1e-3 1e-4; do python3 -B part2_scan.py late $e 16; done
H0=0.04 python3 -B part2_scan.py tau20 1e-3 16; GROWTH=1 MEND=5 python3 -B part2_scan.py tau20 2e-3 16
python3 -B part2_controls.py; python3 -B part2_chaos.py; python3 -B part2_waves.py; python3 -B part2_aggregate.py save
```

The C library compiles automatically (gcc -O3) into `$SP/review/pt_matching/`. Per-run scan files (with traces) are in
`$SP/review/pt_matching/scan/`.
