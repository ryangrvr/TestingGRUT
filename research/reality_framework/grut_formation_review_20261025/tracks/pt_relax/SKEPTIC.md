# pt_relax skeptic: does the "~1.75 X̄ budget relaxes before binding" case hold up? (Codex 5f5513b, issue #13)

Claude adversarial review of track `pt_relax`, 1 November 2026 (started 27 October; resumed after a usage limit).
Brief: [BRIEF_TRANSPORT.md](../../BRIEF_TRANSPORT.md). This is same-session AI work on one supplied, reduced scalar
model. It is not outside replication, and nothing here is observational evidence, a validation, an exclusion or a
GRUT prediction. Nothing was committed and Codex's tree was not edited.

## Labels

- **Reproduction tool:** Codex's FD4/Yoshida engine (`fieldstep.cpp` via `pt_common`) and `prepare.initialize`, run
  from a private `git archive 5f5513b` copy in `$SP/review/pt_relax_skeptic/rf/`. pt_relax cross-checked this engine
  against our own BM6 integrator (2×10⁻⁸ on the windowed action), so I did not repeat that.
- **Own code** (`skeptic/code/`, written fresh, importing nothing from `tracks/pt_relax/code`):
  `scan.py` (cutoff scan driver, Hann and Blackman–Harris estimators, bands, occupations), `analyze_scan.py`,
  `equilibria.py` (self-consistent Rayleigh–Jeans lattice end points, quantum Hartree solver with different
  numerics, onset, massive-phase equilibria), `refit_n32.py`, `mc24.py` and `mc24_checks.py` (tree-level 2↔4
  collision-integral Monte Carlo), `recheck_scalars.py`, `efolds_corrected.py`.
- Results in `skeptic/results/*.json` (largest 80 KB). Scratch and restart files in `$SP/review/pt_relax_skeptic/run/`.

## Summary

pt_relax's verdict survives: the classical ~1.75 budget is a prethermal plateau, it decays even in the continuum
classical theory, and the φ-mediated quantum channels have time to relax it before the φ onset. Two of its numbers do
not survive as stated:

1. **The 2↔4 rate bracket is too high by 10²–10⁵.** An explicit tree-level collision integral (own Monte Carlo,
   three analytic checks) gives c_eff = 3×10⁻⁸ (m/T = 0.36) to 4×10⁻⁷ (m/T = 0.16) in Γ = c V⁴T, against pt_relax's
   assumed c = 3×10⁻⁵–3×10⁻³. The e-fold table changes from (0.4–38, 2×10²–2×10⁴, 2×10⁴–2×10⁶, 10⁶–10⁸) to
   (**5×10⁻³, 0.2–3, 18–250, 1.5×10³–2×10⁴**). "Only pure-X 2↔4 is marginal" becomes: pure-X is negligible, the
   g²-only 2↔4 channel is marginal, and the λ_φ channels carry the relaxation with a margin of 10¹–10⁴, not 10²–10⁸.
2. **The N32 late-time numbers are lattice-specific.** At fixed box and seed, the drop of the all-mode total from
   s = 29 to 249 is 0.185 (N32), 0.073 (N40), 0.060 (N48), 0.056 (N64), extrapolating to ≈ 0.058 in the continuum,
   so N32's early decay is 3.2× the continuum value. But the decay does not vanish with cutoff: at s ≈ 500 the
   totals are 1.523 (N48) and 1.489 (N64) against N32's 1.202, and the late slope on N64 (−1.4×10⁻³ per unit s over
   s 449–499) equals N32's. The classical decay is real; its rate is not converged and is not physical.

## Verdicts (my number next to pt_relax's)

| # | Claim (pt_relax) | pt_relax | Skeptic | Verdict |
|---|---|---|---|---|
| C1 [C] | Codex windows reproduce | 1.7505527/1.7414468; 1.7514644/1.7442332 | 1.7505527/1.7414468 (N32), 1.7514644/1.7442332 (N48): identical; own estimator code, same engine. Blackman–Harris window differs by ≤ 7×10⁻⁵ | CONFIRMED |
| C1b [C] | Budget decays beyond s = 34: N32 1.751 → 0.832 (s 1004); Γ = (1.54 ± 0.02)×10⁻³; N48 −3.4% by s = 254 | as stated | Decay confirmed on four cutoffs. Drop 29→249: 0.185/0.073/0.060/0.056 (N32/40/48/64), continuum ≈ 0.058. At s 494–504: N48 1.523, N64 1.489 (N32: 1.202). Free-fit asymptote on N32 ranges 0.17–0.64 with start point; with A_∞ fixed at the self-consistent RJ value 0.376, Γ = 1.25×10⁻³ (max residual 0.012). Γ is 1.2–2.2×10⁻³ depending on the assumed asymptote | CONFIRMED (decay); CORRECTED (N32 magnitudes are lattice-specific; Γ error bar ±30%, not ±1.3%) |
| C2 [D] | s_onset = 2.5–4.6×10¹²; Γ_crit = 2.2–4.0×10⁻¹³; autonomy of the conformal Hamiltonian | as stated | Engine force is `g2 s² + lx r²` with no a(s) anywhere (fieldstep.cpp lines 23–32): autonomous. κ = 1.44067×10⁻²⁴, a_onset 4.9988×10¹² (r² 3600), 9.1645×10¹² (12100), 5.2211×10¹² (equilibrium 3927); s_onset 2.5173/4.6150/2.6292×10¹²; Γ_crit 3.97/2.17/3.80×10⁻¹³; T_h = 0.448 GeV, T/m_X = 4.85 at a_onset; 0.332 GeV, 3.60 at a_ref; V₀/ρ = 0.122/0.404. All reproduce. Caveat: a = a_m + Cs assumes constant g*; the change in s_onset between g* ≈ 107 and ≈ 86 is O(10%), immaterial at 10¹¹ | CONFIRMED |
| C3 [E] | 2↔4 e-folds: X-only 0.4–38 (marginal); g² 2×10²–2×10⁴; mixed 2×10⁴–2×10⁶; φ self 10⁶–10⁸; XX̄↔φφ 4×10⁵–4×10⁶; 1↔2 10⁶–10⁹ | c = (10–10³)/(4π)⁵ | MC c_eff (Bose, μ = 0, mass on all lines): 3.0×10⁻⁸ (m/T 0.36), 4.1×10⁻⁷ (0.16), 2.7×10⁻⁶ (0.08). Corrected e-folds: X-only **4.9×10⁻³**; g² **0.2–3.0**; mixed **18–250**; φ self **1.5×10³–2.0×10⁴**. XX̄↔φφ by my own route: Γ = 2.0×10⁻⁶ per X̄, 5×10⁶ e-folds (pt_relax 1.6×10⁻⁷–1.4×10⁻⁶: agrees within ×1.4–12) | CORRECTED (2↔4 rows); CONFIRMED (2↔2 row) |
| C4 [A] | Only net charge is protected | stated | Agreed; the cutoff scan shows the decay persists on the finest grid | CONFIRMED |
| C5 [A] | Chemical relaxation is the deciding calculation; the budget very likely relaxes | stated | Agreed, with a thinner margin: the relaxation must go through λ_φ (mixed 2↔4 or XX̄→φφ then φ 2↔4), giving 10¹–10⁴ e-folds; g²-only is marginal and pure-X negligible. In a φ-free sub-model the budget would survive the f ~ 1 phase | CONFIRMED, qualified |
| C6 [E] | End point 0.069–0.110 (pilot); 0.10–0.14 (heavy, R = 0); 1.9–2.9 (R = 40.6–48.5); T = 17.56 | as stated | Own Hartree solver (Gauss–Legendre in k/T, nested roots): T = 17.5561, X̄/net = 0.06914 (φ eq), 19.4026/0.09499 (φ absent): agree to 4 digits. Free BE: 19.403/0.1097 (μ = 0) and 0.0960 with the condensate μ on the X̄ branch; 17.532/0.0809 (3 dof). Heavy release by a different route (massive X at a_ref, free gas, μ_c = 0.861m or m): R = 0: 0.073–0.083; R = 1.14: 0.18–0.19; R = 4.05: 0.39–0.41; R = 40.6: 2.26–2.31 (T/m 8.06); R = 48.5: 2.59–2.65. Self-consistent lattice RJ end points: 0.376/0.333/0.296/0.241 (N32/40/48/64), vs pt_relax's harmonic 0.41–0.61 (N32) and 0.307 (N48) | CONFIRMED (pilot and heavy spread); CORRECTED (N32 RJ end point 0.38, not 0.41–0.61) |
| C7 [C/A] | N32 occupations ≥ 14 (s 200), lattice ends at 10.7 < T; N48 front at occupation 0.1–0.2 by s 200–254 | as stated | N32 min occupation 14.2 (s 200), 10.9 (s 250); ω_max = N/3 = 10.67/13.3/16/21.3; Bose number peak 1.59T = 28, energy peak 2.82T = 49.5. N48 min 0.176 (s 200), 0.102 (s 250), 0.077 (s 500). N64 min 2.7×10⁻³ (s 250, k ≈ 28), 0.022 (s 500); shells with occupation < 1 hold 0.6% (X) + 0.5% (φ) of the energy and 0.0024 of A₋ at s 250 | CONFIRMED |
| C8 | Slow matching, IR width, matrix response, export | not checked | not checked here either | NOT CHECKED |
| C9 | REVIEW_EARLY "conserved through τ 1460" corrected | own correction | Sharpened: −3.2% (N64) to −3.4% (N48) by s = 254; −13% (N48) to −15% (N64) by s = 504 | CONFIRMED |

## 1. Cutoff scan (own driver and estimator; Codex engine as reproduction tool)

Same L = 12, seed 1, width 1 state on N = 32, 40, 48, 64 (Codex's `prepare.initialize` seeds each k ≤ 8 mode by its
own `(seed, L, i, j, k)` SeedSequence, so the infrared realization is identical on every N ≥ 32). h ≤ 0.004, Δs = 0.04,
10-s Hann windows (M = 251), bands: zero mode, 0 < |k| ≤ 8, |k| > 8. Charge sum rule A₊ − A₋ = 1 holds to 4.6×10⁻⁹,
energy drift ≤ 2.0×10⁻⁷, Blackman–Harris minus Hann ≤ 7.4×10⁻⁵. N32 and N48 agree with pt_relax's `long_run.json`
window by window to 0.0 (same engine, independently written estimator). Compute: about 2.4 h on 2 threads.

All-mode total A₋ per net charge (coherent / band / UV in parentheses):

| window | N32 | N40 | N48 | N64 |
|---|---:|---:|---:|---:|
| 24–34 | 1.7510 (0.411/1.330/0.010) | 1.7505 | 1.7501 (0.405/1.339/0.006) | 1.7497 (0.404/1.341/0.005) |
| 94–104 | 1.7173 | 1.7258 | 1.7265 | 1.7260 |
| 194–204 | 1.6507 | 1.6935 | 1.7028 | 1.7044 |
| 244–254 | 1.5655 (0.511/0.982/0.073) | 1.6772 | 1.6900 (0.598/1.013/0.080) | 1.6941 (0.649/0.972/0.073) |
| 294–304 | | | 1.6645 | 1.6730 |
| 394–404 | | | 1.5858 | 1.6007 |
| 494–504 | 1.202 (pt_relax) | | **1.5229** (0.577/0.850/0.095) | **1.4889** (0.568/0.834/0.087) |

Slopes dA/ds (per unit s): 29–249: −7.4/−3.3/−2.6/−2.4×10⁻⁴ (N32/40/48/64); 249–499: −7.1×10⁻⁴ (N48), −7.9×10⁻⁴ (N64);
449–499: −6.4×10⁻⁴ (N48), −13.6×10⁻⁴ (N64). Drop 29→249 by Aitken on N40/48/64: 0.058; drop 29→499: 0.227 (N48),
0.261 (N64), i.e. the late decay is **larger** on the finer grid. φ energy fraction at s ≈ 250: 3.3% (N32), 1.3%
(N40), 1.0% (N48), 0.8% (N64); at s ≈ 500: 8.7% (N48), 8.9% (N64). Quadratic energy above k = 8 at s ≈ 500: 19.5%
(N48), 21.2% (N64).

Reading: the early (s ≲ 250) loss is a converged classical effect of ≈ −2.3×10⁻⁴ per unit s; N32 overstates it 3×
because its cascade front hits the lattice edge (ω_max = 10.7). After s ≈ 250 the loss accelerates on every grid, and
the N64 rate exceeds N48's, which is the classical UV catastrophe opening more phase space. In the quantum theory the
modes at k ≳ T ≈ 18 are Bose-suppressed, so the late classical rate is bounded neither above nor below by these runs,
as pt_relax says. The coherent part does not decay monotonically within s ≤ 504 (N64: 0.40 → 0.65 → 0.57); pt_relax's
"0.411 → 0.152" is an N32, s = 1004 value only.

## 2. 2↔4 collision integral (own Monte Carlo)

Model: one real scalar, L_int = −(λ/4!)φ⁴, mass m on every line (thermal), Bose or Boltzmann statistics at μ = 0,
T = 1. Tree 2 → 4 amplitude: ten 3+3 partitions, M = λ² Σ 1/(q² − m²); with equal masses |q² − m²| ≥ m², so no
on-shell pole. R₂₄ = (1/2!)(1/4!)∫dΠ₁dΠ₂ f₁f₂ ∫dΦ₄ |M|² Π(1 + f_i) by RAMBO (massless, rescaled to mass m with the
standard weight, boosted to the plasma frame for the Bose factors). Linearizing in a number chemical potential,
dn/dt = −4R_eq μ/T and χ = ∫f(1 + f)/T, so Γ_chem = 4R_eq/(Tχ) and c_eff = Γ_chem/(λ⁴T). This holds exactly with
Bose factors (δlnR₂₄ − δlnR₄₂ = −2μ/T).

Checks (`mc24_checks.json`): constant |M|² = 1, massless Boltzmann 2 → 4 against the closed form
(1/48)(2π)⁻⁸(π/2)³(1/12)(12/π⁴): ratio 0.99943 ± 0.001; massive 3-body RAMBO weight against the Dalitz integral:
0.04574/0.04577, 0.30553/0.30559, 0.72994/0.73016 at √s = 3.5/5/10 (m = 1); 2 → 2 massless Boltzmann rate against the
analytic value: 6.3816/6.3824×10⁻⁶; convergence with 4×10⁶ events and a new seed: 3.015×10⁻⁸ (m/T 0.36), 4.136×10⁻⁷
(0.16) versus 3.023×10⁻⁸ and 4.111×10⁻⁷.

| m/T | c_eff Bose | c_eff Boltzmann |
|---:|---:|---:|
| 1.0 | 3.07×10⁻¹⁰ | 1.74×10⁻¹⁰ |
| 0.5 | 8.47×10⁻⁹ | 2.38×10⁻⁹ |
| 0.36 (M_φ/T) | 3.02×10⁻⁸ | 5.56×10⁻⁹ |
| 0.16 (M_X/T) | 4.11×10⁻⁷ | 2.36×10⁻⁸ |
| 0.08 | 2.68×10⁻⁶ | 5.16×10⁻⁸ |

pt_relax's c = 3.2×10⁻⁵–3.2×10⁻³ is 80–10⁵ times larger. The rate is strongly mass dependent (≈ m⁻³ between
0.36 and 0.08, a collinear enhancement cut off by the thermal mass), so the Hartree masses matter. Caveat: the
complex-X and mixed-vertex channels have different diagram combinatorics, O(1–10), which I have not computed; it does
not close a 10²–10⁵ gap. The corrected e-folds (`efolds_corrected.json`, same V⁴T formula and s_onset = 2.63×10¹²):

| channel | pt_relax e-folds | corrected (c_eff 3×10⁻⁸–4×10⁻⁷) |
|---|---:|---:|
| 2↔4 X only (V = 4λ_X) | 0.38–38 | **4.9×10⁻³** |
| 2↔4 two g² vertices | 2.4×10²–2.4×10⁴ | **0.22–3.0** |
| 2↔4 mixed g²·λ_φ | 1.9×10⁴–1.9×10⁶ | **18–250** |
| 2↔4 φ self | 1.5×10⁶–1.5×10⁸ | **1.5×10³–2.0×10⁴** |
| 2↔2 XX̄ ↔ φφ (own route) | 4×10⁵–4×10⁶ | 5×10⁶ |

The chain XX̄ → φφ (5×10⁶ e-folds) followed by φ 2↔4 (≥ 1.5×10³) still equilibrates μ_X̄ = −μ before the onset, so
pt_relax's conclusion stands. What changes is the margin and the identification of the marginal channel.

## 3. Equilibria (own solver)

Quantum Hartree at fixed n_Q = 8112.665, ρ = 105481.25 with a circular condensate: T = 17.5561, μ = 2.0394,
M_X = 2.8355, M_φ = 6.3985, X̄/net = 0.06914, φ share 28.8%, condensate charge 97.6% (φ in equilibrium);
T = 19.4026, X̄/net = 0.09499 (φ absent). These are pt_relax's values to four digits by different numerics, so the
solver agrees; the physics caveats are: (i) the Hartree spectrum has a gap M_X − μ = 0.80 (4.5% of T) at k = 0,
violating Hugenholtz–Pines; a Bogoliubov treatment changes only k ≲ 3, which holds < 1% of the thermal number at
T = 17.6, so the X̄ count moves by ≲ 1%; (ii) the "free BE 2 dof" top of pt_relax's range, 0.110, uses μ = 0 on the
X̄ branch; with the condensate μ it is 0.096, so the bracket is better quoted as 0.069–0.096 (0.110 at μ = 0).

Self-consistent classical (Rayleigh–Jeans) lattice equilibria at the same (n_Q, ρ), FD4 spectrum, X + φ: X̄-action
per net charge 0.376/0.333/0.296/0.241 for N = 32/40/48/64, with T_cl = 1643/843/487/205 and φ share 28–29%. The
N32 value is below pt_relax's harmonic bracket 0.41–0.61 and 26% below the fitted asymptote 0.51, so "the end point
matches the RJ estimate" is weak; "the lattice end point is an artifact that falls with N" is confirmed.

Massive phase after the onset (free X/X̄ at m_conf = 3.62 or 4.88, charge reservoir at μ_c = 0.861m or m, optional
radial scalar at √6 m): pilot X̄/net 0.042–0.096, consistent with pt_relax's 0.056–0.107. Heavy release at a_ref with
R from kin_source: 0.073–0.083 (R = 0), 0.18–0.19 (1.14), 0.39–0.41 (4.05), 2.26–2.31 (40.6, T/m = 8.06), 2.59–2.65
(48.5). The R ≥ 40 cases leave only 1–24% of the charge condensed. pt_relax's massless-phase route gives 1.86–2.92:
same conclusion, 20% apart.

## 4. Attempted counterexamples that failed

- **Estimator artifact?** Two windows (Hann, 4-term Blackman–Harris) agree to 7×10⁻⁵; the charge sum rule holds to
  5×10⁻⁹; the band split moves but the all-mode total decays on all four grids. Not an estimator effect.
- **Pure lattice artifact?** The drop by s = 249 converges (0.073 → 0.060 → 0.056, Aitken 0.058), and the late loss
  grows with N. Not a UV-cutoff artifact in the sense of vanishing in the continuum.
- **Protected number?** All allowed one-condensate-leg processes conserve gas + condensate number in the free
  Hartree picture (X_c + X → XX, X̄ + X_c → X + X̄ is kinematically closed for μ < M_X), as pt_relax argues; number
  change needs 2↔4 or the φ sector. In a φ-free sub-model the f ~ 1 relaxation would be 5×10⁻³ e-folds, so the X̄
  budget would survive the quantum phase; with g² = 0.01 and λ_φ = 0.03 it does not. The conclusion is therefore
  conditional on φ being present and populated, which XX̄ → φφ (5×10⁶ e-folds) ensures.
- **Wrong clock?** Autonomy verified in the engine; Γ_phys/H = Γ_c a/C grows with a, so the conformal count is the
  right one; g* steps change s_onset by O(10%).

## New issues (beyond pt_relax)

1. The 2↔4 bracket c = (10–10³)/(4π)⁵ is unsupported; the computed c_eff is 80–10⁵ smaller (section 2).
2. The late classical decay rate increases with lattice cutoff (N64 > N48 at s ≈ 500): the classical UV catastrophe
   contaminates the s ≳ 250 rates on every grid, so no classical lattice number for Γ should be quoted as a bound in
   either direction.
3. pt_relax's free-fit Γ = (1.54 ± 0.02)×10⁻³ on N32 is overprecise: the asymptote is unconstrained (0.17–0.64 across
   start points) and Γ ranges 1.2–2.2×10⁻³; the self-consistent RJ end point for N32 is 0.376.
4. The Hartree end point has a 0.80 gap at k = 0 (Hugenholtz–Pines violation); immaterial for the count (≲ 1%).
5. The "free BE" 0.110 end point mixes μ = 0 on the X̄ branch with a condensate at μ = 2.01.
6. The coherent zero-mode part does not decay monotonically on N48/N64 through s = 504, so the "coherent 0.41 → 0.15"
   statement is N32-specific.

## Limits

- One seed, one box (L = 12). The N64 run has shells with occupation 2.7×10⁻³ (s 250) to 0.02 (s 500) above k ≈ 18,
  where the classical approximation is wrong; those shells hold ≲ 1% of the energy and 0.0024 of A₋.
- The Monte Carlo is a single real scalar with λφ⁴/4!; complex-X and mixed-vertex combinatorics (O(1–10)) are not
  computed, nor are Bose-enhanced rates in the over-occupied phase (that is what the lattice measures).
- No quantum kinetic evolution was run; nothing here sets n̄, f_sym or the threshold.

## Reproduce

```bash
source $SP/review/pt_relax_skeptic/env.sh      # RF = private git-archive copy of 5f5513b research/reality_framework
cd tracks/pt_relax/skeptic/code
for N in 32 40; do python3 -B scan.py $N 254 3000; done      # restartable: repeat until "SAVED ... s 254.0"
for N in 48 64; do python3 -B scan.py $N 504 540; done       # repeat until s 504 (N64: ~9 min per 60 s of s)
python3 -B analyze_scan.py; python3 -B equilibria.py; python3 -B refit_n32.py
python3 -B mc24.py 1600000; python3 -B mc24_checks.py; python3 -B recheck_scalars.py; python3 -B efolds_corrected.py
```
