# Skeptic of track pt_matching: coherent matching and the φ-condensation transition

27 October 2026. Claude review of Codex 5f5513b, issue #13. Key: `pt_matching_skeptic`. Briefs: BRIEF_TRANSPORT.md
and BRIEF.md. Code and results are in `skeptic/code/` and `skeptic/results/`. Per-run files (about 1,800) are in
`$SP/review/pt_matching_skeptic/runs/` and `.../trackcode/`.

This is one supplied, reduced, classical model and same-session AI work. Nothing here is a validation, an exclusion or
new theory.

## Methods

Two kinds of number appear below. Each one is labelled.

**Independent: my own code and a different method from the track's.**
- **Orbit and ramps:** scipy DOP853 (an adaptive explicit RK, not symplectic) with rtol 10⁻¹²–10⁻¹³.
- **Action:** scipy `quad` in the angle variable.
- **Floquet:** integrates the scalar equation u'' = 4E − 3λ_X u² for u = r², rather than the 2-D orbit. A vector of
  B values is integrated in one solve.
- **Transition:** my own C adaptive DOP853 (`sk_dop.c`; Hairer coefficients taken from scipy). The work integral is an
  extra ODE component.
  - It does not use the track's order-8 symplectic composition.
  - Its window and escape grids differ from the track's.
  - It uses a different phase design: uniform random phases with log-uniform seeds in [10⁻¹², 10⁻²], as well as the
    track's 16-phase × 3-seed grid.
- **Shared input:** only the τ = 20 mean field (`MATCH_L12.npz light_mean_q/p`).

**Reproduction, labelled "track code":** pt_matching's own `hom.c` integrator, imported read-only. It was run on my
random design with its output written to my scratch directory. No Codex code was run.

## Verdicts

| Id | Claim (track) | Track | Skeptic | Verdict |
|---|---|---|---|---|
| 3a [D] | n₋ = m\|B\|²; harmonic E/m = 2J_r + q | harmonic member 1.982382394205 | 1.98238239420 (closed form and spectral); Ω/ω_r = 0.5000000000000 | CONFIRMED |
| 3b [C] | J_r/q, A₋/q of the τ = 20 orbit | 1.7637252793970, 1.7377883531122 | 1.7637252793970, 1.737788353110; per net charge 1.7646827365883 / 1.7387317301624 | CONFIRMED |
| 3c [D] | J_r − A₋ = Σ_{n≥1} nωₙ\|cₙ\|² + Σ_{n≤−2}(\|n\|−1)\|ωₙ\|\|cₙ\|² | 0.025937 = 0.016534 + 0.009403 | 0.0259369 = 0.0165336 + 0.0094033 | CONFIRMED (derivation below) |
| 3d [C] | slow free matching returns J_r | T = 80: mean 1.76372527936, spread 1.99889×10⁻⁵; T = 20: spread 3.5612×10⁻⁴ | 1.76372527937, 1.99890×10⁻⁵; 3.56125×10⁻⁴ | CONFIRMED |
| 3e [C] | mass-only ramp | −4.80, +1.19, +4.79, −1.18 (×10⁻⁷) | −4.796, +1.190, +4.790, −1.181 (×10⁻⁷) | CONFIRMED |
| 3f [C] | sudden range over continuous phase | 0.00535–18.312, mean 8.42 | 0.0053510 (phase 0.77357) – 18.31206 (0.27357); mean 8.42149 = closed form (⟨\|π\|²⟩/2m + m⟨r²⟩/2 − q)/2q | CONFIRMED |
| 3g [C] | T dependence of the ramp error | cosine T^−2.02; C^∞: 3×10⁻³, 1×10⁻⁶, 3×10⁻¹⁰ | cosine slope −2.04 (T 20–160); C^∞: 3.04×10⁻³, 9.65×10⁻⁷, 2.72×10⁻¹⁰ | CONFIRMED |
| Floquet [C] | φ = 0 becomes unstable at 0.732 M_c,peri for both orbits; ν/ω_r 0.05 at 0.87, 0.15 at 1.75 | 0.732 | **0.73376** (τ = 20), **0.73259** (late) by bisection; ν/ω_r 0.0475 and 0.153; band to 3.7926 | CONFIRMED, refined (track's 3000-point grid has resolution 0.0035) |
| 3h [C], τ = 20 orbit | slow-limit transition | n₋/q 0.724 ± 0.026, N_φ/q 0.459 ± 0.017, sum 1.184 ± 0.009 | **0.704 ± 0.017**, 0.475 ± 0.011, 1.178 ± 0.006 (560 runs, ε 10⁻³–3×10⁻⁶). Track code on my random design: 0.712 ± 0.036 | CONFIRMED |
| 3h [C], late orbit | same, late orbit | n₋/q **0.166 ± 0.011**, N_φ/q 0.141 ± 0.006, sum 0.306 | **0.197 ± 0.004**, 0.124 ± 0.002, 0.321 ± 0.002 (736 runs). Track code: 0.189 ± 0.006, 0.131 ± 0.004 (384 runs). Both codes: 0.193 ± 0.004 | **CORRECTED** |
| 3h, ε trend | no ε trend for ε ≤ 3×10⁻³ (down to 3×10⁻⁵) | −0.011 ± 0.037 per decade | +0.006 ± 0.022 per decade (τ = 20, ε to **3×10⁻⁶**); −0.011 ± 0.010 (late) | CONFIRMED, extended one decade |
| 3h, seed | 0.66 / 0.73 / 0.78 for σ = 10⁻² / 10⁻⁶ / 10⁻¹⁰ | (trend read) | slope −0.006 ± 0.006 per decade of seed (416 runs); σ = 10⁻²⁰, 10⁻³⁰, 10⁻⁶⁰ give 0.62 ± 0.12, 0.67 ± 0.11, 0.80 ± 0.10 | no seed dependence; the track's trend is noise |
| 3h, label | "CORRECTED" applied to Codex's "slow matching gives ≈ 1.76" | CORRECTED | Codex stated the result as conditional ("if the actual continuation is adiabatic"; "not calculated alternatives for the binding transition") | label overreaches; see Issue 5 |
| 3i [C] | J_r not conserved; the loss comes from "φ's finite inertia" | slaved mass 4×10⁻⁶; Born–Oppenheimer 5×10⁻⁴ | slaved mass 4.1×10⁻⁶ (ε 10⁻⁴), ∝ ε (own code). **Tilt counterexample** (below): same m_φ, but J_r kept to ≤ 4×10⁻⁴ | non-conservation CONFIRMED; **mechanism REFUTED as stated** |
| 3i, chaos | 10⁻¹² seed change moves n₋ by up to 0.96 at ε ≤ 10⁻³; 2×10⁻¹² at 10⁻² | | the two codes agree run by run to 1.1×10⁻⁶ at ε = 10⁻² (τ = 20 and late), and diverge at ε ≤ 10⁻³ (median \|Δn₋\| 0.29, maximum 1.24) | CONFIRMED |
| 3j [C/E] | waves: actual φ(s) at 2 M_peri, ε = 10⁻³ | k = 6: 1.88; k = 1: 0.0085; smooth reference ≤ 1.1×10⁻⁵ | k = 6: **1.84** (phase range 1.81–1.85); k = 1: 0.0086; smooth reference **5×10⁻⁹** (k = 1) | numbers CONFIRMED; reference CORRECTED (Issue 6) |
| 3k [E] | Γ(φ → X X̄); m_φ/m_X; Γ/H ≈ 2×10¹¹; 1.18 / 0.31 X̄ per charge | | Γ = 1.8757×10⁻⁴ m_X and 2.449 confirmed. Γ/H ≈ 1×10⁸–1×10⁹ (Issue 7). Totals **1.178 / 0.321** | CONFIRMED with caveats |

## Checks and new results

**3c derivation [D].** The θ_r cycle at fixed θ_φ has zero φ-winding, because Ω = Δθ/T_r. So it is homologous to the
radial cycle, and J_r = ⟨Re(ψ̇* ∂ψ/∂θ_r)⟩ = Σ nωₙ|cₙ|².

"Zero only for the harmonic ellipse" is also right. A two-frequency orbit Ae^{iω₁s} + Be^{iω₂s} solves ψ̈ = −f(|ψ|)ψ only
if ω₁² = ω₂². Numerically, the spectral J_r matches the quadrature to 7×10⁻¹³ and the spectral charge to 1.4×10⁻¹².

**New [D]: the φ equation is an exact Lamé n = 4 equation.** For H = |π|²/2 + λ|ψ|⁴/4, u = |ψ|² obeys
u'' = 4E − 3λu² at any charge (checked on the orbit to 9×10⁻¹¹). So u = −(2/λ)℘(s + s₀), and

  φ_k'' = [(2g²/λ_X)℘ + B]φ_k, with B = λ_φM² − k² and 2g²/λ_X = 20 = n(n+1), so n = 4.

The chart in B has exactly four finite instability gaps plus the semi-infinite unstable region:

| Gap | τ = 20 orbit | Late orbit |
|---|---|---|
| 1 | k = 12.18–12.20 at M = 0 (ν/ω_r 8×10⁻⁴) | k = 9.27 |
| 2 | **k = 6.70–7.32** at M = 0 (ν/ω_r 0.027): the early pilot's neutral band | k = 5.00–5.33 |
| 3 | M/M_peri 0.734–3.793 (ν/ω_r ≤ 0.22) | 0.733–1.609 |
| 4 | 4.194–6.842 (≤ 0.70) | 1.833–2.843 |
| semi-infinite | from 6.853 | from 2.851 |

This explains why both orbits put the first edge near 0.733: the gaps follow Lamé algebra. For example, the gap-2 and
gap-4 lower edges sum to b = 0.500 exactly for both orbits.

Two consequences:
- A φ_k mode at mass M behaves like k = 0 at M_eff² = M² − k²/λ_φ, so the homogeneous mode enters gap 3 first.
- Near the slow-limit escape point (0.74–0.76 M_peri), only k ≲ 0.27 is unstable, i.e. wavelengths ≳ 23 code lengths.
  That is larger than the L = 12 box, so the zero-mode treatment is locally reasonable there.

**New [C]: the final state lies on a line.** For the τ = 20 orbit:
- corr(n₋, N_φ) = −0.995;
- n₋/q is consistent with a uniform distribution on [0, 1.443] (KS p = 0.57; mean/max 0.488, sd/max 0.282, against
  uniform 0.5 and 0.289);
- the mean is therefore about half the maximum, and n₋ ≈ 0.40 J_r.

So the exit value of one combination is sharp and the split between X and φ is random, as in a microcanonical exit.
The late orbit is less uniform (p = 2×10⁻⁷).

The spread of the final energy 2n₋ + (m_φ/m_X)N_φ does **not** shrink like √ε: its sd goes 0.18 → 0.15 from ε = 10⁻³
to 3×10⁻⁶. So the phase space is mixed, not ergodic, and the citation of Ott (1979) [R] is an analogy, not a derivation.

With the corrected late value, both orbits give n₋ ≈ 0.40 J_r and N_φ ≈ 0.26 J_r (τ = 20: 0.399 and 0.269; late:
0.410 and 0.258).

**Tilt counterexample [C] (refutes the stated mechanism in 3i).** I added V = −βλ_φM³φ. At large M this moves the
vacuum to φ = xM with x³ − x = β, and the measures were redefined about that vacuum.

| Tilt β | ε = 10⁻³ | 10⁻⁴ | 10⁻⁵ |
|---|---|---|---|
| 0 | 0.709 ± 0.028 | 0.694 ± 0.029 | 0.707 ± 0.056 |
| 0.01 | 0.81 ± 0.07 | – | – |
| 0.1 | 0.84 ± 0.08 | 0.57 ± 0.07 | – |
| **0.3** | **1.7636** (sd < 10⁻⁴) | **1.7636** | **1.7633** |
| 1 | 1.711 ± 0.006 | 0.843 ± 0.007 | 1.457 ± 0.004 |

J_r = 1.7637. What this shows:
- **β = 0.3:** φ has the same mass (comparable to ω_r), yet J_r is conserved to ≤ 4×10⁻⁴ with N_φ ≤ 5×10⁻⁴.
- **So finite φ inertia is not sufficient.** The loss requires φ to grow out of the orbit-pumped (gap-3) instability of
  the symmetric state.
- **β = 1:** the loss is deterministic and non-monotonic in ε, which looks like resonance capture.
- **GRUT's case:** the gauge-point tilt variant (Δ = 10⁻²² GeV⁴, v = 0.924 GeV) corresponds to
  β = 3Δ/(4λ_φv⁴) = 3.4×10⁻²¹ [D, from Q]. That is the unbiased regime, so the track's model is the relevant one.

**Wave Hartree shift [D].** A frozen wave contribution g²W to the φ mass enters only as M² → M² − g²W/λ_φ.
- Codex's s 24–34 export gives W = ⟨δχ_X²⟩ = 5468 against a coherent ⟨|χ̄|²⟩ = 4626, so waves hold 54% of the X
  variance [Q].
- For the late orbit, this moves the onset to ≈ 2.5 M_c,peri and multiplies the effective ε at onset by ≈ 11. That is
  still inside the observed plateau, so the homogeneous result is unchanged.
- The X̄ in waves (1.27 per charge) through the transition is modelled by no one.

**Other checks.**
- **Radiation-era growth M ∝ 1 + κs** (ε = 2×10⁻³, MEND 5, 64 runs): 0.680 ± 0.050, against the track's
  0.618 ± 0.065.
- **Measure convergence:** n₋ at 40 versus 5 M_c,apo differs by ≤ 0.004. Using g|φ| instead of gM changes it by
  < 10⁻⁵ at 40 M_c,apo.
- **Integration quality:** work identity ≤ 1.4×10⁻⁸ and charge ≤ 1.2×10⁻⁸ in every run.

## Issues

1. **Late-orbit numbers [C].**
   - The track's n₋/q = 0.166 ± 0.011 and N_φ/q = 0.141 come from 96 grid runs.
   - Large random-design ensembles in both codes give 0.193–0.197 ± 0.004 and 0.124–0.131. The track's ε = 10⁻³ grid
     was a low sample: both codes give 0.16 on that grid, because the late runs at ε = 10⁻³ are only weakly chaotic
     and partly correlated between the codes.
   - Corrected values: n₋ = **0.41 J_r**, not 0.35 J_r, and the φ → X X̄ total is 0.32, not 0.31.
2. **Mechanism [C].** "φ's finite inertia" is refuted by the tilt counterexample above. The decisive ingredient is
   growth from the Floquet-unstable symmetric state, which is gap 3 of an exact Lamé n = 4 problem.
3. **The seed trend is noise.** There is no seed dependence from 10⁻² to 10⁻⁶⁰. There is also no ε dependence down to
   3×10⁻⁶ (one decade below the track).
4. **The ε → 0 statement is supported, but not by Ott's theory.** The final energy spread levels off (mixed phase
   space). The broad n₋ distribution is a random split along a sharp line.
5. **Verdict label for 3h.** Codex wrote that slow matching gives J_r *if* the continuation is adiabatic. It called its
   ramps "not calculated alternatives for the cosmological binding transition".
   - The homogeneous X + φ result shows that this premise fails in the first dynamical model of the GRUT transition.
     That is a NEW finding.
   - It is not a corrected value of a Codex number. Suggested verdict: Codex's conditional statement is CONFIRMED;
     "premise fails in the homogeneous model (0.70 + 0.47 per charge for the τ = 20 orbit, 0.20 + 0.12 for the late
     orbit)" is NEW.
6. **Wave test (3j).**
   - **Smooth reference:** the 10⁻⁵ is an initial-state mismatch. The modes start in the vacuum of √(k² + 20),
     omitting g²M_start², so δω/ω ≈ 0.6%. The adiabatic value is ≤ 5×10⁻⁹.
   - **Small k is not specific to φ.** A crude isotropic coherent-orbit modulation, 2λ_X r(s)² with φ frozen, gives
     |β|² = 0.026 at k = 1. Only the k ≈ 6 band is specific to φ.
   - **3j and 3k are one channel.** The overflow at 4–8 M_peri is the parametric decay φ → X X̄, i.e. the channel
     counted in 3k. They must not be added.
7. **Scales and λ_φ [E].**
   - The track gives ε ~ 10⁻¹⁵ and Γ/H ~ 2×10¹¹ (from H/m_X ~ 10⁻¹⁵) without a derivation. With REVIEW_EARLY's
     10¹¹–10¹² conformal units to binding, ε ≈ 1/(sω_r) ≈ 10⁻¹²–10⁻¹³ and Γ/H ≈ 10⁸–10⁹. That is still ≫ 1, so the
     conclusion is unchanged, and the extrapolation now spans 6–7 decades from 3×10⁻⁶.
   - φ → X X̄ is open only for λ_φ > 2g² = 0.02. It is open at the current spot point (λ_φ = 0.03, gauge-point packet
     [Q]). It is closed at the Candidate-1 benchmark λ_φ = 0.01 (m_φ/m_X = 1.414).
8. **"Never above J_r" holds only for ε ≤ 3×10⁻³.** At ε = 10⁻², both codes reach 1.806.

## Limits

- I did not rerun the Born–Oppenheimer control.
- The wave spot check covers k = 1 and 6 at ε = 10⁻³ only.
- All limits of pt_matching apply: a homogeneous zero mode, classical dynamics, no gauge/Y fields, and a pure-quartic
  late stand-in.

## Reproduce

```bash
cd research/reality_framework/grut_formation_review_20261025/tracks/pt_matching/skeptic/code
python3 -B skA_orbit_ramps.py; python3 -B skB_floquet.py
for e in 1e-2 1e-3 3e-4 1e-4; do python3 -B skC_scan.py tau20 $e grid16; done
for e in 1e-2 1e-3 1e-4; do python3 -B skC_scan.py late $e grid16; done
python3 -B skC_scan.py tau20 1e-3 rand128; python3 -B skC_scan.py tau20 1e-4 rand128; python3 -B skC_scan.py tau20 3e-5 rand64
MEND=10 python3 -B skC_scan.py tau20 1e-5 rand64; MEND=10 python3 -B skC_scan.py tau20 3e-6 rand32
python3 -B skC_scan.py late 1e-3 rand384; python3 -B skC_scan.py late 3e-4 rand128; python3 -B skC_scan.py late 1e-4 rand128
for s in seed20 seed30 seed60; do python3 -B skC_scan.py tau20 1e-3 $s; done
for b in 0.01 0.1 0.3 1.0; do python3 -B skC_scan.py tau20 1e-3 rand32 0 $b; done
for b in 0.1 0.3 1.0; do python3 -B skC_scan.py tau20 1e-4 rand32 0 $b; done
for b in 0.3 1.0; do MEND=10 python3 -B skC_scan.py tau20 1e-5 rand16 0 $b; done
MEND=5 python3 -B skC_scan.py tau20 2e-3 rand64 1
python3 -B skD_waves.py; python3 -B skF_controls.py; python3 -B skE_final.py
# track-code reproduction on the random design (pt_matching's hom.c, read-only import): python3 -B skG_trackcode_crossrun.py late 1e-3 384; python3 -B skG_trackcode_crossrun.py tau20 1e-3 128
```
