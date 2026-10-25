# Track sp_phaseprep: can inflation prepare the supplied source-phase mismatch of Codex 373d37e?

25 October 2026 · Claude (independent review track, issue #13) · packet `373d37e`
(`origin/codex/grut-inflation-scale-scan-20261025`, `research/reality_framework/grut_inflation_scan_20261025/`),
read only through `git archive 373d37e … | tar -x` into the scratchpad. No Codex code was run. Nothing was committed.
Own code is in `code/`, outputs in `results/`. Reused code of ours: gp_matching `model.py`/`point_full.py` (spectrum),
sp_point `spmodel.py` (potential, radial minimum, output parser), and sp_point `sprel.c`, extended here to
`code/pprel.c` with tagged preparation terms.

**One supplied tree-level working model. Nothing here is evidence that X, Y, Σ or the inflationary inputs exist. The
"phase-preparation interaction" (PPI) below is a hypothetical test term. It is not in the stated action and is not
proposed as physics. This is not a validation, an exclusion or a GRUT prediction. CMB comparisons are retrospective.**
All agents are same-session AI, so this is not outside replication.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled (not re-read this session).

## 1. Inputs [Q]

- Couplings unchanged: M_Y = 1e12, h = 1e11, μ = 5e11 GeV, κ = 1e-3, e = 1e-7, λ_Y = λ_Σ = 1, λ_X = 1e-3. Light
  target g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV, Q = 2e5; Y_Q = 4.8294e-9.
- Point: H_i = 1.1211167304108865e13 GeV, c_H = 1. Aligned inflationary minimum; supplied state Σ = −|Σ| e^{iπ/2}
  with X > 0, Y < 0, zero velocities. Instant removal of the X Hubble mass, radiation, g* = 106.75.
- Tolerance to be met: Codex's Gaussian phase noise σ = 0.03 rad (ε_net = −0.249 ± 0.003). From sp_point: the
  coherent band is δ ≈ 1.46–1.64 rad (half-width 0.09 rad), and σ ≲ 0.08 rad is needed. Radii must stay near the
  aligned values (re-checked in §5).
- Isocurvature limit: S_max = 9.1e-6 for uncorrelated CDM isocurvature [Q: Planck 2018 X, β_iso < 0.038, as used in
  our affleck_dine and gp_matching tracks].

## 2. Method (own code)

- `code/spectrum.py` [C]: the full canonical 6×6 Hessian (from gp_matching `model.py`), the gauge-invariant phases
  ψ₁ = 3θ_X − θ_Y and ψ₂ = 2θ_Y − θ_Σ, and f_eff = 1/|dψ₂/dφ_c| along the light angular eigenvector.
  - Exact Bunch–Davies coarse-grained variance of each mode, σ²(x_c) = ∫₀^{x_c} dx/x P(x), with
    P = (H²/8π) x³ e^{−π Im ν}|H_ν^{(1)}(x)|² (mpmath, complex ν for m > 3H/2).
  - x_c = k/(a_end H) is the cut for modes that are still outside the horizon when the charge forms: x_c = 1/40 for
    the charge-creation epoch a ≈ 25–40 (sp_point), and 1/500 for the end of sensitivity growth.
  - Control at 4.7e13 GeV.
- `code/fokker_planck.py` [C, E]:
  - (A) the overdamped Starobinsky Fokker–Planck equation for compact ψ₂ at fixed radii: equilibrium-preserving finite
    volume, its spectrum, and evolution from uniform, aligned and required starts;
  - (B) the equilibrium P ∝ exp(−8π²V/3H⁴) in the complex Σ plane, with X and Y fixed or re-minimised;
  - (C) a second-order Langevin equation for Σ (noise H/2π per √e-fold) started at the supplied state, 20,000
    realisations.
- `code/pprel.c` + `code/ppi.py` [C]: the hypothetical DRT-type PPI (§5), its inflationary minimum and spectrum (gauge
  mode projected out exactly), and releases. We tested instant and finite switch-off, and ensembles drawn with each
  mode's exact coarse-grained variance (N = 128 each; zero velocities). Integrator: Yoshida-8, c = 0.025 (sp_point
  convention). It reproduces Codex's point at ε = −0.24066123 and N_X = −238.64867, and the Gauss residual stays below
  1e-12 in every run.
- `code/transfer.py` [C]: the long-wavelength transfer T = d ln|Y_X|/dψ_L, local (±1e-7 rad) and smoothed over short-scale
  noise with common random numbers.
- `code/isocurvature.py` [C, E]: observable-scale coherence, S/S_max, and order-of-magnitude new sources.

## 3. Phase spectrum and decay constant at H_i [C]

| Quantity | H_i = 1.1211e13 GeV (ours) | Codex (ENSEMBLE_SUMMARY) | 4.7e13 GeV control (gp_matching) |
|---|---|---|---|
| \|X\|, \|Y\|, \|Σ\| / H | 23.826064, 1.893762, 0.430868 | same | 23.82417, 1.891386, 0.266976 |
| radial masses / H | 1.0546, 1.4096, 4.6516 | m² 1.11218, 1.98688, 21.6376 | 0.6539, 1.4094, 4.6484 |
| angular masses / H | 0 (gauge), **0.60958**, 2.7613 | m² 0.371585, 7.62478 | 0, 0.37809, 2.75063 (0.378, CONFIRMED) |
| light mode (dθ_X, dθ_Y, dθ_Σ per canonical H) | −0.0032, −0.0008, 1.6316 | — | −0.0019, −0.0037, 2.6429 |
| dψ₂, dψ₁ per canonical H | −1.633, −0.0087 | — | −2.650, −0.0020 |
| **f_eff** | **0.6123 H = 6.86e12 GeV** (1.005 √2\|Σ\|) | — | 0.3773 H (0.377, CONFIRMED) |
| **kick per e-fold H/(2π f_eff)** | **0.260 rad** | — | 0.422 rad |

- The light mode is the Σ phase, and therefore the mu source phase [C].
- Σ behaves as a quartic field with a linear tilt μY² [D]. This gives m_t² = 2λ_Σ\|Σ\|² = 0.3712 (Hessian: 0.3716)
  and m_r² = 3 m_t² = 1.114 (Hessian: 1.112). **The phase mass is tied to the radius.** That is the obstruction met
  in §5.

## 4. Stochastic distribution and probability of the required phase [C, E]

*Fokker–Planck, compact ψ₂ at fixed radii* (V = B(1 − cos ψ₂), B = 0.13783 H⁴):

- Equilibrium P ∝ exp(κ cos ψ₂) with κ = 8π²B/3 = 3.628. This gives ⟨cos ψ₂⟩ = 0.847 and rms ψ₂ = 0.58 rad. The
  density at π/2 is 0.027 of the density at the minimum.
- The slowest relaxation rate is 0.1004 per e-fold, i.e. **10.0 e-folds** (3H²/m² = 8.1).
- A state placed exactly at π/2 (σ = 0.03) leaks out:
  - P(\|ψ₂ ∓ π/2\| < 0.09) is 1.00 at the start, 0.36 after 0.5 e-fold, 0.25 after 1, 0.045 after 8, and 0.0073 at
    equilibrium (30–60 e-folds);
  - P(< 0.03) is 0.69 at the start, 0.12 after 0.5 e-fold, and 0.0024 at equilibrium.

*Complex Σ plane* (radius and phase are not hierarchical: 0.61 H against 1.05 H):

- ⟨cos θ⟩ = 0.81 and rms θ = 0.67 rad. Re-minimising X and Y changes these only in the third digit.
- At θ = π/2 the noise-supported radius peaks at **0.31 H**, not at the aligned 0.43 H, nor at the conditional minimum
  h/√2 = 0.006 H. Only 8.6% of that slice lies within 5% of 0.4309 H.
- Exact BD coarse-grained ψ₂ spread: 0.50 rad (x_c = 1), **0.29 rad (x_c = 1/40)**, 0.26 (1/100), 0.21 (1/500).
  The Starobinsky–Yokoyama Gaussian gives 0.52.

*Probability that one Hubble patch is in the band (equilibrium, either CP sign):*

| half-width w | phase band, 1D | phase band, Σ plane | supplied box (band AND \|Σ\| within 5%) |
|---|---|---|---|
| 0.03 rad | 0.24% | 0.37% | **0.031%** |
| 0.09 rad | 0.73% | 1.1% | **0.096%** |
| 0.18 rad | 1.5% | 2.2% | 0.20% |

*The observable universe* [D, E]:

- The light mode relaxes in 10 e-folds, which is much shorter than the N_* ≈ 50–60 e-folds of observable scales. A
  value set when the observable scale left the horizon therefore survives only with the memory factor
  e^{−(3/2−ν)N_*} = 1.6e-3 (N_* = 50) to 4.4e-4 (N_* = 60).
- The long-wavelength rms on the observable scale is σ_L ≈ 7e-4 rad (N_* = 50) to 2e-4 rad (60), centred on the
  aligned minimum.
- Holding the observable-scale mean within 0.09 rad of π/2 is a Gaussian excursion of exponent ~10⁶–10⁷ (see
  isocurvature.json). **The probability is nil.**
- Correlations fall only as r^{−0.26}. The observable universe contains about 10⁶⁶ independent release patches
  (40 H_i⁻¹ at a ≈ 40), and each draws from the equilibrium above. The volume-averaged charge is the CP-symmetric
  ensemble mean, which is zero.
- A single lucky patch (0.03–0.1%) is a region of a few Hubble lengths, not an observable universe.

*Ensemble at the natural spread.* In the "pure A-term" row of §5 (the physical masses are unchanged, and the draws
use x_c = 1/40, ψ₂ rms 0.275) the net fraction is **ε_net = −0.072 ± 0.010**, with 0.58 pairs per net charge. At
x_c = 1/500 it is −0.116 ± 0.011.

*Freshness.* A preparation switched off before the end of inflation decays.
- From noise in Σ alone, P(band 0.03) is 0.75 after 0.01 e-fold and 0.40 after 0.05; P(box 0.09) is 0.54 after 0.05
  and 0.33 after 0.1.
- These noise figures count modes up to the Hubble scale. The modes that act on whole patches (x < 1/40) left the
  horizon at least 3.7 e-folds before the end. For a later switch-off the binding limit is therefore the classical
  drift of the full X, Y, Σ system: sp_point found ε = −0.20 after 0.02 e-fold and −0.17 after 0.05.
- **The preparation must stay on to within about 0.02 e-fold of the end of inflation.**

## 5. A minimal phase-preparation interaction (hypothetical; DRT type) [D, C]

*Requirements from §3–4 and the tolerance scan* (single releases from the supplied state with one radius rescaled;
`ppi.json → radius_tolerance`):

- Phase mass: σ_ψ(x_c = 1/40) ≤ 0.03 rad needs m ≥ 1.15 H, and ≤ 0.08 rad needs m ≥ 0.95 H. At x_c = 1/500 the
  thresholds are 1.0 H and 0.8 H. sp_point's "≳ 1.5 H" used the Starobinsky–Yokoyama Gaussian, which overestimates
  the spread at m ~ H.
- Displacement: ψ₂ = π/2 ± 0.09 rad, with a phase that is the same everywhere. It must therefore be a fixed CP phase
  of the coupling, not a field.
- Radii: |Σ| within about ±1–3% (−3%: −0.218, +3%: −0.023, ±5%: fail); |Y| within about ±0.2–0.5% (±0.5%: −0.199 and
  −0.165); |X| within ±0.1% (all pass).
- Switch-off: below.

*Obstruction [D].* A Hubble-induced A-term alone, V ⊃ μ_P Σ*Y² + c.c., rotates the minimum. But it keeps
m_t = √(2λ_Σ)|Σ|. A mass of 1.15 H would need |Σ| ≥ 0.81 H, about 1.9 times the tolerated radius. A Hubble mass
c_S H²|Σ|² is therefore needed. Its larger source then shifts |Y|: c_S = 0.5, 1 and 2 untuned give δY/Y = 0.25%,
0.51% and 1.0%, and ε = −0.208, −0.073 and −0.014 [C]. That calls for a third term, c_Y H²|Y|².

*Minimal PPI:* V_PPI = c_S H²|Σ|² + c_Y H²|Y|² + [μ_P Σ*Y² + c.c.]. It is gauge invariant, Z2_φ even, and has
dimension ≤ 4. With |μ+μ_P| = [2λ_Σ s(s² − h²/2) + c_S s]/y², arg(μ+μ_P) = π/2 and c_Y = 2(|μ+μ_P| − μ)s, its
inflationary minimum is **exactly** Codex's supplied state (radii to 1e-9, a 6D multistart agrees to 1e-15) [D, C]:

| c_S | c_Y | \|μ+μ_P\|/μ | \|μ_P\| (GeV), arg | light mass / H | f_eff / H | ψ₂ rms (x_c = 1/40) | ε_net (N = 128) | pairs / net charge |
|---|---|---|---|---|---|---|---|---|
| 0 (pure A-term rotation) | 0 | 1 | 7.07e11, 3π/4 | 0.610 | 0.612 | 0.275 | −0.072 ± 0.010 | 0.58 |
| 0.5 | 0.0518 | 2.35 | 1.28e12, 1.974 | 0.927 | 0.622 | 0.081 | −0.192 ± 0.008 | 0.18 |
| **1** | **0.1035** | **3.69** | **1.91e12, 1.835** | **1.153** | 0.635 | **0.026** | **−0.242 ± 0.004** | 0.111 |
| 2 | 0.2071 | 6.39 | 3.23e12, 1.726 | 1.485 (X radial 1.41) | 0.669 | 0.002 | −0.239 ± 0.002 | 0.095 |
| 1, x_c = 1/500 | | | | | | 0.006 | −0.241 ± 0.002 | 0.090 |

- Codex's supplied state is exactly the minimum of a pure A-term rotation, μ → iμ (|a_H| = √2 μ/H = 0.063). That
  preparation leaves the phase as light as before, so it fails on its own fluctuations.
- With c_S ≳ 1, c_Y ≈ 0.1 and |μ_P| ≈ 0.17 H, all O(H) as is generic for Hubble-induced terms [R: DRT], the window
  survives the PPI's own residual spread. Every member has the same sign, and 91–100% have |ε| ≥ 0.19.

*Switch-off speed* (single releases; a chaotic map, so these are patterns, not sharp thresholds). The PPI terms are
multiplied by exp(−(a−1)/Δ), with Δ in units of H_i⁻¹:

| Δ | 0.003 | 0.01 | 0.03 | 0.1 | 0.3 | 1 |
|---|---|---|---|---|---|---|
| PPI alone (X Hubble mass removed instantly), ε for c_S = 0 / 0.5 / 1 / 2 | −0.21/−0.25/−0.28/−0.22 | −0.30/−0.30/−0.25/−0.29 | −0.23/−0.24/−0.27/−0.32 | −0.01/+0.11/+0.04/−0.07 | ≤ 0.12 | ≤ 0.07 |
| PPI and the X Hubble mass together | −0.20/−0.20/−0.19/−0.13 | −0.08/−0.18/−0.20/−0.21 | ≤ 0.08 | ≤ 0.15 | ≤ 0.14 | ≤ 0.11 |
| Codex's X Hubble mass alone (no PPI) | −0.242 | +0.091 | −0.049 | +0.033 | — | — |

- Hubble-tracking persistence (g_A = a⁻², g_m = a⁻⁴ in radiation, as for DRT terms that scale with H) gives
  |ε| ≤ 0.11 in all 8 runs.
- **The window needs the PPI to vanish within ≲ 0.03 H_i⁻¹, and the X Hubble mass to vanish within ≲ 0.003–0.01
  H_i⁻¹.** That is far shorter than one X oscillation (about 5.7 H_i⁻¹).
- Standard DRT/Affleck–Dine Hubble-induced terms persist ∝ H through inflaton oscillation [R]. The window therefore
  needs a sudden end of inflation (instant preheating or a waterfall), and this assumption is new and untested.
- This agrees with sp_yledger, where an inflaton matter era cuts ε to 0.025.

*New sources and perturbations the PPI introduces:*

1. **A new physical CP phase**, arg(μ_P μ*) ≈ 1.84 rad, in the inflaton–dark coupling. The stated action has 0
   physical phases at dimension ≤ 4 with Z2_φ (gp_phase [Q]). This phase is the source of the asymmetry, and it is
   supplied by hand.
2. **Three tuned Hubble-induced coefficients plus a tuned end of inflation.**
   - c_Y must hold |Y| to about 0.2–0.3%. For c_S = 1 that is ±0.04–0.06 out of 0.10 in c_Y, since
     dlnY/dc_Y ≈ −0.049.
   - |μ+μ_P| must be set to ±1.5–5% for |Σ| (d ln|Σ|/d ln|μ+μ_P| = 0.65 at c_S = 1), and arg to ±0.09 rad.
   - H_end must lie within about ±1% of 1.1211e13 GeV (sp_point).
   - The switch-off must take ≲ 0.003–0.03 H_i⁻¹.
3. **Y debris is not cured.** With the PPI, the pairs per net charge are 0.09–0.11 (18–22 times the 0.00493
   allowance). Over all 63 PPI releases the minimum is 0.024 (5 times the allowance) [C]. sp_yledger's leftover
   verdict (×350 or more) stands.
4. **Isocurvature:** negligible for any PPI that meets the spread requirement (§6).
5. **Sudden-quench production** when the PPI masses drop: about 1e-3 to 1e-2 H³ quanta, against N_X = 239, i.e.
   ≲ 1e-4 per net charge [E, flat-space sudden approximation].
6. **Inflaton → Σ*YY through the coupling that produces μ_P** (if it arises from g I Σ*Y² with ⟨I⟩ ~ M_P): pairs per
   net charge ≲ 1e-6, which is model dependent [E].
7. **Not computed:** an oscillating CP-odd A-term during any inflaton oscillation, which could itself source charge
   or Σ/Y quanta; thermal masses after reheating; and sub-horizon (x > 1/40) PPI-era modes as seeds of fragmentation
   (lattice needed).

## 6. CMB isocurvature [C, E]

The X charge (as Q-balls) is the calibrated dark matter. Phase fluctuations are uncorrelated with ζ, so
S = T δψ_CMB, with δψ_CMB = C_ν/(2π f_eff) e^{−(3/2−ν)N_*} per ln k and n_iso − 1 = 3 − 2ν.

- *Transfer* (`transfer.json`):
  - local T = d ln|N_X|/dδ = **1005 per rad** (dε/dδ = 257.2, as in sp_point);
  - smoothed over short-scale noise σ_S = 0.03: **T = 2.3 ± 0.3 per rad** (128 common draws at 5 offsets; ε_net =
    −0.269, −0.261, −0.260, −0.245, −0.220 at δ₀ − π/2 = −0.04 … +0.04);
  - smoothed over σ_S = 0.29: T = 1.4 ± 1.5.
  - The local value applies only if the short-scale spread is below the 1e-3–4e-3 rad scale on which the map is
    jagged (sp_point).

| Case | m / H | δψ_CMB (N_* = 50 / 60) | S/S_max, local T | S/S_max, smoothed T (σ_S = 0.03) |
|---|---|---|---|---|
| Baseline light phase, displacement imposed (not an attractor) | 0.610 | 3.7e-4 / 1.0e-4 | 4.1e4 / 1.1e4 | **95 / 26** (σ_S = 0.29: 58 / 16) |
| PPI c_S = 0.5 | 0.927 | 2.3e-8 / 9.4e-10 | 2.6 / 0.10 | 6e-3 / 2e-4 |
| PPI c_S = 1 | 1.153 | 3.7e-13 / 1.7e-15 | 4e-5 / 2e-7 | 1e-7 / 4e-10 |
| PPI c_S = 2 | 1.485 | ~1e-29 | negligible | negligible |

- **If the phase stayed at 0.61 H with the displacement somehow imposed, the uncorrelated CDM isocurvature would
  exceed the Planck level by about 16–95 times** (smoothed T), or by 1e4 times if the local transfer applied. The
  spectrum would be blue, n_iso ≈ 1.26, whereas the quoted bound is for a near-scale-invariant spectrum [R], so this
  is indicative only.
- S < S_max needs m ≳ 0.70–0.78 H (smoothed T) or 0.88–0.95 H (local T) for N_* = 60–50. **Any PPI that meets the
  spread requirement (m ≥ 0.95–1.15 H) is automatically safe** [C].
- The residual isocurvature then comes only from whatever sets the PPI coefficients. If its CP phase or amplitude
  depends on a light field, there is a new isocurvature channel. If the coefficients follow H_end, the effect is
  adiabatic [A].
- Retrospective note [D, R]: H_CMB ≥ H_end gives r ≥ 2.0e-3 (single-field normalisation, A_s = 2.1e-9). That is
  consistent with BICEP/Keck 2021 (r < 0.036).

## 7. Verdicts

| # | Codex 373d37e statement | Ours | Verdict |
|---|---|---|---|
| 1 | Masses at H_i (ENSEMBLE_SUMMARY): angular m² 0.371585, 7.62478; radial 1.11218, 1.98688, 21.6376 | Same to 5 digits (canonical 6×6 Hessian); light 0.6096 H, heavy 2.761 H [C] | CONFIRMED |
| 2 | CORRECTIONS 2: the light angular mode overlaps the mu source phase; the masses support a linear isocurvature estimate | Light mode = Σ phase (dψ₂ = −1.633, dψ₁ = −0.009 per canonical H), f_eff = 0.612 H, kick 0.26 rad per e-fold; estimate in §6 [C] | CONFIRMED |
| 3 | The supplied mismatch is a favourable assumption, not in the baseline action; no phase selector has been added or tested | The inflationary state of the action gives the box with 0.03–0.1% per Hubble patch and nil per observable universe; at the natural spread ε_net = −0.072 ± 0.010 [C, E] | CONFIRMED |
| 4 | The phase-noise distributions are supplied sensitivity tests, not a computed stochastic state | Computed: FP equilibrium rms 0.58 rad and relaxation 10 e-folds; coarse-grained spread 0.29 rad at x_c = 1/40, ten times Codex's 0.03 [C] | CONFIRMED (now computed) |
| 5 | "A small-noise ensemble preserves it" (σ = 0.03) | True for σ ≲ 0.03–0.08, which needs a phase mass ≥ 0.95–1.15 H against 0.61 H in the action. A DRT-type PPI with c_S ≈ 1 reaches −0.242 ± 0.004 [C] | CONFIRMED conditionally; scope CORRECTED |
| 6 | Instant removal of the X Hubble mass (PROTOCOL assumption) | Single releases: exponential removal over 0.003 H_i⁻¹ keeps −0.242; 0.01, 0.03 and 0.1 give +0.091, −0.049 and +0.033 [C] | NOT CHECKED by Codex; the window needs removal ≲ 0.003–0.01 H_i⁻¹ |
| 7 | A spatial spectrum and cosmic mean need an inflationary fluctuation spectrum and finite-k evolution | Spectrum supplied here (exact BD per mode). Finite-k release not done [C] | CONFIRMED; still open |
| 8 | (sibling sp_point) "transverse Σ mode heavy, ≳ 1.5 H" | Exact coarse-grained BD variance: m ≥ 1.15 H (0.03 rad) or 0.95 H (0.08 rad) at x_c = 1/40; SY overestimates near m ~ H [C] | CORRECTED (self) |
| 9 | (sibling sp_yledger) natural spread 0.52 rad, net ε ≈ 0.05 | 0.52 rad is the x_c = 1 SY value; superhorizon at charge formation it is 0.29 rad, net −0.072 ± 0.010. Conclusion unchanged [C] | CONFIRMED (refined) |

## 8. Limits

- Homogeneous k = 0 releases. Ensembles draw superhorizon modes (x < x_c) as homogeneous patch offsets with zero
  velocities. Sub-horizon modes (x > 1/40) seed inhomogeneity and fragmentation and are not evolved (lattice needed).
- The stochastic treatment is overdamped or light-field at m = 0.6–1 H [E]. Corrections are O(m²/3H²). The exact BD
  mode variances are free-field results at fixed masses.
- The PPI is a hypothetical, tagged test term with instant or exponential switch-off in a radiation background. No
  inflaton dynamics, reheating plasma or thermal masses are modelled. The switch-off results are single chaotic runs.
- Tolerances come from single releases in a jagged map (sp_point: autocorrelation lost at 4e-3 rad).
- The new-source estimates (§5, items 5–6) are order-of-magnitude and model dependent.

## 9. Credit (read levels)

Dine, Randall & Thomas, PRL 75 (1995) 398 and Nucl. Phys. B 458 (1996) 291: Hubble-induced masses and A-terms, and a
phase mismatch from an A-term phase [R]. Affleck & Dine, Nucl. Phys. B 249 (1985) 361 [R]. Starobinsky 1986
(stochastic inflation) and Starobinsky & Yokoyama, PRD 50 (1994) 6357 (equilibrium distribution) [R]. Bunch & Davies
1978 (massive dS modes) [R]. Grain & Vennin 2017 and Tokuda & Tanaka 2018 (stochastic dynamics with momentum) [R].
Linde 1985, Seckel & Turner 1985, Kasuya, Kawasaki & Takahashi 2008, and Enqvist & McDonald 1998–99 (spectator and
Affleck–Dine isocurvature) [R]. Planck 2018 X (Akrami et al., A&A 641 (2020) A10), β_iso [Q via our earlier tracks].
BICEP/Keck 2021, r < 0.036 [R]. Kofman, Linde & Starobinsky 1994/97 and Felder, Kofman & Linde 1999 (preheating and
instant preheating) [R]. Kusenko–Shaposhnikov, Enqvist–McDonald and Kasuya–Kawasaki (Q-ball formation) [R]. None of
these papers was re-read in this session.

## 10. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/sp_phaseprep/code
export PYTHONDONTWRITEBYTECODE=1   # imports gp_matching and sp_point modules; keeps their dirs clean
python3 spectrum.py ../results/spectrum.json            # 25 s
python3 fokker_planck.py ../results/fokker_planck.json  # 40 s
python3 ppi.py ../results/ppi.json                      # ~17 min, 4 cores (binary built in the scratchpad)
python3 transfer.py ../results/transfer.json            # ~13 min, 4 cores
python3 isocurvature.py ../results/isocurvature.json    # 5 s
```
