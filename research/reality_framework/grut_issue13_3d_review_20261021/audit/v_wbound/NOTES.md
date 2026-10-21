# Verifier `v_wbound`: adversarial check of audit track `wbound`

21 October 2026. Verifies `audit/wbound/` (NOTES.md, results/, code), which audits item (4) of Codex `a69bd86`
against Claude `c7b8480`. All code here is new; none of the `wbound` scripts was reused. Codex's packet was read with
`git show` / `git archive` only. Units: v = 1, λ_φ = g_X² = 0.01, λ_X = 0.001 (benchmark couplings). Labels: [T]
theorem, [C] computed, [E] estimate, [A] argued. Nothing here excludes, validates or identifies Candidate 1.

## Method (independent of wbound)

| Script | What it does differently | Output |
|---|---|---|
| `algebra_check.py` | Keeps X complex (X and X̄ treated as independent fields); builds T_μν in Hilbert form; computes Φ_c directly, without the Euler-degree argument. Checks the FRW identity symbolically for arbitrary a(t), with a gradient term up to a spatial divergence. Integrates the boosted kink in the lab frame | results/algebra_check.json |
| `homog_check.py` | Period integrals by event-terminated DOP853 integration (wbound uses θ-substitution quadrature). Targeted φ–X scan at fixed energy (X seeds 0.01–1 v, with or without charge). Homogeneous radiation-FRW runs closing the exact ρa⁴ budget | results/homog_check.json |
| `radial_lump.py` | Spherically symmetric 3D kinetic lumps, R₀ = 30–400. A Cartesian lattice cannot reach these radii. Tests the dilatation theorem and how long the excess persists | results/radial_lump.json |
| `collide_1d.py` | Separate planar φ + X lattice, kink–antikink at γ = 2.29, 7.09 and 22.4. Measures the *local, time-averaged* w of the debris after the walls have left the window | results/collide_1d.json |
| `codex_eos_check.py` | Rebuilds ρ and p from Codex's K, G, U columns and checks them against Codex's ρ and p columns. Integrated continuity check | results/codex_eos_check.json |

## Results

1. **Trace identity and virial integrands [T, sympy].** 3p − ρ − (Σḟ² − Σ(∇f)² − 4U) = 0 and T^μ_μ = 3p − ρ.
   Φ_{+v} = −λv(φ−v)²(φ+v) − 2g²vφ|X|² − D²(m_D0² + h_D²vφ) and Φ₀ = λv²(φ²−v²) − m_D0²D² both have residual 0.
   With D = 0 and φ < 0, the root in |X|² is −λ(φ−v)²(φ+v)/(2g²φ), which reproduces the wbound condition.
2. **FRW identity [T].** 3p − ρ = Φ_c + a⁻³d(a³g_c)/dt has residual 0, homogeneous and with a gradient up to a
   divergence. d(ρa⁴)/dt = −Ha⁴(3p − ρ) has residual 0. Numerically, homogeneous runs with a ∝ t^½ from a = 1 to
   4 (A₀ = 1.2, 2.2, 3.0; H₀ = 0.02, 0.005) close Δ(ρa⁴) = −∫Ha⁴Φ₀dt − [ȧa³g₀] + ∫äa³g₀dt to **1.9×10⁻¹¹–3.0×10⁻⁹**.
   The wbound 3D lattice closure (9.3–9.6×10⁻⁴) is therefore discretization error, not a missing term.
   (2/√3)(v/M_P) = 4.7415×10⁻¹⁶ at v = 1 TeV.
3. **Homogeneous orbits [C].** The independent ODE method gives ε* = 0.587678 λv⁴ (amplitude 1.59160 v) and
   w̄_max = 0.441506 at ε = 2.12703 λv⁴ (amplitude 1.97911 v). Further values: w̄ = 0.39665 at amplitude 3v
   (ε = 16 λv⁴) and 0.33935 at ε = 2500 λv⁴. In-well w̄ runs from −0.00094 down to −0.397 as the barrier is
   approached. The identity holds to ≤ 5×10⁻¹¹. All values agree with zero_d.json to ≥ 5 digits. VERIFY's 0.0225 v⁴
   is the amplitude-2v energy (2.25 λv⁴), as wbound says.
4. **X at fixed energy [C].** 37 φ–X runs at ε = 2.127, 5 and 16 λv⁴. Every run has w̄ below the φ-only value
   (0.4423 / 0.4278 / 0.3966; finite-T): the X runs reach 0.286–0.328, 0.329–0.382 and 0.358–0.388. Even an X seed of
   0.01 v lowers w̄ at 2.127 λv⁴ from 0.442 to 0.309–0.317, because parametric resonance drains the over-barrier
   oscillation. **But wbound's own sampling has X *raising* the maximum below the barrier** (zero_d.json: e = 0.5,
   −0.017 vs −0.066; e = 0.95, −0.063 vs −0.276). Those values are still far below 1/3. "Adding X lowers the maximal
   w̄" therefore holds above the barrier only, and the "no" answer to "does X make it worse" is [C], not [T].
5. **Dilatation theorem [T + C].** Algebra checked: dD/dt = ∫T^ii, |𝒫| ≤ ρ for U ≥ 0, and D(T) ≤ (R₀+T)E, which
   give ∫(3P − E)dt ≤ ER₀ + |D(0)|. Radial runs:

   | R₀ | profile, centre ε/λv⁴ | max ∫(3P−E)dt / (ER₀) | dD/dt = 3P residual | w̄_T last > 1/3 at | w̄_T at end |
   |---|---|---|---|---|---|
   | 30 | cos², 2.0 | 0.140 | 4.3×10⁻⁴ | t = 5 | 0.094 |
   | 100 | cos², 2.0 (dr 0.2 and 0.1 agree) | 0.046 | 2.3×10⁻⁴ | t = 6 | 0.058 |
   | 400 | cos², 2.0 | 0.012 | 3.9×10⁻⁴ | t = 6 | 0.128 |
   | 100 | flat top, 2.128 | 0.065 | 2.7×10⁻⁴ | t = 30 (0.30 R₀) | 0.078 |
   | 400 | flat top, 2.128 | 0.043 | 3.1×10⁻⁴ | t = 196 (0.49 R₀) | 0.123 |
   | 100 | cos², 0.02 (linear) | 0.054 | 1.7×10⁻⁴ | t = 7 | 0.016 |

   The bound is never violated. These spherically symmetric runs contain no X, so they are the most favourable case
   for a lasting excess, and even here it lasts ≤ 0.5 R₀. Spherical focusing drives |φ| up to 17–22 v at the centre
   for R₀ = 400. This is a symmetry effect and does not affect the bound.
6. **Collision debris [C], new.** max |φ| over t > t_col + 20 is 1.26, 1.62 and 2.11 v at γ = 2.29, 7.09 and
   22.4. Over t ≥ t_col − 5, which includes the overlap itself, my run gives at least 2.05, 2.31 and 2.85 v (sampled only
   every 10 time units). wbound's
   1.74–2.83 v are therefore the **collision-instant overshoot**: sims_1d.py takes the maximum from t ≥ t_col − 5.
   The maximum relaxes to ≲ 1.4 v within about 90 time units. **After the walls leave**, the windows |z| < 30 and
   |z| < 100 (t from t_col + 150 to about t_col + 440) have a local time-averaged w (direct, from 3p − ρ) of
   −0.045 to +0.016. Box post-collision w̄ is 0.065, 0.270 and 0.312, against wbound's 0.052, 0.257 and 0.312. That
   rise toward 1/3 is the receding walls (NG u² − 2/3 = 0.143, 0.313 and 0.331), not the debris. **No locally
   time-averaged w̄ > 1/3 state produced by walls was found at γ ≤ 22.**
7. **Codex runs [C].** ρ and p rebuilt from K, G and U match Codex's columns to 4×10⁻¹⁶. In all 12 runs
   w_max ≤ 0.140, ρa⁴ is monotone at every sample, and ln(R_end/R_0) agrees with ∫(1 − 3w)d ln a to ≤ 0.013. ρa⁴
   grows ×5.80–6.83 in 11 runs, but **only ×2.37 in cold_s101**, so wbound's "×5.8–6.8 in all 12" is wrong for the
   cold control.
8. **Finite-T caveat (not an error in wbound).** With origin 0 the boundary term scales as ∝ 1/(A m T) for
   small-amplitude lumps, because a kinetic lump carries zero-mode momentum. In sims_3d_small R6, ⟨Φ₀⟩/ρ = +0.146
   while w̄ = 0.294 < 1/3, and the boundary term closes the gap. So the "exact criterion" is asymptotic. Sign tests on
   finite runs need the boundary term (wbound includes it) or the ±v origin. The cosmological bound already uses the
   first-order worst case |g₀| ≤ (2ρ⟨f²⟩)^½, so it is unaffected.

## Verdicts on the wbound claims

| # | wbound claim | Grade | Reason |
|---|---|---|---|
| V1 | trace identity | CONFIRMED | sympy, complex X, residual 0 |
| V2 | Φ_{+v} ≤ 0 for φ ≥ −v (φ only) | CONFIRMED | residual 0; factored form reproduced |
| V3 | X: pointwise worse; exact criterion X-free | CORRECTED (minor) | Pointwise condition and Φ₀ confirmed. "X lowers max w̄" holds only above the barrier; wbound's own e = 0.5 and 0.95 data show the reverse. "No" is [C], not [T] |
| V4 | Codex: instantaneous w ≤ 1/3 not a theorem | CONFIRMED | Trivial counterexample (w = 1); for radial lumps the running mean is w̄_T(t = 1) = 0.97–0.98 and falls below 1/3 by t = 5–7 (cos²) |
| V5 | Codex: needs a state family; "supplied exactly" | CORRECTED | Identity, budget and 4.74×10⁻¹⁶ bound confirmed (closure ≤ 3×10⁻⁹). But it is an exact *reduction* to one inequality, which remains [A] for the post-parity network. "This is the justified state family" overclaims. The weight is a⁴ d ln a, not just Hubble weighting |
| V6 | Claude w̄ ≤ 1/3 → Hubble-integral iff condition | CONFIRMED | Correct as stated. §6's "it can fail only through coherent over-barrier regions" is [A]: the virial identity is global and cannot be localized (flux terms), so walls mixed with large one-sign waves are not covered by T1/T2 |
| V7 | "only coherent both-vacua oscillations" is too strong | CORRECTED | wbound's counterexamples are themselves coherent over-barrier oscillations (kinetic lumps) or instantaneous collision overshoot; T1 makes "both signs" necessary. Debris after collision has local w̄ ≈ 0. Valid corrections: time-average qualifier; 0.4415 is the homogeneous φ-only maximum, not a proven bound; the incoherent hot-plasma route [E] |
| V8 | "30 orders above parity, so cannot occur" | CONFIRMED (REFUTED as an argument) | The mean-density inference is invalid: wall cores sit at (2.7–4.0)×10³⁰ ρ_h(par). Production of *time-averaged* w̄ > 1/3 states by walls is not shown (item 6) |
| V9 | ε* = 0.58768 λv⁴ = 5.9×10⁹ GeV⁴ | CONFIRMED | independent ODE: 0.587678 |
| V10 | Φ_{+v} ≤ 0 for φ > 0 ⇒ time-averaged w̄ ≤ 1/3 | CONFIRMED | each term ≤ 0 for φ ≥ 0 |
| V11 | NG/field walls | CONFIRMED | lab-frame sympy: ∫(3p−ρ)dz = −3σ/γ, w = u² − 2/3 exactly |
| V12 | von Laue, quanta | CONFIRMED | standard argument (∂_jT^{ij} = 0 ⇒ ∫T^{ij} = 0); not recomputed numerically |
| V13 | ΔN_eff,eq ≥ 7.447 | CONFIRMED | 1/0.227107 + 3.044 = 7.4472 |
| V14 | Codex DERIVATION §3 omits −2G | CONFIRMED | a69bd86 DERIVATION.md line 81 vs lines 30–31; expanding.cpp uses p = (K − G/3 − U)a⁻⁴ |
| V15 | Codex 12 runs | CORRECTED (minor) | w_max ≤ 0.140 and monotone ρa⁴ confirmed; cold_s101 grows ×2.37, not ×5.8–6.8; NOTES says "CONSISTENT", verdicts.json says "CONFIRMED" |
| V16 | "outside F … not produced" | CORRECTED | Local-density theorem (T3) confirmed. "They *are* produced (wall cores, collisions)" is shown only for instantaneous overshoot. Static wall cores have \|φ\| ≤ v and Φ₀ < 0, so they are not outside F |

Also checked: the Cauchy–Schwarz cap w̄ ≤ 1/3 + (2/3)(λv⁴/ρ̄)^½ is CONFIRMED (sympy). The §5 arithmetic is CONFIRMED
(6.5–7.3 ρ_w; R_max ≳ 6t ≈ 3H⁻¹), and its [A] status is agreed. The catalogue row "hot restored plasma … needs mean
density ~ v⁴, irrelevant after parity" uses the mean-density inference that V8 refutes. It should say instead that
local energy densities ≳ v⁴ would be needed (not excluded by the mean; in principle reachable in γ ~ 10³⁰
collisions, where the classical treatment is itself unreliable).

## Limits

- No full 3D rerun of wbound's lattices. The radial runs impose spherical symmetry and contain no X; the collisions are
  planar with γ ≤ 22.4, as in wbound; the homogeneous FRW runs have no gradients. NOT CHECKED: the 1D slab table, the
  3D network and charged-condensate rows, the 3D FRW lattice runs (only their reported closures were read), the
  quantum trace anomaly, metric perturbations and GW back-reaction.
- The φ–X scan is targeted but finite (37 states); it is not a proof that 0.4415 bounds w̄ for all multi-field states.
- Codex diagnostics are sampled every Δτ = 1.2 (41 points), so instantaneous maxima between samples are not seen.

## Reproduce

`python3 algebra_check.py`; `python3 homog_check.py`; `python3 radial_lump.py`;
`python3 collide_1d.py 0.9 0.99 0.999` (each a few minutes); `python3 -I codex_eos_check.py <extract>/raw`, where `<extract>`
comes from `git show a69bd86:research/reality_framework/grut_expanding_charge_20261021/raw_diagnostics.tar.gz | tar -xz`.
All single-threaded.
