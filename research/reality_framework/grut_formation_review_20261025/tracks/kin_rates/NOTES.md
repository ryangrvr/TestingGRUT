# Track kin_rates: collision kernels of the light sector for the kinetic case-L calculation

26 October 2026 · Claude, issue #13 (BRIEF_KINETIC.md) · same-session AI work, not outside replication. Resumed run:
the code of the interrupted attempt (model, amplitudes, phase space, 2→2 quadrature, thermal masses, MC engine) was
re-run and checked; one defect found and fixed on resumption (section 4.2); everything else below was produced after
the resumption. **Third pass (finalisation, after a second usage-limit cut): no recomputation; the deliverables were re-read
against results/, the documentation inconsistencies listed in section 9 were corrected, and the independent spot checks of
`code/finalise_checks.py` (results/finalise_checks.json) were added.**

Codex packets read only through `git show`/`git archive` into the scratchpad: `grut_kinetic_contract_20261026` at
9877518 (KINETIC_SPEC.md, MATHEMATICS.md, THERMODYNAMIC_CONTROL.json, CLAUDE_HANDOFF.md, VERIFICATION_RESULTS.json;
code read, **not run**) and `grut_heavy_vector_budget_20261025` at 225028b (INPUTS.json). Every number here is from
this track's own code unless labelled [Q]. Own earlier tracks reused, with the reuse stated: the Berends–Giele
recursion method of `tracks/gp_leftovers/code/amp42.py` (re-implemented in `code/amp.py` with a batch dimension,
complex masses and per-line propagator switches; its two threshold amplitudes are reproduced to 1e-15 as a check);
the threshold constants A_a, A_b of gp_leftovers and SV_SS of ls_threshold/lsrelic.py are quoted [Q own] only for
the comparison of section 4. Sibling conventions: `tracks/kin_source/NOTES.md` §1 (comoving numbers per H_i³,
f = 2π²(dN/dln q)/q³, one bosonic degree per species) and `tracks/kin_thermo` (Bose state functions). Nothing was
committed. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued.

**Scope.** Kernels, in-medium masses and their checks for one supplied tree-level model. No Boltzmann solution,
no f_sym, no threshold, no viability, no observational statement; nothing here is evidence that the fields exist.

## 1. Model, conventions, units

Fields: complex X = (x₁ + i x₂)/√2 (charge +1), real radial φ = v + s; the D field of the gauge-point packet is
absent from the spot point (parent INPUTS.json lists g_X, λ_X, λ_φ, v only). Potential [Q INPUTS.json convention]

    U = λ_φ/4 (φ² − v²)² + g² φ² |X|² + λ_X |X|⁴,   g = 0.1, λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV.

Tree masses and vertices (exact derivatives of U at the vacuum, sympy, `code/model.py`) [D]: m_X = g v = 0.0924090 GeV,
m_s = √(2λ_φ) v = √6 m_X = 0.2263548 GeV; λ₃ ≡ ∂³U/∂s∂x₁∂x₁ = 2g²v = **0.0184818 GeV** (the s X X̄ trilinear);
λ_s3 = 6λ_φ v = 0.166336 GeV; λ_s4 = 6λ_φ = 0.18; λ_ssXX = 2g² = 0.02; X⁴: ∂⁴U/∂x₁⁴ = 6λ_X, ∂⁴U/∂x₁²∂x₂² = 2λ_X
(X X X̄ X̄ contact amplitude −4λ_X). Vacuum width Γ(s → X X̄) = λ₃² β/(16π m_s) = **1.73328e-5 GeV**, β = (1 − 4m_X²/m_s²)^{1/2}
= 0.8165 (Codex THERMODYNAMIC_CONTROL 1.7333e-5: CONFIRMED). Amplitude convention: all momenta incoming, metric
(+,−,−,−), M = −∏V ∏ 1/(p² − m² + i m Γ) over tree diagrams (so the 2→2 contact amplitude is −V₄); an outgoing
particle is an incoming antiparticle with momentum −k; polarisations X_in = (1, −i, 0)/√2, X̄_in = (1, i, 0)/√2.

Units GeV. Occupation per bosonic degree, n = ∫d³p f/(2π)³ (X, X̄, s one degree each). Rate density γ = reactions
per volume per time (GeV⁴) **with the symmetry factors included**: S_in = ∏ n_i! over identical incoming species,
S_out likewise; a cross section σ (GeV⁻²) includes 1/S_out; ⟨σv⟩ ≡ γ S_in/(n₁n₂) (GeV⁻²), so that the old
rate-density forms A_a n₊²n₋²/4, A_b n₊³n₋/6, ½ σv n_s² are γ with S_in = 4, 6, 2. The species collision term
carries the multiplicity of that species in the channel. Phase space dΦ_n = ∏ d³p_i/((2π)³2E_i) (2π)⁴δ⁴, so
Φ₂ = |p*|/(4π√s). Thermal averages are Maxwell–Boltzmann at μ = 0 unless marked BE; for MB the chemical potentials
enter exactly as e^{Σμ_in/T}; Bose factors (f and 1+f) are evaluated by the lab-frame MC (section 5) and given as
BE/MB ratios at T/m_X = 0.3, 1, 3, 7.344, 12. Energies are lab energies; "T/m" means T/m_X.

Derivations used throughout [D]: for the same |M|² (any n ↔ 2 pair)

    γ_{n→2}(T) = (1/(S_in S_out)) (T/8π³) ∫ds √s K₁(√s/T) G̃(s),      G̃(s) = ∫dΦ_n(s) ∫dΦ₂(s) |M|²  (dimensionless),
    γ_{2→n}(T) = (1/S_in) (T/32π⁴) ∫ds λ(s, m₁², m₂²) σ(s) K₁(√s/T)/√s,   σ = (1/S_out)(1/2λ^{1/2}) ∫dΦ_n |M|²,

from ∏dΠ_i e^{−E_i/T} = ∫d⁴P/(2π)⁴ e^{−P⁰/T} dΦ_n(P) and ∫d⁴P e^{−P⁰/T} F(P²) = ∫ds 2π√s T K₁(√s/T) F(s). The two
are identical by construction, which is the table-level detailed-balance check of section 5.

## 2. In-medium masses and the phase (`code/thermal.py`, `run_thermal.py`, `results/thermal_masses.json`)

One-loop finite-temperature potential (Dolan–Jackiw) for the background φ with the X doublet and the radial
fluctuation in the loop [D]:

    V(φ,T) = λ_φ/4 (φ²−v²)² + (T⁴/2π²)[J_B(m_s²(φ)/T²) + 2 J_B(m_X²(φ)/T²)],  m_s²(φ) = 3λ_φφ² − λ_φv², m_X²(φ) = g²φ²,

J_B by direct quadrature (series for the tachyonic direction, |y| < 0.03 only). The zero-T Coleman–Weinberg piece is
absorbed in the vacuum parameters (its field dependence is of relative size g⁴/16π² ≈ 6e-6, dropped [A]). Thermal
masses are the curvatures at the thermal minimum v_T, with the tadpole J_B′ and the J_B″ pieces:
m_s²(T) = ∂²V/∂φ²|_{v_T}; m_X²(T) = g²v_T² + (T²/2π²)[8λ_X J_B′(y_X) + 2g² J_B′(y_s)] (high-T limit
T²(λ_X/3 + g²/12)); in-medium couplings λ₃(T) = 2g²v_T, λ_s3(T) = 6λ_φv_T. A Parwani-type daisy resummation (the
T² self-energies added to the loop masses) is run as the bracket [E]. μ = 0 for the loop particles (gas chemical
potentials would add to the tadpoles: limit).

| T/m_X | v_T/v | m_X(T)/m_X | m_s(T)/m_s | m_s/(2m_X) | Γ_s,rest(T)/Γ_vac | Bose factor 1+2f(m_s/2) |
|---|---|---|---|---|---|---|
| 0.3 | 1.0000 | 1.0000 | 1.0000 | 1.2247 | 1.034 | 1.034 |
| 1 | 0.9997 | 0.9999 | 0.9994 | 1.2242 | 1.832 | 1.833 |
| 3 | 0.9926 | 0.9955 | 0.9894 | 1.2172 | 4.93 | 5.02 |
| 6 | 0.9585 | 0.9749 | 0.9492 | 1.1924 | 9.46 | 10.36 |
| 7.344 | 0.9329 | 0.9595 | 0.9203 | 1.1747 | 11.22 | 13.06 |
| 8.118 | 0.9149 | 0.9487 | 0.9002 | 1.1621 | 12.10 | 14.75 |
| 10 | 0.8594 | 0.9162 | 0.8391 | 1.1217 | 13.45 | 19.48 |
| 12 | 0.7780 | 0.8700 | 0.7493 | 1.0548 | 11.64 | 26.17 |
| 13 | 0.7253 | 0.8413 | 0.6896 | 1.0039 | 3.59 | 30.8 |
| 15 | 0.5632 | 0.7596 | 0.5666 | 0.914 | closed | – |
| 18 | 0.1282 | 0.6261 | 0.1267 | 0.248 | closed | – |
| 20 | 0 (symmetric) | 0.6818 (thermal) | 0.329 | – | closed | – |

Results [C]: **T_c = 18.14 m_X = 1.676 GeV** (the broken minimum disappears; high-T estimate 18.09 m_X; daisy
18.23 m_X); no degenerate-minima temperature is found (the transition is second-order-like at this order). **The
broken phase holds at every hidden temperature up to 12 m_X and up to T_c**; above T_c there is no v, no s state as
a radial quantum of the broken phase, and no object. **s → X X̄ closes at T = 13.06 m_X (daisy 13.57)**: m_s(T) falls
faster than 2m_X(T) (m_s/2m_X = 1.055 at 12 m_X); the inverse decay closes with it and only s s ↔ X X̄, the 2→2 and
the 3↔2 channels connect s to the pairs between 13.1 m_X and T_c. The in-medium rest-frame width is Γ_vac(m_s(T),
m_X(T), λ₃(T)) × [1 + 2f_B(m_s/2)] (μ = 0; this is the net Weldon damping factor, see section 3.2 for the gross/inverse split): the Bose factor reaches 26 at 12 m_X while the phase-space factor
β(T) and λ₃(T)² reduce it, net 11.6× at 12 m_X (table; the energy-dependent lab width Γ(E,T) = (m_s/E)Γ_rest ⟨1+f₁+f₂⟩
with the daughter energies is tabulated for E/m_s = 1–100 and is the callable `thermal.width_s_medium(E, T, μ₊, μ₋)`).

**Regime of validity of vacuum masses [C]:** within 1 % for T ≤ 2.92 m_X (m_s) and 4.11 m_X (m_X); within 5 % for
T ≤ 5.95 / 8.03 m_X; within 10 % for T ≤ 8.13 / 10.77 m_X. At Codex's matched-energy states T = 7.344 and 8.118 m_X
the vacuum masses are 4–10 % high; at 12 m_X they are 13 % (m_X) and 25 % (m_s) high, and the vacuum width is 12×
too small. The daisy bracket moves these by ≤ 2 % below 12 m_X and by 0.5 m_X in T_close. All kernels of sections 3–4
are tabulated at vacuum masses; the interface carries masses(T) and the validity statement for the solver, and the
m_s(T)/2m_X(T) switch of the decay channel. Limits: one loop; the J_B″ (IR) term is negative in the broken phase
(−4.4 % of m_s² at 12 m_X) and the daisy bracket does not resolve the infrared sector near T_c; our ls_threshold
quoted a wider T_c range (18–42 m_X) from a different loop content [Q own].

## 3. Kernels

### 3.1 2 → 2 (`code/xsec22.py`, `run_22.py`, `results/kernels_22.json`; amplitudes checked against the recursion to ≤ 4e-15)

    X X̄ → X X̄ : M = −4λ_X − λ₃²[1/(s−m_s²) + 1/(t−m_s²)]          s-pole at m_s² (physical, s ≥ 4m_X²)
    X X  → X X  : M = −4λ_X − λ₃²[1/(t−m_s²) + 1/(u−m_s²)]          no pole
    X s  → X s  : M = −2g² − λ₃²[1/(s−m_X²) + 1/(u−m_X²)] − λ₃λ_s3/(t−m_s²)   u-channel X pole inside the angular
                                                                     range for 3.606 < √s/m_X < 5.0 (s → X X̄, X̄ X → s)
    s s  → s s  : M = −λ_s4 − λ_s3²[1/(s−m_s²) + 1/(t−m_s²) + 1/(u−m_s²)]     no pole (s ≥ 4m_s²)
    s s  → X X̄ : M = −2g² − λ_s3λ₃/(s−m_s²) − λ₃²[1/(t−m_X²) + 1/(u−m_X²)]  no pole; X X̄ → s s by crossing

Physical poles are treated by the narrow-width (RIS) subtraction (section 4): with M = A + B D, D = 1/(x + iε),
|M|² = A² + 2AB Re D + B²|D|², and |D|² → Re D² = (x²−ε²)/(x²+ε²)² removes the on-shell piece (π/ε)δ(x)B², which is
the sequential process already in the network (inverse decay × decay for X X̄ → X X̄; decay × inverse decay for
X s → X s). The angular and s integrals are done with explicit finite-part/principal-value quadratures (`fp_quad`,
`pv_quad`: subtraction of h(0) + h′(0)x and analytic remainders), so the physical width Γ_s/m_s = 7.7e-5 is used
exactly for the s pole (the X pole is regulated with Γ_X,reg = 1e-7 m_X: finite-part limit). Checks: fp_quad on
(1+x+x²)/x² to 3e-6 at ε = 1e-6; the on-shell piece of the XX̄ → XX̄ thermal rate equals the MB decay-rate density
n_s Γ_s⟨m/E⟩ to 1e-6 at every T (RIS normalisation); a brute-force quadrature of the full Breit–Wigner at 10–1000×
the physical width (`check_pole.py`, `results/check_pole.json`) reproduces full − RIS = subtracted to 5e-5 (T = 0.3 m,
×1), 1.5e-3 (T = m), 2.5 % (T = 3 m; narrow-width corrections O(Γ/m)).

σ(s) is tabulated from threshold to √s = 10⁶ m_X (110 points; 24 extra points in the X s window). Samples (GeV⁻²):

| √s/m_X | X X̄ (non-res.) | X X | X s (sub) | s s | s s → X X̄ | X X̄ → s s |
|---|---|---|---|---|---|---|
| 2.09 | 3.23e-6 | 2.12e-5 | – | – | – | – |
| 4.0 | 2.22e-7 | 8.09e-7 | 1.08e-4 | – | – | – |
| 5.0 | 2.77e-7 | 1.67e-7 | 3.4e-3 | 1.97e-2 | 2.94e-4 | 6.87e-6 |
| 10 | 2.03e-7 | 4.86e-8 | 3.51e-6 | 1.11e-4 | 1.09e-5 | 4.30e-6 |
| 100 | 3.60e-9 | 1.78e-9 | 8.26e-8 | 3.51e-6 | 9.00e-8 | 4.49e-8 |
| 10³ | 4.17e-11 | 2.09e-11 | 1.10e-9 | 3.56e-8 | 8.79e-10 | 4.40e-10 |
| 10⁴ | 3.07e-13 | 1.53e-13 | 9.17e-12 | 3.20e-10 | 7.91e-12 | 3.95e-12 |
| 10⁵ | 3.38e-15 | 1.69e-15 | 1.13e-13 | 4.23e-12 | 1.05e-13 | 5.23e-14 |

**High-energy behaviour [D, C]:** every 2→2 cross section goes to the contact limit σ → |M_c|²/(16π s S_out), M_c =
4λ_X (X X̄, X X), 2g² (X s, s s ↔ X X̄), λ_s4 (s s); the tables reach it to 1e-9 at √s = 10⁶ m_X with log-log slope
−2.000 in √s (the t-channel pieces fall as 1/s²). For a vector daughter with p = 7.7e4 m_X on a bath X at rest
(√s = 393 m_X) or s (√s = 614 m_X): σ(X X) = 1.3e-10, σ(X X̄) = 2.7e-10, σ(X s) = 2.9e-9 GeV⁻².

MB thermal averages ⟨σv⟩(T) (GeV⁻²) on T/m_X = 0.03–30 (51 points; the transport weight 1 − cosθ also tabulated):

| T/m_X | X X̄ (non-res.) | X X | X s (non-res.) | s s | s s → X X̄ | X X̄ → s s | BE/MB: XX̄, XX, Xs, ss, ss→XX̄ (results/mc_checks.json) |
|---|---|---|---|---|---|---|---|
| 0.3 | 9.16e-7 | 7.43e-6 | 7.39e-5 | 5.65e-3 | 7.80e-5 | 2.11e-8 | 1.05, 1.05, 1.03, 1.00, 1.00 |
| 1 | 3.21e-7 | 1.33e-6 | 2.11e-5 | 1.74e-3 | 3.45e-5 | 4.01e-6 | 1.61, 2.03, 1.47, 1.16, 1.17 |
| 3 | 1.68e-7 | 9.87e-8 | 4.29e-6 | 2.18e-4 | 8.56e-6 | 3.38e-6 | 2.50, 5.81, 2.79, 2.44, 2.24 |
| 7.344 | 5.57e-8 | 2.26e-8 | 1.26e-6 | 4.40e-5 | 1.87e-6 | 8.94e-7 | 4.20, 6.35, 4.37, 4.22, 4.54 |
| 12 | 2.53e-8 | 1.08e-8 | 5.84e-7 | 2.05e-5 | 7.48e-7 | 3.67e-7 | 5.28, 7.19, 5.42, 4.99, 5.64 |

The X X̄ elastic kernel is dominated by the resonant chain: γ_sub/γ_RIS = +1.9 at T = 0.03 m (sub-threshold
s-exchange, positive), −1.6e-4 at 0.3 m, −1.4e-4 at 1 m, +2e-6 at 3 m, +2.2e-3 at 12 m: the subtracted remainder is
sign-indefinite and ≤ 2e-3 of the chain for T ≥ 0.3 m. The X s subtracted kernel is likewise sign-indefinite in the
window 3.6–5 m_X (σ_sub = −1.7e-5 GeV⁻² at 4.5 m_X against σ_nores = 4.2e-5); the recommended transport kernels are
the non-resonant columns with the chain carried by the explicit decay/inverse decay (the table also holds `sub`,
`res`, `int` and the full-width components σ₀, σ₁, σ₂ of X X̄ so that any prescription can be rebuilt).

### 3.2 s ↔ X X̄ (`thermal.width_s_medium`, `thermal_mc.decay_bose`, `inverse_decay_bose`)

Three Bose factors must be kept apart [D]. The collision term of an s mode is f_s(1+f₁)(1+f₂) − f₁f₂(1+f_s) =
f_s(1 + f₁ + f₂) − f₁f₂ (identity of Bose functions at E = E₁ + E₂). (i) Gross forward rate (the Boltzmann gain of
X and X̄, loss of s): factor ⟨(1+f(E₁))(1+f(E₂))⟩ over the decay angle, lab daughter energies E₁, E₂
(`thermal_mc.decay_bose`). (ii) Inverse decay X X̄ → s: f₁f₂(1+f_s) (`thermal_mc.inverse_decay_bose`). (iii) Net
(Weldon) damping rate of the s mode, the coefficient of f_s: Γ_net(E, T, μ₊, μ₋) = (m_s/E) Γ_vac(m_s(T), m_X(T), λ₃(T))
⟨1 + f(E₁) + f(E₂)⟩ (`thermal.width_s_medium`; the Γ(E,T) tabulated in section 2 and `masses_T.width_vs_energy_oneloop`).
A solver uses (i) with (ii), or (iii) with the spontaneous inverse piece ∝ f₁f₂ alone; (iii) together with the full (ii)
double counts. Rate densities MB: γ_dec = n_s Γ_s⟨m_s/E⟩ = m_s²TΓ_s K₁(m_s/T)/2π² = 1.63e-13, 3.27e-10, 1.04e-8,
8.26e-8, 2.33e-7 GeV⁴ at T/m = 0.3, 1, 3, 7.344, 12 (vacuum masses; ⟨m/E⟩ = K₁/K₂ = 0.84, 0.60, 0.32, 0.15, 0.098;
recomputed in closed form to 1e-16, `finalise_checks.py`). Bose, μ = 0: gross forward γ_dec^BE/γ_dec^MB = **1.024, 1.54,
3.78, 9.38, 15.26** (MC; the inverse decay by exact 2-D quadrature agrees to 0.02–0.6 %, 3.2 % at 12 m where the
quadrature tail is the limit); net Γ_net: 1.024, 1.49, 2.98, 5.25, 6.87 (own 2-D quadrature, `finalise_checks.py`;
the difference is the f₁f₂ term, dominated by slow s quanta whose daughters have E ≲ T). Inverse decay X X̄ → s follows
from detailed balance (with Bose factors and μ_s = μ₊ + μ₋, checked in section 5). With in-medium masses the channel
closes at 13.06 m_X (section 2); the interface flags `decay_open(T)` and states the three conventions in
`decay_s_XXb.bose_factors`. (The interrupted attempt's text quoted 1.055–5.86 for the gross factor, a stale run;
corrected in the finalisation pass, section 9.)

### 3.3 n ↔ 2 (`code/kernels_n2.py`, `phasespace.py`, `run_kernels.py`, `run_supplement.py`; `results/kernels_n2.json`, `xsec_2n.json`, `supplement_n2.json`)

Method: Lorentz-invariant n-body phase space by sequential two-body splitting with multichannel importance sampling
(a Breit–Wigner map of the invariant mass of every neutral X X̄ pair at m_s²); checks: massless Φ₃, Φ₄ against the
closed forms to 1.4σ (0.14 %, 0.33 %), massive Φ₃ against a Dalitz quadrature (2e-4), map independence, the
on-shell and momentum-conservation residuals 2e-16. Full tree amplitudes by the recursion (every diagram, all
interference), RIS subtraction of every line that can be on shell (section 4.2), MC with N = 60 000 (6-point) /
100 000 (5-point) per energy, 16 energies √s = √s_th(1+δ), δ = 2e-3 … (400 m_X), two regulators (Γ_reg,s, Γ_reg,X) =
(0.02 m_s, 0.05 m_X) [primary, every point] and (0.04 m_s, 0.10 m_X) [every other point, half the samples]; the 2→n
cross sections likewise to 10⁴ m_X. The engine reproduces the 2→2 quadrature on the control channels (X X → X X
to 0.2 %; X X̄ → X X̄ and X s → X s at enlarged widths) and the threshold constants: G̃/Φ₄ → |M_thr|²Φ₂ = A_old(2m_X)⁴S_out
for the 4→2 channels A and B (section 4.3).

Channels (incoming first; S_in, S_out), the sequential chain that the subtraction removes, and the energy where it opens:

| key | channel | S_in, S_out | √s_th/m_X | chain opens at √s/m_X (physical [D], see below) | removed chain (network processes) |
|---|---|---|---|---|---|
| A | X X X̄ X̄ → X X̄ | 4, 1 | 4 | 4.449 | X X̄ → s then s X X̄ → X X̄ (E); two pairs: twice X X̄ → s then s s → X X̄; nested decay/inverse-decay chains; X X̄ → X X̄* then X̄* X → s |
| B | X X X X̄ → X X | 6, 2 | 4 | 4.449 | X X̄ → s then s X X → X X (E2); elastic + inverse decay |
| C | X X X̄ X̄ → s s | 4, 2 | 4.899 | 4.899 | X X̄ → s then s X X̄ → s s (G); twice X X̄ → s then s s → s s |
| D | X X X̄ → X s | 2, 1 | 3.449 | 3.449 | X X̄ → s then s X → s X; X X̄ → X X̄* then X̄* X → s |
| E | X X̄ s → X X̄ | 1, 1 | 4.449 | 4.899 (= 2m_s) | X X̄ → s then s s → X X̄ (G̃_sub = G̃_full to 1–3 % below 4.899) |
| E2 | X X s → X X | 2, 2 | 4.449 | none | – (no neutral pair, no factorising X-line chain; G̃_sub = G̃_full at every energy) |
| F | s s s → s s | 6, 2 | 7.348 | none | – (s s → s* off shell) |
| G | X X̄ s → s s | 1, 2 | 4.899 | 4.899 | X X̄ → s then s s → s s |
| G2 | X s s → X s | 2, 1 | 5.899 | 5.899 | X s → X s with the u-channel X pole |
| H | X X̄ s s → X X̄ | 2, 1 | 6.899 | 7.348 (= 3m_s) | X X̄ → s then s s s → X X̄ (L): the only tree chain (G̃_sub/G̃_full = 0.99–0.94 up to 7.32) |
| J | X X X̄ s → X s | 2, 1 | 5.449 | 5.449 | X X̄ s → X X̄* (E, real X̄) then X̄* X → s; X X̄ → s then s X s → X s (G2) above 5.899 |
| K | X X̄ s s → s s | 2, 2 | 6.899 | 6.899 (visible from 6.93) | X s s → X* s (G2, real X) then X* X̄ → s; X X̄ → s then s s s → s s (F) above 7.348 |
| L | s s s → X X̄ | 6, 1 | 7.348 | none (X lines at > 14 m_X, 1e-5) | – |
| N | s s s s → s s | 24, 2 | 9.798 | none | – |

The opening energies are the Coleman–Norton condition [D]: an internal line is on shell in the physical region only where the
amplitude factorises into two real processes (a t-channel line between one incoming and one outgoing leg, or any line whose
other side would be a 2→1 or 3→1 process with heavier initial state, never is); they are confirmed by G̃_sub/G̃_full dropping
from ≈ 1 at that energy in the tables. The interface carries them as `chain_threshold_physical_sqrt_s_over_mX`; its
`chain_threshold_sqrt_s_over_mX` is the first tabulated point with a line flagged inside the regulator window (|p² − m²| <
25 mΓ_reg), which for A, B, E, H lies below the physical opening (the flagged line is then off shell and the subtraction changes
G̃ by ≤ 1 %, an O(Γ_reg²/x²) regulator effect).

and the fourteen reverse 2→n channels (X X̄ → X X̄ X X̄, X X → X X X X̄, s s → X X X̄ X̄, X s → X X X̄, X X̄ → X X̄ s, X X → X X s,
s s → s s s, s s → X X̄ s, X X̄ → X X̄ s s, X s → X X X̄ s, s s → s s s s, X X̄ → s s s, X s → X s s, s s → X X̄ s s, X X̄ → s s s s).

MB thermal rate densities γ(T) in GeV⁴ at μ = 0: **subtracted remainder** (★ = resolved: |γ| > 2σ_MC and regulator spread < 50 %) **/ regulated full** (Γ_reg = 0.02 m_s, 0.05 m_X; not a physical rate above the chain threshold, see text). The 2→n reverses follow (same |M|², independent tabulation).

| n→2 | T = 0.3 m_X | T = 1.0 m_X | T = 3.0 m_X | T = 7.344 m_X | T = 12.0 m_X |
|---|---|---|---|---|---|
| A XXX̄X̄ → XX̄ | +1.6e-27★ / 6.6e-26 | +2.1e-21★ / 2.3e-19 | -6.5e-21 / 2.0e-16 | -4.4e-17 / 1.2e-14 | -2.1e-16 / 9.9e-14 |
| B XXXX̄ → XX | +1.9e-30★ / 1.2e-29 | +7.3e-24★ / 3.7e-22 | -1.2e-20 / 1.6e-18 | -2.2e-18 / 1.8e-16 | -3.0e-17 / 1.7e-15 |
| C XXX̄X̄ → ss | +2.0e-25★ / 1.2e-23 | +3.8e-19 / 2.9e-17 | +1.1e-16 / 1.2e-14 | -1.1e-15 / 4.5e-13 | -1.1e-14 / 3.6e-12 |
| D XXX̄ → Xs | -5.4e-22 / 5.9e-20 | -4.0e-18 / 4.6e-15 | +6.6e-16 / 5.2e-13 | -6.2e-14 / 9.8e-12 | -7.1e-13 / 6.0e-11 |
| E XX̄s → XX̄ | -3.9e-25 / 3.4e-23 | -4.3e-18★ / 1.2e-16 | -3.6e-15★ / 1.1e-13 | -1.7e-13★ / 6.6e-12 | -1.0e-12★ / 5.2e-11 |
| E2 XXs → XX | +2.6e-27★ / 2.6e-27 | +9.1e-20★ / 9.1e-20 | +4.1e-16★ / 4.1e-16 | +4.9e-14★ / 4.9e-14 | +4.6e-13★ / 4.6e-13 |
| F sss → ss | +5.5e-28★ / 5.5e-28 | +7.9e-18★ / 7.9e-18 | +1.9e-13★ / 1.9e-13 | +4.3e-11★ / 4.3e-11 | +5.4e-10★ / 5.4e-10 |
| G XX̄s → ss | -9.7e-23★ / 4.5e-21 | +1.1e-17 / 1.1e-14 | +6.1e-14 / 4.8e-12 | +5.2e-13 / 2.0e-10 | -1.2e-11 / 1.7e-09 |
| G2 Xss → Xs | +8.8e-28 / 2.5e-25 | +1.8e-18★ / 2.7e-17 | +1.5e-14★ / 9.4e-14 | +2.9e-12★ / 9.1e-12 | +3.4e-11★ / 8.4e-11 |
| H XX̄ss → XX̄ | -7.8e-33 / 1.1e-30 | -1.4e-22★ / 4.5e-21 | -7.7e-19★ / 4.2e-17 | -6.6e-17 / 5.3e-15 | -2.9e-16 / 5.0e-14 |
| J XXX̄s → Xs | +2.4e-30 / 1.7e-27 | +1.2e-21 / 1.9e-19 | +4.2e-19 / 7.0e-16 | -1.1e-15 / 6.6e-14 | -1.5e-14 / 5.6e-13 |
| K XX̄ss → ss | +4.3e-31 / 1.8e-29 | -3.1e-21 / 1.2e-19 | -8.4e-17★ / 2.4e-15 | -1.3e-14 / 5.5e-13 | -8.6e-14 / 6.9e-12 |
| L sss → XX̄ | +1.1e-28★ / 1.1e-28 | +4.2e-19★ / 4.2e-19 | +3.9e-15★ / 3.9e-15 | +4.7e-13★ / 4.7e-13 | +4.5e-12★ / 4.5e-12 |
| N ssss → ss | +2.5e-36★ / 2.5e-36 | +6.6e-23★ / 6.6e-23 | +3.0e-17★ / 3.0e-17 | +2.5e-14★ / 2.5e-14 | +5.2e-13★ / 5.2e-13 |

| 2→n | T = 0.3 m_X | T = 1.0 m_X | T = 3.0 m_X | T = 7.344 m_X | T = 12.0 m_X |
|---|---|---|---|---|---|
| XX̄ → XX̄XX̄ | -2.2e-28 / 3.2e-26 | -5.8e-22 / 1.7e-19 | -3.6e-19 / 1.8e-16 | -2.4e-17 / 1.2e-14 | -1.1e-16 / 9.5e-14 |
| XX → XXXX̄ | +1.5e-30★ / 1.2e-29 | +7.4e-24★ / 3.2e-22 | -2.2e-21 / 1.4e-18 | -1.6e-18 / 1.7e-16 | -7.5e-18 / 1.6e-15 |
| ss → XXX̄X̄ | +1.9e-25 / 1.2e-23 | +2.2e-19 / 3.0e-17 | +6.8e-17 / 1.2e-14 | -8.2e-16 / 4.7e-13 | -4.8e-14 / 3.7e-12 |
| Xs → XXX̄ | -1.8e-22 / 5.3e-20 | -8.0e-18 / 3.8e-15 | -1.6e-15 / 4.2e-13 | -8.1e-14 / 8.6e-12 | -5.9e-13 / 5.6e-11 |
| XX̄ → XX̄s | +2.5e-24★ / 2.3e-23 | -1.2e-18 / 1.1e-16 | -3.1e-15★ / 1.0e-13 | -1.7e-13★ / 6.4e-12 | -1.1e-12★ / 5.1e-11 |
| XX → XXs | +2.5e-27★ / 2.5e-27 | +8.4e-20★ / 8.4e-20 | +3.9e-16★ / 3.9e-16 | +4.7e-14★ / 4.7e-14 | +4.5e-13★ / 4.5e-13 |
| ss → sss | +5.3e-28★ / 5.3e-28 | +7.6e-18★ / 7.6e-18 | +1.8e-13★ / 1.8e-13 | +4.2e-11★ / 4.2e-11 | +5.3e-10★ / 5.3e-10 |
| ss → XX̄s | -7.2e-23★ / 4.4e-21 | +1.6e-16 / 1.1e-14 | +1.1e-13 / 4.7e-12 | +1.7e-12 / 2.0e-10 | -3.4e-12 / 1.7e-09 |
| Xs → Xss | -2.8e-27 / 2.5e-25 | +1.7e-18★ / 2.6e-17 | +1.4e-14★ / 9.1e-14 | +1.9e-12★ / 9.5e-12 | +2.0e-11★ / 9.2e-11 |
| XX̄ → XX̄ss | +2.4e-32 / 9.6e-31 | -1.3e-22 / 4.2e-21 | -9.6e-19★ / 4.1e-17 | -6.0e-17 / 5.1e-15 | -3.1e-16 / 4.9e-14 |
| Xs → XXX̄s | -8.8e-30 / 1.4e-27 | -2.9e-21 / 1.8e-19 | -2.0e-17 / 6.3e-16 | -9.1e-16 / 6.1e-14 | +5.7e-15 / 5.2e-13 |
| ss → XX̄ss | +1.3e-31 / 1.7e-29 | -4.5e-21 / 1.2e-19 | -1.1e-16 / 2.4e-15 | -1.1e-14 / 5.4e-13 | -6.2e-14 / 6.8e-12 |
| XX̄ → sss | +1.1e-28★ / 1.1e-28 | +4.1e-19★ / 4.1e-19 | +3.7e-15★ / 3.7e-15 | +4.6e-13★ / 4.6e-13 | +4.3e-12★ / 4.3e-12 |
| ss → ssss | +2.4e-36★ / 2.4e-36 | +6.3e-23★ / 6.3e-23 | +2.8e-17★ / 2.8e-17 | +2.4e-14★ / 2.4e-14 | +5.0e-13★ / 5.0e-13 |
| XX̄ → ssss | +3.0e-36★ / 3.0e-36 | +1.3e-23★ / 1.3e-23 | +1.4e-18★ / 1.4e-18 | +4.6e-16★ / 4.6e-16 | +6.0e-15★ / 6.0e-15 |

Reading the table: below the chain threshold the subtracted kernel is the whole kernel and is determined to 0.1–1 %;
above it the regulated full column is dominated by the on-shell chain (it scales as 1/Γ_reg, i.e. ×261 from Γ_reg =
0.02 m_s to the physical width for a single on-shell s line, ×261² for two) and the subtracted remainder is a small
difference, resolved where marked and otherwise bounded. The physical statement for the solver: **above √s =
m_s + 2m_X = 4.45 m_X the 4→2 and 3→2 rates are the sequential chains** (inverse decay, decay, the 2→2 and the 3→2
with a real s), which the network carries with their own energy dependence and Bose factors; the remainder is the
tabulated G̃_sub with its error and regulator bound.

### 3.4 s self-reactions

s s s → s (and s s → s) are kinematically forbidden on shell (equal masses); the number-changing s self-reactions are
s s s → s s (F) and s s s s → s s (N) with their reverses s s → s s s, s s → s s s s, all through the cubic λ_s3 =
6λ_φv and the quartic λ_s4 = 6λ_φ with full interference and no physical pole. F and N are the only n→2 channels
whose kernels are positive and regulator-free at every energy (table); the elastic s s → s s cross section is
2.5e-2 GeV⁻² at threshold (cubic-dominated: the λ_s3²/m_s² ≈ 0.54 "coupling" against λ_s4 = 0.18) falling to
λ_s4²/(32π s) at high energy. MB rate densities (GeV⁴) at T/m_X = 1, 3, 7.344, 12: s s s → s s 7.9e-18, 1.9e-13, 4.3e-11, 5.4e-10; s s s s → s s 6.6e-23, 3.0e-17, 2.5e-14, 5.2e-13; the pair-producing s s s → X X̄ 4.2e-19, 3.9e-15, 4.7e-13, 4.5e-12. For comparison the 2→2 s s → X X̄ rate is 1.7e-14, 1.5e-11, 8.9e-10, 7.0e-9 and the decay-rate density 3.3e-10, 1.0e-8, 8.3e-8, 2.3e-7 at the same T: the 3→2 and 4→2 scalar self-reactions are 10⁻³–10⁻² of the 2→2 ones at T ≥ 3 m_X and negligible below; they do not supply a large number-changing rate by themselves, as Codex asked to be checked.

## 4. The overlap subtraction and the old constants

### 4.1 Statement

Every tree amplitude of the network is a sum over diagrams whose internal lines are s (neutral X X̄ or s s subsets)
or X (charge ±1 subsets). A line is on shell in the physical region exactly when the sequence of real processes it
represents is kinematically allowed (Coleman–Norton); those sequences are themselves channels of the network
(1 ↔ 2, 2 ↔ 2, 3 ↔ 2 with real s or X). Counting both the resonant amplitude and the sequence double counts. The
treatment adopted: **narrow-width (real-intermediate-state) subtraction** — for every line S that can be on shell,
|D_S|² → Re(D_S²) in |M|², i.e. the on-shell delta piece is removed and the Hadamard finite part, all interference
(2AB Re D, the principal value) and all off-shell pieces are kept; the removed pieces are supplied by the explicit
channels with their own in-medium widths, Bose factors and chemical potentials. In the limit Γ_reg → 0 the
subtracted kernel is regulator independent (Hadamard); at finite Γ_reg it carries O(Γ_reg) corrections, measured by
the two regulators.

### 4.2 Implementation and the defect found on resumption

M is multilinear in its propagators. For a set T of lines that can coexist in one tree (pairwise disjoint or nested)
the coefficient N_T of ∏_{S∈T} D_S follows by inclusion–exclusion over amplitudes with subsets of T switched off
(the recursion offers per-line switches), and

    |M|²_sub = Σ_T (−2)^{|T|} ∏_{S∈T} (Im D_S)² |N_T|²        (T = ∅ gives |M|²; |T| ≤ 3 suffices here),

which is the replacement |D_S|² → Re(D_S²) for every flagged line. The interrupted attempt included only disjoint
pairs (the `add` term) and skipped **nested** lines (the chains X X̄ → s → X X̄*, X̄* X → s …): with those omitted the
"remainder" of D, C, G, J, K scaled as 1/Γ_reg (D at 4.0 m_X: −7.5e-8 ± 2e-9 at Γ_reg = 0.02 m_s against −1.8e-8 at
0.04 m_s; C near threshold −3.2e-10 against −2.2e-11). With the nested terms and the triples included these become
D: −0.9 ± 1.1e-9 / −0.2 ± 0.2e-9 (consistent with zero at both regulators), C: +3.7e-11 / +1.7e-11, A at 9.2 m_X:
0 ± 1.7e-13 / 0.6 ± 0.4e-13. The remaining regulator dependence of C near its threshold is the boundary effect: the
double on-shell configuration sits within Γ_reg of the phase-space edge (√s − 2m_s = 0.01 m_X at the first point
against 0.02 m_s² of width), where no finite-regulator MC resolves the finite part; the table flags these points. The
checked identity full − RIS = subtracted (section 3.1) holds for the 2→2 quadrature at the physical width.

### 4.3 The old threshold constants (gp_leftovers/amp42.py; lsrelic.py s route) [C]

At the 4→2 threshold (all four at rest, √s = 4m_X) no internal line is on shell: the pair invariant mass is 4m_X²
against m_s² = 6m_X², three-particle X lines have 9m_X² against m_X². **The old constants A_a = 0.08421 GeV⁻⁸,
A_b = 0.002362 GeV⁻⁸ contain no on-shell piece and no double count** (Codex's statement that none had been
demonstrated: CONFIRMED, now by construction: `amp.py` reproduces M_a = −0.0533994, M_b = 0.0126472 to 1e-15 and the
λ_X-only hand sum 12λ_X²/m_X²). Likewise the old s-route constant SV_SS = 1.2603e-4 GeV⁻² is exactly our
σv_Møller(s s → X X̄) at threshold (|M| = 0.026667, ⟨σv⟩ → 1.19e-4 at T = 0.03 m_X and 1.26e-4 at T → 0), with S_in = 2
correctly carried by the old ½ SV_SS n_s². What the old constants miss is the energy dependence:

| T/m_X | A: A_eff/A_a (subtracted) | B: A_eff/A_b (subtracted) | A, B: regulated-full/old (chain at Γ_reg = 0.02 m_s) |
|---|---|---|---|
| → 0 (threshold, G̃_sub/Φ₄ at δ = 1e-4, 2e-4; `finalise_checks.py`, independent seed) | 0.987 ± 0.003, 0.980 ± 0.003 | 1.005 ± 0.003, 1.008 ± 0.003 | – (kernels_n2.py self-test at δ = 1e-4: 0.986, 1.009) |
| 0.03 | 1.40 | 0.69 | 1.43, 0.69 |
| 0.1 | 3.7 | 0.35 | 25.7, 0.81 |
| 0.3 | 1.5 (0.025 of full) | 0.092 (0.16 of full) | 59.3, 0.58 |
| 1 | 0.006 (0 within error) | 0.001 (0 within error) | 0.61, 0.053 |
| 3, 7.344, 12 | 0 within error | 0 within error | 5.0e-04, 2.1e-04 at 3 m_X (n_eq⁴ outgrows the chain) |

A_eff ≡ S_in γ(T)/n_eq⁴ (A_old is its T → 0 limit; both channels reach it at threshold, with the energy dependence of
|M|² growing toward the s pole for A (M_a and the s-exchange add) and falling for B (M_b is a cancellation).

The old s route (ss → X X̄ with n_s = K(T)n₊n₋ and the threshold σv) therefore **overestimates the pair destruction
through the scalar reservoir by 1.6× at T = 0.3 m_X, 3.6× at m_X, 15× at 3 m_X, 68× at 7.3 m_X and 170× at 12 m_X**
(MB; the Bose factor of the same channel, 1.2–5.7, is a separate correction in the other direction), because
σ(s s → X X̄) falls as 1/s above threshold while the old constant is held fixed. The old 4→2 constants are the T → 0
limit of the subtracted kernel (A_eff/A_a = 1.40, 3.7 at T = 0.03, 0.1 m_X; A_eff/A_b = 0.69, 0.35; results/balance.json `nr_limits`) and
above T ≈ 0.1 m_X they stand in for a rate that is in reality the resonant chain through real s and X, 10³–10⁵ times
larger at the regulator width and ×261 more at the physical width: as Codex wrote, not a kernel for the injected
daughters, and a late limit only for T ≲ 0.03 m_X with the inventory bookkeeping through the explicit s species.

## 5. Validation (`run_mc_checks.py`, `run_balance.py`; `results/mc_checks.json`, `balance.json`)

| check | result |
|---|---|
| Amplitudes: recursion vs the five analytic 2→2 forms; charge violation; s → X X̄ | ≤ 4e-15; 0; M/λ₃ = −1 |
| Phase space Φ₃, Φ₄ (massless exact, massive Dalitz, map independence) | 1.4σ, 1.4σ, 0.4σ; weights vs density 2e-12 |
| MB thermal rates: lab-frame MC vs quadrature, six 2→2 channels × 5 T | 0.1–0.6 % (non-resonant columns) |
| Detailed balance, Bose, point kinematics: 18 forward/reverse pairs × 20 random on-shell configurations (all 14 n↔2 pairs, s ↔ X X̄, s s ↔ X X̄, X X̄ → X X̄, X s → X s), μ_s = μ₊ + μ₋ and μ₊ + μ₋ = 0 for inventory-changing channels | **360 cases, max residual 4.9e-15**, energy residual 3.4e-16 (Codex: 240 cases, 5.3e-16: CONFIRMED and extended) |
| Detailed balance, Bose, integrated: s → X X̄ (MC) vs X X̄ → s (exact 2-D quadrature) at 5 T × 3 (μ₊, μ₋) | 15 cases: ratios 0.9939–1.0021 (0.968 at 12 m, μ = 0: quadrature tail), all within 2.1σ |
| Detailed balance, Bose, integrated: s s → X X̄ vs X X̄ → s s by independent MC samplings, 5 T × 3 μ sets (2μ_s = μ₊ + μ₋) | ratios 0.906–1.046, |σ| ≤ 1.6 over 15 cases |
| Table-level MB detailed balance γ_{n→2}(T) from G̃ vs γ_{2→n}(T) from σ, 14 pairs × 6 T, regulated-full and subtracted | pole-free channels (F, N, L, E2): γ_{n→2}/γ_{2→n} = 1 ± 0.087 at every T (the interpolation systematic of the two 16-point grids); resonant channels, regulated-full column: within 0.37 for T ≥ m_X, 0.95–2.08 at T ≤ 0.3 m_X where the coarse grids straddle the resonance turn-on; subtracted remainders: consistent within their MC errors in 33 of 42 (T ≥ 1 m_X) cases, the rest are the same interpolation systematic on small numbers |
| Vanishing internal energy moment of each elastic kernel at a common T (4 channels × 4 T) | (Q₁+Q₂)/(|Q₁|+|Q₂|) ≤ 2.4e-15 |
| X–s energy exchange at T_X ≠ T_s (closure input; non-resonant kernel) | sign correct, Q_s(T_X = T_s) = 0 to MC noise; e.g. Q_s/(n_X n_s) = −3.5e-6 GeV⁻¹ at (0.5, 7.344), +3.0e-6 at (7.344, 0.5) |
| RIS normalisation (on-shell piece of X X̄ → X X̄ = decay rate density) | 1 − 1e-6 at every T |
| Finalisation spot checks (`finalise_checks.py`: closed forms and own quadratures against the interface, no track module imported) | tree masses, λ₃, Γ_vac exact; γ_dec^MB and ⟨m/E⟩ to 1e-16; Gondolo–Gelmini quadrature of the tabulated σ(s) reproduces ⟨σv⟩ to 0.1–0.3 % (X X at T = 1, 7.3; s s → X X̄ at 3; X s at 1); contact limits to 1e-8 at √s ≈ 1e5 m_X; high-T T_c = 18.09 vs 18.14; gross-forward Bose factors of the decay vs the net ones (3.2); 4→2 threshold limits (4.3) |
| Full Breit–Wigner brute force vs RIS + subtracted (10–1000× width) | section 3.1 |
| Late limit vs old constants | section 4.3 |
| Regulator spread of the subtracted n→2 kernels (0.02 vs 0.04 m_s) | tabulated per point and per T; see 3.3 |

Hard daughters on a thermal bath (`hard_particle`, Bose bath, per bath density; the relaxation input of the solver):
for E = 7.72e4 m_X the collision rate per density ⟨σv_rel⟩ is 8.2e-11 (X), 1.6e-10 (X̄), 2.1e-9 (s) GeV⁻² at T = 0.3 m_X,
and the energy-loss rate per density is 2.9e-7, 5.8e-7, 7.4e-6 GeV⁻¹; at the onset charge density 0.0189 GeV³ the
2→2 energy-loss time of a daughter is ~4e11 GeV⁻¹ ≈ 1e-3 H⁻¹(a_ref) [E, with the solver's densities to be used].

## 6. Interface (`results/rates_interface.json`, built by `build_interface.py`)

`meta` (model, units, symmetry and Bose conventions, the subtraction, provenance); `masses_T` (one-loop and daisy
rows: v_T, m_X(T), m_s(T), λ₃(T), λ_s3(T), decay_open, Γ_s,rest(T), Bose factor; T_c, T_close; validity strings);
`decay_s_XXb` (Γ_vac, formulas, MB rate table, gross-forward Bose/MB at 5 T, `bose_factors` with the gross/inverse/net conventions and the net-width Bose/MB); `channels_22` (σ(s) tables with all prescriptions and
the transport weight, MB γ and ⟨σv⟩ tables on 51 T, Bose/MB factors, high-energy law); `channels_n2` and `channels_2n`
(G̃(s) or σ(s) at both regulators with MC errors and the number of flagged lines, chain threshold and the chain
named, MB thermal tables on 45 T with error, `chain_threshold_physical_sqrt_s_over_mX` and the flagged-point `chain_threshold_sqrt_s_over_mX` with its note, regulated-full, second regulator, `remainder_resolved` flags and
`remainder_bound`); `hard_particle_on_bath`; `energy_exchange_Xs`; `checks`; `usage` (the solver recipe). The solver
can also import `code/tables.py` (interpolation and thermal averages), `xsec22.py` (σ, thermal rates with any
Params), `thermal.py` (masses, widths), `thermal_mc.py` (Bose integrals).

## 7. Limits

- Kernels at vacuum masses and couplings; the in-medium tables give the regime (1 % to 3 m_X, 5 % to 6 m_X) and the
  decay switch at 13.1 m_X; re-tabulation at masses(T) is not done (the code takes any `Params`).
- Thermal tables are MB with Bose factors at five temperatures and μ = 0; a solver with Bose statistics at μ ≠ 0 must
  integrate σ(s)/G̃(s) with its own f's (the MC engine does this for checks; cost is the solver's).
- The subtracted n→2 remainders are resolved only where flagged; near the chain thresholds and where the chain
  dominates they are bounds; the regulated-full columns are not physical rates (they scale with 1/Γ_reg).
- Tree level throughout; no thermal widths in the X s u-pole (regulated as a finite part); no Landau damping in the
  kernels; one loop in the masses; μ = 0 in the loops.
- The hard-particle and energy-exchange tables are for equilibrium baths at μ = 0; the solver's actual densities
  rescale them (rates per density are given).
- The MC tables have statistical errors of 0.1–1 % below the chain thresholds and larger above; the 16-point grids
  are supplemented by 13 extra points for A, B, D, E around 4.25–7.5 m_X; interpolation between points was spot-checked
  to 5–30 % (section 3.3 supplement) before the supplement was added.

## 8. Reproduce

    cd research/reality_framework/grut_formation_review_20261025/tracks/kin_rates/code && ./run_all.sh

(~2.5 h on 4 cores, dominated by `run_kernels.py`; module self-tests `python3 -B model.py | amp.py | phasespace.py |
xsec22.py | thermal.py`). Scratch outputs in
/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/kin_rates/prod/.

Finalisation checks (~40 s): `python3 -B finalise_checks.py` (writes results/finalise_checks.json; the 4→2 threshold rows were
produced by the inline MC call recorded in results/finalise_checks.json `o_threshold_limit_4to2`). The interface alone is
rebuilt from the stored tables in 9 s: `python3 -B build_interface.py ../results/kernels_n2.json ../results/xsec_2n.json
../results/supplement_n2.json ../results/balance.json`.

## 9. Finalisation pass (third session; documentation and consistency only)

Found by re-reading NOTES.md against results/ and by the spot checks; no numerical result changed, and the regenerated
rates_interface.json is numerically identical to the previous one (0 differences over every numeric entry; only strings and
the new fields differ):

1. The s → X X̄ Bose factors quoted in 3.2 (1.055–5.86) were from a stale run; results/mc_checks.json has 1.024–15.26 for the
   gross forward rate. The spot check showed that the track's `width_s_medium` (net Weldon factor 1 + f₁ + f₂) and the MC
   (gross factor (1+f₁)(1+f₂)) are different quantities, and the interface's recipe line ("Γ from formulas + inverse by
   detailed balance") could be read as a double count; the three conventions are now stated (3.2, `decay_s_XXb.bose_factors`,
   `usage`), with the net-width Bose/MB values added (1.024, 1.49, 2.98, 5.25, 6.87).
2. The BE/MB column of the 2→2 table (3.1) was refreshed to results/mc_checks.json (changes up to 30 % at 12 m_X for X X).
3. The "chain opens at" column (3.3) mixed physical openings with the first flagged MC point: E 4.449 → 4.899, E2 5.76 → none,
   H 6.899 → 7.348; the chains named for J and K were not the ones that open first. Both energies are now in the interface.
4. The late-limit table (4.3) was refreshed to results/balance.json (which includes the supplement points): A_eff/A_a at
   T = 0.1, 0.3 m_X 3.6 → 3.7, 1.3 → 1.5; A_eff/A_b at 0.1 m_X 0.29 → 0.35; the unsupported T = 0.01 row was dropped; the
   threshold row was reproduced with an independent seed.
5. Grid counts: 51 (not 53) temperatures in the 2→2 thermal tables, 45 (not 47) in the n ↔ 2 tables; the detailed-balance
   ranges of section 5 were refreshed.

Credits: Dolan–Jackiw finite-temperature potential; Parwani / Arnold–Espinosa resummation; Berends–Giele recursion;
real-intermediate-state subtraction as used in leptogenesis Boltzmann networks (Kolb–Wolfram; Giudice et al.);
Coleman–Norton theorem for physical-region singularities; Kolb–Turner Bessel-series thermodynamics.
