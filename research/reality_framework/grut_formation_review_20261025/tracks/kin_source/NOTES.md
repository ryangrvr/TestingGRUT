# Track kin_source: source terms for the kinetic case-L calculation at the heavy spot point

26 October 2026 · Claude, issue #13 (BRIEF_KINETIC.md) · same-session AI work, not outside replication.

Codex packets read only through `git show`/`git archive` into the scratchpad: `grut_kinetic_contract_20261026` at
9877518 (KINETIC_SPEC.md, MATHEMATICS.md, DECAY_SOURCE.json, THERMODYNAMIC_CONTROL.json, code read, not run) and
`grut_heavy_vector_budget_20261025` at 225028b (INPUTS.json, BACKGROUND_RESULTS.json, VECTOR_RESULTS.json, BUDGET
code read). No Codex code was executed. Every number below is from this track's own code unless labelled [Q].
Own earlier tracks reused, with the reuse stated: hw_gates gate1_vector.json (our 49-mode spectrum, read as data);
the coarse g* table of ls_threshold (copied); the planar-orbit action method of gp_leftovers budget.py
(re-implemented); the certified free-mode fractions and bins of tracks/spectrum (read as data); the Y eigenmode,
width and three-body-spectrum results of sp_yledger [Q own]. Nothing was committed.

**Scope.** Source terms only: no Boltzmann solution, no f_sym, no threshold, no formation kernel, no viability.
One supplied tree-level model; nothing here is evidence that the fields exist. Labels: [D] derived, [C] computed,
[Q] quoted, [E] estimate, [A] argued.

## 1. Conventions and the translation into ours

Codex's convention (DECAY_SOURCE.json `number_convention`): a = 1 at the end of inflation, comoving numbers
N = a³n/H_i³, one bosonic degree per daughter with n = ∫dp p²f/(2π²), comoving momentum q = ap. This is the
convention of our hw_gates/sp_yledger/ls_* releases (N_X = comoving Noether charge/H_i³, E_X = a⁴ρ_X/H_i⁴), so
**no numerical translation is needed between the two vector sources**; the differences are in normalisation
(section 3). The translations the solver needs (all [D], `meta.translation` in the interface):

| from comoving per H_i³ … | multiply by |
|---|---|
| per net charge | 1/\|N_X\| = 3.5866e-4 |
| per visible entropy | Y_Q/\|N_X\| = 1.7321e-12 |
| physical density at a (GeV³) | H_i³/a³; at a_ref: 6.7766e-6 |
| conformal energy (H_i units) → R | 1/E_X = 2.6051e-5 |
| conformal energy → physical ρ at a (GeV⁴) | H_i⁴/a⁴; at a_ref: 1.2823e-7 |
| q/H_i → p/m_X at a_ref | H_i/(a_ref m_X) = 0.20479 |
| pairs per net charge → f (rest mass / DM benchmark) | 2 m_X Y_Q/B = 2.02855 |

Comoving phase-space occupation per degree: f(q) = 2π²(dN/dln q)/q³. For the vector daughters it peaks at
1.5e-11 (e = 1e-7): the 1+f factors on the daughters are 1 to that accuracy.

## 2. The vector source A → X X̄ (own ODE)

**Decay law [D].** m_A = 6eh (charge-6 Σ Higgs), Γ_A = e²m_A/(48π)(1 − 4m_X²/m_A²)^{3/2}, p* = (m_A²/4 − m_X²)^{1/2},
H = H_i/a² (g* constant through every decay: T_vis > 9 TeV). In ln a, with γ = Γ_A/H_i and t H_i = a²/2,

    dN_A/dln a = −γ a² N_A,   dN_±/dln a = +γ a² N_A,   N_A(a) = N_A0 exp[−γ(a² − 1)/2].

The parent is at rest (p_A/m_A < 1e-5 at every decay, section 3), so each daughter is born at q = a p* and keeps
it: **the time profile and the comoving spectrum are the same function**, dN_±/dln q|_{q = ap*} = γa²N_A(a), i.e.
the Rayleigh form N_A0 x² e^{−x²/2} with x = q/q₀ = a/a_unit, q₀ = p*/√γ, a_unit = γ^{−1/2}, ⟨q⟩ = q₀√(π/2),
⟨q²⟩ = 2q₀². The 2D source for the solver is S_±(a, q) = (dN_±/dln a) δ(ln q − ln(ap*/H_i)).

Integrated with DOP853 (rtol 1e-12) from a = 1 to 10 a_unit, state (N_A, N_±, Σq, Σq²); cross-checked by a second
integration in cosmic time. Results (`results/decay_source_own.json`; Codex values from DECAY_SOURCE.json):

| | e = 1e-7 band | e = 1e-7 broad | e = 1e-6 band | e = 1e-6 broad |
|---|---:|---:|---:|---:|
| N_A (per sign) [Q] | 2.0650977 | 2.4691556 | 0.6526450 | 0.7712576 |
| m_A (GeV), Γ_A (GeV) | 3.2003909e6, 2.1223251e-10 | same | 3.2003909e7, 2.1223251e-7 | same |
| a_unit, ⟨a⟩ | 2.50663e10, 3.14159e10 | same | 7.92665e8, 9.93459e8 | same |
| q₀/H_i, ⟨q⟩/H_i | 300795.4, 376991.1 | same | 95119.9, 119215.1 | same |
| pairs per net charge | 7.406624e-4 | 8.855807e-4 | 2.340759e-4 | 2.766172e-4 |
| R (energy moment, ODE) | **40.56262** | **48.49912** | **4.05380** | **4.79054** |
| f_rest if unprocessed | 0.0015025 | 0.0017964 | 0.00047483 | 0.00056112 |
| ⟨p⟩/m_X at a_ref (no collisions) | **77200.86** | 77200.86 | 24413.06 | 24413.06 |
| energy per net charge at a_ref (m_X) | 114.36 | 136.74 | 11.43 | 13.51 |
| T_vis at ⟨a⟩ (GeV) | 9806 | same | 3.10e5 | same |
| parent survival at a_ref | e^{−39515} | same | e^{−3.95e7} | same |

Codex's e = 1e-7 numbers: pairs/net 7.406624e-4, R 40.5626, f_rest 0.0015025, ⟨p⟩/m_X 77200.9: **all CONFIRMED**
(largest relative difference over 112 compared values, including 5 quantiles and the 161 spectrum samples per e:
5.0e-9, which is the trapezoid residual of Codex's own sample grid; the ODE-vs-analytic residuals are ≤ 8e-10).
The finite-start factor e^{γ/2} − 1 ≈ 8e-22 is kept and invisible. The daughters' exact conformal energy at a_ref
(massive) exceeds Σq by 1.3e-10; they become nonrelativistic at a ≈ 5.4e17, 7.7e4 times after a_ref.

**R scan [diagnostic].** For each case the interface lists N_A multipliers for R = 1 … 90 (including 40.5626 and
48.499) at fixed spectral shape; it is a rescaling experiment, not a cosmological point.

## 3. The second (broad) source: our hw_gates grid

Our gate-1 spectrum (own longitudinal-Proca modes, 49 nodes k = 1e-5 … 1 H_i, profile P1 = evolving Higgs norm,
inflationary norm F_inf = 266.06 H from our aligned preparation) integrated in ln k:
N_A = 2.46916 (P0: 2.46849). Split against Codex's band k/k* = 0.2–5 (k* = √(m_A/H_i) = 4.90e-3) [C]:

| | in band | below 0.2 k* | above 5 k* | full |
|---|---:|---:|---:|---:|
| N_A, e = 1e-7 | 2.0204 (0.978 of Codex's 2.0651) | 0.0403 | 0.4081 | 2.4692 (1.196 × Codex) |
| N_A, e = 1e-6 | 0.6384 (0.978 of 0.6526) | 0.0127 | 0.1200 | 0.7713 (1.182 ×) |

Mode-by-mode against Codex's 'actual' modes at the same k: ours/Codex = 0.963–0.981 (number-weighted 0.974). The
2–4% in-band deficit is consistent with the inflationary Higgs norm (Codex's conditional minimum F_inf = 263.45 vs
our 266.06: m_inf² ratio 0.980) [A]. The extra 20% of the broad grid sits almost entirely **above** the band. Codex's
band is a subset, as its packet states; the broad source is what gate 1 used. Normalisations are not mixed: the
broad case uses its own N_A with the fixed point's E_X (R = 48.499); against hw_gates' own E_X = 38398.22 and
H = 1.33350e11 it is R = 48.484 (case `e1e-7_broad_hwnorm`).

Parent motion [C]: number-weighted k_rms = 0.062 H_i (0.11 at e = 1e-6); the daughter comoving momentum is smeared
by δq = (k/2)cos θ*, rms 6e-8 q₀ (max 1.7e-6 q₀ at k = 1); the decayed fraction while any parent of the grid is still
relativistic (k ≤ 1, p/m_A > 0.1) is < 1.4e-10. So the broad source has exactly the Rayleigh spectrum of section 2
with its own N_A; the parent k-spectrum is tabulated in the interface for completeness.

## 4. Mediator debris as late sources (R = 1.14)

Own numbers from INPUTS.json couplings (`results/mediator_source.json`):

- **Y_L → XXX or X̄X̄X̄** (half each; the light real eigenmode) [D]: m_L = (M_Y² − √2μh)^{1/2} = 2.9668e12 GeV
  (22.25 H_i), Γ_L = 3κ²m_L/(256π³) = 1.3414e4 GeV, γ_L = 1.006e-7, a_unit = 3153, ⟨a⟩ = 3952 (T_vis = 7.8e10 GeV).
  Parents N_L = N_YH + N_YL = 0.351718 (Y_H cascades first). Per sign 1.5 daughters per decay: N_± = 0.52758
  (1.892e-4 per net charge). Rest-frame spectrum dN/dE = 8E/m_L² on [0, m_L/2] (contact |M|², massless daughters),
  ⟨E⟩ = m_L/3; q = aE. Comoving spectrum derived and checked numerically (convolution vs closed form 4e-13):
  dN/dq = (4qγ_L/m²)E₁(2γ_Lq²/m²), ⟨q⟩ = m⟨a⟩/3 = 29306 H_i (p/m_X = 6002 at a_ref). **R_Y = 0.80557**
  (Codex 0.805570: CONFIRMED). Time profile from the same ODE form as section 2.
- **Cascade vectors** Y_H → Y_L + A_L: Γ ≈ μ²(m_H² − m_L²)/(8π m_H³) = 4.7e9 GeV (Goldstone equivalence [E]),
  equal to H at a = 5.3; E_A = 1.5e12 GeV at emission, q_A ~ 60 H_i, nonrelativistic from a ~ 2.5e6, i.e. cold for
  four decades before decaying with the primordial vectors: same Rayleigh spectrum with N_A = N_YH = 0.017253.
  R = 0.33888 (e = 1e-7), 0.10716 (e = 1e-6): CONFIRMED against 225028b. Caveat [Q]: N_YH is the late undamped
  classical count (a = 30–40); the population when Γ = H (a ≈ 5) is a bracket, not re-derived here (≤ 0.8% of R).
- **σ → AA**: 6.3e-9 vectors, R = 1.2e-7.
- Totals: R_med = 1.14445 (e = 1e-7), 0.91273 (e = 1e-6); 1.954e-4 pairs per net charge in all.

**Time ordering** (comoving): Y_H → Y_L A at a ≈ 5; Y_L decays at ⟨a⟩ = 4.0e3; vectors at 3.1e10 (e = 1e-7) or
9.9e8 (e = 1e-6); case-L coordinate a_ref = 7.047e12 (T_vis = 45.7 GeV); vector daughters nonrelativistic at 5.4e17.
All sources are complete before onset by 220× (e = 1e-7) to 7000× (e = 1e-6) in a. Every source acts while the X
sector is in the φ = 0 phase (massless X, pure quartic condensate); the s state and the vacuum m_X exist only after
the transition at a_ref [D, phase history otherwise open: ls_threshold one-loop T_c = 18–42 m_X [Q]].

## 5. The no-source control: the release's own population

**Onset coordinate [D, C]:** x_c = (λ_φv²/2g_X²)^{1/2} = 1.13177 GeV, n_c = 2√(2λ_X)x_c³ = 0.129665 GeV³,
ρ_c = 3λ_Xx_c⁴ = 4.9222e-3 GeV⁴, V0 = λ_φv⁴/4 = 5.4691e-3 GeV⁴ (= m_X²x_c²/2 exactly), n_Q = |ε|n_c = 0.0188945 GeV³,
a_ref = H_i(|N_X|/n_Q)^{1/3} = 7.04671e12 (from E_X: 7.04686e12; ratio 0.99998, the window-mean ε vs N_X/N_c(E_X)),
T_vis = 45.7 GeV (46.9 with the g_s node at 86.25), H = 1.037 H_i/a_ref² with g*. Energy per net charge at onset:
ρ_c/n_Q = 2.819 m_X, V0/n_Q = 3.132 m_X, total **5.951 m_X**. Condensate comoving amplitude (E_X/λ_X)^{1/4} = 78.7 H_i,
frequency scale 3.5 H_i, i.e. p ≈ 0.72 m_X at a_ref. Matching the release means these numbers, not a temperature.

Labelled alternatives for the free population (none is measured; the solver should run them as controls):

| label | assumption | X̄ per net charge | quanta per net charge | energy |
|---|---|---:|---:|---|
| R1 coherent | homogeneous condensate until an explicit fragmentation/formation kernel | 0 (gas) | 0 | E_X conformal |
| R2 adiabatic orbit count [D, C] | a³n and a³J_r invariant through a slow transition; planar orbit survives | **2.586** (ε_NR = 0.1620) | 6.17 | 6.17 m_X per net charge (3.7% above the cap: sudden V0 and adiabaticity are not both exact) |
| R3 energy cap [D] | everything cold | ≤ 2.476 (≤ 0.910 without V0) | 5.95 | 5.95 m_X |
| R4 fragmentation template [A] | the release fragments like our certified charge-saturated object (tracks/spectrum, E0/Q0 = 0.9995 m_X) | 0.00657 (⟨ω⟩ 1.69 m_X) + 0.0028 s (⟨ω⟩ 3.24) | 0.66 free X (⟨ω⟩ 1.038) + 0.35 cores | carries only 17% of the release energy: a spectral shape, not a population |
| R5 instant thermalisation [C] | fully thermal X X̄ (s) gas at ρ_c + V0 and n_Q, condensate at ω | 0.185 (no s), 0.134 (with s) at T/m = 3.68 / 3.35; condensed 85–87% | σ/net = 1.85–2.0 | the old proxy's R = 0 start |

R2 reproduces gp_leftovers at ε = 0.1 (f_adiabatic 8.094 vs 8.09) and the harmonic control to 1e-13; the NR count
law gives 2.93. **The release's own symmetric population before number-changing reactions is 2.5–2.6 X̄ per net
charge, 3500 times the vector daughters' 7.4e-4; the vectors supply energy (R = 40.6) not number.** V0 handling:
sudden transition → V0 into s at rest, 1.28 s per net charge; adiabatic zero-pressure coexistence → none
(gp_leftovers); both are listed.

## 6. Equilibrium comparators (own Bose thermodynamics; checks of THERMODYNAMIC_CONTROL)

Bessel-series thermodynamics checked against direct momentum quadrature (16 states, ≤ 9e-15). Zero-μ reservoir
rows at T/m = 0.3, 1, 3, 12, 30: scalar fraction 0.021216, 0.235849, 0.318948, **0.332366**, 0.333177 and n_s/n_-
0.023050 … 0.995928 — Codex's rows to ≤ 3e-14: CONFIRMED (the "33.24% at T = 12 m" claim). Fixed-energy states at
R = 40.5626 (ρ = ρ_c(1+R) + V0 = 0.21005 GeV⁴, μ_c = ω, linear ωn_c): T/m = **8.11769** (no s) → **7.34390** (with s),
condensed 0.004405 → 0.007072, scalar fraction 0.32765, all species n, ρ, P to ≤ 2e-13: CONFIRMED. Added rows
(interface `equilibrium_comparators.states`): R = 0 and 48.499, each with/without s and with/without V0; e.g.
R = 48.5 with s: T/m = 7.665, X̄/net 1.97, s/net 2.18, σ/net 24.6. m_s = 0.226355 GeV, trilinear 2g_X²v = 0.018482,
vacuum Γ(s → XX̄) = 1.73328e-5 GeV, E − Qω = 1148.1075 GeV: CONFIRMED. These are comparators, not populations.

## 7. Verdict table (`results/verdicts.json`, 144 rows)

| Claim (Codex 9877518 / 225028b / own) | Codex | This track | Verdict |
|---|---|---|---|
| e = 1e-7: pairs/net, R, f_rest, ⟨p⟩/m_X at a_ref | 7.406624e-4, 40.5626, 0.0015025, 77200.9 | 7.406624e-4, 40.56262, 0.0015025, 77200.86 (own ODE) | CONFIRMED (≤ 5e-9) |
| e = 1e-7: a_unit, ⟨a⟩, q₀, ⟨q⟩, a_ref (charge, energy), 5 quantiles, early-decay fraction, survival | as DECAY_SOURCE | identical to ≤ 5e-9 | CONFIRMED |
| e = 1e-7 and 1e-6: 161 spectrum samples dN/dlnq and dR/dlnq | samples | ODE on the same grid: max 1.3e-9 | CONFIRMED |
| e = 1e-6: pairs/net, R, ⟨p⟩/m_X | 2.340759e-4, 4.05380, 24413.1 | same | CONFIRMED |
| R_energy vs R_parent_mass agree to 1.6e-15 | stated | 2p*/m_A − 1 = −1.7e-15 | CONFIRMED |
| Reservoir table; 33.24% at T = 12 m; 8.118 → 7.344; condensed charge; species moments | THERMODYNAMIC_CONTROL | reproduced ≤ 2e-13 (own series + quadrature) | CONFIRMED |
| R_Y 0.80557, R_casc 0.33888, R_σ 1.2e-7, mediator total 1.1444 (225028b, hw_gates) | as stated | 0.80557, 0.33888, 1.24e-7, 1.14445 | CONFIRMED |
| hw_gates broad source N_A ≈ 2.47, R ≈ 48.5 [own] | 2.47, 48.5 | 2.46916; 48.499 (fixed E_X) / 48.484 (own E_X) | CONFIRMED |
| Codex band is a conservative subset of the spectrum | stated | band = 0.818 of our grid; in band ours/Codex 0.978 (F_inf 266 vs 263) | CONFIRMED, quantified |
| "Thermal initialisation from R alone assumes the processing" | argued | daughters: f ≤ 1.5e-11, ⟨p⟩ = 7.7e4 m_X at a_ref, 220× before onset; release: 2.5–2.6 X̄/net cold quanta | CONFIRMED as framing (no solution here) |

No CORRECTED items. Not checked: Codex's own numerical residuals (its ODE/quadrature errors), the 240 Bose
detailed-balance samples and the finite-Q derivative test (identities, not inputs).

## 8. The interface file (`results/source_interface.json`, 833 kB)

- `meta`: conventions, translation table, fixed point, background (T_RH = 3.08e14 GeV, s_RH), time ordering.
- `vector_sources.{e1e-7_band, e1e-7_broad, e1e-7_broad_hwnorm, e1e-6_band, e1e-6_broad}`: parent law and
  constants; per-sign totals (N, pairs/net, Σq, R, ⟨q⟩, ⟨q²⟩, ⟨p⟩/m_X at a_ref, survival); quantiles; one 241-point
  grid in x = a/a_unit = q/q₀ with a, tH_i, q/H_i, N_A remaining, dN_±/dln a (= dN_±/dln q), dN_±/dt (per H_i⁻¹),
  dE_conf/dln q, dR/dln q (both signs), f(q), cumulative fraction, p/m_X at a_ref; the analytic forms; the 2D
  injection rule; the R-scan multipliers; for the broad cases the parent k-spectrum and smearing bounds.
- `mediator_debris`: Y_L parents, daughters, time profile (241 points), q spectrum (721 points), 2D rule; cascade
  and σ vectors (their tables under `vector_sources_extra`); totals per e.
- `release_population`: onset numbers; alternatives R1–R5 with their assumptions; V0 handling.
- `equilibrium_comparators`: reservoir rows and the 12 fixed-energy states.
- `sums_for_checks`: per-case N_± and R, mediator totals, release N_X, E_X, n_Q, ρ_c, V0.

## 9. Limits

- The decay law assumes instant reheating with constant g* through the decays (true: T_vis > 9 TeV) and vacuum
  daughter kinematics; collisions are the solver's business. H at a_ref is 3.7% above H_i/a² (g* = 86–94).
- Branching A → XX̄ = 1 (no kinetic mixing, no lighter charged state) [Q gp_leftovers]; width tree level.
- The broad source carries our own inflationary preparation (F_inf = 266.06 H); the 2% in-band offset to Codex is
  attributed, not resolved. Above-band modes k > 5k* are ours alone (16.5% of the number).
- Cascade N_YH and the Y_H → Y_L A width are a late classical count and an equivalence-theorem estimate [Q, E].
- The release population is bracketed, not measured: R2 needs adiabaticity through the φ transition and no
  fragmentation; R4 is a template for a different object; R5 is the proxy's assumption. The φ phase at onset,
  in-medium masses and the fate of V0 are open and handed to the solver as labelled inputs.
- The ε of the fixed point is the a = 30–40 window mean (−0.145717); a = 30–100 gives −0.145669 [Q 225028b].

## 10. Reproduce

    cd research/reality_framework/grut_formation_review_20261025/tracks/kin_source/code && ./run_all.sh

(about 10 s; reads `git show 9877518:…/DECAY_SOURCE.json`, `…/THERMODYNAMIC_CONTROL.json` and
`git show 225028b:…/VECTOR_RESULTS.json` into the scratchpad for the comparisons only). Individual steps:
`decay_source.py`, `broad_spectrum.py`, `mediator_source.py`, `release_population.py`, `build_interface.py`.

Credits: longitudinal-vector production from inflation (Graham–Mardon–Rajendran); Bose thermodynamics by Bessel
series (Kolb–Turner); Q-ball formation and fragmentation (Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki).
