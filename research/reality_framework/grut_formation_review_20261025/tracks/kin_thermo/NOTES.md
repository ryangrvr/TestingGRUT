# Track kin_thermo: thermodynamics with the scalar reservoir, the object sector and the formation alternatives

26 October 2026 · Claude, issue #13 (BRIEF_KINETIC.md) · same-session AI work, not outside replication.

Codex packets read only through `git show`/`git archive` into the scratchpad: `grut_kinetic_contract_20261026` at
9877518 (KINETIC_SPEC.md, MATHEMATICS.md, THERMODYNAMIC_CONTROL.json, DECAY_SOURCE.json; code read, **not run**),
`grut_heavy_vector_budget_20261025` at 225028b (INPUTS.json, BUDGET_RESULTS.json) and `grut_core_conversion_20261023`
at bb211d1 (MATHEMATICS.md, REPORT.md). Every number below is from this track's own code unless labelled [Q]
(quoted) or "reproduction" (our earlier code run unchanged). Labels: [D] derived, [C] computed, [Q] quoted,
[E] estimate, [A] argued, [R] recalled. Nothing was committed.

**Reuse of our earlier tracks, stated.** (i) `code/old_proxies.py` imports `tracks/ls_threshold/code/lsrelic.py` (and its
`lsthermo.py`) unchanged, with no bytecode written there, to produce the reference proxy curves of section 5; those
numbers are labelled "reproduction". (ii) The conversion cross section σv(p) of the capture channel is read as data
from `tracks/conversion/results/rate_scan.json` and `scatter_spot.json` (our own three-channel solver of 25 October,
which CONFIRMED Codex's two bb211d1 spot values to 4e-6) [Q own]. (iii) The profile equations and the quadrature
formulas are the same model equations as `tracks/conversion/code/profile.py` and `tracks/xcharge/code/profile.py`;
the solver here is a third formulation (section 2). (iv) The Bessel-series thermodynamics is the same Kolb–Turner
series as `lsthermo.py`, re-implemented with the s species, the analytic susceptibilities and the tail control, and
cross-checked against an independent direct quadrature; the `lsthermo.py` comparison in `results/checks.json` is
labelled reproduction.

**Scope.** State functions, object energetics, exchange kernels and formation rules for one supplied tree-level model;
callables and tables for the kinetic solver. No Boltzmann solution, no f_sym, no threshold, no viability, no
observational statement. Nothing here is evidence that the fields exist.

## 1. Thermodynamics with the explicit scalar reservoir (`code/thermo.py`)

**Species.** X (charge +1, μ₊), X̄ (−1, μ₋), the neutral radial quantum s (μ_s), one bosonic degree each, Bose
statistics, m_s = √(2λ_φ) v = √6 m_X = 0.226355 GeV. Units m_X = 1 (t = T/m, u = μ/m). For one species of mass r
[D, Kolb–Turner series]:

    n = (r²t/2π²) Σ_k e^{ku/t} K₂(kr/t)/k,   P = (r²t²/2π²) Σ e^{ku/t} K₂/k²,   ρ = (r²t/2π²) Σ e^{ku/t}[3tK₂/k² + rK₁/k],
    s = (ρ + P − Σ_i u_i n_i)/t,
    ∂n/∂u|_t = (1/t) Σ k n_k,   ∂n/∂t|_u = (1/t²) Σ k(ρ_k − u n_k),   ∂ρ/∂u|_t = (1/t) Σ k ρ_k,
    ∂ρ/∂t|_u = (1/t²) Σ k(E2_k − u ρ_k),   E2_k = (r²/2π²) e^{ku/t}[12t³K₂/k³ + 5rt²K₁/k² + r²tK₀/k]

(from f(1+f) = Σ_k k e^{−k(E−u)/t}). Scaled Bessel functions; the series is cut at e^{−40}; within (r−u)/t < 2.7e-5 a
direct momentum quadrature replaces it (∂n/∂u diverges at u → r as it must). **Closures** offered, none imposed:
pairs μ₊ + μ₋ = 0; fast s (s ↔ XX̄ fast) μ_s = μ₊ + μ₋; full equilibrium both (μ_s = 0); caps μ₊ ≤ ω (objects),
μ₋ ≤ ω (anti-objects), μ_s < m_s. **Pair inventory** P_inv = n₋ + n_s (+ anti-condensate), conserved by s ↔ XX̄,
lowered by ss → XX̄ and 4 → 2 (Codex MATHEMATICS.md stoichiometry, agreed).

**State solvers** (inverse maps): `state_equilibrium(n_Q, ρ)` (pairs + fast s; the closure of Codex's control and of
the old proxies, linear object energy e_cap n_c); `state_two_potential(n_Q, P_inv, ρ)` (kinetic equilibrium with the
actual pair inventory, μ_s = μ₊ + μ₋, objects and anti-objects at the caps); `state_three_potential(n_Q, n₋, n_s, ρ)`
(s chemically decoupled; Newton with the analytic 4×4 Jacobian). `susceptibility()` returns
∂(n_Q, n₋, n_s, ρ)/∂(t, u₊, u₋, u_s) and the fast-s 3×3 block.

**Checks [C]** (`results/checks.json`): series vs independent quadrature, 7 quantities × 28 states: 2.0e-15; Gibbs–Duhem
∂P/∂μ = n, ∂P/∂T = s and the analytic ∂n/∂T against finite differences: ≤ 4e-11; the two-potential solver at the
equilibrium inventory reproduces the equilibrium solver to 1e-12; the three-potential Newton round trip to 1e-13.

**Codex's control, reproduced [C]:** reservoir rows at μ = 0, T/m = 0.3, 1, 3, 12, 30: n_s/n₋ = 0.023049850883,
0.450738727920, 0.845563073722, 0.981069806273, 0.995928398062 and scalar energy fraction 0.021215613260,
0.235849436772, 0.318947809560, **0.332365630564**, 0.333177301659 — Codex's rows to ≤ 1e-12. Matched energy at
R = 40.56262 (ρ = ρ_c(1+R) + V₀ = 0.2100491 GeV⁴, n_Q = 0.0188945 GeV³, cap ω = 0.861195 m, linear ωn_c):
T/m = **8.117690604 → 7.343904195**, condensed 0.004404564 → 0.007072206 GeV³, pair inventory 0.04451798 →
0.06850712 GeV³ (s holds 0.03604098), scalar energy 0.327647 of the total (Codex's definition) = 0.328527 of the gas;
all species n, ρ, P to ≤ 1e-9: **CONFIRMED** ("33.24 % at T = 12 m" and "8.118 → 7.344").

**Differences from our old lsthermo.py / thermal_relic.py (gp_leftovers equilibrium.py) [D, C].** The old state
functions carry X and X̄ only; s enters `lsrelic` solely through the rate factor K_MB(T)n₊n₋ (Maxwell–Boltzmann), with
no s number, energy, pressure or entropy in the state, no susceptibilities, and caps for one species at a time. Here
the s reservoir is explicit with Bose statistics and its own potential. Consequences at the old onset: same energy
and charge give T/m 8.118 → 7.344, 1.6× the condensed charge and 1.54× the pair inventory; the entropy per net charge
and the X̄ fraction change at every R (table `equilibrium_comparators` in the interface, R = 0–90, with and without s,
with and without V₀):

| R | s | T/m | condensed fraction | X̄ per net charge | s per net charge | P_inv per net charge | σ per net charge | scalar fraction |
|---|---|---|---|---|---|---|---|---|
| 0 | no / yes | 3.685 / 3.348 | 0.847 / 0.875 | 0.185 / 0.134 | – / 0.158 | 0.185 / 0.293 | 1.85 / 2.02 | 0 / 0.272 |
| 1.144 | no / yes | 4.167 / 3.782 | 0.803 / 0.839 | 0.278 / 0.202 | – / 0.235 | 0.278 / 0.437 | 2.67 / 2.92 | 0 / 0.291 |
| 8 | no / yes | 5.629 / 5.099 | 0.636 / 0.703 | 0.737 / 0.536 | – / 0.611 | 0.737 / 1.147 | 6.55 / 7.21 | 0 / 0.317 |
| 20 | no / yes | 6.875 / 6.223 | 0.453 / 0.553 | 1.395 / 1.016 | – / 1.141 | 1.395 / 2.157 | 11.93 / 13.15 | 0 / 0.324 |
| 40.563 | no / yes | 8.118 / 7.344 | 0.233 / 0.374 | 2.356 / 1.718 | – / 1.907 | 2.356 / 3.626 | 19.63 / 21.66 | 0 / 0.328 |
| 48.499 | no / yes | 8.474 / 7.665 | 0.163 / 0.317 | 2.696 / 1.967 | – / 2.177 | 2.696 / 4.144 | 22.32 / 24.64 | 0 / 0.328 |
| 90 | no / yes | 9.852 / 8.907 | 0 (μ = 0.764) / 0.074 | 4.373 / 3.153 | – / 3.454 | 4.373 / 6.607 | 35.05 / 38.70 | 0 / 0.330 |

With the reservoir the free X̄ per net charge falls by 27 % but the pair inventory that must be destroyed rises by
54 %, the entropy per net charge by 10 % and the condensed fraction by 0.03–0.15; without s the gas at R = 90 is
not even capped (μ = 0.764 m). Codex's fixed-T MB susceptibility 1 + K(n₊+n₋) compared with the exact Bose
∂P_inv/∂n₋|_{T,n_Q} (table `pair_inventory_susceptibility`): identical at T ≲ 0.3 m; the MB form overestimates it by
11 % at T = m (1.835 vs 2.054) and by 18–19 % for T ≥ 3 m (2.744 vs 3.382 in Codex's with-s state at T = 7.34 m,
μ = 0.861 m), because the Bose enhancement of n₊n₋ in the denominator of K exceeds that of n_s: K_Bose/K_MB = 0.98,
0.855, 0.81–0.83 at T/m = 0.3, 1, ≥ 3 (tabulated for T/m = 0.03–100).

**In-medium masses [E].** The state functions take any masses (`species(t, r, u)`); the tadpoles ∂P/∂m² follow from
the same series. The high-T one-loop estimate (Dolan–Jackiw, X doublet + φ in the loop) gives v_T² = v² − T²(λ_φ/4 +
g²/6)/λ_φ, m_X²(T) = g²v_T² + T²(λ_X/3 + g²/12), **T_c = 1.672 GeV = 18.09 m_X**: at T = 12 m, v_T/v = 0.75 and
m_X(T)/m_X = 0.85; at 7.34 m, 0.91 and 0.93. Above T_c there is no broken phase and no object. The exact one-loop
functions are kin_rates' job; all tables here use vacuum masses, as Codex's control does, and the mean-field work is
handled by the standard statement that at a self-consistent mean field the Gibbs–Duhem identities hold with the
explicit (T, μ) derivatives while ρ includes V_eff(v_T, T).

## 2. Object sector (`code/objects.py`, `objects_refine.py`, `results/objects_EQ.json`)

**Profile solver (own, third formulation).** Unknowns (F, F′, S, S′) on [10⁻⁴, L]; inner condition the regular power
series through O(ρ²); outer conditions carrying the Yukawa 1/ρ: F′ + (κ + 1/L)F = 0, κ = √(1−Ω²), and
S′ + (√6 + 1/L)(S−1) = 0; L enlarged until the charge outside L/2 is < 1e-9; collocation (solve_bvp, tol 1e-9);
continuation in Ω with adaptive steps; Simpson quadrature on 60001 points. Checks: Derrick/virial identity
E_grad + 3(E_pot − E_kin) = 0 to ≤ 2.9e-10 over the table; dE/dQ = Ω m_X (exact for stationary Q-balls) by finite
differences at the reference (rel. 3.1e-8 at δΩ = 5e-5) and by cubic splines E(Ω), Q(Ω) along the table (|dE/dQ − Ω| ≤
4.6e-4 m_X on both branches, 130 points, excluding |Ω − Ω_min| < 0.012 where dQ/dΩ → 0).

**Reference Q = 200000 [C]:** Ω = **0.8611951308217** (archived 0.8611951308217225: 1.8e-14), E/(Q m_X) =
**0.923316135372** (archived 0.9233161353726057: 2e-13), E = 17064.5375 GeV, **E − Qω = 1148.10753 GeV** (Codex
1148.1075313764: 2e-12), F₀ = 2.61992, S₀ = 0.006012, r₅₀ = 3.786/m_X, r₉₀ = 5.813/m_X, r₉₉ = 8.127/m_X, energy
split kinetic 46.6 % / gradient 10.1 % / potential 43.3 %. All CONFIRMED.

**The branch [C]** (139 rows, Ω = 0.7601–0.998, Q = 5.20e4–3.51e7; coarse run plus adaptive refinement; the
thin-wall continuation stopped at Ω = 0.7600 when solve_bvp exceeded its mesh-node cap at r₉₀ = 26/m_X, L = 96/m_X):
- stable branch (dQ/dΩ < 0) for Ω ≤ 0.9675, Q from Q_min upward; thin-wall limit ω₀ = min_|X| √(U/|X|²) =
  √(2√(λ_X V₀)) = **0.7401 m_X** [D, from U_eff of section 4, attained at |X| = (V₀/λ_X)^{1/4} = 1.53 GeV > x_c],
  Q → ∞ with E/Q → ω₀; the table approaches it as a surface term, (E/Q − ω₀) Q^{1/3} = 10.7, 10.4, 10.0 m_X at
  Q = 2e5, 1.4e6, 3.5e7 (E/Q = 0.7706 m at Q = 3.5e7);
- **Q_min = 5.196e4 at Ω = 0.9675** (E/Q = 1.0146 m; grid step 0.0025): no stationary object exists below Q_min;
- E/Q = m_X at Ω = 0.9312, **Q = 6.56e4**: objects with 5.2e4 < Q < 6.6e4 are classically stationary but heavier than
  free quanta;
- the thick-wall branch Ω ≥ 0.970 (dQ/dΩ > 0, 13 rows to Ω = 0.998, Q = 1.32e5) is the critical-nucleus branch: its
  E − Qω is the nucleation barrier (section 4).

| Ω = ω/m_X | Q | E/(Q m_X) | E − Qω (GeV) | r₉₀ m_X | branch |
|---|---|---|---|---|---|
| 0.7601 | 3.511e7 | 0.7706 | 33801 | 26.15 | stable |
| 0.7803 | 4.542e6 | 0.8018 | 8997 | 13.43 | stable |
| 0.8000 | 1.433e6 | 0.8323 | 4283 | 9.45 | stable |
| 0.8200 | 6.305e5 | 0.8631 | 2514 | 7.52 | stable |
| 0.8400 | 3.377e5 | 0.8931 | 1658 | 6.46 | stable |
| 0.8612 | 2.000e5 | 0.9233 | 1148 | 5.81 | stable (reference) |
| 0.8800 | 1.368e5 | 0.9481 | 861 | 5.50 | stable |
| 0.9000 | 9.786e4 | 0.9715 | 647 | 5.37 | stable |
| 0.9200 | 7.444e4 | 0.9911 | 489 | 5.43 | stable |
| 0.9400 | 6.025e4 | 1.0058 | 366 | 5.74 | stable |
| 0.9600 | 5.285e4 | 1.0138 | 263 | 6.49 | stable |
| 0.9675 | 5.196e4 | 1.0146 | 226 | 7.00 | Q_min |
| 0.9800 | 5.444e4 | 1.0129 | 166 | 8.53 | thick-wall |
| 0.9900 | 6.617e4 | 1.0082 | 111 | 11.65 | thick-wall |
| 0.9980 | 1.317e5 | 1.0019 | 48 | 25.39 | thick-wall |

**Work term.** dρ_L = ω d(n_L Q) + (E − Qω) dn_L; integrated along the table against the exact ΔE for population
changes (n_L, Q) → (1.3 n_L, Q′) over 20 segments of the stable branch: identity to ≤ 2.7e-5 (path quadrature); the
**linear comparator** ωQn_L with ω = ω(2e5) errs by −14 % (Q 5.2e4 → 5.5e4), −10 % (7.2e4 → 8.8e4), −4 % (1.4e5 →
1.9e5), +0.2 % (2.8e5 → 4.2e5), +5 % (6.9e5 → 1.3e6), +10 % (5.5e6 → 8.5e6), +12 % (2.7e7 → 3.2e7). Per object
at Q_ref, (E − Qω) = 1148 GeV = 12,424 m_X = 6.2 % of E: forming or dissolving a Q = 2e5 object moves 6.2 % more energy
than the linear booking. Both bookings are callable (`expansion.object_work`); the consistent linear comparator is
E_c = μ_c = ω, as ls_threshold's `run_ke_consistent` used.

## 3. Exchange kernels at the hidden temperature (`code/kernels.py`)

**Conversion (bb211d1's channel) and its reverse.** X̄(E₋) + core(Q) → X(E₊) + core(Q−2), E₊ = E₋ + E(Q) − E(Q−2) =
E₋ + 2ω m_X (+O(1/Q)). The hard X: for a cold incoming X̄, E_out = 2.7224 m_X = 251.57 MeV, kinetic energy 1.7224 m_X,
p_out = 2.532 m_X [D; Codex's 251.573 MeV CONFIRMED]. Forward σv(p₋) from our dense table (741 momenta,
p₋/m_X = 0.001–0.708, cold Q = 2e5 profile, frozen-profile weak-wave rate): cold limit 9.683 GeV⁻² = 1.13038e-16
cm³/s; it rises to ~10³ GeV⁻² above p = 0.1 (p-waves and the ℓ = 5 resonance at p = 0.508). Reverse by reciprocity
[Q bb211d1]: (σv)_r(p₊) = (σv)_f(p₋) p₋E₋/(p₊E₊), threshold p₊ = 2.532 m_X. Rates per volume with Bose factors,
R_f = n_L ∫ f₋ σv_f (1+f₊), R_r = n_L ∫ f₊ σv_r (1+f₋); **they balance exactly when μ₊ − μ₋ = E(Q) − E(Q−2)**, i.e.
μ₊ = ω at μ₋ = −μ₊ (checked to ≤ 1e-9, quadrature-limited, four states): the conversion channel has zero net rate in
the capped equilibrium state, and it is the object's "chemical potential" ω that sets the balance, not pair
annihilation. Thermal averages at μ = 0 (T/m = 0.03, 0.1, 0.3, 1, 3, 7.34, 12): ⟨σv⟩_f = 1508, 1247, 838, 726, 770,
830, 860 GeV⁻² with the 'hold' extrapolation; the fraction of the X̄ weight inside the table is 0.998, 0.74, 0.23,
0.033, 0.006, 0.002, 0.001 — **above T ≈ 0.1 m_X the capture rate is an extrapolation**, bracketed by 'hold' and
'geometric' (πr₉₀²v, upper); the cold profile itself is not the object at T ~ 10 m. ⟨E_in⟩, ⟨E_out⟩ are tabulated.

**Surface exchange of X (evaporation/absorption), Langmuir model with detailed balance** [Laine–Shaposhnikov
thermodynamics, R]: dQ/dt = η(A/4)[∫v f₊(T, μ₊) − ∫v f_eq(T, ω)], A = 4πr₉₀² = 4.97e4 GeV⁻²; zero at μ₊ = ω
(checked); growth for μ₊ > ω, evaporation below. η ∈ [0, 1] is a labelled bracket: at the cold-core linear level the
two-body X absorption is closed (no hidden state lighter than m_X − ω), so η → 0 for a cold core, while a thermally
excited surface absorbs with η = O(1). With η = 1 the exchange rate per object is ~5e2 GeV at T = 7.3 m (≈ 10¹⁷ H):
the caps would be enforced instantly, as the old proxies assume; with η → 0 only the conversion channel acts. The
energy bookkeeping offers "cold object" (mismatch returned to the gas) explicitly. Object/anti-object encounters:
geometric bracket σ = π(2r₉₀)² only.

## 4. Formation alternatives (`code/formation.py`)

(a) **gas_only**: no objects; charge beyond the Bose capacity condenses homogeneously at μ = m_X (energy m_X per
charge; case-H comparator). Supplies: the cap u = 1, e = 1.
(b) **instant_condensation** (the old proxies' rule): whenever μ₊ would exceed ω(Q_ref) the excess is condensed into
Q_ref objects at once, and dissolves at once when undersaturated; anti-objects likewise at μ₋ = ω. Variants: e_cap = ω
(consistent linear), e_cap = E/Q (benchmark booking), finite_Q (adds (E − Qω)dn_L). Supplies: no formation time, no Q
distribution. Called on the two-potential state it reproduces Codex's with-s state (T/m 7.3441, μ₋ = −0.8623).
(c) **fragmentation / nucleation bracket** scaled to this model [D for the kernel's scales, E for the brackets]:
- With φ integrated out adiabatically (φ² = v² − 2g²|X|²/λ_φ ≥ 0), **U_eff(|X|) = m_X²|X|² − λ_eff|X|⁴** for |X| ≤ x_c,
  **λ_eff = g⁴/λ_φ − λ_X = 2.3333e-3 > 0**, and U = λ_X|X|⁴ + V₀ beyond (continuous at x_c to 1e-9). The effective
  potential is flatter than quadratic: that is why the two-field Q-balls exist, and it makes the homogeneous
  condensate X = X₀e^{iωt} linearly unstable (Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki, R; own
  re-derivation): for modes e^{αt−ik·x},
      (α² + k²)(α² + k² + U″ − ω²) + 4ω²α² = 0,   band 0 < k² < k_max² = ω² − U″(φ₀) = 4λ_eff|X₀|²,
  ω² = m² − 2λ_eff|X₀|², U″ the second radial derivative at φ₀ = √2|X₀|.
- **At the onset amplitude |X₀| = x_c** (the release before any thermalisation): ω = 0.5477 m (= √0.3 m, continuous with
  the quartic-phase frequency), k_max = 1.183 m_X, α_max = 0.639 m_X = 0.0591 GeV at k = 0.540 m_X, i.e. **α_max/H =
  2.2e13**: fragmentation is instantaneous on Hubble times, t_frag = N_e/α_max = 508 GeV⁻¹ for N_e = 30. Lump charge
  Q_lump = shape × n_c (2π/k_max)³ = **2.46e4 × shape** (shape ∈ [0.1, 1]; Hiramatsu–Kawasaki–Takahashi's 3D peak is
  60 % above the linear estimate, R): **below Q_min = 5.2e4 for every shape** — the direct fragments of the release
  are not stationary objects and must merge or dissolve [flagged by the rule].
- **From a gas condensate** (the gas saturates at μ₊ → m first; see the nucleation barrier): for a condensed fraction
  f_c of n_Q, k_max² = 4λ_eff n_cond/(2ω) and Q_lump = shape × n_cond(2π/k_max)³ ∝ n_cond^{−1/2}: f_c = 0.7 gives
  |X₀| = 0.270 GeV, k_max = 0.283 m, α_max = 0.0204 m = 1.9e-3 GeV (7e11 H), Q_lump = 1.84e5 × shape (5.5e4 at
  shape 0.3, at the edge of Q_min; E_lump = 5167 GeV, ω = 0.951 m, work term 306 GeV per lump). Tabulated for
  f_c = 0.03–1 and shape 0.1/0.3/1.
- **Nucleation from an undersaturated gas** is negligible: F_crit(μ) = (E − Qω) on the thick-wall branch is 214 GeV
  at μ = 0.970 m and still 47.8 GeV at μ = 0.998 m, so F_crit/T ≥ 43 at T = 12 m and ≥ 1700 at T = 0.3 m.
- Constraints supplied with the rule: Q ≥ Q_min, T < T_c (18 m_X high-T estimate), E/Q < m_X only for Q > 6.6e4.
- Anti-objects: for an elliptic release (ε = 0.146) the lattice literature reports Q-balls of both signs [R,
  Kasuya–Kawasaki 2000]; f_anti ∈ [0, 1] of the symmetric number is a bracket the solver must carry and report.
- Brackets: shape [0.1, 1], N_e [20, 40], f_cap [0.5, 1] (KK: "nearly all charge in Q-balls"), f_anti [0, 1].
Credits and read levels: Kusenko–Shaposhnikov PLB 418 (1998) 46 [abstract read; linear analysis recalled];
Enqvist–McDonald PLB 425 (1998) 309 [abstract]; Kasuya–Kawasaki PRD 61 (2000) 041301 and PRD 62 (2000) 023512
[abstracts]; Hiramatsu–Kawasaki–Takahashi JCAP 06 (2010) 008 [abstract]; Laine–Shaposhnikov NPB 532 (1998) 376
[abstract]; Griest–Kolb 1989, Postma 2002 [recalled]. No lattice run; nothing here is a measured formation rate.

## 5. Expansion and the energy equation (`code/expansion.py`)

Visible radiation with g* = 106.75 (task statement), H = √(π²g*/90)T_v²/M_Pl, a ∝ 1/T_v, hidden charge comoving
(n_Q = Y_Q s_vis). At the case-L coordinate a_ref = 7.0467e12, T_vis = 43.72 GeV with g_s = g* = 106.75 (45.67 GeV
with the coarse g* table that the old proxies use and kin_source quotes), H = 2.685e-15 GeV (2.745e-15 with the
table: 2.2 %); the fragmentation rule's α_max/H uses this consistent H. Energy equation

    dρ_gas/dln a = −3(ρ_gas + P_gas) + [S_parent + S_med − W_obj − W_mf]/H,   S_parent = Γ_A m_A n_A(a),
    dn_A/dln a = −3n_A − Γ_A n_A/H,   W_obj = ω d(n_LQ)/dt + (E − Qω) dn_L/dt,

total d[(ρ_gas + ρ_A + ρ_L)a³] = −P_gas d(a³). Parents for the four sources (band/broad, e = 1e-7/1e-6) from
DECAY_SOURCE.json [Q]. Checks [C]: isentropic expansion of the equilibrium gas integrated with the susceptibility
Jacobian conserves s a³ and n_Q a³ to 2e-12 and 1e-11 over three e-folds; the parent ledger reproduces N_A and
R = 40.56262 / 4.053801 to 2e-14; the object-work identity to 2.7e-5.

## 6. The two old proxies at the heavy point, R = 1–90 (`code/old_proxies.py`, reproduction)

|ε| = 0.14571744, yield 1.00 Y_Q, coarse g* table, R grid 1, 1.14445, 2, 3, 5, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40.5626,
41.707, 48.499, 49.643, 60, 75, 90 (the mediator-only, vector + mediator and broad + mediator points included).

| R | isentropic proxy f | isentropic + s route | ke (μ_c = E_c = ω) | ke + s route | ke (μ_c = E_c = 0.923 m) | T_h/m at onset (proxy) | σ/n_Q at onset (ke) |
|---|---|---|---|---|---|---|---|
| 1 | 0 (erased) | 0 (erased) | 0 (erased, 213 fb) | 0 (erased, 213 fb) | – | 4.109 | 2.57 |
| 1.1444 | 0 (erased) | 0 (erased) | 0 (erased, 212 fb) | 0 (erased, 212 fb) | – | 4.161 | 2.67 |
| 2 | 0 (erased) | 0 (erased) | 0 (erased, 211 fb) | 0 (erased, 211 fb) | – | 4.437 | 3.23 |
| 3 | 0 (erased) | 0 (erased) | 0 (erased, 208 fb) | 0 (erased, 209 fb) | 0 (erased, 301 fb) | 4.705 | 3.84 |
| 5 | 0 (erased) | 0 (erased) | 0 (erased, 200 fb) | 0 (erased, 207 fb) | 3.7e-4 | 5.134 | 4.99 |
| 7 | 0 (erased) | 0 (erased) | 1.6e-26 | 0 (erased, 205 fb) | 0.118 | 5.477 | 6.05 |
| 8 | 0 (erased) | 0 (erased) | 0.00568 | 0 (erased, 204 fb) | 0.187 | 5.627 | 6.55 |
| 9 | 0 (erased) | 0 (erased) | 0.050 | 0 (erased, 203 fb) | – | 5.766 | 7.05 |
| 10 | 0 (erased) | 0 (erased) | 0.104 | 0 (erased, 201 fb) | 0.325 | 5.895 | 7.53 |
| 12 | 0 (erased) | 0 (erased) | 0.216 | 0 (erased, 196 fb) | – | 6.131 | 8.47 |
| 15 | 0 (erased) | 0 (erased) | 0.383 | 4.6e-14 | 0.662 | 6.440 | 9.81 |
| 20 | 0 (erased) | 0 (erased) | 0.653 | 0.0931 | – | 6.874 | 11.93 |
| 25 | 0 (erased) | 0 (erased) | 0.910 | 0.241 | – | 7.239 | 13.93 |
| 30 | 0 (erased) | 0 (erased) | 1.16 | 0.389 | – | 7.555 | 15.83 |
| 40.5626 | 0 (erased) | 0 (erased) | 1.65 | 0.691 | – | 8.117 | 19.63 |
| 41.707 | 0 (erased) | 0 (erased) | 1.70 | 0.722 | – | 8.172 | 20.02 |
| 48.499 | 0 (erased) | 0 (erased) | 2.01 | 0.908 | – | 8.474 | 22.32 |
| 49.643 | 0 (erased) | 0 (erased) | 2.06 | 0.939 | – | 8.522 | 22.70 |
| 60 | 0 (erased) | 0 (erased) | 2.49 | 1.21 | – | 8.922 | 26.05 |
| 75 | 4.0e-6 | 1.8e-6 | 3.10 | 1.59 | – | 9.421 | 30.66 |
| 90 | 0.428 | 0.200 | 3.67 | 1.95 | – | 9.852 | 35.05 |

"erased" = p_final at lsrelic's floor (f < 1e-30), "fb" = root fallbacks in the ke step solver (those rows mean only
"erased"); f = 2 m_X Y_Q × pairs per net charge / 4.4e-10 GeV. Proxy: N = 800 ODE steps, s route via the MB K(T)
factor; ke: N = 500, antilumps on, E_c = μ_c. All 92 rows finished, no failures (`results/old_proxies.json`).

Crossings f = 0.01 [C, reproduction]: **isentropic proxy R_L = 82.03–82.19** (bisection, N = 400) and **83.13–83.28
with the s route** (the s route lowers f at R = 75–90 by 2×, which moves the crossing by 1.3); the coarse-grid log
interpolation between R = 75 and 90 would give 84.8 / 85.9 and is not to be used. **ke (μ_c = E_c = ω) R_L = 8.25**
(grid [8, 9]); **ke + s route 19.55** ([15, 20]); ke at μ_c = E_c = 0.923 m: 6.06 ([5, 7]; ls_critic's 7.0 came from
the wider [5, 8.9] bracket, the crossing lies in [5, 7] either way). These reproduce ls_critic's direct values at the
same point (proxy 81.9–82.2; ke 8.2; ke + s between 15 and 20). Both isentropic variants erase everything up to
R = 60 and rise from 4e-6 (R = 75) to 0.43 / 0.20 (R = 90); the ke curves rise smoothly from R ≈ 7 (15 with s) to
f = 1.65 (0.69 with s) at R = 40.5626 and 2.01 (0.91) at 48.499, i.e. the central vector source fails the 1 %
allowance by 69–165× in ke and passes by > 10¹⁸⁰ in the isentropic proxy: that 10× gap in R_L is what the solver
must replace. The ke runs with ~200 root fallbacks
sit at the p floor and mean only "erased"; the cumulative energy-work residual of ke is 1–2 % at N = 500 (as in
ls_threshold). These are equilibrium proxies, not kinetics; they are the comparators the solver must beat.

## 7. The interface (`results/thermo_interface.json`, 1.28 MB)

`meta` (conventions, import recipe, function list, fixed point); `thermodynamics` (series tables for X and s on
87 × 31 and 87 × 19 (t, u) grids with n, P, ρ, s and the four derivatives; K_Bose/K_MB; Codex control; equilibrium
comparators on the R grid; susceptibility examples; thermal-mass estimate; differences from the old code);
`expansion` (background, parents, energy equation, checks); `objects` (E(Q) table with Ω, Q, E, E/Q, E − Qω, radii,
branch, dE/dQ); `kernels` (forward/reverse σv tables, thermal averages at μ = 0 and in the cap state for two
extrapolations, surface-exchange parameters and examples, collision bracket); `formation` (the three rules'
callables and parameters, constraints, nucleation barrier, literature); `old_proxies` (curves and crossings);
`checks`. The solver should import `code/thermo.py` etc. directly; the tables are for use without the code.

## 8. Checks summary (`results/checks.json`, 23/23 passed, residual and threshold stored separately)

| check | residual | threshold |
|---|---|---|
| Bessel series vs direct quadrature, 7 quantities × 28 states | 2.0e-15 | 1e-9 |
| ∂P/∂μ = n, ∂P/∂T = s, analytic ∂n/∂T vs finite differences | 3.7e-11 | 1e-8 |
| Codex reservoir rows (5 T, two ratios) | 2.8e-14 | 1e-9 |
| Codex matched-energy states (T/m 8.118 → 7.344, condensed, species n, ρ, P) | 2.2e-13 | 1e-8 |
| two-potential solver at the equilibrium inventory = equilibrium solver | 7.6e-15 | 1e-9 |
| three-potential Newton round trip | 5.0e-16 | 1e-9 |
| ∂P_inv/∂n₋ Bose exact vs MB 1 + K(n₊+n₋) (reported, no threshold) | 0.19 | – |
| isentropic expansion: s a³ and n_Q a³ drift over 3 e-folds | 9.6e-12 | 1e-9 |
| parent ledgers e = 1e-7, 1e-6 (daughters, R) vs closed forms | 2.1e-14 | 1e-10 |
| Q = 2e5: ω/m, E/(Qm), E − Qω vs archived / Codex | 1.8e-14, 1.4e-13, 1.8e-12 | 1e-10, 1e-10, 1e-8 |
| dE/dQ = Ω at the reference (finite difference) | 3.1e-8 | 1e-6 |
| max \|virial\| over the table | 2.8e-10 | 1e-9 |
| dE/dQ = Ω along the table (splines, both branches) | 4.6e-4 | 2e-3 |
| object-number work identity (path integral vs exact ΔE) | 2.7e-5 | 1e-3 |
| conversion table vs Codex cold / daughter σv | 5.6e-6, 1.0e-6 | 1e-4, 1e-3 |
| conversion forward = reverse at μ₊ − μ₋ = 2ω (4 states) | 3.9e-9 | 1e-7 |
| surface exchange zero at μ₊ = ω | 0 | 1e-12 |
| U_eff continuity at x_c; k(α_max) inside the band | 1.8e-9; pass | 1e-6 |
| reproduction: lsthermo.py series vs this series (labelled, not independent) | 4.4e-16 | – |

## 9. Limits

- Vacuum masses throughout; the thermal phase history (T_c ≈ 18 m at one loop, 42 m with X only [Q ls_threshold])
  and the in-medium s width are open and belong to kin_rates; above T_c no objects exist, which the formation rules
  flag but the state functions do not enforce.
- The conversion kernel is the cold frozen-profile weak-wave rate; above p₋ = 0.708 m_X it is an extrapolation, and a
  hot, newly formed or excited core is not that profile. The surface-exchange η is a bracket, not a computed number.
- The fragmentation kernel is a linear-instability scale with literature O(1) shape factors; the ellipticity
  (ε = 0.146) dependence and the ± lump split are brackets, not results; no lattice run; mergers of sub-Q_min
  fragments are not modelled.
- E(Q) is the zero-temperature ground-state branch; object excitation, finite-T profiles and object–object
  interactions are absent.
- The old-proxy curves inherit lsrelic's limitations (instant caps, MB K(T), 1–2 % energy-work residual, p floor).

## 10. Reproduce

    cd research/reality_framework/grut_formation_review_20261025/tracks/kin_thermo/code
    python3 -B objects.py                                           # ~1 min: coarse E(Q) table, reference profile
    python3 -B objects_refine.py                                    # ~25 min: adaptive continuation, 139 rows
    OMP_NUM_THREADS=1 python3 -B old_proxies.py ../results/old_proxies.json [R1,R2,...]
                                                                    # ~1.5 h on 4 cores (reproduction; resumable:
                                                                    #  finished rows are kept; the optional R list
                                                                    #  runs only those R, for short foreground batches)
    python3 -B checks.py                                            # ~10 s: results/checks.json (23 checks)
    python3 -B build_interface.py                                   # ~4 min: results/thermo_interface.json

Order matters: `objects*.py` before `checks.py`, `formation.py` and `build_interface.py` (they read
`results/objects_EQ.json`); `old_proxies.py` before `build_interface.py` (otherwise the `old_proxies` section is
marked "not yet available"). Dependencies: numpy 2.5, scipy 1.18; the conversion tables are read from
`tracks/conversion/results/`, the old proxies import `tracks/ls_threshold/code/lsrelic.py` unchanged. Module
self-tests: `python3 -B thermo.py`, `kernels.py`, `formation.py`, `expansion.py`. The reference profile is cached
in the scratchpad (`profile_Q200000.npz`) for `objects_refine.py`; `objects.py` regenerates it.

**Status (finalising pass, 26 October, third session).** The only computation added was the completion of the ten
s-route isentropic proxy rows R ≥ 20 that the interrupted run had left out (two foreground batches, same code and
settings, 47–51 s per row), after which `build_interface.py` was re-run so the interface carries the refined 139-row
E(Q) table (Q_min, thin-wall and E/Q = m fields), the complete old-proxy curves and crossings, the 23 checks and the
consistent H at onset (α_max/H = 2.2e13; the first build had used T_vis = 45.67 GeV in H). Spot checks repeated
independently of the track code (direct Bose quadrature): scalar energy fraction 0.33236563 and n_s/n₋ = 0.98106981
at T = 12 m, E − Qω = 1148.1075 GeV, E_out = 251.573 MeV, ω₀ = 0.74008 m, T_vis(a_ref) = 43.717 GeV. Nothing else
was recomputed; nothing was committed.
