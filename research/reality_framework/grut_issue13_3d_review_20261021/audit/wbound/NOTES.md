# Audit track `wbound`: what can be proven about the hidden-sector w after wall–radiation parity

21 October 2026. Claude audit of Codex `a69bd86` item (4), against Claude `c7b8480`
(REPORT.md §0 item 4, §5(b), §6 replies; verification/v_cosmo/VERIFY.md claim 7; code/history.py line 107).
Labels: **[T]** theorem (proof given here), **[C]** computed (numerics in this directory), **[E]** estimate,
**[A]** argued, not proven, **[P]** preparation input. Units in code: v = 1, λ_φ = g_X² = 0.01, λ_X = 0.001,
h_D = 0.01, λ_D = 0.001, m_D0 = 0.1 (the benchmark couplings, unchanged action). Nothing here excludes, validates
or identifies Candidate 1.

## Verdicts

| # | Claim | Verdict | Basis |
|---|---|---|---|
| V1 | (a) 3p − ρ = (φ̇² + 2\|Ẋ\|²) − ((∇φ)² + 2\|∇X\|²) − 4U (+ Ḋ² − (∇D)²) | **CONFIRMED** | sympy from T_μν: residual 0; equals the trace T^μ_μ (algebra.json) |
| V2 | About +v, δU′ − 4U = −λvδ²(φ+v) ≤ 0 for φ ≥ −v (φ only) | **CONFIRMED** | sympy, factored form |
| V3 | X makes it worse? | **CORRECTED/sharpened**: pointwise yes for the +v origin (extra term −2g²vφ\|X\|², ≤ 0 only for φ ≥ 0); for the exact criterion no — with origin 0, X drops out entirely | sympy; origin independence of the average |
| V4 | Codex: instantaneous w ≤ 1/3 is not a theorem even at small mean density | **CONFIRMED** | homogeneous kick at the minimum (w = 1); 1D slabs and 3D lumps start at w = 0.996–1.000 with mean density down to 1.9×10⁻⁸ v⁴ (3D) and 7.5×10⁻⁷ v⁴ (1D) [C]; their time averages fall below 1/3 within t ≈ 3–7 (linear) or about the light-crossing time (over-barrier) |
| V5 | Codex: the cosmological argument needs a justified coarse-grained state family | **CONFIRMED, and supplied** | the dynamics does the coarse-graining exactly (FRW virial budget, §3); the remaining state condition is one inequality, ⟨φ²⟩ ≤ v² + (m_D0²/λv²)⟨D²⟩, Hubble-weighted |
| V6 | Claude: after parity the energy-weighted w̄ of the hidden sector stays ≤ 1/3 | **CORRECTED** | true as a Hubble-time integral (not instantaneously, not "every component") if and only if the inequality in V5 holds; proven for single-sign and \|φ\| ≤ v states; computed for walls, waves, quanta and lumps; for the post-parity network it is [A], not proven |
| V7 | Claude: only coherent oscillations spanning both vacua have w > 1/3, up to 0.44 | **CORRECTED** | homogeneous maximum 0.44151 confirmed [C]; but "only" is too strong: any state with space–time mean φ² > v² qualifies, including transient debris of relativistic wall collisions (max \|φ\| 2.1–2.8 v at γ = 2.3–22) and over-barrier kinetic lumps [C] |
| V8 | Claude: those states need an energy density 30 orders above parity, so they cannot occur after parity | **REFUTED as an argument** | compares a local threshold with the mean density; wall cores themselves sit at λv⁴/2 ≈ (2.7–4.0)×10³⁰ ρ_h(parity); local over-barrier states are produced by the walls. Replaced by §5 (identity + wall deficit + dilatation bound) [A] |
| V9 | VERIFY: w > 1/3 needs ρ ≳ 0.0225 v⁴ ≈ 2×10¹⁰ GeV⁴ | **CORRECTED** | threshold ε* = 0.58768 λv⁴ = 0.0058768 v⁴ = 5.9×10⁹ GeV⁴ (amplitude 1.5916 v); 0.0225 v⁴ is the amplitude-2v energy, next to the w̄ maximum (ε = 2.127 λv⁴) [C] |
| V10 | VERIFY: Φ_{+v} ≤ 0 wherever φ > 0, so any motion within one vacuum has w ≤ 1/3 | **CONFIRMED** (time average) | theorem T1 below; needs the O(1/T) boundary term, not instantaneous |
| V11 | NG walls w = ⟨v²⟩ − 2/3 ≤ 1/3, → 1/3 only as v → 1 | **CONFIRMED** | energy-weighted ⟨u²⟩_E; field theory ∫(3p−ρ)dz = −3σ/γ < 0 (sympy, static); 1D lattice at u = 0–0.999 within 0.0012 of NG [C] |
| V12 | Q-balls, NR quanta w ≈ 0; GW and relativistic quanta w = 1/3 | **CONFIRMED**, with wording | stationary localized lumps: ∫⟨T^ij⟩ = 0 (von Laue), moving: w = ⟨u²⟩_E/3; massive quanta w < 1/3 strictly; only massless w = 1/3 |
| V13 | ΔN_eff ≥ 7.4 | **CORRECTED wording**: ρ_h/ρ_1ν ≥ 7.447 (N_eff,std = 3.044) is a radiation-equivalent energy ratio, conditional on V6; not a Planck N_eff likelihood; "least constraining case" is unsupported | algebra.json; §7 |
| V14 | Codex DERIVATION §3: p − ρ/3 = (2K − 4V_τ)/(3a⁴) | **CORRECTED** (typo-level) | omits the gradient term: p − ρ/3 = (2K − 2G − 4V_τ)/(3a⁴). expanding.cpp line 174 uses the correct p = (K − G/3 − U)a⁻⁴; the homogeneous counterexample (G = 0) is unaffected |
| V15 | Codex's own 3D runs | **CONSISTENT** | in all 12 production runs w ≤ 0.140 at every sample and ρa⁴ grows monotonically (×5.8–6.8 over a = 1→2.2) (codex_eos.json) [C, coarse sampling] |
| V16 | Draft form "states outside F need local densities ≳ λv⁴/4 and are not produced after parity" | **CORRECTED** | the local-density half is a theorem (T3), but such states *are* produced locally after parity (wall cores at λv⁴/2, relativistic collisions). The defensible [A] claim is that their Hubble-averaged energy share stays below what reversing the wall deficit needs (§5) |

## 1. (a) Trace identity

Signature (−,+,+,+), L = −½(∂φ)² − |∂X|² − ½(∂D)² − U. With X = (x₁ + i x₂)/√2 every field is a canonical real
component f_a. Then ρ = Σ(½ḟ_a² + ½(∇f_a)²) + U, p = Σ(½ḟ_a² − (1/6)(∇f_a)²) − U, and

    3p − ρ = T^μ_μ = Σ ḟ_a² − Σ (∇f_a)² − 4U = (φ̇² + 2|Ẋ|² + Ḋ²) − ((∇φ)² + 2|∇X|² + (∇D)²) − 4U.

Verified symbolically from T_μν (algebra.py → results/algebra.json, residual 0).

## 2. (b) Virial identities in flat space

**Identity [T].** Periodic box V (or a lattice: the forward-difference gradient energy is paired with the
(2d+1)-point Laplacian, so the identity is exact there too). For any constant field-space origin c let
G_c(t) = ∫_V (f − c)·ḟ. The equations of motion give dG_c/dt = ∫[ḟ² − (∇f)² − (f − c)·∂U], hence for every T

    (1/T)∫₀ᵀdt ∫_V (3p − ρ) = (1/T)∫₀ᵀdt ∫_V Φ_c + [G_c(T) − G_c(0)]/T,      Φ_c ≡ (f − c)·∂U − 4U.

**No "bounded motion" assumption is needed.** U ≥ (λ/4)(φ² − v²)² + λ_X|X|⁴ + ½m_D0²D² bounds ∫_V|f − c|² by the
conserved energy E, and |G_c| ≤ ‖f − c‖₂ ‖ḟ‖₂ ≤ ‖f − c‖₂ (2E)^½. So the boundary term is O(1/T) for every
finite-energy solution in a periodic box. For an isolated configuration in infinite space the same holds while
G_c grows sublinearly (true for dispersing massive or massless radiation).

**Origin independence [T].** Φ_c − Φ₀ = −c·∂U and ∫_V ∂U = −d/dt ∫_V ḟ, so the long-time average of ∫Φ_c does not
depend on c. Only the pointwise sufficient conditions depend on the origin.

**The two useful origins [T, sympy].**

    Φ_{+v} = −λv(φ − v)²(φ + v) − 2g²vφ|X|² − D²(m_D0² + h_D²vφ)
    Φ₀    =  λv²(φ² − v²) − m_D0²D²          (X drops out: every X term is quartic about the origin)

Equivalently Φ = Σ_k (k − 4)U_k for the expansion of U in homogeneous degrees about the origin. About 0 the only
positive piece is the tachyonic quadratic term −½λv²φ², which is the single source of w̄ > 1/3 in this action.

*Does X make it worse?* For the +v origin, yes pointwise: the φ-only condition φ ≥ −v becomes φ ≥ 0 wherever
X ≠ 0 (the D term only needs φ ≥ −m_D0²/(h_D²v) = −100v). Exactly, Φ_{+v} ≤ 0 at a point iff
λ(φ − v)²(φ + v) + 2g²φ|X|² + D²(m_D0²/v + h_D²φ) ≥ 0; for −v ≤ φ < 0 (D = 0) this needs
|X|² ≤ λ(φ − v)²(φ + v)/(2g²|φ|). So "bounded motion that never enters φ < −v" (the brief's wording) is sufficient
only where X = 0; with X the clean sufficient condition is φ ≥ 0 (theorem T1). For the exact criterion, no: X never enters Φ₀. Through
the dynamics X raises φ's mass² by 2g²|X|² and pulls φ² below v²: quanta of X about the vacuum shift
⟨s⟩ ≈ −(3/2v)⟨s²⟩ − (g²/λv)⟨|X|²⟩, giving λv²⟨φ² − v²⟩ = −m_φ²⟨s²⟩ − 2m_X²⟨|X|²⟩ (the correct massive-quanta
trace). In 0D sampling, adding X at fixed energy *lowers* the maximal w̄ (e = 9 barrier units: 0.325 with X vs
0.4416 φ-only; e = 6: 0.305 vs 0.438) [C].

**Exact criterion [T].** Long-time flat-space average:

    w̄ − 1/3 = [λv²⟨φ² − v²⟩ − m_D0²⟨D²⟩] / (3ρ̄),       ⟨·⟩ = space–time mean.

**Sufficient conditions [T].**
- T1: φ keeps one sign everywhere at all times ⇒ Φ_{±v} ≤ 0 pointwise ⇒ w̄ ≤ 1/3 (X and D included).
- T2: |φ| ≤ v everywhere at all times ⇒ Φ₀ ≤ 0 ⇒ w̄ ≤ 1/3.
- T3: a sign change of φ at (x,t) needs ρ(x,t) ≥ U(0, X, D) ≥ λv⁴/4. So a single-sign start whose local energy
  density stays below λv⁴/4 everywhere obeys T1.
- General cap (Cauchy–Schwarz with (λ/4)(φ² − v²)² ≤ ρ): w̄ ≤ 1/3 + (2/3)(λv⁴/ρ̄)^½. No sustained stiff phase at
  ρ̄ ≫ λv⁴; the cap is weak at low mean density, which is why §5 is needed.

**Which states can have w̄ > 1/3 (exactly):** those whose space–time mean of φ² exceeds v² + (m_D0²/λv²)⟨D²⟩.
Necessary (pointwise): φ takes both signs (else T1) and |φ| > v somewhere (else T2). Catalogue:

| State | w̄ (time-averaged) | Source |
|---|---|---|
| homogeneous orbit inside one well (ε < λv⁴/4) | ∈ (−1, 0], → 0 at small amplitude | quadrature [C] |
| homogeneous over-barrier orbit, λv⁴/4 < ε < ε* | < 1/3 (time near φ = 0 dominates) | [C] |
| homogeneous over-barrier orbit, ε > ε* = 0.58768 λv⁴ (amplitude > 1.5916 v) | > 1/3; maximum **0.44151** at ε = 2.1270 λv⁴ (amplitude 1.9791 v); → 1/3⁺ as ε → ∞ (0.339 at ε = 2500 λv⁴) | quadrature, identity to 10⁻¹² [C] |
| homogeneous φ, X, D (384 random orbits, e = 0.5–100 barrier units) | max 0.4416 (φ-only); 0 one-sign violations; identity to ~10⁻⁷ after the boundary term | zero_d.json [C] |
| NG walls | ⟨u²⟩_E − 2/3 < 1/3 (u < 1) | [T] |
| field walls, steady at speed u | ∫(3p − ρ)dz = −3σ/γ: 1D lattice u = 0, 0.5, 0.9, 0.99, 0.999 → −0.6655, −0.4159, 0.1432, 0.3130, 0.3311 (NG −0.6667, −0.4167, 0.1433, 0.3134, 0.3313) | sympy + sims_1d [C] |
| charged stationary walls with \|φ\| ≤ v | λv²∫⟨φ² − v²⟩dz < 0 ⇒ w < 1/3 (X-rich cores push φ toward 0) | [T] |
| waves/quanta about either vacuum | −m_φ²⟨s²⟩ − 2m_X²⟨\|X\|²⟩ − m_D²⟨D²⟩ (mode w_k = k²/3ω_k² < 1/3) | leading order [T]; 1D random waves 0.046–0.307 [C] |
| stationary localized lumps (Q-balls, oscillons): "Q-matter" | ∫⟨T^ij⟩ = 0 (von Laue) ⇒ w = 0; a gas moving at u_i: w = ⟨u²⟩_E/3 | [T] |
| restored rotating X condensate (φ ≈ 0, \|X\| = 1.06 v) | −0.167 (3D lattice) | sims_3d_charged.json [C] |
| 1D random states spanning both vacua (\|φ\| up to 1.8–3.2 v, with X) | 0.060, 0.297, 0.322 | sims_1d A [C] |
| 1D quench to kink networks (32 and 1270 sign changes at end) | 0.233, 0.322 | sims_1d C [C] |
| 3D quench network with X noise (96³, flat; forms walls that annihilate in the small box) | 0.260 (0.165 for t ≤ 150; 0.283 late); max instantaneous ⟨Φ₀⟩/ρ = +0.013 | sims_3d_net.json [C] |
| hot restored plasma at T ≳ √12 v (quantum, ⟨φ²⟩_T ≈ T²/12) | slightly > 1/3 [E]; needs mean density ~ v⁴, irrelevant after parity | [E] |

**Minimum local energy density.** Any sign change: ρ_loc ≥ λv⁴/4 [T]. A coherent over-barrier region with
w̄ > 1/3: ρ_loc > ε* = 0.588 λv⁴ (homogeneous threshold) [C]; inhomogeneous lumps need more because gradient
energy adds to ρ without adding to ⟨φ² − v²⟩ — e.g. a 3D cos² lump of radius 10 at local ε = 2λv⁴ gives
w̄_T ≈ 0.26 [C, smoke test]. **But local density does not separate F from its complement**: wall cores in the
admissible family already have ρ_loc = λv⁴/2 (static) and more when moving.

## 3. FRW with the adiabatic correction (exact budget)

Closed hidden sector, classical minimally coupled fields, statistically homogeneous (periodic comoving box), any
a(t). With g_c = ⟨(f − c)·ḟ⟩ (spatial mean) the equations of motion give **exactly**

    3p − ρ = ⟨Φ_c⟩ + a⁻³ d(a³ g_c)/dt.

With the continuity equation d(ρa⁴)/dt = −Ha⁴(3p − ρ) and Ha = ȧ:

    Δ(ρa⁴) = −∫ H a⁴ ⟨Φ₀⟩ dt − [ȧ a³ g₀] + ∫ ä a³ g₀ dt.

The instantaneous kinetic-minus-gradient imbalance (Codex's counterexample) appears only in the total-derivative
terms. Bound: |g₀| ≤ (2ρ⟨f²⟩)^½ and H² = ρ_tot/(3M_P²), so |ȧa³g₀|/(ρa⁴) ≤ (2⟨f²⟩ρ_tot/(3ρM_P²))^½ ≈
(2/√3)(v/M_P) = **4.7×10⁻¹⁶** at parity (ρ_tot = 2ρ_h, ⟨f²⟩ ≈ v²); the ä term is of the same order per e-fold
(|ä| ≤ aH² for w_tot ≤ 1/3). Over the ≈ 10 e-folds from 6 keV to recombination: < 10⁻¹⁴ (cosmo_numbers.json) [T+E].
For coherent oscillation about a vacuum the same terms are O(H/m_φ) ≈ 10⁻³¹.

**Consequence [T].** ρ_h a⁴(rec) ≥ ρ_h a⁴(par) (to < 10⁻¹⁴) **if and only if**

    ∫_par^rec d ln a · a⁴ ⟨ λv²(φ² − v²) − m_D0² D² ⟩ ≤ 0.

This is the justified state family Codex asked for: no oscillation averaging needs to be assumed; the
Hubble-weighted integral is exactly what the background equation uses.

Lattice check (sims_3d.py `frw`, 64³, dx = 1.5, prescribed a ∝ t^½) [C]:

| state | H₀/m_φ | a range | ρa⁴ end/start | budget closure | boundary term / ρa⁴ | max instantaneous ⟨Φ₀⟩/ρ |
|---|---|---|---|---|---|---|
| quench network (both signs) + X noise | 0.1 | 1 → 2 | 2.658 | 9.3×10⁻⁴ | 4.8×10⁻² | −0.45 |
| same | 0.02 | 1 → 1.41 | 1.283 | 9.4×10⁻⁴ | 2.0×10⁻³ | **+0.12** |
| one-sign waves about +v | 0.1 | 1 → 2 | 1.082 | 9.6×10⁻⁴ | 7.1×10⁻³ | **+0.18** |
| one-sign waves about +v | 0.02 | 1 → 1.41 | 1.028 | 9.6×10⁻⁴ | 7.2×10⁻⁴ | +0.057 |

The instantaneous spatial mean of Φ₀ can be positive (oscillation phase); only the Hubble-weighted integral enters,
and ρa⁴ grew in every run. The continuity equation itself closes to 1.4–1.8×10⁻³ (KDK, coarse sampling).

## 4. (c) The counterexample family and its persistence

**Instantaneous violations are real [C].** Homogeneous kick at the minimum: w = 1 at t = 0, period average → 0
(quadrature: w̄ = −0.0009 at ε = 0.0025 λv⁴). Sparse kinetic lumps (φ = v, π = A·cos² bump, seed noise confined to
the lumps):

| run | geometry | mean ρ (v⁴) | local ε/λv⁴ | w(0) | long-time w̄ | w̄_T last above 1/3 at |
|---|---|---|---|---|---|---|
| lin_A0.02_R5 | 1D slabs, filling 0.005 | 7.5×10⁻⁷ | 0.02 | 0.996 | 0.214 | t = 5 |
| lin_A0.02_R50 | 1D slabs, filling 0.05 | 7.5×10⁻⁶ | 0.02 | 0.997 | 0.018 | t = 7 |
| over_A0.2_R20 | 1D slabs | 3.0×10⁻⁴ | 2.0 | 1.000 | −0.196 | t = 5 |
| over_A0.2_R100 | 1D slabs | 1.5×10⁻³ | 2.0 | 1.000 | −0.013 | t = 57 |
| over_A0.2_R400 | 1D slabs | 3.0×10⁻³ | 2.0 | 1.000 | 0.086 | t = 264 |
| over_A0.12_R100 | 1D slabs | 5.4×10⁻⁴ | 0.72 | 1.000 | −0.233 | t = 6 |
| 3D lump A0.2 R30 | isolated, 128³ | — | 2.0 | 1.000 | 0.115 at t = 66 | t ≈ 6 (peak 0.428) |
| 3D flat-top lump A0.2 R40 (edge 8) | isolated, 128³ | — | 2.0 | 1.000 | 0.024 at t = 56 | t ≈ 5.4 (max after t = 5: 0.379) |
| 3D small lump A0.02 R6 | 64³ periodic (one lump per box) | 1.9×10⁻⁸ | 0.02 | 1.000 | 0.294 | t = 2.7 |
| 3D small lump A0.02 R20 | 64³ periodic | 6.8×10⁻⁷ | 0.02 | 1.000 | 0.172 | t = 5.4 |

The excess lasts about the light-crossing time of the coherent region (0.25–0.66 R₀ in 1D) and the long-time
averages are all below 1/3. Every run satisfies the virial identity (relative residuals 2×10⁻⁵–2.4×10⁻⁴ in 1D).

**Dilatation theorem (isolated configurations, any U ≥ 0 with U(vacuum) = 0) [T].** Let D(t) = ∫x·𝒫 d³x with
𝒫_i = T^{0i} = −Σ ḟ ∂_i f. Then dD/dt = ∫T^{ii} = 3P_tot(t). Canonical scalars with U ≥ 0 obey |𝒫| ≤ ρ, and the
energy stays inside radius R₀ + t, so

    ∫₀ᵀ (3P_tot − E) dt = D(T) − D(0) − ET ≤ E R₀ + |D(0)| ≤ 2E R₀,   i.e.   w̄_T ≤ 1/3 + (R₀ + |D(0)|/E)/(3T).

A kinetic lump has D(0) = 0 (∇f = 0 where ḟ ≠ 0), so w̄_T ≤ 1/3 + R₀/(3T): **any localized excess, whatever the
potential, is transient on the light-crossing time.** 3D check (A = 0.2, R₀ = 30, 128³, dx = 1.5): dD/dt = 3P_tot
holds to 8.7×10⁻⁴ of E·T (flat-top R₀ = 40: 3.3×10⁻³), the bound is never violated, and w̄_T falls below 1/3 at
t ≈ 6 (sims_3d_lumpdil_*.json). The homogeneous over-barrier period at ε = 2λv⁴ is ≈ 49, comparable to the 56–66
isolated window, and edge disturbances reach most of the lump volume within one half-period. So the 3D runs
test the bound but do not show a sustained 3D excess; the 1D slabs do (R₀ = 400: excess until t = 264). The theorem does not apply to planar (1D) or homogeneous configurations; there
the scalar virial identity governs.

## 5. Relativistic wall collisions and the complement of F

The walls themselves make locally over-barrier states. 1D kink–antikink collisions (sims_1d B/E) [C]:

| γ | max \|φ\| after collision | max instantaneous ⟨Φ₀⟩/ρ after | post-collision w̄ | late post-collision w̄ |
|---|---|---|---|---|
| 1.15 | 1.74 v | +1.88 | −0.479 | — |
| 2.29 | 2.08 v | +1.23 | 0.052 | — |
| 7.09 | 2.60 v | +0.35 | 0.257 | — |
| 22.4 | 2.83 v | +0.10 | 0.312 | 0.309 |
| 70.7 | not executed (exceeded the 30-min cost cap) | | | |

So debris with |φ| > v and instantaneous ⟨φ²⟩ > v² exists; its time average rises toward 1/3 from below with γ in
these planar runs. A collapsing horizon-size closed wall reaches γ ~ H⁻¹/δ ≈ (4.9–5.9)×10³⁰ (cosmo_numbers.json), far
outside the tested range.

**Energy-budget argument [A].** In a scaling network the walls contribute (3p − ρ)_w = −3(1 − ⟨u²⟩_E)ρ_w ≈ −(2.1–2.4)ρ_w
(⟨u²⟩_E ≈ 0.21–0.3 [E]). Two routes show sub-horizon over-barrier debris cannot reverse the sign:
1. *Equation-of-state route.* With w_s ≤ 0.4415 (homogeneous maximum), they would need ρ_s > 6.5–7.3 ρ_w on
   Hubble average.
2. *Dilatation route.* If every joule released by the network (Q ≈ 0.188 ρ_w/t, Codex's corrected value) passes
   once through isolated lumps of size ≤ R_max, their total excess is ≤ 2R_max Q ≈ 0.38 (R_max/t)ρ_w. Reversing the
   sign needs R_max ≳ 6t ≈ 3H⁻¹: super-horizon.

Not established: isolation of debris from the network; FRW corrections O(HR₀) to the dilatation identity for
R₀ ~ H⁻¹; repeated re-concentration of radiated energy; γ ≫ 70.

## 6. (d) Strongest honest statement

> Within the unchanged action treated as classical, minimally coupled fields (the ξRφ² effect is ~10⁻³³ [E]),
> with gravity the only coupling to visible matter and a statistically homogeneous FRW background, the hidden
> energy obeys **exactly** Δ(ρ_h a⁴) = −∫H a⁴⟨λv²(φ² − v²) − m_D0²D²⟩dt + B, with |B| ≲ (2/√3)(v/M_P) ρ_h a⁴ per
> e-fold (≈ 5×10⁻¹⁶). Hence ρ_h a⁴ does not decrease between parity and recombination if and only if the
> Hubble-weighted spatial mean of φ² stays ≤ v² + (m_D0²/λv²)⟨D²⟩. No averaging over oscillations has to be
> assumed: instantaneous w > 1/3 states (Codex's counterexample) enter only through B.
>
> This condition holds **[T]** for every configuration in which φ keeps one sign or |φ| ≤ v. It holds
> **[computed or leading order]** for walls at any speed (−3σ/γ per unit area; NG w = ⟨u²⟩_E − 2/3), for waves and
> quanta of φ, X, D about either vacuum, for stationary or moving localized lumps including Q-balls (von Laue,
> w = ⟨u²⟩_E/3), for the restored charged condensate, for every lattice state run here, and for Codex's 12 3D runs.
>
> It can fail only through coherent over-barrier regions (local ρ ≳ 0.59 λv⁴; homogeneous w̄ ≤ 0.4415) that hold,
> on Hubble average, more than ≈ 7 times the wall energy. Such states are created locally by the walls themselves
> (relativistic collisions; 1D debris reaches |φ| ≈ 2.8 v at γ = 22), so the mean-density argument does not
> exclude them. Isolated sub-horizon regions lose any pressure excess within their light-crossing time
> (dilatation theorem), and on the scaling energy budget they cannot outweigh the wall deficit unless they are
> super-horizon or re-concentrated **[A]**. Under this condition, and if a persistent network reaches parity,
> ρ_h ≥ ρ_γ+ν at recombination.

Labels: the identity and T1–T3 are theorems; the wall/lump catalogue is computed; the post-parity network
statement is [A].

## 7. (e) ΔN_eff-equivalent: meaning and wording

ΔN_eff,eq ≡ ρ_h/ρ_1ν with ρ_1ν = (7/8)(4/11)^{4/3} ρ_γ = 0.22711 ρ_γ. Then ρ_h ≥ ρ_γ + ρ_ν(N_eff = 3.044) ⟺
ΔN_eff,eq ≥ 1/0.22711 + 3.044 = **7.447** (7.403 with N_eff = 3).

Suggested wording: "If a persistent network reaches parity, the hidden energy density at recombination is at
least the standard radiation density (ρ_h/ρ_rad ≥ 1). In energy this equals ≥ 7.4 extra neutrino species. This is a
radiation-equivalent energy ratio, conditional on §6; it is not an N_eff likelihood. The hidden component would
be walls (w ≈ ⟨u²⟩_E − 2/3), non-relativistic lumps and radiation, inhomogeneous on horizon scales. A CMB test needs
its own background, perturbation and likelihood model, which neither packet computes. Planck's N_eff is quoted only
as context." Drop "the least constraining case": w̄ = 1/3 gives the *minimum energy ratio*, but whether it is
least constraining for a CMB likelihood was not computed.

## 8. Corrections to Claude's packet (c7b8480)

1. REPORT §0 item 4, §5(b), §6 reply and history.py line 107: replace "the energy-weighted mean w̄ of the hidden
   sector stays ≤ 1/3" by the Hubble-integral statement of §6, including its condition.
2. REPORT §0 item 4: "Only coherent oscillations spanning both vacua have w > 1/3 (up to 0.44), and they need an
   energy density about 30 orders of magnitude above that at parity, so they cannot occur after parity" → the
   §2 criterion (space–time mean φ² > v²), the homogeneous threshold 0.588 λv⁴ and maximum 0.4415, and the fact
   that walls create such states locally; drop the mean-density inference.
3. VERIFY claim 7: "needs ρ ≳ 0.0225 v⁴ ≈ 2×10¹⁰ GeV⁴" → ε* = 0.00588 v⁴ = 5.9×10⁹ GeV⁴; "every hidden component
   has w ≤ 1/3" → "time-averaged, under §6's condition"; keep "Verdict: sound" only for single-sign/|φ| ≤ v states.
4. "ΔN_eff ≥ 7.4 … (the least constraining case; w̄ < 1/3 is more constraining)" → §7 wording.
5. "Particle gases … p ≤ ρ/3" stays, for on-shell m² ≥ 0; "relativistic quanta have w = 1/3" → "< 1/3 if
   massive, = 1/3 if massless".

## 9. Limits

- Classical fields. The quantum trace anomaly and renormalized ⟨φ²⟩ are omitted; effects of order β_λ/λ ~ 10⁻²
  of the wall term are expected [E], not computed.
- FRW background assumed. Claude's report itself expects O(1) horizon-scale inhomogeneity after parity. The
  ρa⁴ bookkeeping and the recombination comparison then need a non-perturbative GR treatment. This limits the
  conclusion itself, not only the w bound.
- The complement-of-F argument (§5) is [A]: debris isolation, γ ≫ 70, FRW corrections to the dilatation bound for
  R₀ ~ H⁻¹ and re-concentration are not computed.
- Reduced grids: 1D planar N ≤ 22800, 3D 64³–128³ (identity residuals up to 0.8% in 3D from the KDK step); flat or prescribed radiation background; λ = 0.01 benchmark;
  D = 0 on lattices (D tested only in 0D); seed noise ≤ 10⁻³; no gravitational waves or metric perturbations.
- The 1D planar runs are not covered by the dilatation theorem.
- Codex's diagnostics have 41 samples per run (Δτ = 1.2), so codex_eos.json is coarse.
- Planck 2018 N_eff value: NOT CHECKED (not this track).

## Files and reproduction

All in `audit/wbound/` (numpy/scipy/sympy; single-threaded):
- `algebra.py` → results/algebra.json: (a), Φ_{+v}, Φ₀, origin shift, Euler decomposition, static-kink identity,
  ΔN_eff,eq. Runs in under 1 min.
- `zero_d.py` → results/zero_d.json: homogeneous quadrature (threshold, maximum) and 384 random φ/X/D orbits
  (≈ 11 min).
- `lattice.py`: KDK lattice (1D planar/3D, φ + complex X, flat or radiation FRW, identity diagnostics,
  dilatation moment).
- `sims_1d.py A B C D E` → results/sims_1d_*.json: one-sign waves, boosted walls, quench, kinetic slabs,
  ultra-relativistic collisions (`WB_E_SPEEDS=0.999` for E; γ = 70.7 exceeded the cost cap, log kept as
  results/sims_1d_E_killed_at_cap.log). Families A–C come from the first pass (sims_1d_ABCD.json); D was rerun
  with seed noise confined to the slabs (sims_1d_D.json supersedes the D block of sims_1d_ABCD.json, whose
  linear-lump rows were dominated by white lattice noise).
- `sims_3d.py lumpdil tophat net small charged frw` → results/sims_3d_*.json (small rerun with noise confined to
  the lumps).
- `codex_eos.py <extracted raw dir>` → results/codex_eos.json (read-only use of a69bd86 raw_diagnostics).
- `cosmo_numbers.py` → results/cosmo_numbers.json.
