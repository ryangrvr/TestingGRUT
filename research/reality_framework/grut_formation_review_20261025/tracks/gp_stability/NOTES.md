# Track `gp_stability`: gauge stability, screening and the discrete Z6 charge at the stated gauge point

Claude review of Codex b5174e5 (`grut_gauge_point_test_20261025`, handoff item 3) and its parent 7737487
(`grut_z6_gauge_bridge_20261025`). 25 October 2026. Packets read only through `git archive` into the
scratchpad. No Codex code was run (no reproduction check was done); every [C] number below is from this
track's own code.

Labels: [D] derived, [C] computed here, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

Stated point [Q]: M_Y = 1e12 GeV, h = 1e11 GeV (|Sigma| = h/sqrt2), mu = 5e11 GeV, y (kappa) = 0.001,
e = 1e-7, lambda_Y = lambda_Sigma = 1; light target g_X = 0.1, lambda_X = 0.001, lambda_phi = 0.03,
v = 0.924089637 GeV, Q = 200000. Charges X : Y : Sigma = 1 : 3 : 6.

This is one supplied working model. Nothing here is evidence that the fields exist.

## 1. Method

1. **Own object profile** (`code/profile_fixed_charge.py`). The frequency is an unknown parameter of the
   boundary-value problem and the charge Q = 200000 is imposed as a boundary condition on the running
   charge. There is no root finding. This differs from Codex's Background class and from track
   `conversion`. Two box sizes were used (L = 70 and L = 90 in units of 1/m_X).
2. **Gauge energies** (`code/gauge_energy.py`).
   - The screened (Yukawa) energy uses a real-space radial Green's function, exact for a piecewise-linear
     source. Codex used a Fourier/Parseval route.
   - The local limit and its leading gradient correction are computed.
   - The unscreened Coulomb energy is computed two ways: Green's function, and enclosed charge.
   - Also computed: central potentials, thin-wall Q_max scalings, heavy masses and the Z6 algebra.
3. **Z6 string tension** (`code/z6_strings.py`). An own Nielsen-Olesen solver in s = ln(m_v r), which
   resolves beta = m_s^2/m_v^2 up to the stated point's 5.6e12.
   - Checks: eps(1) = 1 (BPS) to 6e-11; the magnetic and potential parts are equal (Derrick); the result
     agrees with the Hill-Hodges-Turner fit 1.04 beta^0.195 to within 4% for beta between 0.5 and 100.

Reproduce: `cd code && python3 profile_fixed_charge.py && python3 gauge_energy.py && python3 z6_strings.py`
(about 1 + 1 + 15 minutes). Results: `results/*.json`.

## 2. Results

### 2.1 Object profile (input to everything below)

- Omega = omega/m_X = 0.86119513082 [C]. E = 17064.5374544 GeV [C] (Codex 17064.537454786 [Q]).
- Binding margin 200000 m_X - E = 1417.2552931 GeV [C] (Codex 1417.2552927 [Q]; relative difference 3e-10).
- L = 70 and L = 90 agree to 1e-13.
- Charge radii [C]: half-charge radius 40.97 GeV^-1 (8.1e-15 m), rms radius 45.1 GeV^-1, 99% radius 87.9 GeV^-1.

### 2.2 Vector mass and screening

- **m_A = 6 e h = 6.000e4 GeV** [D]. From |D Sigma|^2 with D = d + 6ieA and |Sigma| = h/sqrt2,
  m_A^2 = 2(6e)^2 h^2/2 = 36 e^2 h^2. CONFIRMED. The X condensate in the core adds e^2 f0^2 = 5.9e-14 GeV^2
  to m_A^2 [C], which is negligible.
- Screening length 1/m_A = 1.67e-5 GeV^-1 = 3.3e-21 m [C]. m_A times the half-charge radius = 2.46e6 [C].
  The object is deep in the locally screened (Higgs) regime.

### 2.3 Screened energy

- **Local bound.**
  - int n^2 d^3x = 47241.833 GeV^3 [C] (Codex 47241.832 [Q]).
  - Delta E_local = int n^2/(72 h^2) = **6.56137e-20 GeV** [C] (Codex 6.56137e-20 [Q]).
  - This equals the Higgs-phase kinetic energy of an E = 0, n_Sigma = -n/6 slice [D].
  - It is independent of e.
  - In effect it adds a quartic term of strength delta lambda = omega^2/(18 h^2) = 3.5e-26, which is
    3.5e-23 of lambda_X.
- **Variational minimum** (fixed Higgs modulus) [D]:
  - (-lap + m_A^2) Phi = e n, giving Delta E_min = (e^2/2) int |n_k|^2/(k^2 + m_A^2).
  - Its relative distance below the local bound is at most int|grad n|^2/(m_A^2 int n^2) = 5.7165e-13 [C]
    (Codex 5.7165e-13 [Q]). The bound is correct [D].
- **Real-space Yukawa check against Codex's Fourier ratios** [C vs Q], for m_A/m_X = 0.1, 1, 10 and 1000:
  - mine: 0.08165947, 0.82486394, 0.99760116, 0.99999975;
  - Codex: 0.08165947, 0.82486394, 0.99760117, 0.99999976.
  - The two agree to at most 6e-9, which is my grid error.
- **Against the binding margin:** 6.56e-20/1417.26 = **4.63e-23** [C]. CONFIRMED.

### 2.4 Unscreened upper bound (new; Codex gave none)

- **Coulomb energy** with m_A set to 0 (unbroken U(1), same e):
  - E_C = **3.4005e-7 GeV** [C], by Green's function and by enclosed charge; the two agree to 2e-9.
  - A uniform sphere with the same rms radius gives 3.28e-7 GeV.
  - E_C/margin = 2.4e-10.
  - Higgs screening lowers the energy by 5.2e12 [C].
- **Central potentials.** e Phi(0) = 3.9e-14 GeV unscreened and 1.3e-24 GeV screened [C]. The gap
  m_X - omega is 0.0128 GeV.
- **Thin-wall scaling** [E] (Lee-Stein-Schabes-Watkins-Widrow 1989 type argument, E_C/Q proportional to
  Q^(2/3) at fixed density). For an unbroken U(1) the object would stop being bound at about Q ~ 5e19.
  The central potential would reach the gap at about Q ~ 4e22.
- LSSWW-type existence condition [E]: Coulomb energy per charge must stay below the binding per charge at
  the thin-wall onset radius. This gives alpha <~ b (m_X^2 - omega^2)/n0 = 3.4e-5, i.e. e_max ~ 2e-2;
  from the Q_max scaling, e_max ~ 6e-3. The stated e = 1e-7 is 1e5 times smaller.
- So even unscreened, the gauge coupling does not affect classical binding at Q = 2e5.

### 2.5 Gauged-Q-ball stability criteria [A + C]

- m_A R ~ 2e6, so this is the heavy-vector limit of Proca and Higgs-gauged Q-balls:
  - Heeck-Rajaraman-Riley-Verhaaren, JHEP 10 (2021) 103, arXiv:2107.10280;
  - charge screening by a Higgs: Ishihara-Ogawa arXiv:1811.10894; Forgacs-Lukacs PRD 102 076017.
  - There the object behaves like the global-symmetry Q-ball. The literature identifiers were checked by
    web search [Q]. The results are quoted from abstracts, not re-derived.
- The global criteria therefore carry over with corrections of order 1e-23:
  - E < Q m_X (margin 1417 GeV);
  - dE/dQ = omega < m_X (0.861);
  - dQ/domega < 0. Track `conversion` (our earlier, independent profile) found -4.42e6 [Q-own], which means
    classical stability and stability against fission.
- Not certified: a full coupled fluctuation spectrum, quantum metastability and thermal histories.
  Codex says the same.

### 2.6 Heavy vacuum masses vs the object's energy

- Masses [C]: m_A = 6.0e4 GeV; Higgs radial sqrt(2 lambda_Sigma) h = 1.414e11 GeV; mediators
  M_Y sqrt(1 -/+ sqrt2 mu h/M_Y^2) = 9.640e11 and 1.0348e12 GeV.
- The minimum is 6.0e4 GeV, which is **3.52 times** E = 17064.5 GeV. So an isolated cold object cannot
  produce any on-shell heavy quantum. CONFIRMED.
- Margins and side numbers:
  - The vector margin is only a factor 3.5. Four objects' rest energy exceeds m_A, which matters only
    for mergers.
  - The vector decays quickly: Gamma(A -> X Xbar) = e^2 m_A/(48 pi) = 4.0e-12 GeV, so tau = 1.7e-13 s [D].
  - The Y induced in the core is about kappa |X0|^3/M_Y^2 = 5e-27 GeV [C], so the Y = 0 trial is accurate.

### 2.7 The mod-6 statement

- Arithmetic [C]: 200000 = 6 x 33333 + 2, so Q is 2 mod 6. Three cores are the smallest neutral multiple
  (600000), with Sigma number -100000. CONFIRMED.
- Correct invariant [D]: on a compact slice, Q_X + 3Q_Y + 6Q_Sigma = 0. The residual Z6 charge is
  therefore (Q_X + 3Q_Y) mod 6, not Q_X alone. Codex's code and parent state the full identity; the b5174e5 prose speaks of Q_X alone.
  At late times Y = 0, so this changes nothing in practice.
- Scope [A]:
  - The constraint applies to the **whole compact slice**: the closed universe GRUT assumes. It does not
    apply to single objects or to bound groups, and it does not apply at all to non-compact space.
  - Compensation can sit anywhere: four extra X quanta or two anti-X quanta anywhere in the universe
    suffice. No grouping into triples is required.
  - Locally, the Higgs condensate screens the full charge in expectation, including a fractional
    -Q/6 = -33333.3. The residual 2 mod 6 is visible only as an Aharonov-Bohm phase (quantum hair:
    Krauss-Wilczek 1989; Preskill-Krauss 1990; Coleman-Preskill-Wilczek 1992 [R]).
  - Verdict: CONFIRMED, with that clarification. As a population constraint it is trivially satisfiable
    and selects nothing.

### 2.8 Observable consequences of the discrete Z6

1. **Selection rule** [A]. In vacuum, only Delta Q_X = 0 mod 6 processes occur. X^5 is forbidden exactly
   if the Z6 is a genuine gauge remnant (Krauss-Wilczek; Banks-Seiberg [R]). The pair channel removes 4
   units from the core and 2 escape, 6 in total, consistent with this. This is the protection itself; it
   is seen only through the leakage rates.
2. **Short-range force.** The force mediated by the 60-TeV vector has range 3.3e-21 m. With zero kinetic
   mixing (assumed), it has no visible-sector signature.
3. **Aharonov-Bohm phase on Z6 strings.** The phase is exp(2 pi i Q/6) = exp(2 pi i/3) for the object
   (sin^2(pi nu) = 3/4). Transport cross-section per unit length is about 2 sin^2(pi nu)/k: about
   0.09 GeV^-1 for v ~ 1e-3 [E]. That is negligible at cosmic string densities.
4. **Cosmic strings** (only if a network forms; see 2.9). This is the one potentially testable signature.

### 2.9 Z6 string tension and G mu (retrospective comparison)

- beta = m_s^2/m_v^2 = 2 lambda_Sigma/(36 e^2) = 5.56e12, independent of h [D]. These are very strongly
  type-II strings: flux tube 1/m_A = 1.7e-5 GeV^-1, scalar core 7e-12 GeV^-1.
- eps(beta) at large beta: eps - (1/2) ln beta -> 0.15024 [C] (converged to 1e-6 between beta = 1e8 and 1e10).
- At the point: eps = (1/2) ln beta + 0.150240 = **14.823** [C], with an extrapolation uncertainty of about 1e-6.
  Direct collocation at beta = 5.6e12 is too stiff for solve_bvp.
- mu = pi h^2 eps = 4.657e23 GeV^2. **G mu = 3.12e-15** [C], with M_Pl = 1.2209e19 GeV.
  The naive BPS value pi h^2 would give 2.1e-16, so the logarithm matters by a factor of about 15.
- Limits [Q] (retrospective, Nambu-Goto-like loops assumed):
  - Planck 2013 XXV: 1.5e-7 (Nambu-Goto), 3.2e-7 (Abelian-Higgs). Passes by about 5e7.
  - LIGO-Virgo-KAGRA O3 stochastic background (arXiv:2101.12248), excluded above:
    - model A (Blanco-Pillado-Olum-Shlaer): (9.6e-9 to 1e-6);
    - model B (Lorenz-Ringeval-Sakellariadou): (4.0 to 6.3)e-15;
    - model C-1: (2.1 to 4.5)e-15;
    - model C-2: (4.2 to 7.0)e-15.
  - **So the point lies just below model B's bound (ratio 0.78) and inside the C-1 excluded range for
    some kink numbers (ratio 1.49 against its lower end).**
  - At fixed beta, G mu is proportional to h^2. h <~ 1.13e11 GeV clears model B; h <~ 8.2e10 GeV clears
    all of C-1 [C]. Pulsar timing (about 1e-10 [R]) is irrelevant. LISA (about 1e-17 [R]) would probe it.
- Conditions [A]:
  - A network forms only if the U(1) is restored after inflation and re-broken at T ~ h. Thermal
    restoration needs T_rh above roughly 1e11 to 3e11 GeV; this is plausible for the review's high-scale,
    fast-reheating example but was not computed.
  - The Affleck-Dine condensate itself carries U(1) charge and may keep the symmetry broken early. Its
    amplitude during inflation is |Sigma| ~ 1.25e13 GeV [Q].
  - If strings form, this is a near-edge gravitational-wave signature that Codex's packet does not mention.

## 3. Verdicts

| Claim (Codex) | Codex | This track | Verdict |
|---|---|---|---|
| m_A = 6 e h | 6.0e4 GeV | 6.0e4 GeV [D] | CONFIRMED |
| Screening length << core | (m_A/m_X = 6.49e5) | 1/m_A = 1.7e-5 GeV^-1 vs r_1/2 = 41 GeV^-1; ratio 2.5e6 [C] | CONFIRMED |
| Screened energy bound | <= 6.56137e-20 GeV | 6.56137e-20 GeV [C]; Yukawa ratios agree to 6e-9 | CONFIRMED |
| Margin and ratio | 1417.26 GeV; 4.63e-23 | 1417.2553 GeV; 4.63e-23 [C] | CONFIRMED |
| Unscreened bound | not given | 3.40e-7 GeV = 2.4e-10 of margin [C] | new (passes) |
| Classical binding survives gauging; full spectrum open | PARTIAL | Heavy-vector limit inherits global criteria; corrections 1e-23 [A,C] | CONFIRMED (scope as stated) |
| Heavy masses exceed object energy | yes | min 6.0e4 = 3.52 E [C] | CONFIRMED |
| Q = 2 mod 6; three cores neutral | yes | yes [C] | CONFIRMED |
| "Compact population needs compensating sectors" | stated | True for the whole closed slice only; invariant Q_X + 3Q_Y; compensation anywhere, trivial | CONFIRMED with clarification |
| Observable Z6 consequences | not discussed | strings, G mu = 3.1e-15, near LVK O3 model-B/C bounds (conditional on formation) | NOT CHECKED by Codex; flagged |

## 4. Limits and open items

- No fluctuation/Floquet spectrum of the coupled gauged object, no quantum metastability, no thermal history.
- String network formation is conditional. Gravitational-wave bounds depend strongly on the loop model.
  Abelian-Higgs simulations radiate more into particles, which weakens the bounds [R].
- Weak-gravity caution [A/R, conjectural]. The magnetic WGC (Arkani-Hamed-Motl-Nicolis-Vafa 2007) suggests
  a cutoff Lambda <~ e M_Pl ~ 2.4e11 GeV (reduced Planck mass), which is below M_Y = 1e12 GeV and below
  H_inf = 4.7e13 GeV. This is not a test, but e = 1e-7 sits in tension with it by a factor of about 4
  to 200, depending on the scale used.
- Kinetic mixing is assumed zero. A nonzero value would turn the 60-TeV vector into a visible Z'. Not checked here.
