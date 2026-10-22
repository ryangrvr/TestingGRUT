# Track `scenario`: the solved example of construction QV1

21 October 2026. Stage 2 of the q-theory construction pass (issue #13). Briefs: `../../BRIEF.md`,
`../../BRIEF_STAGE2.md`. Stage-1 records used: `../equations/NOTES.md` (cited "equations §n"),
`../literature/NOTES.md` ("literature §n"), `../ledger/NOTES.md` and `../ledger/results/ledger.json` ("ledger").

**What this is.** One closed set of equations (QV1: the Klinkhamer-Volovik action with the state condition P_q),
integrated in full and in its validated reduction, from a Planck-scale displaced vacuum to the far future, with the
vacuum sector's gravitational response computed throughout. Every physical element is borrowed and credited
(Klinkhamer-Volovik for the action, the equations, the relaxation and the state condition; Saikawa-Shirai for the
Standard Model equation of state; Starobinsky / Vilenkin / Gorbunov-Panin for scalaron decay). Nothing here is
GRUT-native. No new theory, empirical validation or confirmatory test is claimed. Every comparison with data is
**retrospective** (all data exposed).

**Labels.** [D] derived here; [C] computed by `code/` (output in `results/`); [Q] stated by the cited source (read
level in §8); [E] estimate; [A] argued; [R] recalled, not checked against a source here; [P] supplied input.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdicts in one screen

| # | Result | Label |
|---|---|---|
| 1 | Closed system: action -> (q), (g), (m) -> flat-FRW first-order system with q, H, SM plasma, dust, radiation. The FRW Bianchi identity dC/dt = -4HC holds identically (also with the dissipative extension, whose matter source it fixes); numerically \|C\|/scale <= 4e-13 in every QV1 run (<= 1e-11 in the strongly misaligned light-q runs of §4.2) and C a^4 drifts <= 1e-11. | [D][C] |
| 2 | **Reduction validated.** The full q dynamics at m_q/H_c = 10...1e4 reproduces "GR + V_m(t) - c*(mu)" with error eps0 = 0.030 (H_c/m_q)^2 (fitted slope -2.007 over 1e2...1e4; prediction and measurement agree to < 1% for m/H >= 300). Adding the adiabatic q-correction (R1) removes it to <= 7e-8 (m/H = 100), <= 2e-10 (300), <= 1e-11 (>= 1000, numerical floor). | [C] |
| 3 | **Deposited oscillation energy** = (4/3) Delta^2/(K m_q^2) (x/sinh x)^2, x = pi m_q tau/2, confirmed: measured/formula = 1.0000 +- 0.0004 for x <= 2 and m/H >= 100; at x = 8 the ratio is 3.35, 1.028, 1.012, 1.002 for m/H = 30, 100, 300, 1000 (finite-H correction, vanishing as m/H grows). Fitted E ~ exp(-kappa x) with kappa = 2.000-2.004 for m/H >= 100 (formula 2.005 on the same points). Scales as Delta^2 to 0.1%. <rho_V> = S^2/(12 K m_q^2) (= rho_m^2/(12Km^2) in a dust era) and = KV's (9/8) kappa^2 H^4 K'^2/eps'' to O(H^2/m^2). | [C] |
| 4 | **Early relaxation.** KV 0806.2805's <rho_V> = (2/3) d_M^2 M^2/t^2 reproduced (t^2<rho_V> = 0.6620 vs 0.6667, d_M = 1; 0.6497 vs 0.6495, d_M = 0.987; 0.2335 vs 0.2332, d_M = 0.591). With scalaron decay (effective source, Gamma = 4 m^3/(192 pi M^2) = 4.4e-3 m_q): the oscillation energy decays at exactly Gamma and turns into radiation; a q displaced at rest hands 0.96-1.06 of its initial vacuum energy (per comoving volume) to radiation for \|y_i\| <= 0.1 and 0.85-1.22 for -0.25 <= y_i <= +1 (larger, frame-dependent values for y_i <= -0.5, where K -> 0); reheating T_rh = (4.8-5.2)e16 GeV for every displacement (Gorbunov-Panin formula: 5.03e16). EFT control at Planck densities is absent (§3.4). | [C][Q/E] |
| 5 | **Realistic history, natural m_q = 1.988e18 GeV.** From 1 TeV to the far future: H/H_GR - 1 <= 2.6e-64 (1 TeV), 1.4e-65 (EW), 4e-77 (QCD), 3e-92 (BBN), 1.8e-121 (today); rho_V(q) = 6.9e-55 GeV^4 at 1 TeV, 6.8e-58 at the EW crossover, 4.7e-169 today; delta G/G = -1.6e-121 and Gdot/G = 1.9e-131 /yr today. No constant remnant; asymptotic state Minkowski (t H -> 2/3, rho_V -> 0 as a^-6). | [C] |
| 6 | **Bracket m_q = 4 meV** (lowest unscreened value): the reduction holds below T ~ 100 GeV (m/H >= 300), but the q-mode is light above ~1 TeV, gets misaligned by the SM trace anomaly and leaves a stable coherent oscillation: Omega_q h^2 = 378 (adiabatic start at 1 TeV) to 4e5 (q at q0 from 10-1000 TeV), i.e. 3e3-3e6 times Omega_c h^2 (115 to 2e6 across the trace-anomaly uncertainty). Same overproduction for every m_q from 4 meV to 1 GeV computed (Omega_q h^2 = 2e3 to 6e9). **Retrospective failure of the light-q bracket for untuned initial data**; viable unscreened m_q needs decay before BBN, m_q >~ 2e5 GeV [E] (moduli-problem bound, credited). | [C][E] |
| 7 | **QV1 (P_q) late universe has no vacuum term**: H0 = 37.74 km/s/Mpc (from T0 and the CMB physical densities), age 17.27 Gyr, q0 = +0.500. Pantheon+ Delta chi^2 = 699 against flat LambdaCDM (shape, M free); DESI DR2 BAO shape Delta chi^2 = 1412, absolute chi^2 = 15333 for 13 points; 100 theta_* = 0.943 vs 1.0411 +- 0.0003 (-327 sigma); H0 vs SH0ES -34 sigma; q0 vs SH0ES +42 sigma. **QV1 fails retrospectively.** The failure is reported, not repaired inside QV1. | [C], retrospective |
| 8 | **QV1-cal** (mu = mu0[eps_inf] - rho_Lambda/q0, calibrated to rho_Lambda,obs = 2.52e-47 GeV^4) reproduces flat LambdaCDM: H0 = 67.36, age 13.81 Gyr, q0 = -0.529, chi^2 equal to Planck-LambdaCDM's within 0.01 (SN) and 0.5 (BAO absolute). **The calibration fixes the residual**; second-order correction (4/3) rho_L/(K m_q^2) = 1.4e-120. | [C] |
| 9 | **Rejected readings.** Initial-equilibrium mu leaves rho_V,inf = -Delta_EW = -4.7e54 rho_Lambda: turnaround at T = 44.6 GeV, t = 2.02e-10 s, crunch at 4.05e-10 s (never reaches the QCD crossover). Fixed-number mu: a heavy threshold c moves rho_V,inf one-for-one (rho_V,inf/c - 1 = (4/3)c/(K m_q^2) = 1.4e-120 ... 1.4e-60); a 10 TeV scalar (c = -2.37e13 GeV^4) gives a crunch at 8e-13 s (or de Sitter domination below 0.91 TeV for the opposite sign). Under P_q the same c is re-absorbed: rho_V,inf = 0 exactly; only K(q0) shifts, delta G/G = -chi c = -3.4e-66 (EW) ... -6.8e-61 (10 TeV), absorbed in the measured G. Confirmed with the full equations at m/H = 100 (§5.4). | [C][D] |
| 10 | **What emerges / what does not.** Emerges from the equations: cancellation of every constant under P_q (any order, any sign within c >= -M^4/8), adiabatic tracking through transitions, no constant remnant (second-order terms only as decaying transients), Planck-era relaxation into radiation (with the decay source). Does not emerge: the observed residual (QV1 gives exactly zero; QV1-cal inserts it by calibration), dark matter (the q-oscillation decays at natural m_q and overcloses when light). | [D][C][A] |

---------------------------------------------------------------------------------------------------------------------

## 1. The closed system that is solved

### 1.1 Action and field equations (credited to Klinkhamer-Volovik)

    S = ∫ d^4x sqrt(-g) [ (1/2) K(q) R - eps(q) ] + S_m[psi, g],      F = dA,  q = F~/sqrt(-g)

KV 0806.2805 eqs (2.1)-(2.8); 0905.1919 eqs (2.1)-(2.2). The field equations are those confirmed in stage 1
(equations §1; literature §2):

    (q)  eps'(q) - (1/2) K'(q) R = mu                       mu: one global integration constant (the state)
    (g)  K G_mn = T_mn - (eps - mu q) g_mn - (g_mn Box - nabla_m nabla_n) K
    (m)  nabla^m T_mn = 0
    rho_V = eps(q) - mu q     (gravitating vacuum term)

Concrete functions (QV1: generic natural-scale functions, KV forms where a concrete choice is needed) [P, credited]:

    eps(q) = (1/2chi) ( -q^2/q0^2 + q^4/(3 q0^4) )     (KV 0806.2805 eq (4.7))
    K(q)   = q/3                                        (KV (4.14) with sigma = 1, so K0 = K(q0) = q0/3 = M^2)

with the natural choice chi = 1/M^4 (q0^2 eps''(q0) = M^4, q0 K' = K0 = M^2). Then [D, checked in s1]:

    q0 = 3 M^2,  mu0 = eps'(q0) = -1/(3 chi q0),  eps(q0) - mu0 q0 = 0,
    m_q^2 = (2/3) K0 eps''(q0)/K'^2 = 2/(3 chi) = (2/3) M^2   =>   m_q = 1.988e18 GeV.

### 1.2 Degrees of freedom [D/R]

Graviton (2); the q-mode, a scalar of mass m_q (1; it exists because K' != 0; omega_BD = 0 scalar-tensor theory,
equations §2.2); the 3-form carries no local degree of freedom; one global constant mu per connected spacetime
(changed only by membrane nucleation, outside the action). Matter: SM plasma as a perfect fluid with a tabulated
equation of state, dust (baryons + CDM), radiation (inside the SM table; decay products in §3). Ghost-free (K > 0)
and tachyon-free (eps'' > 0) on the KV branch.

### 1.3 State condition and the matter normalisation it dictates [P, credited; D]

    (P_q)  mu = mu0[eps_inf],   eps_inf(q) = eps(q) + c_inf,   eps_inf'(q0) = mu0,  eps_inf(q0) - mu0 q0 = 0,

c_inf = every q-independent vacuum energy of the final T -> 0 matter vacuum, to all orders (KV 1102.3152 §6 for the
final-state reading; 0905.1919 §II for the form; ledger §1.4). Consequence [D]: c_inf moves into eps_inf, so the
matter stress tensor in (g) is normalised by **p_m(T = 0) = 0**. This is exactly how Saikawa-Shirai normalise their
tables (they renormalise the pressure so that it vanishes at T = 0, their §3.3 [Q]). Under P_q the gravitating
vacuum term is rho_V = eps_inf(q) - mu0 q (zero at q0, positive nearby), and every pre-transition plateau (Higgs
symmetric phase, QCD bag) sits inside rho_m - 3 p_m, where it gravitates as in GR until it is released.

### 1.4 Flat-FRW first-order system (the equations integrated) [D; s1]

Variables (units K0 = 1): Y = (K(q) - K0)/K0 = (q - q0)/q0 (so Y = delta K/K = -delta G/G), V = dY/dt, H, N = ln a,
J = ∫H^2 dt (bookkeeping for N-space averages), matter state. With S = rho_m - 3 p_m (minus the trace of T^m):

    dY/dt = V
    dV/dt = S/3 + (c/3)[(4 - 2 beta) + 2 beta Y] - m^2 F(Y) - 3 H V - Gamma V/(1+Y)
    dH/dt = R/6 - 2 H^2,    R = 3 m^2 Y (1 + 1.5 b3 Y + 2 b4 Y^2) + 2 beta c + 3 Gamma V/(1+Y)^2
    rho_V(Y) = (3/4) m^2 Y^2 (1 + b3 Y + b4 Y^2) + c (1 + beta Y),   F(Y) = Y + 1.5 b3 Y^2 + (b3/2 + 2 b4) Y^3 + b4 Y^4
    matter: d rho_d/dt = -3H rho_d;  d ln rho_pl/dt = -3H (1 + w(ln rho_pl));  toy: d rho_r/dt = -4H rho_r - dDeltaV/dt
    constraint: C = 3 H^2 (1+Y) + 3 H V - rho_m - rho_V = 0

- KV shape: b3 = 4/3, b4 = 1/3, beta = K0/(K'q0) = 1 [D, s1]; exact, since eps is quartic.
- c = -(mu - mu0) q0 = 0 under P_q. Non-zero only for QV1-cal (c = rho_Lambda) and the rejected readings. A constant
  added to eps at fixed mu (a heavy threshold) is equivalent to a constant in rho_m and S ("type i").
- Gamma = 0 in QV1; Gamma > 0 only in §3 (QV1+Gamma, §1.6).
- **Bianchi identity in FRW form** [D, s1]: dC/dt = -4 H C identically, for generic eps(q), K(q) (K'' kept), with or
  without dissipation provided its energy enters the matter equation (residual 0; without the source the residual is
  Gamma_q qdot^2). Numerically C and C a^4 are monitored in every run.
- **Adiabatic reduction** [D, s1]: inserting Y = S/(3 K0 m^2) + O(m^-4) into C:

      R0:  3 K0 H^2 = rho_m                       (GR + V_m(t) - c*(mu); under P_q c* = c_inf, already removed)
      R1:  3 K0 H^2 = rho_m + delta rho_q,   delta rho_q = [ S^2/(12 K0) - H^2 S - H dS/dt ] / m^2 + O(m^-4)

  i.e. 3 K0 H^2 (1 + (S + S_N)/(3 K0 m^2)) = rho_m + S^2/(12 K0 m^2), S_N = dS/dN. The q-sector pieces: rho_V =
  S^2/(12 K0 m^2) (equations track's rho_m^2/(12 K m^2) in a dust era; identical to KV's (9/8) kappa^2 H^4 K'^2/eps''
  with kappa = S/rho_m when 3K0H^2 = rho_m) and the G-variation terms -3 K0 H^2 Y - 3 K0 H Ydot. The q-mode alone:
  Ydd + 3H Yd + m^2 F(Y) = S/(3 K0), an oscillator of frequency m_q driven by the matter trace.

### 1.5 Matter content

- **SM plasma** [Q]: Saikawa-Shirai 1803.01038 tables (g_*rho, g_*s with errors; w = p/rho for T >= 10 MeV),
  p(T=0) = 0. Trace S_pl = rho_pl (1 - 3w) for T >= 10 MeV (quoted error carried: +-1.5e-3 absolute in S/rho above
  300 GeV, larger than the value 2-4e-4; 1-10% elsewhere); below 10 MeV (neutrinos decoupled) the e+- and mu+-
  ideal-gas trace at the photon temperature [C], 16% below the table at the 10 MeV junction [E]. In full-equation
  runs the plasma is evolved as d ln rho/dt = -3H(1 + w(ln rho)), so nabla.T = 0 holds exactly.
- **Dust** [Q, retrospective]: omega_b = 0.0224, omega_c = 0.120 (Planck 2018 VI abstract, a LambdaCDM fit), scaled
  with the entropy density; neutrinos massless (Sum m_nu = 0.06 eV would add 0.4% to omega_m [E]).
- **Toy transition** [P] (compressed runs, §2): radiation + dust (6% of rho at N = 0) + DeltaV(N) =
  Delta (1 - tanh(N/w))/2, Delta = 0.15 = 5% of 3 H_c^2 (also 20% and 1.25%), latent heat to radiation; duration
  tau = w/H_c.

### 1.6 Dissipative extension, used only in §3 (QV1+Gamma) [Q/E]

Scalaron decay into N_s minimally coupled scalars, Gamma = N_s m_q^3/(192 pi M^2) [Q: Gorbunov-Panin 1009.2448 eq.
(T4+) with reduced M_P, stated there to agree with Starobinsky 1980 and Vilenkin 1985]. Implemented as an effective
Rayleigh-type source (q-dot^2 dissipation of the kind KV use in 1601.04676, per literature §3.2):

    eps'(q) - K'R/2 = mu - Gamma_q qdot,   Gamma_q = (3/2) Gamma K0 K'^2/K^2;  matter gains Q = Gamma_q qdot^2;
    the metric equation keeps the constant mu.

The Bianchi identity fixes Q (s1). Linearised: Ydd + (3H + Gamma) Yd + m^2 Y = S/3, so the Einstein-frame oscillation
energy decays at exactly Gamma (measured 1.0000000 Gamma, §3). mu itself does not move: P_q and the Minkowski end
state are kept. N_s = 4 (Higgs doublet, assumed minimally coupled, xi = 0 [A]); conformal fields do not contribute at
tree level. This source is an effective description of a quantum process, not part of the classical action.

### 1.7 Initial data

- §2: quasi-static track to O(m^-4) (Y_qs = Y_alg(S) - (Ydd_1 + 3H Yd_1)/(m^2 F'), V_qs its time derivative), H from
  C = 0. The O((H_i/m)^6) residual leaves a small oscillation, measured by Delta = 0 control runs (§2.3 floor).
- §3: KV's boundary data at tau = 1, (a, h, f, r_M) = (1, 1/2, 1/2, 1/20), w_M = 1/3, i.e. Y_i = -0.5, H_i = 0.289 M,
  rho_M = 0.05 M^4, rho_V,i = 0.052 M^4, V_i from C = 0; and a displaced q **at rest** (V_i = 0, H_i from C = 0) with
  Y_i from -0.9 to +1 (Planck-scale initial vacuum energy 4e-4 ... 1.3 M^4).
- §4, natural m_q: at T = 1 TeV, q on the quasi-static track. Justified by §3: the Planck-era oscillation has decayed
  (e^(-Gamma t)), and m_q/H >= 570 from reheating on, so the SM trace cannot re-excite it.
- §4, bracket m_q = 4 meV: q on the quasi-static track at T = 1 TeV, or q = q0 at rest at T_i = 1e3...1e6 GeV. The
  outcome depends on this choice, which the construction does not fix.

### 1.8 Inputs table

| Input | Value | Provenance |
|---|---|---|
| M = K0^(1/2) | 2.435e18 GeV | brief [P]; reduced Planck mass |
| eps(q), K(q) | KV 0806.2805 (4.7), (4.14), sigma = 1 | [Q] literature §4.1; choice of form [P] |
| chi | 1/M^4 | natural-scale assumption [P] (QV1 definition); chi^-1 = E_P^4 (non-reduced) would raise m_q by 8 pi [D] |
| q0, mu0, m_q, b3, b4, beta | 3M^2; -1/(3 chi q0); 1.988e18 GeV; 4/3, 1/3, 1 | derived [D] |
| mu (QV1) | mu0[eps_inf] | state condition P_q [P], credited KV 1102.3152 §6 |
| mu (QV1-cal) | mu0[eps_inf] - rho_Lambda/q0 | **calibrated** to rho_Lambda,obs = 2.52e-47 GeV^4 (brief; Planck 2018 inputs) |
| m_q bracket | 4 meV | ledger (Eot-Wash, alpha = 1/3) [E] |
| SM EoS | Saikawa-Shirai 2018 tables | [Q] (section level, §8) |
| omega_b, omega_c | 0.0224, 0.120 | [Q] Planck 2018 VI abstract (LambdaCDM fit; retrospective use) |
| T0 | 2.7255 K | [Q] Fixsen 2009, 0911.1955 abstract (2.72548 +- 0.00057 K) |
| radiation today | g_*rho = 3.3830836 (N_eff = 3.044 equivalent) | [Q] SS table |
| Gamma | 4 m_q^3/(192 pi M^2) = 4.42e-3 m_q | [Q] Gorbunov-Panin; N_s = 4, xi = 0 [A] |
| Delta_EW | m_H^2 v^2/8 = 1.19e8 GeV^4 (m_H = 125.25, v = 246.22 GeV) | tree level [D/C] equations §7; inputs [R]; loop corrections O(10%) [E] |
| Delta_QCD | (0.145-0.235 GeV)^4 | [R] equations §7 |
| heavy threshold | 10 TeV real scalar, c = -(3/2) m_h^4/(64 pi^2) = -2.37e13 GeV^4 | illustrative [P]; one-loop MS-bar at mu_RG = m_h [R] |
| toy step | Delta/3H_c^2 = 0.05 (0.2, 0.0125), dust 0.02, widths w | [P] |
| KV boundary data | (h, f, r_M) = (1/2, 1/2, 1/20) | [Q] KV figure captions (literature §4.1) |
| Data (retrospective) | DESI DR2 BAO, Pantheon+, SH0ES H0 and q0, Planck 100 theta_*, globular-cluster age | §5, §8 |

Nothing is calibrated in QV1. The only calibrated quantity in this track is the QV1-cal offset (datum rho_Lambda,obs).

---------------------------------------------------------------------------------------------------------------------

## 2. Validated reduction (compressed hierarchies) [C] — `code/s2_reduction_validation.py`, `results/reduction_validation.json`

### 2.1 Set-up

Units K0 = 1, time 1/H_c (H_c = H at the transition centre N = 0); m = m_q/H_c = r. Full system of §1.4 (KV shape,
Gamma = 0, c = 0) with the toy matter of §1.5. Diagnostics: (i) H averaged in N over windows of k = max(1, 0.002 r)
oscillation periods, compared at equal N with R0 and R1 (eps0 = <H>/<H_R0> - 1, eps1 = <H>/<H_R1> - 1, eps2 = with R1
plus the measured oscillation energy as dust), plus pointwise values; (ii) the Einstein-frame oscillation energy about
the O(m^-4) quasi-static track, E_osc a^3 = (3/4)[(V - V_qs)^2 + m^2 (Y - Y_qs)^2] a^3, averaged over k periods;
(iii) <rho_V> against S^2/(12 m^2) and KV's kick formula; (iv) C and C a^4. Fortran DOP853 (scipy `ode`), rtol 1e-12.
Oscillation-averaged quantities only: a snapshot of the oscillating rho_V is not used.

### 2.2 Reduction error versus m_q/H (adiabatic transition, w = 1, x = pi r/2 >= 16)

| m/H_c | max eps0 (window) | R1 correction predicted | max eps1 (window) | max eps1 (pointwise) | \|C\|/scale | C a^4 drift |
|---|---|---|---|---|---|---|
| 10.1 | 1.7e-4 | 1.7e-4 | 6.5e-6 | 1.4e-4 (osc.) | 1.5e-14 | 3e-15 |
| 30.3 | 2.5e-5 | 2.4e-5 | 1.2e-6 | 2.4e-5 (osc.) | 1.1e-14 | 3e-15 |
| 101 | 2.83e-6 | 2.77e-6 | 7.0e-8 | 3.2e-6 (osc.) | 4.7e-15 | 3e-14 |
| 303 | 3.129e-7 | 3.131e-7 | 2.0e-10 | 1.6e-8 | 4.2e-15 | 1e-14 |
| 1010 | 2.743e-8 | 2.743e-8 | 1.2e-11 | 3.9e-11 | 4.6e-15 | 2e-14 |
| 3029 | 3.048e-9 | 3.048e-9 | 2.0e-11 (floor) | 2.1e-13 | 8.0e-15 | 2e-14 |
| 10098 | 2.745e-10 | 2.765e-10 | 8.9e-11 (floor) | 1.5e-13 | 8.6e-14 | 5e-13 |

- **eps0 = 0.0298 (H_c/m)^(2.007)** (log-log fit, 5 points, m/H = 1e2...1e4). The error of "GR + V_m(t) - c*(mu)" is
  exactly the O(H^2/m^2) adiabatic q-correction: measured and predicted agree to 1e-3 relative for m/H >= 300.
- **R1 removes it.** Residual eps1 <= 7e-8 at m/H = 100 (2.5% of eps0, limited by leakage of the initial-data
  oscillation through the window), 2e-10 at 300 (6e-4 of eps0), and the floor 1e-13 (pointwise) to 1e-11 (window, a
  J-difference round-off) for m/H >= 1000. Consistent with eps1 = O((H/m)^4) (<= 7 (H/m)^4); a clean slope is not
  measurable because the floor is reached at m/H ~ 1000. eps2 = eps1: no oscillation energy is deposited at w = 1.
- At m/H = 10-30 the residual is set by the non-adiabatic start (m/H_i = 3.7 and 9): R1's error is 4% (10) and 5%
  (30) of eps0. With w = 0.3 at m/H = 10 (x = 4.7) a measurable oscillation is deposited and eps1 = 3.6e-5.
- **Window and tolerance checks.** k = 1 vs 3 periods (m/H = 300): eps0 3.129e-7 vs 3.071e-7 (the window bias of the
  smooth part), eps1 2.0e-10 vs 2.1e-10. rtol 1e-9...1e-13 (m/H = 100): eps0 identical to 1e-8 relative, eps1 to 1e-6
  relative; RHS evaluations 2.7e4 -> 9.0e4.
- **Limit.** m/H >= 3e4 was not integrated: in these variables dV/dt is a difference of O(1) terms (S/3 - m^2 F(Y))
  that cancel to O(H^2/m^2), so round-off exceeds the requested accuracy and DOP853 stalls (at m/H = 3000 already a
  rtol of 1e-13 multiplies the step count by 8 with unchanged results). A deviation-variable formulation would remove
  this; it is not needed for any conclusion below.

### 2.3 Deposited oscillation energy versus transition duration

Formula (equations §5.3): E_osc a^3 = (4/3) Delta^2/(K m^2) a_c^3 (x/sinh x)^2, x = pi m tau/2, tau = w/H_c.

| x | m/H = 30 | m/H = 100 | m/H = 300 | m/H = 1000 | formula (x/sinh x)^2 |
|---|---|---|---|---|---|
| 0.1 | 1.0052 | 1.0000 | 1.0000 | 1.0000 | 0.99668 |
| 1.0 | 1.0054 | 0.9997 | 0.9999 | 1.0000 | 0.72409 |
| 2.0 | 1.0101 | 0.9998 | 0.9999 | 1.0000 | 0.30408 |
| 4.0 | 1.111 | 1.0072 | 1.0013 | 1.0002 | 0.02148 |
| 6.0 | 1.516 | 1.0246 | 1.0049 | 1.0009 | 8.84e-4 |
| 8.0 | 3.349 | 1.0277 | 1.0117 | 1.0020 | 2.87e-5 |

(entries: measured / formula; full table with 11 x-values per mass in the JSON)

- The formula is confirmed. Its error is the variation of H during the transition, O(x H/m) at large x: it vanishes
  as m/H grows (ratio - 1 at x = 8: 2.35, 0.028, 0.012, 0.002).
- **Scaling fit** for x >= 2: ln(E/x^2) = ln A - kappa x gives kappa = 1.999, 2.003, 2.004 (m/H = 100, 300, 1000);
  the same fit to the formula's own values gives kappa = 2.005, A = 4.116 (asymptotically kappa = 2, A = 4:
  E ~ 4 x^2 e^(-pi m tau)). At m/H = 30: kappa = 1.81 (finite-H regime).
- **Delta^2 scaling**: measured/formula = 0.9989, 0.9997, 0.9985 for Delta/3H_c^2 = 0.0125, 0.05, 0.2 (m/H = 100, x = 1).
- **Floor** (Delta = 0 controls): E/E_sudden = 5.8e-6, 2.5e-8, 4.4e-12, 2.9e-16 at m/H = 30, 100, 300, 1000; always
  below the smallest deposit measured. Plateau spread over late windows <= 1.5e-4.
- **Tolerance**: E at x = 2, m/H = 100 changes by 1e-9 relative from rtol 1e-9 to 1e-13.
- At the physical hierarchy (electroweak: m_q/H = 5.6e31, tau ~ 0.1/H) x ~ 1e31, so (x/sinh x)^2 = exp(-1e31): no
  oscillation energy is deposited by any SM transition [C].

### 2.4 rho_V and the KV kick formula

Window-averaged <rho_V> / [S^2/(12 K0 m^2)] = 0.9971-0.9990 (m/H = 10), 0.99993-1.00007 (100), 1.000000 (>= 300),
in the dust-dominated part and through the transition. Against KV's (9/8) kappa^2 H^4 K'^2/eps'' with kappa = S/rho_m:
0.9964-0.9987 (10), 0.9994-1.0001 (100), 1.0000 (>= 1000). Both formulas describe the same quasi-static rho_V and are
confirmed to O(H^2/m^2). Neither is a constant: rho_V follows S^2, i.e. a^-6 in a dust era.

### 2.5 Error budget of the reduction at realistic scales [C/E]

From §2.2: eps0 = 0.03 (H/m)^2 x (S/rho)-dependence (here S/rho ~ 0.2), eps1 <= 7 (H/m)^4. At natural m_q:
(H/m)^2 <= 5e-61 at T <= 1 TeV, so R0 is exact to 1e-62 and R1 to 1e-121; the R1 correction itself is computed in §4.
At m_q = 4 meV, (H/m)^2 = 1e-4 at the EW crossover: R1 needed, error <= 1e-7 for T <~ 100 GeV.

---------------------------------------------------------------------------------------------------------------------

## 3. Early relaxation: Planck-scale displaced vacuum, mu = mu0 [C] — `code/s3_early_relaxation.py`, `results/early_relaxation.json`

Natural scales (M = 1, chi = 1, m = 0.8165 M); full nonlinear equations.

### 3.1 KV's dust-like regime reproduced (Gamma = 0, QV1 proper)

| KV case | t<H> (KV: 2/3) | t^2 <rho_V>/M^2 | KV (2/3) d_M^2 | d_M |
|---|---|---|---|---|
| w_M = 1/3, r_M = 1/20 (Figs 1-2) | 0.6655 | 0.6620 | 0.6667 | 1 |
| w_M = 0, r_M = 1/200 (Figs 3-4) | 0.6670 | 0.6497 | 0.6495 | 0.987 |
| w_M = 0, r_M = 1/2 (Fig 5) | 0.6667 | 0.2335 | 0.2332 | 0.591 |

Window tau = 1800-3000. **<rho_V> = (2/3) d_M^2 M^2/t^2** [Q KV (5.15)-(5.17); D literature §4.4] is reproduced
(0.7% residual for w_M = 1/3 is the expected O(1/tau) term). The same numbers were obtained by the literature track
from KV's own form of the equations (5.2a,b); this track's independent Y-variable implementation agrees to 4 digits,
a cross-check of the code. <rho_V> is half the critical density (ratio 0.496): the oscillating q-sector acts as dust
(a ~ t^(2/3)) for any initial matter. **Without a decay channel QV1 stays dust-dominated from the Planck time on**
(KV's own caveat: conflicts with BBN, literature §4.4).

### 3.2 With scalaron decay (QV1+Gamma, effective source [Q/E])

Gamma = 4 m^3/(192 pi M^2) = 3.61e-3 M = 4.42e-3 m. Integrated to t = 14/Gamma.

Displaced q at rest (rho_M,i = 0; the only initial energy is the vacuum energy rho_V(y_i)):

| y_i | rho_V,i/M^4 | H_i/m | f_inj | f_com | T_rh (GeV) |
|---|---|---|---|---|---|
| -0.90 | 0.0284 | 0.376 | 155 | 8.0e3 | 5.18e16 |
| -0.75 | 0.0527 | 0.325 | 8.92 | 220 | 5.18e16 |
| -0.50 | 0.0521 | 0.228 | 2.04 | 31.1 | 5.15e16 |
| -0.25 | 0.0215 | 0.120 | 1.22 | 12.0 | 5.12e16 |
| -0.10 | 0.00435 | 0.049 | 1.062 | 6.04 | 5.05e16 |
| -0.03 | 4.3e-4 | 0.015 | 1.016 | 2.99 | 4.82e16 |
| +0.03 | 4.7e-4 | 0.015 | 0.986 | 2.95 | 4.83e16 |
| +0.10 | 0.00568 | 0.051 | 0.959 | 5.79 | 5.06e16 |
| +0.25 | 0.0423 | 0.130 | 0.917 | 10.3 | 5.12e16 |
| +0.50 | 0.219 | 0.270 | 0.878 | 16.5 | 5.15e16 |
| +1.00 | 1.333 | 0.577 | 0.847 | 28.5 | 5.16e16 |

- f_inj = ∫Q a^3 dt / (rho_V,i a_i^3): the energy handed to radiation, counted at injection, per comoving volume.
  **For small displacement (\|y_i\| <= 0.1, H_i << m) the whole initial vacuum energy becomes radiation**:
  f_inj = 0.96-1.06. For positive y_i up to 1 a part (8-15%) is lost to Hubble friction during the first oscillation.
  For y_i -> -1 the Jordan-frame accounting is not conserved (K = (1+Y) K0 -> 0, G -> infinity; the Einstein-frame
  energy rho_V/(1+Y)^2 is 4-100 times larger): f_inj > 1 there is a strong-coupling artefact of the frame, outside any
  control.
- f_com = rho_r a^4 (end)/(rho_V,i a_i^4): comoving radiation; it exceeds 1 because the dust-like phase lasts until
  t ~ 1/Gamma. All final radiation comes from the vacuum (no initial matter).
- **T_rh is set by Gamma, not by the displacement**: (4.8-5.2)e16 GeV for every y_i (g_* = 106.75), from rho_r at
  t = 1/Gamma; Gorbunov-Panin's formula [Q] gives 5.03e16. T_rh scales as Gamma^(1/2): 1.63e16 (Gamma/10), 1.63e17
  (10 Gamma) against 1.59e16, 1.59e17 from the formula.
- The oscillation energy decays at Gamma_fit/Gamma = 1.0000000-1.0000007; residual E_osc/rho_r = 4e-6 at
  t = 14/Gamma (e^-14); late t H = 0.51 -> 1/2 (radiation era). |C|/3H^2 <= 2e-14, C a^4 drift <= 1e-11.
- With KV's boundary data (h_i = 1/2 imposed, so C = 0 forces a large initial qdot when rho_M + rho_V < 3H_i^2), the
  q-kinetic energy dominates for small y_i and f_inj = 0.8-98 is not a measure of the vacuum energy; T_rh is again
  (4.9-5.2)e16 GeV and 90-99.9% of the final radiation comes from the q-sector (the rest is KV's initial r_M).

### 3.3 What §3 shows [C/A]

The relaxation of a Planck-scale vacuum energy at mu = mu0 is dynamical (KV) and, once the standard decay channel
of the scalaron is included, ends in a radiation-dominated universe at T_rh ~ 5e16 GeV with the q-field at q0 and
no oscillation left (e^(-Gamma t)). The fraction of the initial vacuum energy that becomes radiation is ~1 for
small displacements and O(1) for large ones; the reheating temperature does not depend on it. This supplies the
initial data of §4 (q on the quasi-static track, m_q/H >= 570 from T_rh on).

### 3.4 EFT-validity caveat [A]

The initial state has rho ~ 0.05 M^4 and H ~ 0.2-0.6 m ~ 0.3 M: curvature R ~ M^2. The classical action with two
functions K(q), eps(q) is not a controlled description there: higher-curvature terms (R^2, Weyl^2) with O(1)
coefficients, quantum-gravity corrections and the q-dependence of loop corrections are all O(1). The decay rate is a
perturbative formula (valid for Gamma << m, H << m, \|Y\| << 1, and particles that do not backreact); at \|Y\| ~ 0.5
non-perturbative production (parametric resonance) may dominate. The numbers of §3 illustrate the mechanism and its
scaling; they are not predictions.

---------------------------------------------------------------------------------------------------------------------

## 4. Realistic history from T = 1 TeV to the far future [C] — `code/s4_realistic_history.py`, `results/realistic_history.json`

GR + SM plasma (SS EoS, p(0) = 0) + dust, with the R1 q-correction; t from entropy conservation (t = 1/(2H) at
1 TeV, error < 1e-3). Today is defined by T0 = 2.7255 K.

### 4.1 Natural m_q = 1.988e18 GeV (QV1 proper)

| epoch | T (GeV) | t (s) | m/H | H/H_GR - 1 | rho_V(q) (GeV^4) | delta G/G | Gdot/G (1/yr) |
|---|---|---|---|---|---|---|---|
| T = 1 TeV | 1e3 | 2.4e-13 | 1.4e30 | 2.6e-64 | 6.9e-55 | -2.0e-64 | 4.8e-44 |
| EW crossover | 160 | 9.3e-12 | 5.6e31 | 1.4e-65 | 6.8e-58 | -6.2e-66 | 5.8e-47 |
| T = 100 GeV | 100 | 2.4e-11 | 1.5e32 | 1.8e-66 | 4.3e-59 | -1.6e-66 | 3.4e-48 |
| QCD crossover | 0.155 | 2.0e-5 | 1.1e38 | 4.1e-77 | 7.6e-81 | -2.1e-77 | 8.5e-65 |
| T = 10 MeV | 1e-2 | 7.4e-3 | 4.5e40 | 3.7e-84 | 3.2e-96 | -4.2e-85 | 1.7e-74 |
| e+e- annihilation | 2e-4 | 22.7 | 1.4e44 | 6.5e-90 | 4.7e-106 | -5.2e-90 | 1.3e-83 |
| BBN (70 keV) | 7e-5 | 264 | 1.6e45 | 3.2e-92 | 1.4e-111 | -8.9e-93 | 4.4e-87 |
| matter-radiation equality (z = 3403) | 8.0e-10 | 1.6e12 | 8.8e54 | 6.9e-111 | 7.4e-148 | -6.5e-111 | 2.1e-115 |
| recombination (z = 1090) | 2.6e-10 | 1.2e13 | 6.0e55 | 2.3e-112 | 8.0e-151 | -2.1e-112 | 1.0e-117 |
| today | 2.35e-13 | 5.45e17 (17.3 Gyr) | 2.5e60 | 1.8e-121 | 4.7e-169 | -1.6e-121 | 1.9e-131 |
| a = 10 a0 | | 1.7e19 | 7.8e61 | 1.8e-124 | 4.7e-175 | -1.6e-124 | 6.0e-136 |
| a = 1000 a0 | | 1.7e22 | 7.8e64 | 1.9e-130 | 4.7e-187 | -1.6e-130 | 6.0e-145 |

- **rho_V(t) of the q-sector** = S^2/(12 K0 m^2): it follows the square of the matter trace (SM trace anomaly,
  S/rho = 4e-4 at 1 TeV, 1.9% at the EW crossover, 27% at the QCD crossover, 8e-4 at 10 MeV, then e+- and dust).
  It is never constant and never larger than 1e-54 GeV^4 = 3e-8 rho_Lambda after 1 TeV.
- **The q-correction** delta rho_q = [S^2/(12K0) - H^2(S + S_N)]/m^2: its G-variation part (-H^2(S+S_N)/m^2)
  dominates rho_V by 4 rho/S (in a dust era delta rho_q = (3/4) rho_d^2/(K0 m^2) = 9 rho_V). H(T) equals the GR
  history to 1e-63 or better at every epoch.
- **G_eff** = 1/(8 pi K(q)): delta G/G = -S/(3K0 m^2), Gdot/G = -H S_N/(3K0 m^2): 1.9e-131/yr today (LLR bound
  ~1e-14/yr: met trivially, ledger row 1).
- **Transitions** (EW, QCD, e+e-): the plateau energy sits in rho - 3p and gravitates as in GR until released; m/H
  >= 1e30 so the q-mode follows adiabatically and deposits nothing ((x/sinh x)^2 ~ exp(-1e31)).
- **Asymptotic state**: no vacuum term; matter domination forever, t H -> 2/3 (0.66667 at a = 1000 a0), H -> 0,
  rho_V(q) ~ a^-6 -> 0: the end state is Minkowski, as P_q demands. Under QV1-cal: de Sitter, H -> H_Lambda
  (H/H_Lambda = 1.0000000002 at a = 1000 a0).
- Today in QV1: H0 = 37.72 km/s/Mpc and t0 = 17.28 Gyr (GR + no vacuum term with the CMB physical densities); see §5.1.

### 4.2 Bracket m_q = 4 meV

Reduction (valid where m/H >~ 300): H/H_GR - 1 = 3.5e-6 at the EW crossover (m/H = 112), 4.6e-7 at 100 GeV, 7e-12 at
10 GeV, 1e-17 at the QCD crossover; delta G/G = -1.5e-6 (EW), -5e-12 (10 GeV), -4e-62 today; rho_V(q) = 170 GeV^4 at
the EW crossover (7e-9 of rho), 1e-109 GeV^4 today. At 1 TeV m/H = 2.9: the reduction does not apply there.

Full equations (SM plasma table, dust), integrated from T_i to m/H = 1e5 (T ~ 5 GeV):

| start | Omega_q h^2 (today, stable oscillation) | / Omega_c h^2 |
|---|---|---|
| quasi-static track at 1 TeV (m/H = 2.9) | 378 | 3.2e3 |
| q = q0 at rest at 1 TeV | 502 | 4.2e3 |
| q = q0 at rest at 10 TeV | 4.1e4 | 3.4e5 |
| q = q0 at rest at 100 TeV | 1.9e5 | 1.6e6 |
| q = q0 at rest at 1000 TeV | 4.1e5 | 3.4e6 |
| 100 TeV, S/rho = 0 above 300 GeV | 115 | 960 |
| 100 TeV, S/rho + 1 sigma above 300 GeV | 2.1e6 | 1.7e7 |

- Mechanism [C/D]: above T ~ 1 TeV the q-mode is lighter than H and is pushed by the SM trace anomaly,
  Ydd + 3H Yd = S/3 => Y ~ (S/rho) x (e-folds); at m ~ H (T_osc = 0.98 TeV) it starts to oscillate about the
  quasi-static track with amplitude \|delta Y\| = \|delta G/G\| ~ 3e-5 (adiabatic start) to 7e-4 (q0 at 100 TeV), and
  the oscillation redshifts as dust. Its ratio to radiation grows as a from T_osc ~ 1 TeV, so an oscillation carrying
  2e-9 of the energy density at 1 TeV already gives Omega_q h^2 = 378 today.
- The oscillation energy is constant (E a^3 spread <= 1e-8) from m/H ~ 300 down to the end; tolerance check (rtol
  1e-12) and a longer run (to m/H = 3e5) give the same Omega h^2 (378.4).
- Decay is negligible: Gamma (N_s = 4, an upper bound) / H(1 MeV) = 1.6e-49.
- **Mass scan** (same full equations, start on the quasi-static track at m/H = 3 / q = q0 at m/H = 0.01): Omega_q h^2 =
  2.2e3 / 6.2e5 (1 eV), 4.3e4 / 1.2e7 (1 keV), 8.8e5 / 2.6e8 (1 MeV), 1.9e7 / 5.6e9 (1 GeV); roughly ~ m_q^(0.3-0.45).
- The reduction holds where it should: max \|eps2\| (R1 + measured dust) = 2.6e-6 for m/H = 100-3e4 in the adiabatic
  4 meV run, limited by window leakage of the relic oscillation itself; \|C\|/scale <= 4e-13 (<= 1e-11 in the
  misaligned runs); E a^3 constant to <= 4e-9 over the late windows in every run.
- **Verdict (retrospective)**: unless the oscillation amplitude at the onset (T ~ 1 TeV) is below \|delta G/G\| ~ 5e-7
  (1% of the trace-driven displacement there), an unscreened 4 meV q-mode overcloses the universe by 3e3-3e6
  (Omega_c h^2 = 0.120, ledger row 13); the SM trace anomaly alone produces more. The window 4 meV < m_q < (decay before BBN) is excluded for generic initial data: Gamma(N_s = 4) =
  H(T = 4 MeV) needs m_q >~ 2e5 GeV [E]. This is the cosmological moduli problem for the q-mode (Coughlan et al. 1983
  [R]; Banks-Kaplan-Nelson hep-ph/9308292 and de Carlos-Casas-Quevedo-Roulet hep-ph/9308325 [Q abstract]); R^2
  scalaron dark matter (Cembranos 0809.1653 [Q abstract]) is the related literature. It strengthens the ledger's
  lower bound (~4 meV from Eot-Wash) for the light branch; QV1 at natural m_q is unaffected.

---------------------------------------------------------------------------------------------------------------------

## 5. Variants and comparison (all retrospective) [C] — `code/s5_variants_data.py`, `results/variants_comparison.json`

Data used (downloaded to the scratchpad): DESI DR2 BAO Gaussian likelihood (13 numbers + covariance; arXiv:2503.14738);
Pantheon+ (1590 SNe with z_HD > 0.01, STAT+SYS covariance, absolute magnitude marginalised; arXiv:2202.04077); SH0ES
H0 = 73.04 +- 1.04 and q0 = -0.51 +- 0.024 (2112.04510 abstract); Planck 100 theta_* = 1.0411 +- 0.0003 (1807.06209
abstract); globular-cluster age of the universe 13.5 +0.16/-0.14 (stat) +- 0.5 (sys) Gyr (2007.06594 abstract).
Code checks: flat-LambdaCDM fits reproduce the published SN-only Omega_m = 0.3316 +- 0.018 (published 0.334 +- 0.018)
and BAO-only Omega_m = 0.2975 +- 0.0086 (published 0.2975 +- 0.0086). r_d from the Aubourg et al. 2015 fitting formula
[R] (147.04 Mpc). theta_* computed with the Hu-Sugiyama z_* fit [R] is 0.27% low for Planck LambdaCDM (1.0383); the
calibrated values below divide out this code systematic.

### 5.1 (a) QV1 (P_q): no late vacuum term — **retrospective failure**

| quantity | QV1 | data / best fit | tension |
|---|---|---|---|
| H0 (from T0, omega_b, omega_c, flatness) | 37.74 km/s/Mpc | SH0ES 73.04 +- 1.04 | -34 sigma |
| q0 | +0.500 | SH0ES+Pantheon+ -0.51 +- 0.024 | +42 sigma |
| SN shape (Pantheon+) | chi^2 = 2101.4 | flat LCDM 1402.9 (Omega_m = 0.332) | Delta chi^2 = 699 |
| BAO shape (DESI DR2, H0 r_d free) | chi^2 = 1421.9 | flat LCDM 10.3 (Omega_m = 0.2975) | Delta chi^2 = 1412 |
| BAO absolute (r_d = 147.0 Mpc) | chi^2 = 15333 / 13 | | |
| 100 theta_* (calibrated) | 0.943 | 1.0411 +- 0.0003 | -327 sigma |
| age | 17.27 Gyr | t_U = 13.5 +- 0.5 (lower bound) | passes |
| d_L (z = 0.5, 1, 2) | 4373, 9306, 20141 Mpc | QV1-cal 2930, 6810, 15958 | |

QV1 is an Einstein-de Sitter universe (Omega_m = 0.9997) at late times, i.e. the Omega_m = 1 point of the flat
LambdaCDM family; the SN and BAO shapes exclude it at ~26 and ~38 sigma (sqrt Delta chi^2). Forcing H0 = 73.04 instead
requires omega_m = 0.533 (3.7 x Planck) and gives an age of 8.92 Gyr, 8.8 sigma below the globular-cluster age. The
theta_* comparison uses LambdaCDM-fitted omega's and is a consistency statement, not a CMB fit. **QV1 fails
retrospectively; the failure is not repaired inside QV1.** Any repair is a new labelled version (QV1-cal or QV2).

### 5.2 (b) QV1-cal: reproduces flat LambdaCDM; the calibration fixes the residual

mu = mu0[eps_inf] - rho_Lambda/q0 (c = rho_Lambda = 2.52e-47 GeV^4, relative size of the offset 7.2e-121, calibrated
to the observed value): H0 = 67.36, Omega_m = 0.314, age 13.81 Gyr, q0 = -0.529, 100 theta_* = 1.0410 (calibrated,
-0.2 sigma). chi^2: SN 1403.88 (Planck-LCDM 1403.89; Delta = 0.96 from the SN best fit), BAO shape 13.80 (Planck-LCDM
13.79), BAO absolute 33.4 (Planck-LCDM 33.0; the known 2.3 sigma DESI-CMB tension), H0 vs SH0ES -5.5 sigma (the
known Hubble tension). Late attractor: de Sitter with rho_V = c + c^2/(K m^2) + ... (type ii offset; relative
correction 1.4e-120 at natural m_q). Confirmed at m/H = 100 with the full equations: rho_V,late = 0.03000012 against
the static root 0.03000012 (fixed-mu case of §5.4, same algebra). **The residual is put in by the calibration; nothing
in QV1 produces it.** DESI's w0wa preference (ledger row 8) is a retrospective tension shared with LambdaCDM.

### 5.3 (c) The rejected readings

**Initial-equilibrium mu** (mu = mu0 of the vacuum before the electroweak step): rho_V,inf = -Delta_EW =
-1.19e8 GeV^4 = -4.7e54 rho_Lambda. With the SM EoS from 1 TeV: H = 0 at **T = 44.6 GeV, t_turn = 2.02e-10 s**
(radiation-era analytic (pi/4) sqrt(3K/Delta) = 2.00e-10 s); recollapse to a crunch at **t ~ 4.0e-10 s** (time-
symmetric, entropy-conserving recollapse [A]). The universe never reaches the QCD crossover; adding the QCD bag
changes nothing (turnaround precedes it). Excluded (retrospective: we exist after T = 44.6 GeV).

**Fixed-number mu plus a heavy-threshold constant c** (eps -> eps + c at fixed mu, type i): exact static roots
(mpmath, 300 digits) give rho_V,inf = c [1 + (4/3) c/(K0 m^2) + ...]: rho_V,inf/c - 1 = 1.4e-120 (c = rho_Lambda),
6.8e-66 (c = Delta_EW), 1.35e-60 (c = 2.37e13 GeV^4). **The history moves one-for-one**:
- c = +rho_Lambda is QV1-cal (H0 67.36, age 13.81 Gyr): a fixed mu is a free Lambda.
- c = -rho_Lambda: turnaround when rho_m = rho_Lambda at a = 0.770 a0, photon temperature 3.54 K (2.7255 K is never
  reached); t_turn = 18.4 Gyr, crunch at 36.7 Gyr.
- 10 TeV scalar threshold, c = -2.37e13 GeV^4: turnaround at T = 912 GeV, t = 4.0e-13 s, crunch ~8e-13 s; with the
  opposite sign, de Sitter domination below T = 912 GeV with H = 1.16e-12 GeV (no electroweak plasma era).

**P_q with the same c re-absorbed into mu0[eps_inf + c]**: rho_V,inf = 0 exactly for every c (and every order). With
the KV quartic the new equilibrium is exact: (q0'/q0)^2 = (1 + sqrt(1 + 8 chi c))/2 [D]; delta K/K = chi c (1 + ...):
3.4e-66 (Delta_EW), 6.8e-61 (10 TeV), 7.2e-121 (rho_Lambda); m_q^2/K0 shifts by 4 chi c (m_q^2 itself by 5 chi c). This is a renormalisation of
the measured G, not an observable change. Existence requires c >= -M^4/8 (stable branch). **The history does not
move.** The tuning is carried by the state (dmu0/dc = 1/q0), not removed (equations §6.3).

### 5.4 (d) Full-equation confirmation at m/H = 100 (toy matter, KV shape)

| reading | full equations | analytic |
|---|---|---|
| P_q | late rho_V,tot/3H^2 = 3.4e-10 (Minkowski attractor; quasi-static S^2/12m^2 only) | 0 + O(S^2/m^2) |
| initial-equilibrium (mu = mu0[eps_inf + Delta], exact KV) | turnaround t = 3.51 H_c^-1; rho_V at turnaround -0.1499971 | static AdS root -0.1499970; series -Delta + (4/3)Delta^2/m^2 = -0.149997 |
| fixed mu + c (c = +0.2 Delta) | rho_V,late = 0.03000012, H^2 = 0.01000000 | root 0.03000012, 0.0100000; rho_V/c = 1.000004 = 1 + (4/3)c/m^2 |
| fixed mu - c | turnaround t = 7.92 H_c^-1, crunch | AdS |
| P_q, c re-absorbed | binned H in own-G units differs by 5.6e-10; late rho_V,tot/3H^2 = 3.4e-10 | physical shift = G renormalisation only, delta H/H = -(chi c)/2 = -1.0e-6 |

---------------------------------------------------------------------------------------------------------------------

## 6. What the solved example shows

**Emerges from the equations** (no input beyond P_q and the stated functions):
1. **Cancellation of constants under P_q** [D][C]: every q-independent constant (heavy thresholds, loops, the SM
   T = 0 vacuum) is absorbed by the state; rho_V,inf = 0 exactly; the only trace is a renormalisation of G by
   chi c <= 1e-60 (§5.3, §5.4). With a fixed number mu the same constants gravitate one-for-one; with an initial-state
   mu the electroweak step alone crunches the universe at 4e-10 s.
2. **Adiabatic tracking through transitions** [C]: the full q dynamics equals GR + V_m(t) - c*(mu) up to
   0.03 (H/m)^2, and R1 removes that to O((H/m)^4) (§2.2). Transitions deposit (4/3)Delta^2/(Km^2)(x/sinh x)^2 of dust-
   like oscillation energy, confirmed to 1e-4 (§2.3); at the physical hierarchy this is exp(-1e31).
3. **No constant remnant** [C]: second-order terms appear only as rho_V = S^2/(12K0m^2) (follows the trace, a^-6 in a
   dust era) and as decaying oscillations; at natural m_q the q-sector never exceeds 1e-54 GeV^4 after 1 TeV and is
   4.7e-169 GeV^4 today (§4.1).
4. **Relaxation of a Planck-scale vacuum energy** [C, with the effective decay]: dust-like (KV) then radiation, T_rh ~
   5e16 GeV independent of the displacement (§3).

**Does not emerge**:
1. **The observed residual** rho_Lambda: QV1 gives exactly zero and fails the late-universe data at 26-327 sigma
   (§5.1); QV1-cal matches LambdaCDM because rho_Lambda is inserted by calibration (§5.2). No relation between
   measurable quantities beyond the calibration target follows from QV1 at this level.
2. **Dark matter**: at natural m_q the q-oscillation decays at the Planck era; a light q-mode overproduces it by
   >= 1e3 (§4.2). The CDM density is an input (omega_c).
3. **A radiation era without the decay source**: classical QV1 alone stays dust-like from the Planck time (§3.1).

**Error budget** (largest first; none changes a verdict):

| source | size | effect |
|---|---|---|
| EFT at Planck densities (§3) | O(1) | §3 numbers illustrative only |
| SM trace anomaly above 300 GeV (SS error +-1.5e-3 abs) | O(100%) of S/rho there | 4 meV relic: factor 0.002-10 (115 to 2e6); natural m_q: none |
| initial data of the light q (not fixed by QV1) | factor 1e3 in Omega_q | §4.2 verdict needs only Omega_q > 0.12 |
| decay rate formula (factor ~2, GP footnote) and N_s (xi) | T_rh ~ Gamma^(1/2) | +-40% in T_rh |
| omega_b, omega_c from a LambdaCDM fit; massless nu | ~1% (0.4% nu) | QV1 H0 37.7 +- 0.3 [E] |
| theta_* code systematic | 0.27% (calibrated out) | -327 sigma -> -336 sigma raw |
| r_d formula | 0.02% [R] | BAO absolute chi^2 |
| e+- trace below 10 MeV | 16% at the junction | q-correction at BBN (1e-92) |
| R0/R1 reduction | 0.03 (H/m)^2 / <= 7 (H/m)^4 | 1e-62 / 1e-121 at natural m_q |
| Delta_EW loop corrections | O(10%) | crunch time +-5% |
| numerical integration | rtol 1e-11..1e-12; C <= 1e-13 | none |

**GRUT's contribution in this track**: none of the physics; the scenario uses only KV's construction, standard SM
thermodynamics and standard scalaron decay. P_q is used as stated in the ledger (§1.4 there).

**Feeds for the gate** (for the relation agent): Q2 — QV1 fixes the residual to exactly zero (sign and scale: zero),
which fails the exposed late-universe data; QV1-cal fixes it by calibration. Q3 — stable under heavy thresholds
without retuning any Lagrangian parameter, because the state is re-selected (§5.3). Q4 — no independent consequence
at natural m_q (all q-sector effects <= 1e-54 GeV^4 after 1 TeV); the light branch is excluded retrospectively by
Omega_c h^2.

---------------------------------------------------------------------------------------------------------------------

## 7. Reproduction, files, limits

**Reproduction** (from this directory; downloads go to the scratchpad, never into the repository):

    python3 -I code/run_all.py        # ~10 min on 4 cores: fetch inputs, then s1..s5

| file | content | NOTES |
|---|---|---|
| `code/qv1_core.py` | model, full FRW system, reductions, quasi-static track, integrator wrapper (Fortran DOP853) | §1 |
| `code/sm_eos.py` | Saikawa-Shirai EoS loader, trace S(T), e+- trace below 10 MeV | §1.5 |
| `code/s1_symbolic_checks.py` -> `results/symbolic_checks.json` | Bianchi propagation (with/without dissipation), Y-variable reduction, KV shape, R1 formula | §1.4 |
| `code/s2_reduction_validation.py` -> `results/reduction_validation.json` | §2 (all tables, fits, controls, convergence, series) | §2 |
| `code/s3_early_relaxation.py` -> `results/early_relaxation.json` | §3 | §3 |
| `code/s4_realistic_history.py` -> `results/realistic_history.json` | §4 (epoch tables, 220-point series for figures, 4 meV and mass-scan full runs) | §4 |
| `code/s5_variants_data.py` -> `results/variants_comparison.json` | §5 | §5 |
| `code/run_all.py` | reproduction driver | §7 |

Figures are provided as data series in the JSON files (`series` keys): eps0/eps1/eps2 and rho_V/S^2 versus N
(m/H = 10, 100, 1000), E_osc a^3 versus N (x = 2), the decay run (E a^3, rho_r a^4, tH versus t), H(T), rho_V(T),
delta G/G(T), Gdot/G(T) for both masses.

**Limits.**
- Classical field equations plus one effective decay source; no loop corrections to K(q), eps(q) beyond the
  constant c; no inhomogeneities (bubbles, walls, q-perturbations); membrane nucleation outside the action.
- The KV quartic/linear forms are one concrete choice; at realistic scales only m_q enters (Y <= 1e-64), the shape
  coefficients b3, b4 matter only in §3 (O(1) amplitudes).
- §2 stops at m/H = 1e4 (round-off cancellation, §2.2); the realistic-scale claims use the analytic reduction whose
  error law is validated over 1e2-1e4.
- §3: Jordan-frame energy accounting at large \|Y\| is frame dependent; EFT caveat §3.4.
- §4.2 relic depends on pre-TeV initial data and on the poorly known high-T trace; the verdict needs only a factor
  > 1 and holds in every variant computed, but the precise excluded mass window above 1 GeV is an extrapolation [E].
- §5: comparisons use LambdaCDM-fitted early-universe parameters (consistency, not a refit of QV1 to the CMB); SN and
  BAO shape tests do not depend on them.
- Read levels below; nothing in §5 is held out — all retrospective.

## 8. Sources and read levels (this track)

| source | identifier | read level | used for |
|---|---|---|---|
| Klinkhamer-Volovik q-theory papers | 0711.3170, 0806.2805, 0905.1919, 0907.4887, 1102.3152, 1601.04676 | via literature track (full/section there); not re-read | §1, §3 |
| Saikawa, Shirai, JCAP 05 (2018) 035 | arXiv:1803.01038 | ancillary data tables (full); TeX §3 intro, §3.3 (normalisation p(0) = 0), §3.4 (full EoS, uncertainty) | §1.5, §4 |
| Gorbunov, Panin, PLB 700, 157 (2011) | arXiv:1009.2448 | section (TeX: decay rates eqs (T4+), (T4++), T_reh formula (T5+), footnote on Starobinsky/Vilenkin) | §1.6, §3 |
| Starobinsky 1980; Vilenkin PRD 32, 2511 (1985) | | recalled via Gorbunov-Panin [R] | §1.6 |
| DESI DR2 II | arXiv:2503.14738 (PRD 112, 083515) | abstract + released Gaussian BAO likelihood (CobayaSampler/bao_data desi_bao_dr2) | §5 |
| Pantheon+ (Brout et al., ApJ 938, 110 (2022); Scolnic et al.) | arXiv:2202.04077 | abstract + data release (distances, STAT+SYS covariance) | §5 |
| SH0ES, Riess et al. | arXiv:2112.04510 | abstract (H0 = 73.04 +- 1.04; q0 = -0.51 +- 0.024) | §5 |
| Planck 2018 VI | arXiv:1807.06209 (A&A 641, A6) | abstract (omega_b, omega_c, 100 theta_*) | §1.8, §5 |
| Fixsen 2009 | arXiv:0911.1955 (ApJ 707, 916) | abstract (T0) | §1.8 |
| Valcin et al. 2020 | arXiv:2007.06594 (JCAP 12, 002) | abstract (t_U) | §5 |
| Bond et al. 2013 (HD 140283, 14.46 +- 0.8 Gyr) | arXiv:1302.3180 | abstract | not used for a verdict |
| Banks, Kaplan, Nelson | hep-ph/9308292 (PRD 49, 779) | abstract | §4.2 credit |
| de Carlos, Casas, Quevedo, Roulet | hep-ph/9308325 (PLB 318, 447) | abstract | §4.2 credit |
| Coughlan, Fischler, Kolb, Raby, Ross 1983 | PLB 131, 59 | recalled [R] | §4.2 credit |
| Cembranos, "Dark matter from R^2 gravity" | arXiv:0809.1653 (PRL 102, 141301) | abstract | §4.2 related work |
| Laine, Meyer; Borsanyi et al. | arXiv:1503.04935; 1606.07494 | abstract (inputs of the SS tables) | §1.5 credit |
| Aubourg et al. 2015 (r_d fit); Hu-Sugiyama 1996 (z_*) | arXiv:1411.1074; ApJ 471, 542 | recalled formulas [R] | §5 |

Downloads (scratchpad only): `.../scratchpad/qtheory/scenario_dl/{ss,desi,pplus,abs,gp}`.
