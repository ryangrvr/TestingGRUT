# Module cosmo_observable: cosmology of a persistent Z2 wall network in unchanged Candidate 1

21 October 2026. Claude calculation module for the issue #13 wall-history packet (BRIEF.md). The action is
unchanged: no bias, no portal and no inflaton coupling. Whether a persistent, unbiased, percolating network forms
at all is decided by the history modules (thermal_potential, charge_phase, hidden_heating, inflation_spectator,
network_sim). Inflation and reheating are outside-action preparation data. **This module takes the existence of
such a network as its premise.** It then computes the background expansion that follows, with expansion and
backreaction included, and one observable whose magnitude is fixed before comparison.

Nothing here validates or excludes Candidate 1. The results are **conditional history constraints**. Every
comparison with data (H0, age, acceleration, Omega_c h^2, N_eff, CMB bounds) uses numbers that were already known.
These comparisons are **retrospective**: each is a consistency test of a stated history, not a new prediction.

Labels: **[D]** derived in closed form here; **[C]** computed by the code here; **[Q]** quoted from literature
or from the repository; **[E]** order-of-magnitude estimate. Every headline number is in
`results/headline.json`, and the full outputs are in the other `results/*.json` files.

## Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/cosmo_observable/code
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py      # about 27 s on one core
```

`run_all.py` runs these steps in order:
- `normalization.py` (part 3)
- `history_benchmark.py` (parts 1-2)
- `scale_family.py` (part 4)
- `bubbles_gw.py` (part 5)
- `checks.py` (35 checks)

Core physics is in `core.py` (inputs, units, g_*, background solver) and `observables.py` (H0, age, acceleration,
fractions, tuned vacuum). The code uses only numpy and scipy, with no random numbers. Two consecutive runs gave
identical JSON apart from the runtime fields.

## 0. Result in brief

1. **Domination [C].** Take the tree tension sigma = 9.428090416e7 GeV^3 (v = 1 TeV, lambda_phi = 0.01) and the
   HKS scaling law rho_w = A sigma/t with A = 0.8. Including backreaction, the walls overtake radiation at
   **T_dom = 6.367 keV** (z = 2.71e7, t = 2.71e4 s). The range over all variants is 5.45-7.34 keV:
   - A = 0.7 to 0.9;
   - backreaction on or off;
   - sigma/t versus the 2A sigma H convention.

   Walls dominate after BBN and long before recombination.
2. **Observable [C].** Define R_H = H0(predicted from T_CMB and this history) / H0(measured). In the history the walls
   dominate from T_dom until today. The H0 range [67, 74] km/s/Mpc is known retrospective data.
   - **Frozen network** (rho_w ∝ a^-1, w = -2/3): **R_H = 1.30e9** against 70 km/s/Mpc (1.23e9-1.36e9 across
     74-67). Varying A, the switch point and the scaling convention gives 1.14e9-2.06e9. The predicted age is
     21.4 yr. The universe accelerates for all z < 2.15e7.
   - **Bracket** rho_w ∝ a^-2 (w = -1/3): **R_H = 2.50e5** (range 2.32e5-2.76e5). The age is 5.6e4 yr and the
     universe never accelerates.

   Either way the history differs from the known H0 range by at least **5.4 orders of magnitude**. The conclusion
   does not depend on the post-domination law. Suppose the walls instead convert at T_dom into stable matter (n=3)
   or into radiation (n=4). Then **Omega h^2 ≈ 1.1e3** (n=3) or **Delta N_eff ≈ 7.4 at recombination** (n=4).
3. **Wall fractions [C].** BBN is only mildly affected:
   - at T = 1 MeV, Omega_w = 1.92e-5 (Delta N_eff-equivalent 1.18e-4);
   - at T = 0.1 MeV, Omega_w = 2.93e-3 (Delta N_eff-equivalent 0.0214).

   At recombination (T = 0.2562 eV) the walls dominate completely: 1 - Omega_w = 2.7e-13 (n=1) or 6.7e-9 (n=2).
   A subdominant-network CMB template therefore does not apply to this history, and none is used.
4. **Tuning [C].** The vacuum energy of the action is not supplied. Setting rho_Lambda = 0 gives R_H = 0.539 even
   without walls. A negative rho_Lambda can force H0 = 70, but only by cancelling rho_w0 to 1 part in 2.4e18 (n=1)
   or 8.8e10 (n=2). The cancellation leaves q0 = 8.5e17 or 6.3e10 and an age of 34 yr or 8.8e4 yr.
5. **Normalization [C, Q].**
   - Conversion: 1 GeV^3 = 4.114693e21 J/m^2 = 4.578214e4 kg/m^2 (mass per area).
   - Literal reading: the Lazanu-Martins-Shellard (LMS) printed bounds are 4.22e-9 kg/m^2 (R+M) and 3.85e-9 kg/m^2
     (R+M+Lambda). Read literally they give sigma^(1/3) = 45 keV, which is **inconsistent with their own
     0.96/0.93 MeV labels**.
   - Reconciled reading: if the unit is read as kg/cm^2, sigma = 9.22e-10 and 8.41e-10 GeV^3, with
     sigma^(1/3) = 0.973 and 0.944 MeV. This matches the labels within 1.5%. It also matches their stated agreement
     with deltaT/T ~ 1e-5, since G sigma t0 = 4e-6 rather than 4e-10. Their "energy scale" is then sigma^(1/3)
     in natural units, not a vev.
   - Comparison: Codex's dimensional value G sigma/H0 = 1e-5 is reproduced as 2.2257e-9 GeV^3, with the benchmark
     4.236e16 times above it. The benchmark is 1.12e17 times the reconciled LMS (R+M+Lambda) bound.
   - These are bounds for a subdominant network. They are used only for part 4.
6. **Scale family [C].** Hold lambda_phi = 0.01 and g_X = 0.1 fixed. A persistent network does not dominate before
   today only if **v < 47.3-68.1 MeV** (m_X < 4.7-6.8 MeV). The reconciled CMB bound for a subdominant network needs
   **v < 2.07 MeV** (m_X < 0.207 MeV). Read literally, the LMS bound gives v < 96 keV. At m_X = 0.207 MeV the
   product bound m_X|r| < 63.93 GeV only requires |r| < 3.1e5. No viable point is claimed, because the other
   consequences of an MeV-scale hidden sector were not computed.
7. **Collapsing minority bubbles [C, E].** In the unchanged action the products h, X and D are stable:
   - h -> X Xbar and h -> DD are kinematically closed;
   - nothing in the hidden sector is lighter than X;
   - there is no portal.

   The injected rest mass therefore remains as relic matter. Requiring Omega h^2 <= 0.12 means the collapse must end
   by roughly **T_end ≳ 60 MeV × (f c_b p_-)**. The scaling is linear within a factor of 2 down to T_eq. The values are:

   | f c_b p_- | T_end must be at least |
   |---|---|
   | 1 | 60.3 MeV |
   | 0.1 | 7.8 MeV |
   | 1e-2 | 0.82 MeV |
   | 1e-3 | 0.107 MeV |
   | 1e-6 | 126 eV |
   | 1e-10 | 0.029 eV |

   This is equivalent to rho/s <= 4.35e-10 GeV.
8. **Gravitational waves at domination [E].** The HKS amplitude gives Omega_gw(T_dom) = 7.3e-3 relative to the total
   density (range 2.4e-3 to 1.45e-2). The present-day rho_gw/rho_c,100 is 6.1e-7 (2.0e-7 to 1.2e-6), at a peak
   frequency of 8.0e-13 Hz. In the wall-dominated n=1 history this is only 7e-25 of the total density today.
   **This is an estimate only.**

## 1. Conventions and inputs

Natural units with hbar = c = k_B = 1. Energies are in GeV, and the metric signature is (-,+,+,+) as in BRIEF.
The reduced Planck mass is M_P = 2.435323e18 GeV and G = 1/(8 pi M_P^2). Today is defined by the photon temperature
T_CMB.

| Input | Value | Provenance |
|---|---|---|
| sigma | (2 sqrt(2 lambda_phi)/3) v^3 = 9.428090416e7 GeV^3 (sigma^(1/3) = 455.136 GeV) | ACTION / BRIEF; Codex DERIVATION.md at 135d535 [Q] |
| v, lambda_phi, g_X, m_h, m_X, m_D | 1000 GeV, 0.01, 0.1, 141.42, 100, 100.50 GeV | ACTION benchmark [Q] |
| A = rho_w t / sigma | 0.8 ± 0.1 | HKS arXiv:1309.5001 eqs. (4), (29); radiation-era Z2 lattice [Q]. LMS matter-era value 0.643 used as a variant |
| eps_tilde (GW) | 0.7 ± 0.4 | HKS eq. (30) [Q] |
| g_*rho(T), g_*s(T) | Saikawa-Shirai fitting functions, App. C (C.1)-(C.9); g_rho0 = 3.383, g_s0 = 3.931 | arXiv:1803.01038 [Q] |
| T_CMB | 2.7255 K | BRIEF; Fixsen arXiv:0911.1955 (2.72548 ± 0.00057 K) [Q], calibration |
| omega_b, omega_c | 0.02237, 0.1200 | Planck 2018 arXiv:1807.06209 Table 2 (TT,TE,EE+lowE+lensing) [Q], calibration |
| z_* | 1089.92 | same table; fixes T_rec = T0(1+z_*) = 0.2562 eV as a temperature [Q] |
| rho_Lambda | 0 (base); Planck-like 0.3113 rho_c100 (closes a no-wall universe at h = 0.6736); tuned negative | FREE (vacuum matching not supplied, BRIEF) |
| H0 data | [67, 74] km/s/Mpc: Planck 67.36 ± 0.54 (LCDM-inferred), Riess et al. 73.04 ± 1.04 | retrospective, known [Q] |
| Omega_DM h^2 limit | 0.12 | Planck 2018, retrospective [Q] |
| Neutrinos | massless, N_eff = 3.044 (implicit in the g_* fit) | simplification. Sum m_nu = 0.06 eV would add 6.4e-4 to omega_m and is irrelevant here |

**What A multiplies [Q, D].** HKS eq. (4) defines A ≡ rho_wall t / sigma_wall, where t is cosmic time and
sigma_wall = 2 sqrt(2 lambda) eta^3/3 is the same tree tension used here. The simulations are radiation-era, with
H = 1/(2t). So rho_w = A sigma/t = 2 A sigma H there. A convention of the form rho_w = A_H sigma H has
A_H = 2A = 1.6 in the radiation era and A_H = 1.5A in the matter era. LMS (arXiv:1505.03673 eq. 27) write
rho = (1 - lambda)(A eta/V) sigma/t for a ∝ t^lambda. Their matter-era value A eta/V = 1.93 gives
rho t/sigma = 0.643 in the HKS convention. The production runs use rho_w = A sigma/t, with t integrated
self-consistently. A variant uses rho_w = 2 A sigma H with H solved self-consistently. The two agree exactly in a
pure radiation era. With backreaction they differ by 9% in T_dom, and this spread is reported.

## 2. Equations and derivations

**Background [D].**

```
H^2 = (rho_r + rho_m + rho_w + rho_Lambda) / (3 M_P^2)
rho_r = (pi^2/30) g_*rho(T) T^4
s = (2 pi^2/45) g_*s T^3
a/a0 = (g_s0/g_*s(T))^(1/3) T0/T
rho_m = omega_m rho_c100 (a0/a)^3
dt/dlnT = -(1 + (1/3) dln g_*s/dlnT) / H
```

Here rho_c100 = 3 M_P^2 (100 km/s/Mpc)^2 = 8.095894e-47 GeV^4. The walls do not interact with the SM plasma
(there is no portal), so SM entropy is conserved and the photon temperature maps to a/a0 exactly. The code
integrates ln t in x = ln T with DOP853 at rtol 1e-10 (convergence 7e-10, see §8).

**Scaling phase with backreaction [D, C].** The wall density is rho_w = A sigma/t, with t(T) solved together with
H(T, rho_w). From d ln rho_w/dt = -1/t = -3H(1 + w), the scaling walls have 1 + w = 1/(3Ht). In a radiation era
this gives w = -1/3 (rho_w ∝ a^-2); in a matter era it gives w = -1/2 (a^-3/2).

The scaling phase starts at T_start = 100 GeV. Changing it to 1000 GeV changes T_dom by less than 1e-6. The
assumption is that friction is negligible: the SM does not couple to the walls, and hidden quanta of mass
≳ 100 GeV are Boltzmann suppressed at T ≲ 100 GeV.

**T_dom closed form (radiation era, no backreaction) [D].** Set rho_w = rho_r with t = 1/(2H). This gives

```
H_dom = 2 A sigma / (3 M_P^2)
T_dom = [40 A^2 sigma^2 / (pi^2 g_* M_P^2)]^(1/4)
t_dom = 3 M_P^2 / (4 A sigma)          (Saikawa 1703.02576 eq. 2.19)
```

The numbers are: closed form 5.8223 keV; numerical integration without backreaction 5.8221 keV; t_dom = 3.8823e4 s
numerically versus 3.8817e4 s from eq. (2.19) [C]. With backreaction (production) T_dom = 6.3669 keV and
t_dom = 2.7145e4 s, because the walls themselves speed up the expansion. With the 2 A sigma H convention,
T_dom = 6.924 keV, which is exactly 2^(1/4) times the no-backreaction value [D, C].

**Equation of state of a wall gas; post-domination laws [D].** Consider a planar wall moving at speed v normal to
itself. Its energy per area is gamma sigma, the normal stress is gamma sigma v^2, and the two tangential stresses
are -sigma/gamma. Averaging isotropically,

```
w = (v^2 - 2/gamma^2)/3 = v^2 - 2/3
rho_w ∝ a^-n,  n = 3(1 + w) = 1 + 3 v^2        (no energy loss to radiation)
```

- A frozen network that is stretched by the Hubble flow (v -> 0) has **n = 1, w = -2/3**. This is the standard
  isotropic non-relativistic wall gas: Saikawa 1703.02576 eq. (2.20) and its ref. [36] (Vilenkin-Shellard);
  Bucher-Spergel astro-ph/9812022; Battye-Bucher-Spergel astro-ph/9908047; Kawasaki-Takahashi hep-ph/0410158
  ("decreases as slow as the inverse of the scale factor").
- Once n = 1 walls dominate, the universe accelerates and the comoving Hubble radius shrinks. Walls are then
  stretched outside the horizon, and causal processes cannot remove them, so the frozen law is self-consistent.
- The bracket **n = 2** (v^2 = 1/3, w = -1/3) is the radiation-era scaling exponent continued past domination.
- Saikawa (1703.02576 sec. 2.2, footnote 4) states that dynamics in a wall-dominated universe are "uncertain".
  That is why n is treated as a bracket, not a result.
- A Z2 network has no junctions, so it is not a "frustrated" network in the Bucher-Spergel sense. The n = 1 law
  rests on super-horizon freezing during acceleration, not on frustration.
- n = 3 and n = 4 are **stress tests**, not networks. They represent all wall energy converted at T_dom into stable
  non-relativistic matter (for example horizon-mass black holes, §3c) or into radiation.

The scaling phase ends where rho_w = f_switch rho_r, with f_switch in {0.1, 0.3, 1, 3}; production uses 1.
Afterwards rho_w = rho_w,sw (a_sw/a)^n.

**H0 given T_CMB [D].** For T_dom ≪ m_e, g_*s(T_dom) = g_s0, and the wall density today is

```
rho_w0 = rho_r(T_dom) (T0/T_dom)^n
n = 1:  H0 = [pi^2 g_* T_dom^3 T0 / 90]^(1/2) / M_P
```

The ratio R_H = H0/H0_meas therefore scales as T_dom^(3/2) ∝ v^(9/4) for n = 1, and as T_dom ∝ v^(3/2) for n = 2.
The code checks rho_w0 = rho_w,dom (a_dom/a0) to 1e-9. The v^(9/4) and v^(3/2) slopes appear in the scan
(factor 11.8 and 5.2 per factor 3 in v).

**Acceleration, age, tuned vacuum [D].** The deceleration equation is

```
a''/a = -(rho + 3p)/(6 M_P^2),    rho + 3p = 2 rho_r + rho_m + (n - 2) rho_w - 2 rho_Lambda
```

For n = 1, acceleration starts where rho_w = 2 rho_r + rho_m, at T ≈ T_dom/2^(1/3) = 5.05 keV. The age today is
t0 = 2/H0 for a ∝ t^2 (n = 1) and t0 = 1/H0 for a ∝ t (n = 2).

If rho_Lambda < 0 is tuned so that H(T0) = H0_meas, then

```
q0 = (n rho_w0 + 4 rho_r0 + 3 rho_m0 - 2 rho_c) / (2 rho_c) ≈ n R_H^2 / 2
age ≈ pi/(R_H H0)        (n = 1)
age ≈ (pi/2)/(R_H H0)    (n = 2)
```

The ages follow from the integrals ∫dx/sqrt(x(1-x)) = pi and ∫dx/sqrt(1-x^2) = pi/2. The code evaluates H with
expm1 to avoid the cancellation at 1 part in R_H^2.

**Delta N_eff-equivalent [D].** Delta N_eff,eq(T) = rho_w(T)/rho_1nu(T), with rho_1nu = (7/8)(pi^2/15) T_nu^4.
T_nu/T comes from instantaneous decoupling: (T_nu/T)^3 = (2 + g_e(m_e/T))/(11/2), where g_e is the exact e±
entropy weight obtained by quadrature. T_nu/T = 0.99407 at 1 MeV and 0.76325 at 0.1 MeV.

## 3. Results for the benchmark (v = 1 TeV): parts 1-2

**(a) Domination [C].**

| Variant | T_dom (keV) | z_dom | t_dom (s) |
|---|---:|---:|---:|
| A = 0.8, sigma/t, backreaction (production) | 6.3669 | 2.711e7 | 2.7145e4 |
| A = 0.7 / 0.9, same | 5.9557 / 6.7530 | 2.536e7 / 2.875e7 | 3.10e4 / 2.41e4 |
| A = 0.8, no backreaction | 5.8221 | 2.479e7 | 3.8823e4 |
| A = 0.8, 2 A sigma H, backreaction | 6.9240 | 2.948e7 | 2.316e4 |
| All variants | 5.45 - 7.34 | | |

**(b) Full history to T_CMB, rho_Lambda = 0, A = 0.8, f_switch = 1 [C].**

| | n = 1 (frozen, w = -2/3) | n = 2 (w = -1/3) |
|---|---:|---:|
| rho_w0 (GeV^4) | 6.746457e-29 | 2.488682e-36 |
| H0 predicted (km/s/Mpc) | 9.1286e10 | 1.7533e7 |
| **R_H** vs 67 / 70 / 74 | 1.362e9 / **1.304e9** / 1.234e9 | 2.617e5 / **2.505e5** / 2.369e5 |
| R_H (vs 70) over A in [0.7, 0.9], f_switch in [0.1, 3], both conventions | 1.145e9 - 2.063e9 | 2.317e5 - 2.758e5 |
| Age | 21.4 yr (13.5-24.4) | 5.58e4 yr (5.07e4-6.03e4) |
| Acceleration | z in [2.15e7, 0] (T from 5.05 keV to T0) | none; q -> 0+ (q0 = 2.3e-12) |
| q0 | -0.5000 | +2.3e-12 |
| 1 - Omega_w at T_rec = 0.2562 eV | 2.68e-13 | 6.67e-9 |
| Tuned rho_Lambda < 0 forcing H0 = 70: tuning 1 - abs(rho_Lambda)/rho_w0 | 4.2e-19 | 1.1e-11 |
| Tuned: q0, age | 8.50e17, 33.6 yr | 6.27e10, 8.76e4 yr |

With the Planck-like rho_Lambda the n = 1 and n = 2 numbers are identical to the digits shown, because rho_w0
dwarfs rho_Lambda. The no-wall reference is H0 = 37.74 km/s/Mpc (R_H = 0.539) for rho_Lambda = 0, and
67.36 km/s/Mpc for the Planck-like value. The action alone therefore does not give R_H ≈ 1 even without walls,
because vacuum matching is not supplied.

The value q0 = -0.5 for n = 1 is close to the LCDM value of about -0.53. That closeness is coincidental and is not
a test. The acceleration begins at z ≈ 2e7, not z ≈ 0.6, and H0 is off by a factor 1e9.

**(c) Wall fractions [C].** These are set before T_dom, so they are the same for n = 1 and n = 2.

| T | Omega_w | rho_w/rho_r | Delta N_eff-equiv (A = 0.7-0.9) |
|---|---:|---:|---:|
| 1 MeV | 1.917e-5 | 1.917e-5 | 1.18e-4 (1.04e-4 - 1.33e-4) |
| 0.1 MeV | 2.931e-3 | 2.940e-3 | 0.0214 (0.0187 - 0.0240) |
| 0.2562 eV (recombination) | 1 - 2.7e-13 (n=1); 1 - 6.7e-9 (n=2) | 1.53e13 (n=1) | not meaningful; walls are not radiation |

The wall fraction grows as T^-2 in the radiation era. In this history BBN alone is not a strong test: the
equivalent Delta N_eff is at most 0.024 by T = 0.1 MeV. The walls dominate before z_* by a factor of about 1e13.

**Stress tests at T_dom [C].**
- n = 3, all wall energy into stable matter: Omega h^2-equivalent 1134 (332-2162 over A and f_switch), R_H = 48
  (26-66).
- n = 4, all wall energy into radiation: R_H stays at the no-wall value 0.539, but Delta N_eff at T_rec = 7.45.

Both conversions are rejected by known data (Omega_c h^2 = 0.12; N_eff ≈ 3), retrospectively.

**Gravitational scales [E].**
- The planar-wall repulsion time t_sigma = 1/(2 pi G sigma) is 1.66e5 s.
- A closed spherical wall reaches r_s(R) = R at R = M_P^2/sigma = 4.1e4 s × c.
- r_s exceeds the wall width once R > 4.2 mm.
- At t_dom the horizon mass is 4.7e9 M_sun, and a horizon-size closed wall has mass 4 pi sigma t^2 = 1.8e9 M_sun.

Horizon-size closed walls near t_dom are therefore self-gravitating. Black-hole formation would act as matter,
which is the n = 3 stress test. No black-hole fraction is computed.

## 4. The observable: R_H

**Definition.** R_H = H0(predicted | T_CMB, persistent-network history) / H0(measured).

**Inputs.** R_H is fixed by the following, none of which is fitted to H0 or to any expansion data:
- sigma (action);
- A (literature);
- the post-domination law n (bracketed);
- T_CMB (calibration);
- rho_Lambda (free; base 0);
- the SM g_*.

omega_b and omega_c do not enter at any digit shown.

**Data.** The comparison range H0 ∈ [67, 74] km/s/Mpc is known retrospective data. It was not held out, and R_H is
not a new prediction. A consistent history would need R_H within about [0.905, 1.105].

**Result.** The computed R_H lies in **2.3e5 - 2.1e9** for every persistent-network variant with n in [1, 2].

**Conclusion (conditional history constraint).** Suppose an unbiased persistent Z2 network with the tree tension at
v = 1 TeV, lambda_phi = 0.01 formed in our past and survived (no bias or portal, as the unchanged action requires).
Then the expansion history is inconsistent with the known value of H0 by 5.4-9.3 orders of magnitude. With
T_rec = 0.2562 eV, the wall fraction at recombination is 1 - O(1e-13 to 1e-8). Tuning rho_Lambda < 0 does not
avoid the inconsistency; it moves it to q0 = 6.3e10-8.5e17 and an age of 34 yr-8.8e4 yr.

This rejects a history, not the action. Whether such a network forms is the open question for the other modules.

## 5. Normalization reconciliation: part 3

**Conversion [D].** 1 GeV^3 = (1 GeV)/(hbar c/1 GeV)^2:
- energy per area: 1.602176634e-10 J / (1.973269804e-16 m)^2 = **4.114693e21 J/m^2**;
- mass per area (divide by c^2): **4.578214e4 kg/m^2**;
- inverse: 1 kg/m^2 = 2.184258e-5 GeV^3.

An independent MeV/fm^2 route agrees to 1e-12. The benchmark is sigma = 4.316e12 kg/m^2 = 3.879e29 J/m^2.

| Bound | Printed | sigma (GeV^3) | sigma^(1/3) | benchmark/bound | v_max at lambda_phi=0.01 | m_X,max (g_X=0.1) |
|---|---|---:|---:|---:|---:|---:|
| LMS R+M, literal | 4.22e-9 kg/m^2, "0.96 MeV" | 9.218e-14 | 45.2 keV | 1.02e21 | 99.3 keV | 9.9 keV |
| LMS R+M+Lambda, literal | 3.85e-9 kg/m^2, "0.93 MeV" | 8.409e-14 | 43.8 keV | 1.12e21 | 96.3 keV | 9.6 keV |
| LMS R+M, reconciled (kg/cm^2) | | 9.218e-10 | 0.9732 MeV | 1.02e17 | 2.138 MeV | 0.214 MeV |
| **LMS R+M+Lambda, reconciled** | | **8.409e-10** | **0.9439 MeV** | **1.12e17** | **2.074 MeV** | **0.207 MeV** |
| ZKO-type, G sigma t0 < 1e-5 (Saikawa eq. 2.21; t0 = 13.797 Gyr) | "sigma^(1/3) ≲ O(MeV)" | 2.253e-9 | 1.311 MeV | 4.18e16 | 2.881 MeV | 0.288 MeV |
| Codex, G sigma/H0 = 1e-5, H0 = 70 | 2.23e-9, ratio 4.24e16 | 2.2257e-9 | 1.306 MeV | 4.236e16 | 2.869 MeV | 0.287 MeV |

**LMS reconciliation [C].** The printed SI values convert to sigma^(1/3) = 45.2 and 43.8 keV. These are a factor
21.25 and 21.23 below the printed MeV labels, which means sigma is off by a factor 9.6e3 ≈ (100)^2.

Read with the unit kg/cm^2 (equivalently 4.22e-5 and 3.85e-5 kg/m^2), the values give 0.9732 and 0.9439 MeV. That
is 1.4% and 1.5% above the labels; the residual is a common conversion-constant offset.

The LMS amplitude coefficient is 6 pi G sigma t0 (eq. 28). The kg/cm^2 reading gives G sigma t0 = 4.1e-6 and a
matter-era rho_w/rho_bar = 0.643 × 6 pi G sigma t0 = 5.0e-5. This is consistent with their statement that the
result agrees with the deltaT/T ≲ 1e-5 estimate. The literal reading gives G sigma t0 = 4.1e-10, which would be
far below cosmic variance.

Conclusion: the "energy scale of formation" in LMS is sigma^(1/3) in natural units (ħ = c = 1), not a vev. Reading
it as a vev at lambda = 0.01 would give sigma = 8.3e-11 GeV^3, 11× smaller. This is our reconciliation, not an
erratum by the authors. The literal reading is carried as a stricter alternative.

**Scope.** LMS adopt a LambdaCDM background and treat walls as first-order perturbations (sec. 2: "the domain walls
always have a small contribution"). ZKO's argument also assumes subdominance today. Neither applies to the
wall-dominated benchmark history of §3, and neither is used there. They enter only part 4.

## 6. Scale dependence at fixed lambda_phi = 0.01, g_X = 0.1: part 4

sigma = 0.0942809 v^3, so v_max = (3 sigma_max/(2 sqrt(2 lambda_phi)))^(1/3) = 2.1972 sigma_max^(1/3) [D].

**No domination before today [C].** This is the root of rho_w(T0) = rho_r0 + rho_m0 + rho_Lambda, for a network
that scales from early times to today with backreaction.

| Variant | sigma_max (GeV^3) | v_max (MeV) | m_X,max (MeV) |
|---|---:|---:|---:|
| A = 0.8, Planck-like rho_Lambda | 2.393e-5 | 63.32 | 6.33 |
| A = 0.643 (LMS matter era) / 0.7 / 0.9, Planck-like | 2.976e-5 / 2.735e-5 / 2.127e-5 | 68.09 / 66.20 / 60.88 | 6.81 / 6.62 / 6.09 |
| A = 0.8, frozen (a^-1) after the onset of acceleration | 2.230e-5 | 61.85 | 6.18 |
| A = 0.8, rho_Lambda = 0 | 9.976e-6 | 47.30 | 4.73 |
| Analytic 3 M_P^2 H0^2 t0/A (no backreaction) | 3.04e-5 | | |

The range is **v_max = 47.3-68.1 MeV**. The rho_Lambda = 0 value is confirmed independently: the scan gives
T_dom = T0 at v = 47.298 MeV. Note that rho_Lambda = 0 together with v ≈ 59 MeV (n=1) or 74 MeV (n=2) gives
R_H = 1. A wall network could therefore act as "dark energy" with w = -2/3, but sigma would be 2.3e4-4.5e4 times
the reconciled CMB bound. This is noted, not pursued.

**CMB (subdominant) [C].** v < 2.074 MeV (m_X < 0.207 MeV) from the reconciled LMS R+M+Lambda bound, or v < 96 keV
read literally.

**Characteristic scales [C].** At fixed lambda_phi the persistent-network T_dom equals:

| T_dom | v |
|---|---:|
| 1 MeV (domination before BBN) | 35.2 TeV |
| 0.1 MeV | 6.58 TeV |
| T_rec | 1.571 GeV |
| T0 (rho_Lambda = 0) | 47.3 MeV |

From these points, T_dom ∝ v^(3/2) [D].

**Product bound [C].** m_X |r| < 63.93 GeV (primary; 74.64 GeV stress) gives:
- |r| < 3.08e5 (3.60e5) at m_X,max(CMB) = 0.207 MeV;
- |r| < 1.0e4 at m_X,max(dom) = 6.33 MeV.

At fixed couplings the wall-tolerability conditions and the product bound do not conflict on paper. **No viable
point is claimed.** At v ~ MeV, the relic abundances and N_eff of the MeV-mass hidden states h, X and D, the charge
floor and the original distance-envelope test have not been recomputed. Nor has any other consequence of rescaling
v. Lowering lambda_phi changes sigma only as sqrt(lambda_phi); the formula above gives v_max(lambda_phi).

## 7. Collapsing minority bubbles and gravitational waves: part 5

**Population model (stated order of magnitude) [E].** The minority-vacuum domains have volume fraction p_- and
form closed walls. Domains of size ~ t at time t have number density n_b ~ p_-/t^3, and each releases E ~ sigma t^2
on collapse within about one Hubble time:

```
d rho_inj / d ln t = c_b p_- sigma / t,       c_b in [1, 3]
```

c_b = 3 corresponds to spheres with R = t (area/volume = 3/R), and c_b = 1 to the literal n_b E. Here p_- means the
minority fraction in domains whose size is about the horizon at that time. A fraction f of the injected energy ends
as rest mass of stable quanta; f ∈ [0.01, 1] is a free input. f -> 1 if collapse forms black holes, which are
stable over the age at these masses.

Injection per Hubble time relative to rho_tot is c_b p_- sigma/(t rho_tot) [C]. It reaches 1 at:

| c_b p_- | T where injection per Hubble time = rho_tot |
|---|---:|
| 1 | 6.51 keV (≈ T_dom, as it must) |
| 0.1 | 2.06 keV |
| 0.01 | 0.65 keV |
| 1e-3 | 0.21 keV |

**Products in the unchanged action [D].**
- m_h = 141.42 GeV < 2 m_X = 200 GeV and < 2 m_D = 201.0 GeV, so h -> X Xbar and h -> DD are kinematically forbidden.
- h couples to the SM only through gravity, so it is stable except gravitationally.
- hh -> X Xbar and hh -> DD are open (282.8 GeV), as is DD -> X Xbar. X Xbar -> DD is closed at rest.
- X is the lightest hidden state and carries the U(1)_X charge, and D is Z2_D-odd; both are stable.
- No hidden state is lighter than X, and the action has no portal.

The injected rest mass therefore cannot annihilate into anything lighter, except through gravity. It survives as
relic matter. If every h converts by hh -> X Xbar, the retained rest-mass fraction is 2 m_X/(2 m_h) = 0.707.

**Relic [C].** The injection continues at constant p_- until T_end, with the energy kept as non-relativistic rest
mass (fraction f) from injection on:

```
rho_relic(T0) = f c_b p_- ∫ sigma/t^2 (a/a0)^3 dt
```

In a radiation era this equals 2 f c_b p_- sigma/(t_end s(T_end)) × s0, which agrees with the numerical integral
to 0.7%. The integral is dominated by the latest collapses. The requirement Omega h^2 <= 0.12 is rho/s <= 4.349e-10 GeV.

| f c_b p_- | Minimum T_end for Omega h^2 <= 0.12 |
|---:|---:|
| 1 | 60.3 MeV |
| 0.1 | 7.80 MeV |
| 1e-2 | 0.823 MeV |
| 1e-3 | 0.107 MeV |
| 1e-6 | 126 eV |
| 1e-10 | 0.029 eV (matter era) |

In a single-epoch variant every minority domain collapses at t_c, with rho_inj = 3 p_- sigma/t_c. Its maximum f p_-
is a factor 1.5-1.7 below the continuous maximum f c_b p_-, so the two variants agree in order of magnitude
(`relic_vs_T_end` in `results/bubbles_gw.json`). Below T_eq the bound tightens faster than linearly. This is a conditional constraint linking the minority-domain size spectrum to the
history. It needs p_-(t) from the inflation_spectator or network_sim modules, and is not a predicted abundance.

**Gravitational waves of a scaling network at domination [E].** HKS eqs. (32)-(33) and (38) give:

```
(d rho_gw/d ln k)_peak = eps_tilde G A^2 sigma^2
f_peak(t0) = (a_dom/a0) H_dom
```

- In the analytic radiation-era limit Omega_gw(t_dom) = 3 eps_tilde/(32 pi) = 0.0209 for eps_tilde = 0.7,
  independent of A and sigma [D].
- Computed with backreaction and normalized to the total density at T_dom: 7.30e-3 (range 2.40e-3 to 1.45e-2 over
  eps_tilde 0.3-1.1 and A 0.7-0.9).
- Redshifting as a^-4: rho_gw0/rho_c,100 = 6.11e-7 (2.0e-7 to 1.2e-6). This agrees with HKS eq. (35) to 1%.
- Peak frequency today: 8.0e-13 Hz. The f^-1 tail extends up to about 1e18 Hz (the wall-width scale; order of
  magnitude only).
- The peak is equivalent to about 0.11 in Delta N_eff per ln k at T_rec.
- In the n = 1 history it is 7e-25 of the wall-dominated total today.

HKS amplitudes are measured in radiation-era scaling, and T_dom is the edge of their validity. Emission after
domination is not modelled.

## 8. Checks, and what failed

`results/checks.json` holds 35 checks; **33 pass and 2 fail.** The passes include:
- g_*(T0) against Saikawa-Shirai eqs. (3.18)-(3.19);
- continuity at 120 MeV;
- Omega_gamma h^2;
- the BRIEF tension;
- two independent unit routes;
- reproduction of Codex's 2.23e-9 and 4.24e16;
- Saikawa eq. (2.19);
- the closed-form T_dom;
- insensitivity to T_start;
- production ODE convergence (7e-10);
- q0 = -1/2 and 0 for n = 1 and 2;
- age = 2/H0 and 1/H0;
- the tuned-vacuum ages pi/(R_H H0) and (pi/2)/(R_H H0);
- two code paths for v_max;
- the analytic relic;
- HKS eq. (35).

Failures and corrections, all retained:
1. **FAIL (expected, informative):** "LMS literal SI reading reproduces the 0.96 MeV label" fails by a factor 21 in
   scale. It documents the unit inconsistency in arXiv:1505.03673v1 (§5).
2. **FAIL:** "T_dom at ODE tolerance 1e-7 vs 1e-11 within 1e-6". The check as first written gave 1.39e-6. It is
   kept as a record. Production uses rtol 1e-10, and a production-vs-1e-12 check passes at 7.1e-10.
3. **Corrected modelling error.** The first scale-family scan switched from scaling to a^-n where rho_w = rho_r
   (switch_on = 'rad'). For v ≲ 1 GeV that switch fell in the matter era, when the walls did not yet dominate the
   total, and the network froze too early. The scan now switches where rho_w = rho_r + rho_m + rho_Lambda. The
   benchmark is unaffected to 1e-4, because rho_m/rho_r ≈ 1e-4 at 6 keV.
4. **Corrected range error.** The relic saturation search first stopped at T = 10 eV and returned "none" for
   f c_b p_- = 1e-10. It now extends to T0.
5. Husdal arXiv:1609.04979 was downloaded but not used; the Saikawa-Shirai fits were used instead. A PDF text
   extraction printed an xref warning for one file, but the text was recovered.

## 9. What was not done

- **No network dynamics into wall domination.** There was no lattice or VOS evolution through domination. The
  dilution law n is a bracket. The A values are radiation-era (HKS) or matter-era (LMS) measurements, extrapolated
  with backreaction.
- **No CMB anisotropy calculation** (C_l, likelihood) for any history. No wall-dominated CMB template is used.
- **Physics left out.** Not computed: PBH formation fraction or PBH constraints, friction from a hidden bath,
  X-condensate-decorated or charged walls, and thermal or quantum corrections to sigma. The planar tree wall with
  X = D = 0 is used throughout.
- **Background simplifications.**
  - Massless neutrinos.
  - The radiation term in rho + 3p uses w = 1/3.
  - Recombination physics in a wall-dominated universe is not computed: T_rec is the standard T0(1+z_*), used as a
    temperature.
  - Hidden relics from network formation and any Candidate 1 charged component are omitted from the background.
    Being positive densities, they would only raise H0.
- **Not in this module.** Whether the network forms and percolates (other modules); the A in a Lambda era; GW
  emission after domination; and the GW spectral shape.
- **No viability study at MeV-scale v** (§6).

## 10. Free inputs, and the conditional constraints they support

**Free inputs** (the result depends on them as stated):
1. existence and survival of an unbiased, percolating Z2 network (preparation, outside this module and outside the
   action);
2. v and lambda_phi (benchmark 1 TeV, 0.01);
3. A (0.7-0.9; 0.643 variant);
4. the post-domination law n (1-2; 3-4 as stress tests);
5. f_switch (0.1-3);
6. the scaling convention (sigma/t or 2 A sigma H);
7. rho_Lambda (0, Planck-like, tuned negative);
8. calibrations T_CMB, omega_b, omega_c (retrospective);
9. for bubbles: p_-(t), c_b, f, T_end;
10. for GWs: eps_tilde.

**Conditional history constraints:**
- **CHC-1.** A persistent unbiased network at the benchmark tension is incompatible with the known H0 by a factor of
  at least 2.3e5 (n = 2) and typically 1.3e9 (n = 1). This holds for every A, switch point and convention examined,
  and for rho_Lambda >= 0. A negative rho_Lambda moves the inconsistency to q0 and the age (§4). Conversion of the
  wall energy at domination into matter or radiation is also incompatible (Omega h^2 ~ 1e3; Delta N_eff ~ 7).
  Admissible histories at v = 1 TeV therefore have no such network: either it never forms or percolates, or it is
  removed before domination by something the unchanged action does not contain. That last option requires a
  preparation or extension statement outside this module.
- **CHC-2.** At fixed lambda_phi = 0.01, a persistent network is tolerable against domination only for
  v < 47-68 MeV, and against the reconciled CMB bound only for v < 2.07 MeV (96 keV for the literal LMS reading).
- **CHC-3.** A non-percolating minority population must finish collapsing by roughly T_end ≳ 60 MeV × f c_b p_-
  (within a factor of 2; exact values in the §7 table), or else stable hidden relics exceed Omega h^2 = 0.12.
- **CHC-4.** At v = 1 TeV the network does not affect BBN appreciably (Delta N_eff-equiv ≤ 0.024 at T ≥ 0.1 MeV).
  Domination before BBN would need v ≳ 35 TeV.

## 11. Literature (read level)

| Source | Identifier | Read level | Used for |
|---|---|---|---|
| Hiramatsu, Kawasaki, Saikawa 2013 | arXiv:1309.5001 | sections II, III.B, IV.A, V (full text of those) | A definition (eq. 4), A = 0.8 ± 0.1 (eq. 29), eps_tilde = 0.7 ± 0.4 (eq. 30), GW eqs. (32)-(39) |
| Saikawa 2017 (review) | arXiv:1703.02576 | sections 1, 2.1-2.3, footnote 4 | scaling (2.17)-(2.18), t_dom (2.19), w = -2/3 and R ∝ t^2 (2.20), ZKO estimate (2.21)-(2.22), "wall-dominated dynamics uncertain" |
| Lazanu, Martins, Shellard 2015 | arXiv:1505.03673v1 | sections 1, 2, 4 (eqs. 27-28), 6, 7 | CMB bound 4.22e-9 / 3.85e-9 kg/m^2, 0.96 / 0.93 MeV; matter-era A eta/V = 1.93; subdominance assumption |
| Saikawa, Shirai 2018 | arXiv:1803.01038 | abstract, eqs. (3.18)-(3.19), Appendix C | g_*rho, g_*s fitting functions |
| Husdal 2016 | arXiv:1609.04979 | skimmed (structure only) | not used |
| Planck 2018 VI | arXiv:1807.06209 | Tables 1-2 | omega_b, omega_c, z_*, H0, age, Omega_Lambda (calibration and retrospective comparison) |
| Riess et al. 2022 | arXiv:2112.04510 | abstract | H0 = 73.04 ± 1.04 (retrospective upper end) |
| Fixsen 2009 | arXiv:0911.1955 | abstract | T_CMB |
| Bucher, Spergel 1999 | astro-ph/9812022 | abstract | frozen wall network stretched by the Hubble flow, w = -2/3 |
| Battye, Bucher, Spergel 1999 | astro-ph/9908047 | abstract | wall-dominated cosmology, p = -(2/3) rho |
| Kawasaki, Takahashi 2005 | hep-ph/0410158 | abstract | frustrated-wall density ∝ 1/a |
| Hattori et al. 2015 | arXiv:1510.03595 | abstract | cited by Saikawa for wall domination; not used quantitatively |
| Zel'dovich, Kobzarev, Okun 1974 | Zh. Eksp. Teor. Fiz. 67, 3 [Sov. Phys. JETP 40, 1 (1975)] | **not read** (not on arXiv); secondary via 1703.02576 sec. 2.2 and 1505.03673 sec. 1 | MeV-scale bound, in the form G sigma t0 ≲ 1e-5 |
| Vilenkin 1981 | Phys. Rev. D 23, 852 | **not read**; bibliographic entry from 1703.02576 ref. [24] | wall gravity, context |
| Vilenkin 1983; Ipser, Sikivie 1984 | Phys. Lett. B 133, 177; Phys. Rev. D 30, 712 | **not read**; bibliographic data from memory, not verified | repulsive planar-wall gravity, t_sigma = 1/(2 pi G sigma) (context for §3 gravitational scales) |
| Sikivie 1982 | Phys. Rev. Lett. 48, 1156 | **not read**; cited in 1703.02576 ref. [81] | wall problem, context |
| Kibble 1976 | J. Phys. A 9, 1387 | **not read**; cited in 1703.02576 ref. [22] and 1505.03673 ref. [1] | defect formation, context |
| Vilenkin, Shellard (book) | CUP 1994 | **not read**; cited via 1703.02576 ref. [36] | w = -2/3 of a non-relativistic wall gas (also derived here, §2) |
| Codex packet | git show 135d535: README.md, DERIVATION.md, FEASIBILITY.md, SOURCES.md | sections | tension, the G sigma/H0 estimate, the scope of the LMS use |

No paper text or PDF is copied into the repository. The downloads are in the session scratchpad.

## 12. Files

- `code/core.py`: inputs and provenance, unit conversions, g_* fits, neutrino temperature, background solver.
- `code/observables.py`: H0, age, acceleration intervals, wall fractions, tuned negative vacuum.
- `code/history_benchmark.py` → `results/history_benchmark.json` (parts 1-2, all variants).
- `code/normalization.py` → `results/normalization.json` (part 3).
- `code/scale_family.py` → `results/scale_family.json` (part 4; v scan from 3 MeV to 30 TeV).
- `code/bubbles_gw.py` → `results/bubbles_gw.json` (part 5).
- `code/checks.py` → `results/checks.json`.
- `code/run_all.py` → `results/headline.json`.
