# Module minority_relic: collapsing minority domains from inflation regimes B and C

21 October 2026. Claude module for the issue #13 wall-history packet (BRIEF.md). It closes gap 5 of
`verification/v_critic/CRITIC.md` §4: the consequence of the closed opposite-vacuum (minority, phi = -v) domains that
an inflationary spectator leaves when H_inf lies between about 140 GeV and the percolation threshold (1.6-1.8 TeV).

The action is **unchanged**. **Inflation, H_inf, N_\*, the start condition of the observable patch, the reheating
history (instant reheating at T_RH(H_inf) or any non-instant variant), and the domain-size and debris conventions are
outside-action preparation data [P].** Nothing here adds a bias, portal or inflaton coupling. Nothing here says that
Candidate 1 is excluded or validated. All results are **conditional history constraints**. Every comparison with data
(Omega_DM h^2 = 0.12, recombination, BBN epoch) uses known numbers and is **retrospective**.

Labels: **[D]** derived here, **[C]** computed by the code here, **[Q]** quoted (sibling module, verifier, BRIEF),
**[E]** order-of-magnitude estimate, **[P]** outside-action preparation datum.

## Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/minority_relic/code
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 run_all.py   # 7.4 min wall, 2 workers
```

Order: `collapse.py` → `production.py` (2 single-threaded workers, 6.8 min) → `thresholds.py` → `headline.py` →
`checks_mr.py`. numpy/scipy only, no random numbers. The code imports, read-only, the Fokker-Planck chain of
`inflation_spectator/code` (`common.Spectator`, `fp1d.make_grid`, `fp1d.odd_chain`, `fp1d.lambda1_mfpt`) and the
Saikawa-Shirai g_* fits of `cosmo_observable/code/core.py`; it writes no bytecode into those directories.
This module's shared file is `mr_common.py` (not `common.py`, which would shadow the inflation_spectator module).

## 0. Result in brief

1. **Production [C].** One forward solve of the inflation_spectator chain from the barrier cell gives the exact
   per-e-fold minority production dp_-/dn for both start conditions (§1). It reproduces the verifier's corrected p_-
   to ≤ 0.12 % (250, 335, 1000, 2000 GeV, N_* = 50 and 60) and the corrected Lambda_1 eigenvalue column (6.316e-3 at
   1 TeV, 2.130e-3 at 501 GeV). The fresh start is placed exactly at y = v/H (verifier correction).
2. **e-fold → size → re-entry [D].** A domain born n e-folds after the start of the window, with radius c_R/H_inf,
   re-enters at T_re given by `a(T_re) H(T_re) = a0 H0 exp(n')`, `n' = n - DeltaN - ln c_R`, `DeltaN = N_* - N_0`,
   `N_0 = ln(a_end H_inf/(a0 H0))`. With instant reheating [P], N_0 = 48.05-49.41 for H_inf = 130-2000 GeV, and in
   the radiation era `T_re = T_RH e^{-m}/c_R x [g_rho^{1/2} g_s^{-1/3}](T_RH)/[same](T_re)` (m = N_* - n e-folds
   before the end; the g factor is 1.87 below m_e). N_* = 60 is not N_0 for any w ≤ 1 history ending above 4 MeV
   (N_0 ≤ 57.8-59.4), so "N_* = 60" can only mean a fresh start 11-12 e-folds before today's Hubble scale exits (§2).
3. **Collapse energy [D, C].** Thin-wall sphere frozen outside the horizon, released at re-entry: it collapses at
   tau_c = 1.8670 r0 (network_sim: 1.867) and releases **E_c = kappa_E sigma R_re^2, kappa_E = 7.787 = 0.620 x 4 pi**
   (radiation era), with R_re = 1/H_re. Matter era: kappa_E = 1.888 (§3).
4. **Relic [D].** `rho/s = f C_E sigma ∫ (dp/dn) H(T_re)/s(T_re) dn'`, with C_E = 3 kappa_E c_A/(4 pi) = 1.859 c_A
   (c_A ≥ 1 irregular-shape area factor). This is identical to cosmo_observable's CHC-3 with c_b = C_E/2. The bound is
   rho/s ≤ 4.373e-10 GeV (Omega h^2 ≤ 0.12; Y_N ≤ 4.373e-12 at 100 GeV, 3.092e-12 at 141.4 GeV). The relic is
   dominated by domains that re-enter just above T_dom; collapses before 1 MeV alone never exceed it below 2 TeV.
5. **Thresholds (instant reheating, 1D X = 0 baseline) [C].** H_inf must lie below:

   | criterion | fresh N_*=50 | fresh N_*=60 | relaxed N_*=50 | relaxed N_*=60 |
   |---|---:|---:|---:|---:|
   | relic Omega h^2 ≤ 0.12, f = 1 | 505 | 438 | 362 | 367 |
   | relic, f = 0.1 | 672 | 605 | 489 | 522 |
   | relic, f = 0.01 | none below 2 TeV | none | none | none |
   | no domain re-entering after T_dom (8.0 keV) in today's Hubble volume (N < 1) | 242 | 221 | 212 | 212 |
   | late-domain mass ≤ Omega_DM (upper estimate) | 484 | 409 | 337 | 340 |
   | late-domain mass ≤ 1e-6 Omega_DM | 296 | 262 | 241.5 | 241.8 |
   | no domain re-entering after recombination (N < 1) | 917 | 560 | 371 | 387 |
   | no domain re-entering after 1 MeV (N < 1, count only) | 202 | 195 | 192 | 192 |
   | regime A (no domain of any size, N < 1e-3) | 145.6 | 145.6 | 145.6 | 145.6 |

   Systematic ranges, the "exit" convention and a DeltaN scan are in §6. The X-fluctuation correction (verifier-
   calibrated) lowers the relic and mass thresholds by 5-12 % and the count thresholds by 1-3 %.
6. **Critic's estimate H_inf ≲ 240 GeV [C].** Its arithmetic is reproduced (Lambda_1/2 = 1e-11 at 237.4 GeV; a
   wall-energy integral continued to 0.03 eV gives 245-248 GeV for the relaxed start). The integral it rests on,
   however, is dominated by domains that re-enter after T_dom, where the walls dominate their own patch and the wall
   energy is not a hidden-quanta relic. Separated: the hidden-quanta relic gives only ≲ 362-367 GeV (relaxed, f = 1);
   "no locally wall-dominated domain in our Hubble volume" gives ≲ 212 GeV (relaxed) or 221-250 GeV (fresh); a
   late-domain mass bound at 1e-6 Omega_DM gives 241-242 GeV (relaxed). So 240 GeV is the right order, but it is a
   bound on locally wall-dominated (black-hole-candidate) domains, not on the relic; a fresh start loosens all of them.

## 1. Per-e-fold production dp_-/dn [D, C]

**Chain.** inflation_spectator's square-root-approximation chain on y = phi/H > 0 with couplings and xi run to mu = H
(X = D = 0 path), killing rate kappa_0 = 2D/h^2 in the cell next to y = 0 (the phi-odd sector). For a start in cell
m0 the minority fraction is p_- = F/2, with F the absorbed probability.

**Identity [D].** Reversibility, pi_i Q_ij = pi_j Q_ji, gives

```
dp_-/dn (fresh: delta at y = v/H)       = (kappa_0/2) (pi_0/pi_m0) P(X_n = m0 | X_0 = 0)
dp_-/dn (relaxed: pi restricted to y>0)  = (kappa_0/2) (pi_0/Z_+) sum_m P(X_n = m | X_0 = 0)
```

One forward solve from the barrier cell gives both. The factor e^{-B} sits in the closed-form Boltzmann ratio, and
P(X_n = m0 | X_0 = 0) is the downhill probability. That probability is O(1e-2) after the roll time and is
exponentially small (the fresh-start lag) before it. BDF with component-wise error control (rtol 1e-9, atol 1e-120)
resolves the lag tail relative to itself: atol 1e-80 and 1e-140 and dn = 0.025 agree to 4 decimals in log10 dp/dn at
200 and 280 GeV (e.g. log10 dp/dn = -32.1854 at n = 10, 200 GeV, in all three). atol ≤ 1e-180 overflows scipy's
error norm (failed, §11).

**Verifier corrections applied.** (i) The grid spacing is adjusted so that a cell centre sits exactly at y = v/H
(the module's offset start biased p_- up by ≤ 2 %). (ii) Lambda_1 is the odd eigenvalue where it is resolvable
(> 1e-12), else the MFPT value; they agree at large barrier. (iii) H_perc = 1792/1617 GeV (N_* = 50/60) and
1424/1348 GeV after trapping are quoted from v_inflation as the upper end of regime C.

**Validation [C]** (`results/checks.json`, 31/32 pass): fresh p_- = 2.4148e-9/5.2814e-9 (250 GeV, N50/N60; verifier
2.412e-9/5.278e-9), 1.8205e-4/3.4028e-4 (335; 1.820e-4/3.402e-4), 0.18380/0.20754 (1 TeV; 0.18381/0.20755),
0.33011/0.34625 (2 TeV; 0.33011/0.34626). The thresholds p_- = 1e-3 at 367.8/355.1 GeV and 1e-10 at 239.0/236.7 GeV
match the verifier (367.8/355.07) and the module (239.1/236.6).

**Shape [C]** (`results/production.json`, 57 H values, 130-2000 GeV, n = 0-62 in steps of 0.05):
- Fresh start: dp/dn is suppressed by the lag. At 241 GeV log10 dp/dn = -40.7, -23.6, -17.9, -15.2, -13.6, -12.0 at
  n = 4, 8, 12, 16, 20, 28, and -10.6 at n = 60. Lambda_1/2 = 3.4e-11 there.
- Relaxed start: flat at about Lambda_1/2 after a transient. At 241 GeV dp/dn = 10^-8.51 at n = 0, 10^-10.25 at n = 4,
  and ≈ 10^-10.45 from n ≈ 20 on. The transient comes from near-barrier weight in pi restricted to y > 0. A
  quasi-stationary start would have no transient; the closed forms of §6.4 are that variant.

**What "a domain" means here [D, E].** Each crossing makes the patch that crosses the top a minority seed of comoving
radius ~ 1/(a_n H_inf) (c_R = 1). Near the top the drift is weak: mu = |V''(0)|/(3H^2) = lambda v^2/(3H^2) =
0.05-0.08 at 200-250 GeV. The sign pattern inside a crossing patch is therefore set over ~1/(2 mu) ≈ 6-10 further
e-folds, and the true domain-size distribution is spread over that range. This module carries the spread as the
offset s_off = -ln c_R ∈ [-2, +3] and through the DeltaN scan (only DeltaN - s_off enters). A two-point
(correlated) stochastic calculation was not done (§11).

## 2. Creation e-fold, comoving size and re-entry temperature [D, P]

Let n count e-folds from the start of the last-N_* window (where the fresh start is imposed), and let m = N_* - n be
the e-folds before the end. A domain of physical radius c_R/H_inf at creation has comoving wavenumber
k = a_n H_inf/c_R. It stays frozen in comoving size while super-horizon (network_sim NG and field theory:
R(tau = R0/2)/R0 = 0.94-0.97). It re-enters when a R = 1/H:

```
a(T_re) H(T_re) = a_end H_inf e^{-m}/c_R = a0 H0 e^{n'},     n' = N_0 - m - ln c_R = n - DeltaN + s_off
N_0 = ln(a_end H_inf / a0 H0),   DeltaN = N_* - N_0,   s_off = -ln c_R
```

n' is the number of e-folds after the exit of today's Hubble scale (comoving radius today c_R e^{-n'} x 4451 Mpc). The
map T_re(n') uses only the standard post-reheating history (SM entropy conservation, Saikawa-Shirai g_*, matter and a
Planck-like Lambda as calibration). It is exact for T_re below the end of any non-standard era. Reference values [C]:
n'(T_rec = 0.256 eV) = 3.06, n'(T_eq = 0.80 eV) = 3.83, n'(T_dom,net = 6.37 keV) = 12.47,
n'(T_dom,iso = 7.97 keV) = 12.69, n'(70 keV) = 14.88, n'(1 MeV) = 17.77 (inflation_spectator: 17.766 and 12.692).

**Instant reheating [P] (primary).** 3 M_P^2 H_inf^2 = (pi^2/30) g_rho T_RH^4, a_end = a(T_RH). This gives
T_RH = 0.97e10, 1.19e10, 1.32e10, 1.47e10, 1.89e10 and 3.79e10 GeV at H_inf = 130, 198, 241, 302, 500 and 2000 GeV
(v_inflation: 1.46e10 at 300 GeV). It also gives N_0 = 48.05-49.41 (module 48.115, 48.462, 49.064 at 150, 300,
1000 GeV; mine 48.119, 48.466, 49.068). In the radiation era

```
T_re = T_RH e^{-m} c_R^{-1} [g_rho^{1/2} g_s^{-1/3}](T_RH) / [g_rho^{1/2} g_s^{-1/3}](T_re)     (factor 1.87 below m_e)
```

Example at 250 GeV: T_re = 6.09e5 GeV, 1.83 MeV and 106 eV for m = 10, 30 and 40. The closed form agrees with the
numerical map to 4e-3, and with inflation_spectator's table (6.1e5 GeV, 1.8 MeV, 106 eV).

With instant reheating, DeltaN = N_* - N_0 = 1.0-1.9 for N_* = 50 and 10.6-11.9 for N_* = 60. For N_* = 60 the fresh
start occurs ~11.6 e-folds before today's Hubble scale exits, so part of the lag is spent on scales that are still
super-horizon.

**Reach of N_0 [D, C].** aH ∝ rho^{(1+3w)/(6(1+w))}, so N_0 - N_0(instant) = (1/12) ln(rho_inf/rho_x) for a kination
era (w = 1) down to rho_x, and -(1/12) ln(rho_inf/rho_RH) for matter-dominated reheating to T_RH:

| H_inf (GeV) | N_0 instant | N_0 max (kination until 4 MeV) | N_0 min (MD reheating, T_RH = 10 MeV) |
|---:|---:|---:|---:|
| 150 | 48.12 | 57.84 | 38.71 |
| 250 | 48.37 | 58.18 | 38.88 |
| 1000 | 49.07 | 59.10 | 39.34 |
| 1600 | 49.30 | 59.42 | 39.50 |

So N_0 = 60 cannot be reached for H_inf ≤ 1.6 TeV with any w ≤ 1 era ending above 4 MeV.

**Alternative ("exit" convention) [P].** DeltaN = 0: the fresh start is imposed when today's Hubble scale exits, as in
inflation_spectator. It needs N_0 = N_*, i.e. a non-instant, stiffer-than-radiation era. That is possible for N_* = 50
(+1-2 e-folds of kination), but not for N_* = 60. For late re-entry the results do not depend on the label N_*.
General non-instant histories are covered by the DeltaN scan (-8 to +21) of §6.3.

**Consequence for inflation_spectator's regime A [C].** With instant reheating today's Hubble volume holds only N_0
e-folds of scales. The "no domain of any size" boundary (N < 1e-3) is therefore 145.6 GeV for both N_* = 50 and 60.
The module's 138.1 GeV for N_* = 60 requires N_0 = 60, which is unreachable (my exit-convention value is 137.8).

## 3. Energy released by a collapsing super-horizon domain [D, C]

**Equation [D].** Nambu-Goto sphere in a ∝ tau^alpha, comoving radius r(tau). The action is
S = -4 pi sigma ∫ dtau a^3 r^2 sqrt(1 - r'^2), which gives

```
r'' = -2 (1 - r'^2)/r - 3 (a'/a) r' (1 - r'^2)
```

With u = gamma r' this becomes du/dtau = -2 sqrt(1+u^2)/r - 3(a'/a) u. The conformal Hamiltonian is
H_c = 4 pi sigma a^3 r^2 gamma = a E_phys. Its explicit time dependence gives, for eps = r^2 gamma,

```
d eps/dtau = -3 (a'/a) eps r'^2,     r' = -sqrt(1 - r^4/eps^2)
```

This form is regular at r → 0 and is used for the final plunge. Start: at rest at tau_i = r0/x, with x = 1e2-1e4
(results independent of x to 2e-5). Re-entry is a r0 = 1/H, i.e. tau_re = alpha r0. Released energy:
E_c = 4 pi sigma a_c^2 eps_c ≡ kappa_E sigma R_re^2, with R_re = a(tau_re) r0 = 1/H_re.

**Results [C]** (`results/collapse.json`):

| | radiation (alpha = 1) | matter (alpha = 2) |
|---|---:|---:|
| tau_c / r0 | **1.8670** (network_sim 1.867; v_network 1.8670) | 2.2663 |
| t_c / t_re (cosmic time) | 3.486 | 1.455 |
| r(tau_re)/r0 | 0.734 (network_sim NG 0.734) | 0.262 |
| eps_c / r0^2 | 0.1778 | 0.0911 |
| **kappa_E = E_c/(sigma R_re^2)** | **7.787 (= 0.620 x 4 pi)** | 1.888 (= 0.150 x 4 pi) |
| relative to the actual radius at tau_re | 14.47 | 27.6 |
| maximum wall energy before the plunge | 0.774 x 4 pi sigma R_re^2 at tau = 1.34 r0 | 0.197 x 4 pi at 1.85 r0 |

**Flat-space checks:** tau_c = 1.3110287771 r0 (exact Gamma(1/4)^2/(4 sqrt(2 pi))), and eps is conserved to 3e-13.

**Reading.** The wall starts to fall at tau ~ r0/2. By re-entry it has shrunk to 0.73 r0, and Hubble friction on the
infall partly offsets the stretching. The released energy is therefore 0.62 of the static sphere energy at re-entry.
The ansatz "E ~ sigma R^2" holds with factor 7.79. This number is for thin-wall NG without gravity. It is invalid for
supercritical domains (sigma/(M_P^2 H_re) > 1, T_re < 7.97 keV), which are treated separately in §5.

## 4. Relic of stable hidden quanta [D, C]

**Stable products [Q].** h (141.42 GeV), X (100 GeV) and D (100.50 GeV) are stable apart from gravity in the unchanged
action: h → X Xbar and h → DD are closed, there are no portals, and the gravitational lifetime is ≥ 4e38 s
(v_relic §1, CONFIRMED). The collapse debris starts relativistic in a cold hidden sector. A fraction f of the released
energy ends as rest mass of stable quanta [P: f ∈ [0.01, 1] is a free input; hh → X Xbar alone keeps 0.707 of the
rest mass (cosmo_observable)]. Number is conserved after the debris dilutes, so

```
rho_relic/s = m Y_N = f E_c n_dom/s = f C_E sigma ∫ dn' (dp/dn) H(T_re)/s(T_re),     C_E = 3 kappa_E c_A/(4 pi)
```

This is independent of the quantum mass. The mass enters only through Y_N = (rho/s)/m. The bound Omega h^2 ≤ 0.12
(retrospective) is rho/s ≤ 4.373e-10 GeV, i.e. Y_N ≤ 4.373e-12 (m = 100 GeV, the strict v_relic value) or 3.092e-12
(m = 141.42 GeV). C_E(sphere) = 1.859; irregular domains have c_A ≥ 1 (variant c_A = 3).

**Equivalence with CHC-3 [D, C].** cosmo_observable writes d rho_inj/d ln t = c_b p_- sigma/t, with p_- the minority
fraction per ln t. In a radiation era d ln t = 2 dn' and H = 1/(2t), so the two coincide for
dp/dn = 2 p_-(per ln t) and C_E = 2 c_b. Their c_b ∈ [1, 3] is c_A = 1.08-3.2 here. Check [C]: f c_b p_- = 1e-6 gives
T_end = 126.17 eV (theirs 126 eV). At 0.1 MeV and 60 MeV the offsets are 5.5 % and 14 %, because their t(T) carries
the g_* history and t ≠ 1/(2H) at those steps (one check FAILS at a 3 % tolerance; §11).

**Where the relic comes from [C].** For a flat dp/dn the integrand goes as 1/T_re, so the relic is dominated by the
latest collapses still counted, those re-entering just above T_cut. The cut is T_cut = T_dom,iso = 7.97 keV; later
domains are supercritical (§5). Collapses before 1 MeV alone never reach Omega h^2 = 0.12 below 2 TeV for f = 1 (at
1 MeV sigma/(M_P^2 H) = 1.1e-5). **The hidden-relic constraint is therefore set by collapses between BBN and T_dom.**
Closed form for the relaxed start [D]:

```
rho/s ≈ f C_E (Lambda_1/2) sigma H/s |_{T_dom,iso} = f C_E (Lambda_1/2) x 1.716e-6 GeV
=> Lambda_1 ≤ 2.74e-4/f   (C_E sphere)
```

At T_dom,iso, H = sigma/M_P^2, so sigma H/s = rho_rad/(3s) ≈ 0.215 T_dom. The closed form gives 372 GeV against the
full 362-367 GeV; the difference is the relaxed transient. It is checked against the integral at 198 GeV (5 %).

## 5. Late domains: self-domination, black-hole candidates, BBN [D, C, E]

- **Self-domination.** At re-entry the wall energy of a sphere of radius 1/H over the radiation energy inside it is
  sigma/(M_P^2 H_re). This equals 1 at **T_dom,iso = 7.972 keV** (critic I2: 7.97). The network value 6.37 keV
  (cosmo_observable) is the variant; it moves the count thresholds by ≤ 1.3 %.
- **Domains re-entering after T_dom** dominate their own patch. They are black-hole or baby-universe candidates
  (GVZ 1512.01819, quoted via inflation_spectator at abstract level). Their mass is ~ the horizon mass,
  4 pi M_P^4/sigma = 4.2e9 M_sun at T_dom,iso and more for later re-entry. Two measures:
  - **counts** N = ∫_0^{n'_X} (dp/dn) e^{3n'} dn' in today's Hubble volume (x 32.3 for the observable universe,
    chi_p a0 H0 = 3.185);
  - an **upper estimate of their mass density**. Each such domain is taken to collapse with all its enclosed energy
    plus its wall energy, the latter capped at the enclosed energy:
    rho_late/s ≤ ∫ (dp/dn) (rho_tot/s)(T_re) [1 + min(kappa_E sigma/(4 pi M_P^2 H_re), 1)] dn'.
  
  N ~ 1 such object of ≥ 4e9 M_sun is not by itself in conflict with data. Literature PBH-abundance limits were not
  imported. Because Omega_late rises by ~1 decade per 1-3 % in H, the thresholds for Omega_late ≤ (1, 1e-3, 1e-6)
  Omega_DM bracket any such limit.
- **Re-entry after recombination (or T_eq).** A sigma = 9.4e7 GeV^3 wall of horizon size inside the last-scattering
  volume would be grossly visible. N < 1 is the natural criterion.
- **Earlier re-entry.** For a sphere the Schwarzschild radius of the debris exceeds the wall width once
  T_re ≲ 180-200 GeV (`spectra.json`, r_s/delta). The relative precision of sphericity needed for a black hole is
  sigma/(M_P^2 H_re), about 1e-5 at 1 MeV. Irregular stochastic domains are not expected to meet it except near
  T_dom [E]. If some do, the energy goes into small black holes. These evaporate mostly into SM states (lowering f)
  and carry evaporation constraints that were not computed (§11).
- **After BBN (1 MeV; 70 keV variant).** Counts only. The debris is hidden and has no portal, so no photodissociation
  or n/p effect arises at this order. The local energy fraction at re-entry is sigma/(M_P^2 H_re) = 1.1e-5 at 1 MeV,
  0.036 at 42 keV and 1 at 7.97 keV.

## 6. Thresholds [C]

`results/thresholds.json`, `results/headline.json`. Thresholds are the lowest H at which a criterion is crossed. They
come from log-log interpolation on the 57-point grid (2.8 % spacing below 500 GeV). Every quantity is monotone below
the crossing, with one exception. The R_X-bound variant of the fresh start switches from R_X to the 2D ratio at
B = 3 (H ≈ 390 GeV), so its quantities step down there and are not monotone. Its thresholds (first crossing) are kept
only as a bound; for a fresh start, where X starts at zero, the factor-2 calibration is the relevant one.

### 6.1 Primary (baseline: 1D X = 0, c_R = 1, C_E sphere, T_cut = 7.97 keV, today's Hubble volume)

H_inf thresholds in GeV. i = instant reheating; x = exit convention (DeltaN = 0).

| criterion | F50 i | F50 x | F60 i | F60 x | R50 i | R50 x | R60 i | R60 x |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| relic f = 1 | 505.4 | 520.9 | 437.8 | 520.9 | 362.5 | 361.6 | 366.6 | 361.6 |
| relic f = 0.1 | 672.1 | 689.3 | 605.0 | 689.3 | 489.2 | 484.4 | 522.2 | 484.4 |
| relic f = 0.01 | n.r. | n.r. | n.r. | n.r. | n.r. | n.r. | n.r. | n.r. |
| relic f = 1, c_A = 3 | 462.9 | 476.8 | 401.6 | 476.8 | 338.6 | 337.9 | 341.3 | 337.9 |
| relic f = 1, T_cut = 6.37 keV | 497.5 | 513.5 | 430.2 | 513.5 | 356.6 | 355.8 | 360.5 | 355.8 |
| relic f = 1, only collapses before 1 MeV | n.r. | n.r. | n.r. | n.r. | n.r. | n.r. | n.r. | n.r. |
| N(T_re < 7.97 keV) < 1 | 242.5 | 249.9 | 221.3 | 249.9 | 211.9 | 211.9 | 212.0 | 211.9 |
| N(T_re < 7.97 keV) < 1e-3 | 226.5 | 233.0 | 208.9 | 233.0 | 201.9 | 201.8 | 201.9 | 201.8 |
| N(T_re < 6.37 keV) < 1 | 245.2 | 253.0 | 223.0 | 253.0 | 213.0 | 213.0 | 213.1 | 213.0 |
| same, observable universe (x 32.3) | 233.8 | 240.7 | 214.6 | 240.7 | 206.5 | 206.5 | 206.6 | 206.5 |
| Omega_late,max ≤ Omega_DM | 484.4 | 502.8 | 408.7 | 502.8 | 337.0 | 336.6 | 340.0 | 336.6 |
| Omega_late,max ≤ 1e-3 Omega_DM | 351.5 | 365.2 | 303.8 | 365.2 | 270.5 | 270.4 | 271.2 | 270.4 |
| Omega_late,max ≤ 1e-6 Omega_DM | 296.3 | 307.0 | 261.7 | 307.0 | 241.5 | 241.4 | 241.8 | 241.4 |
| N(T_re < T_eq) < 1 | 710.7 | 790.6 | 452.2 | 790.6 | 325.1 | 322.9 | 331.4 | 322.9 |
| N(T_re < T_rec) < 1 | 917.1 | 1038 | 559.9 | 1038 | 370.5 | 364.9 | 387.1 | 364.9 |
| same, observable universe | 719.4 | 826.0 | 429.6 | 826.0 | 309.6 | 307.4 | 315.0 | 307.4 |
| N(T_re < 1 MeV) < 1 | 202.2 | 204.7 | 195.1 | 204.7 | 192.4 | 192.4 | 192.4 | 192.4 |
| N(T_re < 1 MeV) < 1e-3 | 194.0 | 196.3 | 188.0 | 196.3 | 186.0 | 186.0 | 186.0 | 186.0 |
| N(T_re < 70 keV) < 1 | 220.9 | 225.6 | 207.9 | 225.6 | 202.3 | 202.3 | 202.4 | 202.3 |
| super-horizon minority wall crossing today's Hubble volume, N < 1 | n.r. | – | n.r. | – | n.r. | – | n.r. | – |
| N(all sizes) < 1e-3 (regime A) | 145.6 | 144.2 | 145.6 | 137.8 | 145.6 | 144.2 | 145.6 | 137.8 |

n.r. means the criterion is not crossed below 2 TeV. It is then not limiting below percolation (1.6-1.8 TeV), where a
network forms and the other modules take over. For comparison with inflation_spectator (exit convention): after
1 MeV < 1e-3 gives 196.3/186.0 (module 196.5/186.2); after 1 MeV < 1 gives 204.7/192.4 (module 205.1/192.7); after
T_dom < 1e-3 gives 233.0/201.8 (module 235.3/202.1). The module's 235.3 is ~1 % high because its count integral stops
at the grid point n = 12.5 instead of n'_X = 12.69 (a factor e^{0.57} in N). That is a small error in the module,
not in this one.

### 6.2 Systematics (instant reheating)

Ranges cover: the exit convention; the X factor 2 (fresh, calibrated); the 2D/1D rate ratio vs B (relaxed,
calibrated, v_inflation) and the Langer R_X bound; and s_off = +3, -2. The X-calibrated values use F_X = 2 (fresh) and
the 2D rate ratio (relaxed).

| criterion | fresh N50 | fresh N60 | relaxed N50 | relaxed N60 |
|---|---|---|---|---|
| relic f = 1 | 505 (X-cal 477; 465-545) | 438 (414; 396-521) | 362 (325; 314-364) | 367 (328; 316-367) |
| relic f = 0.1 | 672 (608; 597-722) | 605 (529; 513-689) | 489 (395; 358-496) | 522 (406; 362-528) |
| relic f = 0.01 | n.r. (X-cal 1163) | n.r. | n.r. | n.r. |
| N after T_dom < 1 | 242 (241; 234-258) | 221 (220; 215-250) | 212 (207; 206-212) | 212 (207; 206-212) |
| Omega_late ≤ Omega_DM | 484 (462; 452-533) | 409 (391; 356-503) | 337 (308; 301-338) | 340 (310; 302-340) |
| Omega_late ≤ 1e-6 Omega_DM | 296 (292; 280-320) | 262 (259; 249-307) | 241.5 (233; 232-242) | 241.8 (233; 232-242) |
| N after T_rec < 1 | 917 (866; 769-1038, one n.r.) | 560 (523; 513-1038) | 371 (331; 319-376) | 387 (341; 326-389) |
| N after 1 MeV < 1 | 202 (201; 198-207) | 195 (194; 191-205) | 192 (189; 189-192) | 192 (189; 189-192) |

### 6.3 DeltaN scan (fresh start; DeltaN - s_off is the only combination that enters)

DeltaN = -8 corresponds to kination; 1.0-1.9 to N_* = 50 instant; 10.6-11.9 to N_* = 60 instant; up to 21 to N_* = 60
with matter-dominated reheating to 10 MeV.

| DeltaN | -8 | -4 | 0 | 2 | 4 | 8 | 12 | 16 | 21 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| relic f = 1 | 716 | 590 | 521 | 498 | 480 | 454 | 435 | 422 | 410 |
| N after T_dom < 1 | 357 | 280 | 250 | 241 | 235 | 226 | 221 | 218 | 215 |
| N after 1 MeV < 1 | 233 | 214 | 205 | 202 | 200 | 197 | 195 | 194 | 193 |
| N after T_rec < 1 | n.r. | n.r. | 1038 | 833 | 729 | 615 | 549 | 512 | 481 |
| Omega_late ≤ Omega_DM | 774 | 590 | 503 | 476 | 456 | 427 | 407 | 393 | 380 |

As DeltaN → ∞ the fresh thresholds approach the relaxed ones (212, 192, 340 GeV).

### 6.4 Closed forms (quasi-stationary relaxed start, dp/dn = Lambda_1/2) [D, C]

- N(T_re < T_X) ≈ (Lambda_1/2) e^{3 n'_X}/3 < 1 requires Lambda_1 < 6 e^{-3 n'_X}:
  - 4.3e-23 for 1 MeV, giving 192.4 GeV;
  - 1.7e-16 for T_dom, giving 212.0 GeV;
  - 6.2e-4 for T_rec, giving 403.9 GeV.
- Relic: Lambda_1 ≤ 2.74e-4/f, giving 372.1 GeV at f = 1.

The full results (192.4, 211.9, 365-387, 362-367) differ only where the restricted-equilibrium transient matters (T_rec
count, relic).

## 7. Spectrum example: H_inf = 241 GeV, instant reheating [C] (`results/spectra.json`)

Lambda_1 = 6.84e-11 and B = 19.7 here. The table is the relaxed start, N_* = 50 (DeltaN = 1.64). Per e-fold of creation:

| n (from start) | m (before end) | n' | radius today | T_re | log10 dp/dn | domains per Hubble volume per e-fold | E_c per domain | sigma/(M_P^2 H_re) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 46 | 2.36 | 420 Mpc | 0.076 eV | -10.25 | 10^-7.2 | supercritical | 3e9 |
| 12 | 38 | 10.36 | 0.14 Mpc | 0.77 keV | -10.40 | 10^3.1 | supercritical | 110 |
| 16 | 34 | 14.36 | 2.6 kpc | 42 keV | -10.43 | 10^8.3 | 3.4e6 M_sun | 0.036 |
| 20 | 30 | 18.36 | 47 pc | 1.8 MeV | -10.44 | 10^13.5 | 0.32 M_sun | 1.1e-5 |
| 28 | 22 | 26.36 | 3.3e3 AU | 3.8 GeV | -10.46 | 10^23.9 | 2e-15 M_sun | 9e-13 |
| 48 | 2 | 46.36 | 1e3 km | 1.8e9 GeV | -10.46 | 10^49.9 | 3e-50 M_sun | 4e-30 |

The relic at this H is Omega h^2 = 10^-7.6 (f = 1). N(after T_dom) = 10^5.5, N(after 1 MeV) = 10^12.1,
N(after T_rec) = 10^-6.9 and Omega_late,max/Omega_DM = 10^-6.2. For the fresh start the corresponding values are
10^-12.3/10^-9.6 (relic), 10^-0.4/10^3.3 (after T_dom), 10^8.6/10^10.6 (after 1 MeV) and 10^-33/10^-12.7 (after
T_rec), each pair for N_* = 50/60 (`headline.json` → quantities_at_fixed_H).

## 8. Comparison with the critic's estimate (H_inf ≲ 240 GeV) [C]

- The critic's reading "relaxed start, f c_b ~ 1, injection continuing to late times needs Lambda_1/2 ≲ 1e-11" is
  reproduced: Lambda_1/2 = 1e-11 at **237.4 GeV**. A wall-energy relic integral continued to 0.03 eV (f = 1, C_E
  sphere) gives 245-248 GeV (relaxed), 304 GeV (fresh N_* = 60) and 420-444 GeV (fresh N_* = 50, exit).
- That integral is dominated by domains re-entering after T_dom. These are supercritical: the NG energy does not
  apply, the energy is not hidden-quanta rest mass, and they dominate their patch.
- The **hidden-quanta relic** alone (collapses before 7.97 keV) is weaker: ≲ 362-367 GeV (relaxed, f = 1), 325-328
  with the X correction, and no constraint for f ≤ 0.01.
- The **late population** gives:
  - "none in our Hubble volume": ≲ 212 GeV (relaxed), 221-250 GeV (fresh);
  - mass ≤ 1e-6 Omega_DM: 241-242 GeV (relaxed);
  - mass ≤ Omega_DM: 337-340 GeV (relaxed);
  - no horizon-scale wall after recombination: 365-387 GeV (relaxed).
- **Verdict.** The critic's 240 GeV is the right scale for a conservative "no locally wall-dominated domain / tiny
  PBH-like mass fraction" condition with a relaxed start. It is not a relic bound. A fresh start (inflation_spectator's
  default) raises every threshold, to 221-250 GeV for the count and 410-520 GeV for relic and mass. The weakest
  physically robust conditions (relic, Omega_late ≤ Omega_DM, no wall after recombination) sit at 337-387 GeV
  (relaxed). That is the regime B/C boundary (p_- = 1e-3 at 355-368 GeV), not 240 GeV.

## 9. Free inputs (all outside the action unless marked)

- H_inf [P]; N_* = 50/60 [P]; start condition (fresh at +v, or + restricted equilibrium) [P]
- reheating: instant at T_RH(H_inf), or DeltaN = N_* - N_0 for any other history [P]
- domain size at creation c_R (s_off) [P/E]; shape area factor c_A [E]
- debris rest-mass fraction f ∈ [0.01, 1] [P]
- criteria conventions: T_dom definition (7.97 vs 6.37 keV); count target (1 vs 1e-3; Hubble vs observable volume);
  T_BBN (1 MeV vs 70 keV)
- X-fluctuation correction (1D baseline vs verifier-calibrated 2D factors)
- Calibration and retrospective data: T_CMB, omega_m, Planck-like Lambda (re-entry bookkeeping), Omega_DM h^2 = 0.12,
  z_*
- Action/BRIEF: v = 1 TeV, lambda_phi = 0.01, sigma = 9.428090416e7 GeV^3, M_P

## 10. Conditional history constraints (not predictions)

Within the class "inflation precedes the hidden history, the hidden sector is never heated above T_c, no
non-thermal restoration, our patch at phi = +v at the start of the last N_* e-folds" (inflation_spectator regimes
B/C), with the tree tension at the benchmark:

1. **CHC-MR1 (relic, retrospective).** Stable hidden debris from collapses before T_dom stays below Omega h^2 = 0.12
   only if H_inf ≲ 362-367 GeV (relaxed) or 438-521 GeV (fresh) for f = 1. With f = 0.1 the bounds are 484-522 and
   605-689 GeV. For f ≤ 0.01 there is no constraint below the percolation threshold. Including the X correction and
   c_A ≤ 3 the f = 1 range is 314-545 GeV.
2. **CHC-MR2 (late domains).** Domains re-entering after T_dom = 7.97 keV dominate their patch (black-hole candidates
   of ≥ 4e9 M_sun). None is expected in today's Hubble volume only if H_inf ≲ 206-212 GeV (relaxed) or 215-258 GeV
   (fresh). Their mass density stays ≤ Omega_DM only if ≲ 301-340 (relaxed) or 356-533 GeV (fresh). Specific PBH
   limits were not imported.
3. **CHC-MR3 (recombination).** No horizon-scale minority wall re-enters after recombination in today's Hubble volume
   only if H_inf ≲ 319-389 GeV (relaxed) or ≲ 513-1038 GeV (fresh; one variant is not crossed below 2 TeV).
4. **CHC-MR4 (BBN, count only).** No minority domain re-enters after 1 MeV in today's Hubble volume only if
   H_inf ≲ 189-192 GeV (relaxed) or 191-207 GeV (fresh). No BBN effect was computed: the debris is hidden.
5. All thresholds lie inside regimes B/C, between 145.6 GeV (regime A, instant reheating) and 1.35-1.79 TeV
   (percolation). The numbers constrain a preparation history. They do not predict H_inf. They do not exclude or
   validate Candidate 1. GRUT remains empirically unvalidated.

## 11. Failures, corrections, limits, not done

**Failed or corrected during the work.**
- BDF with atol ≤ 1e-180 overflowed scipy's RMS error norm ("step size less than spacing"). Production uses
  atol = 1e-120; 1e-80 and 1e-140 agree to 4 decimals.
- The first shared file was named `common.py`. It would have shadowed inflation_spectator's `common.py` on the import
  path, so it was renamed `mr_common.py` before any result was produced.
- The first edit of the late-domain estimate had the wall/enclosed energy ratio low by a factor 3 (it used
  1/(3 M_P^2 H) instead of 1/(M_P^2 H)). This was caught on review and fixed before the reported numbers.
- Check FAILS at 3 % tolerance: CHC-3 T_end for f c_b p_- = 1e-3 (0.1128 MeV vs 0.107 MeV) and, with a 20 %
  tolerance, passes at 60 MeV (68.96 vs 60.3 MeV). Both come from the e±/QCD g_* steps, where cosmo_observable's
  integrated t(T) differs from 1/(2H); at 126 eV the agreement is 0.14 %. Not an error in either module, but recorded.
- Found in a sibling module: inflation_spectator's count integral stops at the grid point below n'_X (12.5 vs 12.69),
  which raises its "after T_dom < 1e-3" fresh threshold by ~1 % (235.3 vs 233.0). Its N_* = 60 regime-A boundary
  (138 GeV) needs an unreachable N_0 = 60; with instant reheating it is 145.6 GeV (§2).

**Limits (not done).**
- **1D X = 0 Fokker-Planck.** X enters only through the verifier-calibrated factors: fresh x 2 (measured at 330 GeV);
  relaxed 2D/1D rate ratio vs B; the Langer R_X as an upper bound. No 2D production curve dp/dn(n) was computed.
- **Overdamped stochastic picture.** It is marginal at m_phi ~ H (m_phi/H = 0.71 at 200 GeV, 0.59 at 241 GeV).
  Exponents are Hawking-Moss, prefactor errors move thresholds by < 1 %.
- **Domain geometry.** Domain size and shape from a 1-point FP: c_R and c_A are carried as conventions. The
  correlated sign field (two-point stochastic statistics, percolation of the correlated field) was not computed.
  Within regime C (p_- up to 0.3) the isolated-domain counting is only order of magnitude.
- **Collapse.** Thin-wall NG sphere, no gravity. Supercritical domains (T_re < 7.97 keV) are only bounded (Omega_late)
  and counted. There is no field-theory collapse of irregular super-horizon domains in an expanding background. Black
  holes from near-spherical collapse (r_s > wall width for T_re ≲ 200 GeV), their evaporation into SM states and the
  evaporating-PBH constraints were not computed. f, the debris spectrum and number-changing in dense debris were not
  computed (f is a free input).
- **Non-instant reheating** is covered only through DeltaN. Scales re-entering during a non-standard era are mapped
  with the radiation-era table; they carry negligible relic weight and do not enter the counts.
- **Matter era.** kappa_E(MD) was computed but is used only for counts. Re-entry after T_eq uses the full LambdaCDM
  H(T).
- **Not used:** charged walls, X/D condensates on walls (charge_phase/charged_exit_sim), and friction from a hidden
  bath (the hidden sector is cold in this class).

## 12. Literature

No new paper was read for this module. Physics taken from the literature is quoted through sibling modules and
verifiers with their stated read levels:
- Starobinsky-Yokoyama astro-ph/9407016 (FP formalism; abstract, via inflation_spectator)
- Garriga-Vilenkin-Zhang arXiv:1512.01819 (inflationary walls → black holes; abstract, via inflation_spectator)
- Wang et al. arXiv:1302.0421 (p_c; abstract, via inflation_spectator)
- Saikawa-Shirai arXiv:1803.01038 (g_* fits; sections, via cosmo_observable)
- Planck 2018 arXiv:1807.06209 (calibration and Omega_DM h^2; tables, via cosmo_observable)
- HKS arXiv:1309.5001 (T_dom,net through A; via cosmo_observable)

## 13. Files

- `code/mr_common.py`: inputs, sibling imports, background, re-entry map n' ↔ T, T_RH, N_0, T_dom,iso.
- `code/collapse.py` → `results/collapse.json`: NG collapse, kappa_E (RD, MD), flat-space check.
- `code/production.py` → `results/production.json`: dp_-/dn (fresh, relaxed) on n = 0-62 for 57 H values; p_-,
  Lambda_1.
- `code/relic.py`: counts, super-horizon walls, relic integral, late-domain mass, threshold finder.
- `code/thresholds.py` → `results/thresholds.json` (all scenarios, systematics, DeltaN scan, closed forms, critic
  comparison, N_0 reach), `results/per_H_quantities.json`, `results/spectra.json`.
- `code/headline.py` → `results/headline.json`: summary table with systematic ranges and values at fixed H.
- `code/checks_mr.py` → `results/checks.json`: 32 checks (31 pass; the failure is explained in §11).
- `code/run_all.py` → `results/runtimes.json`, `results/run_all.log`.
