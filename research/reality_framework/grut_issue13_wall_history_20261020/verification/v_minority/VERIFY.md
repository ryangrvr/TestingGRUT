# Verification v_minority: module minority_relic

21 October 2026. Adversarial check of `modules/minority_relic` (NOTES.md, results/*.json). I wrote all the code in
this directory from scratch before I opened the module code. I did not import or copy the module's code or the
inflation_spectator Fokker-Planck chain: the 1D Fokker-Planck solver, the RG running, the g_* model and the
Nambu-Goto ODE are all my own. Only after I had my numbers did I read `mr_common.py`, `thresholds.py`
(`Nstar_reach`), the `collapse.py` header and the module's JSON files, to locate discrepancies. For one diagnostic
(the source of a 0.8 % T_RH offset) I evaluated cosmo_observable's `core.gstar_rho/gstar_s` once. None of my
numbers uses it.

Inflation, H_inf, N_*, the start condition, reheating and the domain conventions are outside-action preparation data
[P], as in the module. Python 3, numpy and scipy only, OMP_NUM_THREADS=2. Every job ran in under 80 s.

Labels: [V] computed here; [D] derived here; [Q] quoted.

## Files

| file | content |
|---|---|
| `vm_collapse.py` → `results_collapse.json` | NG sphere in a ∝ tau^alpha. Two independent integrations: (A) (r, u = gamma r') then (tau, eps = r^2 gamma)(r) for the plunge; (B) brute second-order (r, r') to r = 1e-3. Flat-space check against Gamma(1/4)^2/(4 sqrt(2 pi)). |
| `vm_cosmo.py` → `results_cosmo.json` | Own SM g_rho, g_s (massive ideal gases, QGP/pion blend, instantaneous nu decoupling at 2 MeV). Instant reheating T_RH, N_0, n'(T), T_dom,iso, relic bound, re-entry examples, N_0 reach. |
| `vm_fp.py` → `results_fp.json` | Own one-loop running (generic real-scalar tensor betas for lambda, m^2, xi) to mu = H. Own 1D FP on y = phi/H with the exact 1D steady conductance c = D/∫e^s on each link, absorbing at y = 0, exact time evolution by eigendecomposition, and Lambda_1 from the MFPT formula in log space. Relic integral and threshold root-finding. Regime A. |
| `vm_sensitivity.py` → `results_sensitivity.json` | Where the relic comes from (wall/radiation ratio at re-entry), the threshold with a gravity cut, r(tau = r0/2), regime A for extreme N_0, N_0 reach with the g_* term. |
| `vm_extra.py` → `results_extra.json` | Count closed forms, relic at 241 GeV, closed-form relic threshold. |
| `results_summary.json` | Claim, module value, my value, verdict. |

## 1. Collapse energy factor kappa_E, tau_c

**Method [D, V].** S = -4 pi sigma ∫ dtau a^3 r^2 sqrt(1 - r'^2) gives
r'' = -2(1 - r'^2)/r - 3(a'/a) r'(1 - r'^2) (derived independently; same as the module). The conformal Hamiltonian
H_c = 4 pi sigma a^3 r^2 gamma has explicit time dependence. That gives d(r^2 gamma)/dtau = -3(a'/a) r^2 gamma r'^2,
which is regular at r → 0. Initial state: at rest at tau_i = r0/x. Re-entry: a r0 = 1/H, so tau_re = alpha r0.

**Released energy, definition.** The module (collapse.py header and NOTES §3) uses E_c = 4 pi sigma a_c^2 lim r^2 gamma,
the wall energy at r → 0. There it is all kinetic: the rest part 4 pi sigma a^2 r^2 → 0. kappa_E = E_c/(sigma R_re^2)
with R_re = 1/H_re. This is the physically relevant energy converted into debris, and I use the same definition.

**Numbers [V]** (x = 1e4; x = 1e2 to 1e4 change kappa_E by 2e-4; brute method B agrees to 1e-3):

| | module | mine |
|---|---:|---:|
| tau_c/r0 (RD) | 1.8670 | 1.866979 |
| eps_c/r0^2 | 0.1778 | 0.177771 |
| kappa_E (RD) | 7.787 | **7.78661** |
| kappa_E/(4 pi) | 0.620 | 0.61964 |
| r(tau_re)/r0 | 0.734 | 0.73364 |
| t_c/t_re | 3.486 | 3.4856 |
| max wall energy /(4 pi sigma R_re^2) | 0.774 at 1.34 r0 | 0.77396 at 1.3404 r0 |
| kappa_E (MD), tau_c (MD) | 1.888, 2.2663 | 1.88810, 2.266323 |
| flat tau_c | 1.3110287771 | 1.31102877715 (exact to 1e-13), eps conserved to 2e-12 |
| C_E = 3 kappa_E/(4 pi) | 1.859 | 1.85892 |

**Alternative definitions [V].** Wall energy at tau_re: 0.653 x 4 pi sigma R_re^2. Maximum ("collapse onset"):
0.774. At r = r0/2: 0.774. Using the onset energy instead would raise kappa_E by 25 % and move the relic threshold
by about -1.6 % (d ln rho/d ln H ≈ 14 at threshold). The statement "0.62 of the static sphere energy at re-entry"
is correct for the energy at collapse.

Verdict: **CONFIRMED** (kappa_E = 7.787, tau_c = 1.867 r0, 0.62).
Wording, minor: "the wall starts to fall at tau ~ r0/2" is loose. r' = -0.1 already at tau = 0.199 r0, and
r(r0/2) = 0.937 r0. The ODE handles the motion correctly, so this does not affect any number.

## 2. Re-entry map, instant reheating, N_0

**Method [D, V].** a_end = a(T_RH) with 3 M_P^2 H_inf^2 = (pi^2/30) g_rho T_RH^4. By entropy conservation
a/a0 = (g_s0/g_s)^{1/3} T0/T. N_0 = ln(a_end H_inf/(a0 H0)), h = 0.6736. n'(T) = ln(a H/(a0 H0)), with H including
matter. In RD, aH ∝ g_rho^{1/2} g_s^{-1/3} T, so T_re = T_RH e^{-m} [g_rho^{1/2} g_s^{-1/3}](T_RH)/[...](T_re)
(c_R = 1).

| quantity | module | mine |
|---|---|---|
| T_RH at 130/198/241/302/500/2000 GeV | 0.97/1.19/1.32/1.47/1.89/3.79 e10 | 0.962/1.187/1.310/1.466/1.887/3.773 e10 |
| N_0 at 130/150/300/1000/2000 GeV | 48.05/48.115/48.462/49.064/49.41 | 48.045/48.117/48.463/49.065/49.412 |
| g factor high/low | 1.87 | 1.870 |
| T_re at 250 GeV, m = 10/30/40 | 6.09e5 GeV / 1.83 MeV / 106 eV | 6.06e5 GeV / 1.830 MeV / 105.6 eV |
| n'(T_rec, T_eq, 6.37 keV, T_dom,iso, 70 keV, 1 MeV) | 3.06, 3.83, 12.47, 12.69, 14.88, 17.77 | 3.058, 3.834, 12.468, 12.694, 14.879, 17.765 |
| T_dom,iso; sigma H/s there | 7.972 keV; 1.716e-6 GeV | 7.979 keV; 1.7166e-6 GeV |
| DeltaN = N_* - N_0 | 1.0-1.9 (N50), 10.6-11.9 (N60) | 0.59-1.95, 10.59-11.95 (130-2000 GeV) |

The T_RH offset (≤ 0.8 %) is located. The module takes the Saikawa-Shirai fit of cosmo_observable, which gives
g_rho(1e10 GeV) = 104.9; my ideal-gas SM gives 106.75. This is a convention, not an error. N_0 is insensitive to it
(≤ 0.005).

Verdict: **CONFIRMED**.

## 3. Relic bound and the threshold (relaxed start, N_* = 50, f = 1)

**Dimensional form [D].** A domain of re-entry radius R = 1/H_re takes up volume fraction dp per e-fold. So
n_dom = 3 dp/(4 pi R^3) and E_c = kappa_E c_A sigma R^2. Hence
rho/s = f (3 kappa_E c_A/4 pi) sigma H_re (dp/dn)/s_re dn', in GeV. dn = dn' because n' = n - DeltaN - ln c_R.
Y_dom = n_dom/s is conserved, so evaluating at re-entry or at collapse is the same. The bound is
0.12 rho_c/(h^2 s0) = **4.3733e-10 GeV** (g_s0 = 3.909, T0 = 2.7255 K), which is h-independent. **CONFIRMED.**
C_E = 1.859 c_A: **CONFIRMED**. The CHC-3 mapping C_E = 2 c_b, dp/dn = 2 p_-(per ln t) follows from d ln t = 2|dn'|
and H = 1/(2t) in RD: **CONFIRMED** (algebra only).

**Own FP [V].** V = 1/2 (m_phi^2(H) + 12 xi_phi(H) H^2) phi^2 + lambda_phi(H) phi^4/4, with my own one-loop running
from mu_0 = 100 GeV. At 1 TeV: lambda = 0.0100322, xi_phi = -2.4397e-4; this matches the packet's xi table to all
printed digits. dp_-/dn = (absorbed flux at y = 0)/2, by the reflection principle. Relaxed start:
e^{-s} restricted to y > 0. Validation:
- fresh p_- at 335 GeV (N50/N60): 1.8197e-4 / 3.4018e-4 (module 1.8205e-4 / 3.4028e-4);
- at 1 TeV: 0.18381 / 0.20755; at 2 TeV: 0.33011 / 0.34626;
- Lambda_1 = 6.3158e-3 at 1 TeV (module 6.316e-3) and 2.1303e-3 at 501 GeV (2.130e-3);
- at 250 GeV, p_- = 2.404e-9 extrapolated in h (module 2.4148e-9, v_inflation 2.412e-9; 0.4 %).

The relaxed dp/dn shape at 241 GeV matches NOTES §1: log10 = -10.25 at n = 4 and -10.45 at n = 20.

**Thresholds (H_inf where rho/s = 4.373e-10 GeV, instant reheating, T_cut = T_dom,iso) [V]:**

| | module | mine |
|---|---:|---:|
| **relaxed N50, f = 1** | **362.5** | **362.35** (grid h = 0.008/0.004/0.002: rho/s changes 4e-5) |
| relaxed N60, f = 1 | 366.6 | 366.6 |
| fresh N50 / fresh N60, f = 1 | 505.4 / 437.8 | 505.0 / 437.9 |
| relaxed N50, f = 0.1 | 489.2 | 489.0 |
| closed form (dp/dn = Lambda_1/2): Lambda_1 ≤ 2.74e-4/f | 372.1 | 372.0 (Lambda_1 ≤ 2.741e-4) |
| tree potential, no running or xi | – | 358.6 (running + xi: +1 %) |

The relaxed transient matters at this low barrier (B = 4.0 at 362 GeV): dp/dn at n = 14.3 (T_dom) is 1.42 Lambda_1/2.
That is the source of the 362 vs 372 GeV difference.

Verdict: **CONFIRMED** (362 GeV).

**Where the relic comes from, and a missing systematic [V].** The relic is dominated by domains that re-enter just
above T_cut, where the integrand goes as 1/T_re. At the threshold, **70 % of rho/s comes from domains with
sigma/(M_P^2 H_re) ≥ 0.1 (T_re < 25 keV), and 91 % from ≥ 0.01 (T_re < 78 keV)**. In this window the wall energy is
10-100 % of the horizon radiation energy. The no-gravity NG kappa_E is uncontrolled there, and so is the assumption
that the debris becomes hidden quanta rather than black holes. If only domains with ratio < 0.1 are counted, the
relaxed N50 f = 1 threshold moves to **404 GeV**; with < 0.01 it moves to **507 GeV** (fresh N50, ratio < 0.1:
556 GeV). The module's systematic range (314-545 GeV, covering X, c_A, s_off and the 6.37 keV cut) does not include
this. The cut variant it does include, 6.37 keV, only tightens the bound. If those domains form PBHs instead, the
energy stays as non-relativistic matter, so a comparable abundance constraint may still apply. That was not
computed, and the direction of the correction is therefore open.

## 4. "N_* = 60 is unreachable with instant reheating for H_inf ~ 100s GeV"

**[V]** N_0(instant) = 48.05 (130 GeV) to 48.72 (500 GeV), and 49.41 at 2 TeV. So N_* = 60 means DeltaN = 10.6-12.0:
the fresh start lies 11-12 e-folds before today's Hubble scale exits. **CONFIRMED.**

**N_0 reach table: CORRECTED (minor).** The module uses N_0 ± (1/12) ln(rho_inf/rho_x)
(`thresholds.py:Nstar_reach`). That formula holds only at constant g. Done properly,
N_0,kin - N_0,inst = (1/12)[ln(rho_inf/rho_x) + ln(g_s(T_RH)/g_s(T_x))], and the same g term (+0.19) enters the MD
case. I checked this directly by integrating ln(aH) through a w = 1 era (59.604 at 1.6 TeV).

| H (GeV) | module kin 4 MeV | mine | module MD 10 MeV | mine | mine, kination to 1 MeV |
|---:|---:|---:|---:|---:|---:|
| 150 | 57.84 | 58.03 | 38.71 | 38.90 | 58.49 |
| 250 | 58.18 | 58.37 | 38.88 | 39.07 | 58.83 |
| 1000 | 59.10 | 59.29 | 39.34 | 39.53 | 59.76 |
| 1600 | 59.42 | 59.60 | 39.50 | 39.68 | 60.07 |

The conclusion "N_0 = 60 not reachable for H ≤ 1.6 TeV with any w ≤ 1 era ending above 4 MeV" survives, with a margin
of 0.4 e-fold at 1.6 TeV. If kination is allowed to last until 1 MeV, N_0 = 60 is reached at about 1.45 TeV. For
H ~ 100s GeV the claim is robust either way (≥ 1.2 e-folds of margin).

## 5. Regime-A boundary 145.6 GeV

**[V]** N(all sizes in today's Hubble volume) = ∫ (dp/dn) e^{3n'} dn' ≈ (Lambda_1/2) e^{3 N_0}/3. It is dominated by
the last e-folds, where both starts are quasi-stationary, so it is independent of N_* and of the start condition. I
get Lambda_1 by MFPT; it agrees with the eigenvalue to 0.03 % at 250-300 GeV.

| N_0 used | threshold (N < 1e-3) | module |
|---|---:|---:|
| instant reheating | **145.57 GeV** (Lambda_1 = 1.42e-65, B = 146.7) | 145.6 |
| exit convention N_0 = 60 | 137.83 | 137.8 (inflation_spectator 138.1) |
| exit convention N_0 = 50 | 144.19 | 144.2 |
| kination to 4 MeV (corrected N_0) | 139.0 | – |
| MD reheating to 10 MeV | 153.4 | – |
| tree potential, instant | 145.51 | – |

Verdict: **CONFIRMED**, for instant reheating as labelled. Over w ≤ 1 histories the boundary spans 139-153 GeV.

## 6. Other claims spot-checked

- Lambda_1/2 = 1e-11 at **237.29 GeV** (module 237.4): CONFIRMED.
- Count closed forms Lambda_1 < 6 e^{-3 n'_X}: T_dom 211.93, 1 MeV 192.40, T_rec 402.6 GeV (module 212.0, 192.4,
  403.9): CONFIRMED.
- **"sigma/(M_P^2 H) = 1.1e-5 at 1 MeV" (NOTES §4, §5): CORRECTED.** It is 3.6e-5 at 1 MeV; 1.1e-5 is the 1.8 MeV
  value (§7 table). This is text only and enters no threshold.
- **§7 "relic at 241 GeV = 10^-7.6 (relaxed), 10^-12.3/10^-9.6 (fresh)": CORRECTED.** These are the H = 240 GeV rows
  of headline.json. The module's own per_H_quantities at 241.04 GeV gives -7.48 / -12.18 / -9.46; mine is -7.488.
- §5: "if some [near-T_dom domains] form black holes, these are small and evaporate into SM states (lowering f)".
  This is an overclaim for the population that dominates the relic. At T_re ~ 8-80 keV a horizon-scale collapse gives
  roughly 1e5-1e9 M_sun objects, which do not evaporate. Black holes that evaporate by today need M ≲ 1e15 g, so
  re-entry above about 30 GeV.

## 7. Overclaims and gaps

1. The relic thresholds (362-367 relaxed, 438-521 fresh at f = 1; the f = 0.1 values) are quoted without the
   near-criticality systematic. 70 % (91 %) of the integral comes from domains with wall/radiation ≥ 0.1 (≥ 0.01),
   where no-gravity NG and the hidden-quanta fate are not established. Excluding them raises the relaxed threshold to
   404 (507) GeV. NOTES §3 marks the NG result as invalid only above ratio 1.
2. The PBH alternative for near-critical domains is not computed, and the "small black holes, lower f" text is wrong
   for them (item 6).
3. Minor numerical text errors: the N_0 reach is missing a +0.19 g_* term; the 1 MeV wall fraction is off by 3.3x; the
   241 GeV relic row is the 240 GeV row.
4. Regime A at 145.6 GeV holds for instant reheating only (139-153 GeV over w ≤ 1 histories). The module labels it as
   such.
5. Not checked here: late-domain mass thresholds (Omega_late), the full T_rec/T_eq count thresholds with transients
   (only the closed forms), the X-correction factors, the DeltaN scan, the CHC-3 numerical offsets, and the
   percolation inputs.

## Summary

The five requested items reproduce independently:
- kappa_E = 7.7866, tau_c = 1.86698 r0, E_c = 0.6196 x 4 pi sigma R_re^2 (energy at collapse);
- N_0 = 48.05-49.41;
- the relaxed N50 f = 1 relic threshold is 362.35 GeV;
- N_* = 60 is unreachable with instant reheating;
- regime A is at 145.57 GeV.

Three numerical text slips are corrected: the N_0 reach (+0.19), the 1 MeV wall fraction, and the 241/240 GeV relic
row. The main overclaim is robustness. The relic bound is set by near-critical domains, where the no-gravity NG
energy and the hidden-quanta fate are not controlled. A gravity-motivated cut moves the f = 1 relaxed threshold from
362 to 404-507 GeV. All results remain conditional history constraints on outside-action preparation data.
