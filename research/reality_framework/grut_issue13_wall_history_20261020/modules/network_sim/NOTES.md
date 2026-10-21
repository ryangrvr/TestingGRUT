# Module network_sim: Z2 wall population of the unchanged phi potential (X = D = 0), radiation era

21 October 2026. Calculation module of the issue-13 wall-history packet (BRIEF.md). This is a **classical lattice
(physical-equation and Press-Ryden-Spergel) simulation of the X = D = 0 sector only**. It does not include the
charged field X, the reporter D, thermal masses after formation, quantum fluctuations, wall gravity/backreaction, or
any term outside the unchanged action. It does not exclude or validate Candidate 1; it supplies conditional
statements about what the classical phi sector does once a restored-phase preparation is given.

Reproduce everything (about 2.5 CPU-hours; two single-threaded worker processes; each job capped at 24 min):

    bash code/run_all.sh
    # = validate.py ; run_networks.py run (128 lattice jobs) ; bubbles.py ; scales.py ; analyze.py

Existing `results/runs/*.json` are skipped; delete `results/runs/` to recompute. Every number quoted below is in
`results/summary.json` (aggregates), `results/runs/<name>.json` (time series of one run), `results/bubbles.json`,
`results/validation.json` or `results/scales.json`. Status labels: **[computed]** = output of the code here,
**[derived]** = closed-form from the action/conventions, **[quoted]** = literature, with read level in §9.

## 1. Headline

| Quantity | Value | Status / provenance |
|---|---|---|
| A (3D PRS, alpha=3, beta=0), 256^3, dx = delta/2, 3 seeds, window h = tau/L in [1/8,1/4] | 0.777 +- 0.013 (seed std); 384^3 (1 seed): 0.770 | computed; summary.json groups.prs3d_N256_L256.A_win |
| A (3D PRS), box-affected window h in [1/4,1/2] | 0.849 +- 0.043 (still rising) | computed; A_box |
| A (3D physical eq.), 256^3, L = 24, tau 2 -> 12, post-formation tau in [6.4,12] (h up to 1/2) | 0.884 +- 0.062 (3 seeds) | computed; A_postform (box-affected, H^-1/delta <= 72) |
| A (2D physical eq., z-invariant 3+1 walls), 2048^2, L = 64, window tau in [8,16] | 0.709 +- 0.014 (3 seeds) | computed |
| A (2D PRS alpha=2), 2048^2, L = 2048, window tau in [256,512] | 0.652 +- 0.002 (3 seeds); 4096^2, L = 4096: 0.725 | computed |
| A (2D PRS alpha=3 = damping of z-invariant 3+1 walls) | 0.673 +- 0.016 (3 seeds) | computed |
| Resolution checks, same continuum initial field, A(tau) rel. difference (late mean / max) | 3D PRS dx 1 vs 0.5: 0.9 % / 1.3 %; 3D phys 128^3 vs 256^3: 1.8 % / 2.8 %; 2D phys 2048 vs 4096: 0.4 % / 3.0 % (1024 vs 2048: 3.5 % / 6.6 %); 2D PRS dx 1 vs 0.5: 1.9 % / 3.5 % | computed; summary.json resolution_checks |
| Scaling exponent mu (area density ~ tau^mu; scaling = -1) | 3D PRS: -0.896 +- 0.032 (256^3, tau 25-62), -0.975 (384^3, tau 28-93); 2D PRS alpha=2: -0.991 +- 0.010, alpha=3: -0.988 +- 0.018; 2D phys: -0.91 +- 0.14 (tau 10-15.6 only) | computed |
| Fraction of wall area in box-spanning ("infinite") wall clusters, window h in [1/8,1/4] | 3D PRS 0.994 +- 0.001; 3D phys 0.996; 2D phys 0.71 +- 0.06; 2D PRS 0.81 +- 0.14 (alpha=2), 0.84 +- 0.14 (alpha=3) | computed |
| Persistence verdict p = 1/2 | Network formed and walls were present at every output up to tau = L/2 in all 23 unbiased runs; in every 3D run a box-spanning wall carried >= 97.4 % of the area at every post-formation output; in 2 of 14 2D runs no spanning wall remained before h = 1/4 (finite-box 2D critical percolation) | computed |
| Wall velocity v_rms (estimator sum pi^2 / sum grad^2 over abs(phi) < 1/2) | 3D PRS 0.456 +- 0.006; 2D PRS 0.487 +- 0.011; 2D phys 0.478 +- 0.015 | computed |
| Biased decay, 2D PRS alpha=2: coarsening factor K10 = xi_ref(tau_10)/xi_ref(tau_form) | K10 = 0.262 (1/2 - p)^(-1.255) (early formation); 0.247 (1/2 - p)^(-1.226) (formation deep inside horizon, damped preparation); K10 ~ 9-11 at p = 0.45, 4.3-4.9 at 0.40, 2.9-3.0 at 0.35 | computed; fits over 0.05 <= 1/2 - p <= 0.3 only |
| Biased decay, 3D PRS (128^3) | tau_10 = 40.7, 22.3, 16.3 (p = 0.45, 0.40, 0.35; formation 12.3); K10 = 2.73, 1.63, 1.27; p = 0.45: only 1 of 3 seeds wall-free by tau = 64 = L/2 | computed |
| Biased, 2D physical eq. | p = 0.25, 0.20: no network forms (minority gone before formation, 3/3 seeds); p = 0.30-0.35: all walls gone by tau ~ 10-12 | computed |
| Closed-bubble collapse (thin-wall NG, radiation era, from rest at tau_i), x = R0/tau_i | x << 1: Delta tau_c = 1.5708 R0 (circle/cylinder), 1.3110 R0 (sphere); x >> 1: tau_c -> 2.459 R0 (cylinder, D_h = 3), 2.214 R0 (circle, D_h = 2), 1.867 R0 (sphere, D_h = 3) | derived/computed (ODE); bubbles.json |
| Field theory vs NG collapse time tau_10 (R = R0/10) | physical radial: max abs dev 2.5 % (22 runs); PRS radial: max 4.0 % (21 runs); Cartesian 2D/3D vs radial: <= 2.1 % | computed |
| Super-horizon freezing (field, physical eq., x = 8) | R(tau = R0)/R0 = 0.875 (cylinder), 0.744 (sphere); R(tau = R0/2)/R0 = 0.972, 0.945 | computed |
| Wall-domination time if an unbiased scaling network persisted with the tree tension | t_w = 3 M_P^2/(4 A sigma) = 3.9e4 s (A = 0.8; 6.2e4 s for A = 0.5) | derived (conditional; validity limit of these no-backreaction runs, not a prediction) |

Interpretation in one line: **given an unbiased restored phase, the classical X = D = 0 sector forms a percolating
wall network that persists and approaches scaling with A ~ 0.8 in 3D; given an initial population bias (an
outside-action preparation datum), the network disappears after a finite, p-dependent amount of coarsening; isolated
closed walls always collapse within ~1.3-2.5 comoving radii of conformal time after release/horizon re-entry.** No
statement about which preparation Candidate 1's cosmology actually supplies is made here.

## 2. Equations, conventions and derivations

### 2.1 Units and the wall

Sector X = D = 0 of the unchanged potential: U = (lambda_phi/4)(phi^2 - v^2)^2. Field in units of v, lengths and times
in units of 1/m_phi, m_phi^2 = 2 lambda_phi v^2. Then lambda = 1/2,

    U(phi) = (1/8)(phi^2 - 1)^2,   U'(phi) = (1/2) phi (phi^2 - 1),   U''(+-1) = 1.

Planar wall phi = tanh(z/delta), delta = sqrt(2/lambda_phi)/v = 2/m_phi; tension sigma = int dz phi_z^2 =
(2/3) m_phi v^2 = (2 sqrt(2 lambda_phi)/3) v^3. Benchmark [derived, scales.json]: m_phi = 141.42136 GeV,
sigma = 9.4280904e7 GeV^3, sigma^(1/3) = 455.136 GeV, delta = 0.0141421 GeV^-1 = 2.79e-18 m. The dynamics depend
only on the dimensionless ratios; sigma enters only through rho_w = A sigma / t.

### 2.2 Radiation-era field equation

ds^2 = a^2(tau)(-dtau^2 + dx^2). From S = int sqrt(-g)[-1/2 (dphi)^2 - U]:

    phi'' + 2 (a'/a) phi' - lap(phi) + a^2 U'(phi) = 0.

Radiation era a proportional to tau; we set a = tau (tau in 1/m_phi; a = 1 at tau = 1). Then a'/a = 1/tau,
H = a'/a^2 = 1/tau^2, comoving Hubble radius (aH)^-1 = tau, comoving wall half-width delta/a = 2/tau, cosmic time
t = int a dtau = tau^2/2 = a tau/2, and H^-1/delta = tau^2/2. The background is supplied (no wall backreaction).
The 2D "physical" runs solve this 3+1 equation for fields independent of z (straight-in-z walls; area per volume =
line length per area).

### 2.3 PRS family and the wall damping coefficient

PRS (ApJ 347 (1989) 590) use phi'' + alpha (dln a/dln tau) phi'/tau - lap(phi) + a^beta U'(phi) = 0, which follows
from L = a^alpha [phi'^2/2 - (grad phi)^2/2 - a^beta U]. A static kink then has comoving width proportional to
a^(-beta/2) and comoving energy per comoving area proportional to a^(alpha + beta/2), so the thin-wall action is
-sigma_0 int d^(p+1)xi a^(D_h) sqrt(-gamma_flat) with

    D_h = alpha + beta/2.

For a Nambu-Goto p-brane the induced metric gives a^(p+1), so the physical damping is D_h = p + 1.
Physical 3+1 equation: alpha = 2, beta = 2, D_h = 3 (p = 2). PRS 3D: alpha = 3, beta = 0, D_h = 3 (standard).
PRS 2D with the literature choice alpha = 2 gives D_h = 2 (line walls of a 2+1 universe), whereas our 2D physical
runs (z-invariant 3+1 walls) have D_h = 3; the matching 2D PRS choice is alpha = 3. Both 2D PRS values were run.
Discretisation (PRS 1989): delta_n = alpha dt/(2 tau_n), pi_{n+1/2} = [(1 - delta_n) pi_{n-1/2} + dt(lap phi_n -
a_n^beta U'(phi_n))]/(1 + delta_n), phi_{n+1} = phi_n + dt pi_{n+1/2}; (2d+1)-point Laplacian; periodic box;
dt = dx/3.

### 2.4 Area parameter and estimator

Comoving area A_c in comoving volume V_c: rho_w = sigma a^2 A_c/(a^3 V_c) = sigma A_c/(a V_c). With t = a tau/2:

    A = rho_w t / sigma = (A_c/V_c) tau / 2 = tau / (2 xi_c),   xi_c = V_c/A_c.

This is the Hiramatsu-Kawasaki-Saikawa definition (their eqs. 4, 12, 13). **Convention note:** "comoving area
density times conformal time" equals 2A in radiation era, not A. Estimator (stated): PRS link count with gradient
weight, A_c = dx^(d-1) sum_{sign-changing links} |grad phi| / sum_i |d_i phi| (link-direction derivative = link
difference; transverse derivatives = mean of central differences at the two ends, as in Garagounis-Hindmarsh).
Also recorded: raw count x isotropic factor (pi/4 in 2D, 2/3 in 3D), `A_iso`; it agrees with A to <= 1.2 % at every
output in the main windows.

Windows: h = tau/L. Main window h in [1/8,1/4] (box at least 4 comoving horizons across). Box-affected window
h in [1/4,1/2] (HKS also stop at h = 1/2 and note a late rise from the box). Slope fit over
[max(L/32, 2 tau_form), L/4]. Formation time tau_form = first output where the volume fraction with abs(phi) < 1/2
outside walls (fraction minus 2 atanh(1/2) delta_com A_c/V_c) is < 0.05 and both vacua are reached.

### 2.5 Thin-wall collapse of an isolated closed wall

Comoving radius R(tau) of a circle/cylinder (kappa = 1) or sphere (kappa = 2): L = -a^(D_h) R^kappa sqrt(1 - R'^2),
whose Euler-Lagrange equation is

    R'' = -(1 - R'^2) [kappa/R + D_h (a'/a) R'].

Flat space from rest: R'^2 = 1 - (R/R0)^(2 kappa), collapse time int_0^R0 dR/sqrt(1 - (R/R0)^(2 kappa)) =
(pi/2) R0 (kappa = 1) and R0 int_0^1 dx/sqrt(1 - x^4) = 1.3110288 R0 (kappa = 2). Radiation era: the equation is
invariant under (tau, R) -> (s tau, s R), so tau_c/tau_i = F(x), x = R0/tau_i = comoving radius over comoving
Hubble radius at release. For x >> 1, R' ~ -kappa tau/((D_h + 1) R0) at early times (frozen), and tau_c -> c_inf R0.
In cosmic time a super-horizon sphere is gone at t_c = c_inf^2 t_R = 3.49 t_R, t_R = time at which the comoving Hubble
radius equals R0.

### 2.6 Unit translation

H^-1/delta = tau^2/2: 512 at the end of the 2D physical runs, 72 for 3D physical runs [derived]. At an illustrative
formation temperature T = 1.5476 v with g_* = 106.75 the real hierarchy is H^-1/delta = 2.1e13 [derived, scales.json;
T and g_* are free/illustrative inputs, T is the leading curvature-zero CHECK of the brief, not a transition result].
The extrapolation from simulated to physical hierarchy rests on the scaling attractor being independent of
delta/H^-1, which the PRS-vs-physical agreement supports but does not prove.

## 3. Configurations

Initial condition (all runs): Gaussian random field, white noise smoothed by a Gaussian kernel of std s_c
(C(r) proportional to exp(-r^2/(4 s_c^2))), zero mean, std 0.1 v, zero velocity, deterministic seeds 1, 2, 3.
Resolution pairs share the same continuum field (white noise drawn on a base grid N_base and band-limited
zero-padding; validated to 1.7e-16). Biased runs: constant shift chosen so the positive-volume fraction is exactly p
on the lattice.

| Group | Eq. | d | N | L | dx | tau_i -> tau_end | s_c | seeds | min 2delta/(a dx) |
|---|---|---|---|---|---|---|---|---|---|
| phys2d_N2048_L64 | alpha=2, beta=2 | 2 | 2048 | 64 | 1/32 | 2 -> 32 | 0.25 | 3 | 4.0 |
| phys2d_N1024 / N4096_L64 | same | 2 | 1024 / 4096 | 64 | 1/16, 1/64 | 2 -> 32 | 0.25 | 1 | 2.0 / 8.0 |
| phys3d_N256_L24 (+N128) | alpha=2, beta=2 | 3 | 256 (128) | 24 | 0.094 (0.19) | 2 -> 12 | 0.3 | 3 (1) | 3.56 (1.78) |
| prs2d_a2 / a3_N2048_L2048 | beta=0, alpha=2/3 | 2 | 2048 | 2048 | 1 | 2 -> 1024 | 1.5 | 3 + 3 | 4 |
| prs2d_a2_N4096_L2048 / L4096 | beta=0, alpha=2 | 2 | 4096 | 2048 / 4096 | 0.5 / 1 | 2 -> 1024 / 2048 | 1.5 | 1 + 1 | 8 / 4 |
| prs2d_a2_N2048_L2048_late | beta=0, alpha=2 | 2 | 2048 | 2048 | 1 | 20 -> 1024 | 1.5 | 1 | 4 |
| prs3d_N256_L256 (+N384_L384, N256_L128) | beta=0, alpha=3 | 3 | 256 (384, 256) | 256 (384, 128) | 1 (1, 0.5) | 2 -> L/2 | 1.5 | 3 (1, 1) | 4 (4, 8) |
| phys2d_bias_N1024_L32 | phys | 2 | 1024 | 32 | 1/32 | 2 -> 16 | 0.25 | 3 x 7 p | 4 |
| prs2d_bias_a2_N1024_L1024 | PRS a2 | 2 | 1024 | 1024 | 1 | 2 -> 512 | 1.5 | 3 x 7 p | 4 |
| prs2d_biaslate (failed setup) / biasdamp | PRS a2 | 2 | 1024 | 1024 | 1 | 100 -> 512 (biasdamp: extra friction 0.5 phi' for tau < 120) | 1.5 | 3 x 7 p | 4 |
| prs3d_bias_N128_L128 | PRS a3 | 3 | 128 | 128 | 1 | 2 -> 64 | 1.5 | 3 x 7 p | 4 |

"2delta/(a dx)" is the number of cells across the core abs(phi) < tanh(1) = 0.76 of a planar wall at the end (physical
runs) or always (PRS). HKS required delta/(a dx) > 1 (i.e. > 2 here); all main runs have >= 3.56, the low-resolution
partners 1.78-2.0. p in {0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2}.

## 4. Validation [computed, results/validation.json]

- Area estimator on exact tanh circles/spheres (R = 8-24) and tilted periodic planes (k = (1,0), (1,1), (1,2),
  (1,2,3)): max abs error 3.6 % (sphere R = 8 at dx = delta/2), <= 1.0 % for R >= 16; planes <= 0.06 %. The raw
  isotropic count errs by up to 33 % for single orientations (it is only an average correction).
- Periodic labelling: stripe wrapping, boundary-crossing blob, corner connectivity in 3D: pass.
- IC embedding: max difference on common points 1.7e-16.
- Integrator energy budget E(tau) - E(tau_i) - int dE/dtau|_{friction, expansion} (float64, small boxes):
  physical 2D 1.4e-3 of E (production dt), 3.4e-4 at dt/2; 3D PRS 1.7 % of instantaneous E (0.66 % of max E),
  0.44 % at dt/2; **2D PRS: 21 % of the instantaneous E at late times (1.8 % of max E), 5.3 % at dt/2.** The residual
  converges as dt^2 (factor 4.0-4.1 per halving), so it is time-discretisation error dominated by the
  post-formation scalar oscillations (omega dt = 1/3), not a code error; it nevertheless fails a strict 1 %
  instantaneous-energy criterion for 2D PRS at the production step. The wall observable A changes by <= 3.5 % when
  dx and dt are both halved (resolution pair), which is the operative check.

## 5. Results

### 5.1 Unbiased formation, persistence and approach to scaling

All 23 unbiased runs (p = 1/2) formed walls (tau_form = 4.9 physical 2D, 5.1 physical 3D, 10-14 PRS with tau_i = 2, 29 for tau_i = 20) and had sign
changes at every output to tau = L/2. Time series are in `results/runs/*.json` (`A`, `area_density`,
`n_domains`, `spanning_cluster_area_frac`, `largest_cluster_area_frac`, `closed_wall_area_frac`, `v2_wall`, energies).

3D: after formation one box-spanning wall cluster carries 99.4 +- 0.1 % (PRS 256^3) of the wall area; closed walls
(non-spanning clusters) carry ~0.6 %; the number of sign domains is 125 +- 14 in the window (256^3), most of them
small closed islands inside the two percolating phases. A rises slowly through the window (3D PRS 256^3:
0.77 at tau = 32-64, 0.85 at 64-128); mu = -0.90 +- 0.03 (256^3) and -0.975 (384^3) - consistent with a slow
approach to scaling and with the literature, not a clean plateau.

Physical 3D: the dynamic range is small (tau 2 -> 12, H^-1/delta <= 72). The main window h in [1/8,1/4] = tau 3-6
straddles formation (A = 1.04, v_rms = 0.99 there are rolling-field values, not wall properties): **this window
fails as a scaling measurement**. After formation A = 0.82, 0.94, 0.89 (seeds) over tau 6.4-12, but this interval is
already in the box-affected range h > 1/4. No slope fit was possible.

2D: A = 0.65-0.73 (PRS) and 0.71 (physical) with mu = -0.99 (PRS). Box-spanning walls carry 66-94 % of the area;
in 2 of 14 2D runs (phys2d s2, prs2d late-start) all walls had become closed before h = 1/4 and one phase began
to take over the box (f_plus -> 0.95 at tau = L/2 in the late-start run). 2D at p = 1/2 sits exactly at the
continuum percolation threshold, so this is expected finite-box behaviour; the 3D runs never showed it.

Approach from either side: early formation (tau_i = 2) approaches A from below (2D PRS: 0.47 at tau = 10 -> 0.65-0.73);
formation deep inside the horizon (damped late start, tau_i = 100, horizon/s_c ~ 70) approaches from above
(A = 3.1 at tau = 123 -> 1.05-1.19 at tau ~ 256 = L/4, 3 seeds) and **has not converged by the box limit**.
IC sensitivity tau_i = 20 vs 2 (same seed): +7.2 % in A at h <= 1/4.

### 5.2 Literature comparison (read levels in §9)

| Source | Their quantity | Converted | Ours |
|---|---|---|---|
| HKS 1309.5001 (3D physical, RD, 512^3/1024^3, tau 2 -> 40/60) | A = 0.8 +- 0.1 (same definition) | - | 3D PRS 0.777 +- 0.013 (box-safe), 0.849 (h <= 1/2); 3D phys 0.88 +- 0.06 (h <= 1/2, small range): consistent |
| Leite-Martins 1110.3486 (3D PRS alpha=3, RD, 512^3/1024^3) | xi_c/eta = 0.60 +- 0.05 / 0.64 +- 0.03; gamma v = 0.46-0.48 | A = 1/(2 xi_c/eta) = 0.83 / 0.78 if xi_c = V/A_c (they write A/V proportional to 1/xi_c, so indicative) | xi_c/tau = 0.644 +- 0.010 (256^3), 0.649 (384^3); v_rms = 0.456 (gamma v ~ 0.51; estimators differ) |
| Leite-Martins 2D (matter era only) | xi_c/eta = 0.58-0.91 | - | not comparable (radiation era here) |
| Oliveira-Martins-Avelino hep-ph/0410356 | rho_w ~ eta^(lambda-1): 2D rad 8192^2 lambda = 0.0285 +- 0.0067; 3D rad 512^3 lambda = 0.0343 +- 0.0023 | mu = -0.97 / -0.966 | 2D PRS mu = -0.99; 3D mu = -0.90 (256^3) / -0.975 (384^3) |
| Garagounis-Hindmarsh hep-ph/0212359 (2D and 3D, RD) | comoving area density = 0.928 t^-0.996 (2D), 0.925 t^-0.994 (3D), t conformal | 2A = 0.93 -> A = 0.46 | **not reconciled**: factor ~1.7 below HKS and below ours; their conformally rescaled, initially dissipated equation and time origin were not reproduced |

### 5.3 Biased initial populations (p = initial positive-volume fraction)

tau_10 = first time A < 0.1 A_ref(tau) (A_ref = seed mean at p = 1/2, same group); tau_0 = first output with no
sign change (an upper bound on wall lifetime: residual scalar oscillations can create non-wall sign changes,
identified by v2_wall > 1); K10 = xi_ref(tau_10)/xi_ref(tau_form), the factor by which the unbiased network's length
scale grows between formation and tau_10. Seed means (std), 3 seeds each.

2D PRS alpha=2, early formation (tau_form = 10):

| p | 0.45 | 0.40 | 0.35 | 0.30 | 0.25 | 0.20 |
|---|---|---|---|---|---|---|
| tau_10 | 155.9 (16.6) | 63.9 (6.7) | 37.3 (4.0) | 24.3 (0) | 16.1 (1.5) | 12.0 (0) |
| tau_0 | 250.3 | 138.3 (13.9) | 78.0 (21.9) | 53.3 (5.8) | 44.7 (4.6) | 33.1 (3.3) |
| K10 | 10.8 (1.0) | 4.85 (0.48) | 2.99 (0.28) | 2.06 | 1.46 (0.12) | 1.13 |

2D PRS alpha=2, formation deep inside the horizon (tau_i = 100, damped preparation to tau = 120, tau_form = 111):
tau_10 = 258.3 (13.5), 169.3, 144.7, 123.3, 117.0, 111.0; K10 = 9.24 (0.71), 4.32, 2.85, 1.64, 1.31, 1.00.
Fits over the tested range: K10 = 0.262 (1/2 - p)^-1.255 (early), 0.247 (1/2 - p)^-1.226 (damped late).
The two preparations give the same K10(p) within 5-20 % (p = 0.45-0.25), although their ratio of comoving horizon
to initial correlation length at formation differs by a factor ~11 (10/1.5 vs 111/1.5): K10, not tau_10, is the
transferable quantity. (Observation only: the exponent is close
to the 2D percolation exponent nu = 4/3, as expected if the lifetime is set by the largest minority clusters; not
derived here.)

3D PRS alpha=3 (128^3, tau_form = 12.3): tau_10 = 40.7, 22.3, 16.3 for p = 0.45, 0.40, 0.35 and = tau_form for
p <= 0.30; K10 = 2.73, 1.63, 1.27; tau_0 = 47.3, 36.9, 26.0, 23.6, 17.0 for p = 0.40 ... 0.20; at p = 0.45 only 1 of 3
seeds was wall-free at tau = 64 = L/2 (box-limited). The minority phase still percolated after formation for
p >= 0.40 (and initially for all p >= 0.20), yet the network decayed: percolation of the minority at formation does
not guarantee persistence. Fit K10 = 0.331 (1/2 - p)^-0.70 (3 points only).

2D physical equation (L = 32, tau <= 16): p = 0.45: 1 of 3 seeds wall-free by tau = 16; p = 0.40: 2 of 3 by
tau = 10.5 +- 0.6; p = 0.35, 0.30: all gone by tau = 11.5, 10.4; p = 0.25, 0.20: the minority vanished before
formation (tau_0 = 4.4, 3.6 < tau_form = 5.0): no network formed.

### 5.4 Isolated closed bubbles [computed, results/bubbles.json]

NG ODE (validated: flat-space collapse 1.5707962 R0 and 1.3110287 R0 vs exact 1.5707963, 1.3110288). Radiation-era
tau_c/tau_i at x = R0/tau_i = 0.001, 0.178, 1, 5.6, 31.6, 1000: sphere (D_h = 3) 1.001, 1.238, 2.454, 10.66, 59.07,
1867; cylinder (D_h = 3) 1.002, 1.292, 2.917, 13.94, 77.78, 2459; circle (D_h = 2) 1.002, 1.288, 2.787, 12.62, 70.04,
2214. Law: tau_c = tau_i + c_flat R0 for R0 << tau_i (sub-horizon bubbles collapse in a light-crossing time) and
tau_c = c_inf R0 for R0 >> tau_i with c_inf = 1.867 (sphere), 2.459 (cylinder), 2.214 (2D circle, D_h = 2): the
bubble is frozen in comoving size until the horizon reaches it, then collapses. NG freezing: R(tau = R0/2)/R0 =
0.937 (sphere), 0.968 (cylinder); R(tau = R0)/R0 = 0.734, 0.870.

Field theory (radial finite-volume solver, physical equation, tau_i = 4 and 8, x = 0.125-8, 5 cells per delta_com
at the end; PRS radial solver tau_i = 50, x = 0.2-16): tau_10 agrees with NG to max 2.5 % (physical, mean +0.26 %)
and 4.0 % (PRS, mean +1.2 %). Radial resolution (2.5 vs 10 cells/delta): tau_10 20.07 vs 19.99 (cylinder), 15.82 vs
15.71 (sphere), NG 19.96, 15.68. Cartesian lattice checks with the network integrator (2D R0 = 4, 16; 3D R0 = 2, 4):
within 0.6-2.1 % of the radial solver. Super-horizon freezing in the field: x = 8: R(R0)/R0 = 0.875 (cylinder),
0.744 (sphere); R(R0/2)/R0 = 0.972, 0.945.
After the first collapse a transient opposite-vacuum region re-forms at the centre (debris/bounce) with radius up to
0.43 R0 (small x) and lives at most 9.2 R0 (sphere, x = 0.125, R0 = 1); in every run without numerical failure it
had vanished by the end (no negative cells). Cartesian 3D: last negative cell at tau = 11.8 (R0 = 4,
tau_c^NG = 9.8) and 8.8 (R0 = 2, tau_c^NG = 6.8).

## 6. Conditional history constraints (for the packet)

C1 (unbiased preparation). If the hidden Z2 is restored and the transition starts from a statistically unbiased
restored phase (p = 1/2; what a Z2-symmetric bath gives), the classical X = D = 0 sector forms a macroscopic wall
network with an infinite (box-spanning) wall carrying 99.4 % of the area in 3D (>= 97.4 % at every post-formation output), and it persists and approaches
scaling with A ~ 0.78-0.88 (3D; HKS 0.8 +- 0.1). Isolated-bubble collapse is not the fate of this network. If this
scaling persisted with the tree tension, walls would reach radiation density at t_w = 3 M_P^2/(4 A sigma) ~ 4e4 s
(A = 0.8) [derived, conditional]; beyond that the no-backreaction treatment here is invalid. This module does not
show that the hidden sector is ever restored, nor at what temperature, nor with which T_hidden/T_visible.

C2 (biased preparation; outside-action datum). If the initial vacuum population is biased, |p - 1/2| >= 0.05, the
network disappeared in every tested case after its length scale grew by K10 <= 11 (2D) or <= 2.7 (3D, box-limited)
relative to formation; at |p - 1/2| >= 0.25 no network formed in the physical 2D runs. For a formation hierarchy
H^-1/xi_f ~ 1e13 this occurs long before the network reaches horizon size. The behaviour for |p - 1/2| < 0.05 is not
tested; the fitted K10(p) laws are not extrapolated here. Any such bias requires a preparation outside the unchanged
action (e.g. an inflationary mean field), which must be supplied and justified separately.

C3 (closed regions). A closed opposite-vacuum region of comoving radius R collapses by conformal time ~1.87 R
(sphere) after release from rest, i.e. within ~3.5 times the cosmic time at which the Hubble radius reaches it.
Only a region larger than ~ the present comoving Hubble radius/1.87 could survive to today; that is a single
horizon-scale domain, not a network.

## 7. Free inputs the results depend on

- p, the initial positive-volume fraction (outside-action preparation datum).
- The restored-phase initial state: Gaussian spectrum (s_c = 0.25 comoving at tau_i = 2 physical; 1.5 = 0.75 delta in
  PRS), amplitude 0.1 v, zero velocity, instantaneous quench to the T = 0 potential; formation time relative to the
  horizon (tau_i = 2, 20, 100). A time-dependent thermal mass, latent heat and plasma friction are not modelled (the
  damped late start is a numerical preparation device, not physics of the hidden bath).
- Radiation-dominated background a proportional to tau, no wall backreaction; X = D = 0 throughout.
- PRS modification (beta = 0) for the large-range runs; alpha choice in 2D (2 vs 3).
- v and lambda_phi set the units only; sigma = 9.428e7 GeV^3 enters only the derived t_w.
- Illustrative cosmological inputs in scales.json (g_* = 106.75/110.75/3.36 and T = 1.5476 v, 1 TeV, 100 GeV) are not
  results.

## 8. What failed, limitations, not done

1. Physical 3D: box-safe window lies on formation; A there (1.04, v_rms ~ 1) is not a wall measurement; post-formation
   A only available in the box-affected interval; no slope fit. Physical 3D dynamic range (H^-1/delta <= 72) is far
   below HKS's tau range.
2. Undamped late-start biased runs (prs2d_biaslate, tau_i = 100): the weakly damped post-quench oscillations
   (amplitude ~1 about the vacuum, alpha/tau ~ 0.02) produce non-wall sign changes; their tau_10/tau_0 (e.g. 188 for
   all p <= 0.40) measure oscillation damping, not wall decay. **Failed setup**; kept in results/runs and summary.json
   for the record and replaced by the damped variant.
3. Radial spherical solver: 7 runs (3 physical with x >= 4, 4 PRS with x >= 2) overflowed after the first collapse
   (energy focusing at r = 0 with a fixed step); tau_10 precedes the failure (unaffected), debris diagnostics of those
   runs are invalid (listed in summary.json bubble_summary.*.numerical_failures_after_collapse).
4. Energy-budget check fails a strict 1 % instantaneous criterion for 2D PRS at the production step (§4).
5. Scaling is approached, not established: A still drifts by ~5-10 % per factor 2 in tau in 3D PRS, mu = -0.90 to
   -0.975; the damped late-start network has not converged (A ~ 1.1 at h = 1/4).
6. Finite-box 2D: spanning walls lost before h = 1/4 in 2 of 14 2D runs.
7. tau_10/tau_0 are quantised by the log-spaced output grid (~5-15 % spacing); 3D biased runs are box-limited for
   p = 0.45.
8. Literature: the task's identifier arXiv:1104.5925 returned HTTP 404 on arxiv.org/abs; the Leite-Martins paper
   "Scaling properties of domain wall networks" (PRD 84, 103523) is arXiv:1110.3486 and was used. The
   Garagounis-Hindmarsh amplitude is not reconciled (§5.2).
9. Not done: the optional charged-X run (X != 0, subcritical/supercritical charge relative to F^2 = 1); any thermal
   or charged formation calculation; wall gravity/backreaction; quantum fluctuations; decay through X/D; GW or CMB
   observables; 4096^3-class 3D runs; extrapolation of K10(p) outside 0.05 <= 1/2 - p <= 0.3.

## 9. Literature (PDFs kept in the scratch directory only)

| Id | Title | Read level | Used for |
|---|---|---|---|
| arXiv:1309.5001 (Hiramatsu, Kawasaki, Saikawa) | On the estimation of gravitational wave spectrum from cosmic domain walls | full text of §§II-IV (definition of A, simulation setup, Table I, A = 0.8 +- 0.1 eq. 29) | A definition, resolution criterion, comparison |
| arXiv:1110.3486 (Leite, Martins), PRD 84 103523 | Scaling properties of domain wall networks | §§II-IV and Tables I-III | PRS alpha/beta, xi_c/eta, gamma v, mu |
| hep-ph/0212359 (Garagounis, Hindmarsh), PRD 68 103506 | Scaling in numerical simulations of domain walls | §§II-IV and appendix (estimator, isotropic factors, Tables I-II) | area estimator; amplitude (not reconciled) |
| hep-ph/0410356 (Oliveira, Martins, Avelino) | The cosmological evolution of domain wall networks | abstract, eqs. 6, 10-12 and Table I | slow approach to scaling, mu |
| arXiv:1407.3905 (Correia, Leite, Martins), PRD 90 023521 | Effects of biases in domain wall network evolution | abstract and §IV | biased initial conditions decay; Hindmarsh law fails for biased IC |
| Press, Ryden, Spergel, ApJ 347 (1989) 590 | Dynamical evolution of domain walls in an expanding universe | not read directly (prescription as stated in the four papers above) | alpha/beta prescription, estimator |
| Coulson-Lalak-Ovrut; Larsson-Sarkar-White; Hindmarsh PRL 77 (1996) 4495 | biased walls | not read (cited via Correia et al.) | context only |

## 10. Files

- `code/wallsim.py` lattice evolution and diagnostics; `code/run_networks.py` job list and 2-worker runner;
  `code/bubbles.py` NG ODE, radial and Cartesian bubble solvers; `code/validate.py`; `code/scales.py`;
  `code/analyze.py`; `code/run_all.sh`.
- `results/headline.json` (compact machine-readable headline), `results/runs/*.json` (128 runs: config, dx, dt,
  time series), `results/summary.json` (aggregates, fits, resolution checks, bubble summary), `results/bubbles.json`, `results/validation.json`, `results/scales.json`.
  Total < 5 MB; no field dumps.
