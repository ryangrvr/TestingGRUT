# Verification v_network: module network_sim (Z2 wall network, X = D = 0, radiation era)

21 October 2026. Adversarial verifier. I wrote all code here myself (`lattice_v.py`, `run_network.py`, `validate_v.py`,
`ng_v.py`, `bubble_v.py`, `cosmo_v.py`, `analyze_v.py`, `run_queue*.sh`). I did not import or copy the module code. I
read the module code (`wallsim.py`, `bubbles.py`) only after I had my own 2D numbers, to find where the results
differed. Libraries: numpy and scipy. At most 2 threads were used (two single-thread workers). Longest job: 7.3 min.
Total CPU about 32 min over about 30 min wall. This is above the "25 min total CPU" reading of the task, because I
added two independent-seed 2048^2 runs once I saw that seeds 1 and 2 copy the module's realizations (see V-a2).

Results: `results/summary_v.json` (aggregates), `results/runs/*.json` (time series), `results/bubble_v.json`,
`results/ng_v.json`, `results/validation_v.json`, `results/cosmo_v.json`, `results/claims_v.json`.

## Method (independent)

- Units: v = 1, lam = 2, so the kink is phi = tanh(z/w) with w = 1, m = 2, sigma = 4/3. The module uses lam = 1/2
  with length unit 1/m_phi = w/2, so tau_module = 2 tau_here. Same physical grid as the module: dx = w/2.
- Equation: phi'' + alpha phi'/tau - lap phi + U'(phi) = 0 (beta = 0), a proportional to tau.
- Integrator (not the PRS one): the friction is handled exactly through P = tau^alpha phi'. phi' = P/tau^alpha and
  P' = tau^alpha F(phi), stepped with a staggered leapfrog. dt = 0.3 dx in 2D and 0.25 dx in 3D. 5-point (2D) and
  7-point (3D) Laplacian, periodic box, float64.
- IC: white noise with a Gaussian Fourier filter exp(-k^2 s_c^2/2), s_c = 0.75 w, std 0.1 v, zero mean and zero
  velocity at tau_i = 1 w. In physical units this matches the module's IC (s_c = 1.5, tau_i = 2 in module units).
  **Consequence:** with the same `default_rng(seed)` call, my seeds 1 and 2 reproduce the module's seed-1 and seed-2
  realizations on the same grid. My time series match theirs to about 1 % (code agreement), but those two seeds
  are not new statistics. Seeds 11 and 12 are independent realizations.
- Wall measure, four estimators: (i) raw count of sign-changing links x dx^(d-1) (Manhattan); (ii) raw x pi/4 (2D)
  or x 2/3 (3D); (iii) gradient-weighted link estimator (HKS eq. 14 / Garagounis-Hindmarsh); (iv, 2D only)
  marching-squares contour length with linear interpolation. Validation (`validation_v.json`): on tanh circles
  R = 8-32 the ms estimator is within 2e-4 and grad within 1.3e-3; tilted lines are exact (ms) or within 1.1e-3
  (grad); on spheres grad is within 1.8e-3. For a circle the raw count is exactly 4/pi x 2 pi R. The single-
  orientation isotropic errors are 0.79-1.11. Flat-space energy error: 1.0 % (dt = 0.15) and 0.26 % (dt = 0.075),
  converging as dt^2.
- **A definition used:** rho_w = A sigma / t, with t cosmic time. The physical wall energy per physical volume (3D)
  or area (2D) is sigma a^(d-1) S_c / (a^d V_c) = sigma S_c/(a V_c). In RD, a = tau/tau_1 and
  t = int a dtau = tau^2/(2 tau_1) = a tau/2. So **A = (S_c/V_c) tau / 2**. That is: comoving length per comoving
  area (2D), or comoving area per comoving volume (3D), times conformal time, divided by 2. The tau_1 normalization
  drops out. The same formula is used in 2D PRS as a convention: a 2+1 PRS run with a ~ tau has no physical
  radiation-era t.

## Claims and verdicts

### Definitions and conversions

| # | Claim (module) | My method | My number | Verdict |
|---|---|---|---|---|
| D1 | A = rho_w t/sigma = (A_c/V_c) tau/2 = tau/(2 xi_c). This is the HKS definition (their eqs. 4, 12, 13). | Own derivation. Read HKS 1309.5001v2 §II eq. 4 (A = rho_wall t/sigma_wall), §III.B eqs. 12-13 (rho = sigma A/(a V), A = A t/(a V)), §III.D (a ~ tau, a(tau_i = 2) = 1). | t = a tau/2 holds for any a = tau/tau_1. A = (A_c/V_c) tau/2 (checked numerically in cosmo_v.json). | CONFIRMED |
| D2 | "Comoving area density x conformal time = 2A in RD" | same | 2A | CONFIRMED |
| D3 | a = tau gives H = 1/tau^2, (aH)^-1 = tau, t = tau^2/2, H^-1/delta = tau^2/2 (512 at tau = 32, 72 at tau = 12) | arithmetic; delta = 2/m_phi; module code `a_of(t) = t/tau_a` | 512, 72 | CONFIRMED |
| D4 | sigma = 9.428e7 GeV^3, delta = 0.014142 GeV^-1, H^-1/delta = 2.1e13 at T = 1.5476 v, g* = 106.75 | cosmo_v.py, reduced M_P | 9.4281e7; 0.0141421; 2.101e13 | CONFIRMED (T and g* are illustrative inputs) |
| D5 | t_w = 3 M_P^2/(4 A sigma) = 3.9e4 s (A = 0.8), 6.2e4 s (A = 0.5) | rho_rad = 3 M_P^2 H^2 with H = 1/(2t) and reduced M_P | 3.88e4 s, 6.21e4 s | CONFIRMED (conditional, as the module labels it) |
| D6 | D_h = alpha + beta/2. PRS 2D alpha = 2 gives D_h = 2; physical z-invariant gives D_h = 3 | re-derived: comoving kink width a^(-beta/2), energy per comoving area a^(alpha + beta/2) | as stated | CONFIRMED |
| D7 | Isotropic factors pi/4 (2D), 2/3 (3D) agree with the gradient estimator to <= 1.2 % | own estimators in the networks | 2D: A_iso/A_ms within 0.6 % in the window; 3D: A_iso/A_grad within 0.3 % | CONFIRMED |

### (a) 2D PRS alpha = 2, beta = 0, 2048^2, dx = w/2, L = 1024 w, window h = tau/L in [1/8, 1/4]

| # | Claim | My number | Verdict |
|---|---|---|---|
| a1 | Code reproduces the 2D PRS network | Seeds 1 and 2 (same realizations as the module's): A_ms = 0.6490, 0.6499 (grad 0.6508, 0.6516). Module seeds 1 and 2: 0.6525, 0.6506. The time series agree to about 1 %, including f_plus = 0.830 at tau = L/2 in seed 1. | CONFIRMED (independent integrator and estimator) |
| a2 | **A = 0.652 +- 0.002 (3 seeds)** | Independent seeds 11 and 12: A_ms = **0.748, 0.760** (grad 0.750, 0.762). Over all 5 distinct realizations (module 1, 2, 3 + mine 11, 12): **A = 0.693, seed std 0.056, s.e.m. 0.025.** Per-output values inside the window range from 0.58 to 0.80. | **CORRECTED: A(2D PRS alpha=2) = 0.69 +- 0.06 (realization std). The +-0.002 is not a valid uncertainty: the module's three seeds happen to cluster, and its window holds only 3 output times.** The module's single 4096^2/L=4096 value 0.725 fits this spread. |
| a3 | mu = -0.991 +- 0.010 (2D PRS alpha=2) | Fit of S_ms ~ tau^mu over [L/32, L/4]: seeds 1, 2: -1.002, -0.993; seeds 11, 12: -0.873, -0.907. Mean -0.944 +- 0.064. | CORRECTED: mu = -0.94 +- 0.06. It is realization-dependent; A is still rising through the window in seeds 11 and 12. |
| a4 | v_rms 2D PRS = 0.487 +- 0.011 | sqrt(sum phi'^2/sum grad^2) over abs(phi) < 1/2: 0.498, 0.477, 0.454, 0.454 | CONFIRMED within the scatter (about 0.47 +- 0.02) |
| a5 | Manhattan counting correction | Raw count / A_ms = 1.27 = 4/pi in the network | CONFIRMED |

### (b) 3D PRS alpha = 3, beta = 0, 160^3, dx = w/2, L = 80 w, 2 seeds

| # | Claim | My number | Verdict |
|---|---|---|---|
| b1 | A = 0.777 +- 0.013 (256^3, L = 128 w, window tau in [16, 32] w) | My box is smaller, so its window h in [1/8, 1/4] is tau in [10, 20] w: A_grad = 0.716, 0.737 (mean **0.726 +- 0.015**). At matched tau in [16, 20] w: **0.747 +- 0.034**. The module's 256^3 runs at the same tau (their time series): 0.748-0.763. Box-affected h in [1/4, 1/2]: 0.735, 0.875 (module A_box 0.849 +- 0.043). | CONFIRMED to the extent my box allows: A(tau) agrees with the module's within 1-3 % at matched tau, and A is still rising. I could not test the value 0.777 itself (that window needs tau up to 32 w with h <= 1/4, so L >= 128 w). |
| b2 | 3D mu = -0.896 +- 0.032 (256^3), approach to scaling, not a plateau | Fit over tau in [5, 20] w: -0.886, -0.822 | CONFIRMED qualitatively (slow approach, mu above -1). The fit windows differ. |
| b3 | v_rms 3D PRS = 0.456 +- 0.006 | 0.476, 0.457 | CONFIRMED (about 0.47) |
| b4 | Comparison with HKS A = 0.8 +- 0.1 | Read HKS §IV.A. Their eq. 29 gives A ~= 0.8 +- 0.1, where +-0.1 is the realization scatter of their 10 x 512^3 runs. They used the physical equation (10), 4th-order stencils, tau 2 -> 40/60 (h up to 1/2), and note a late rise from the box. | CONFIRMED: same definition (rho_wall t/sigma with cosmic t, comoving area via eq. 13). The module's 3D PRS 0.777 and my 0.73-0.75 (earlier tau) lie within 0.8 +- 0.1. HKS is a physical-equation result, not PRS; the module states this. |

### (c) Biased initial population, 2D PRS alpha = 2, 1024^2, L = 512 w

| # | Claim | My number | Verdict |
|---|---|---|---|
| c1 | At p = 0.35 the network disappears. tau_10 = 37.3 +- 4.0, tau_0 = 78 +- 22, K10 = 2.99 (module units) | Three seeds (1, 2 shared with the module; 11 independent). Walls entirely gone (no sign change) at tau = 35.0, 44.9, 44.9 w (= 70, 90, 90 module units). tau_10 (A < 0.1 A_ref) = 18.85 w (= 37.7 module) in all three. K10 = 2.9-3.1. Independent seed 11: A = 0.29 at tau = 4 -> 0.0023 at tau = 31 -> 0 by tau = 50.8, with f_plus -> 0. | CONFIRMED |
| c2 | p = 0.45: tau_10 = 155.9, tau_0 = 250, K10 = 10.8 (module units) | seed 1: tau_10 = 73.9 w (= 148 module), tau_0 = 137.5 w (= 275 module), K10 = 9.6-10.2 | CONFIRMED (same realization as the module's seed 1) |
| c3 | Reference p = 1/2 run keeps walls | 1024^2 p = 0.5 seed 1: walls present to tau = 256 = L/2, A = 0.86, f_plus = 0.67 at the end | CONFIRMED |

### (d) Isolated circular bubble, 2D, thin wall, from rest

| # | Claim | My method | My number | Verdict |
|---|---|---|---|---|
| d1 | x << 1: Delta tau_c = 1.5708 R0 (circle). Sphere 1.3110 R0. | Own NG ODE (Euler-Lagrange re-derived: R'' = -(1-R'^2)[kappa/R + D_h R'/tau]), DOP853 at rtol 1e-11. Closed forms pi/2 and Gamma(1/4)^2/(4 sqrt(2 pi)). | 1.570796, 1.311029. x = 1e-3: 1.5711 (circle D_h = 2). x = 0.01: 1.5735. | CONFIRMED |
| d2 | Field theory agrees with NG | My Cartesian lattice, R0 = 20-80 w, x = 0.01 and flat space, dx = w/2, w/4, w/8 | Delta tau_10/R0 (NG flat 1.4706, x = 0.01: 1.4733): R0 = 40, x = 0.01: 1.5706 / 1.4955 / 1.4786 at dx = 0.5 / 0.25 / 0.125 (dev +6.6 / +1.5 / +0.36 %). Flat: +6.7 / +1.6 / +0.38 %. R0 = 80, dx = 0.5: +8.2 %. Convergence goes as dx^2. | CONFIRMED in the continuum limit. **Caveat:** at the module's network resolution (dx = delta/2) a large ultra-relativistic collapse runs 5-8 % slow, from lattice dispersion of the Lorentz-contracted wall. The module's own Cartesian 3D check (R0 = 4, x = 1, `bubbles.json` cartesian_checks) gives tau_10 = 10.0 vs NG 9.415 (+6.2 %), and its paired radial run gives +4.0 %. That exceeds the headline "max 2.5 %", which holds only for the 22-run radial set. Its tau_c there is 11.6 vs NG 9.81. |
| d3 | x >> 1: tau_c -> 2.459 R0 (cylinder), 2.214 R0 (circle, D_h = 2), 1.867 R0 (sphere). t_c = 3.49 t_R. | NG ODE at x = 1000 | 2.4589, 2.2139, 1.8670; 1.8670^2 = 3.486 | CONFIRMED |
| d4 | tau_c/tau_i table (x = 0.178, 1, 5.6, 31.6) | NG ODE | sphere 1.238, 2.454, 10.62 (x = 5.6; 10.66 at x = 10^0.75 = 5.623), 59.03; cylinder 1.292, 2.917, 13.88, 77.72; circle 1.288, 2.787, 12.57, 69.99 | CONFIRMED (the module's x values are 10^0.75 and 10^1.5) |

### Persistence and decay statements, literature, wording

| # | Claim | Evidence | Verdict |
|---|---|---|---|
| P1 | Persistence at p = 1/2: a network forms and walls are present at every output to tau = L/2 | All 7 of my unbiased runs (2D: 4 x 2048^2, 1 x 1024^2; 3D: 2 x 160^3) had walls at every output to the end (L/2, or L/4 for seeds 11 and 12). The 2D boxes drift toward one phase (f_plus = 0.83, 0.67 at L/2). | CONFIRMED as a statement limited to the box and the simulated hierarchy. Persistence to H^-1/delta ~ 1e13 is an extrapolation (scaling assumption), which the module states. |
| P2 | 3D: a box-spanning wall carries >= 97.4 % of the area | not computed (no cluster labelling in my code) | NOT CHECKED |
| P3 | Biased networks decay (C2): "the network disappeared in every tested case after its length scale grew by K10 <= 11 (2D) or <= 2.7 (3D, box-limited)" | My runs: p = 0.35 and 0.45 decay completely. **The module's own summary.json contradicts the literal wording:** n_seeds_gone_by_end = 1/3 at p = 0.45 in 3D PRS and in 2D physical, and 2/3 at p = 0.40 in 2D physical. | PARTLY SUPPORTED, with a WORDING OVERCLAIM. Every tested case fell below 0.1 A_ref (tau_10 criterion), but complete disappearance was not seen within the box time in 2/3 seeds at p = 0.45 (3D PRS, 2D physical) and in 1/3 at p = 0.40 (2D physical). |
| P4 | "For a formation hierarchy H^-1/xi_f ~ 1e13 this occurs long before the network reaches horizon size" | K10 was transferred between formation hierarchies of about 7 and about 70 only (module §5.3). | PLAUSIBLE but extrapolated over about 11 decades. It should be labelled an extrapolation, not a result. |
| L1 | HKS A = 0.8 +- 0.1, same definition | HKS 1309.5001v2, read level: §II (eqs. 3-4), §III intro, §III.A-B (eqs. 11-16), §III.D (Table I), §IV.A (eq. 29, Fig. 1 caption). PDF text extraction in my scratch directory. | CONFIRMED |
| L2 | Leite-Martins 1110.3486: xi_c/eta = 0.60 +- 0.05 / 0.64 +- 0.03, so A ~ 0.83 / 0.78 "indicative" | Read level: §III-IV text and Table III (512^3 and 1024^3 rad. rows). LM write only A/V proportional to 1/xi_c. | CONFIRMED (the module correctly calls the conversion indicative) |
| L3 | Garagounis-Hindmarsh hep-ph/0212359: 0.928 t^-0.996 (2D rad.), 0.925 t^-0.994 (3D), giving A ~ 0.46, "not reconciled" | Read level: §II-IV and Tables I-II. The quoted numbers are present, and the tables show the isotropically corrected A. My extraction of their first-order eq. (7) shows friction (R'/R + eta) pi with no visible factor alpha. That could explain the factor ~1.7 (weaker damping), but I did not establish it. | Numbers CONFIRMED. Reason for the discrepancy NOT CHECKED. The module's "not reconciled" is honest. |
| L4 | Correia-Leite-Martins 1407.3905: biased initial conditions decay; Hindmarsh law fails for biased IC | Read level: abstract | CONFIRMED at abstract level |

## Not checked

2D and 3D physical-equation runs (A = 0.709, 0.884); the spanning-cluster fractions; the 3D biased runs; K10 fits and
exponents; the damped late-start preparation; the radial solver's debris lifetimes; the energy-budget failure in 2D
PRS (my integrator's flat-space energy check passes, but I did not reproduce theirs).

## Bottom line

- The definitions and unit conversions are correct: A vs rho_w, comoving vs physical, t = a tau/2, reduced M_P in
  t_w. No factor-2 or time-convention error found.
- 2D PRS A = 0.652: the central value is low and the error is understated by about 30x. The corrected value is
  0.69 +- 0.06 (5 realizations).
- 3D PRS: consistent with the module at matched tau, and with HKS 0.8 +- 0.1.
- Biased p = 0.35 decay: confirmed, including with an independent seed.
- NG bubble collapse numbers: confirmed exactly. The field theory converges to NG, but is 5-8 % slow at dx = delta/2.
- Wording overclaim: C2's "disappeared in every tested case" is contradicted by the module's own n_seeds_gone_by_end
  counts at p = 0.45 and 0.40.
- Extrapolation: K10 is extrapolated to a hierarchy of 1e13.
