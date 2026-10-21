# Verification v_charge: module charge_phase (+ relic_yield BEC claim)

21 October 2026. Adversarial verifier. My own code is in this directory, written from the action in BRIEF.md. No
module code was imported or copied. I read `z2_memory.py` only after my toy results existed, and only to compare
seed conventions; the module also uses `-eps*w(r)`. Codex data were read with `git show 135d535:...` and
`git show 479f111:...` (RESULTS.json, README.md, MATHEMATICS.md, REPORT.md, protocol.json). Numerical jobs used at
most 2 threads and each ran in under 5 minutes.

**Files.**

| File | Contents |
|---|---|
| `t0_check.py` → `results_t0.json` | T = 0 phase structure. Brute-force 2D Nelder-Mead, independent reduced analysis, and a real-field rotating-frame linear dispersion matrix |
| `tension.py` → `results_tension.json` | Q-matter / vacuum tension at mu = m*. (A) `solve_bvp` on the static field equations; (B) L-BFGS minimisation of the discretised 1D grand potential |
| `cosmo_check.py` → `results_cosmo.json` | Y_B from scratch (G_N from CODATA 2018), T_n(r), H/m_X, KZ length, n_Q,max, xi_BEC |
| `thermal_check.py` → `results_thermal.json` | The NOTES §3 leading-order (T, n) scheme implemented from its verbal definition. Analytic dF/dphi^2 and a full-F finite difference, plus a frozen-partition sensitivity variant |
| `memory_toy.py` → `results_memory.json` | My own spherical sign-memory toy |
| `archived_quench.py` → `results_archived_quench.json` | Archived soliton with its charge removed |
| `results_summary.json` | Combined numbers |

Units: v = 1. Densities are in v^3. ell = 1/(sqrt(lam_phi) v) = 1/(100 GeV).

## Claims

| # | Claim (module) | Independent method | My number | Verdict |
|---|---|---|---|---|
| 1 | \|X\|_c^2 = 0.5 v^2, n_c = 0.0316228 v^3 | closed form + 2D brute force (max rel. dev 2e-16) | 0.5; 0.031622777 | CONFIRMED |
| 2 | n_s1 = 0.042767 v^3 (x_s1 = 0.740741) | sign change of the branch-1 radial Hessian 2W'/y + W'' | 0.04276669; 0.7407407 | CONFIRMED |
| 3 | n_t = 0.040103 v^3; S jumps 0.6480 → 0, x 0.5801 → 1.1716 | root of rho_1(n) = rho_2(n) | 0.04010262; 0.64800, 0.58009 → 1.17161 | CONFIRMED |
| 4 | n_* = 0.177828 v^3, rho/n = m* = 56.234 GeV, x_* = 3.1623 | numerical min of rho/n on branch 2 and closed form | 0.17782794; 56.23413 GeV; 3.162278 | CONFIRMED |
| 5 | c_s^2 < 0 on the phi != 0 branch for all n < n_s1 | finite-difference n dmu/dn / mu, and the small-k eigenvalue of the 6×6 linear system | -4.4e-4 (n = 1e-4) … -0.170 (Codex state) … -4.7 (n = 0.0425); the two methods agree | CONFIRMED. This holds on the radially stable part. On the radially unstable segment (n from n_c up to n_s1) c_s^2 > 0, but that segment is a saddle. |
| 6 | c_s^2 = 1/3 on phi = 0; branch 2 linearly stable, S mass^2 = a x - 1 | the same two methods | 0.33333 at all n tested; max growth rate 0 | CONFIRMED |
| 7 | Growth rates 0.0045 / 0.045 / 0.112 / 0.226 at x = 0.01 / 0.1 / 0.25 / 0.5 | my own dispersion matrix | 0.004498 / 0.04485 / 0.11214 / 0.22640 | CONFIRMED |
| 8 | Codex production: n = 0.022009 v^3 = 0.696 n_c, c_s^2 = -0.1698, k < 0.595, gamma = 0.1121 at k = 0.404, 6 harmonics, ln(200)/gamma = 47 | same | 0.0220085; 0.69597; -0.16981; 0.5953; 0.11214 at 0.4035; 6; 47.2 | CONFIRMED |
| 9 | sigma_Q = 0.6800 sqrt(lam_phi) v^3 = 6.80e7 GeV^3, sigma^(1/3) = 408 GeV, virial 0.34/0.34 | BVP (L = 150 and L = 250 agree to 1e-10) and L-BFGS (h = 0.125) | 0.680013 (BVP), 0.680013 (L-BFGS); grad/pot = 0.340007/0.340007; 6.800e7 GeV^3; 408.2 GeV | CONFIRMED. The module's 0.679976 is 5e-5 low (discretisation). |
| 10 | Thin wall: E/Q = 0.8665 m_X at Q = 39,944 (R = 3.77); binding at Q ≈ 1.3e4; S_4 = 2.5e5 at n_c | hand recomputation with my sigma_Q | 0.8666; 1.34e4; 2.5e5 | CONFIRMED (arithmetic). Thin wall is not valid near n_c, as the module states. |
| 11 | Y_B = 8.674e-11 | omega_b rho_c/h^2/m_p/s_0, with rho_c from CODATA 2018 G | 8.6731e-11 (m_p); 8.6756e-11 (m_b = 0.938, Codex: 8.675814e-11); 8.7362e-11 (m_u) | CONFIRMED. The module's rho_c/h^2 = 1.05375e-5 is 7e-5 above the CODATA-2018 value 1.053672e-5. The baryon-mass convention moves Y_B by up to 0.7%. |
| 12 | T_n(r = 1) = 198.2 TeV; scaling r^(-1/3) v; n_t/n_s1/n_* rows (214.5/219.2/352.5) | n = r Y_B (2pi^2/45) g*s T^3 with g*s = 106.75 | 198.20 / 214.54 / 219.19 / 352.46 TeV; every r row matches | CONFIRMED. Assumes g*s = 106.75 (no thermalised nu_R, no hidden dof) and entropy conservation since T_n. |
| 13 | H/m_X = 5.5e-10 at T_n; KZ length 850 ell (r = 1) and 300 ell (r = 0.01); R_H = 1.8e9 ell; about 1e19 domains | H = sqrt(pi^2 g*/90) T^2/M_P | 5.52e-10; 845, 304; 1.81e9; 9.8e18 | CONFIRMED as arithmetic. KZ is a mean-field estimate. |
| 14 | "fragments in ~50–100 ell … which is 1e-9 of a Hubble time" (NOTES §5) | 100 ell / R_H | **5.5e-8** (at T_n); Codex: H dt(240) = 1.04e-7 | CORRECTED. It is (3–6)e-8 Hubble times, not 1e-9. Qualitatively harmless. |
| 15 | T_c0 = 1.5575 v (+0.64% over 1.54765); continuous, with phi_min = 0.057 at T_c0 - 0.0025 | root of dF/dphi^2(n = 0) using exact Bose integrals | 1.557496; 0.0573 | CONFIRMED (within the stated leading-order scheme) |
| 16 | n_c(T_h) rises to about 0.042 v^3 and falls to 0 at T_c0 (table: 0.03670, 0.03994, 0.04141, 0.04208, 0.04218, 0.04077; 0.0278 at 1.555; 0.0131 at 1.557) | my implementation of the stated scheme at 8 temperatures; analytic and full-F finite difference agree to ≤ 2e-5 | 0.036703, 0.039943, 0.041411, 0.042076, 0.042180, 0.040772, 0.027776, 0.013063 | CONFIRMED for the scheme. Sensitivity: freezing the gas/condensate partition gives 0.03228 (T = 0.5), 0.03574 (T = 1.0), 0.04155 (T = 1.5). The low-T rise is ~10% scheme-dependent; the ≈ 0.042 plateau is robust. |
| 17 | Codex: negative cells have F^2 > 1.07 / 1.28 / 1.78; quartiles 3.15/3.65/4.00, 3.04/3.47/3.93, 2.99/3.44/3.99; medians 1.09–1.15 x_*; volumes 1228/1301/1291 ell^3; every-step min of mean F^2 2.12/1.91/1.89 | read RESULTS.json at 135d535 | identical; volumes 1228×1, 549×(4/3)^3 = 1301, 1291×1 | CONFIRMED (quoted data) |
| 18 | Negative phi are oscillations about phi ≈ 0 in charged fragments, not opposite-vacuum domains | Analytic: at fixed F^2 > 1/a = 1, V(S) has a unique minimum at S = 0 with mass^2 = F^2 - 1 (median 2.4–2.65, so oscillation frequency ≈ 1.6 >> growth rate 0.11). Codex data: no cell below -0.5 at any of 4801 steps on 3 grids; zero negative low-X plateaus; min S -0.30 / -0.24 / -0.26 versus vacuum -1 | — | CONFIRMED for "not opposite-vacuum domains". "Oscillations" is an inference from the local potential and transient minima; per-cell time series were not archived. The every-step record gives only the min of the *mean* F^2 on negative cells, so per-cell F^2 < 1 at intermediate times is not excluded by the summaries. |
| 19 | Codex reports 248–315 components at the 2% threshold, so negative patches are smaller than a merged 6.7-ell ball | REPORT.md at 479f111 | 315 / 281 / 248 are components at 2% of peak **F^2** (charge network) | CORRECTED (attribution). These counts are not negative-phi patches. The size inference is indirect and Codex warns the count is not a soliton count. The merged-ball radius 6.64–6.77 is correct. |
| 20 | Sign-memory radius: R_first_flip = 8 / 14 / 20 ell (x_after = 0; eps = 1e-2 / 1e-4 / 1e-6), 12 / 20 / > 28 (x_after = 0.5); R_* ≈ 8 ell = 0.016 fm; eps_q ≈ 0.011 | my own spherical toy (same toy class), R step 1 | 8 / 14 / 20 exactly; x_after = 0.5: 10 / 18 / 25; eps_q = 0.01125; 8 ell = 0.0158 fm | CONFIRMED (order of magnitude; the x_after = 0.5 values differ by ≤ 20% because of grid and convergence) |
| 21 | Archived soliton with charge removed inherits +v for coherent seeds up to 0.3 | own leapfrog from the archived CSV | +v for eps ≤ 0.3; flips at eps = 0.5 | CONFIRMED |
| 22 | Archived soliton: F^2(0) = 3.917, S(0) = +0.049, F^2 > 1 to r = 3.38, S < 0.5 to 3.26 | read CSV | 3.9166, 0.04908, 3.376, 3.255 | CONFIRMED |
| 23 | relic_yield: n_Q,max = n(mu = M) = 0.011386 T^3, with M = sqrt(PiX) T = 0.034157 T; agrees with M T^2/3 to 3e-5 | exact Bose integral for particle + antiparticle (1 dof each); normalisation checked: massless small-mu limit gives mu T^2/3 to 5e-8 | 0.0113858; 2.96e-5 | CONFIRMED |
| 24 | relic_yield: xi_BEC = 7.09e-3 r^(1/3) | (r Y_B s_SM/T^3 / n_max)^(1/3) | 7.092e-3 | CONFIRMED. It uses visible-only entropy and assumes a kinetically equilibrated hidden bath, which is needed for PiX. |
| 25 | Z2 proof p(+) = p(-) | checked the logic | — | CONFIRMED (standard). It needs a P-invariant initial state, which the module states. |
| 26 | Literature: Felder et al. hep-ph/0012142, Laine–Shaposhnikov hep-ph/9804237, Kasuya–Kawasaki hep-ph/9909509 | fetched arXiv abstracts (abstract read) | — | CONFIRMED as used. Felder: "completes within a single oscillation". L–S: an infinite system at finite density is inhomogeneous with Q-balls at any T. K–K: size set by the "most developed" linear mode. |
| 27 | Critical endpoint between T_h = 1.545 and 1.55 v; n_global and binodal columns of the thermal table | — | — | NOT CHECKED. I computed only n_local. |
| 28 | Charged-interior control (x frozen at 3.4, kick -0.3/-0.6, S oscillates about 0); O(4) S_4 near n_*; xi-trajectory table | — | — | NOT CHECKED numerically. The control is consistent with the unique S = 0 minimum at x = 3.4. |

## Skeptical assessment

**Normalizations.** These are all right:

- The complex-field factors check out: n = 2 w |X|^2, |dX|^2 = F'^2/2, and the X thermal mass PiX = gX^2/12 + lam_X/3.
- The phi thermal mass is 0.004175. The X loop gives gX^2 T^2/6 through the Bose gas.
- The Bose-gas density for 2 dof is mu T^2/3.
- The Planck mass is the reduced M_P in H.
- The tension uses the canonical F = sqrt2 |X| and the grand potential at mu = m*. Both bulk phases have omega = 0, and V >= 0 on a grid.

**Sign-memory radius R_*.** I reproduce the toy numbers, but their physical relevance is limited.

- The toy assumes sudden, complete removal of charge from a patch. It also assumes a frozen X background, exact spherical symmetry and a coherent seed.
- With exact U(1)_X, Q-balls above Q ≈ 1.3e4 are classically bound (E/Q < m_X). In the cold H2 case their charge does not leave, so their phi ≈ 0 cores never "decide" a sign. R_* matters only for evaporating, fissioning or diluted lumps.
- For dilution, tau_Q ~ 1/H is about 1e9 ell. The module's own trend (memory grows with tau_Q) then pushes the relevant scale up towards xi_KZ ~ 1e2–1e3 ell. So "R_* ≈ 8 ell" is a sudden-quench lower bound, not a cosmological value.

**H1 sentence.** "A random-sign domain pattern spanning the horizon forms" is stated in bold. It rests on several unproven steps:

1. outside-action homogeneous super-horizon supply above n_s1;
2. the coherent phi oscillation decohering, which the module only says is "expected" and did not compute;
3. a mean-field KZ estimate;
4. the absence of phase separation in n_c < n < n_*, which uses a thin-wall estimate that fails near n_c.

Persistence and backreaction are not computed. It should read as a conditional formation statement ("would form if …").

**Thermal part.** The (T, n) scheme leaves the X gas and the condensate at different chemical potentials (mu_gas = M < mu_cond). That mismatch drives a partition term worth ~10% of n_c(T) at T ≤ 1 v. The module states the approximation but does not quantify it. I quantify it above.
