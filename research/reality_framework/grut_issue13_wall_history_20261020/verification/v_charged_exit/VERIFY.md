# VERIFY: modules/charged_exit_sim (verifier key v_charged_exit)

21 October 2026. Adversarial verification of `modules/charged_exit_sim` (NOTES.md, results/headline.json,
results/summary.json).

Rules of BRIEF.md were followed:
- no git;
- numpy and scipy only;
- every job single-threaded, with at most two jobs at once;
- every job under 25 min wall (the longest took 19 min).

All files are in this directory. Field states went to the scratch directory, not here.

The module's code was read only after my numbers existed. I read only two parts of it, the decharge procedure and
the L_vor routine, to locate differences. No module code was imported or copied.

## 1. Method (independent code)

**`vlat.py`** is my lattice code. It uses a different formulation from the module. The module rescales conformally
(chi = a f). Mine keeps the comoving fields f in {S = phi/v, q1, q2}, with X = v (q1 + i q2)/sqrt(2), and uses
their canonical momenta p_f = a^2 f'.

- Units: lengths in ell = 1/(sqrt(lambda_phi) v) = 1/m_X. With the benchmark couplings,
  A = g_X^2/lambda_phi = 1 and B = lambda_X/lambda_phi = 0.1.
- Equations: p_f' = a^2 lap f - a^4 dU/df and f' = p_f/a^2, with a = tau/tau_i (radiation era, conformal time).
- Integrator: Strang splitting in extended phase space. The kick is taken at fixed tau. The drift is exact:
  q += p tau_i^2 dt/(tau (tau+dt)).
- Grid: periodic 5-point (2D) or 7-point (3D) stencil. dt = 0.05 in 2D and 0.06 in 3D.
- Grids used: 2D at dx = 0.25 with 1024^2 (L = 256) and 2048^2 (L = 512); 3D at dx = 0.3 with 128^3 and 192^3.
- Noise: my own Gaussian random field. White noise is filtered by exp(-k^2 s_c^2/2) with s_c = 1 ell and has
  std 1e-2 in phi. X carries a relative complex perturbation of 1e-3 that co-rotates with the condensate.
- Optional frozen expansion (a fixed) gives the H -> 0 continuation.

**Other code.**
- `checks.py` writes `results/checks.json`: closed-form, convergence, Z2-mirror and diagnostic tests.
- `run_job.py` and `analyze.py` write `results/<run>.json` and `results/v_charged_exit_summary.json`.

**My wall criterion** was fixed before I read the module's code.
- Vacuum seeds are cells with S > 0.5 or S < -0.5. Seeds within 1 ell (physical) of Q-matter (F^2 = 2|X|^2/v^2 >= 1)
  are excluded.
- Every other cell is assigned to the nearest seed sign (periodic Euclidean distance transform).
- The wall is the interface between the + and - basins. Its length per area is (pi/4) x (sign-changing links) x dx /
  area. This Crofton estimate is unbiased for isotropic curves: a circle of R = 16 is reproduced exactly, and an
  axis-aligned straight wall reads pi/4 low.
- Two variants are also stored:
  - the same basin interface without Q exclusion (closest to the module's L_vor);
  - L_raw, the raw count of every S sign change (closest to the module's contaminated L_all / L_wall).
- Rates are d xi_c/d tau with xi_c = 1/L. The fit runs from tau_form + 5 to xi_c <= L_box/8 ("box-safe"), or over
  stated windows.

**Charge near walls**, measured two ways:
- (a) charge in "wall components": connected components of core cells (|S| <= 0.5 or F^2 >= 1) that touch both + and
  - seeds;
- (b) charge within 2, 4 or 6 ell (physical) of the basin interface, together with the fraction of area inside that
  distance.

## 2. Claims, numbers and verdicts

| # | Claim (module) | Module value | My value | Verdict |
|---|---|---|---|---|
| 1 | Conserved comoving charge Q = a^2 int (q2 q1' - q1 q2') d^3x; the code conserves it to round-off | max dQ/Q = 3.9e-16 | derivation in §3. My code: max \|dQ/Q\| <= 3.3e-16 in every production run (2D, 3D, frozen) and 3.8e-15 on the exact orbit | CONFIRMED |
| 2 | Homogeneous rotating phi = 0 condensate, \|X\|^2 = (n^2/(8 lambda_X))^{1/3} with n = n_i a^-3, is exact in radiation domination | V3: 8e-5 | derivation in §3. Lattice at dt = 0.01: the relation holds to 2.6e-6 (kappa = 2) and 4.2e-6 (kappa = 4); n a^3 is constant to 1.2e-13; S = 0 exactly | CONFIRMED |
| 3 | phi turns tachyonic at n = n_c, a = kappa^{1/3}; n_c = sqrt(8 lambda_X) \|X_c\|^3 = 0.0316 v^3 | a_c = 1.26 (kappa = 2); V4 growth 3.52e5 | m_phi^2(phi = 0)/(lambda_phi v^2) = kappa^{2/3}/a^2 - 1: +0.020 at 0.99 a_c and -0.020 at 1.01 a_c. n_c = 0.0316228 v^3 and n_s1 = 0.0427667 v^3. k = 0 lattice growth from a = 1 to 2 (tau_i = 40): 1.7082e12 against 1.7086e12 from an independent ODE (max deviation 2.3e-4 along the path) | CONFIRMED |
| 4 | Formation, kappa = 2: a_form = 1.55 (tau_i = 40) and 1.383 (tau_i = 150); n/n_c = 0.537 / 0.756 | as stated | a_form = 1.55 in all 7 tau_i = 40 runs (6 at 1024^2, 1 at 2048^2) and 1.383 in both tau_i = 150 runs (2048^2); n/n_c = 0.537 / 0.756. Identical to the output spacing | CONFIRMED |
| 5 | p(+) at formation = 0.5003 +- 0.0046 (11 2D runs); no charge-induced sign bias | 0.5003 +- 0.0046 | 9 charged 2D runs (6 at 1024^2, 3 at 2048^2): 0.4919, 0.4910, 0.4997, 0.5173, 0.5144, 0.4975, 0.5083, 0.5265, 0.4968 -> 0.505 +- 0.012 (std), standard error 0.004, so mean - 0.5 is 1.2 standard errors. Neutral: 0.502, 0.496. The exact Z2 mirror reproduces the mirrored evolution bit for bit (max \|S + S_mirror\| = 0) | CONFIRMED; my scatter is larger because my boxes are smaller |
| 6 | Walls present at every output of every charged run to the end; the network coarsens and is neither pinned nor decaying | true in all runs | true in all 21 of my 2D runs and both 3D runs. After formation the minority vacuum fraction never fell below 0.28 in 2D. At the end both signs are macroscopic. xi_c grows overall; isolated dips of at most 2 comoving units occur, mostly in the transient right after formation | CONFIRMED over the simulated range (§4 discusses what "persists" means) |
| 7 | 15-37% of the charge sits on walls (bare cores + attached Q-balls); 60-84% of the lump charge is attached | t40 15-17%, t150 25-27%, 3D 28%, frozen 28 -> 37% | wall-component charge at the end: t40 14.8 +- 1.0% (6 seeds, 1024^2) and 16.1% (2048^2); t150 27.3 / 28.0% at a = 1.6 and 24.6% at a = 2.5; 3D 192^3 30.2%; frozen 27-28% -> 32-33% over 130 tau (the module's 37% is at 260 tau). Charge within 2 ell of the interface: 21-25% of Q on 2.3-3.5% of the area, an 8-11x over-density. Attached share of the Q-matter charge: 66% (t150) to 86% (t40) | CONFIRMED |
| 8 | Decharged control at tau_i = 150 from a = 1.6: d xi_c/d tau = 0.42 vs 0.82, ratio 0.51 (0.58 box-safe) | 0.416 / 0.820 / 0.51 / 0.58 | own seed, 2048^2: 0.389 / 0.702, ratio 0.554 over the whole continuation. Box-safe (xi_c <= 64): 0.420 / 0.734, ratio 0.572. Common window tau = 240-282.5: ratio 0.537 | CONFIRMED |
| 9 | Decharged control at tau_i = 40 from a = 1.8: 0.53 vs 0.69, ratio 0.76 | 0.525 / 0.688 / 0.763 | own seed, 2048^2: 0.543 / 0.714, ratio 0.761 (first 40 tau: 0.770) | CONFIRMED |
| 10 | The charge slows coarsening by 1.3-2x (rate ratio 0.5-0.76) | 0.51-0.76 | 0.54-0.76 in the expanding pairs; 0.53 / 0.53 in two frozen pairs (box-safe) | CONFIRMED |
| 11 | Frozen (H -> 0) continuation: L_wall plateaus because of sign oscillations inside wall-attached Q-matter (artifact) while L_vor keeps falling; rate 0.32 vs 0.62, ratio 0.52 | L_wall about 0.027-0.028 flat; L_vor falls 3.7x over 260 tau | 2 seeds (301, 303), 2048^2, frozen at a = 1.6, tau 240 -> 370. L_raw stays at 0.0276-0.0319 from tau = 260 (L_raw-based rate 0.02 / 0.01), while L_basin falls 2.15x / 2.23x. Whole-run rates: 0.310 / 0.644 (ratio 0.48) and 0.306 / 0.720 (ratio 0.43). Box-safe ratios: 0.53 and 0.53 (common windows 0.53 and 0.58) | CONFIRMED (artifact diagnosis and ratio about 0.5) |
| 12 | "The slowdown does not grow within a run"; in the frozen run "the rate is the same in the first and last thirds" and "xi_vor grows linearly at a constant rate" while the charge per wall length rises | constant 0.32 over 416 ell | t40 pair: ratio 0.76 -> 0.74 (first / second half): constant. **Charge-rich cases are not constant.** Frozen, both seeds: the charged rate falls from 0.347 / 0.367 (first 65 tau) to 0.243 / 0.255 (last 65 tau) while the charge per wall length rises 2.5x / 2.7x. The decharged rate rises from 0.646 / 0.644 to 0.705 / 0.854, so the ratio goes 0.54 -> 0.34 and 0.57 -> 0.30. Expanding tau_i = 150: charged 0.417 -> 0.341, decharged 0.722 -> 0.752. The decline starts at tau of about 300-310 (xi_c about 55-59, just inside L/8 = 64) and continues into the box-affected range | REFUTED as a general statement. It holds for the t40 pair and inside the early box-safe windows. In 3 of my 3 charge-rich realizations the charge-induced slowdown grows within the run, while the module's single frozen realization shows a constant rate. Box effects cannot be excluded for the late halves, so "does not grow" is unsupported, not disproven |
| 13 | Sub-critical control C1 (n = 0.9 n_c on the +v branch, S0 = 0.811, F^2 = 0.342) keeps the sign | min f+ = 0.988, f- = 0 | branch state recomputed: F^2 = 0.34210, S0 = 0.81111, omega = 0.83193. 1024^2, tau_i = 150, a = 1 -> 2.2: f- = 0 at every output; min f+ = 0.992; mean S rises monotonically 0.811 -> 0.986; transient Q-matter at most 1.0% of the cells (a = 1.6), dissolving afterwards | CONFIRMED |
| 14 | Neutral control: a_form = 1.15; rate 0.73 (L_vor) / 0.75 (L_wall); the charged-from-scratch rate is lower | 0.727 / 0.753 +- 0.042; charged 0.536 / 0.565 | a_form = 1.15 (2 seeds). Box-safe rate (L = 256, short window): neutral 0.80 +- 0.09; charged 0.563 +- 0.062 (6 seeds) | CONFIRMED |
| 15 | Physical domain size at formation is 21-29 ell (charged) and 9-13 ell (neutral) | 21-29 | the module's own run files give 21.3 (3D) to 31.0 (A_k2_t150_s1) with L_wall and 26.5-32.8 with L_vor. Mine (basin): 29.3-29.5 (t40), 35.0 (t150); neutral 10.0-10.4 | CORRECTED: about 21-33 ell; the t150 runs exceed 29. Immaterial to the conclusions |
| 16 | 3D (one 256^3 run, L = 76.8): p(+) = 0.514; both vacua percolate at every output after formation | 0.514; 54% / 41% at the end | 192^3 (L = 57.6): p(+) = 0.469 at formation; both signs span all three axes at every output after formation; at the end f+/f- = 0.395/0.548 (p drifts to 0.419); 30% of Q on walls and 21% in Q-matter. 128^3 (L = 38.4): p(+) = 0.546, drifting to 0.744; at the last output (a = 2.2) the minority no longer spans all axes. L/xi_c at the end is 3.0 and 2.4 | CONFIRMED only as a small-box qualitative check. The percolation statement is box-size dependent at these sizes (§4) |
| 17 | The L_wall/L_vor artifact handling is sound; in expanding runs L_vor/L_wall = 0.91-1.00 | 0.91 (t150), 0.98-1.00 (t40) | L_basin/L_raw = 0.78 (t150 at a = 1.6) -> 0.90 (a = 2.5) and 0.98-0.99 (t40 end). Excluding Q-matter from the seeds changes the basin length by <= 0.6% in the t150 continuation and frozen runs, and by <= 6% just after formation, so the module's L_vor (no exclusion) is equivalent to mine. The frozen artifact is reproduced (claim 11) | CONFIRMED |
| 18 | The integrator converges as dt^2 (V5) | 8.8e-4 -> 2.2e-4 | my canonical KDK against an independent friction-form DOP853 solution (24^2, strongly nonlinear): max \|dS\| = 3.25e-3, 8.1e-4, 2.0e-4 and 5.1e-5 at dt = 0.02, 0.01, 0.005 and 0.0025 (ratio 4.0 per halving). Statistics at dt = 0.05 vs 0.025 (512^2): L_basin within 0.6%, f+ within 3e-4, wall charge within 0.2%; a_form and p_form identical | CONFIRMED for my code; the module's own value NOT CHECKED |
| 19 | Physical vs simulated H: H/m_X = 5.6e-10 at T_SM = 200 TeV; simulated 1.6e-2 to 2.6e-3; ratio 4.7e6 to 2.8e7 | as stated | H = sqrt(pi^2 g*/90) T^2/M_P = 5.62e-8 GeV. 1/(tau_i a_c^2) = 1.575e-2 (tau_i = 40), 4.20e-3 (tau_i = 150, kappa = 2) and 2.6e-3 (kappa = 4). Ratios 2.8e7 and 4.7e6 | CONFIRMED (arithmetic) |

## 3. Derivations

**Charge (claim 1).** Take L = -g^{mu nu} d_mu X* d_nu X - U on ds^2 = a^2(-dtau^2 + dx^2). The U(1) Noether current
j^mu = i g^{mu nu}(X* d_nu X - X d_nu X*) obeys d_mu(sqrt(-g) j^mu) = 0. With sqrt(-g) = a^4 and g^{tau tau} = -a^-2,
the conserved comoving charge is Q = a^2 int i(X* X' - X X*') d^3x.

With X = (q1 + i q2)/sqrt(2): i(X* X' - X X*') = i (i)(q1 q2' - q2 q1') = q2 q1' - q1 q2'. Hence
Q = a^2 int (q2 q1' - q1 q2') d^3x, the brief's expression.

In canonical variables Q = int (q2 p1 - q1 p2) d^3x, and the physical density is n = Q-density/a^3. Each sub-flow of
the splitting is U(1) invariant, so Q is conserved to round-off:
- The kick adds a^2 lap q - a^4 m^2(x) q. Its contribution to q2 p1 - q1 p2 sums to zero by summation by parts on a
  symmetric periodic stencil.
- The drift q += w p leaves q2 p1 - q1 p2 unchanged.

The module's Q is -1 times this: with chi = aY, 2(chi_r pi_i - chi_i pi_r) = a^2(q1 q2' - q2 q1'). This is a sign
convention only, as NOTES §2 states.

The module's Q0 values match n_i = kappa n_c: 165795/512^2 = 0.63246 = 2 n_c, and 74608.4/512^2 = 0.28461 =
0.9 n_c (in units of sqrt(lambda_phi) v^3).

**Exact condensate (claim 2).** At phi = 0 the X sector is a massless |X|^4 theory, which is classically conformal.
For conformal coupling the curvature term would be xi R |X|^2, and R = 6 a''/a^3 = 0 in radiation domination. The
brief also sets the xi couplings to zero at mu_0.

In conformal variables chi = aX the equation is chi'' - (a''/a) chi - lap chi + 2 lambda_X |chi|^2 chi = 0. With
a'' = 0 this is the flat-space equation, so chi = chi0 e^{i Omega tau} with Omega^2 = 2 lambda_X chi0^2 is exact for
all tau.

Physically:
- |X| = chi0/a and omega = Omega/a, so omega^2 = 2 lambda_X |X|^2.
- n = 2 omega |X|^2 = n_i a^-3, which gives |X|^2 = (n^2/(8 lambda_X))^{1/3}.
- Equivalently, in cosmic time with H = 1/(2t): (1/a)'' + 3H (1/a)' = (3H^2 - 3H^2)/a = 0, so the radial equation
  reduces to the rotating-orbit condition.
- phi = 0 is consistent because dU/dphi is proportional to phi (Z2).

The solution fails for a'' != 0 (a non-radiation era), with a thermal mass or bath, or with phi != 0. The module states
all three: classical, no bath, radiation background.

**Tachyonic onset (claim 3).** At phi = 0 the phi mass is m^2 = -lambda_phi v^2 + 2 g_X^2 |X|^2. It vanishes at
|X_c|^2 = lambda_phi v^2/(2 g_X^2) = v^2/2. On the orbit 2 g_X^2 |X|^2/(lambda_phi v^2) = kappa^{2/3}/a^2, so the sign
changes at a = kappa^{1/3}, exactly where n = kappa n_c a^-3 = n_c. In addition,
n_c = 2 omega_c |X_c|^2 = sqrt(8 lambda_X) |X_c|^3 = sqrt(0.008) x (1/2)^{3/2} v^3 = 0.0894427 x 0.353553 v^3 =
0.0316228 v^3.

## 4. Assessment of the qualitative claims

**"Persists".** This is supported only as "persists over the simulated range".
- In my charged runs the physical domain scale grows by 4.8 +- 0.5x at tau_i = 40 (xi_c by 2.5x) and by 6.1x at
  tau_i = 150 (from formation at a = 1.38 to a = 2.5). The module's "factor 2-8" is consistent.
- At the end the networks are sub-horizon (xi_c/tau = 0.23-0.40) and not in scaling.
- What the runs do show: no decay, no sign bias, and a coarsening speed of 0.2-0.6 c (comoving) that never drops to
  zero.
- Persistence to a wall-domination or late epoch is an extrapolation. It rests on the exact Z2, the absence of bias
  (claim 5) and 3D percolation, not on the simulated dynamic range. The physical range from formation to the
  horizon is about 1e6 in xi.
- NOTES §7, §8.3 and §9 E1 ("over the simulated range") say this. The §0 bullet "Persistence [C] ... It is not pinned
  and it does not decay" should carry "over a factor 2-8 in xi".
- "Not pinned" holds on the simulated time scale. "Does not grow" does not hold in my charge-rich realizations
  (claim 12). There, the friction from the Q-matter lenses increases as the charge per wall length grows. The
  module's own open item ("whether the Q-matter lenses eventually pin or detach at much larger charge per wall area")
  is therefore the live question, not a remote one. A slower (physical) exit puts more charge on walls.

**3D claim.** NOTES qualifies it in §8 item 2 ("one 256^3 run with L/xi of about 4-7 at formation; that is a small
box") and §8 item 4 (no L_vor in 3D). The §0 headline and §9 E1 ("both vacua percolate in 3D") state it without
that caveat.
- My two smaller 3D boxes show the box dependence directly:
  - at L = 57.6 both signs span at every output, which reproduces the claim;
  - at L = 38.4 one phase takes over (p goes 0.55 -> 0.74) and the minority stops spanning at the last output.
- The module's own 3D run also drifts (p goes 0.514 -> 0.568 at the end, with L/xi_c = 3.7).
- So the percolation statement is carried by the standard 3D result (p = 1/2 is far above the site threshold of about
  0.31), not by the single run.
- The run is fine as a qualitative 3D sanity check. Q-balls are real, and the charge on walls (module 28%, mine 30%)
  and in Q-matter (module 22%, mine 21%) agree.

**L_wall / L_vor artifact handling.** Sound.
- A basin (Voronoi) interface is the right observable once Q-matter sits on walls. My independent basin estimator,
  with explicit Q-matter exclusion, agrees with the no-exclusion version (the module's choice) to <= 0.6% late and
  <= 6% just after formation.
- The raw sign-change count overestimates the wall length by up to 53% just after formation (t150) and by 11% at the
  end of the expanding t150 continuation, in line with the module's L_vor/L_wall = 0.91.
- In both of my frozen runs L_raw is flat from tau = 260 on (rate 0.01-0.02), while the basin length falls 2.2x.
  Reading "freezing" from L_wall would be wrong, exactly as the module found.
- The pre-fix A, C2 and D3 rates are L_wall-based. In expanding runs they lie within about 10% of the basin rates (my
  L_raw/L_basin at the end of expanding runs is 1.00-1.11).

**Decharged control.** It removes the whole X sector, not only the charge.
- The phi vacuum value jumps from sqrt(1 - F^2) to 1, a small quench (F^2 of about 0.05-0.2 in the vacuum at
  a = 1.6-1.8). The measured factor is therefore "the effect of the charged X sector (condensate + Q-matter)".
  Attributing it to charge is reasonable but not strictly isolated.
- The quench transient does not drive the ratio: the first-40-tau ratios (0.537 at t150, 0.770 at t40, 0.53-0.58
  frozen) equal the box-safe ratios.
- The late-time ratios are below the early ones in the t150 and frozen runs (claim 12). The headline 0.5-0.76 range is
  an early-window value.

**Scope.**
- All runs are classical, 2D z-invariant (plus small 3D), with D = 0, no bath, and H/m_X 5e6 to 3e7 times the
  physical value.
- The trend with tau_i is monotonic: a slower exit puts more charge on walls and in Q-matter. The physical slowdown
  factor is therefore not bounded by these runs and could be larger. The module states this (§7, §8), and the
  conditional phrasing of §9 is appropriate.
- Not checked by me:
  - the kappa = 4 runs;
  - the C2 forced-restoration controls (p(+) = 0.47 / 0.50 / 0.71);
  - the module's resolution pair;
  - the 2D topology-class statistics;
  - the frozen continuation beyond tau = 370 (the module ran to 500);
  - lump statistics (Q = 2.9e5 quanta per ell).

## 5. Files

- Code: `vlat.py`, `checks.py`, `run_job.py`, `analyze.py`.
- `results/checks.json`: closed-form, convergence, Z2-mirror and synthetic-diagnostic checks.
- `results/<run>.json`: per-run configuration and time series.
  - 2D: A_k2_t40_N1024_s101-106, B_t40_N1024_s201-202, C1_k0.9_t150_N1024_s401,
    A_k2_t150_N2048_s301_to1.6, CC/DC_k2_t150_N2048_s301_from1.6, A_k2_t40_N2048_s302_to1.8,
    CC/DC_k2_t40_N2048_s302_from1.8, A_k2_t150_N2048_s303_to1.6, F_k2/F_dc_k2_t150_N2048_s301_freeze1.6,
    F_k2/F_dc_k2_t150_N2048_s303_freeze1.6, A_k2_t40_N512_s501_dt05/dt025.
  - 3D: D3_k2_t40_N128_s601, D3_k2_t40_N192_s602.
- `results/v_charged_exit_summary.json`: per-run metrics, groups, decharged and frozen comparisons (including
  half-run rates).
- Reproduce, one job per call: `OMP_NUM_THREADS=1 python3 run_job.py <name>`. The CC/DC/F jobs need the
  corresponding `*_to1.x` state, written to `$VCE_SCRATCH`. Then run `python3 analyze.py` and `python3 checks.py`.
