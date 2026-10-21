# Module charged_exit_sim: lattice simulation of the charged exit (phi, X; D = 0) in a radiation era

21 October 2026. Claude module of the issue #13 wall-history packet (BRIEF.md). The action is unchanged. This
module adds no bias, portal or inflaton coupling. The charged preparation (a homogeneous X condensate, or charge put
into a +v vacuum) is **outside-action preparation data**. The reporter D is set to zero (allowed by the task; see §2).
Nothing here validates or excludes Candidate 1. The statements are conditional history constraints for a classical
field system whose expansion rate is far faster than the physical one (§7).

Labels: **[C]** computed by the code here, **[D]** derived in closed form, **[Q]** quoted from sibling modules
(charge_phase, network_sim, relic_yield) or from the brief, **[E]** estimate.

## 0. Result in one paragraph

A classical lattice simulation solved the unchanged phi-X equations in a radiation background. The X sector starts
as the exact homogeneous phi = 0 rotating condensate at n_i = 2 or 4 n_c. Expansion dilutes the charge through n_c,
and the field then breaks into domains of both signs.

- **Formation [C].** In 11 2D runs (2048^2) the formation fraction is p(+) = 0.5003 +- 0.0046 (mean +- seed std).
  The single 3D run (256^3) gives 0.514. Both vacua percolate in 3D at every output after formation. In 2D, p = 1/2 sits exactly at
  the percolation threshold, so which sign wraps the torus is a lottery, as in the neutral runs.
- **Persistence [C].** Walls were present at every output of every charged-exit run up to the end. The network
  coarsens. It is not pinned and it does not decay.
- **Effect of the charge [C].** 15-37% of the charge sits on walls, as bare wall cores and attached Q-balls. The charge
  slows the coarsening by a factor of about 1.3 to 2, and the slowdown does not grow within a run.
  - The cleanest test starts from the same phi configuration at a = 1.6 and compares runs with and without the X field
    ("decharged").
  - Expanding background, tau_i = 150: d xi/d tau = 0.42 with charge, 0.82 decharged (ratio 0.51; 0.58 in the
    box-safe window).
  - Expansion switched off at a = 1.6 (Minkowski, the H -> 0 limit): 0.32 with charge, 0.62 decharged (ratio 0.52).
  - Expanding background, tau_i = 40, where only about 16% of the charge is on walls: 0.53 with charge, 0.69 decharged
    (ratio 0.76).
  - In the Minkowski run the charge on walls grows from 28% to 37% and Q-lumps reach Q = 2.9e5 quanta per ell. The rate
    stays constant over 416 ell of physical time.
- **Sign bias [C].** No charge-induced sign bias is seen: p(+) at formation has mean 0.5003 over 11 charged 2D runs. This
  is expected from the exact Z2 symmetry (§4 of charge_phase), and an exact mirror test reproduces the mirrored
  evolution bit for bit.

The sign-memory controls go both ways:

- **Charge below n_c on a +v background** (C1, 0.9 n_c): not a single -v cell or sign change appears. Memory is kept.
- **Charge forced above n_s1 into a +v vacuum** (C2): phi oscillates coherently through zero and decoheres within about
  40 ell. With 4 n_c the sign memory is lost (p(+) = 0.47 and 0.50 at formation, 2 seeds) and a persistent
  network results. With 2 n_c the oscillation leaves a bias (p(+) = 0.71) and the minority phase disappears.

Conditional answer to Codex: within these classical simulations, a charged preparation that exits through the
n_c spinodal produces the same kind of macroscopic, persistent, unbiased wall network as the neutral quench, only
slowed by a factor of about 1.3 to 2. Its relevance to a charged preparation therefore cannot be dismissed. This holds within
§8's limits (classical, D = 0, no bath, mostly 2D z-invariant, simulated H/m about 1e7 times the physical value, finite
box and dynamic range).

## 1. Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/charged_exit_sim/code
OMP_NUM_THREADS=1 python3 validate.py            # 12 s, writes results/validation.json
# 1) runs that need a saved state, in this order (the S_* run writes the state), one worker each:
OMP_NUM_THREADS=1 python3 runs.py worker 1 2 S_k2_t150_s1_to1.6,F_k2_t150_freeze1.6,F_dc_k2_t150_freeze1.6,CC_k2_t150_s1_from1.6,DC_k2_t150_s1_from1.6 &
OMP_NUM_THREADS=1 python3 runs.py worker 0 2 S_k2_t40_s1_to1.8,CC_k2_t40_s1_from1.8,DC_k2_t40_s1_from1.8 &
# 2) then everything else: two single-threaded workers share the list (finished runs are skipped);
#    each job checkpoints after 22 min wall and is resumed later by the same worker
OMP_NUM_THREADS=1 python3 runs.py worker 0 2 &
OMP_NUM_THREADS=1 python3 runs.py worker 1 2 &
python3 analyze.py                                # results/summary.json, results/headline.json
```

Checkpoints, saved states and snapshots go to `$CES_CKPT`, which defaults to this session's scratch directory.
The total cost was about 6.3 CPU-hours on 2 cores, of which 1.2 h went to runs superseded by the L_vor re-runs. A 2048^2 run takes 0.14 s per step (6-17 min per run), and every
job stayed under 22 min wall. Only numpy and scipy are used, with deterministic seeds. A determinism check
(S_k2_t150_s1_to1.6 against A_k2_t150_s1) gives a difference of exactly 0. No field dumps are in results/.

## 2. Equations and conventions [D]

**Action and units.** This is the unchanged Candidate 1 scalar sector with D = 0. D enters its own equation linearly, so
D = 0 is an exact classical solution; D fluctuations are not included. Fields are in units of v. Lengths and times are
in units of ell = 1/(sqrt(lambda_phi) v) = 1/(g_X v) = 1/(100 GeV) = 1/m_X; comoving lengths equal physical ones at
a = 1. With S = phi/v and Y = X/v:

    Ubar = (S^2 - 1)^2/4 + A_G S^2 |Y|^2 + B_X |Y|^4,   A_G = g_X^2/lambda_phi = 1,  B_X = lambda_X/lambda_phi = 0.1

The energy density unit is lambda_phi v^4. The kinetic terms are -1/2 (dS)^2 - |dY|^2.

The charge density unit is sqrt(lambda_phi) v^3 = 0.1 v^3, so the charge_phase thresholds are:

- nbar_c = 0.316228 (n_c = 0.0316 v^3);
- nbar_s1 = 0.427667 (n_s1 = 0.0428 v^3);
- nbar_* = 1.77828 (Q-matter density n_* = 0.178 v^3).

Vacuum masses: m_phi = sqrt(2), i.e. 141 GeV, and m_X = 1, i.e. 100 GeV. The planar wall is S = tanh(z/delta) with
delta = sqrt(2) ell. Its tension is 2 sqrt(2)/3 = 0.9428 in units of sqrt(lambda_phi) v^3, i.e. 9.43e7 GeV^3 at v = 1 TeV.

**Background.** ds^2 = a^2(tau)(-dtau^2 + dx^2), radiation era a = tau/tau_i (a = 1 at the start; a'' = 0). The physical
Hubble rate is H = a'/a^2 = 1/(tau_i a^2) in units of m_X. Cosmic time is t = a tau/2.

**Field equations.** The requested comoving equations f'' + 2(a'/a) f' - lap f + a^2 dU/df* = 0 (for f = S and Y) are
solved exactly in conformal variables chi_S = a S, chi_Y = a Y. Because a'' = 0 in radiation domination:

    chi_S'' = lap chi_S - chi_S (chi_S^2 - a^2 + 2 A_G |chi_Y|^2)
    chi_Y'' = lap chi_Y - chi_Y (A_G chi_S^2 + 2 B_X |chi_Y|^2)

This is flat-space dynamics with Hamiltonian

    H(tau) = sum dx^d [ p_S^2/2 + |grad chi_S|^2/2 + |p_Y|^2 + |grad chi_Y|^2 + Ut ],
    Ut = (chi_S^2 - a^2)^2/4 + A_G chi_S^2 |chi_Y|^2 + B_X |chi_Y|^4.

Only the vev term depends on time. Energy budget: dH/dtau = -(a a') sum dx^d (chi_S^2 - a^2).

The two forms are equivalent because f'' + 2(a'/a) f' = (chi'' - (a''/a) chi)/a and a^2 dU/df = (1/a) dUt/dchi.
Validation V5 checks this numerically against an independent RK4 integrator of the friction form. The maximum field
difference is 8.8e-4 at dt = 0.02 and 2.2e-4 at dt = 0.01, i.e. dt^2 convergence.

**Charge and its a-power [D].** With j^mu = i g^{mu nu}(X* d_nu X - X d_nu X*), d_mu(sqrt(-g) j^mu) = 0 and
sqrt(-g) g^{tau tau} = -a^2, the conserved comoving charge is

    Q = a^2 int i (X* X' - X X*') d^3x = a^2 int (q2 q1' - q1 q2') d^3x,   X = (q1 + i q2)/sqrt(2).

The physical density is n = Q/(a^3 V_comoving). In conformal variables Q = 2 int (chi_r pi_i - chi_i pi_r) d^dx. The
code's Q is minus the brief's form, which is only a sign convention; the code makes the initial charge positive. The 2D
fields are z-invariant, so Q is per unit comoving z-length and n dilutes as a^-3.

**Integrator.** Kick-drift-kick (velocity Verlet) for the time-dependent Hamiltonian, i.e. Strang splitting. The
Laplacian is the periodic 5-point (2D) or 7-point (3D) stencil, in float64, with dt = dx/5.

- Each sub-flow is U(1)-invariant: the kick adds (lap - m^2(x)) chi_Y, whose contribution to sum chi* p is real, and the
  drift adds p to chi. Q is therefore conserved to round-off.
- The energy-budget residual is O(dt^2).

**Exact homogeneous charged state [D].** At S = 0 the X sector is conformal (pure quartic). The rotating condensate
chi_Y = Y0 e^{i omega tau}, omega^2 = 2 B_X Y0^2, is therefore an exact solution at every a.

- Physically, |Y|^2 = Y0^2/a^2 and n = n_i a^-3. This is exactly the charge_phase phi = 0 branch,
  |X|^2 = (n^2/(8 lambda_X))^{1/3} and omega = n/(2|X|^2), at every a.
- The phi mode obeys delta chi'' + (k^2 + kappa^{2/3} - a^2) delta chi = 0. phi = 0 is stable until a_c = kappa^{1/3},
  where n = n_c, and tachyonic after that.
- V3 checks the orbit (to 3e-5, the integration error). V4 checks the mode equation: the growth factor is 3.52e5 against
  3.52e5 from the ODE, with 1.8% maximum deviation through the crossing.

## 3. Runs [C]

Common settings: 2D, N = 2048, dx = 0.25 (L = 512 comoving), dt = 0.05.

Initial noise (stated spectrum): white noise smoothed by a Gaussian kernel of std s_c = 1 ell (comoving at a = 1),
i.e. P(k) ~ exp(-k^2 s_c^2) and correlation exp(-r^2/(4 s_c^2)). It is normalised to std 1e-2 v in phi (zero velocity)
and 1e-3 v per real component of X. The X noise is given co-rotating velocity.

Resolution of the wall core 2 delta/(a dx) at the end: 3.77 cells (a = 3), 4.53 (a = 2.5), 4.04 (a = 2.8),
7.07 (a = 1.6 frozen), 7.54 (resolution partner). Lump resolution:
- A Q-ball core of radius 3 ell covers at least 4 cells of radius at a = 3.
- Bound-size lumps (R > 1.75 ell) cover at least 2.3 cells of radius.
- The median lump radius is 2.0-2.5 ell at tau_i = 150 and 1.3-1.5 ell at tau_i = 40. Those small, mostly unbound lumps
  are marginally resolved and not converged (§5).

| Name | Preparation | tau_i | a range | seeds | H/m_X at n_c crossing |
|---|---|---|---|---|---|
| A_k2_t40 | (A) homogeneous phi = 0 condensate, n_i = 2 n_c | 40 | 1 -> 3 | 6 | 1.57e-2 |
| A_k4_t40 | (A) n_i = 4 n_c | 40 | 1 -> 3 | 2 | 9.9e-3 |
| A_k2_t150 | (A) n_i = 2 n_c | 150 | 1 -> 2.5 | 2 | 4.2e-3 |
| A_k4_t150 | (A) n_i = 4 n_c | 150 | 1 -> 2.8 | 1 | 2.6e-3 |
| B_t40, B_t150 | (B) neutral X = 0, phi noise about 0 (instantaneous quench at a = 1) | 40, 150 | as A | 6, 2 | - |
| C1_sub_k0.9_t150 | (C) phi = +v side: branch-1 state S0 = 0.811, F^2 = 0.342 at n = 0.9 n_c | 150 | 1 -> 2.5 | 1 | - |
| C2_forced_k4_t150 | (C) phi = +1 with X on its own orbit in the fixed +v vacuum, n = 4 n_c > n_s1 (F^2 = 1.195) | 150 | 1 -> 2.5 | 2 | - |
| C2_forced_k2_t150 | (C) same, n = 2 n_c > n_s1 (F^2 = 0.614 on that orbit) | 150 | 1 -> 2.5 | 1 | - |
| S/CC/DC_k2_t150_s1 | state of A_k2_t150_s1 at a = 1.6, continued with charge (CC) or with X and its momentum set to 0 (DC, "decharged") | 150 | 1.6 -> 2.5 | 1 | - |
| S/CC/DC_k2_t40_s1 | same from A_k2_t40_s1 at a = 1.8 | 40 | 1.8 -> 3 | 1 | - |
| F_k2 / F_dc_k2_t150_freeze1.6 | same state at a = 1.6, expansion then switched off (Minkowski, H -> 0 limit), charged / decharged | 150 | a = 1.6, tau 240 -> 500 | 1 | 0 after freeze |
| F_neutral_t150_freeze1.6 | neutral quench, expansion stopped at 1.6 | 150 | | 1 | |
| R_lo / R_hi_k2_t40 | resolution pair, L = 128: 512^2 dx = 0.25 vs 1024^2 dx = 0.125, same continuum initial field | 40 | 1 -> 3 | 1 | |
| D3_k2_t40 / D3_neutral_t40 | 3D 256^3, dx = 0.3 (L = 76.8), dt = 0.06 | 40 | 1 -> 2.2 | 1 | 1.57e-2 |

The exit times were chosen so that the n_c crossing is slow compared with the tachyonic time: 1/(3H) at the crossing
is 21-126 ell, against the tachyonic time of about 1 ell. They were also chosen so that the exit fits in the run.
Formation happens at a = 1.38-1.88, when n has dropped to 0.49-0.79 n_c because of the finite quench rate.

## 4. Criteria (fixed before the production analysis)

All criteria use the physical fields S = chi_S/a, F^2 = 2|chi_Y|^2/a^2 and n = 2(chi_r pi_i - chi_i pi_r)/a^3.

- **Cells.** Vacuum cells have S > +0.5 ("+") or S < -0.5 ("-"). Core cells have |S| <= 0.5.
- **Core components.** These are connected components of core cells (nearest neighbour, periodic). A component is:
  - **wall-type** if it is adjacent to at least one "+" cell AND at least one "-" cell;
  - **single-vacuum** if it touches only one sign: an isolated Q-ball or Q-matter lump, or a collapsing bubble interior;
  - **restored** if it touches no vacuum cell (before formation the whole box is one restored component).

  A Q-ball sitting on a wall merges with the wall's core and is counted as charge on walls.
- **Wall vs Q-ball surface.** A nearest-neighbour link across which S changes sign is a **wall link** if either holds:
  - at least one endpoint lies in a wall-type component;
  - both endpoints are vacuum cells of opposite sign.

  Sign changes inside single-vacuum or restored components are **internal**, not walls. A Q-ball surface separates a
  vacuum from Q-matter and is never a vacuum sign change. L_wall is the gradient-weighted (PRS/Garagounis-Hindmarsh)
  wall-link length per area. **Bare** length counts only wall links with F^2 < 1 at both ends.
- **Vacuum-basin interface L_vor** (added during the work; see §8, item 1). Every core cell is assigned to the nearer
  vacuum sign (Euclidean distance, periodic). L_vor is the length of the interface between the + and - basins. It counts
  each wall exactly once, including stretches covered by Q-matter, and it ignores sign oscillations inside Q-matter.
  Validation:
  - two walls with an oscillating Q-ball on one of them give 127.7 for exact 128, where L_wall gives 145;
  - circles R = 8 and 16 are reproduced to 2.9e-3 / 1.8e-3 (L_wall: 5e-4 / 1e-4).

  The decharged and frozen comparisons and the re-run charged runs use L_vor.
- **Area parameter** (network_sim convention): A = L tau/2 (= rho_wall t/sigma in the radiation era). The comoving domain
  scale is xi_c = 1/L, the physical one a xi_c. The **coarsening rate** d xi_c/d tau is the peculiar coarsening speed in
  units of c. It is fitted linearly from tau_form + 5 up to xi_c = L_box/8 (box-safe window).
- **Percolation.** A vacuum cluster (periodic) **wraps** if it winds around the torus. This is detected with integer shift
  potentials on the boundary-link graph; network_sim's projection test is also stored. In 2D at p = 1/2 both signs can
  wrap only as stripes. Topology classes:
  - both_wrap (stripes);
  - one_sign_wraps_both_axes (the other sign forms islands only);
  - one_sign_wraps_one_axis;
  - no_wrap.
- **Q-lumps.** These are connected components of F^2 > 1, the cells where phi = 0 is the local minimum at fixed |X|.
  - Lump charge in 2D is per physical ell of z-length, in units of sqrt(lambda_phi) v^3 ell^3 = 1/lambda_phi =
    100 X quanta. In 3D it is the total charge.
  - Lumps with Q > 1 are counted. Lumps with Q > 17 are called **bound-size**: a 2D tube of radius > 1.75 ell, with
    thin-wall E/Q = m* + 2 sigma_Q/(n_* R) < m_X using charge_phase's m* = 0.5623, sigma_Q = 0.68, n_* = 1.778.
  - A lump is **attached** if it overlaps a wall-type component.
- **Formation.** For A and B runs: the first output with f+ + f- >= 1/2. For memory controls: the first output after the
  field has been restored (core >= 1/2) and both vacua again hold >= 5%. **p(+)** = f+/(f+ + f-).

## 5. Validation and conservation [C] (results/validation.json, summary.json)

- **Charge.** The largest |Q(tau)/Q0 - 1| over all 36 production runs is 3.9e-16, i.e. round-off. This includes 3D,
  nonlinear phases, the Minkowski continuations and 1e4-step runs. Small validation runs at strong nonlinearity give
  0 and 1.8e-16.
- **Energy budget.** The residual H(tau) - H0 - int dH/dtau is at most 4.1e-3 of the kinetic + gradient energy, over
  all runs and outputs. In 2D validation it converges as dt^2: 1.82e-3 -> 4.54e-4 when dt is halved, ratio 4.007.
- **Exact branch-2 orbit (V3).** It stays on the branch to 3e-5 with S = 0 exactly. The relation
  |Y|^2 = (n^2/8B)^{1/3} with n = n_i a^-3 holds to 8e-5 for kappa = 2 and 4.
- **Linear tachyonic growth through the n_c crossing (V4).** The growth factor is 3.521e5 against 3.523e5 from the
  ODE.
- **Z2 mirror (V7).** Initial data with phi -> -phi gives exactly the mirrored evolution: max |S + S_mirror| = 0, and
  f+ and f- are exchanged.
- **Labelling.** Stripes, a blob through a corner and diagonal winding are classified correctly. The synthetic
  wall + Q-ball configuration is classified correctly: 2 walls, one attached and one free lump.
- **Resolution (R_lo vs R_hi, dx halved, same continuum field, L = 128).** Mean / max relative differences after
  formation:

  | Quantity | Mean | Max |
  |---|---|---|
  | L_wall | 1.2% | 4.4% |
  | f+ | 0.6% | 1.8% |
  | f- | 0.3% | 0.9% |
  | charge on walls | 3.8% | 11% |
  | total core charge | 2.7% | 10% |
  | Q-lump charge fraction | 10% | 67% |
  | Q-lump number | 9% | 37% |

  The Q-lump rows are small numbers: 3% of the charge and about 20 lumps. a_form is identical. Walls, signs and charge
  on walls are converged at the few-% level. Individual small Q-lumps are not.

## 6. Results [C]

### 6.1 Formation in the charged exit (question 1)

| Group (2D, 2048^2) | runs | a_form | n at formation / n_c | p(+) at formation (mean +- std) | Topology at end (2D torus) |
|---|---|---|---|---|---|
| A_k2_t40 | 6 | 1.55 | 0.54 | 0.5005 +- 0.0037 (0.497-0.505) | 4 both_wrap, 2 one-sign-wraps-both-axes |
| A_k4_t40 | 2 | 1.85-1.88 | 0.61-0.63 | 0.4997 +- 0.0025 | 1 / 1 |
| A_k2_t150 | 2 | 1.383 | 0.76 | 0.501 +- 0.012 (0.493, 0.510) | 2 both_wrap |
| A_k4_t150 | 1 | 1.717 | 0.79 | 0.4986 | one-sign-wraps-both-axes |
| B_t40 (neutral) | 6 | 1.15 | - | 0.4998 +- 0.0035 | 2 both_wrap, 4 one-sign |
| B_t150 (neutral) | 2 | 1.05 | - | 0.504 +- 0.006 | 2 one-sign |
| D3_k2_t40 (3D) | 1 | 1.594 | 0.49 | 0.514 | both signs wrap at 100% of outputs after formation |
| D3_neutral_t40 (3D) | 1 | 1.198 | - | 0.500 | both signs wrap at 100% of outputs |

- **Formation.** The charged exit forms a random-sign domain pattern with p(+) = 1/2 within the seed scatter. The
  physical domain size at formation is xi = 21-29 ell, larger than the neutral quench's 9-13 ell because of the slow
  spinodal.
- **3D percolation.** In the 3D run both vacua percolate: each sign is a single wrapping cluster holding 54% / 41% of
  the volume at the end. This is the 3D statement at p = 1/2, above the site threshold of about 0.31.
- **2D percolation.** In 2D, p = 1/2 is the critical point. Both signs wrap only as stripes, and which configuration
  occurs is a lottery: charged 7 of 11 runs, neutral 2 of 8. The charged runs end more often in the stripe class,
  plausibly because they coarsen more slowly and stay nearer the critical pattern. This is not a 3D statement.
- **Charge at formation.** Most of the charge is not in Q-matter when phi breaks. The local condensate drops onto the
  phi != 0 branch, F^2 of about 0.05-0.4. In the simulation that branch only fragments partly before expansion dilutes
  it.
  - Vacuum fraction at the end: tau_i = 40 about 81%; tau_i = 150 about 65%; Minkowski continuation 44% and still
    falling.
  - Charge on walls at the end: 15-17% (tau_i = 40), 25-27% (tau_i = 150), 23-27% (kappa = 4), 28% (3D).
  - Charge in Q-lumps at the end: 5% (k2 t40), 15% (k4 t40), 30% (k2 t150), 40% (k4 t150), 22% (3D).
  - 60-84% of the lump charge sits in lumps attached to walls.
  - The slower exit (larger tau_i, closer to the physical rate) puts more charge into Q-matter and onto walls.

### 6.2 Persistence: does the charged network persist, decay or freeze? (question 2)

- **No decay.** Every charged-exit run, and every decharged and frozen run, had wall links at every output after
  formation. The end-of-run comoving domain scale was xi_c/L_box = 0.06-0.16 (2D, still box-safe). p(+) at the end
  ranges 0.453-0.587.
- **No freezing.** xi grows steadily in all charged runs. Box-safe coarsening rates d xi_c/d tau (L_vor where
  available, L_wall otherwise):

  | Group | Charged rate | Neutral control |
  |---|---|---|
  | A_k2_t40 | 0.536 (L_vor, 4 runs); 0.565 +- 0.036 (L_wall, 6 runs) | B_t40: 0.727 (L_vor, 2); 0.753 +- 0.042 (L_wall, 6) |
  | A_k4_t40 | 0.420 (L_vor); 0.436 (L_wall, 2) | |
  | A_k2_t150 | 0.396 (L_vor, s2); 0.365 (L_wall, 2) | B_t150: 0.707 (L_wall, 2) |
  | A_k4_t150 | 0.369 (L_wall) | |

  There is no slowing within a run. The first-half / second-half rates are:
  - k2 t150: 0.33/0.40 and 0.36/0.43;
  - k2 t40: 0.53-0.60 / 0.50-0.67;
  - k4 t40: 0.43/0.43.

  L_vor/L_wall at the end is 0.98-1.00 for tau_i = 40 and 0.91 for tau_i = 150 (expanding). It is 0.51 in C2, where
  walls carry most of the charge, and 0.31 in the frozen run.
- **Clean charge effect (same phi configuration, with or without X).**

  | Pair (start state) | Charged d xi_c/d tau | Decharged d xi_c/d tau | Ratio | xi_c at end, charged / decharged |
  |---|---|---|---|---|
  | CC/DC_k2_t150 from a = 1.6, expanding to 2.5 | 0.416 | 0.820 | 0.51 | 87.8 / 138.3 |
  | F_k2/F_dc frozen at a = 1.6, 416 ell of physical time | 0.320 | 0.616 | 0.52 | 113.9 / 180.5 |
  | CC/DC_k2_t40 from a = 1.8, expanding to 3 | 0.525 | 0.688 | 0.76 | 43.9 / 51.7 |

  Fit windows: the whole continuation. In the box-safe window (xi_c <= 64) the tau_i = 150 pair gives 0.388 / 0.673
  (ratio 0.58). The slowdown grows with the charge carried by walls: ratio 0.76 at about 16% on walls
  (tau_i = 40), 0.51-0.58 at 27-37% (tau_i = 150, frozen).

- **Minkowski (H -> 0) continuation.** In the frozen charged run the charge keeps condensing, because nothing dilutes
  it any more:
  - Q-lumps: 35% -> 54% of the charge; attached lumps 21% -> 36%; charge on walls 28% -> 37%.
  - The largest lump reaches Q = 2888 (2.9e5 quanta per ell of length) and the median lump radius is 2.5 ell.
  - Charge per unit wall length rises about 4.8-fold. Even so, xi_vor grows linearly (31 -> 114 comoving) at a constant
    rate.
  - Both signs keep wrapping, with f+/f- = 0.43/0.50.

  The walls are dressed with Q-matter lenses. This is visible in the field maps and agrees with the partial-wetting
  estimate: cos theta = sigma_wall/(2 sigma_Q) = 0.693 [E].
  - The dressing halves the coarsening speed.
  - It does not pin the network on the simulated time scale.
  - Its effect does not grow as the charge per wall length grows: the rate is the same in the first and last thirds of
    the run.
- **Measurement artifact found and corrected.** In the frozen charged run, L_wall stops falling (0.035 -> 0.028) while
  L_vor falls 3.7-fold. L_wall counts the sign oscillations of phi inside the Q-matter attached to walls, so it alone
  would have suggested "freezing". The bare wall length (no Q-matter) fell 7-fold. The freezing reading is wrong; see
  §8.
- **Area parameter.** At the end, A = 1.25-1.39 (k2 t40), 2.0 (k4 t40), 2.3-2.4 (k2 t150) and 3.1 (k4 t150), against
  0.93-1.06 (B_t40) and 0.68-1.17 (B_t150, box-affected). network_sim's 2D physical-equation scaling value is 0.71.
  - These are not scaling values. Both networks are still sub-horizon: their comoving domain scale is 0.2-0.7 of the
    comoving horizon tau.
  - The charged ones form later and coarsen more slowly, so they carry more wall area per Hubble volume at equal tau.
- **Neutral comparison.** The neutral networks reproduce network_sim's 2D behaviour [Q], including the finite-box
  take-over of one phase in 2D: B_t150_s1 reaches p(+) = 0.87 at xi_c/L = 0.53. The frozen neutral run reaches 0.96.

### 6.3 Charge bias (question 3)

- **At formation.** Over 11 charged 2D runs p(+) = 0.5003 +- 0.0046, against 0.5009 +- 0.0042 over 8 neutral runs. The
  3D run gives 0.514.
- **At the end.** Charged runs give 0.501 +- 0.043. Neutral runs give 0.586 +- 0.141, which is 2D take-over in a finite
  box, not bias.
- The charge does not select a sign. The rotation sense of X is irrelevant (exact Z2 x U(1) commutation; charge_phase
  §4). The exact mirror test is V7.

### 6.4 Sign-memory controls (C)

- **C1 (n = 0.9 n_c, phi on the +v branch, Codex-like H2).**
  - min f+ = 0.988, max f- = 0, and no sign change at any output.
  - The charge fragments transiently: 23% in lumps at a = 1.5, falling to 2% by a = 2.5 as small lumps dissolve in the
    +v vacuum.
  - **Memory is kept.** This agrees with charge_phase H2 and with Codex's 3D runs [Q].
- **C2, n = 4 n_c > n_s1 put into the +v vacuum** (X on its own orbit at S = 1; the phi != 0 homogeneous state does not
  exist).
  - S leaves +1 and oscillates coherently through zero with amplitude about 1: mean S = -0.61 and f- = 1.0 at some
    outputs. The output spacing of 2.5 conformal units does not resolve the period, so oscillations were not counted.
  - The oscillation decoheres by a = 1.25 (about 40 ell). Domains of both signs then form
    by phase separation into Q-matter and vacuum, while n is still above n_c.
  - p(+) at formation: 0.466 (s1) and 0.498 (s2). At the end: 0.38 and 0.55.
  - This network is the most charge-laden one: 54-60% of the charge on walls and 76% in Q-lumps. It coarsens most
    slowly: d xi/d tau = 0.17-0.19 from L_wall, 0.29 / 0.30 from L_vor.
  - **Memory is lost.**
- **C2, n = 2 n_c > n_s1** (F^2 = 0.61 on the +v orbit).
  - The same coherent oscillation leaves a bias: p(+) = 0.71 at formation.
  - The minority (-) phase then disappears, as in network_sim's biased runs: f- = 0.16 -> 0.02 by a = 2.5.
  - **Memory is partly kept**, through the phase of the coherent oscillation at decoherence.

  The C2 outcome is therefore not universal. It depends on n/n_s1 and on how the oscillation decoheres, which in turn
  depends on the noise and the expansion rate. charge_phase's "expected to decohere by parametric resonance" holds; the
  sign it leaves can be biased or unbiased.

## 7. Physical scales vs simulated [D]

At v = 1 TeV and T_SM = 200 TeV with g* = 106.75:

- H = sqrt(pi^2 g*/90) T^2/M_P = 5.62e-8 GeV;
- H/m_X = 5.6e-10 and H/m_phi = 4.0e-10;
- at T_n(r = 1) = 198.2 TeV, H/m_X = 5.5e-10, which matches charge_phase.

The simulated H/m_X at the n_c crossing is 1.6e-2 down to 2.6e-3. The simulation therefore expands
**4.7e6 to 2.8e7 times faster** than the physical exit.

- The quench time tau_Q = 1/(2H) is 9e8 ell physically, against 32-190 ell simulated.
- The KZ estimate xi_KZ ~ tau_Q^{1/3} is about 960 ell physically [E], against 3-6 (x log factors) simulated. Measured
  physical domain scales at formation are 21-29 ell.
- The physical dynamic range from formation to the horizon is about 1e9 ell / 1e3 ell. The simulated range is a factor
  2-8 in xi.

Consequences of the fast expansion:

- In the simulation the n_c crossing overshoots to n = 0.54-0.79 n_c before phi becomes nonlinear. Physically the
  overshoot is negligible.
- Dilution competes with branch-1 fragmentation. Physically the dilute vacuum condensate has about 1e9 ell to fragment
  completely into Q-balls, which the Minkowski continuation approximates.

The trend from tau_i = 40 to 150, and then to frozen, is consistent throughout: more charge goes into Q-matter and onto
walls, and the coarsening speed falls from about 0.57 to 0.37-0.42 and then 0.32, against decharged values of
0.62-0.82. No trend toward zero speed is seen within a run.

## 8. What failed, was corrected, or is limited

1. **L_wall is contaminated by Q-matter on walls (corrected).**
   - In the frozen charged run, L_wall plateaued at about 0.027 for 250 conformal time units. This looked like pinning.
   - Snapshots showed walls covered with Q-matter whose interior S oscillates through zero; those zero crossings sit in
     wall-type components and were counted as wall.
   - Fix: the L_vor diagnostic was added (validated, §5). The frozen run, and new decharged, continuation, C2 and A
     replicas, were rerun with it. A, C2 and D3 runs computed before the fix carry L_wall only.
   - In the expanding runs L_vor/L_wall is 0.91-1.00, so L_wall-based rates there are within about 8% of the L_vor
     rates.
     L_wall-based quantities from the frozen run, such as "xi_phys" and beta in summary.json per_run.F_k2...freeze,
     are not used.
   - The first frozen run (run from scratch, L_wall only) was replaced by the run from the saved state. On the 105
     common outputs their time series agree exactly (difference 0).
2. **The 2D percolation question is ill-posed at p = 1/2**: it is the critical point. Both-sign percolation is
   established only in 3D, with one 256^3 run with L/xi of about 4-7 at formation; that is a small box.
3. **Finite box and range.**
   - 2D: L = 512, xi_c up to 0.16 L for charged runs and up to 0.53 L for neutral tau_i = 150 runs, which are
     box-affected late.
   - 3D: L = 76.8, xi_c/L = 0.27 at the end.
   - The coarsening range is a factor 2-8 in xi; scaling is not reached.
4. **2D z-invariant geometry.** Q-balls are infinite Q-tubes (no Rayleigh-Plateau or 3D lens physics) and walls are
   straight in z. The 3D run (true Q-balls with Q up to 1.7e5 quanta) behaves the same way qualitatively. Its 3D
   coarsening rate was not compared with L_vor (the 3D re-run was dropped for time).
5. **Q-lump statistics are not converged** (resolution check: 10% mean, 67% max on a 3% quantity). Lump counts depend
   on the F^2 > 1 threshold. Charge on walls and wall measures are converged.
6. **The A-run area parameter is not a scaling A.** Comparing it with network_sim's A = 0.71 (2D) / 0.8 (3D) is not
   meaningful, because the domains are far sub-horizon.
7. **Classical fields only.** There is no thermal hidden bath: relic_yield's restoring history leaves T_h of order
   0.15 v at the exit, which this module does not model (charge_phase: n_local(T_h = 0.2 v) = 0.034 v^3). There are no
   quantum fluctuations beyond the 1e-2 seed, no D, no wall gravity or backreaction, and no SM coupling (none exists in
   the action).
8. **Not done:**
   - long-time behaviour beyond about 400 ell of Minkowski time;
   - whether the Q-matter lenses eventually pin or detach at much larger charge per wall area;
   - 3D L_vor coarsening;
   - more C2 seeds and n values;
   - initial charge supplied as a thermal plasma or as pre-existing Q-balls rather than a condensate.

## 9. Conditional history constraints (no prediction, no exclusion, no validation)

- **E1 (charged spinodal exit forms a network).**
  - Condition: the hidden charge is a homogeneous condensate on the phi = 0 branch (n > n_c, e.g. relic_yield R3 /
    charge_phase H1 preparation; outside-action data) and is diluted through n_c by expansion.
  - Consequence (classical lattice, D = 0, H/m about 1e7 times the physical value): random-sign domains form with
    p(+) = 1/2 within 0.004 (2D seed scatter), and both vacua percolate in 3D.
  - The network persists over the simulated range with walls dressed by 15-37% of the charge (bare cores plus attached
    Q-balls). It coarsens 1.3-2 times more slowly than the same configuration without charge, in the expanding runs and
    in the H -> 0 continuation. It neither decays nor freezes.
  - Hence the wall-network consequences analysed in network_sim C1 (unbiased preparation) are not removed by the charge.
    The charged network carries more wall area per Hubble volume at a given time than the neutral one in these runs.
- **E2 (sub-critical charge on a selected vacuum).** Charge supplied below n_c onto a +v background (charge_phase H2;
  Codex's preparation) creates no opposite-vacuum region. This was computed once with n = 0.9 n_c.
- **E3 (forced restoration of a selected vacuum).** Charge forced above n_s1 into a +v vacuum does not guarantee sign
  memory.
  - At 4 n_c the memory was lost (2 seeds, persistent network).
  - At 2 n_c a coherent-oscillation bias of p(+) = 0.71 survived and the network decayed.
  - Whether a given forced preparation keeps its sign is a property of that preparation, not of the action.
- **Free inputs.**
  - Charge form and spatial homogeneity, n_i/n_c and its timing (outside action).
  - Initial noise spectrum and amplitude (1e-2 phi, 1e-3 X, s_c = 1 ell).
  - The expansion rate (simulated, not physical).
  - D = 0 and no bath.
  - Benchmark couplings and v (dynamics scale-free in ell units; v enters only through the physical H/m).

## 10. Files

- `code/cesim.py`: integrator, initial states, periodic labelling with winding detection, all diagnostics including
  L_vor, checkpoint/restart and state save/load.
- `code/runs.py`: the run list and worker.
- `code/validate.py`: V1-V7.
- `code/analyze.py`: per-run metrics, comparisons, group aggregates.
- `results/runs/*.json`: per-run config and compact time series (floats to 6 significant digits).
- `results/validation.json`.
- `results/summary.json`: per-run metrics, charged-vs-neutral and decharged comparisons, resolution, sign statistics,
  physical scales, compact tables.
- `results/headline.json`: group aggregates.
- No PNGs: matplotlib is not installed. No field dumps.
