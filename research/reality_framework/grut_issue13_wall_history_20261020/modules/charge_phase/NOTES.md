# Module charge_phase: phase structure of the unchanged action at nonzero X charge

21 October 2026. Claude module for the issue #13 wall-history packet. The action is unchanged (BRIEF.md;
GRUT_CANDIDATE_1.md §2). This module adds no bias, portal or inflaton coupling. Inflation, reheating and charge
generation are **outside-action preparation data**. Nothing here validates or excludes Candidate 1. The statements
below are conditional history constraints.

Every number is marked **derived** (closed form), **computed** (numerical, from the code here), **estimate**
(order-of-magnitude or model-dependent), or **quoted** (from the brief, an archived packet or Codex's commits,
read with `git show`).

## Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/charge_phase/code
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py      # about 5.3 minutes total, deterministic
```

This runs `t0_phase.py` (22 s), `thermal_phase.py` (129 s), `z2_memory.py` (160 s), `cosmo_map.py` (8 s) and
`codex_consistency.py` (1 s). Each writes `results/<name>.json`. Only numpy and scipy are used, with no random
numbers. `codex_consistency.py` reads the archived soliton profile at
`grut_candidate_1_20261018/code/results_profile.csv` without modifying it. Codex's numbers come from
`git show 135d535:.../grut_issue13_response_20261020/{README.md,RESULTS.json}` and
`git show 479f111:.../grut_3d_empirical_20261020/{MATHEMATICS.md,REPORT.md,protocol.json}`. Codex's solver
was not rerun.

## Conventions

These follow Codex's 3D variables, so the results can be compared with them directly:
`S = phi/v`, `x = F^2 = 2|X|^2/v^2`, `a = gX^2/lam_phi = 1`, `b = lam_X/lam_phi = 0.1`.

| Quantity | Unit |
|---|---|
| Length and time | `ell = 1/(sqrt(lam_phi) v)` = 1/(100 GeV) = 1.973e-3 fm |
| Energy density | `lam_phi v^4` |
| Charge density | `sqrt(lam_phi) v^3` = 1e8 GeV^3 at v = 1 TeV |
| Frequency | `sqrt(lam_phi) v` = 100 GeV |

At the benchmark, `m_X = g_X v = sqrt(a) sqrt(lam_phi) v` = 100 GeV, so `ell = 1/m_X`. The dimensionless
potential is

`Vbar(S,x) = (S^2-1)^2/4 + a S^2 x/2 + b x^2/4`.

A homogeneous rotating condensate `X = |X| e^{i w t}` has `n = 2 w |X|^2`. In the units above this is
`nbar = Omega x`, with `Omega = w/(sqrt(lam_phi) v)`. Physical values convert as `n = 0.1 nbar v^3`,
`rho = 0.01 rhobar v^4` and `rho/n = rhobar/nbar x m_X`. Part 3 works in units of v.

---------------------------------------------------------------------------------------------------------------

## 1. T = 0 homogeneous charged condensate (D = 0): exact phase structure

**Derivation.** At fixed charge density the kinetic energy is `n^2/(4|X|^2)`, so

`rhobar(nbar) = min_{S,x} [ nbar^2/(2x) + Vbar(S,x) ]`.

First minimise over S at fixed x (fixed |X|): `dV/dS^2 = 0` gives `S^2 = 1 - a x` when `a x < 1`, and `S = 0`
otherwise. Writing `W(x) = min_S Vbar` gives:

- **Branch 1** (`phi != 0`, `x < 1/a`): `W = a x/2 - (a^2-b) x^2/4`. W is concave because `a^2 > b`
  (binding class `gX^4 > lam_phi lam_X`).
- **Branch 2** (`phi = 0`, `x >= 1/a`): `W = 1/4 + b x^2/4`. W is convex.

W is C^1 at `x = 1/a`. Stationarity in x gives `Omega^2 = 2 W'(x)`. This reproduces the X equation of motion.

| | `Omega^2` | `nbar` | `rhobar` |
|---|---|---|---|
| Branch 1 | `a - (a^2-b) x` | `x sqrt(a - (a^2-b) x)` | `a x - (3/4)(a^2-b) x^2` |
| Branch 2 | `b x` | `sqrt(b) x^{3/2}` | `1/4 + (3/4) b^{1/3} nbar^{4/3}` |

The second x-derivative at fixed n is `2W'/x + W''`. On branch 1 this equals `a/x - (3/2)(a^2-b)`, so branch 1
is a local minimum only for `x < x_s1 = 2a/(3(a^2-b))`. On branch 2 it equals `3b/2 > 0`.

**Results** (derived; a brute-force 2D Nelder-Mead minimisation agrees to a maximum relative deviation of
1.2e-10 over 35 densities):

| Quantity | Value | Meaning |
|---|---|---|
| `\|X\|_c^2` | **0.5 v^2** (`F^2 = 1`) | Brief's estimate **confirmed**. At fixed \|X\| the phi transition is continuous: `S^2 = 1 - F^2`. |
| `n_c = sqrt(8 lam_X) (lam_phi/(2 gX^2))^{3/2} v^3` | **0.0316228 v^3 = 3.16228e7 GeV^3** at v = 1 TeV | Brief's 0.0316 **confirmed**. phi = 0 is a local minimum only for n >= n_c, so this is the **lower spinodal** of the phi = 0 phase. |
| Branch-1 radial spinodal | `x_s1 = 0.740741`, `n_s1 = gX v^3 x_s1/sqrt(3)` = **0.0427667 v^3** | The phi != 0 homogeneous state ceases to exist for n > n_s1. |
| Homogeneous global switch | **n_t = 0.0401026 v^3** | `S` jumps from **0.6480** (x = 0.5801) to **0** (x = 1.1716) |
| Q-matter density (binodal) | `n_* = sqrt(lam_phi) v^3 (lam_phi/lam_X)^{1/4}` = **0.177828 v^3** | `x_* = 1/sqrt(b) = 3.16228`, `\|X\|_* = 1.2574 v`. The numerical convex hull gives 0.17783. |
| Grand-canonical T = 0 transition | at `mu = m* = (lam_phi lam_X)^{1/4} v` = **56.234 GeV** | First order: n jumps from 0 to `n_*`. Computed `p(mu)` matches `max(0, Omega^4/(4b) - 1/4)` on 21 points. |

**Order of the transition.** In n, the transition is **first order**, not continuous. The point where phi
reaches 0 along the stationary branch (`n_c`) lies on the radially unstable part of branch 1 (`x = 1 > x_s1`).
This leaves a metastability window `n_c < n < n_s1` in which both phases are local minima. phi = 0 becomes the
global homogeneous minimum at `n_t = 0.04010 v^3`, *before* the phi != 0 state loses local stability at
`n_s1 = 0.04277 v^3`. Once inhomogeneous states are allowed, the convex hull of `rho(n)` is the straight line
`rho = m* n` for `0 < n < n_*`. That line is phase coexistence of the vacuum with phi = 0 "Q-matter". The
homogeneous phi = 0 phase is the ground state only for `n > n_*`.

**Energy per charge, `rho/n`** (computed, `results/t0_phase.json: rho_over_n_table`):

| n / v^3 | State | `rho/n` |
|---|---|---|
| → 0 | — | 1.000 m_X (free carriers) |
| 0.0220 (Codex production) | phi != 0 | 0.9442 m_X |
| n_c (stable branch) | phi != 0 | 0.9138 m_X |
| n_c | phi = 0 | 1.0277 m_X |
| n_t | — | 0.8801 m_X |
| 0.1 | phi = 0 | 0.5981 m_X |
| n_* | phi = 0 | 0.56234 m_X = **m*** (minimum) |
| 1.0 | phi = 0 | 0.775 m_X |

`m* = 0.56234 m_X` is the homogeneous minimum; above n_* the energy per charge rises again. Even the
homogeneous phi != 0 condensate has `rho/n < m_X`, because integrating out phi leaves an attractive effective
quartic `b - a^2 < 0`.

## 2. Stability of the homogeneous condensate and fragmentation (FLS / Q-ball criterion)

**Linear analysis** (derived, computed). Perturb about the rotating background in the co-rotating frame, with
`psi = sqrt(x) + p + i q` and `S = S0 + s`, as `e^{lambda t + i k z}`. Write `A = lambda^2 + k^2` and
`m_s^2 = 3 S0^2 - 1 + a x`. The dispersion relation is the cubic

`(A + 2bx) A (A + m_s^2) + 4 Omega^2 (A - k^2)(A + m_s^2) - 4 a^2 x S0^2 A = 0`.

As k → 0 on branch 1, this gives the sound speed

`c_s^2 = -x(a^2-b) / (2 Omega^2 - x(a^2-b))`.

This equals the thermodynamic `c_s^2 = (n/mu) dmu/dn`. It is **negative for every** `0 < x < x_s1`. The
unstable band edge is `(k^2 + 2bx)(k^2 + 2S0^2) = 4 a^2 x S0^2`; the analytic and numerical band edges agree
to 1e-3. On branch 2, `c_s^2 = 1/3` and the phi-fluctuation mass^2 is `a x - 1 >= 0`, so it is linearly
stable.

**Verdict:**

- **Branch 1** (phi != 0, all `n < n_s1`, including all n < n_c) is fragmentation-unstable (FLS/Q-ball type,
  `d mu/dn < 0`). Example growth rates: `gamma_max` = 0.0045, 0.045, 0.112 and 0.226 (in units of
  `sqrt(lam_phi) v`) at x = 0.01, 0.1, 0.25 and 0.5.
- **Branch 2** (phi = 0, `n >= n_c`) is linearly stable.
  - For `n_c < n < n_*` it is metastable to phase separation. Thin-wall bubbles of either vacuum sign would be
    needed.
  - For `n > n_*` it is the ground state.

**Thin-wall vacuum-bubble nucleation in the stretched phi = 0 phase** (**estimate**, T = 0):

- The Q-matter/vacuum interface tension at `mu = m*` is `sigma_Q = 0.6800 sqrt(lam_phi) v^3` =
  **6.80e7 GeV^3**, so `sigma_Q^(1/3)` = 408 GeV. This was computed by 1D relaxation, which satisfies the
  virial check: gradient part 0.33998 vs potential part 0.33999.
- With an O(4) bounce, `S_4 = 27 pi^2 sigma_Q^4/(2 eps^3)` ranges from 2.5e5 at n_c to 7.5e7 near n_*. One
  bubble per Hubble time needs `S_4 ≈ 85`.
- So the stretched homogeneous phi = 0 phase does not nucleate before the barrier disappears near n_c. Its exit
  is spinodal (tachyonic) at ≈ n_c.

Thin wall is not valid close to n_c, so this conclusion is a qualitative estimate.

**Which state does dense charge prefer?**

- For mean density `n < n_*` the minimum-energy state is phase-separated Q-matter. Its energy per charge tends
  to m* = 0.5623 m_X for large lumps, plus surface energy.
  - The thin-wall estimate `E/Q ≈ m* + 3 sigma_Q/(n_* R)` gives 0.8665 m_X at the archived charge
    Q = 39,944 (R = 3.77 ell).
  - The archived solved soliton gives **0.91108652 m_X** (quoted). So thin wall is 4.9% low at that size, as
    expected for a thick wall.
  - The thin-wall binding threshold (`E/Q < m_X`) is Q ≈ 1.3e4 (estimate).
- The homogeneous state has *higher* energy than the archived soliton for n < 0.0324 v^3 (phi != 0 branch). It
  has lower energy only in narrow homogeneous windows (phi != 0 branch for 0.0324 < n < n_s1; phi = 0 branch for
  0.0379 < n < 1.71 v^3). Those homogeneous states are themselves above the Maxwell line `m* n`.
- Conclusion: dense charge (n > n_*) prefers the homogeneous phi = 0 phase. Dilute charge prefers Q-balls
  (localized phi ≈ 0 lumps) in a phi = ±v vacuum.

## 3. Finite hidden temperature with chemical potential (LEADING ORDER)

**Scheme** (`thermal_phase.py`). At fixed `(T, n)`:

`f(phi; T, n) = U_T(phi) + f_gas(n_g; T, M) + E_cond(n - n_g; M)`, with

- `U_T = (lam_phi/4)(phi^2 - v^2)^2 + 1/2 (lam_phi/4 + hD^2/12) T^2 phi^2`. These are the phi and D loops at
  O(T^2) only.
- `M^2 = gX^2 phi^2 + Pi_X T^2`, with `Pi_X = gX^2/12 + lam_X/3 = 1.1667e-3`.
- `f_gas = -P + mu n_g` is the **exact one-loop ideal Bose gas** of X and anti-X with chemical potential
  `|mu| <= M`. Its phi dependence supplies the X-loop thermal mass (`gX^2 T^2/6` at mu = 0) and its mu
  dependence.
- `E_cond(q; M) = min_c [q^2/(4c) + M^2 c + lam_X c^2]` is the homogeneous condensate.
- `n_g = min(n, n_max)`, with `n_max = n(T, mu = M)`. Bose-Einstein condensation occurs when mu reaches the X
  effective mass, and the excess charge then goes into the condensate.

**Approximations, stated:**

1. Zero-temperature (Coleman-Weinberg) pieces are dropped. The archived quantum packet found them at the
   1e-5 level in m* (+1.86e-5, quoted).
2. phi and D cubic/daisy terms are dropped.
3. Gas-condensate interactions beyond `Pi_X` are dropped. This makes `mu_gas = M` while the condensate has
   `mu_c = sqrt(M^2 + 2 lam_X c)`. The mismatch is `O(lam_X c/M)`. The partition is chosen variationally.
4. The full φ-even structure is kept.

**Checks.**

- The T → 0 limit reproduces part 1 exactly (computed): n_local = 0.031623, n_global = 0.04015 on a 1e-4 grid,
  broken spinodal 0.0427, and coexistence [0, 0.1778] at `mu_coex` = 0.056234 v.
- The Bose quadrature agrees with scipy `quad` to ≤ 6e-10 relative. For the `mu = 0` charge the reference
  value is 0 and the agreement is ≤ 2e-17 absolute. The integrand is written without cancellation.

**Zero-charge transition** (computed):

- `T_c0 = 1.55750 v`. The leading-T^2-only value is 1.54765 v.
- The +0.64% shift comes from evaluating the X loop with its thermal mass (`-M T/(4 pi)` term).
- `phi_min(T)` goes to zero continuously: 0.057 at `T_c0 - 0.0025 v`, with a single minimum. So the zero-charge
  transition is **second order (continuous) at this order**.

**Phase diagram in `(T_h, n)`.** Units: T/v and n/v^3. Computed; `results/thermal_phase.json: rows`.

| T_h | n_local (phi = 0 is a local min above this) | n_global (phi = 0 global among homogeneous) | phi != 0 spinodal | Coexistence [n_low, n_high] (binodal) | gas `n_max` at phi = 0 |
|---|---|---|---|---|---|
| 0 | 0.03162 | 0.04015 | 0.0427 | [0, 0.1778] | 0 |
| 0.2 | 0.03409 | 0.04035 | 0.0427 | [0.0006, 0.1759] | 9e-5 |
| 0.5 | 0.03670 | 0.04105 | 0.0428 | [0.0044, 0.1640] | 0.00142 |
| 1.0 | 0.03994 | 0.04215 | 0.0429 | [0.0180, 0.1193] | 0.01139 |
| 1.3 | 0.04141 | 0.04245 | 0.0427 | [0.0302, 0.0797] | 0.02501 |
| 1.5 | 0.04208 | 0.04245 | 0.0426 | [0.0397, 0.0498] | 0.03843 |
| 1.545 | 0.04218 | 0.04225 | 0.0422 | [0.0420, 0.0427] | 0.04199 |
| 1.55 | 0.04077 | 0.04075 | 0.0407 | none | 0.0424 |
| 1.5575 (`T_c0`) | 0 | 0 | — | — | — |

**Reading the diagram:**

- **The critical line `n_c(T_h)` is not monotonic.** It rises from 0.0316 v^3 at T_h = 0 to a plateau of
  ≈ 0.042 v^3 near T_h ≈ 1.5 v.
  - Charge stored in a relativistic thermal X gas barely raises phi's mass, because `<|X|^2>_gas <= T^2/12`
    while `mu <= M`. Only condensed charge (`2 gX^2 c`) restores the symmetry efficiently.
  - The gas's thermal X mass also raises the condensate's cost per amplitude.
- **Near `T_c0`.** n_local falls steeply from ≈ 0.042 to 0 between T_h = 1.545 and 1.5575 v. Some values:
  0.0336 at 1.5535, 0.0278 at 1.5550, 0.0131 at 1.5570.
- **Critical endpoint.** First-order coexistence ends between T_h = 1.545 v and 1.550 v, at n ≈ 0.042 v^3.
  Above it, the n-transition is continuous (n_local = n_global within the grid).
- **Whole picture.** For T_h < T_c0, φ = 0 needs a *condensed* charge density of at least about 0.032–0.042 v^3.
  For T_h > T_c0 the symmetry is restored thermally whatever the charge.

## 4. Z2 sign symmetry and sign memory

**Proof that `p(+) = p(-)`** (derived). Let P be the unitary implementing `phi → -phi`, `pi_phi → -pi_phi`,
which leaves X, D, the SM fields and the metric invariant.

1. The potential U depends on `phi^2` only, and the kinetic and metric couplings are even in phi. So
   `P H P^† = H`, including the time-dependent FRW Hamiltonian.
2. The U(1)_X charge `Q = ∫ j^0_X` involves only X, so `P Q P^† = Q`. Z2_phi and U(1)_X act on different fields
   and therefore commute.
3. Hence `rho_GC = e^{-beta(H - mu Q)}/Z` satisfies `P rho P^† = rho`. Any state built from it by P-invariant
   unitary evolution stays P-invariant: `rho(t) = U rho0 U^†` with `P U P^† = U`.
4. For any Z2-odd observable O, `Tr(rho O) = Tr(P rho P^† P O P^†) = -Tr(rho O) = 0`. P maps the event "sign of
   coarse-grained phi in region A is +" onto the event "... is -". So `p(+) = p(-)`, and the joint sign
   distribution over any set of regions is invariant under a global flip.

A charged symmetric phase (phi = 0) that later breaks therefore yields equal sign probabilities, unless the
initial state is not P-invariant. Two examples of a non-invariant state:

- a +v boundary or background that percolates;
- a coherent phi oscillation remaining from a time when phi = +v had been selected.

The leading-order free energy here is even in phi by construction. The numerical check
`max|f(phi) - f(-phi)| = 0` is trivial.

Equal probabilities do not by themselves give domains. A network needs a sign correlation length smaller than
the restored region (here the Hubble volume). Part 5 estimates that correlation length.

**Sign memory of a finite restored patch** (`z2_memory.py`). This is a **toy model**:

- spherically symmetric 3D phi equation;
- X background frozen as `x(r,t) = x_in(t) w(r)`, with w a tanh ball of radius R and edge width 0.6;
- pre-quench `x_in = 3.4`, a typical Codex fragment interior; the initial S is the static positive solution;
- then x_in drops below 1, either suddenly or linearly over tau_Q;
- "decision" is the sign when `|S(0,t)|` first exceeds 0.5.

All lengths below are in ell = 1/(100 GeV) = 1.973e-3 fm.

- **Dirichlet estimate (derived).** The lowest Dirichlet mode of a sharp ball is unstable iff `R > R_D = pi/m_tach`,
  with `m_tach^2 = lam_phi v^2 - 2 gX^2 |X|^2`, i.e. `(1 - a x)` in units. This gives
  R_D = 3.14 / 4.44 / 7.02 / 9.93 / 31.4 ell for x_after = 0 / 0.5 / 0.8 / 0.9 / 0.99. For full charge loss,
  R_D = 6.2e-3 fm.
- **The boundary keeps imposing its sign even for R > R_D**, because its forcing projects positively onto the
  unstable mode. A coherent negative seed of amplitude eps flips the centre only beyond `R_first_flip`
  (computed, sudden quench):

  | x_after | eps = 1e-2 | eps = 1e-4 | eps = 1e-6 |
  |---|---|---|---|
  | 0 | **8** | 14 | 20 |
  | 0.5 | 12 | 20 | 28 |
  | 0.8 | 16 | 28 | > 28 |

  `R_D` is 3.14, 4.44 and 7.02 for these three rows. A causal heuristic, `R_heur ≈ 2 + ln(0.5/eps)/m_in + ln(1/eps)/m_tach`
  with `m_in = sqrt(a x_before - 1)`, gives 9.1 / 16.7 / 24.3 for x_after = 0. It overestimates by about 15–20%.
- **Phase-randomised seeds** (8 phases, eps = 1e-4, x_after = 0). Fraction of negative centre decisions:

  | tau_Q | R = 4 | R = 8 | R = 16 | R = 32 |
  |---|---|---|---|---|
  | 0 | 0 | 0 | 0.5 | 0.5 |
  | 10 | 0 | 0 | 0 | 0.5 |
  | 100 | 0 | 0 | 0 | 0 |

  Sign memory grows with the quench time, a Kibble-Zurek-like trend. 0 means the boundary sign is always
  inherited; 0.5 means the outcome follows the seed (no memory).
- **Noise amplitude (estimate).** Quantum fluctuations of S coarse-grained on `1/m_tach` give
  `eps ≈ sqrt(lam_phi)/(2 sqrt(2) pi) ≈ 0.011`. So **R_* ≈ 8 ell (sudden charge removal)**, i.e. 0.016 fm.
  Slower removal raises it. A thermal hidden bath raises eps and lowers R_*.
- **Flipped interiors are closed bubbles** (computed). For R = 16, eps = 1e-2 the centre decides negative at
  t = 4.8, spends 13.4 time units below -0.5 (inner-half volume), and ends at S = +1.000: the closed wall
  collapsed. A −v region created inside a connected +v background cannot become a network unless something
  stabilises it. Charged walls are not computed here.
- **Charged-interior control** (computed). With x kept at 3.4 (no quench) and a coherent kick of -0.3 or -0.6,
  S oscillates about 0. The final centre mean is 0.0017 / 0.0004, and S never settles at -v.

## 5. Mapping to cosmology (homogeneous charge with yield Y_Q = r Y_B)

**Quoted inputs:** omega_b = 0.02233, rho_crit/h^2 = 1.05375e-5 GeV cm^-3, m_p = 0.938272 GeV,
s_0 = 2891.2 cm^-3, g_*s = 106.75.

These give `n_b = 2.50783e-7 cm^-3` and **Y_B = 8.6740e-11**. With `s_SM = 46.823 T^3`,
`n_Q = r Y_B s_SM = 4.0614e-9 r T^3`.

Codex's production density 0.022009 v^3 at r = 1 corresponds to T_vis = 175.64 TeV. Codex reports 175.63 TeV,
which agrees to 7e-5.

**T_n(r) for a HOMOGENEOUS charge** (derived/computed, T_h = 0, v = 1 TeV). Scaling:
`T_n = 198.2 r^{-1/3} (v/TeV) TeV` at fixed couplings; in general
`T_n = [sqrt(8 lam_X)(lam_phi/(2gX^2))^{3/2}/(r Y_B 2pi^2 g_*s/45)]^{1/3} v`.

| Threshold | r = 1 | r = 0.64 | r = 0.1 | r = 0.01 |
|---|---|---|---|---|
| n_c (phi = 0 locally stable) | **198.2 TeV** | 230.0 | 427.0 | 920.0 |
| n_t (phi = 0 homogeneous global) | 214.5 | 248.9 | 462.2 | 995.8 |
| n_s1 (phi != 0 branch gone: forced restoration) | 219.2 | 254.3 | 472.2 | 1017 |
| n_* (homogeneous phi = 0 is the ground state) | 352.5 | 409.0 | 759.3 | 1636 |

At T_n(r = 1), `H/m_X = 5.5e-10`.

**Hidden temperature ratio.** With `T_h = xi T_vis` held constant (an outside-action input), the symmetric
phase is exited through the charge line only for small xi. On the sampled xi grid the charge route holds at
xi <= 0.005 and the exit is thermal from xi = 0.0075 (r = 1, 0.64). The corresponding pairs are 0.003 and 0.005
for r = 0.1, and 0.001 and 0.003 for r = 0.01. For larger xi the exit is thermal, at `T_h ≈ T_c0 = 1.5575 v`, where the charge is
negligible: for xi = 0.1, r = 1, n_exit = 1.5e-5 v^3. Full table in `results/cosmo_map.json: xi_trajectories`.

**Kibble-Zurek domain size (estimate, mean-field nu = 1/2, z = 1).** This assumes a homogeneous phi = 0 phase
diluted through n_c by expansion, with `tau_Q = 1/(3H)`. Then
`xi_KZ ≈ (m tau_Q)^{1/3}/m` ≈ 850 ell (r = 1), 300 ell (r = 0.01). The Hubble radius is 1.8e9 ell (r = 1),
giving about 1e19 sign domains per Hubble volume (2e16 for r = 0.01). Both signs then percolate in 3D
(`p = 1/2` exceeds the site-percolation threshold ≈ 0.31). This is a formation statement only;
persistence and backreaction are not computed here.

**When is "homogeneous" the right assumption?**

- Only if an **outside-action** mechanism supplies the charge as a spatially homogeneous condensate or plasma,
  over super-horizon regions, at `T_vis > T_n(r)`. For n >= n_c, the phi = 0 condensate does not fragment
  (`c_s^2 = 1/3`).
- If the charge is supplied below n_c, the condensate is on branch 1. It fragments in ~50–100 ell (Codex: first
  phi < 0 at tau ≈ 72), which is 1e-9 of a Hubble time. The charge is then localized in Q-balls of fixed
  interior density (x ≈ 3–4). Further dilution of the *mean* density no longer re-breaks anything
  homogeneously. Patches of radius ~3–5 ell are inside R_*, and the vacuum between them keeps its pre-existing
  sign.
- Charge already localized in Q-balls before the phi vacuum is selected would require Q-balls to exist in the
  restored phase, where they are not bound. Laine-Shaposhnikov's finite-T equilibrium-with-Q-balls statement
  (abstract read) concerns an infinite system at fixed charge density. It was not applied quantitatively here.

## 6. Consistency with Codex's production replays (135d535, 479f111; not rerun)

**Preparation.** `F^2 = 0.25`, `S = sqrt(0.75)`, `Omega = sqrt(0.775)`, side-64 periodic box (length unit ell,
checked from MATHEMATICS.md: `a = 1` and the S equation fix the unit at `1/(sqrt(lam_phi) v)`), 3D seed of RMS
amplitude 0.005. This is exactly the branch-1 homogeneous state at `nbar = 0.220085`, matching Codex's
0.2200852108. Its density is **n = 0.022009 v^3 = 0.696 n_c**.

- It is fragmentation-unstable: `c_s^2 = -0.1698`, unstable band k < 0.595, fastest growth
  `gamma = 0.1121` at k = 0.404 (wavelength 15.6), 6 unstable harmonics per axis.
- The time to non-linearity, `ln(1/0.005)/gamma = 47`, is consistent in order with Codex's first phi < 0 at
  tau = 72.35–73.35. The seed is spread over modes, and the compression to `F^2 > 1` takes extra time.

**Is F^2 > 1 exactly the region where phi = 0 is the local minimum at fixed |X|?** **Yes** (derived). At fixed
x, `Vbar(S,x)` has minima `S = ±sqrt(1 - a x)` for `a x < 1`, and the **single** minimum S = 0 for
`a x >= 1`. The threshold is `F^2 = 1/a = lam_phi/gX^2 = 1` exactly.

Codex's retained negative endpoint cells (quoted) all satisfy this:

| Grid | Minimum F^2 | F^2 quartiles | Local phi mass^2 at the median (`a F^2 - 1`) |
|---|---|---|---|
| Coarse | 1.072 | 3.15 / 3.65 / 4.00 | 2.65 |
| Medium | 1.287 | 3.04 / 3.47 / 3.93 | 2.47 |
| Fine | 1.788 | 2.99 / 3.44 / 3.99 | 2.44 |

- The every-step minimum of the mean F^2 over negative cells is 2.12 / 1.91 / 1.89.
- The quartiles match Q-matter, `x_* = 3.162`; the medians are 1.09–1.15 x_*. They also match the archived
  soliton centre, F^2 = 3.917.
- **Verdict:** the shallow negative phi (minima -0.2975, -0.2377, -0.2618) are oscillations about the
  phi ≈ 0 minimum inside charge-restored fragments, not opposite-vacuum domains. The charged-interior control
  in part 4 reproduces this behaviour.
- Collisions keeping phi > 0 is also consistent. Static Q-ball centres sit at small positive S (archived
  S(0) = +0.049) because the +v exterior penetrates as `e^{-sqrt(aF^2-1) R}`.

**Patch size vs R_*.**

- The archived Q = 39,944 soliton has its `F^2 > 1` region out to r = 3.38 and S < 0.5 out to r = 3.26
  (R90 = 5.10).
- The negative-cell volumes (1228, 1301 and 1291 ell^3) would make a ball of radius about 6.7 only if merged.
  Codex reports 248–315 components at the 2% threshold, so the real patches are smaller.
- The patches are comparable to the bare Dirichlet `R_D = pi` and **below the toy sign-memory radius**
  (first flip at 8 ell for eps = 1e-2, 14 for 1e-4).
- Directly: the archived profile with its charge removed instantly inherits +v for coherent negative seeds up to
  0.3 (computed).
- These finite patches in a connected +v background could only produce closed bubbles, which collapse.

**Codex-consistency verdict: CONSISTENT.**

---------------------------------------------------------------------------------------------------------------

## Conditional history constraints from this module (no prediction, no exclusion)

- **H1 (charge-forced restoration).**
  - **Condition:** the hidden X charge existed as a homogeneous condensate or plasma with `n_Q >= n_s1 = 0.0428 v^3`
    (cold hidden sector) over super-horizon regions. Equivalently, it was supplied homogeneously at
    `T_vis >= 219 r^{-1/3} (v/TeV) TeV`.
  - **Consequence:** the phi != 0 homogeneous state does not exist, and phi = 0 holds there.
  - **Re-breaking:** spinodal, at `n ≈ n_c` (`T_vis ≈ 198 r^{-1/3}` TeV, cold hidden sector). By the Z2 proof,
    `p(+) = p(-)`, with correlation length ≈ 1e2–1e3 ell, far below the Hubble radius. **A random-sign domain
    pattern spanning the horizon forms.**
  - **Exception:** a coherent sign memory in phi would prevent this. One example is a phi oscillation remaining
    from an earlier +v selection, if charge was added later. This was not computed; such oscillations are
    expected to decohere by parametric resonance at `H/m ~ 1e-9`.
  - The brief's exclusion therefore cannot assume charged preparations avoid walls: those above T_n(r) do not.
- **H2 (sub-critical supply).**
  - **Condition:** charge supplied homogeneously at `n < n_c` (`T_vis < 198 r^{-1/3}` TeV) on a background
    already of one sign.
  - **Consequence:** the charge sits on branch 1, fragments within about 100 ell into Q-balls whose phi ≈ 0
    cores (radius ~3–5 ell) lie inside the sign-memory radius. The charge creates no opposite-vacuum domain.
  - Codex's production preparation is this case. The open inputs are the existence of the single-sign
    background (outside-action: inflationary selection, no hidden reheating above T_c0) and the charge
    production mechanism.
- **H3 (thermal hidden sector).** If `T_h,max > T_c0 = 1.5575 v` (leading order, continuous at mu ≈ 0), the
  symmetry is restored independently of charge. With `T_h/T_vis ≳ 0.0075 r^{1/3}`-ish (table), the exit is
  thermal and the charge plays no role in the sign pattern.

## Free inputs on which these results depend

- **Couplings and scale.** The benchmark `(gX, lam_phi, lam_X, hD, lam_D) = (0.1, 0.01, 0.001, 0.01, 0.001)` and
  `v = 1 TeV` are inputs. All thresholds scale as v^3 (densities) or v (temperatures) at fixed couplings.
- **Charge yield.** `r = Y_Q/Y_B`. The primordial charge yield is not supplied by the action.
- **Charge form and timing** (outside action):
  - how the charge was supplied: homogeneous condensate, thermal plasma, or localized;
  - the visible temperature at supply.
- **Hidden temperature.** `xi = T_h/T_vis`, assumed constant, and `T_h,max` (outside action).
- **Initial vacuum state.** Selection of phi = +v (outside action) and the absence of coherent phi oscillations.
- **Noise.** The noise amplitude eps (quantum estimate ≈ 1e-2) and the charge-loss time tau_Q, for R_*.

## What failed, was corrected, or is limited

1. **First toy criterion failed and was replaced.** The `z2_memory` toy first scored a "flip" by the *final*
   centre sign. It never flagged flips, because flipped interiors collapse within the run. It was replaced by
   the first-decision criterion; the collapse is reported separately.
2. **First ramp test was confounded and was replaced.** It used a coherent zero-velocity seed. For `x > 1` the
   seed oscillates before the tachyon appears, so the outcome depended on its phase. It was replaced by an
   8-phase randomised test.
3. **Path bug.** The first run could not find the archived profile (wrong relative path). This was fixed and
   rerun; no results from the failed attempt are used.
4. **Bose-integrand cancellation.** `scipy quad` in the check diverged at `mu = M` (`eps - mu` evaluates to 0 in
   floating point). Both the check and the production integrand now use `(y^2 + M^2 - mu^2)/(eps + mu)`. The
   full pipeline was rerun from cleared results after the fix. Tc0 changed by 2e-16.
5. **Spherical focusing.** Exact spherical symmetry produces large transient negative S at r = 0 (down to -9)
   after the + front converges. This is an artifact of the toy, not domain formation. Decisions and final
   states are reported instead. The toy has no X dynamics, no non-spherical modes, no expansion and no
   gravity.
6. **Thermal part is leading order only.** It has no phi/D daisy cubic terms, no zero-T loops, and no
   gas-condensate Hartree terms beyond Pi_X. Treat the critical endpoint (1.545–1.55 v) and the plateau shape
   as leading-order. n_global and the spinodal use a 1e-4 v^3 grid in n and a 0.004 v grid in phi.
7. **Not done:**
   - nucleation beyond thin wall;
   - finite-T bounce;
   - charged-wall (Q-balls on walls) dynamics and whether charge stabilises closed bubbles;
   - wall-network persistence and scaling;
   - coherent-oscillation sign memory;
   - a quantitative Laine-Shaposhnikov-type finite-T Q-ball equilibrium;
   - 3D runs.

## Literature (read level stated; no text copied)

| Reference | Read level | Used for |
|---|---|---|
| Felder, Garcia-Bellido, Greene, Kofman, Linde, Tkachev, PRL 87, 011601 (2001), arXiv:hep-ph/0012142 | abstract (fetched) | Spinodal/tachyonic breaking completes within about one oscillation. This supports the random-sign spinodal exit. |
| Laine & Shaposhnikov, NPB 532, 376 (1998), arXiv:hep-ph/9804237 | abstract (fetched) | At finite T, an infinite system at finite charge density is inhomogeneous, with Q-balls. Qualitative support for coexistence. |
| Kasuya & Kawasaki, PRD 61, 041301 (2000), arXiv:hep-ph/9909509 | abstract (fetched) | Condensate fragmentation, with Q-ball size set by the dominant growing linear mode. Supports the linear-growth comparison with Codex's timing. |
| Enqvist & McDonald, NPB 538, 321 (1999), arXiv:hep-ph/9803380 | abstract (fetched) | Context only: condensate-to-B-ball cosmology. No number used. |
| Benson, Bernstein, Dodelson, PRD 44, 2480 (1991), doi:10.1103/PhysRevD.44.2480 | bibliographic only (search result); not read | Context: effective potential at fixed charge; SSB as BEC. |
| Coleman, "Q-balls", NPB 262, 263 (1985); Friedberg, Lee, Sirlin, PRD 13, 2739 (1976); Haber & Weldon, PRL 46, 1497 (1981); Kibble, J. Phys. A 9, 1387 (1976); Zurek, Nature 317, 505 (1985) | not read in this session; identifiers from background knowledge, unverified | Thin-wall `omega_0 = min sqrt(2U/\|phi\|^2)` (rederived here); FLS two-field class; ideal relativistic Bose gas with mu (computed here exactly); KZ scaling (estimate only). |
| Mossa et al., Nature 587, 210 (2020) | quoted via BRIEF.md | omega_b |
