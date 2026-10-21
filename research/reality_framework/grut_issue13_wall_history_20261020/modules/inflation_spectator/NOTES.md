# Module inflation_spectator: initial Z2 vacuum selection of phi by an inflationary era

21 October 2026. Claude module for the issue #13 wall-history packet (BRIEF.md). The action is **unchanged**.
**Inflation, H_inf, N_\*, the inflaton, reheating and the start condition of the observable patch are
outside-action preparation data.** Nothing here adds a bias, portal or inflaton coupling. This module does not
say that Candidate 1 is excluded or validated. It gives conditional history constraints.

Labels: **[D]** derived here, **[C]** computed by the code here, **[Q]** quoted (literature, BRIEF, sibling
module), **[E]** order-of-magnitude estimate.

## Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/inflation_spectator/code
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py     # about 4-5 min, deterministic
```

`run_all.py` runs `xi_running.py`, `checks_1d.py` (cross-checks quoted in section 2), `stochastic_scan.py`
(1D Fokker-Planck scan), `fp2d_checks.py` (phi-|X| Fokker-Planck), `post_inflation.py`, `final_sign.py`,
`classification.py` and `plot_svg.py`. Outputs go to
`results/*.json` and `results/p_minus_vs_Hinf.svg`. Only numpy and scipy are used, with no random draws (one
fixed-seed coupling sample in an algebra check). Wall times are in `results/runtimes.json`.

## 0. Result in brief

1. **Nonminimal running [C].** From xi(mu_0 = 100 GeV) = 0, one-loop running gives xi_phi(mu = H) = 0 at
   H = 1e2 GeV and -2.74e-3 at H = 1e13 GeV. The scale band mu in [H/2, sqrt(12) H] gives -2.66e-3 to -2.87e-3.
   The Hubble-induced mass^2 is **12 xi_phi H^2 = -0.0328 H^2 at 1e13 GeV. It is tachyonic for every H > mu_0**,
   so it deepens the double well. Its magnitude is small (phi stays light).
2. **Stochastic selection [C].** Use the X = D = 0 path with couplings run to mu = H, and the fresh start
   phi = +v at the beginning of the last N_* e-folds. The expected opposite-vacuum volume fraction p_- then falls
   below 1e-3 for H_inf < 368 GeV (N_* = 50) or 355 GeV (N_* = 60). It falls below 1e-10 for H_inf < 239 GeV or
   237 GeV. The expected number of opposite-vacuum regions of any size in the observable patch is below 1e-3 for
   H_inf < 144 GeV or 138 GeV. Two variants move these thresholds by at most 20%. A relaxed + well start gives
   324, 233 and 144 GeV (N_* = 50). The X-fluctuation (entropic) correction gives 295, 225 and 143 GeV. The
   threshold for no opposite-vacuum domain that re-enters after T = 1 MeV is 183-197 GeV.
3. **H_inf >> v [C, D].** P(phi) = P(-phi) exactly. For a fresh start at +v, p_- -> 1/2: it is 0.4996 at
   1e6 GeV and 1/2 - O(1e-10) at 1e13 GeV. For equilibrium draws the ensemble mean is exactly 1/2. A single
   observable patch has a minority fraction with median 0.07-0.20 when X is frozen at zero, and 0.38-0.41 when X
   fluctuates (2D). After inflation the oscillating condensate randomizes the sign wherever H_inf >~ 14 v. For
   H_inf >= 5 TeV the post-trapping fraction is 0.4997-0.4999.
4. **Classification.** A fresh +v start gives four regimes: (A) no opposite region expected (H_inf < ~140 GeV);
   (B) rare closed minority bubbles; (C) closed minority domains with 1e-3 < p_- < p_c; and (D) both vacua
   percolate, so a macroscopic wall network forms at horizon re-entry. Regime D starts at H_inf > 1.78 TeV
   (N_* = 50) or 1.61 TeV (N_* = 60) at the end of inflation, and at 1.43 TeV or 1.35 TeV after post-inflation
   trapping. For H_inf >> v the result is p_- = 1/2: a percolating network forms, unless the history is changed
   by data outside this module. Collapse of minority domains that re-enter late carries real consequences; see
   section 3.
5. **Post-inflation [C, E].** The condensate makes N_cross = 0.31 phi_end/v zero crossings [C]. Preheating of X
   (q = g_X^2/lambda_phi = 1) has Floquet mu_max = 0.1470 [C]. Its back-reaction needs only about 7-9 crossings
   [E], so for H_inf >~ 20 v nonthermal restoration is plausible. The criterion is <dphi^2> > v^2/3 or
   g_X^2 <|X|^2> > lambda_phi v^2/2.

The inflationary selection matters **only** if the hidden sector is never heated above T_c ~ 1.56 v afterwards
and suffers no nonthermal restoration. With all portals zero at mu_0, a visible-only reheating leaves the hidden
sector cold (section 4.5). That cold-hidden-sector chronology is a preparation assumption, not a result.

## 1. Nonminimal couplings: matching, running, Hubble mass

**Convention [D].** Metric (-,+,+,+). The potential gains (1/2) xi_phi R phi^2 + xi_X R |X|^2 + (1/2) xi_D R D^2.
With X = (q1 + i q2)/sqrt2, every real component carries (1/2) xi R z^2. GRUT_CANDIDATE_1.md §3 writes
"xi_phi R phi^2" without the 1/2. Its coefficient is xi_phi/2 here; both vanish at mu_0. In de Sitter space
R = +12 H^2, so V gains 6 xi_phi H^2 phi^2 and the Hubble mass^2 is m_H^2 = 12 xi_phi H^2. The conformal value is
xi = 1/6.

**Derivation [D].** Write V4 = Lambda_abcd z_a z_b z_c z_d/4!. On a curved background the one-loop
heat-kernel/MSbar divergence of the scalar loop depends on the field-dependent mass matrix
H_ab(z) + (xi_ab - delta_ab/6) R (cf. 1407.3141 eq. 2, theta_i = xi - 1/6 for scalars). The R-linear,
field-dependent pole is therefore the mass pole of rg_beta_functions.txt with m^2_cd -> (xi_cd - delta_cd/6) R:

    16 pi^2 beta_{xi_ab} = sum_cd Lambda_abcd (xi_cd - delta_cd/6).

For the Candidate 1 tensor (RG_DERIVATION.md §2) this gives

    16 pi^2 beta_xi_phi = 6 lp (xi_phi - 1/6) + 4 g_X^2 (xi_X - 1/6) + 2 h_D^2 (xi_D - 1/6)
    16 pi^2 beta_xi_X   = 2 g_X^2 (xi_phi - 1/6) + 8 lX (xi_X - 1/6) + 2 alpha (xi_D - 1/6)
    16 pi^2 beta_xi_D   = 2 h_D^2 (xi_phi - 1/6) + 4 alpha (xi_X - 1/6) + 6 lD (xi_D - 1/6)

These are the same coefficients as B_{m^2} in rg_beta_functions.txt. The anomalous dimension vanishes at one
loop (RG_DERIVATION.md §2).

**Checks [C]** (`results/xi_running.json` → `checks`):
- The closed forms agree with the tensor contraction to 2.8e-16 (relative), with Lambda taken from exact
  finite-difference fourth derivatives of the quartic polynomial. This holds at the benchmark and at five random
  coupling/xi points.
- The curvature-linear RG identity sum_ab beta_xi_ab z_a z_b/2 = Tr[H4(z) Xi]/(16 pi^2) holds to 3.2e-16 at 20
  random field/coupling points. (The first version of this check failed by an exact factor 2. The bug was the
  normalization inside the check, not in the betas; see section 7.)
- For an O(4) multiplet with V = lambda (h^2/2)^2 the formula gives 12 lambda (xi - 1/6). This is the scalar term
  of Herranen-Markkanen-Nurmi-Rajantie 1407.3141 eq. (7), 16 pi^2 beta_xi = (xi - 1/6)(12 lambda + 6 y_t^2 - ...)
  [Q, read in the text of the paper].

**Numbers [C].** At mu_0: B_xi = (-0.0167, -0.004667, -0.001033) for (phi, X, D), so
beta_xi_phi = -1.0575e-4 per e-fold of mu (leading log: xi_phi ~= -1.0575e-4 ln(mu/mu_0)). Full one-loop running
of all twelve parameters, mu = H:

| H (GeV) | xi_phi | xi_X | xi_D | 12 xi_phi = m_H^2/H^2 | lambda_phi(H) |
|---:|---:|---:|---:|---:|---:|
| 1e2 | 0 | 0 | 0 | 0 | 0.01000 |
| 1e3 | -2.440e-4 | -6.82e-5 | -1.51e-5 | -2.93e-3 | 0.01003 |
| 1e4 | -4.889e-4 | -1.366e-4 | -3.01e-5 | -5.87e-3 | 0.01006 |
| 1e6 | -9.815e-4 | -2.742e-4 | -6.03e-5 | -1.178e-2 | 0.01013 |
| 1e8 | -1.478e-3 | -4.128e-4 | -9.05e-5 | -1.774e-2 | 0.01020 |
| 1e10 | -1.978e-3 | -5.524e-4 | -1.208e-4 | -2.374e-2 | 0.01026 |
| 1e13 | -2.736e-3 | -7.637e-4 | -1.661e-4 | -3.284e-2 | 0.01036 |

**Sign.** For mu > mu_0, xi_phi < 0, so the curvature term is tachyonic for phi. During inflation the minima
move out to phi_m^2 = v^2 + 12|xi_phi| H^2/lambda_phi (phi_m ~= 1.78 H at 1e13 GeV on the X = 0 slice). The
barrier at phi = 0 grows by the same amount. In the threshold region H ~ 140-370 GeV, |12 xi_phi| <= 1.7e-3. That
raises the exponent B by 2 x 12|xi_phi| H^2/(lambda v^2) = 0.2% (140 GeV) to 4.5% (370 GeV) [C, D]; it is included
in every number below. xi_X < 0 makes X slightly tachyonic at phi = 0 above
H ~ 101 GeV (|m^2| <= 0.0092 H^2, stabilized by lambda_X). xi_D < 0 makes D tachyonic at phi = 0 above
H ~ 5.6 TeV. Z2_D is unbroken in flat space, so any D sign domains disappear when R -> 0 after inflation, and
there are no D walls [D].

**Scope.** The SM does not enter: every portal vanishes at mu_0 and stays zero at one loop. Graviton-loop
contributions (~(H/M_P)^2 <= 2e-11) are not included. The non-logarithmic one-loop de Sitter terms
(~g^2 H^2/(16 pi^2) without the large log) are O(1/ln(H/mu_0)) relative to the RG term [E], and were not
computed. The task's H range starts at H = mu_0 = 100 GeV, so no running below mu_0 is needed.

## 2. Stochastic spectator dynamics

### 2.0 Set-up [D]

The phi Langevin equation (Starobinsky-Yokoyama) for the coarse-grained field, per e-fold N, units H = 1,
y = phi/H, is dy/dN = -w'(y)/3 + eta/(2 pi). Here w = V/H^4 and V = U(phi, X=0, D=0; couplings(mu=H)) +
6 xi_phi(H) H^2 phi^2. Fokker-Planck:

    dP/dN = d/dy [w'/3 P] + D d^2P/dy^2,   D = 1/(8 pi^2),   P_eq ~ exp(-(8 pi^2/3) w).

**Discretization.** A square-root-approximation (SqRA) Markov chain on a cell-centered grid, with
k_{i->j} = (D/h^2) exp(-(s_j - s_i)/2) and s = (8 pi^2/3) w. It obeys detailed balance with exp(-s) exactly.

**Odd-sector identity [D].** w is even, so for a start z > 0, E_z[sgn phi_N] equals the survival probability
S(z,N) of the half chain (y > 0) with killing rate 2D/h^2 in the cell next to y = 0. Then p_-(z,N) = (1 - S)/2,
and the smallest phi-odd eigenvalue Lambda_1 is the ground state of that killed chain.

**Exact large-barrier tools [D].**
(i) Lambda_1 = 1/E_qsd[tau], computed from the birth-death mean absorption time
T_m = sum_{j<=m} (sum_{l>=j} pi_l)/(pi_j r_j) in log space. There is no cancellation, and it works at any
barrier.
(ii) A delta start in cell m0 uses reversibility, p_j(t)/pi_j = [e^{At} e_m0]_j/pi_m0. Then
p_- = (kappa/2)(pi_0/pi_m0) int_0^N [e^{At} e_m0]_0 dt. The integrand is the O(1e-2) probability of going
downhill from the barrier cell to m0, so p_- comes with relative precision even at 1e-286.

**Validation [C]** (`results/checks_1d.json`, `stochastic_1d.json` → `validation_pure_quartic`, `grid_convergence`).
- Pure quartic, lambda = 0.01: Lambda_1/sqrt(lambda) = 0.088921 (Starobinsky-Yokoyama value 0.08892, recalled,
  not read) and 2 Lambda_2/sqrt(lambda) = 0.57870 (MRST 1904.11917 abstract: 0.579).
- Spectral sum vs delta-start identity: 1.0329e-5 vs 1.0329e-5 at H = 300 GeV; 2.39999e-9 vs 2.40022e-9 at
  250 GeV.
- Kramers asymptotics vs MFPT Lambda_1: agree to 0.05% at B = 658 and 5% at B = 8.5.
- Grid halving changes Lambda_1 by <= 7.5e-5 (relative) and log10 p_- by <= 0.003.

**Validity of the overdamped stochastic picture.** It needs m << H. m_phi/H = 0.38 at 370 GeV, 0.56 at 250 GeV
and 1.0 at 141 GeV. So the lowest thresholds sit where the slow-roll Langevin approach is **marginal**. The
exponent itself is the Hawking-Moss exponent B = 8 pi^2 DeltaV/(3 H^4), and it does not rely on the overdamped
approximation. H_* scales as B^{1/4}, so a factor-10 prefactor error moves H_* by less than 1% at B ~ 150 [D].
The competing wall-bubble (Basu-Guth-Vilenkin) instanton has B_wall = 2 pi^2 sigma/H^3 [D]: a tree wall on the
equatorial S^3 of Euclidean de Sitter, area 2 pi^2/H^3, with gravitational back-reaction sigma/(M_P^2 H) < 1e-30
neglected. B_wall exceeds B_HM by a factor 4H/m_phi > 1 throughout the threshold region (e.g. 660 vs 166 at
141 GeV). Hawking-Moss dominates, and it exists because |V''(0)| = lambda v^2 < 4 H^2 (H > 50 GeV) [D].

### 2.1 (a) Equilibrium for H_inf >> v

**Symmetry [D, C].** P_eq(phi) = P_eq(-phi) exactly, because U and the xi term are even and X enters through
phi^2. The code solves the odd sector explicitly, and the equilibrium mean p_- is 1/2 identically.

**X = 0 slice [C].** At 1e13 GeV the distribution is bimodal: peaks at y = +-1.77, P_peak/P(0) = 1.98,
<phi^2> = 2.71 H^2, Lambda_1 = 3.70e-3 per e-fold. At 1e4 GeV: peaks at +-0.77, ratio 1.02, <phi^2> = 1.49 H^2,
Lambda_1 = 9.48e-3.

**With X fluctuations (2D phi-|X| FP) [C].** X is light wherever g_X|phi| < H, and its stochastic variance gives
phi a positive mass. At 1e13 GeV, <rho^2> = <2|X|^2> = 8.0 H^2, so 2 g_X^2 <|X|^2> ~= 0.08 H^2. That beats the
tachyonic 12 xi_phi H^2 = -0.033 H^2, and the phi marginal becomes **unimodal**: peak at 0, <phi^2> = 0.97 H^2
(0.73 H^2 at 1e4 GeV), Lambda_1(2D) = 8.45e-3 (1.37e-2 at 1e4 GeV). The equilibrium remains exactly symmetric.

### 2.2 (b) p_- after the last N_* e-folds

There are three start conditions:
- **fresh start**: phi = +v exactly, with X = 0 in 2D. The patch was in the +v vacuum, unfluctuated, when the
  last N_* e-folds began.
- **relaxed start**: P_eq restricted to phi > 0. This represents a long earlier inflation conditioned on the
  + sign.
- **equilibrium draws**: phi_* drawn from P_eq. The ensemble mean of p_- is 1/2, and the table reports the
  distribution of the patch's minority fraction.

Selected rows (`results/stochastic_1d.json`, 71 H values; X = D = 0 path):

| H (GeV) | B_HM | Lambda_1 | p_- fresh N50 | p_- fresh N60 | p_- relaxed N50 | eq-draw minority quantiles 10/50/90 (N50) |
|---:|---:|---:|---:|---:|---:|---|
| 100 | 658 | 2.6e-287 | 4.4e-286 | 5.7e-286 | 6.6e-286 | - |
| 141 | 165.5 | 9.7e-74 | 1.1e-72 | 1.6e-72 | 2.4e-72 | - |
| 200 | 41.8 | 2.7e-20 | 1.3e-19 | 2.5e-19 | 6.7e-19 | - |
| 251 | 16.8 | 1.23e-9 | 3.27e-9 | 7.15e-9 | 3.07e-8 | - |
| 298 | 8.48 | 3.35e-6 | 8.8e-6 | 1.79e-5 | 1.16e-4 | 3.5e-6 / 9.5e-6 / 5.6e-5 |
| 335 | 5.40 | 5.6e-5 | 1.79e-4 | 3.35e-4 | 2.1e-3 | 7.1e-5 / 2.1e-4 / 1.5e-3 |
| 398 | 2.77 | 5.6e-4 | 2.98e-3 | 4.87e-3 | 2.22e-2 | 1.1e-3 / 3.8e-3 / 0.043 |
| 501 | 1.15 | 2.3e-3 | 0.0224 | 0.0316 | 0.083 | 0.008 / 0.032 / 0.26 |
| 1000 | 0.110 | 7.6e-3 | 0.186 | 0.210 | 0.199 | 0.050 / 0.161 / 0.415 |
| 2037 | 0.026 | 9.4e-3 | 0.333 | 0.349 | 0.225 | 0.068 / 0.199 / 0.426 |
| 1e4 | 0.023 | 9.5e-3 | 0.461 | 0.465 | 0.227 | 0.070 / 0.200 / 0.433 |
| 1e6 | 0.090 | 8.0e-3 | 0.4996 | 0.4997 | 0.204 | 0.054 / 0.171 / 0.416 |
| 1e13 | 0.684 | 3.7e-3 | 0.5000 | 0.5000 | 0.121 | 0.017 / 0.066 / 0.336 |

**Fresh-start lag [C, D].** From a delta at +v, crossings are suppressed for the first 20-45 e-folds. In the
barrier regime, p_-(N) ~= (Lambda_1/2)(N - N_lag), with N_lag ~= 45 at 250-300 GeV and ~22 at 120 GeV. The most
probable uphill path is the time reverse of the roll down from the barrier top. That roll takes
~(3H^2/(lambda v^2)) ln(...) e-folds, because the top curvature |V''(0)|/(3H^2) = 0.053 at 250 GeV is small.
So a fresh start lowers p_- by up to 10x relative to the relaxed start. (A first MFPT-difference estimate of
N_lag did not match. It neglected fast-mode terms; see section 7.)

**2D check of the fresh start [C].** At 330 GeV, p_-(N50) = 2.55e-4 with X fluctuating (dt-halving change 1e-5
relative), against 1.33e-4 with X = 0. The factor is 1.9; X has to build up from zero first.

**H >> v, fresh start.** v/H is far below the grid spacing there. p_- = 1/2 - (1/2 - p(h/2)) (v/H)/(h/2) by
oddness of S, which gives 1/2 - O(v/H). Equilibrium draws give mean 1/2 by symmetry. Each patch's minority
fraction is broad. With X = 0 the probability that the minority is below p_c is 0.72-0.88. With X fluctuating
(2D) the medians are 0.38-0.41 (N50-N60) and P(minority < p_c) = 0.24-0.37 (`stochastic_2d.json`).

### 2.3 (c) Barrier crossing, counts in the observable patch, and thresholds

For H << DeltaV^{1/4} = 224 GeV, Lambda_1 ~= 2 x Kramers rate = (sqrt2 lambda v^2/(3 pi H^2)) exp(-B), with
B = (8 pi^2/3)(lambda_phi/4) phi_m^4/H^4 ~= (2 pi^2 lambda_phi/3)(v/H)^4 [D]. Each crossing creates one minority
region of the Hubble size at its time. The expected number inside the observable patch is

    N_reg = int_0^{N_*} (dp_-/dn) e^{3n} dn   (~ (Lambda_1/2) e^{3 N_*}/3 for a relaxed start) [D].

Only crossings in the first n'_X e-folds after the patch's horizon exit re-enter after T_X. Here
n'_X = ln(a_X H_X/a_0 H_0) is 17.77 for T_X = 1 MeV and 12.69 for T_dom = 8.0 keV [C]. At T_dom a tree wall per
Hubble patch carries sigma/(M_P^2 H) = 1.

**Thresholds H_* (GeV) [C]** (`results/stochastic_1d.json` → `thresholds_GeV`; interpolation on the 71-point
grid):

| quantity < target | fresh N50 | fresh N60 | relaxed N50 | relaxed N60 | relaxed + X,D entropic N50 |
|---|---:|---:|---:|---:|---:|
| p_- < 1e-3 | 368.5 | 355.0 | 323.9 | 321.9 | 295.4 |
| p_- < 1e-10 | 239.1 | 236.6 | 233.1 | 232.7 | 224.7 |
| N_reg (all sizes) < 1e-3 | 144.4 | 138.1 | 144.4 | 138.1 | 143.4 |
| N_reg (all sizes) < 1 | 146.1 | 139.3 | 146.1 | 139.3 | 145.1 |
| N(re-enter after 1 MeV) < 1e-3 | 196.5 | 196.5 | 186.2 | 186.2 | 183.2 |
| N(re-enter after T_dom) < 1e-3 | 235.3 | 235.3 | 202.1 | 202.1 | - |
| p_- = p_c (end of inflation) | 1784 | 1607 | not reached | not reached | - |

The relaxed-start p_- never reaches p_c in 1D; it peaks at 0.23 at H ~ 3-5 TeV.

**X and D corrections [C].** Multi-field Langer theory gives the rate x R_X R_D, with R = Z_perp(saddle)/Z_perp
(well) [D]. R_D = 1.005 throughout. R_X is 13.6 at 335 GeV, 23.7 at 251 GeV, 73 at 141 GeV and 144 at 100 GeV.
The 2D FP confirms the trend. Lambda_1(2D)/Lambda_1(1D) = 1.54, 2.31, 3.48, 6.40, 11.7 and 16.7 at
B = 0.11, 1.16, 2.71, 5.73, 10.9 and 17.1. As a fraction of R_X that is 0.62, 0.35, 0.36, 0.46, 0.61 and 0.70,
approaching R_X as B grows. Grid refinement changes the 330 GeV value by 4e-5. The entropic column above is
therefore a mild over-correction. The exponent is unchanged because the only saddle of V(phi, |X|) is (0, 0) [D].

**Consistency of N_\* [E].** With instant reheating into the visible SM bath, N_* for today's Hubble scale is
48.1, 48.5 and 49.1 at H_inf = 150, 300 and 1000 GeV, and 60.6 at 1e13 GeV. N_* = 60 together with
H_inf ~ 100s GeV therefore needs extra post-inflation expansion. N_* and the reheating history are preparation
data.

**Requirement for a +v start [D, E].** If the hidden field sat near phi = 0 when inflation began (for example a
hot hidden sector before inflation), the roll to +-v takes N_roll ~ (3H^2/(lambda v^2)) ln(2 pi v/H) e-folds.
That is ~41 e-folds at H = 200 GeV and ~82 at 300 GeV [D]. The fresh or relaxed start then needs N_tot - N_* >~ N_roll + N_lag.
Otherwise the roll happens inside the observable e-folds and creates domains at those scales. The
time-dependent roll was not simulated.

**Curve.** `results/p_minus_vs_Hinf.svg` shows log p_- against log H_inf for the fresh start (N_* = 50, 60), the
relaxed start (N_* = 50), the reference levels 1e-3, 1e-10 and p_c, and post-trapping points.

## 3. Classification

For a fresh start at +v, N_* = 50 (N_* = 60 in parentheses):

| regime | H_inf range | what forms |
|---|---|---|
| A | < 144 GeV (138) | no opposite-vacuum region expected anywhere in the observable patch (N_reg < 1e-3) |
| B | 144-368 GeV (138-355) | p_- < 1e-3: rare closed minority bubbles, isolated in a +v sea; they collapse after re-entry |
| C | 368-1784 GeV (355-1607) | 1e-3 < p_- < p_c: closed minority domains of all sizes (site-percolation proxy) |
| D | > 1784 GeV (1607); 1.43 TeV (1.35) after trapping | both vacua percolate: a macroscopic wall network forms at re-entry |

p_c = 0.3116077 (simple-cubic site percolation, Wang et al. 1302.0421, abstract) [Q]. The inflationary sign
field is correlated, with a power-law two-point function of tilt 2 Lambda_1 (Starobinsky-Yokoyama; MRST
spectral expansion). The uncorrelated-site p_c is therefore only a proxy [E].

**Minority-domain size distribution (small p_-) [C, D].** Stationary crossing makes the opposite-sign volume
fraction per unit ln k flat: dp_-/dn = Lambda_1/2. The number of minority domains of comoving size ~k^-1 per
Hubble volume at their re-entry is likewise ~Lambda_1/2 per e-fold of scale. A fresh start suppresses the largest
scales, those that exited in the first ~20-45 e-folds. With instant reheating and RD at constant g_*,
T_re = T_reh e^{-n} for a scale that exited n e-folds before the end. Example, H_inf = 250 GeV, N_* = 50
(`classification.json`):

| e-folds before end | T_re | dp/dn fresh | dp/dn relaxed | sigma/(M_P^2 H(T_re)) |
|---:|---:|---:|---:|---:|
| 0 | 1.3e10 GeV | 2.3e-10 | 4.5e-10 | 6e-32 |
| 10 | 6.1e5 GeV | 1.1e-10 | 4.5e-10 | 3e-23 |
| 20 | 28 GeV | 2.2e-11 | 4.5e-10 | 1.6e-14 |
| 30 | 1.8 MeV | 3.9e-13 | 4.5e-10 | 1.1e-5 |
| 40 | 106 eV | 3.4e-19 | 4.5e-10 | 5.7e3 |

Domains that re-enter below T_dom = 8.0 keV have walls that dominate their own Hubble patch before collapse. They
become candidate black-hole/strong-inhomogeneity sources (Garriga-Vilenkin-Zhang 1512.01819, abstract) [E].
Domains that re-enter early collapse into hidden scalar radiation and carry a tiny energy fraction. This is why
the "re-enter after 1 MeV" and "after T_dom" thresholds are the physically relevant ones. The all-sizes count is
dominated by bubbles of size ~1/H_inf that collapse at T ~ T_reh.

**Gaussian-field approximation (as requested) [C].** Linearize about +v. The smoothed-field variance is then
sigma^2(k) = (3H^4/(8 pi^2 m^2))(e^{-nu n_k} - e^{-nu N_*}) with nu = 2 m^2/(3H^2), and the minority fraction is
erfc(v/(sqrt2 sigma))/2. The exponent is exactly 4x the Hawking-Moss exponent, because the quadratic well
overestimates the barrier by a factor 4. Trapping in the other well is missed, and all opposite-sign volume sits
on scales that exited in the last ~3H^2/(2m^2) e-folds. At 250 GeV it gives 7e-32 against 2.4e-9 from the FP,
and at 300 GeV 1.6e-16 against 1.0e-5. **The Gaussian approximation fails for small p_- and is not used for any
threshold.**

**Literature on the fate of the walls [Q, abstracts].** If one vacuum is favored by initial conditions it
eventually dominates, and the network disappears after a scaling stage (Larsson-Sarkar-White hep-ph/9608319).
Walls in biased non-equilibrium transitions decay exponentially, unlike equilibrium ones (Coulson-Lalak-Ovrut
hep-ph/9508226; the task's identifier hep-th/9508226 returns 404). Wall survival depends on bias, tension and
volume pressure (Gelmini-Gleiser-Kolb, PRD 39 (1989) 1558). Here a fresh start in regimes B/C is an
"initial-condition bias" with no energetic bias. Whether the minority domains of regime C shrink away, and on what
time scale, is a network-evolution question that this module does not compute.

## 4. Post-inflation evolution when |phi_end| >> v (estimates)

### 4.1 Zero crossings [C, D]

Assume an outside-action instant radiation era. With a = tau and chi = a phi, the equation is
chi'' + lambda chi^3 - lambda v^2 tau^2 chi = 0, which is exact in RD because a'' = 0. phi is frozen until
H ~ sqrt(lambda) phi. The ODE gives N_cross = (0.300-0.320) phi_end/v (mean 0.313) for phi_end/v >= 30. The
derived bounds are 0.381 (trapping when the amplitude reaches v) and 0.270 (trapping at amplitude sqrt2 v, where
E = lambda v^4/4). Below phi_end ~= 2.6-3.0 v there is no crossing at all, because Hubble friction wins.

### 4.2 Sign sensitivity [C, D]

N_cross is proportional to phi_end, so dN/dy_end = 0.313 H/v. Two Hubble patches one e-fold apart differ by
delta y ~ sqrt2/(2 pi) = 0.225 in isocurvature. The expected number of sign flips between them is 1.1, 2.3, 6.8
and 21 at H/v = 10, 30, 100 and 300 [C]. **Signs decorrelate between neighboring patches for H_inf >~ 14 v**
(threshold 14.2). Purely adiabatic zeta only shifts the local clock. On uniform-density slices it does not
change the crossing count, so only the isocurvature delta phi ~ H/(2 pi) matters [D].

### 4.3 Combined with the FP end distribution [C]

`final_sign.py` maps every phi_end of the 1D fresh-start distribution through the RD trapping ODE:

| H (GeV) | 500 | 700 | 1000 | 1400 | 2000 | 3000 | 5000 | 1e4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| p_- end of inflation, N50 | 0.0222 | 0.0902 | 0.186 | 0.265 | 0.333 | 0.386 | 0.427 | 0.461 |
| p_- after trapping, N50 | 0.0222 | 0.0902 | 0.188 | 0.304 | 0.445 | 0.514 | 0.4999 | 0.4997 |

H_perc after trapping is 1.43 TeV (N50) and 1.35 TeV (N60). Below ~500 GeV no end value exceeds 2.4 v with
appreciable probability (1.7e-12 at 500 GeV), so the thresholds of section 2.3 are unaffected.

### 4.4 Preheating and nonthermal restoration [C for Floquet, E otherwise]

Lame chart for x'' + (kappa^2 + q cn^2(T, 1/2)) x = 0 with Floquet period 2K(1/2) = 3.7081 (one zero crossing):
- q = g_X^2/lambda_phi = 1 (X): band 0 < kappa^2 < 0.499, mu_max = 0.1470 at kappa^2 = 0.228.
- q = 3 (phi self-resonance): band 1.500 < kappa^2 < 1.731 (= sqrt3), mu_max = 0.03598 at kappa^2 = 1.614.

GKLS hep-ph/9705347 (abstract) states that analytic solutions exist for g^2/lambda = 1 and 3. The numbers above
agree with values recalled from that literature, but the band data were not read from the paper.

Back-reaction from conformal-vacuum seeds (an upper estimate; superhorizon X modes are much larger after
inflation) is reached after about 7-9 zero crossings for X, at eps = 0.1-1 of lambda phi^2, and after 30-39 for
phi self-resonance. With N_cross = 0.313 phi_end/v, X preheating back-reacts for phi_end >~ 22-29 v. For a typical
phi_end ~ 1.5 H_inf that means H_inf >~ 15-20 v.

**Restoration criterion (Hartree) [D].** m_eff^2(0) = -lambda v^2 + 3 lambda <dphi^2> + 2 g_X^2 <|X|^2> > 0, that
is <dphi^2> > v^2/3 or g_X^2 <|X|^2> > lambda_phi v^2/2. After back-reaction both variances and phi_amp^2 redshift
as a^-2. Write eps_fl = (3 lambda <dphi^2> + 2 g^2 <|X|^2>)/(lambda phi_amp^2). Restoration then lasts until
phi_amp^2 ~ v^2/eps_fl, and the final breaking starts from a fluctuation-dominated state with uncorrelated signs
[E]. Supporting literature: Kofman-Linde-Starobinsky hep-th/9510119; Tkachev hep-th/9510146; and
Khlebnikov-Kofman-Linde-Tkachev hep-ph/9804425, whose abstract reports temporary restoration of phi -> -phi by X
fluctuations and a first-order breaking for g^2 >> lambda. All were read at abstract level. Here g^2 = lambda
(q = 1), and no lattice run was made.

### 4.5 When the inflationary selection survives [E]

- **Hidden thermalization.** If the condensate energy thermalizes into the 4 hidden degrees of freedom,
  T_h = 0.209 phi_amp. That exceeds T_c = 1.557-1.564 v (thermal_potential module [Q]) while phi_amp > 7.46 v. In
  that case a thermal transition, not inflation, sets the domains.
- **Gravitational particle production.** The dimensional bound rho <~ H_inf^4 gives T_eq <= 0.93 H_inf. That is
  below T_c for H_inf < 1.67 v, which includes every threshold above.
- **Curvature oscillations after inflation.** Particle production from oscillating R needs xi >~ 1 (Herranen et
  al. 1506.04065, abstract). Here |xi| <= 2.7e-3, so it is negligible.
- **Visible reheating.** Even at H_inf ~ 300 GeV the visible T_reh ~ 1.5e10 GeV, but with zero portals that heat
  does not reach the hidden sector at this order.

## 5. Free inputs and conditional constraints

**Free inputs (all outside the action unless marked):**
- H_inf (constant during the last N_* e-folds)
- N_* = 50 or 60
- the start condition of the observable patch (fresh +v, relaxed + well, or equilibrium draw), and
  N_tot - N_* relative to N_roll + N_lag
- the post-inflation history (instant RD assumed; non-instant reheating changes the e-fold-to-temperature map and
  the crossing count)
- whether the inflaton or anything else heats the hidden sector above T_c
- RG scale choice mu = H (band H/2 to sqrt12 H). This is inside the EFT prescription but is a choice.
- the benchmark parameters (v = 1 TeV etc., from the action as specified)
- the 3D site-percolation p_c used as a proxy

**Conditional history constraints (not predictions):**
1. Within the class "inflation precedes the hidden history, hidden sector never re-heated above T_c, no
   nonthermal restoration, observable patch in +v at the start of the last N_* e-folds", the expected
   opposite-vacuum volume fraction is < 1e-3 iff H_inf <~ 295-368 GeV and < 1e-10 iff H_inf <~ 225-239 GeV. The
   expected number of opposite-vacuum regions of any size is < 1e-3 iff H_inf <~ 137-144 GeV. No opposite domain
   re-enters after 1 MeV (expected number < 1e-3) iff H_inf <~ 183-197 GeV. Ranges cover N_* = 50/60, the start
   conditions and the X correction.
2. For H_inf >~ 1.4-1.8 TeV the inflationary preparation **does not select a vacuum**: both vacua percolate. For
   H_inf >~ 14 v, post-inflation oscillation alone randomizes the sign between neighboring Hubble patches, and
   for H_inf >~ 15-20 v preheating plausibly restores the symmetry. In this class the wall network is generic at
   formation. Its persistence and cosmological consequence are outside this module.
3. None of these numbers is a prediction of H_inf or of a relic abundance. GRUT remains empirically unvalidated,
   and the physical scale v and the primordial charge yield still need independent justification.

## 6. Literature (read level stated)

| Reference | Read level | Used for |
|---|---|---|
| Starobinsky & Yokoyama, astro-ph/9407016 (PRD 50, 6357) | abstract | stochastic FP formalism; equilibrium and correlation radius. Lambda_1 = 0.08892 sqrt(lambda) is recalled, then reproduced numerically |
| Herranen, Markkanen, Nurmi, Rajantie, arXiv:1407.3141 (PRL 113, 211102) | sections (pdftotext of v2: eqs. 1-7, 10, 12) | beta_xi normalization (12 lambda (xi - 1/6)); M_i^2 = ... + (xi - 1/6) R; scale choice mu^2 = phi^2 + R |
| Herranen et al., arXiv:1506.04065 | abstract | curvature-oscillation production needs xi >~ 1 |
| Markkanen, Rajantie, Stopyra, Tenkanen, arXiv:1904.11917 | abstract | check 2 Lambda_2 = 0.579 sqrt(lambda); power-law sign correlations |
| Wang, Zhou, Zhang, Garoni, Deng, arXiv:1302.0421 | abstract | p_c(site, sc) = 0.3116077 |
| Larsson, Sarkar, White, hep-ph/9608319 (PRD 55, 5129) | abstract | initial-condition-favored vacuum eventually dominates |
| Coulson, Lalak, Ovrut, hep-ph/9508226 | abstract | biased walls decay exponentially (the task's hep-th/9508226 is a 404) |
| Gelmini, Gleiser, Kolb, PRD 39 (1989) 1558, doi:10.1103/PhysRevD.39.1558 | abstract (INSPIRE) | bias/tension regimes |
| Greene, Kofman, Linde, Starobinsky, hep-ph/9705347 | abstract | Lame structure; q = 1, 3 special cases (band numbers computed here) |
| Kofman, Linde, Starobinsky, hep-th/9510119 | abstract | nonthermal phase transitions after preheating |
| Tkachev, hep-th/9510146 | abstract | symmetry restoration at preheating |
| Khlebnikov, Kofman, Linde, Tkachev, hep-ph/9804425 | abstract | temporary restoration via X fluctuations; first-order breaking |
| Basu, Guth, Vilenkin, PRD 44 (1991) 340, doi:10.1103/PhysRevD.44.340 | abstract (INSPIRE) | wall bubbles nucleated during inflation. The exponent 2 pi^2 sigma/H^3 is derived here, not quoted |
| Garriga, Vilenkin, Zhang, arXiv:1512.01819 | abstract | inflation-nucleated walls can form black holes after re-entry |
| Hawking & Moss, PLB 110 (1982) 35 | not read | standard HM exponent (identical to the stochastic one, derived here) |
| Buchbinder, Odintsov, Shapiro (book) | not read | beta_xi derived here instead |

Downloads are kept in the scratch directory only. The 1407.3141 PDF was text-extracted there.

## 7. What failed, what was corrected, what was not done

**Failed or corrected during the work:**
- The first RG-identity check failed by an exact factor 2. The check, not the betas, had the wrong
  normalization. It was fixed and passes to 3e-16.
- The first scan used the two-state formula (1 - e^{-Lambda_1 N})/2 for the relaxed start even where no barrier
  exists (B < 1), where it is meaningless. It was replaced by the exact spectral average over P_eq restricted to
  phi > 0. Region counts are now set to null where p_- > 0.05.
- The fresh start at H >> v initially gave p_- = 0.4956, a grid artifact: the start cell sat at h/2 instead of
  v/H. It was replaced by the odd-symmetry interpolation.
- The tridiagonal eigenvalue loses Lambda_1 below ~1e-13 (it returned 1e-14 and negative noise at 200 and
  150 GeV). Every small-Lambda_1 number therefore uses the MFPT identity, which matches the eigenvalue to 1e-5
  where both work.
- An MFPT-difference estimate of the fresh-start lag disagreed with the direct calculation: 77 vs 45 e-folds at
  300 GeV, 12 vs 22 at 120 GeV. It neglects fast-mode terms and is not used.
- The 2D BDF and scipy expm_multiply runs were too slow (expm_multiply also emitted a NaN warning in
  start/stop mode). Both were replaced by L-stable BDF2, checked against expm_multiply (2.55361e-4 vs
  2.55358e-4) and by dt halving (6e-7 absolute in S).
- The first sign-sensitivity window was too coarse (zero flips at H/v = 10). It was replaced by a numerical
  dN/dy_end.
- The analytic crossing coefficient 0.381 overestimates the computed 0.313 by 22%, because trapping happens at
  amplitude ~sqrt2 v. Both are reported.

**Not done:**
- full 3-field (phi, X, D) Fokker-Planck; D enters only through the entropic factor R_D = 1.005
- the time-dependent roll from phi ~ 0 at the start of inflation
- non-overdamped (m ~ H) stochastic dynamics or a QFT de Sitter tunneling prefactor near H ~ 140 GeV
- non-logarithmic one-loop de Sitter terms and graviton loops in xi
- non-instant reheating
- lattice preheating and restoration
- percolation of the correlated sign field
- network evolution, persistence, scaling, backreaction, and any observable or CMB comparison
- black-hole formation from late-re-entering domains (only flagged)

## 8. Files

- `code/common.py`: parameters, one-loop running of the 12 parameters including the xi's, spectator potential
- `code/xi_running.py` → `results/xi_running.json`
- `code/fp1d.py`: SqRA chain, odd-sector identities, MFPT, delta-start relative-precision solver
- `code/checks_1d.py` → `results/checks_1d.json` (spectral vs delta identity, eigenvalue failure, N_lag)
- `code/stochastic_scan.py` → `results/stochastic_1d.json` (all rows, thresholds, validation, convergence)
- `code/fp2d_checks.py` → `results/stochastic_2d.json`
- `code/post_inflation.py` → `results/post_inflation.json`
- `code/final_sign.py` → `results/final_sign.json`
- `code/cosmo.py`: re-entry and N_* bookkeeping (estimates)
- `code/classification.py` → `results/classification.json`
- `code/plot_svg.py` → `results/p_minus_vs_Hinf.svg`
- `code/run_all.py` → `results/runtimes.json`
