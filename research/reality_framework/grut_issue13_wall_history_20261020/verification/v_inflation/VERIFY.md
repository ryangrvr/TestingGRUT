# Verification v_inflation: module inflation_spectator

21 October 2026. Adversarial check of `modules/inflation_spectator` (NOTES.md and results). All numbers below come
from code in this directory. I wrote that code from scratch before I opened the module code. I read the module code
afterwards, only to locate the discrepancies listed here. Inflation, H_inf, N_* and the start condition are
outside-action preparation data, as in the module. Python 3, numpy and scipy only, 2 threads; every job ran in under
2 minutes.

Labels: [V] computed here, [D] derived here, [Q] quoted (literature, with read level).

## Files

| file | content |
|---|---|
| `rg_xi.py` → `results_rg.json` | One-loop running of the full quartic tensor, the mass matrix and the xi matrix for 4 real fields (phi, x, y, D). Exact symbolic 4th derivatives. Tests: O(4) normalization against HMNR, and closed forms. |
| `fp1d_verify.py` → `results_fp1d.json`, `results_fp1d_aux.json` | 1D Starobinsky-Yokoyama FP. Backward-Kolmogorov solve for f = 1 - E[sgn phi_N] with f(0) = 1, central differences, implicit Euler, Richardson in dt and h. Lambda_1 by power iteration on the explicit Green function. Spectral cross-check and threshold root finding. |
| `hm_exponent.py` → `results_hm.json` | Lambda_1 against B_HM and the Kramers form, from B = 2.8 to 658. |
| `start_offset_check.py` → `results_start_offset.json` | Reproduces the module's mid-H p_- when the start is moved to the module's start cell. |
| `fp2d_verify.py` → `results_fp2d.json` | 2D (phi, rho = sqrt2 abs(X)) FP, finite-volume, with face conductances. Lambda_1(2D) and a fresh-start p_- at 330 GeV. Decoupling test (g = 0 reproduces 1D to 1e-5). |
| `entropic.py` → `results_entropic.json` | R_X and R_D (Langer transverse partition ratios), and diagnostics of their validity. |
| `postinf_lame.py` → `results_postinf_lame.json` | RD zero-crossing ODE in conformal variables (chi = a phi, a = eta), and the Lame/Floquet chart. |
| `final_sign_verify.py` → `results_final_sign.json` | Forward FP end distribution mapped through the trapping ODE: p_- after trapping and H_perc. |

Convention used throughout (the module uses the same one): V_R = (1/2) xi R z^2 for each real field, R = +12 H^2 in de
Sitter, conformal value xi = 1/6. GRUT_CANDIDATE_1.md §3 writes `xi_phi R phi^2` without the 1/2, so in the action-text
normalization every xi value below is halved. The module states this.

## Claim-by-claim

### 1. beta_xi and the running of xi_phi from xi(mu_0 = 100 GeV) = 0

**Method [D, V].** I use the general one-loop curved-space result 16 pi^2 beta_{xi_ab} = Lambda_abcd (xi_cd - delta_cd/6),
with Lambda_abcd = d^4V/dz^4 (Buchbinder-Odintsov-Shapiro, **not read**). Its normalization is fixed by reproducing
HMNR 1407.3141 eq. (7), 16 pi^2 beta_xi = (xi - 1/6)(12 lambda + ...), for V = lambda (H^dag H)^2. I read the text of
eqs. (1)-(7) and Table I from the arXiv PDF, which also gives the convention xi R phi^2/2 and theta = xi - 1/6. My O(4)
test returns 12 lambda (and 24 lambda^2 for beta_lambda).

**Hand derivation.** g_X^2 phi^2 abs(X)^2 = (g_X^2/2) phi^2 (x^2 + y^2) gives Lambda_{phi phi x x} = Lambda_{phi phi y y}
= 2 g_X^2. h_D^2 phi^2 D^2/2 gives Lambda_{phi phi D D} = 2 h_D^2, and Lambda_{phi phi phi phi} = 6 lambda_phi. Therefore

    16 pi^2 beta_xi_phi = 6 lambda_phi (xi_phi - 1/6) + 4 g_X^2 (xi_X - 1/6) + 2 h_D^2 (xi_D - 1/6).

This equals the module's formula. Both real components of X are counted, so no factor-2 error is present.

**Numbers [V].** B_xi_phi(mu_0) = -0.0167, so beta = -1.05754e-4 per e-fold. The full tensor running of all couplings,
evaluated at mu = H:

| H (GeV) | xi_phi | xi_X | xi_D | lambda_phi |
|---:|---:|---:|---:|---:|
| 1e3 | -2.4397e-4 | -6.817e-5 | -1.507e-5 | 0.010032 |
| 1e13 | **-2.7363e-3** | -7.637e-4 | -1.661e-4 | 0.010365 |

- Band at 1e13 GeV: mu = H/2 gives -2.660e-3, and mu = sqrt12 H gives -2.874e-3.
- Leading log: -2.679e-3.
- 12 xi_phi = -0.03284.
- phi_m/H = 1.780.
- beta_alpha(mu_0) = 1.26651e-8 (BRIEF value reproduced).

**Verdict: CONFIRMED.** The value, the sign (tachyonic for every mu > mu_0, because Lambda > 0 and xi - 1/6 < 0) and
the whole module table agree to every printed digit.

**Caveat (overclaim, minor).** NOTES §1 says the non-logarithmic one-loop de Sitter terms are O(1/ln(H/mu_0)) relative
to the RG term. In the threshold region ln(H/mu_0) is only 0.9-1.3, so these terms are O(1) relative to the
RG-generated xi there. The whole xi effect raises B by only 0.2-5.6% [V: B with and without xi at 141-398 GeV], which
moves H_* by about 1% or less. The caveat does not change the thresholds, but the justification as written is wrong in
that region.

### 2. Fresh-start p_- after N_* e-folds, and the thresholds

**Method [V].** I use V = (1/2) m^2(H) phi^2 + lambda(H) phi^4/4 + 6 xi(H) H^2 phi^2, with m^2(mu_0) = -lambda v^2 and
all couplings run to mu = H. The Langevin equation is dy/dN = -w'/3 + eta/(2 pi), so D = 1/(8 pi^2) and
P_eq ∝ exp(-8 pi^2 w/3). I solve the backward equation for f(y, N) = 2 p_-(start y). This is a different method from
the module's SqRA chain and reversibility identity.

Accuracy:
- The dt and h Richardson extrapolations change p_- by less than 1e-5 relative.
- A spectral check by eigen-expansion agrees: 5.2727e-9 against 5.2778e-9 at 250 GeV, and 7.82165e-4 against
  7.82170e-4 at 350 GeV.
- For the pure quartic, Lambda_1/sqrt(lambda) = 0.088924, against 0.08892 (SY; the value is in my memory, and SY was not
  read). 2 Lambda_2/sqrt(lambda) = 0.57876, against 0.579 (MRST 1904.11917; abstract not re-read).

| H (GeV) | p_- N50 | p_- N60 | module (N50/N60) |
|---:|---:|---:|---|
| 250 | 2.412e-9 | 5.278e-9 | (251.19: 3.27e-9 / 7.15e-9; mine at 251.19: 3.260e-9 / 7.132e-9) |
| 335 | 1.820e-4 | 3.402e-4 | 1.79e-4 / 3.35e-4 |
| 350 | 4.330e-4 | 7.822e-4 | - |
| 1000 | 0.18381 | 0.20755 | 0.186 / 0.210 |
| 2000 | 0.33011 | 0.34626 | - (2037: 0.333 / 0.349; mine 0.33304 / 0.34894) |

Thresholds [V] (root finding on the converged p_-):

| quantity | mine | module | verdict |
|---|---:|---:|---|
| p_- = 1e-3, N60 | 355.07 GeV | 355.0 | CONFIRMED |
| p_- = 1e-3, N50 | 367.8 | 368.5 | CONFIRMED (0.2%) |
| p_- = 1e-10, N60 | 236.6 | 236.6 | CONFIRMED |
| p_- = p_c = 0.3116077, N60 | **1616.6** | 1607 | CORRECTED (+0.6%) |
| p_- = p_c, N50 | **1792** | 1784 | CORRECTED (+0.45%) |

**Cause of the small mid-H offset (located in the module output after I had my numbers).** For H ≳ 375 GeV the module
uses a fixed grid with h = 0.02 and starts in the cell at y = 0.99, 0.73 or 0.59 instead of v/H = 1.0, 0.7305 or 0.5988.
When I move my start to the module's cell, I reproduce its p_- to 4 digits (0.209667 against 0.209661 at 1000 GeV N60;
0.319736 against 0.319731 at 1670 GeV). A start too close to phi = 0 biases p_- upward by up to about 1.2% at
1-1.7 TeV. This shifts H_perc, and the regime C/D boundary, down by 0.5-0.6%. At 330-335 GeV the same effect is about
2% in p_- (1D fresh p_- at 330 GeV: 1.309e-4 mine against 1.334e-4 module). It has a negligible effect on the 1e-3
threshold.

### 3. Hawking-Moss / stochastic exponent

**Method [V].** I computed Lambda_1 (odd sector, killed half-line) with a Green-function power iteration and
Collatz-Wielandt brackets, from B = 2.8 to 658.
- The local slope d ln Lambda_1/dB is -0.865, -0.917, -0.956, -0.974, -0.984, -0.991, -0.996, -0.998 and -0.999
  between neighbouring H from 398 GeV down to 100 GeV. It tends to -1.
- Lambda_1/Kramers (Kramers = sqrt2 a/(3 pi) e^{-B}) is 0.90 (B = 2.8), 0.95 (8.5), 0.976 (16.8), 0.991 (41) and 0.998
  (165, 317).
- At B = 658 my central-difference grid is not converged: Richardson gives 2.60e-287, against 2.6e-287 in the module,
  and a ratio of 0.984 to Kramers.
- The B_HM values agree with the module table: 658.0, 164.8 at 141.42 GeV (165.5 at the module's 141.25), and 41.4 at
  200 GeV (41.8 at 199.5).
- H_* ∝ B^{1/4} [D]: a factor 10 in the prefactor moves H_* by ln10/(4B), which is 0.4% at B = 150.

**Verdict: CONFIRMED.** The stochastic exponent equals B_HM = (8 pi^2/3) DeltaV/H^4, and the Kramers prefactor is
confirmed. "Agree to 0.05% at B = 658" is NOT CHECKED (my grid is not converged at that B).

**Lambda_1 column (NOTES table and §2.1): CORRECTED.** The module's Lambda1_mfpt is 1/E_pi[tau], with pi the
equilibrium restricted to phi > 0. NOTES calls it 1/E_qsd[tau]. It equals Lambda_1 only at large barrier. Correct
values from the eigenvalue, confirmed independently three ways (power iteration, dense tridiagonal eigensolver, and an
ODE shooting method with Y = 5, 6 and 7):

| H (GeV) | module | correct | overstatement |
|---:|---:|---:|---:|
| 501 | 2.3e-3 | **2.13e-3** | 8% |
| 1000 | 7.6e-3 | **6.32e-3** | 20% |
| 2037 | 9.4e-3 | **7.61e-3** | 24% |
| 1e4 | 9.48e-3 | **7.68e-3** | 23% |
| 1e6 | 8.0e-3 | **6.57e-3** | 22% |
| 1e13 | 3.70e-3 | **3.30e-3** | 12% |

The module's own `Lambda1_eig` column holds the correct values. Lambda_1 must lie below the pure-quartic value
0.0889 sqrt(lambda) = 8.9e-3 when the quadratic term is tachyonic, so 9.4-9.5e-3 is impossible.

Thresholds are not affected. The p_- values use the delta-start identity, and the large-B rates are accurate.

The 2D/1D ratios in NOTES §2.3 also use this denominator:
- B = 0.11: 1.54 → **1.86**.
- B = 1.16: 2.31 → **2.49**.
- B = 2.71: 3.48 → **3.55**.
- 6.40 (B = 5.73), 11.7 (10.9) and 16.7 (17.1) are CONFIRMED.

### 4. Post-inflation zero crossings in RD

**Method [V].** I integrate chi'' + chi^3 - eta^2 chi = 0 (chi = a phi/v, a = eta, which is exact in RD), with phi
frozen at phi_end at the start. The count stops when E < 0 (trapped). The results do not depend on eta_0 from 1e-5 to
1e-2.
- N_cross(phi_end/v) is 31, 47, 62, 94, 157 and 314 at 100, 150, 200, 300, 500 and 1000.
- A linear fit on [30, 300] gives a slope of **0.3149** (intercept -0.49).
- The mean ratio on a uniform 30-300 grid is 0.3107. The ratio swings from 0.289 to 0.315 because N is an integer.
- There is no crossing at phi_end ≤ 2.8 v, and one crossing at 3.0 v.
- The analytic bounds are 0.38138 and 0.26968 [D]: chi_A^2 = sqrt2 Phi_end and N = chi_A^2/(2 K(1/2) c), with c = 1
  or sqrt2.

**Verdict: CONFIRMED.** N_cross = 0.313 phi_end/v holds to within 1% (asymptotic slope 0.315). The bounds are
confirmed, and so is "no crossing below 2.6-3.0 v". The stated scatter "0.300-0.320" is slightly too narrow (0.289 at
phi_end ≈ 38 v).

**Sign-flip numbers (§4.2): CORRECTED.** The flip count is 0.315 x 0.225 x H/v, which gives 0.71, 2.1, 7.1 and 21 at
H/v = 10, 30, 100 and 300. The module's 1.1 at H/v = 10 comes from an integer finite difference over Delta y = 0.2.
The quoted threshold of 14.2 is right (mine 14.1).

**Post-trapping [V]** (forward FP end distribution mapped through the ODE with eta_0 = sqrt(sqrt(lambda) v/H)):
- p_- after trapping, N50: 0.1856, 0.3044, 0.4457 and 0.5135 at 1, 1.4, 2 and 3 TeV. The module has 0.188, 0.304,
  0.445 and 0.514.
- H_perc after trapping: **1424 GeV (N50) and 1348 GeV (N60)**. The module has 1.43 and 1.35 TeV.

Verdict: CONFIRMED.

### 5. Lame chart

**Method [V].** Monodromy of x'' + (kappa^2 + q cn^2(T, m = 1/2)) x = 0 over 2K(1/2) = 3.7081, with
mu = arccosh(abs(tr)/2)/(2K).

| q | band | mu_max | at kappa^2 |
|---:|---|---:|---:|
| 1 | 0 < kappa^2 < 0.500000 | 0.146993 | 0.22847 |
| 3 | 1.500000 < kappa^2 < 1.7320508 (= sqrt3) | 0.035979 | 1.6145 |

Spot values:
- q = 1: 0.1238, 0.1460, 0.1470, 0.1411 and 0.1117 at kappa^2 = 0.1, 0.2, 0.228, 0.3 and 0.4.
- q = 3: 0.02981, 0.03570, 0.03598, 0.03427 and 0.02459 at kappa^2 = 1.55, 1.6, 1.615, 1.65 and 1.7.

**Literature [Q, section level].** I read GKLS hep-ph/9705347 eqs. (18)-(19), (33), (34), (39), (40) and the table in
§VI-VII, from the arXiv PDF text in scratch. The paper gives:
- for g^2/lambda = 1: band 0 < kappa^2 < 1/2 and mu_max ≈ 0.1470 at kappa^2 ≈ 0.228;
- for g^2/lambda = 3: band 3/2 < kappa^2 < sqrt3 and mu_max ≈ 0.03598 at kappa^2 ≈ 1.615;
- a table that matches every spot value above to 4 digits.

GKLS eq. (34), transcribed, reproduces my q = 1 values to 1e-12. The GKLS model is lambda phi^4/4 + g^2 phi^2 chi^2/2
with real chi, and g_X^2 phi^2 abs(X)^2 = (g_X^2/2) phi^2 (x^2 + y^2), so q = g_X^2/lambda_phi = 1 with no factor 2.

**Verdict: CONFIRMED.** The module's band edge "0.499" is exactly 1/2 (GKLS eq. 33). Its literature row says
"abstract; band data not read". I confirm that its numbers are those of the paper. The abstract statement that analytic
solutions exist for g^2 = lambda and 3 lambda is also in the abstract [Q, read].

The back-reaction "7-9 crossings" [E] is CONFIRMED at the order-of-magnitude level. My own estimate is 8-10, from
n_k growth e^{2 mu_max 2K} ≈ 2.97 per crossing and g^2 <X^2>_vac/(lambda phi^2) ~ 2e-5.

## X entropic correction: is it justified?

**Arithmetic [V]: CONFIRMED.** R_X = Z_perp(phi = 0)/Z_perp(phi = phi_m) over both real X components, and R_D is the
same ratio over one component.
- R_X = 143.8, 72.9, 23.70 and 13.62 at 100, 141.42, 251.19 and 335 GeV. The module has 144, 73, 23.7 and 13.6.
- R_D = 1.0050 (Gaussian value sqrt(1 + h^2 v^2/m_D0^2) = 1.00499).

**Physical justification: weak. The factor is an over-correction, not a controlled correction.**
1. **The saddle is degenerate in X.** At phi = 0, X is massless at tree level (m_X0 = 0) and is held only by
   lambda_X = 0.001. Langer's formula assumes Gaussian transverse modes and asymptotically small noise, and neither
   holds here.
2. **X fluctuations at the saddle leave the region where phi = 0 is a saddle.** The saddle equilibrium has
   <rho^2> ≈ 7.1 H^2. The saddle condition fails for rho^2 > rho_c^2 = abs(m_phi,eff^2)/g^2, which is 9.1 H^2 at
   335 GeV and 16 H^2 at 251 GeV. P(rho^2 > rho_c^2) is 0.31 at 335 GeV, 0.21 at 300 GeV, 0.07 at 251 GeV and 1e-8 at
   141 GeV [V].
3. **X never reaches the saddle equilibrium during the window.** The X relaxation time at the saddle is
   ~1/(lambda_X <rho^2>/H^2) ≈ 140 e-folds, far longer than N_* = 50-60. The quasi-equilibrium that the factor assumes
   is not reached within the observable window. It is reached only for the asymptotic stationary rate.
4. **The 2D FP shows the overshoot.** My independent 2D FP reproduces the module's Lambda_1(2D) to 1e-4:
   1.1750e-2, 5.271e-3, 2.0457e-3, 2.6977e-4, 3.939e-6 and 1.5076e-8 at 1000, 500, 400, 330, 280 and 250 GeV. The true
   2D/1D ratio is only 0.75, 0.38, 0.36, 0.46, 0.61 and 0.70 of R_X. It approaches R_X only as B grows.
5. **The fresh start is far from R_X.** My 2D fresh-start p_- at 330 GeV is 2.501e-4 (N50) and 5.37e-4 (N60), against
   1.309e-4 and 2.476e-4 in 1D. The factor is 1.91 (N50) and 2.17 (N60), not R_X = 14.
6. **Effect on the threshold.** Rescaling the "relaxed + X,D entropic" p_- < 1e-3 threshold (295 GeV) by the 2D/R_X
   fraction (~0.4 at B ≈ 3-6) moves it to about 307 GeV [E: d ln p/d ln H ≈ 23 there].

The direction of the X effect is confirmed: it lowers the thresholds. The entropic column overstates it. Use it as a
bound, not as a central value. Below the band, at H ~ 140 GeV with B ~ 160, Langer is closer to valid.

## Other observations

- **Thin-wall wall bubble (BGV) comparison (§2.0): overclaim, minor.** B_wall = 2 pi^2 sigma/H^3 and
  B_wall/B_HM = 4H/m_phi are correct algebra [D]. The thin-wall instanton requires a wall thickness ~2/m_phi much less
  than 1/H. Here m_phi/H ≤ 1 in the threshold band, so the thin-wall formula is not valid there and the comparison
  "660 vs 166" is not meaningful. The conclusion survives, for the module's second reason: there is no CdL/wall
  instanton when abs(V''(0)) = lambda v^2 < 4H^2 (H > 50 GeV), so Hawking-Moss is the relevant exponent.
- **Typical phi_end (§4.4) [E].** "phi_end ~ 1.5 H_inf" uses the X = 0 slice. The module's own 2D marginal gives an
  rms phi of 0.85-0.98 H (<phi^2> = 0.73-0.97 H^2). The preheating back-reaction condition H_inf ≳ 15-20 v therefore
  becomes about 22-34 v. Order-of-magnitude statement.
- **Section 4.5 arithmetic [D]: CONFIRMED.**
  - T_h = (30 lambda/(4 pi^2 g))^{1/4} phi_amp with g = 4 gives 0.2088 phi_amp, so phi_amp > 7.46 v for T_h > 1.557 v.
  - T_eq ≤ 0.934 H_inf, so H_inf < 1.67 v.
  - Visible T_reh(H = 300 GeV) = 1.46e10 GeV (g_* = 106.75, reduced M_P).
- **Re-entry bookkeeping [D]: CONFIRMED.** n'_X(1 MeV) = 17.76 (module 17.77). T_dom = 8.0 keV for the stated
  criterion sigma/(M_P^2 H) = 1. That criterion is an O(1) convention: rho_wall = A sigma/t = 2A sigma H against
  3 M_P^2 H^2 gives T_dom × (2A/3)^{1/2}, about 0.7x for A = 0.8.
- **Gaussian-approximation exponent = 4 B_HM [D]: CONFIRMED.**
- **N_* ≈ 48.5 at H = 300 GeV and ≈ 61 at 1e13 GeV for instant reheating [D, rough]: CONFIRMED.** The use of
  N_* = 60 with H ~ 100s GeV therefore needs extra expansion, as the module says. It is preparation data.
- **NOT CHECKED:**
  - N_reg counts and their thresholds (138-146 GeV and 183-197 GeV);
  - the relaxed-start and equilibrium-draw columns;
  - the 2D high-H equilibrium quantiles;
  - Kramers agreement to 0.05% at B = 658;
  - literature statements read only at abstract level by the module: SY, MRST, Wang et al. p_c, LSW, CLO, GGK, KLS,
    Tkachev, KKLT, BGV, GVZ and HMNR 1506.04065. The values 0.08892, 0.579 and p_c = 0.3116077 agree with what I
    reproduce or know. They are not re-read.

## Overall

The module's headline numbers survive independent recomputation:
- xi_phi(1e13) = -2.736e-3;
- the 355/368 GeV thresholds;
- the 236.6 GeV 1e-10 threshold;
- the HM exponent;
- N_cross = 0.313 phi_end/v;
- the Lame chart;
- the post-trapping H_perc of 1.42/1.35 TeV.

Corrections:
- H_perc at the end of inflation is 1617 GeV (N60) and 1792 GeV (N50), against 1607 and 1784;
- the small-B Lambda_1 values are overstated by up to 24% (wrong estimator, mislabelled as QSD);
- three small-B 2D/1D ratios;
- the H/v = 10 flip count.

The X entropic correction is an over-estimate, by 1.3-2.8x in the rate at B ≤ 17 and about 7x for the fresh start. It
should be presented as a bound, not as one of the "ranges". None of these changes alters the module's qualitative
classification. The conditional, preparation-data framing is appropriate.
