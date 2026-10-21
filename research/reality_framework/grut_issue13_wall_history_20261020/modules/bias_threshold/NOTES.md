# Module bias_threshold: how much systematic population bias removes the Z2 wall network, and when

21 October 2026. Claude module of the issue #13 wall-history packet (BRIEF.md). The action is unchanged. Nothing here
validates, excludes or identifies Candidate 1. The results are **conditional thresholds on a preparation datum** (the
systematic bias of the vacuum population at formation), and they rest on an extrapolation of the phase-ordering
scaling law over 11-17 decades in L/xi0 beyond any run.

Labels: **[C]** computed by code here, **[D]** derived in closed form here, **[Q]** quoted (source given), **[P]**
preparation input outside the action, **[E]** estimate.

Reproduce (about 4 min wall on 2 single-threaded workers, 7.1 CPU-min; then 3 s):

    cd modules/bias_threshold && OMP_NUM_THREADS=1 python3 code/run_all.py

`run_all.py` runs `lattice.py` (31 deterministic jobs, the largest 16 s; existing `results/lattice_runs/*.json` are
reused), then `analyze_lattice.py`, `fit_K10.py` (reads network_sim's stored runs, read only) and `thresholds.py`.
Outputs: `results/lattice_response.json`, `results/fit_K10.json`, `results/thresholds.json`, `results/headline.json`,
`results/runtimes.json`, `results/lattice_runs/` (time series only, 0.7 MB). A re-run of one job reproduced its
series exactly (difference 0.0).

## 0. Result

Notation: m = 2p - 1 is the systematic magnetization (population bias) of the sign field. m_f is its value when the
domains have formed (the lattice formation time; physically the Kibble-Zurek freeze-out). K is the comoving growth
factor of the network length V/Area since formation.

1. **Law [D, C].** A small systematic bias grows as m(K) = erf(y_f K^gamma), with y_f = erfinv(m_f) and
   gamma = d - lambda_bar. Here lambda_bar is the Fisher-Huse autocorrelation exponent. The network falls below 10% of
   the unbiased wall area when y_f K^gamma = x10.
   - 3D relativistic calibration (network_sim's biased PRS runs, collapsed on the formation bias):
     **gamma = 1.41 +- 0.05, x10 = 1.54** [C].
   - The required coarsening is K10 = (x10/erfinv m_f)^(1/gamma).
2. **Relativistic walls follow the model-A exponent within the errors [C].** These numbers come from the 3D wave
   equation in the radiation era (PRS).
   - lambda_bar = 1.61 +- 0.02 (box-safe) and 1.69 +- 0.07 (to L/2), from the unbiased autocorrelation.
   - gamma = 1.41 +- 0.05 from the bias collapse.
   - Model A gives lambda_bar = 1.60-1.67 (gamma = 1.33-1.40) [Q]. Our own overdamped 96^3 control gives 1.64.
   - In 2D at equal dynamic range: relativistic gamma = 0.84-0.85, our model A 0.85. The model-A asymptote is 0.75
     [Q, Bray], so both runs sit above it at short range.
   - The fit to the 3D data in network_sim's convention (m0 = initial bias) is **K10 = c m0^-alpha with
     alpha = 0.72 +- 0.05, c = 0.49 (+0.045/-0.041)** [C]. Theory predicts 1/(3 - lambda_bar) = 0.71-0.75, and the
     percolation exponent nu_3D = 0.88 is disfavoured.
   - In 2D, alpha = 1.27 +- 0.05 against 1/(2 - 1.25) = 1.33.
3. **Thresholds for H-star** [P: xi = 5e-4; formation at T_SM = 3.1e6 GeV, KZ domains] [C, conditional].
   - Parity needs K_par = 10^16.7 (band 10^15.8-10^17.1).
   - The relic-safe decay temperature is T_safe = 31 MeV (f = 1). The band is 0.31-62 MeV over f in [0.01, 1] [P],
     A = 0.7-0.9 and the emission factor kappa = 1-2 [E]. It needs K_safe = 10^12.9 (band 10^11.7-10^15.3).
   - **eps_persist = 6e-24** (band 10^-25.3 to 10^-19.8).
   - **eps_safe = 1.5e-18** (band 10^-22.7 to 10^-14.9).
   - The fate classes by m_f:
     - m_f < eps_persist: the network survives to parity (rejected by the parity/ΔN_eff test).
     - eps_persist < m_f < eps_safe: the network decays after T_safe and overproduces stable h/X/D (rejected
       retrospectively by Omega h^2 <= 0.12).
     - m_f > eps_safe: the network decays early. That is harmless for the parity and relic tests (other signals, e.g.
       gravitational waves, not computed).
   - The exponent gamma dominates the band. At central inputs, gamma = 1.30 / 1.41 / 1.50 gives eps_persist =
     3e-22 / 6e-24 / 2e-25 and eps_safe = 3e-17 / 1.5e-18 / 9e-20.
   - H-star-Q and the inflationary-randomization histories give similar numbers (§4).
4. **Consequence for the report's classification [C, conditional].** Every residual bias the review lists is far
   above eps_safe, if it is systematic across the observable patch:
   - a thermal-remnant zero mode of 1e-3 of the soft-mode level gives m_f = 8e-4;
   - post-inflation trapping gives 2-6e-4;
   - forced charge at 2 n_c gives 0.42;
   - a relaxed start gives 0.54.
   - Each such network decays at K10 = 10^0.3-10^2.8, before it reaches horizon size (K = 10^4.7 for H-star). That
     is case (iii).
   - A hot history therefore ends with a persistent network **only if the restored phase erases the sign memory to
     better than about 1e-23 (band 1e-25 to 1e-20)**. For a 1e-3 remnant that means 46 amplitude e-folds of damping
     (band 39-51); 34 e-folds (27-45) suffice to reach the relic-overproduction class.
   - The draft's sorting "unbiased network: rejected" is therefore incomplete. Between "exactly unbiased" and
     "memory kept" lie three distinct outcomes, set by a residual bias that has not been computed.
5. **Statistical versus systematic [D].**
   - An exact-Z2 ensemble has m = 0 systematically. A region of size R has a random imbalance (xi0/R)^(3/2), and the
     memory of it decays relative to the domain-scale pattern because lambda_bar >= d/2.
   - So finite-volume fluctuations cannot remove the network. A systematic bias is the same on every scale and grows
     as K^gamma.
   - The statistical floor at parity, K_par^-3/2 = 1e-25 (H-star), lies 1.8 decades below eps_persist for
     gamma = 1.41.
6. **Extrapolation [D].**
   - The 3D lattice reaches K <= 2.7 (network_sim biased runs) and K <= 5.4 (here).
   - Parity needs 15.1-16.7 more decades and the relic-safe point 10.9-14.9. The inflationary cases need up to 21.
   - The thresholds are therefore theory-based extrapolations of the scaling form, calibrated on runs. They are not
     run results.

## 1. Theory [D, Q]

**Response identity.** Let the initial field be a Gaussian ensemble with a uniform shift. Then
S(k,t,0) = Delta G(k,t) holds exactly, with G(k,t) = <d phi_k(t)/d phi_k(0)> (Bray, Adv. Phys. 43 (1994) 357,
arXiv:cond-mat/9501089, §4.2 eq. (93)).

- The proof is an integration by parts over the Gaussian initial distribution, so it holds for any deterministic
  dynamics, including second-order relativistic dynamics.
- With dynamical scaling, C(r,t,0) = L^-lambda_bar f(r/L) and G(k,t) = L^lambda g_R(kL), so lambda = d - lambda_bar
  (Bray eq. (95); Bray's lambda is the response exponent, and his lambda_bar is the Fisher-Huse lambda).
- The k = 0 response therefore gives a mean order parameter that grows as m0 L^(d - lambda_bar) for small m0.
- In Bray's §5.2 ("External fields/initial bias", eq. (175)) the Gaussian theory gives the full time dependence
  <phi> = erf(m0(0) (8 pi t)^(d/4)/sqrt(2 Delta)). The exact statement replaces t^(d/4) = L^(d/2) by L^lambda
  (Bray-Kissner 1992, quoted there).
- We use m(K) = erf(y_f K^gamma), gamma = d - lambda_bar, normalized at formation.

**Bound [D].** Cauchy-Schwarz gives S(k,t,0) <= sqrt(S(k,t) S(k,0)). With S(k->0,t) ~ L^d g(0) and a short-range
initial S(k,0), this yields L^(d - lambda_bar) <= const L^(d/2). Hence lambda_bar >= d/2 and gamma <= d/2.

- This is the Fisher-Huse lower bound (PRB 38 (1988) 373; via Bray ref. [19], not read). It needs only scaling and
  short-range initial correlations, so it applies to relativistic networks.
- Gaussian/OJK theory, and Hindmarsh's relativistic Gaussian ansatz (hep-ph/9605332), saturate it: gamma = d/2. That
  is why Hindmarsh's biased law decays as exp(-kappa eps^2 eta^2) in 2D (quoted from arXiv:1407.3905 §IV).
- Correia-Leite-Martins (arXiv:1407.3905, 1804.10761) found that this Gaussian law fails for biased initial
  conditions. That is consistent with gamma < d/2 in the nonlinear regime.

**Decay criterion and thresholds [D].**
- The network is gone in the tau_10 sense when y_f K^gamma = x10, giving K10(m_f) = (x10/erfinv m_f)^(1/gamma).
- Given a target growth K_target: eps(K_target) = erf(x10 K_target^-gamma) ~ 1.128 x10 K_target^-gamma.
- Survival to parity needs K10 (or K0 = tau_0-based, 1-2.3 x K10 [C network_sim]) > K_par.
- Relic safety needs K10 < K_safe.
- If the bias is quoted as a coherent zero-mode offset phibar over the fluctuation level sigma at freeze-out, then
  m_f = erf(phibar/(sqrt2 sigma)), and phibar/sigma = 1.2533 m_f for small m_f.

**Statistical versus systematic [D].**
- A region of R/L_f formation domains per side has a random initial imbalance ~ (L_f/R)^(3/2).
- Amplified by K^gamma while L < R, it reaches (R/L_f)^(gamma - 3/2) <= O(1) when L ~ R. It never exceeds what the
  scaling state already has.
- The systematic part m_f K^gamma is coherent on all scales and exceeds the random part on the final domain scale
  once m_f > K^-3/2.

## 2. Lattice measurements of gamma [C] (`code/lattice.py`, `code/analyze_lattice.py`)

**Setup.** This is an independent implementation of network_sim's PRS scheme:

- equation phi'' + (alpha/tau) phi' - lap phi + U'(phi) = 0, with U = (1/8)(phi^2 - 1)^2 and units v = m_phi = 1;
- dx = 1, dt = 1/3, tau_i = 2, 3D alpha = 3, 2D alpha = 2;
- Gaussian initial field smoothed with s_c = 1.5, amplitude 0.1, zero velocity.

The overdamped model-A control is phi_t = lap phi - U'(phi) with dt = 0.1.

**Bias and response.**
- Each seed has three runs from one noise realization: the exactly zero-mean field, and the same field shifted by
  +-mu, with m0 = erf(mu/sqrt2).
- Response estimators:
  - R_y = (erfinv m_+ - erfinv m_-)/(2 erfinv m0), erf-linearized;
  - the raw R_m;
  - the unbiased autocorrelation C(t) = <sign phi(t) phi0/std>, which goes as K^-lambda_bar.
- Formation is the first output where the bulk |phi| < 1/2 fraction (minus the wall cores) is below 5% and both
  vacua are reached.
- K = xi_ref/xi_ref(formation), with xi = V/A and an isotropic link count.
- Windows: box-safe tau <= L/4 (PRS) or xi <= N/6 (model A); extended tau <= L/2 or xi <= N/3. Fits use K >= 1.5.

| Run group | Seeds | K reached (safe / ext) | lambda_bar (C, safe / ext) | gamma from response (safe) | Pre-formation gamma | Comment |
|---|---|---|---|---|---|---|
| PRS 3D 160^3, m0 = 0.002 (+ 0.02 on 2 seeds) | 3 | 2.5-2.6 / 4.4-5.4 | **1.610 +- 0.019 / 1.687 +- 0.067** (C_ss: 1.618 / 1.66) | 1.58 +- 0.12 (1.67, 1.45, 1.63); m0 = 0.02: 1.56 +- 0.12 | 1.457 +- 0.038 (theory d/2 = 1.5) | R_y(form) = 4.39 +- 0.09 |
| PRS 2D 1024^2, m0 = 0.001 | 3 | 15-18 / 25-44 | 1.38 +- 0.13 / 1.56 +- 0.12 | 0.43 +- 0.48: **failed** (see below) | 0.968 +- 0.024 (theory 1) | |
| Model A 2D 512^2, m0 = 0.002 | 2 | 6.4-6.6 | 1.146 +- 0.045 (literature asymptote 1.25) | 0.66 +- 0.26 | 1.02 +- 0.14 | |
| Model A 3D 96^3, m0 = 0.002 | 1 | - / 2.8 | ext 1.643 | 3.2 (box-dominated, discarded) | 1.51 | |

**What these show.**

- **Linear (pre-formation) stage.** The growth exponent is d/2, as the Gaussian theory predicts: 1.46 (3D) and
  0.97-1.02 (2D). This checks the estimator. The bias is amplified 4.4-fold in 3D between the quench start and
  formation, consistent with (xi_form/xi_i)^1.41 = 2.76^1.41 = 4.2.
- **Post-formation (nonlinear) stage.**
  - The relativistic 3D autocorrelation gives lambda_bar = 1.61-1.69, i.e. gamma = 1.31-1.39.
  - Model A gives 1.6-1.67 [Q: Liu-Mazenko theory 1.673, Bray Table 1 (lambda = 1.327); brief: 1.6-1.7]. Our
    small model-A box gives 1.64.
  - The direct small-bias response is higher, 1.58 +- 0.12. It covers only K <= 2.6 right after formation, while A is
    still rising, and it agrees with the bound 1.5 at 0.7 sigma.
- **2D model A.** At K <= 7 our model A gives lambda_bar = 1.15, 0.1 below the asymptote 1.25. That is the known slow
  convergence, and it is the size of short-range bias we allow in 3D (band lower edge 1.30).
- **Failure (2D direct response).** In two of three 2D PRS seeds the unbiased reference drifted to m = 0.70 and 0.77
  by tau = L/2. This is the finite-box 2D critical-percolation takeover also seen in network_sim. It saturates the
  response, so the 2D response numbers are not used.
- **Finite-box drift in 3D.** The unbiased 3D runs drift to m = -0.12, 0.24 and -0.02 by tau = L/2. That is the
  statistical imbalance of a 160^3 box amplified by K^gamma, which is why paired runs are needed.

## 3. Fits to network_sim's biased-decay data [C] (`code/fit_K10.py`)

Read-only re-analysis of `network_sim/results/runs/*bias*`. tau_10 is the first time A < 0.1 A_ref(tau). It is now
interpolated in log(A/A_ref) between outputs, and the output-grid bracket is kept as the error. K10 =
xi_ref(tau_10)/xi_ref(tau_form) with network_sim's tau_form. m0 = 2p - 1 is the initial bias (|m0| = 2(1/2 - p)).

| Group | Points | K10 = c m0^-alpha: alpha | c | Notes |
|---|---|---|---|---|
| 3D PRS 128^3 | m0 = 0.1, 0.2, 0.3 (K10 = 2.55, 1.59, 1.15; first-output 2.73, 1.63, 1.27) | **0.722 +- 0.053** | **0.490** (sigma_lnc 0.088) | first-output definition 0.695 / 0.544 (= network_sim's 0.331 x 2^0.70); jackknife alpha 0.68-0.79; local slope of the two smallest biases 0.68 +- 0.08 |
| 2D PRS 1024^2 | m0 = 0.1-0.6 | 1.273 +- 0.048 | 0.569 | local slope of the two smallest biases 1.15 +- 0.11 |
| 2D PRS, damped late start | m0 = 0.1-0.5 | 1.228 +- 0.060 | 0.547 | |

- Theory, 1/(d - lambda_bar):
  - 3D: 0.714 (lambda_bar = 1.60), 0.754 (1.673) and 0.667 (Gaussian bound). With alpha fixed at those values, c =
    0.50, 0.47, 0.55, and every point is fit to <= 3.4% (5.5% for the bound).
  - 2D: 1.333 (lambda_bar = 1.25).
  - The 3D percolation exponent nu = 0.88 is 3.0 sigma from the fit, so the 2D coincidence 1/(2 - 5/4) = 4/3 = nu_2D
    does not carry to 3D.
- **Formation-bias collapse.** All of network_sim's 3D biased runs are already nonlinear at formation: |m| =
  0.44/0.76/0.92 at tau_form for p = 0.45/0.40/0.35. Using y_f = erfinv|m(tau_form)| per run, ln K10 = (ln x10 -
  ln y_f)/gamma gives:
  - 3D: **gamma = 1.406 +- 0.053, x10 = 1.54** (sigma_ln 0.03), chi^2/dof = 5.9/7, 9 runs.
  - 2D: gamma = 0.841 +- 0.011 (x10 = 1.55); damped late start 0.851 +- 0.018 (x10 = 1.62).
  - The 2D agreement between early and deep-subhorizon formation supports gamma not depending on how far inside the
    horizon the walls form.
  - x10 = 1.54 corresponds to m = erf(1.54) = 0.97 at tau_10, as observed (0.977 at p = 0.45).

## 4. Cosmological coarsening factors and thresholds [D, C] (`code/common.py`, `code/thresholds.py`)

**Network length.**
- In scaling, rho_w = A sigma/t, so V/Area = t/A (physical). For a charged network with slowdown factor sQ, A_eff =
  sQ A.
- K(T) = [t(T)/(sQ A)] (g_s(T)/g_s(T_f))^(1/3) (T/T_f) / L_f, with t = 1/(2H) and H from the visible radiation.
- g from the Saikawa-Shirai fits (arXiv:1803.01038 App. C, as in cosmo_observable).
- L_f = g_geo xi0, with g_geo in [1, 3] [E] converting a correlation length into V/Area.

**Parity [D].** 2 A sigma H = rho_rad gives T_par = 5.45/5.82/6.18 keV for A = 0.7/0.8/0.9 without backreaction. With
backreaction T_par = 6.37 keV [Q cosmo_observable]; both are scanned.

**Relic-safe temperature [D].**
- The network's wall energy at decay, plus the emission of the scaling era (kappa = 1-2 [E]: the energy lost per
  Hubble time is ~ rho_w, and rho/s grows as 1/T, so the earlier emission adds ~ rho_w(T_d) once more), becomes stable
  rest mass with fraction f [P]. This must satisfy f kappa (2 A sigma H/s)(T_d) <= 4.373e-10 GeV (Omega h^2 <= 0.12;
  minority_relic's value, retrospective). With 4.35e-10, T_safe moves by 0.5%.
- Closed form: T_safe = f kappa A sigma (45/pi^2) sqrt(pi^2 g_rho/90)/(g_s M_P Y_max).
- T_safe = 30.6 MeV (f = kappa = 1, A = 0.8); 56 MeV (kappa = 2); 3.3 MeV (f = 0.1); 0.35-0.67 MeV (f = 0.01).
  Band 0.31-62 MeV.
- H-star's thermal relic omega_sym = 0.035 already uses part of the budget. Using the remaining 0.085 instead of 0.12
  raises T_safe by 1.4 and eps_safe by 0.2 decades (not applied).

**Histories.** Band = all combinations of: gamma in {1.30, 1.35, 1.41, 1.46, 1.50} with x10(gamma) = 1.49-1.61 refit
on the 3D runs; A = 0.7-0.9; g_geo = 1-3; f = 0.01-1; kappa = 1-2; K0/K10 = 1-2.3 for persistence; both T_par values;
and the listed formation scales.

| History | Formation T_f [inputs] | Formation scale xi0 | log10 K_hor (subhorizon phase) | log10 K_par (band) | log10 K_safe (band) | eps_persist central (log10 band) | eps_safe central (log10 band) |
|---|---|---|---|---|---|---|---|
| H-star | 3.11e6 GeV [P: xi = 5e-4; Q thermal_potential] | xi_KZ = 0.42-2.04 GeV^-1 (z = 2 / z = 1, Ginzburg-matched; log10 N_H = 13.7-15.7) [Q] | 4.70 | 16.70 (15.79-17.13) | 12.86 (11.66-15.28) | **5.8e-24** (-25.3 to -19.8) | **1.5e-18** (-22.7 to -14.9) |
| H-star, all six KZ variants | same | 0.40-3.64 GeV^-1 | | | | (-25.4 to -19.5) | (-22.7 to -14.6) |
| H-star-Q | 5.8e5 GeV [Q REPORT: T_h = 0.29 v] | 21-960 ell (ell = 1/m_X) [Q/E charged_exit_sim]; sQ = 1-2 | 5.97 | 17.24 (15.41-18.16) | 13.40 (11.14-16.31) | 1.0e-24 (-26.9 to -19.5) | 2.5e-19 (-24.2 to -14.3) |
| Inflationary, H_inf = 2.4e4 GeV [P] | 1.4e9 GeV (H_trap = 2.8 GeV) [E] | 1/sqrt(H_inf H_trap), 1-10/m_phi [E] | 1.50 | 16.16 (14.59-17.01) | 12.32 (10.46-15.16) | 3.4e-23 (-25.1 to -18.3) | 8.4e-18 (-22.5 to -13.4) |
| Inflationary, 1e6 GeV | 2.2e8 GeV | same | 3.12 | 16.97 | 13.13 | 2.4e-24 (-26.4 to -19.3) | 6.1e-19 (-23.7 to -14.4) |
| Inflationary, 1e8 GeV | 2.2e7 GeV | same | 5.12 | 17.97 | 14.13 | 9.6e-26 (-27.9 to -20.6) | 2.4e-20 (-25.2 to -15.7) |
| Inflationary, 1e10 GeV | 2.2e6 GeV | same | 7.12 | 18.97 | 15.13 | 3.8e-27 (-29.4 to -21.9) | 9.4e-22 (-26.7 to -17.0) |
| Inflationary, 1e13 GeV | 6.9e4 GeV | same | 10.12 | 20.47 | 16.63 | 2.9e-29 (-31.6 to -23.9) | 7.4e-24 (-28.9 to -19.0) |

Central inputs: gamma = 1.406, x10 = 1.528, A = 0.8, g_geo = 1, f = kappa = 1, the geometric-mean formation scale,
and T_par without backreaction.

**Inflationary-case inputs [P, E].**
- The spectator starts oscillating at H ~ sqrt(lambda_phi) phi_end, with phi_end ~ 1.5 H_inf, and is trapped when
  its amplitude reaches v: H_trap = sqrt(lambda_phi) v^2/phi_end (inflation_spectator §4).
- Instant visible reheating is assumed.
- A preheating restoration may delay formation to T_trap sqrt(eps_fl), with eps_fl in [0.1, 1]; this is scanned.

**One-at-a-time sensitivity for H-star [C]** (other inputs central):

| Input | Values | eps_persist | eps_safe |
|---|---|---|---|
| gamma | 1.30 / 1.327 / 1.35 / 1.39 / 1.41 / 1.46 / 1.50 | 3.3e-22 / 1.2e-22 / 4.4e-23 / 1.1e-23 / 5.8e-24 / 7.8e-25 / 1.6e-25 | 3.2e-17 / 1.5e-17 / 6.9e-18 / 2.3e-18 / 1.5e-18 / 3.1e-19 / 9.4e-20 |
| f, kappa | f = 1, 0.1, 0.01 (kappa = 1) | - | 1.5e-18 / 6.0e-20 / 2.5e-21 (kappa = 2: x 2.4-2.6) |
| A | 0.7 / 0.8 / 0.9 | 4.4e-24 / 5.8e-24 / 7.5e-24 | 1.0e-18 / 1.5e-18 / 2.0e-18 |
| KZ variant | z1 naive ... z2 mean-field (log10 N_H 12.9-15.8) | 4.0e-23 ... 1.8e-24 | 1.0e-17 ... 4.5e-19 |

**Stress and comparison values (not part of the band).**
- The review's naive extrapolation of the 2D fit gives |p - 1/2| = 9.8e-15 at K = 1e17 (m = 3.4e-14 at H-star's
  K_par).
- If 3D had the 2D relativistic exponent gamma = 0.84: eps_persist = 1.6e-14 and eps_safe = 2.8e-11.
- network_sim's 3D fit in its own convention (initial m0, c = 0.49, alpha = 0.722) gives m0 = 2.9e-24 at H-star's
  K_par. That is consistent with eps_persist = 5.8e-24 at formation, given the lattice's quench amplification of 4.4
  and the slightly different exponent.

## 5. Residual biases from the review, classified [C, conditional on the bias being systematic over the patch]

| Candidate (source) | m_f | log10 K10 | Fate (H-star central) |
|---|---|---|---|
| Thermal-restoration remnant, 1e-3 of the soft-mode level, undamped [E, review] | 8.0e-4 | 2.37 | (iii) decays at K ~ 230, before horizon size (K_hor = 10^4.7) |
| Post-inflation trapping p_- = 0.4999 / 0.4997 [Q inflation_spectator] | 2e-4 / 6e-4 | 2.80 / 2.46 | (iii), before horizon size |
| Forced charge at 2 n_c, p(+) = 0.71 [Q charged_exit_sim] | 0.42 | 0.42 | (iii); the lattice already saw the minority vanish |
| Relaxed start p_- = 0.23 at end of inflation [Q] | 0.54 | 0.33 | (iii) |

**Damping requirement [D].** For a remnant of 1e-3 to fall below eps_persist (the persistent-network class) it must be
damped by 46 amplitude e-folds (band 39-51) before freeze-out. Reaching only eps_safe (the relic-overproduction
class) takes 34 e-folds (27-45). Whether the hidden bath supplies this (zero-mode damping rate times time above T_c)
is not computed here. It is the open input that decides the class of a hot history.

**What "systematic" means here.** In the inflationary case the residual p - 1/2 is a property of the end-of-inflation
distribution, and the long-wavelength spectator modes make the local bias vary on super-horizon scales. A spatially
varying (random-field-like, long-range correlated) bias is not the uniform m_f of this module. Its treatment
(Bray-Humayun-Newman-type correlated initial conditions) is not done.

## 6. Conditional statements for the report

- (C1) With the scaling law calibrated in 3D (gamma = 1.41 +- 0.05, band 1.30-1.50), a hot history whose exit is
  unbiased to better than eps_persist ~ 1e-23 (band 1e-25 to 1e-20) at formation ends with a persistent network. The
  parity rejection then applies.
- (C2) A systematic formation bias between eps_persist and eps_safe ~ 1e-18 (band 1e-23 to 1e-15) gives a network
  that decays after T_safe = 0.3-62 MeV. Its stable h/X/D relic exceeds Omega h^2 = 0.12 (retrospective; f in
  [0.01, 1] [P]), so that history is rejected too.
- (C3) A systematic bias above eps_safe removes the network early. For every residual bias in the review (>= 2e-4)
  that happens before the network reaches horizon size. The wall-parity and relic tests then do not reject the
  history; other consequences (GW, BBN details) are not computed.
- (C4) The draft's identification "hot history => persistent network => rejected" holds only in sub-case (C1). It
  requires erasure of the sign memory to ~1e-23, i.e. tens of e-folds of zero-mode damping in the restored phase.
  That was argued for a kinetically thermal bath but not computed.
- (C5) An exact-Z2 (statistically unbiased) ensemble is case (C1): finite-volume fluctuations do not act as a bias.
- None of this is a prediction, validation or exclusion of Candidate 1. The bias, the hidden temperature, the
  inflationary scale and v are free preparation data.

## 7. Limits, failures, not done

1. **Extrapolation.** The 3D exponent is measured over K <= 2.7-5.4. Using it at K ~ 1e11-1e21 assumes:
   - dynamical scaling with a single length;
   - no change of universality through the sub-horizon (flat-space or friction-dominated) phase and the scaling
     phase;
   - no late-time crossover.

   Friction-dominated, overdamped motion is model-A-like (gamma = 1.33-1.40), inside the band. The bound
   gamma <= d/2 = 1.5 caps how fast a bias can grow, under scaling. Nothing caps how slowly it grows: a late
   crossover to smaller gamma would raise both thresholds (stress case gamma = 0.84: eps_persist 1.6e-14).
2. **PRS only.** The physical 3+1 equation (shrinking comoving wall width) was not run for the response. network_sim
   found PRS and physical A in agreement for unbiased runs.
3. **The direct 3D small-bias response** gives gamma = 1.58 +- 0.12, above the collapse and autocorrelation values. It
   covers only K <= 2.6 just after formation, with 3 seeds, and is not used for the central value.
4. **Failed.** The 2D PRS direct response (finite-box percolation drift in 2 of 3 seeds) and the 3D model-A response
   (96^3 box too small) are recorded but unused.
5. **The KZ-to-V/A conversion** g_geo = 1-3 and the emission factor kappa = 1-2 are estimates. The formation scales of
   H-star-Q (21-960 ell) and of the inflationary case are wide [E].
6. **Not computed:**
   - the actual residual bias of any history (thermal zero-mode damping, the spatial structure of the inflationary
     bias);
   - wall dynamics with charge on the walls beyond the A_eff = sQ A rescaling;
   - D;
   - black-hole formation from late large domains;
   - GW signals of early-decaying networks;
   - the relic budget reduced by H-star's thermal omega_sym (would move eps_safe by +0.2 decades).
7. **Literature.**
   - Bray's review was read at section level from the arXiv PDF (kept in scratch, not in the repository).
   - The 3D simulation values of lambda_bar are not from a read source: Liu-Mazenko 1.673 is from Bray's Table 1, and
     1.6-1.7 is the task brief's range.
   - Fisher-Huse 1988, Huse 1989, Liu-Mazenko 1991, Bray-Kissner 1992 and Humayun-Bray 1991 were not read; they are
     cited via Bray.

## 8. Literature

| Id | Read level | Used for |
|---|---|---|
| Bray, Adv. Phys. 43 (1994) 357, arXiv:cond-mat/9501089 | sections: §2.1 eq. (9) (lambda_bar), §4.2 eqs. (85)-(95) (response identity S = Delta G, lambda_bar = d - lambda), §5.1 Table 1 and text (lambda: 2D sims 0.75, Mazenko 0.711 / 1.327 for d = 2 / 3), §5.2 "External fields/initial bias" eqs. (172)-(175) | growth law, erf form, exponent values |
| Correia, Leite, Martins, PRD 90 (2014) 023521, arXiv:1407.3905 | abstract + §IV (via ar5iv HTML) | biased initial conditions; Hindmarsh's law exp(-kappa eps^2 eta^2) fails for biased IC |
| Correia, Leite, Martins, PRD 97 (2018) 083521, arXiv:1804.10761 | abstract | the failure is attributed to the Gaussian approximation |
| Hindmarsh, PRL 77 (1996) 4495, arXiv:hep-ph/9605332 | abstract (search snippet) | relativistic Gaussian ansatz (gamma = d/2) |
| Larsson, Sarkar, White, PRD 55 (1997) 5129, arXiv:hep-ph/9608319 | abstract | biased-IC decay context |
| Fisher, Huse PRB 38 (1988) 373; Liu, Mazenko PRB 44 (1991) 9185; Bray, Kissner J. Phys. A 25 (1992) 31; Humayun, Bray J. Phys. A 24 (1991) 1915 | not read (cited via Bray) | lambda_bar, its bound, response in the initial regime |
| Saikawa, Shirai arXiv:1803.01038 | not read here; fit coefficients copied from cosmo_observable/code/core.py | g_rho, g_s |
| Hiramatsu, Kawasaki, Saikawa arXiv:1309.5001 | via network_sim (A = 0.8 +- 0.1) | A |

## 9. Files

- `code/lattice.py`: lattice runs (PRS and model A, paired bias runs). `code/analyze_lattice.py`: gamma and
  lambda_bar.
- `code/fit_K10.py`: re-analysis of network_sim's biased runs, the c, alpha fits and the formation-bias collapse.
- `code/common.py`: constants, g fits, H, t, parity, T_safe, K(T). `code/thresholds.py`: K targets, thresholds,
  bands, candidates, sensitivity. `code/run_all.py`: orchestration and the headline.
- `results/lattice_runs/*.json` (31 runs), `results/lattice_response.json`, `results/fit_K10.json`,
  `results/thresholds.json`, `results/headline.json`, `results/runtimes.json`.
