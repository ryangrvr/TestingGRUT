# Verification v_bias_threshold: module bias_threshold

21 October 2026. Adversarial verifier. I wrote all code here myself: `fit_ns_v.py`, `lattice_v.py`,
`analyze_lattice_v.py`, `cosmo_v.py` and `run_queue_v.sh`. I did not import or copy module code.

- I read the module code (`common.py`, `fit_K10.py`, `thresholds.py`) to understand definitions. My
  implementations differ on purpose (see Method).
- Libraries: numpy and scipy only. At most 2 single-thread workers.
- Longest job: one 256^3 lattice run, 309 s. Total lattice CPU about 27 min. Everything else took under 1 min.
- No git was used. Nothing outside this directory was written.
- Bray's review PDF was downloaded to the scratchpad only.

Results:
- `results/fit_ns_v.json`: re-analysis of network_sim's biased runs.
- `results/lattice_v/*.json`: 12 own runs.
- `results/lattice_v_analysis.json`.
- `results/cosmo_v.json`.

**Bottom line.**
- Every headline number I could recompute is reproduced, mostly to 3 digits.
- The theory the module relies on is quoted correctly from Bray (read at section level here). This covers the
  response identity S = Delta G, lambda_bar = d - lambda, the k = 0 response growing as L^(d - lambda_bar), and
  Table 1 with Mazenko's lambda = 1.327 in 3D. The Fisher-Huse bound gamma <= d/2 follows correctly from
  Cauchy-Schwarz under scaling.
- My own 3D lattice independently gives d - lambda_bar = 1.43 +- 0.04 (box-safe) and 1.40 (to L/2). The direct
  small-bias response is 1.48-1.55, consistent with the module.
- The central gamma = 1.41 is therefore confirmed **as a value in the lattice window K <= 3-7**.
- The thresholds eps_persist ~ 6e-24 and eps_safe ~ 1.5e-18 are correct arithmetic. They remain conditional
  extrapolations over 15-17 and 11-15 decades.
  - The ±0.053 on gamma is statistical only. The spread between estimators is about ±0.08.
  - The band edge gamma = 1.30 is a judgment, not a bound. If gamma is pushed below 1.30 (see V-G4),
    eps_persist rises to about 1e-20 to 1e-18.
- The qualitative conclusions are robust to any gamma >= 0.32: three fate classes, and every listed candidate
  bias falls in class (iii).
- Corrections are minor:
  - The trapping bias values are |m| = 1.4e-4 to 6.5e-4, not 2e-4 to 6e-4.
  - They come from an H_inf = 5-10 TeV inflationary history, but they are classified under H-star.
  - The lambda_bar error bar understates the seed scatter.
  - The erf crossover form is an ansatz, not a derivation.

## Method (independent)

1. **network_sim re-analysis (`fit_ns_v.py`).** These choices differ from the module:
   - The unbiased reference is **paired per seed**. network_sim uses the same white noise for every p, shifted by
     a quantile. The module uses the seed-mean reference instead.
   - tau_10 is found by log-log interpolation of A/A_ref.
   - gamma is computed by three estimators:
     - G1: regression of ln K10 on ln erfinv|m(t_ref)|, plus the inverse regression and a seed bootstrap;
     - G2: in-run trajectory slopes d ln erfinv|m| / d ln K up to |m| = 0.97;
     - G3: the reference output varied over 5 outputs around network_sim's tau_form.
   - Four length measures: tau/2A (PRS area), tau/2A_iso, 1/(raw link count), and conformal time tau.
2. **Own 3D PRS lattice (`lattice_v.py`).** Same physics as the module, different code:
   - equation phi'' + (3/tau) phi' - lap phi + U'(phi) = 0, U = (1/8)(phi^2 - 1)^2;
   - friction integrated exactly with P = tau^3 phi'; kick-drift-kick, dt = 0.25, dx = 1;
   - Gaussian IC with s = 1.5 and std 0.1, own seeds 101-104, tau_i = 2;
   - runs: 3 seeds at 192^3 and 1 seed at 256^3. Each seed has a paired triple mu = 0, +-2.507e-4, i.e.
     m0 = 0.002.
   - K = nl(form)/nl(t), with nl the raw sign-changing link count of the unbiased run.
   - Formation is the time when the fraction of cells with |phi| < 1/2 falls below 0.2. Variants 0.3 and 0.12 were
     also run.
   - Measured quantities: gamma from R = (erfinv m+ - erfinv m-)/2, and lambda_bar from C = <sign phi phi0>/std.
   - Windows: box-safe tau <= N/4 and extended tau <= N/2.
3. **Cosmology (`cosmo_v.py`).** My own g_rho and g_s come from ideal-gas FD/BE integrals of all SM species, not from
   the Saikawa-Shirai fit.
   - The QCD crossover is a tanh switch between a quark-gluon gas and a pion gas.
   - Neutrino decoupling and e+e- annihilation are handled by entropy conservation (N_eff = 3 or 3.044).
   - Formulas, all written from scratch:
     - rho_w = A sigma/t and t = 1/2H;
     - parity: rho_w = rho_rad;
     - relic-safe decay: f kappa rho_w/s <= Y_max, with Y_max = 0.12 rho_c/(h^2 s0) = 4.3733e-10 GeV;
     - K = (t/A)(a_f/a)/L_f.
   - Thresholds: eps = erf(x10 K^-gamma).
4. **Literature (read level).**
   - Bray, arXiv:cond-mat/9501089 (section level, PDF in scratch): §2 eq. (9); §4.2 eqs. (85)-(95); §5.1 Table 1 and
     the text after eq. (143), "d = 2 simulations give lambda ≃ 0.75"; §5.2 eqs. (173)-(175) and the following text
     (initial-bias erf form, and "t^(d/4) gets replaced by L^lambda" in the initial regime).
   - Midya-Majumder-Das, J. Phys. Condens. Matter 26 (2014) 452202 (abstract via repository record): 3D
     nonconserved lambda_bar agrees with the Gaussian-auxiliary-field (Liu-Mazenko) value.
   - Chakraborty-Das, EPL 119 (2017) 50005 (abstract): Fisher-Huse bound lambda_bar >= d/2.
   - Parwani, arXiv:hep-ph/9204216 (eq. 3.29 via ar5iv; only for the residual_bias flag below).

## Claims and verdicts

### Theory

| # | Claim | My check | Verdict |
|---|---|---|---|
| T1 | S(k,t,0) = Delta G(k,t) exactly (Bray §4.2 eq. 93), valid for any deterministic dynamics | Read eq. (93) and the text "exact result ... integration by parts on the gaussian distribution". The proof uses only the Gaussian initial distribution, so it covers second-order dynamics with fixed initial velocity. | CONFIRMED |
| T2 | G = L^lambda g_R(kL) and lambda_bar = d - lambda (eq. 95); the k = 0 response grows as L^(d - lambda_bar) | Read eqs. (94)-(95). Own lattice: gamma_direct + lambda_bar = 1.48 + 1.57 = 3.05 ± 0.07 (box-safe, K measure), consistent with d = 3. | CONFIRMED |
| T3 | Literature values: Mazenko lambda = 1.327 (3D), so lambda_bar = 1.673; 2D simulations lambda ≈ 0.75 | Table 1 (d = 3, n = 1: 1.327; d = 2: 0.711) and the text after eq. (143) | CONFIRMED (read level: section) |
| T4 | 3D simulation lambda_bar ≈ 1.6-1.7 | Not from a read table. The Midya et al. 2014 abstract says the 3D values agree with the Gaussian-auxiliary-field theory (about 1.67). The module labels its own source honestly. | CONFIRMED at abstract level |
| T5 | Bound lambda_bar >= d/2, i.e. gamma <= 3/2 (Cauchy-Schwarz) | Derivation re-done: S(0,t,0)^2 <= S(0,t) S(0,0) with S(0,t) ~ L^d gives L^(2(d - lambda_bar)) <= L^d. It needs scaling and short-range initial correlations, and holds only asymptotically. The early-time lattice response exceeds it (1.48-1.58). | CONFIRMED (asymptotic only) |
| T6 | Law m(K) = erf(y_f K^gamma) "derived" | The erf form is Bray eq. (175) for the **Gaussian** theory. The exact statement covers only the small-m initial regime: L^(d/2) is replaced by L^lambda. The full crossover function is an ansatz. In network_sim's 3D runs the in-run slope d ln erfinv|m|/d ln K (G2) is 1.51, while the cross-run endpoint relation (G1) gives 1.39. A single-gamma erf form does not fit both. This affects x10 (factor ≲ 1.5, i.e. ≲ 0.2 decades in eps), not the power law. | CORRECTED (label: ansatz calibrated on data, not derived) |
| T7 | Statistical imbalance (xi0/R)^(3/2) amplified to (R/L_f)^(gamma - 3/2) <= O(1); finite-volume fluctuations cannot remove the network | Correct counting (N^(-1/2)), and the response at kL < 1 is ~ L^gamma. Valid given gamma <= 3/2. | CONFIRMED (conditional on scaling) |

### Fits to network_sim's biased runs (3D PRS 128^3, 2D PRS 1024^2)

| # | Claim (module) | My value | Verdict |
|---|---|---|---|
| F1 | 3D formation-bias collapse: gamma = 1.406 ± 0.053, x10 = 1.54 (9 runs, chi2/dof 5.9/7) | **1.392** with paired reference (bootstrap 16-84%: 1.27-1.52), x10 = 1.535. Inverse regression 1.349. With the reference output moved over tau = 9.0-14.3: 1.377-1.394. Area estimator (PRS, iso, raw links): 1.372-1.392. The 9 runs give only **3 distinct bias values**, and seeds are nearly identical, so "dof = 7" overstates the information. | CONFIRMED (value); the error is statistical only |
| F2 | x10 = 1.54 means m = erf(1.54) = 0.97 at tau_10, as observed | \|m\|(tau_10) = 0.968-0.978 over the runs | CONFIRMED |
| F3 | K10 = c m0^-alpha, 3D: alpha = 0.722 ± 0.053, c = 0.49; local slope of the two smallest biases 0.68 | alpha = 0.727, c = 0.486 (paired, unweighted); smallest two 0.694 | CONFIRMED |
| F4 | 2D PRS: alpha = 1.273 ± 0.048, c = 0.57; damped late start 1.228 | 1.247 / 0.590; damped 1.263 / 0.510 | CONFIRMED (within 1 sigma) |
| F5 | 2D collapse gamma = 0.841 ± 0.011 (x10 1.55); damped 0.851 ± 0.018 | 0.846 (x10 1.56); damped 0.825 (paired ref) | CONFIRMED (range 0.82-0.85) |
| F6 | nu_3D = 0.88 lies 3.0 sigma from alpha | (0.88 - 0.722)/0.053 = 2.98 | CONFIRMED |
| F7 | (not claimed) length-measure systematic | With conformal time instead of xi = tau/2A, gamma becomes 1.14-1.19 (3D) and 0.74 (2D). The reference network has not reached scaling: A_ref rises from 0.65 to 0.75 over the window, so L ~ tau^0.89. In my own lattice the response identity gamma + lambda_bar = d holds in the K (area) measure (3.05) and fails in tau (2.67). So the module's K measure is the self-consistent one, and the tau-based values are not a valid lower edge. The real caveat is that **the whole calibration window is pre-asymptotic**. | NOTE (supports the module's choice; flags pre-asymptotic window) |
| F8 | prs2d_biaslate group | Not usable: x10 = 8-650, K10 flat at about 4 for m0 >= 0.3. The module did not use it either. | n/a |

### Own 3D lattice (4 seeds; module values in brackets)

| # | Quantity | Mine | Module | Verdict |
|---|---|---|---|---|
| L1 | lambda_bar box-safe (K <= 3.3) | **1.571 ± 0.040** (per seed 1.50, 1.60, 1.67, 1.51) | 1.610 ± 0.019 | CONFIRMED (0.9 sigma). The module's ±0.019 from 3 seeds understates the seed scatter: my per-seed std is 0.08. |
| L2 | lambda_bar to L/2 (K <= 6-7) | 1.602 ± 0.037 | 1.687 ± 0.067 | CONFIRMED (1.1 sigma) |
| L3 | direct small-bias gamma (m0 = 0.002), box-safe | 1.483 ± 0.056 (extended 1.554 ± 0.060) | 1.58 ± 0.12 | CONFIRMED. It sits at or above the bound 1.5, which shows the window is pre-asymptotic. |
| L4 | pre-formation (linear) gamma, theory d/2 = 1.5 | 1.469 ± 0.014 | 1.457 ± 0.038 | CONFIRMED |
| L5 | 3D quench amplification of the bias, initial to formation | 4.76 ± 0.06 (formation at frac < 0.2); 4.24 (at 0.3); 11.4 (at 0.12) | 4.39 ± 0.09 | CONFIRMED. The value depends on the formation definition (factor 2.7 across 0.3-0.12). |
| L6 | "relativistic walls follow the model-A exponent within the errors" | d - lambda_bar = 1.40-1.43 against model A 1.33-1.40 | | CONFIRMED in the window K <= 7 only |

### Coarsening factors, temperatures, thresholds (H-star central: T_f = 3.111e6 GeV, xi0 = 0.9256 GeV^-1, A = 0.8)

| # | Claim | Mine | Verdict |
|---|---|---|---|
| C1 | KZ inputs xi_KZ = 0.42-2.04 GeV^-1 (z2/z1 Ginzburg-matched), log10 N_H = 13.7-15.7 | thermal_potential row xi = 5.01e-4: 0.4203 / 2.0385; 15.73 / 13.67. H(T_f) = 1.3599e-5 GeV in both codes. | CONFIRMED (inputs quoted correctly) |
| C2 | T_par = 5.45 / 5.82 / 6.18 keV (A = 0.7 / 0.8 / 0.9) | 5.454 / 5.831 / 6.185 (N_eff 3); 5.446 / 5.822 / 6.176 (N_eff 3.044) | CONFIRMED |
| C3 | T_safe = 30.6 MeV (f = kappa = 1); 56 MeV (kappa = 2); 3.3 MeV (f = 0.1); 0.35-0.67 MeV (f = 0.01); band 0.31-62 MeV | 30.59; 56.3; 3.26; 0.348-0.665; band 0.309-62.6. The closed form matches the root (my g = 12.65 / 12.43 at T_safe). Y_max = 4.3733e-10 GeV reproduced from rho_c and s0. | CONFIRMED |
| C4 | log10 K_par = 16.70 (band 15.79-17.13); log10 K_safe = 12.86 (11.66-15.28); log10 K_hor = 4.70 | 16.695 (15.786-17.125); 12.855 (11.659-15.281); 4.696 | CONFIRMED |
| C5 | eps_persist = 5.8e-24; eps_safe = 1.5e-18 (central) | 5.80e-24; 1.45e-18 | CONFIRMED (arithmetic; conditional extrapolation) |
| C6 | Sensitivity to gamma: 1.30: 3.3e-22 / 3.2e-17; 1.41: 5.8e-24 / 1.5e-18; 1.50: 1.6e-25 / 9.4e-20 | 3.33e-22 / 3.27e-17; 5.80e-24 / 1.45e-18; 1.64e-25 / 9.47e-20 | CONFIRMED |
| C7 | Band of eps_persist 10^-25.3 to 10^-19.8; eps_safe 10^-22.7 to 10^-14.9 (gamma 1.30-1.50 plus nuisance inputs) | The edges are reproduced by my nuisance scan. Widening gamma to 1.15-1.50 (stress) gives -25.4 to -17.5 and -22.7 to -13.2. Lattice estimators in the self-consistent K measure span 1.39-1.55. Model A asymptotically gives 1.33. **No bound protects the low edge.** | CONFIRMED as computed. The band is a judgment over gamma, not a bound. |
| C8 | H-star-Q: log10 K_par 17.24, K_safe 13.40, K_hor 5.97; eps 1.0e-24 / 2.5e-19 | 17.239, 13.399, 5.969; 9.97e-25 / 2.50e-19 | CONFIRMED |
| C9 | Inflationary cases: H_trap = sqrt(lambda) v^2/phi_end (amplitude ~ 1/a, H ~ a^-2); eps_persist 3.4e-23 (2.4e4 GeV), 9.6e-26 (1e8), 2.9e-29 (1e13) | H_trap formula re-derived. 3.31e-23, 9.44e-26, 2.88e-29; T_trap(1e8) = 2.18e7 GeV. Note that 1/sqrt(H_inf H_trap) = 3.87e-3 GeV^-1 is independent of H_inf. | CONFIRMED (inputs are [P]/[E]) |
| C10 | Extrapolation span: parity 15.1-16.7, relic-safe 10.9-14.9 decades | 15.79 - log10 5.4 = 15.06; 17.13 - log10 2.69 = 16.70; 11.66 - 0.73 = 10.93; 15.28 - 0.43 = 14.85 | CONFIRMED |
| C11 | Statistical floor K_par^-3/2 = 9e-26, 1.8 decades below eps_persist | 9.06e-26; 1.81 decades | CONFIRMED |
| C12 | Damping for a 1e-3 remnant: 46 (39-51) e-folds to the persistent class; 34 (27-45) to the relic class | 46.37; 33.94. With gamma = 1.19 (stress): 38.1 / 27.6. | CONFIRMED |
| C13 | Review's naive 2D extrapolation 9.8e-15 at K = 1e17; stress gamma = 0.84: 1.6e-14 / 2.8e-11 | (1e17/0.262)^(-1/1.255) = 9.79e-15; 1.66e-14 / 2.78e-11 | CONFIRMED |
| C14 | phibar/sigma = sqrt2 erfinv(m) = 1.2533 m (small m) | P(phi > 0) for a Gaussian with mean phibar gives m = erf(phibar/(sqrt2 sigma)) | CONFIRMED |

### Candidate residual biases (module §5)

| # | Claim | Mine | Verdict |
|---|---|---|---|
| R1 | Thermal remnant 1e-3: m_f = 8.0e-4, log10 K10 = 2.37 | 7.98e-4, 2.372. If "1e-3" refers to the start of the quench rather than freeze-out, the lattice amplification adds ×4-5. That moves log10 K10 by -0.5 and the damping requirement by +1.5 e-folds. | CONFIRMED |
| R2 | Trapping p_- = 0.4999 / 0.4997 gives m_f = 2e-4 / 6e-4 | inflation_spectator `final_sign.json`: p_- after trapping = 0.49993 / 0.49969 (N50, H_inf = 5e3 / 1e4 GeV) and 0.49970 / 0.49967 (N60). So **\|m\| = 1.4e-4 to 6.5e-4**, log10 K10 = 2.42-2.91. These are fresh-start (+v delta at N_*) ensemble means with instant RD and no preheating [P]. They belong to an H_inf = 5-10 TeV history, not to H-star, under which the module classifies them. Equilibrium draws give exactly 1/2. Locally the bias varies on super-horizon scales: the never-crossing core \|phi_end\| < 2.4 v carries the asymmetry. The module itself says a spatially varying bias is not treated. | CORRECTED (values 1.4-6.5e-4; history label). The fate (iii) is unchanged. |
| R3 | Forced charge p(+) = 0.71 gives 0.42; relaxed start p_- = 0.23 gives 0.54 | arithmetic 0.42 / 0.54; log10 K10 0.42 / 0.33. The relaxed-start 0.23 is an end-of-inflation value at H_inf ~ 3-5 TeV, before the trapping randomization of inflation_spectator §4.3. | CONFIRMED (arithmetic; provenance as stated) |
| R4 | All candidates are class (iii), before horizon size (K_hor = 10^4.7) | Using the smallest corrected candidate, \|m\| = 1.4e-4: for gamma = 1.406, log10 K10 <= 2.91; for gamma = 1.15, <= 3.52; for gamma = 0.84, <= 4.88. Class (iii) (K10 < K_safe = 10^12.86) holds for any gamma > 0.32. "Before horizon size" fails only for gamma <~ 0.87. | CONFIRMED (robust) |

### Not checked

- **residual_bias** (soft-mode damping rate, the integral of Gamma dt, the coherent-residual-to-m0 conversion): NOT
  CHECKED. While I verified, the module had no NOTES.md and its files were still being written (08:55-08:59). It is
  not the module assigned here.
  - One flag for its verifier: Parwani (hep-ph/9204216 eq. 3.29) gives the **pole (amplitude)** damping
    gamma = g^2 m/(64 pi) = lambda^2 T^2/(1536 pi m), with lambda/4! normalization.
  - residual_bias's `parwani_check` value Gamma/T = 2.98e-5 (g^2 = 0.06, m/T = 0.05) equals lambda^2 T/(768 pi m),
    i.e. 2 gamma (the width).
  - If that Gamma is used as an amplitude rate, the e-folds are doubled. That matters for its "hard-rate floor"
    branch (D = 27.5, i.e. about 13.8 amplitude e-folds), not for D ~ 3000.
- The K0/K10 = 1-2.3 input from network_sim, g_geo = 1-3 [E] and kappa = 1-2 [E] were taken as given. They are
  inside the scanned band.

## What the module overstates or omits

1. "gamma = 1.406 ± 0.053": the error is the fit error of one estimator on 3 distinct bias values. The spread across
   estimators in the self-consistent K measure (collapse 1.39, d - lambda_bar 1.40-1.43, direct response 1.48-1.58,
   in-run slope 1.51) is about ±0.08. All of it comes from a pre-asymptotic window (L ~ tau^0.89; local exponents
   noisy for K > 3). The thresholds band uses gamma 1.30-1.50, which covers this. The headline ± should be read as
   statistical.
2. "Law ... derived": only the small-bias power law with exponent d - lambda_bar is exact (under scaling). The erf
   crossover is a Gaussian-theory ansatz (T6).
3. The candidates table mixes histories: the trapping biases come from H_inf = 5-10 TeV but are classified with
   H-star's K targets. This is harmless for the fate.
4. The extrapolation also assumes things not tested on any lattice: the same gamma through the 4.7-decade
   sub-horizon phase (thermal friction from the hidden bath, deep-subhorizon formation in 3D), and through the QCD
   and e+e- epochs (only g changes there, which is included in K). The 2D damped late-start test supports
   insensitivity to the formation depth, but only in 2D.
