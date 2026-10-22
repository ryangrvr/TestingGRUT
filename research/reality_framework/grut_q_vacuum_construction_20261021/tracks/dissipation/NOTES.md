# Track `dissipation`: can mu be formed by dynamics? Klinkhamer-Volovik dissipative q-theory

21 October 2026. Stage-2 track of the q-theory construction (`../../BRIEF.md`, `../../BRIEF_STAGE2.md`). This is the
"formed" side of the protect-or-form dichotomy (`../ledger/NOTES.md` §1.5): if mu is conserved it must be prepared
(QV1 with P_q); here the question is whether a dissipative extension forms mu, so that P_q is not needed.

**Labels.** [D] derived here; [C] computed here (code in `code/`, output in `results/`); [Q] stated in the cited source
(paraphrased, with the source's equation numbers); [E] estimate; [A] argued; [R] recalled, not checked here.
No new theory, empirical validation or confirmatory test is claimed. Every comparison with data is **retrospective**.
No parameter is calibrated to rho_Lambda,obs anywhere in this track.

**Sources and read levels** (TeX from arXiv, kept in the scratchpad
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/dissipation_dl/`; nothing
copied into the repository):

| Source | Identifier | Read level |
|---|---|---|
| Klinkhamer, Volovik, "Dynamic cancellation of a cosmological constant and approach to the Minkowski vacuum" (KV16a) | arXiv:1601.00601v6, MPLA 31, 1650160 (2016) | **full** (TeX, incl. App. A and App. B; the literature track had not read App. B) |
| Klinkhamer, Savelainen, Volovik, "Relaxation of vacuum energy in q-theory" (KSV16) | arXiv:1601.04676v4, JETP 125, 268 (2017) | **full** (TeX, incl. appendix) |
| Klinkhamer, "Q-field from a 4D-brane: cosmological constant cancellation and Minkowski attractor" (K22) | arXiv:2207.03453v4, LHEP 2022, 312 | **full** (TeX) |
| Klinkhamer, "Extension of unimodular gravity and the cosmological constant" | arXiv:2207.02826v6, PRD 106, 124015 | section (§§I, VII, by search in the TeX) |
| Stage-1 records | `../literature/NOTES.md`, `../equations/NOTES.md` (+ code), `../ledger/NOTES.md` (+ `results/ledger.json`) | full |

Figures of the papers were not viewed as images; captions, tables and text were used, and every figure was
re-computed from its stated parameters (§3).

---------------------------------------------------------------------------------------------------------------------

## 0. Result in one screen

1. **Common structure** [D, C]. All three KV models replace the conserved Gauss constant by an effective
   mu_eff(t) with d mu_eff/dt = S. The off-shell identity of stage 1 (`../equations/NOTES.md` §1.4) then forces the
   matter equation rho_m-dot + 3H(rho_m + p_m) = q S: the dissipated vacuum energy goes into matter (massless quanta in
   KV16a). With that source the Bianchi identity and the Friedmann constraint hold exactly (`s1` checks A, C:
   C-dot = -3hC; without the source C is violated at rate f s).
2. **What is derived and what is not** [Q, D]. KV16a's rate (gamma/36)|H|R^2 rests on the Zel'dovich-Starobinsky R^2
   production law (derived, global, for asymptotically static in/out regions). The local form, E ~ |H|, the
   coefficient gamma = 1/(32 pi^2) as a local coefficient, and the assignment of the produced energy to rho_V (rather
   than to the geometry's backreaction in general) are heuristic or phenomenological. KSV16's S is a pure ansatz. K22's
   dissipative equation is not a deformation of its own action's field equation (§2.4).
3. **All three papers are reproduced** [C] (§3): KV16a r_V(1e5) = 3.0285e-3 (paper 3e-3) and App. B exact solution to
   3e-12; KSV16 u_eff(98) = -0.333276 (paper -0.333276), Table I Y/N grid identical in all three slices, Table II to the
   printed digits; K22 Table III to the printed digits.
4. **The "conflict" is not one** [D, C] (§4.1). Both dissipation functions are non-negative, so one quantity is
   monotone (rho_V non-increasing in KV16a; u_eff non-decreasing in KSV16). Whether Minkowski attracts depends only on
   whether de Sitter states are in the zero set of S: |H|R^2 is nonzero in de Sitter (de Sitter is assumed to radiate),
   so for conformal matter R = 0 <=> rho_V = 0 is the only fixed point and **the basin is exactly the half-space
   r_V,i >= 0, H_i > 0** (finite measure: 40.2% of the stated (lambda, f_i) box; the boundary is r_V,i = 0 for every
   lambda in [-1, 1]; dust or any trace moves it up, to r_V,i ~ gamma r_d,i^2/24). The KSV16 S
   vanishes on every static state, so the de Sitter family is a continuum of fixed points and Minkowski is a
   codimension-1 separatrix (measure zero). The approach is one-sided in both: no trajectory that starts below (or
   overshoots) the Minkowski value comes back.
5. **Too slow by 121 orders** [D, C] (§5). For radiation plus vacuum the KV16a law is exact and closed:
   d rho_V/dN = -(4 gamma/9) rho_V^2/M^4, so rho_V(N) = rho_Vi/(1 + (4 gamma/9)(rho_Vi/M^4) N) >= (1/2) min(rho_Vi,
   9M^4/(4 gamma N)). QV1's K'(q) != 0 changes nothing: the deviation from this law falls as ~(H/m_q)^6 (§4.5). Reaching rho_Lambda needs **N = 9.9e122 e-folds**; a standard history from T = M to T_CMB
   has 72.5. On the attractor the universe stays vacuum dominated (Omega_M = 1/(4N)): at the observed age
   rho_V(t0) = 2.24e35 GeV^4 = **8.9e81 rho_Lambda**, H(t0) = 7.8e40 H0, H t0 = 7.4e40 (observed 0.95), Omega_m = 2e-42
   (observed 0.315). Retrospective failure of the attractor history.
6. **One-sided, so later vacuum releases are fatal** [D, C]. The drain can only lower rho_V. Anything released after
   relaxation leaves rho_V < 0 and a crunch: EW + QCD plateaux -4.7e54 rho_Lambda (computed history ends in turnaround);
   an inflaton plateau V^(1/4) = 1e16 GeV over 60 e-folds drains -9.5e99 rho_Lambda; the trace drain of the hot plasma
   is -gamma kappa^2 H_RH^4/32 = -(1.5e93 to 1.5e95) rho_Lambda for T_RH = 1e16 GeV, kappa = 1e-3 to 1e-2.
7. **Thermal-history compatibility only if prepared** [D, C]. A standard radiation era with BBN (H age = 1/2,
   H = H_std(T)) is obtained only when rho_V,i is put in by hand at its final value; the dissipation is then inert
   (relative change 1e-86 after BBN, 8e-13 over the whole history in the run). Late-time effects are unobservable:
   FIRAS-type injection 7e-109 of rho_gamma per e-fold at z = 2e6 (and into gravitons, not photons), graviton
   (N_eff) fraction <= 1e-16 per e-fold, no G-variation at constant G.
8. **Verdict: QV2-diss** (§6). It does **not** remove the need to prepare mu in any history with a standard thermal era.
   It forms mu_eff only asymptotically, only from above, only in an empty vacuum-dominated universe, after ~1e123
   e-folds. The residual is fixed (positive, rho_V = 9M^4(6 gamma sqrt(3) M t)^(-2/3)) only on the attractor, where it
   is wrong by 1e82; otherwise it is initial-data and history dependent, and the set of initial data compatible with
   BBN and rho_V(t0) ~ rho_Lambda has relative measure ~rho_Lambda/rho_scale (2e-55 against the EW release alone).
   **Single next dependency:** a first-principles derivation of S[mu_eff] (one-loop matter effective action in the
   q-dependent background, with the 3-form constraint) that decides whether it contains a two-sided term linear in
   rho_V with rate >= H; the KV term's drain per e-fold is O(gamma R^2) ~ gamma H^4 and cannot form mu (§6.4).

---------------------------------------------------------------------------------------------------------------------

## 1. Conventions and the dissipative system in BRIEF.md notation

Signature (-,+,+,+), hbar = c = 1, K = M^2 = 1/(8 pi G), M = 2.435e18 GeV, R = 6(Hdot + 2H^2) in flat FRW (MTW sign).
KV use Weinberg's curvature sign (R_W = -R; literature §2.1); every equation below is translated. Only R^2 and
Hdot + 2H^2 enter the dissipative terms, so the sign convention cancels there [D].

**Generic dissipative q-system** (all three papers are special cases) [D]:

    (q-diss)  d/dt [ eps'(q) - (1/2) K'(q) R ] = S(t)          i.e.  d mu_eff/dt = S,  mu_eff := eps' - K'R/2
    (g)       K G_mn = T^m_mn - (eps - mu_eff q) g_mn - (g Box - nabla nabla) K      (unchanged, mu -> mu_eff(t))
    (m-diss)  nabla^m T^m_mn = q S u_n      (u^m the cosmic rest frame; u_n = (-1,0,0,0))
              <=>  rho_m-dot + 3H (rho_m + p_m) = q S
    rho_V := eps(q) - mu_eff q,   rho_V-dot = (1/2) R K-dot - q S

(m-diss) is not an extra assumption: the identity nabla^m E_mn = -q d_n mu(x) holds off shell for any metric and q(x)
(`../equations/NOTES.md` §1.4; FRW form re-verified symbolically, `results/symbolic_checks.json` check A). Once mu_eff
varies, the Einstein equation can hold only if matter absorbs q mu_eff-dot. Two consequences [D]:
- The modified system is **energy-conserving overall**: total rho_m + rho_V obeys the standard continuity equation,
  so the Bianchi identity is satisfied. The dissipated energy is not lost; it is in matter.
- The modification needs a **preferred frame** u^m (the rest frame of the produced quanta / the FRW frame), and S is
  defined only along it. KV give no covariant form; a covariant completion would be
  nabla_m mu_eff = -S u_m with S built from the expansion theta = nabla.u and R. This is outside the action: the
  q-equation is no longer a variational equation [A].

---------------------------------------------------------------------------------------------------------------------

## 2. Extraction and translation of the three papers

### 2.1 KV16a, arXiv:1601.00601 (constant G, particle-production backreaction)

Equations (paper numbers) and translation [Q, D]:
- (2.1) action as 0806.2805 with G(q); §4 sets **G = G_N** (eq 4.1). Then K' = 0, the q-equation reads
  eps'(q) = mu_eff(t), and rho_V = eps - q eps' (eq 4.4, "the equilibrium form", exact at K' = 0).
- (3.1) d mu_eff/dt = S; (4.6a) rho_V-dot = -q S; (4.6b) rho_M-dot + 3H(rho_M + P_M) = +q S. This matches (m-diss).
- (5.3) rho_M-dot + 4 H rho_M = (gamma/36)|H| R^2 =: S_M; (5.6) rho_V-dot = -S_M = -gamma |H| (Hdot + 2H^2)^2.
  **Radiation equation** (where the energy goes), in BRIEF notation: rho_r-dot + 4 H rho_r = (gamma/36) |H| R^2,
  with q S = (gamma/36)|H| R^2.
- (5.4) gamma ~ 1/(8 pi) (naive), (5.5) gamma = 1/(32 pi^2) = 3.1663e-3 for gravitons (App. A).
- Dimensionless system (6.1): f fdot eps'' = gamma |h| (hdot + 2h^2)^2, hdot + 2h^2 = 2 r_V (radiation only),
  eps = lambda + f^2 + 1/f^2, r_V = lambda - f^2 + 3/f^2 (6.2). **Unit identification [D]** (`s1` check B): (6.1)
  equals the dimensionful (5.6) with the same gamma only in units E_* = (3/(8 pi G))^(1/2) = sqrt(3) M = 4.22e18 GeV
  (r = rho/E_*^4, h = H/E_*, tau = E_* t). The text says "rescaling with powers of E_P"; with E_P (G = 1) units, as
  K22 eq (7c) uses, the same printed equations mean gamma_printed = (8 pi/3) gamma_phys. So "gamma = 1/100" in the
  figures is gamma_phys = 1/100 (E_* reading) or 1.19e-3 (E_P reading). This changes no conclusion; §5 uses the
  dimensionful (5.6) with gamma = 1/(32 pi^2), which is unambiguous.
- App. B [Q, D verified]: with xi = ln a, dr_V/dxi = -4 gamma r_V^2, dk/dxi = -2 gamma k^2, k = hdot + 2h^2;
  r_V = 1/(4 gamma (xi - xi0)), h^2 = (C1 + e^{4 xi0} Ei(4 xi - 4 xi0)) e^{-4 xi}/gamma. Verified symbolically (`s1`).
  Asymptotically r_M/r_V = 1/(4(xi - xi0)) + O(xi^-2) [D, from the Ei series], and xi = (3 tau/(4 sqrt(gamma)))^(2/3),
  which gives KV's (7.1) r_V = (6 gamma tau)^(-2/3), h = (6 gamma tau)^(-1/3) **exactly at leading order** [D].
- App. B, last paragraph [Q]: the attractor domain is "f(1) such that r_V(1) >= 0, and h(1) >= 0". That is the exact
  basin statement (§4.2); the paper's "finite domain" means this half-space, not a compact region.

**Derived versus phenomenological** [Q, A]:

| Element | Status |
|---|---|
| Massless non-conformal fields are produced at a rate ~R^2 (Zel'dovich-Starobinsky 1977; Dobado-Maroto gr-qc/9803076) | derived in QFT in curved spacetime, perturbatively in V = a''/a, as a **global** in-out Bogoliubov result for asymptotically static past and future (KV16a §5, first remark) |
| Energy per quantum ~abs(H) (eq 5.2) | assumed ("coefficient of unity") |
| Local rate from the global integral (App. A: Taylor expansion with Delta eta V' ~ V, keep one term of V', drop higher-derivative and nonlocal terms, then move the upper limit to t) | heuristic; O(1) uncertainty in the local coefficient |
| gamma = 1/(32 pi^2) | computed: the prefactor of the exact double integral (A.5) for gravitons; as a **local** coefficient it inherits the heuristic step |
| The produced energy is taken from rho_V (the q-field is the "agent", Unruh/Schwinger analogy) | phenomenological; KV: "the simplest possible ... other changes are certainly to be expected" |
| Use in de Sitter and on the attractor (H never -> 0 in the past) | outside the stated validity of the ZS result. In exact de Sitter the law gives rho_V-dot = -4 gamma H^5, i.e. it **assumes** de Sitter decays (KV16a and KSV16 say so). For free massive or conformal fields in the Bunch-Davies state <T_mn> is proportional to g_mn and there is no steady loss [R]; the massless minimally coupled case (gravitons) has infrared subtleties and the question is disputed (Polyakov vs Starobinsky-Yokoyama, cited by KSV16) [R] |
| G = G_N, eps = lambda + f^2 + 1/f^2 | choices (eps chosen so r_V spans R and chi > 0) |

### 2.2 KSV16, arXiv:1601.04676 (nonconstant G, condensed-matter dissipation function)

- Units of 0806.2805 with sigma = 1: K = f q0/3 (G^-1 = (8 pi/3) f q0), eps = (1/2chi0)(-q^2/q0^2 + q^4/3q0^4)
  (eq 5a), G^-1 = sqrt(q^2)/(G_N |q0|) (eq 7), mu0 = -1/(3 chi0 q0) (eq 6b), u0 = -1/3.
- (11) d mu_eff/dt = S; (12a,b) rho_V-dot = -q S, rho_M-dot = -3H(P + rho) + q S (constant-G illustration).
- (13) q S = Gamma_q qdot^2 + Gamma_H Hdot^2; the Gamma_H term "has a parallel with" the ZS R^2 term:
  (R^2 - 6E)/36 = Hdot^2 in flat FRW (E = Gauss-Bonnet) (eq 14) [Q]. **Pure ansatz**, modelled on
  condensed-matter dissipation functions; "we follow the lead of condensed matter physics".
- Dimensionless ODEs (17a-c): hdd + 4h hdot = fdot eps'' - s~; h fdot + f h^2 = r~_V + r_M; r_M-dot + 3h(1+w)r_M = f s~,
  with r~_V = eps - f eps' + f(hdot + 2h^2) = eps - u_eff f, u_eff = eps' - hdot - 2h^2, du_eff/dtau = s~ (eqs 18a,
  19-21), s~ = gamma_f fdot^2 + gamma_h hdot^2 (18b). **All verified from the BRIEF equations** with K = f/3 (`s1` check C). The
  de Sitter family (23): f, h constant, h^2 = eps' - eps/f, r_M = 0, exists for every u_eff [Q, D].
- Notational point [D]: eq (13) puts the quadratic form in q S, while (18b), (21) put it in s~ = du_eff/dtau and (17c)
  injects f s~ into matter. The two differ by a factor f = q/q0. The numerics follow (17)-(21).
- Energy: total conserved; the dissipated energy goes into the matter fluid with w_M (= 1/3 in all runs) [Q].
- Paper's result [Q]: no Minkowski attractor; de Sitter or collapse to a curvature singularity (App.: h ~ ln(tau_sing - tau),
  r_M ~ 1/(tau_sing - tau)); Minkowski only for initial data on a 3-dimensional separatrix in the 4-dimensional space
  {h, hdot, f, r_M}(1); "the decay of the de-Sitter vacuum is a necessary condition".

### 2.3 Reading of the two 2016 papers together

KV16a states that the KSV16 source "does not discriminate between de Sitter and Minkowski vacua" and that its own
source "does prefer the Minkowski vacuum" [Q, KV16a §3]. So the authors present them as two hypotheses about one
open question (does de Sitter radiate?), not as conflicting results. §4.1 makes this exact.

### 2.4 K22, arXiv:2207.03453 (4D-brane realization, mu in the action)

- Action (1): S = -∫sqrt(-g)[R/16 pi G_N + eps(n/sqrt(-g)) + L^M] + mu ∫ n. q = n/sqrt(-g) has mass dimension 4,
  **mu is a dimensionless parameter of the action**, not an integration constant [Q]. Varying n gives
  eps'(q) = mu (2a), an **algebraic** equation; rho_V = eps - q eps' = eps - mu q (4); the Einstein equation is
  standard GR with rho_V (3). Then d rho_V/dq = 0 (5) and, for homogeneous q, d rho_V/dt = 0 (6).
- Dissipative model (7a-d) / (15a-e): r_M-dot + 3h(1+w)r_M = +Gamma, r_V-dot = -Gamma, with
  Gamma = gamma|h|(hdot+2h^2)^2 (8, from KV16a) or the "ad hoc" Gamma = g~ |h| h^4 (17) used for all figures;
  xi-dot dr_V/dxi = -Gamma (15b); r_V(xi) = lambda + xi^2 - mu xi (13, 15e, 16a).
- **Finding [D]** (`s1` check D): at Gamma = 0, (15b) xi-dot r_V'(xi) = 0 admits every constant xi. So (15b) is a
  dissipative deformation of the weaker statement (6), not of the field equation (2a). The late-time endpoint xi_bar
  (r_V(xi_bar) = 0) violates (2a) by eps'(xi_bar) - mu = sqrt(mu^2 - 4 lambda) = 3.46 (lambda = 1) and 4.47
  (lambda = -1) at mu = 4, unless lambda = mu^2/4. This is the content of K22's "puzzle" (§III.E): the endpoint is not
  the equilibrium of the action. The brane model therefore does not form mu either; it abandons the n-equation and
  drains rho_V directly, exactly as KV16a at constant G ("exactly the same" ODEs for w = 1/3 and Gamma (8), K22 §III.B;
  verified up to the unit factor 8 pi/3 of §2.1).
- Energy: total conserved; Gamma goes into the w_M fluid [Q].
- The ad hoc Gamma does not vanish at r_V = 0 when radiation is present (it vanishes only at h = 0). With radiation
  domination it overshoots r_V below zero; §4.4 quantifies this.
- 2207.02826 §VII (section read) uses the same structure, Gamma = g~ |r-dot/r| r_vac^2 with a switch-on function, and
  needs gamma = 2e11 to show relaxation from 6e-6 to 1e-10 by tau = 100 [Q]: fourteen orders above the graviton value.

### 2.5 Energy conservation and Bianchi identity, summary [D, C]

| Model | Bianchi / energy | Where the energy goes | Preferred frame |
|---|---|---|---|
| KV16a | holds; rho_V-dot = -qS, rho_r-dot + 4H rho_r = +qS | massless quanta (gravitons), w = 1/3 | yes |
| KSV16 | holds (constraint propagates as C-dot = -3hC, `s1` C) | matter fluid w_M | yes |
| K22 | holds (GR + rho_V(t) + sourced matter) | matter fluid w_M | yes |
| QV2-diss (§6) | holds (generic identity, `s1` A) | massless quanta | yes |

---------------------------------------------------------------------------------------------------------------------

## 3. Reproduction of each paper's headline result [C] (`code/s2_reproduce.py`, `results/reproductions.json`)

| Paper, item | Parameters and data (paper) | Paper value [Q] | This track [C] |
|---|---|---|---|
| KV16a Fig. 1 | lambda = 1, gamma = 1/100, {h, f}(1) = {6, 3/5} | r_V(1), r_M(1) = 8.97333, 27.0267; r_V(1e5) = 3e-3 | 8.97333, 27.0267; **3.0285e-3** |
| KV16a Fig. 2 | lambda = -1 | 6.97333, 29.0267; 3e-3 | 6.97333, 29.0267; **3.0284e-3** |
| KV16a lambda = 0 | "similar" | | 3.0284e-3 |
| KV16a Fig. 3 | 8 corners {6 +- 1, 3/5 +- 1/5}, lambda = +-1 | all -> r_V = 0 | all 8: r_V(1e5) = 3.0284e-3 to 3.0285e-3; f -> f_bar(lambda) (1.51685 vs 1.51749 at lambda = 1; 1.14091 vs 1.14139 at -1) |
| KV16a App. B | exact r_V(xi) | | agreement 3.0e-12 (relative) |
| KV16a (7.1) | r_V ~ (6 gamma tau)^(-2/3) | | ratio 1 - 2.8e-5 at tau = 1e5; gamma = 1, 1e-3, 1/(32 pi^2) all within 5e-4 |
| KSV16 Fig. 1 | no dissipation, {h, hdot, f, r_M}(1) = {0.55, 0, 0.2953, 0} | u_eff = -0.883133 constant; de Sitter | -0.8831328 constant; h -> 0.795511, h^2 = eps' - eps/f to 2e-15 |
| KSV16 Fig. 2 | gamma_f = gamma_h = 1 | u_eff(100) = -0.333499; de Sitter | **-0.3334995**; de Sitter, h -> 0.01243 |
| KSV16 Fig. 3 | f(1) = 0.2952 | u_eff(98) = -0.333276; recollapse; a -> O(10) at tau ~ 198 | **-0.333276**; turnaround, singularity at tau = 197.58 with **a = 10.7** |
| KSV16 Table I | Y/N on 6 x 6 grid, three slices | pattern NNNYYY / NNNYYY / NYYYYY / YYYYYY x3 | **identical in all three slices** |
| KSV16 Table II | h(1) = 0.525, f(1) = 0.352 ... 0.39 | 0.352: none; 0.353: h, f, r_V, r_M = 0.0362, 1.00131, 0.00132, 5e-12; ... 0.39: 0.187, 1.0334, 0.0363 | 0.352: turnaround at tau = 68.6; 0.353: 0.0363, 1.00131, 0.00132, 3e-12; ... 0.39: 0.1874, 1.03344, 0.03631 |
| K22 Figs. 1-2 | mu = 4, lambda = +-1, g~ = 1e-3, {h, xi}(1) = {9, 4.55} | h ~ (g~ tau)^(-1/3), xi -> xi_bar | xi(1e6) = 3.7321 (xi_bar 3.7321), 4.2361 (4.2361) |
| K22 Table III | lambda = 4 | tau = 1e2...1e6: h = 0.926, 0.430, 0.200, 0.0927, 0.0430; xi = 2.32 ... 2.01; r_M = 1.85e-4 ... 8.56e-10; r_V = 0.102 ... 0.000221 | 0.926, 0.4299, 0.1996, 0.09266, 0.04301; 2.3197 ... 2.0149; 1.845e-4 ... 8.557e-10; 0.1022 ... 2.208e-4 |
| K22 rectangle (20) | h in [8.9, 9.1], xi in [4.54, 4.56] | attractor | all 8 corners -> r_V(1e5) = 1.0249e-3, xi -> xi_bar |

Exact asymptote of K22 [D] (`s1` D): r_V = (3/2)^(-2/3)(8 pi/3)^(-5/3)(g~ tau)^(-2/3) = 0.022084 (g~ tau)^(-2/3),
h = 0.430127 (g~ tau)^(-1/3), r_M = g~ h^4/4. This reproduces Table III to its 3 digits; K22's "~" in (18) hides the
coefficient 0.0221.

Integration: DOP853, rtol 1e-10 to 1e-11; KSV16 in a regular first-order form (`diss_models.py`) that is equivalent on
the constraint and does not divide by h (the paper's form does, which matters at turnaround). Constraint residuals
<= 1e-12 throughout.

---------------------------------------------------------------------------------------------------------------------

## 4. Why one paper finds an attractor and the other none; basin maps

### 4.1 Structural statement [D]

Let S >= 0 (both papers). Then:
- KV16a: rho_V-dot = -qS <= 0 (q > 0 branch). rho_V is a Lyapunov-type function: **non-increasing on every
  trajectory**, expanding or contracting (|H|).
- KSV16: u_eff-dot = s~ >= 0: **u_eff is non-decreasing**.
- A trajectory can end only where S = 0 persistently. KV16a: S = 0 iff H = 0 or R = 0. With radiation as the only
  matter, R = 4 rho_V/K exactly, so the only expanding end state is rho_V = 0 (Minkowski). De Sitter (R = 12H^2) is
  not a fixed point. KSV16: s~ = 0 iff fdot = hdot = 0, i.e. **every static state**: Minkowski and the whole de Sitter
  family, one member per value of u_eff. A trajectory freezes when the motion damps out, at whatever u_eff it has
  reached.
- Therefore KV16a has an attractor whose basin is "everything that starts on the correct side", and KSV16 has a
  continuum of attractors (de Sitter) with Minkowski as the boundary between "froze too early" (de Sitter,
  u_eff < u0) and "overshot" (u_eff > u0: anti-de Sitter-like, collapse). Both are **one-sided**: a monotone quantity
  cannot return. The two papers do not conflict; they differ in one assumption, whether de Sitter radiates
  (KV16a yes, KSV16 no), which KSV16 state as the conclusion ("decay of de Sitter is a necessary condition").

### 4.2 KV16a basin (`code/s3_basin.py`, `results/basin_maps.json`)

**Exact statement [D].** Radiation only: dr_V/dxi = -4 gamma r_V^2 and h^2 = r_V + r_r >= 0. If r_V,i > 0 and h_i > 0,
r_V stays positive and decreases as 1/(4 gamma xi): Minkowski is reached. If r_V,i < 0, r_V decreases further and
h^2 = r_V + r_r reaches zero as r_r dilutes: turnaround and crunch. r_V,i = 0 stays (measure zero). In terms of the
paper's variables, for any lambda the basin is f_i <= f_*(lambda) = sqrt((lambda + sqrt(lambda^2 + 12))/2), H_i > 0.

**Grid [C].** lambda in [-1, 1] step 0.05 (41 values; units E_*^4 = 9 M^4; since the boundary is r_V,i = 0 the picture
is the same for lambda in [-1, 1] E_P^4 up to the f-scale) x f_i in [0.20, 3.00] step 0.05 (57); r_r,i = 10;
dust r_d,i in {0, 10, 100, 1000}; gamma = 1/(32 pi^2) and 0.01. Rigorous classification by monotonicity (crunch once
r_V < 0; Minkowski once r_V > 0 and r_d < 1e-6 r_V, after which the remaining dust drain is < 1e-6 r_V).

| gamma | r_d,i | basin fraction of the box | basin boundary in r_V,i (bisection, lambda-independent) | estimate gamma r_d^2/24 [D] |
|---|---|---|---|---|
| 1/(32 pi^2) | 0 | **0.4018** | 0 (bisection f vs analytic f_*: agree to 1e-9 at lambda = -1, -0.5, 0, 0.5, 1) | 0 |
| 1/(32 pi^2) | 10 | 0.4014 | 0.01329 | 0.0132 |
| 1/(32 pi^2) | 100 | 0.3188 | 1.419 | 1.32 |
| 1/(32 pi^2) | 1000 | 0 | above the box | 132 |
| 0.01 | 0 | 0.4018 | 0 | 0 |
| 0.01 | 10 | 0.3992 | 0.04261 | 0.0417 |
| 0.01 | 100 | 0.1904 | 5.368 | 4.17 (nonlinear regime) |
| 0.01 | 1000 | 0 | above the box | 417 |

Convergence [C]: radiation-only grids agree with the analytic prediction at all 2337/2337 points; rtol 1e-9 vs 1e-11
changes no classification on any grid; no point is undecided at tau_max (1e3 radiation only, 1e6 with dust). One
bisection point (lambda = 1, gamma_phys, r_d = 10) landed within the tau_max resolution of the boundary and is reported
as undecided there.

Dust (any trace, rho - 3p > 0) moves the fixed point of the drain to rho_V = -(rho_m - 3p_m)/4 < 0, so with enough
non-conformal matter early the drain overshoots and the trajectory crunches. Leading order [D]: the boundary is
r_V,i ~ gamma r_d,i^2/24 when dust dilutes first (`s1` E: matter-era drain gamma H_1^4/24).

### 4.3 KSV16 outcome map and separatrix

Grid [C]: lambda in {-0.12, -0.1, -0.05, 0, 0.1, ..., 1.0} (eps + lambda; a Minkowski equilibrium with f > 0 exists only
for lambda >= -1/8 [D]) x f_i in [0.10, 1.50] step 0.05; h_i = 0.55, hdot_i = 0, r_M,i = 0 (paper slice);
gamma_f = gamma_h = 1; tau_max = 1500. Outcomes: de Sitter 33.0%, collapse 67.0%, undecided 0; rtol 1e-9 vs 1e-11:
no change. No grid point ends near Minkowski: the smallest final r~_V among de Sitter points is 6.4e-3, the median 0.38.
(h_i, f_i) map at lambda = 0 (h in [0.45, 0.90] step 0.025, f in [0.10, 0.60] step 0.025): de Sitter 80%, collapse 20%.

Separatrices [C] (bisection to 1e-10 in f_i at every lambda; de Sitter side sampled at delta f = 1e-3, 1e-4, 1e-5
with tau_max = 1500 and 6000, which agree to the printed digits):
- For **every** lambda in [-0.12, 1] at least one collapse/de Sitter boundary is a Minkowski separatrix: approaching
  it from the de Sitter side, u_eff,inf -> u0(lambda) from below and r~_V,inf -> 0. Example lambda = 0: f_sep = 0.295245
  (KSV's Table I/II boundary) and 1.286909; lambda = 1: f_sep = 1.400694.
- Fine-tuning law: r~_V,inf = c (delta f)^p with local p = 0.77-0.85 between delta f = 1e-4 and 1e-5 (c ~ 1 to 10).
  Examples: lambda = 0, f_sep = 0.295245: r_V = 1.69e-3, 2.51e-4, 4.29e-5; lambda = 1: 2.16e-2, 2.98e-3, 4.46e-4.
  So |rho_V,inf| <= epsilon needs |delta f| <~ (epsilon/c)^(1/p), which for epsilon = rho_Lambda in natural units
  (~1e-121) is ~1e-140 or smaller: a **measure-zero condition, tuned more finely than epsilon itself**.
- For lambda < 0 a second boundary at small f (towards the q = 0, G^-1 = 0 vacuum) is **not** through Minkowski:
  r~_V on its de Sitter side stays at 0.22 (lambda = -0.12), 0.17 (-0.10), 6.8e-4 (-0.05).
- With lambda the de Sitter end states carry r~_V from 6e-3 to O(1) across the grid. KSV's conclusion is confirmed and
  extended to lambda in [-0.12, 1]: no attractor, Minkowski on a codimension-1 set.

### 4.4 K22 brane basin

Grid [C]: mu = 4; lambda in [-1, 5] step 0.25 x xi_i in [0.1, 6.0] step 0.1; h_i = 9; r_M,i from Friedmann (points
with r_M,i < 0 invalid); tau_max = 1e5 (undecided points re-run to 1e6: none).

| Gamma | basin fraction (valid points) | lambda <= mu^2/4 = 4 | lambda > 4 |
|---|---|---|---|
| ad hoc g~abs(h)h^4, g~ = 1e-3 | 0.309 | 0.363 | 0 (all 'singular': xi reaches mu/2 with r_V > 0 and xi-dot = -Gamma/(2xi - mu) diverges) |
| particle production (8), gamma = 1e-3 | 0.392 | 0.460 | 0 (singular) |

The ad hoc Gamma loses basin because it does not vanish at r_V = 0 while radiation is present: e.g. lambda = 1,
xi_i = 3.8 (r_V,i = 0.24 > 0, h = 9) crunches at tau = 1.005. For lambda > mu^2/4 no zero of r_V exists and the model
has no continuation; KV's condition (16b) lambda <= mu^2/4 is a hard boundary [D, C].

### 4.5 With the QV1 structure K'(q) != 0

Our transcription [D] (KV16a §7, second remark, states the result for nonconstant G without printing equations):
units of 0806.2805/KSV16 (sigma = 1, chi0 q0^2 = 1), K = f/3, eps = lambda + (1/2)(-f^2 + f^4/3), radiation, and the
drain du_eff/dtau = s with f s = gamma |h| (hdot + 2h^2)^2 (energy into radiation). The q-mode is massive,
m_q^2 = 2 f0 eps''(f0) (m_q = sqrt(2) at lambda = 0), so the hierarchy is compressed to m_q/H = 2.4 to 53
(`code/s5_qv1diss.py`, `results/qv1diss.json`).

1. **Reduction to the constant-G law** [C]. Constant-G prediction with K0 = f0/3: k = 2 r_V/f0, so
   r_V(N) = r_V1/(1 + 4 gamma r_V1 N/f0^2). Runs start on the static de Sitter branch at u = u0 + du, gamma = 0.01,
   6 e-folds:

   | lambda | du | H_i/m_q | max relative deviation from the law over 6 e-folds |
   |---|---|---|---|
   | 0 | -0.3 | 0.41 | 1.1e-2 |
   | 0 | -0.1 | 0.23 | 7.3e-4 |
   | 0 | -0.03 | 0.12 | 2.5e-5 |
   | 0 | -0.01 | 0.071 | 1.0e-6 |
   | 0 | -0.003 | 0.039 | 2.8e-8 |
   | 1 | -0.3 ... -0.003 | 0.19 ... 0.019 | 6.2e-4 ... 2.3e-9 |

   The deviation falls roughly as (H/m_q)^6 (local slope 5.8-5.9) [C]. At the physical hierarchy (m_q/H ~ 1e32 at the
   EW scale, `../equations/NOTES.md` §7) the constant-G law of §5.1 holds to any accuracy that matters. QV1's K'(q)
   does not change any conclusion of §5.
2. **KV's nonconstant-G remark** [C]. From {h, hdot, f, r_M}(1) = {1/2, 0, 1/2, 1/20} (0806.2805 data, u_eff free):
   tau^(2/3) r_V = 6.54 (lambda = 0) and 11.64 (lambda = 1) at tau = 1e5, against the constant-G attractor prediction
   f0^(5/3)(6 gamma)^(-2/3) = 6.53 and 11.64; tau^(1/3) h = 2.554 and 2.867 against f0^(1/3)(6 gamma)^(-1/3) = 2.554 and
   2.866; f -> f0 (1.003, 1.4155 vs 1, 1.4142). So the remark is confirmed for lambda = 0, +1, and the K' != 0 attractor
   is the constant-G attractor with K0 = K(q0).
3. **lambda range** [D]: for this eps, r_V(f) = lambda + f^2/2 - f^4/2 has maximum lambda + 1/8, so no lambda < -1/8 can
   be cancelled. KV16a's "arbitrary lambda" for the nonconstant-G model needs a different eps (not printed).

---------------------------------------------------------------------------------------------------------------------

## 5. Late-time behaviour and compatibility with a standard thermal history (`code/s4_history.py`, `results/history.json`)

### 5.1 The relaxation law per e-fold [D]

In BRIEF conventions at constant K (or K' != 0 with m_q >> H, §4.5), for H > 0:

    d rho_V / dN = -(gamma/36) R^2,     R = (rho_m - 3 p_m + 4 rho_V)/M^2,     N = ln a.

- Radiation only (exact): rho_V(N) = rho_Vi / (1 + (4 gamma/9)(rho_Vi/M^4) N). The amount of radiation does not enter.
  Hence rho_V(N) >= (1/2) min(rho_Vi, 9M^4/(4 gamma N)): **a vacuum energy is reduced by a factor 2 or more only if
  (4 gamma/9)(rho_Vi/M^4) N >= 1**.
- Reaching rho_Lambda needs N_Lambda = 9M^4/(4 gamma rho_Lambda) = **9.9e122 e-folds** [C].
- A standard history from T = M to T_CMB has N_std = ln(M (106.75/3.91)^(1/3)/T_CMB) = **72.5** e-folds [C]. After it,
  the floor is 9M^4/(4 gamma N_std) = 3.4e74 GeV^4 = 9.8 M^4 = 1.4e121 rho_Lambda.
- Non-conformal matter: the drain stops at R = 0, i.e. rho_V = -(rho_m - 3p_m)/4 < 0. Its contribution per e-fold
  is (gamma/4) kappa^2 H^4 for a plasma with trace fraction kappa; integrated over a radiation era from T_RH it is
  D = gamma kappa^2 H_RH^4/32 [D]; over a matter era from H_1 it is gamma H_1^4/24 [D].

### 5.2 On the attractor at the observed age [D, C]

Vacuum plus radiation, vacuum dominated (KV16a (7.1) translated):
rho_V = 9M^4 (6 gamma sqrt(3) M t)^(-2/3), H = sqrt(3) M (6 gamma sqrt(3) M t)^(-1/3),
N = (3 sqrt(3) M t/(4 sqrt(gamma)))^(2/3), Omega_M = 1/(4N), and t(H) = M^2/(2 gamma H^3).

| Quantity | Attractor value at t0 = 13.787 Gyr | Observed [R], retrospective |
|---|---|---|
| rho_V(t0) | 2.24e35 GeV^4 = **8.9e81 rho_Lambda** | rho_Lambda = 2.52e-47 GeV^4 |
| H(t0) | 0.112 GeV = 7.8e40 H0 | H0 = 1.44e-42 GeV |
| H(t0) t0 | 7.4e40 | 0.950 (H0 t0 from H0 = 67.4, t0 = 13.787 Gyr; LCDM integral 0.951) |
| e-folds since tau = 1 | 1.1e41 | |
| Omega_M | 2.2e-42 (produced gravitons only; no photons, no CMB) | Omega_m = 0.315 |
| age when H = H0 | 4.5e122 / H0 = 6.6e123 Gyr, with Omega_V = 1 - O(1e-123) | |

Verdict [C]: a universe on the KV attractor has neither the observed H0-age relation nor a matter or radiation
component. **Retrospective failure** of the attractor history. KV16a §8 present the attractor as an early,
inflation-like phase, "possibly" followed by a standard FRW phase [Q]; §5.4 examines that reading.

### 5.3 Does the q-sector become subdominant? [C]

Integrations in N from T = 1e18 GeV with SM radiation (g* steps), EW and QCD plateaux (Delta_EW = 1.19e8 GeV^4,
Delta_QCD = 4.4e-4 GeV^4) released into radiation, dust (Omega_m) appearing at T = 100 GeV, and the KV16a drain with
gamma = 1/(32 pi^2):

| Initial rho_V (q-sector) | Outcome |
|---|---|
| M^4 (kappa = 0 or 1e-2) | vacuum dominated throughout; rho_V falls only by 9% in 73.5 e-folds (to 1.26e120 rho_Lambda); the produced quanta keep the radiation "temperature" above 1.3e17 GeV, so no epoch with T = 1 MeV is ever reached |
| 1e-60 M^4 = 3.5e13 GeV^4 | radiation era down to ~1 TeV, then vacuum dominated: at T = 1 MeV, H age = 14.8 and H/H_std = 3.2e12 (BBN needs 1/2 and 1); rho_V unchanged |
| rho_Lambda (prepared) | standard history: H age = 0.500 and H/H_std = 1.000 at T = 1e16 GeV, 1 TeV, 1 GeV, 1 MeV; final total vacuum = 0.9999999999992 rho_Lambda (the drain changed it by 8e-13) |
| -(Delta_EW + Delta_QCD) ("relaxed to zero total before the EW crossover") | the plateau release leaves -4.7e54 rho_Lambda; turnaround at T = 49 GeV, crunch |
| rho_Lambda with kappa = 1e-3 from 1e18 GeV | the trace drain gives -3.9e62 GeV^4 by N = 4.6; turnaround, crunch |

So the q-sector is subdominant at BBN **only if rho_V was already below ~rho_r(1 MeV) = 3.5e-12 GeV^4 by preparation**;
the relative change the drain can then make after BBN is 1.0e-86 [C].

### 5.4 The attractor as an early phase followed by reheating [D, C]

The final vacuum energy would be rho_V,end - Delta_EW - Delta_QCD - D_trace - D_inflaton, with
rho_V,end = 9M^4/(4 gamma N_end):
- To end at rho_Lambda after the EW and QCD releases alone: N_end = 2.1e68 e-folds, tuned to a fractional precision
  rho_Lambda/(Delta_EW + Delta_QCD) = **2.1e-55**: the preparation of mu is replaced by a preparation of N_end.
- At that point H = 2.6e-15 GeV and the radiation is 1/(4 N_end) = 1.2e-69 of the total (1.4e-61 GeV^4); a hot plasma at
  T = 160 GeV needs 2.3e10 GeV^4. The model has no reheating source: any reservoir that survives ~1e68 e-folds must be
  vacuum-like, and a vacuum-like reservoir is drained with the q-vacuum, because the drain sees only R (total trace).
  Example [D]: an inflaton plateau V >> rho_V drains Delta rho_V = -(4 gamma/9)(V^2/M^4) N_inf =
  -2.4e53 GeV^4 = **-9.5e99 rho_Lambda** for V^(1/4) = 1e16 GeV, N_inf = 60, which is then left behind as a negative
  vacuum energy when the inflaton decays.
- Trace drain of a hot plasma (kappa = rho - 3p over rho, [E] 1e-3 to 1e-2 from running couplings and thresholds):

| T_RH | kappa = 1e-3 | kappa = 1e-2 |
|---|---|---|
| 1e16 GeV | -1.5e93 rho_Lambda | -1.5e95 rho_Lambda |
| 1e9 GeV | -1.5e37 rho_Lambda | -1.5e39 rho_Lambda |
| 1e3 GeV | -1.5e-11 rho_Lambda | -1.5e-9 rho_Lambda |
| 160 GeV | -6.6e-18 rho_Lambda | -6.6e-16 rho_Lambda |

  QCD epoch (kappa = 0.3, T = 0.2 GeV): -1.2e-36 rho_Lambda. Matter era from equality: -1.4e-104 rho_Lambda.
  Note: the drain carries the KV/AHHKM scale H^4 (H(1.26 TeV)^4 = 0.98 rho_Lambda), but multiplied by
  gamma kappa^2/32 = 1e-10 to 1e-8 and with **negative** sign. It does not reproduce the KV electroweak "kick" (which needs the
  separate phenomenological freeze law of 0905.1919 eq (4.1); literature §5).
- Every term is negative or zero after relaxation. The residual is therefore set by N_end and the whole later
  history, with either sign, and is not fixed by the dynamics.

### 5.5 Late-time observables [C]

| Channel | Size | Bound (ledger row) | Verdict |
|---|---|---|---|
| Energy injection today | 1.6e-119 rho_gamma per e-fold; w + 1 = 4e-124 | FIRAS (row 15) | unobservable |
| Injection at z = 2e6 (mu-era) | 6.7e-109 rho_gamma per e-fold, and into gravitons: photons are conformally coupled and not produced by conformally flat FRW at tree level [R, Parker 1969] | FIRAS abs(mu) < 9e-5 | unobservable |
| Graviton energy fraction (N_eff) | 8.8e-94 per e-fold at BBN; 8.8e-17 per e-fold at T = 1e16 GeV (kappa = 1e-2) | N_nu = 2.898 +- 0.141 (row 12) | unobservable |
| G variation | none at constant G; with K' != 0, dK/K ~ H^2/m_q^2 = 1.6e-121 today as in QV1 | rows 1-5 | unobservable |

Late-time dissipation does not touch thermal-history observables. Its cosmological consequences are confined to the
early universe, where §§5.2-5.4 exclude the formed route.

---------------------------------------------------------------------------------------------------------------------

## 6. Verdict: QV2-diss

### 6.1 Definition

**QV2-diss** (a QV2 candidate in the sense of BRIEF_STAGE2: an extension of the equations, not a state condition, so
not a QV1 variant). Credited to Klinkhamer and Volovik (1601.00601) for every element:
- action, matter and generic eps(q), K(q) as QV1; **P_q dropped**;
- q-equation replaced by d mu_eff/dt = S along the cosmic frame, with q S = (gamma/36)|H| R^2, gamma = 1/(32 pi^2);
- matter equation nabla^m T^m_mn = q S u_n (energy into massless quanta); Einstein equation unchanged with mu_eff(t).
QV1 and its stage-1 results (`../equations`, `../ledger`) are kept beside it unchanged.

### 6.2 Does it remove the need to prepare mu?

| Question | Answer |
|---|---|
| Is mu_eff formed by the dynamics? | Yes, asymptotically: mu_eff -> mu0 of the eps present during relaxation, from every initial datum with rho_V,i > 0 (finite measure, a half-space; 40% of the stated box). Not from below. |
| How fast? | rho_V = 9M^4/(4 gamma N): 9.9e122 e-folds to reach rho_Lambda; 72.5 available in a standard history. |
| Compatible with a radiation era, BBN, H0, T_CMB, t0? | Only if rho_V,i is prepared at its final value; then the dissipation is inert (relative effect 1e-86 after BBN). On the attractor: rho_V(t0) = 8.9e81 rho_Lambda, H t0 = 7.4e40, Omega_m = 2e-42: retrospective failure. |
| Stable under later thresholds (gate Q3)? | No. One-sided: EW + QCD releases after relaxation give -4.7e54 rho_Lambda (crunch); the trace drain and any inflaton plateau add further negative terms. |
| Residual: fixed or initial-data dependent? | Fixed only on the attractor (positive, scale 9M^4(6 gamma sqrt(3) M t)^(-2/3) ~ gamma^(-2/3) M^(10/3) t^(-2/3)), and wrong there by 1e82. In any history with a thermal era it is initial-data and history dependent, of either sign. |
| Measure of initial data giving BBN and rho_V(t0) ~ rho_Lambda | relative measure ~ rho_Lambda/rho_scale of the initial range: 2e-55 against the EW release, 7e-121 against a Planck-scale range. Same tuning as QV1-cal, moved from mu to the initial mu_eff (or to N_end). |

**Verdict.** QV2-diss forms mu only in the formal limit of an empty, vacuum-dominated universe, from one side, after
~1e123 e-folds. In every history with a standard thermal era, mu (equivalently the initial mu_eff, or the duration of
the attractor phase) must still be prepared, to a precision of 1e-55 or better. **The need to prepare mu is not
removed.** The protect-or-form dichotomy (`../ledger` §1.5) sharpens to: QV1 protects mu and needs it prepared;
QV2-diss forms it, but too slowly and one-sidedly to matter, so it needs it prepared too.

### 6.3 Four-question gate for QV2-diss

1. Extra structure: a non-variational source in the mu_eff equation, a preferred frame, and massless quanta that
   receive the energy (KV16a).
2. Residual: removes a large positive term only asymptotically; fixes a residual (positive, (gamma t)^(-2/3)-scaling)
   only on the attractor, where it fails by 82 orders. Otherwise no.
3. Stable without retuning: no (one-sided; post-relaxation thresholds leave -Sum Delta and a crunch).
4. Independent consequence: none observable at late times (§5.5). The attractor's consequences (H t >> 1, Omega_m -> 0)
   are already excluded by exposed data (retrospective). No held-out quantity.

### 6.4 Single next mathematical dependency

**A first-principles form of S[mu_eff].** Derive the dissipative term in the mu_eff equation from the one-loop effective
action of quantum matter on the q-dependent background (4-form plus metric, with the Gauss constraint), and determine
(i) whether particle-production backreaction acts on mu_eff at all (rather than only through <T_mn> on the geometry),
and (ii) whether it contains a **two-sided term linear in rho_V with rate >= H** during the radiation era. Without such a
term no dissipative q-model can form mu: a drain that is quadratic in curvature, O(gamma R^2) ~ gamma H^4 per e-fold
(the KV16a and K22 forms), changes rho_V by less than (4 gamma/9)(rho_V/M^4) per e-fold (§5.1), and a dissipation
function that vanishes on static states (KSV16) leaves Minkowski on a separatrix. With such a term, the open question
becomes why it would stop at rho_Lambda rather than at zero.

GRUT contribution [A]: none beyond the protect-or-form bookkeeping already recorded in `../ledger` §1.5 and its
quantitative completion here. Nothing in this track is GRUT-native.

---------------------------------------------------------------------------------------------------------------------

## 7. Inputs and provenance

| Input | Value | Provenance | Role |
|---|---|---|---|
| M | 2.435e18 GeV | BRIEF.md | units, K = M^2 |
| rho_Lambda,obs | 2.52e-47 GeV^4 | BRIEF.md (Planck 2018; exposed) | comparison only, retrospective; never calibrated |
| gamma | 1/(32 pi^2) = 3.1663e-3 | KV16a eq (5.5), App. A [Q] (heuristic local reduction, §2.1) | drain coefficient |
| gamma = 1/100; eps = lambda + f^2 + 1/f^2; lambda = +-1, 0; {h, f}(1) | | KV16a figures [Q] | reproduction, basin |
| gamma_f = gamma_h = 1; eps of 0806.2805; boundary data | | KSV16 figures and tables [Q] | reproduction, basin |
| mu = 4, g~ = 1e-3, {h, xi}(1) = {9, 4.55} | | K22 figures [Q] | reproduction, basin |
| H0 = 67.4 km/s/Mpc, Omega_m = 0.315, t0 = 13.787 Gyr, T_CMB = 2.7255 K, Omega_gamma h^2 = 2.473e-5, N_eff = 3.046 | | Planck 2018 [R]; exposed | §5.2, §5.5, retrospective |
| g* steps (106.75, 61.75, 17.25, 10.75, 3.36), g*s0 = 3.91 | | [R] | radiation era |
| Delta_EW = 1.19e8 GeV^4; Delta_QCD = 4.4e-4 (to 3.0e-3) GeV^4 | | `../equations/NOTES.md` §7 | plateaux |
| kappa (plasma trace fraction) | 1e-3, 1e-2 (hot), 0.3 (QCD) | [E] | trace drain |
| V^(1/4) = 1e16 GeV, N_inf = 60 | | supplied [P], illustrative | inflaton drain |
| T_nr = 100 GeV (dust switch-on), z_eq = 3400 | | [P] / [R] | dust drain |
| FIRAS abs(mu) < 9e-5, abs(y) < 1.5e-5; N_nu = 2.898 +- 0.141 | | `../ledger/results/ledger.json` rows 15, 12 [Q] | §5.5 |
| hbar = 6.582e-25 GeV s, 1 Gyr = 3.156e16 s, 1 Mpc = 3.0857e19 km | | [R] | units |

Calibrated quantities: **none**.

---------------------------------------------------------------------------------------------------------------------

## 8. Limits

- The dissipative laws are the papers' (phenomenological or heuristic, §2). This track tests their consequences; it
  does not derive S. The §6.4 dependency is exactly this gap.
- Classical homogeneous FRW only. No perturbations, no spatial variation of mu_eff (KSV16 mention local catastrophes
  with mu fixed at infinity), no membrane nucleation.
- QV2-diss at K' != 0 is tested only at compressed hierarchies (m_q/H = 2.4 to 53, §4.5); the reduction to the
  constant-G law is controlled by H/m_q, which is ~1e-32 at the EW scale, so the extrapolation is safe for the
  smooth drain. Oscillation-driven particle production (parametric resonance at m_q) is outside KV16a's |H| R^2 law
  and is not modelled; it would add a further one-way drain [A].
- §5.3 histories use step g*, smoothed plateaux, a constant trace fraction and dust switched on at 100 GeV. These set
  orders of magnitude only; the verdicts rest on the exact radiation law (§5.1) and on monotonicity, not on these
  details.
- KSV16 classification uses tau_max = 1500 (paper: 200); near the separatrix the de Sitter endpoint converges slowly,
  so the fine-tuning exponents in §4.3 carry that resolution.
- The ZS/KV rate is applied in de Sitter-like phases outside its stated domain of validity; this is KV's assumption,
  kept here to test their scenario on its own terms.
- Read levels: the three main papers full; 2207.02826 section level; observational inputs recalled from Planck 2018
  and from the ledger track.

## 9. Files and reproduction

- `code/diss_models.py`: right-hand sides (KV16a with dust extension; KSV16 regular form, with the 'pp' source for
  K' != 0; K22 brane) and the integrator wrapper.
- `code/s1_symbolic.py` -> `results/symbolic_checks.json`: identity, unit map, App. B, KSV16 equations and constraint
  propagation, K22 endpoint violation and asymptote, per-e-fold laws (14/14 checks pass).
- `code/s2_reproduce.py` -> `results/reproductions.json`: §3.
- `code/s3_basin.py` -> `results/basin_maps.json`: §4.2-4.4 (about 10 min on 4 cores).
- `code/s4_history.py` -> `results/history.json`: §5.
- `code/s5_qv1diss.py` -> `results/qv1diss.json`: §4.5.

Reproduction (numpy, scipy, sympy, mpmath):

    cd research/reality_framework/grut_q_vacuum_construction_20261021/tracks/dissipation/code
    python3 -I s1_symbolic.py && python3 -I s2_reproduce.py && python3 -I s3_basin.py \
      && python3 -I s4_history.py && python3 -I s5_qv1diss.py
