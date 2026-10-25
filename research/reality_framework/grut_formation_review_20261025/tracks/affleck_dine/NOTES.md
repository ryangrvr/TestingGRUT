# Track affleck_dine: can one Z_6-allowed X^6 coefficient both generate and preserve the X charge?

25 October 2026 · Claude (independent review track) · reads Codex 8986d34 (origin/codex/grut-cstar-audit-20261024)
through `git archive` only. Own code in `code/`, outputs in `results/`. No commits.

**This is a conditional calculation inside one model with supplied fields. It is not a GRUT prediction.**
Nothing here is evidence that X, phi or any inflationary input exists.

Claim labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## 1. Question and the supplied inputs

The hypothesis is that the U(1)_X-breaking operator dV = c_d X^d / M^(d-4) + h.c. could generate the
primordial X asymmetry by the Affleck–Dine (AD) mechanism, as well as threaten late charge survival.
We test it for d = 6 (first allowed by a residual gauge Z_6), with d = 5 (Z_5) for comparison.

| Input | Value | Source |
|---|---|---|
| g_X, lambda_X | 0.1, 1e-3 | Candidate 1 / C1-HP-03/v1 [Q] |
| lambda_phi, v_phi (formation point) | 0.03, 0.924089637375362 GeV (Q=200000, omega/m_X=0.86119513) | 59a1a71/ff49e5d/8986d34 [Q] |
| Old benchmark | v=1 GeV, lambda_phi=0.01, Q=39944 | [Q] |
| M (operator scale) | reduced Planck, 2.435e18 GeV | 8986d34 [Q] |
| Required net yield | Y_Q = n_Q/s_vis = 4.8294009076988175e-9 (calibrated to B_DM = 4.4e-10 GeV) | portal-debris packet [Q] |
| Late survival (1 % core-charge loss in 13.8 Gyr) | c_6 <= 7.654e14, c_5 <= 1.2227e-4 (full response, R=40, tol 2e-8) | 8986d34 CHARGE_LEAKAGE_RESULTS.json [Q] |
| g* during AD | 106.75 | assumed |
| H_inf maximum | 4.7e13 GeV (r<0.036, P_R=2.1e-9) | [E] from BICEP/Keck 2021 [R] |
| Isocurvature limit | beta_iso < 0.038 (uncorrelated CDM, n_iso=1) → S_max = 9.1e-6 | Planck 2018 X, search-level [Q] |
| CMB-scale e-folds before end of inflation | N = 50 | assumed |
| Hubble-induced mass | -c_H H^2 abs(X)^2, c_H in {1, 10, 100}; switched off at release | supplied (scanned) |
| Release epoch | end of inflation (W1/W2), reheating (W3), or a frozen light field (W4) | supplied |
| Initial radial velocity | Xdot/X = -beta H_rel, beta in {0, 1.5} | supplied |
| Reheating | T_R in {1e9, 1e12, 1e14 GeV} or instant | supplied |
| CP angle | phi0 = d·theta0 + arg(c_d), scanned pi/8..pi; set by a mismatch between an inflaton-dependent and the low-energy A-term (DRT-type) | supplied |
| Excluded | stabilizing higher operators, thermal masses (hidden sector not in equilibrium with the SM at kappa ≈ 9e-10), all other breaking operators | assumption |

## 2. Method [D, C]

During the AD epoch abs(X) >> v, so phi has mass sqrt(2) g abs(X) >> H, omega and sits at phi = 0 [A]. X then
has V = lambda_X abs(X)^4 + 2 abs(lambda_d) R^d cos(phi), with phi = d·theta + arg c_d and lambda_d = c_d/M^(d-4).
With chi = aX, conformal time and a = 1 at release, the dimensionless equation solved by
`code/ad_condensate.py` (DOP853, rtol 1e-10) is

  u'' = (a''/a) u − abs(u)^2 u − (d/2) eta_d e^{−i arg c} conj(u)^(d−1) a^(4−d),   eta_d = (c_d/lambda_X)(X_0/M)^(d−4).

n a^3 = q X_0^2 omega_0 with q = −2 Im(conj(u) u'), omega_0 = sqrt(2 lambda_X) X_0 and h = H_rel/omega_0. The
yields follow [D]:

- release into radiation: **Y = q h^−3 (H_rel/M)^(3/2) / (2 lambda_X D)**, with D = (2π²/45) g* (90/π² g*)^(3/4) = 7.40;
- release into inflaton matter domination: **Y = q h^−3 T_R H_rel / (8 lambda_X M²)**.

For a Hubble-held field, X_0 = sqrt(c_H/2lambda_X) H_rel, h = 1/sqrt(c_H). In the linear regime q ≈ A_d eta sin(phi0), so for
d = 6 and release into radiation, **Y ≈ A_6 sin(phi0) c_6 c_H^(5/2) (H_rel/M)^(7/2) / (4 lambda_X³ D)**. Matter era:
Y ≈ A_6 sin(phi0) c_6 c_H^(5/2) T_R H_rel³/(16 lambda_X³ M⁴).

Linear coefficients A_d = q/(eta sin phi0) [C]:

| Release | A_6 | A_5 |
|---|---:|---:|
| radiation, c_H = 1, beta = 0 | 12.97 | 7.89 |
| radiation, c_H = 1, beta = 1.5 | 2.20 | 0.352 |
| radiation, c_H = 10 / 100 | 8.66 / 17.5 | 2.54 / 0.891 |
| matter, c_H = 1 / 10 / 100 | 13.0 / 6.47 / 11.9 | 11.1 / 2.84 / 0.912 |
| frozen field, radiation (q h^−3/2) / matter (q h^−2) | 8.47 / 7.96 | 6.07 / 7.92 |

The charge fraction is epsilon ≈ 4.5 eta in the linear regime. It saturates at 0.8–0.9 for eta ≈ 0.3 at the best
phase. **Boundedness:** for d = 6, 4 of 8 phases run away along the A-term valley at eta = 0.4 and none at
eta = 0.3, so eta_safe = 0.3 (d = 5: 0.3–0.4). The valley turnover is at eta = 1/3 for d = 6 and 0.4 for d = 5 [D].

**Controls [C].** eta = 0 gives q = 1e-15. Halving the tolerance changes q by 2e-8 (d = 6) and 5e-7 (d = 5).
Extending to a = 1e4 changes q by 2e-8 (d = 6) and 4e-6 (d = 5). The mirror phase flips the sign exactly. Frozen
release scales as h^(3/2) in radiation and h² in matter (0.8 % and 0.6 % between h = 10 and 30). The
physical-time control (`code/control_physical_time.py`) integrates Xddot + 3H Xdot + V' directly and computes
n/s from T(H). It agrees with the conformal yield formula to 3.0e-8 (c_H = 1, eta = 0.1, phi0 = π/2).

**Window per input point.** c_req is the smallest c_d giving Y_Q at the best scanned phase. c_cap is the
boundedness limit, eta ≤ 0.3 at the largest amplitude of the history. c_heavy is the angular mass needed for
isocurvature safety, d² eta c_H/2 ≥ m²/H². Its relaxed form uses the e-fold decay of a light-field S ≈ d δθ to below
S_max (m²/H² ≈ 0.46 for d = 6, c_H = 1). The conservative form is m ≥ 3H/2. c_surv is Codex's late bound. The window
is [max(c_req, c_heavy), min(c_cap, c_surv)].

## 3. Results

### 3a. Required asymmetry [Q, C]

Y_Q = 4.83e-9 per visible entropy, which is n/n_γ = 3.40e-8 (≈ 55 × Y_B) [C]. It gives n_Q = 1.82e-13 GeV³ at 20 MeV [Q]. The
abundance calibration assumes all net charge ends in reference objects plus their debris. A surviving
charge-symmetric remainder would lower the required net yield to Y_Q/(1+f_sym), with v recalibrated (not done).

### 3b. Windows for c_6 (Z_6) [C]

| Inputs (end-of-inflation release) | Generation only | Isocurvature-safe (relaxed) | Upper edge |
|---|---|---|---|
| c_H=1, H_inf=4.7e13, instant reheating | [0.35, 1.61e3] | [137, 1.61e3] (1.1 decades) | boundedness |
| c_H=1, H_inf=4.7e13, T_R=1e14 | [20, 1.61e3] | [137, 1.61e3] | boundedness |
| c_H=10, H_inf=1e13, instant | [0.37, 3.56e3] | [26, 3.56e3] (2.1) | boundedness |
| c_H=10, H_inf=4.7e13, T_R=1e12 | [15, 161] | [15, 161] (1.0) | boundedness |
| c_H=100, H_inf=4.7e13, instant | [2.6e-6, 16.1] | [0.010, 16.1] (3.2) | boundedness |
| c_H=100, H_inf=1e12, instant | [1.85, 3.56e4] | [22, 3.56e4] (3.2) | boundedness |
| any c_H ≤ 100, T_R = 1e9 GeV | empty | empty | — |
| Hubble mass persisting until reheating (W3), all T_R tested | empty | empty | — |

Minimum H_inf for a window (d = 6; the relaxed isocurvature-safe value is the same):

| | instant | T_R=1e14 | T_R=1e12 |
|---|---:|---:|---:|
| c_H = 1 | 2.8e11 | 1.2e12 | none |
| c_H = 10 | 4.6e10 | 7.7e10 | 7.8e12 |
| c_H = 100 | 4.5e9 | 1.4e10 | 2.7e11 |

**The late-survival bound never binds for d = 6.** Inside every window found, c_cap ≤ 1.8e9, far below 7.65e14. c_cap
scales as 3.6e30/(c_H H_inf²) with H in GeV.

### 3c. X^5 (Z_5) [C]

A window exists only for c_H ≥ 10 at H_inf ≳ 1.6e13 GeV (c_H=10, instant) or ≳ 3.8e12 GeV (c_H=100). Isocurvature
safety requires c_H = 100 and H_inf ≳ 2.4e13 GeV. At 4.7e13 GeV, c_5 ∈ [6.0e-5, 1.22e-4] (0.31 decades). There the
**late-survival bound is the binding upper edge**.

### 3d. Further constraints and consequences [A, E]

- **Random phase is excluded [E].** A light angle during inflation gives S ≈ 0.030 cot(phi0)/sqrt(c_H) for d=6. That
  is 3,300 × S_max at c_H = 1. A frozen light field (W4) gives radial S ≈ 0.028. Only the heavy-angle DRT branch
  survives, and its CP phase is a supplied mismatch. This is a retrospective comparison.
- **Charge-symmetric remainder [E, not computed].** At c_req, epsilon ~ 1e-5–0.1. The Codex hidden sector has no
  lighter receiver (visible branch ~1e-27). If the symmetric part survives as cold X X* or anti-objects, the extra
  matter is up to ≈ 0.9 (epsilon^(−4/3) − 1) of the calibrated amount. The window then collapses to eta ≈ 0.3 with a
  tuned phase: H_inf ≈ H_min and c_6 ≈ 4.4e7 (c_H=1), 1.7e8 (10) or 1.8e9 (100). Asymmetric cannibal number-changing
  could instead erase it; that is uncomputed.
- **Washout by the same operator** in a hidden plasma: Γ/H ≲ 1e-8 inside the windows [E].
- **Fragmentation is delayed [C].** No growing mode exists while R > X_c = 1.132 GeV (new point) or 0.707 GeV (old
  point), checked at R/X_c = 1.001–10. This differs from flat-direction AD (Kusenko–Shaposhnikov; Enqvist–McDonald;
  Kasuya–Kawasaki). Onset is at T ≈ 86, 69 and 40 GeV for epsilon = 1, 0.5 and 0.1 (old point: 54, 43, 25 GeV).
  Near X_c (R ≳ 0.97 X_c; old point ≳ 0.9 X_c) the homogeneous k → 0 mode is unstable. The transition is therefore
  first-order-like, and linear theory cannot fix the first fragment scale.
- **Fragment scale [E].** In the developed band (R = 0.95–0.3 X_c, two-field), k_peak = 0.023–0.045 GeV and
  gamma_max = 3e-3–2.7e-2 GeV ≫ H ≈ 1e-14 GeV. The wavelength cube holds **2.6e5–8.6e5 charges** (old point
  8.3e4–1.1e5 at R ≤ 0.7 X_c). A slow-quench estimate at onset gives 5e14–6e15. Codex's supplied 20 MeV dilute condensate
  has a diagnostic of 1.45e11. The reference Q = 2e5 lies near the low end of the bracket but is not selected by it.
- **Z_2 domain walls [E].** AD pins phi = 0 everywhere, so the sign of phi is chosen at random when R reaches X_c. The
  wall tension is σ^(1/3) = 505 MeV (new point) or 455 MeV (old point, = Codex). The existing wall audit finds a
  surviving network incompatible. Annihilation before 1 MeV needs a bias ΔV ≳ 5.8e-26 GeV⁴; a Planck-suppressed odd
  term v⁵/M ≈ 2.8e-19 GeV⁴ would suffice but is absent from the action.

## 4. Verdicts on checked claims

| Claim (8986d34 or hypothesis) | Codex | Mine | Verdict |
|---|---|---|---|
| Rate normalization 4π m_X abs(ε_d)² k abs(A)²/(g²ħ), ε_d formula | — | independent derivation | CONFIRMED [D] |
| X^6 unit-coefficient core loss 3.92e-50 /s; bound c_6 ≤ 7.65e14 (8.8e10 GeV) | full response | free channel 5.590e-50 (= Codex free control), bound 6.41e14 | CONFIRMED (free channel; ×0.70 response factor quoted) |
| X^5 1.54e-12 /s; c_5 ≲ 1.23e-4 (GLOBAL_CHARGE.md, summary) | JSON 1.2227e-4 | free 1.005e-4 | CORRECTED: 1.22e-4 (rounding) |
| Z_6 passes at any c_6 below the survival bound | survival only | boundedness needs c_6 ≲ 0.3 λ_X (M/X_max)² or a stabilizer | CORRECTED (scope) |
| One c_6 can generate Y_Q and survive | hypothesis | window exists for d=6 under the listed inputs; not for T_R ≤ 1e9 or tracking until reheating | CONFIRMED, conditional |
| Same for c_5 | hypothesis | only c_H ≥ 10 at H_inf ≳ 1e13; isocurvature-safe only at c_H = 100 | REFUTED except a marginal edge |
| "AD condensates fragment into Q-balls naturally" | hypothesis | stable until R ≈ X_c; first-order-like onset; walls | CORRECTED |
| Hidden-sector washout, Z_6 matching, other channels | — | — | NOT CHECKED |

## 5. Limits

Classical homogeneous condensate; no lattice or two-field nonlinear evolution. No stabilizing operator. No thermal
masses or hidden-plasma evolution. The isocurvature and H_inf limits are retrospective. The Hubble-induced mass is
assumed to switch off at release, and if it persists, which is common, W3 applies and there is no window. Q_cell is a length-scale diagnostic,
not dN/dQ. The full-response survival factor (0.70) is quoted, not recomputed.

## 6. Credit and read levels

Affleck & Dine, NPB 249 (1985) 361 [R]. Dine, Randall & Thomas, PRL 75 (1995) 398 and NPB 458 (1996) 291 [R].
Kusenko & Shaposhnikov, PLB 418 (1998) 46 [R]. Enqvist & McDonald, PLB 425 (1998) 309 [R]. Kasuya & Kawasaki, PRD 61
(2000) 041301 and PRD 62 (2000) 023512 [R]. Kasuya, Kawasaki & Takahashi, JCAP 0810:017 (2008), AD isocurvature [R].
Frieman, Gelmini, Gleiser & Kolb, PRL 60 (1988) 2101 [R]. Krylov, Levin & Rubakov, PRD 87 (2013) 083528 [R]. Lloyd-Stubbs &
McDonald, arXiv:2008.04339, non-SUSY AD with a quartic [abstract via search]. Planck 2018 X, arXiv:1807.06211, beta_iso
(search-level; table not opened) [Q]. BICEP/Keck, PRL 127 (2021) 151301 [R]. Kibble and Zel'dovich–Kobzarev–Okun on walls [R].

## 7. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/affleck_dine
python3 code/ad_condensate.py          # ~12 min on 4 cores -> results/ad_scan.json
python3 code/window.py                 # -> results/window.json
python3 code/fragmentation.py          # ~2.5 min -> results/fragmentation.json
python3 code/control_physical_time.py  # -> results/control_physical_time.json
git archive 8986d34 research/reality_framework/grut_portal_debris_window_20261023/data | tar -x -C <scratch>
python3 -I code/survival_check.py <scratch>/.../profile_Q200000_lambda0.03_refined.csv 0.8611951308217225 results/survival_free_check.json
```
