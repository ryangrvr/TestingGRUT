# Track gp_matching: early matching at the 7737487 gauge point (handoff item 1 of b5174e5)

25 October 2026 · Claude (independent review track) · Codex packets read only through `git archive` of b5174e5
(`grut_gauge_point_test_20261025/`) and 7737487 (`grut_z6_gauge_bridge_20261025/`) into the scratchpad. Own code in
`code/`, outputs in `results/`. No Codex code was run. No commits.

**This is a conditional, classical, tree-level-plus-estimate calculation inside one supplied model. It is not a
GRUT prediction and nothing here is evidence that X, Y, Sigma or the inflationary inputs exist.** Comparisons with
CMB limits are retrospective. Same-session AI agents are not outside replication.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## 1. Inputs

Stated point (7737487): M_Y=1e12, h=1e11, mu=5e11 GeV, kappa(=y)=1e-3, e=1e-7, lambda_Y=lambda_Sigma=1, lambda_X=1e-3.
Light target g_X=0.1, lambda_phi=0.03, v=0.924089637 GeV, Q=200000 (enter only through phi=0 at high field and the
yield target Y_Q=4.83e-9 [Q, affleck_dine]). Review inputs: H_inf=4.7e13 GeV, -c_H H^2|X|^2 with c_H=1 on X only,
release into radiation with instant reheating (T_R=5.78e15 GeV, g*=106.75), N=50 (and 60) CMB e-folds,
S_max=9.1e-6 [Q, Planck 2018]. Potential as in 7737487 CONSTRUCTION.md; phi is pinned at 0 (mass^2=2g_X^2|X|^2=11.4H^2).

## 2. Method [C]

- `code/model.py`: V in the six canonical real components of X, Y, Sigma (sympy-exact gradient/Hessian). The
  Hessian is directly the mass matrix; no polar ansatz.
- `code/point_full.py`: 48 random *complex* starts, BFGS plus Newton polish; full 6x6 spectrum; gauge-invariant
  reduction G=(W K^-1 W^T)^-1 with W rows psi1=3thX-thY, psi2=2thY-thS, K=diag(2R_i^2); light-mode dictionary;
  stochastic equilibrium; massive-field CMB amplitude (H/2pi) 2^(nu-3/2) Gamma(nu)/Gamma(3/2) e^-(3/2-nu)N.
- `code/boundedness.py`: H=0 global search (400 starts, amplitudes 1e-3..1e4 M_Y, V/(1+|z|^4)), sphere sampling,
  Young and Sigma>M_Y bounds, exact Y,Sigma-eliminated V_eff(|X|) versus the truncated sextic EFT.
- `code/release.py`, `code/instability.py`, `code/variants.py`: homogeneous three-field release in conformal
  variables (A0=0 gauge, zero Gauss charge, e drops out), DOP853; linear-regime runs at delta0=1e-14; tolerance
  and 1e-9-rad shift pairs; 40-point delta0 ensemble; MD and H=1e13 variants.
- `code/running.py`: one-loop scalar betas from the general tensor formula (8 real fields, closed 16-operator
  basis incl. the portals 7737487 sets to zero); control beta(lambda)=20 lambda^2/16pi^2 exact; projection residual
  3e-16. Flat-space Coleman-Weinberg at the point, RG-improved couplings.

## 3. Results

### 3a. Inflationary stationary point and validity [C]

| Quantity | Codex b5174e5 | Mine | |
|---|---|---|---|
| abs(X), abs(Y), abs(Sigma) (GeV) | 1.1197e15, 8.8895e13, 1.2548e13 | 1.119736e15, 8.889516e13, 1.254787e13 (48/48 starts) | CONFIRMED |
| V_min/H^4 | -283.80 | -283.7998 (floor -500) | CONFIRMED |
| Higgs / vacuum | 177.5 | 177.45 | CONFIRMED |
| M_Y/H | 0.0213 | 0.02128 | CONFIRMED |
| kappa abs(X)^2/M_Y^2 | ~1100 (single-field start) | 1104.5 (start); 1254 (full min) | CONFIRMED |
| radial mass^2/H^2 | 0.428, 1.986, 21.61 | same to 1e-12 (lightest is Sigma-radial) | CONFIRMED |
| physical angular masses | 0.378H, 2.751H | 0.37809H, 2.75063H; FD 2e-11; re/im blocks 3e-15 | CONFIRMED |

Analytic check [D]: for abs(X) >> M_Y/sqrt(kappa), Y is held by its quartic, Y = (kappa/2lambda_Y)^(1/3) X, giving
lambda_eff = lambda_X - (3/2) kappa^(4/3) (2lambda_Y)^(-1/3) = 8.8094e-4. Then abs(X)=sqrt(c_H/2lambda_eff)H=23.824H,
abs(Y)=1.8909H, abs(Sigma)=(mu Y^2/2lambda_Sigma)^(1/3)=0.2669H (matches to 1e-4). Quadratic elimination would put
Y 1.58e4 times too high; lambda_Y abs(Y)^2/M_Y^2 = 7.9e3. The matched EFT torque 2|alpha6||X|^6 exceeds the actual
psi2 torque B = 2 mu |Sigma||Y|^2 = 0.0203 H^4 by 1.4e6.

### 3b. Boundedness [C, D]

- Young: kappa=1e-3 < 2.467e-3; numerical margin >= 0 on a 2001^2 grid. Sigma>M_Y branch bound minimum 0.5 >= 0.
- H=0 global search: min V = -1e-20 (round-off at the vacuum); min V/r^4 on spheres r=1e2..1e5 M_Y: +2.3e-4.
- **Induced sextic.** Truncated EFT lambda_X x^4 - (|beta6|+2|alpha6|) x^6 turns negative at abs(X)=30.5 M_Y =
  3.05e13 GeV, where kappa abs(X)^2/M_Y^2 = 0.93: the truncation fails exactly there. The exact eliminated V_eff >= 0
  for all abs(X), with V_eff/abs(X)^4 -> 8.80945e-4 = lambda_eff (1e-7). At the core (abs(X0)~1.7 GeV) the sextic is
  3e-27 of the quartic. Not a physical instability at this point.
- With -c_H H^2|X|^2: floor -c_H^2H^4/(2 lambda_X) = -500 H^4; found -283.80.

### 3c. Quantum running [E]

| | lambda_Y=lambda_Sigma | kappa | lambda_X | lambda_XY (generated) | Young ratio | min V4 on sphere |
|---|---|---|---|---|---|---|
| mu0=M_Y=1e12 | 1 | 1.0000e-3 | 1.0000e-3 | 0 | 37 | 1.16e-4 |
| 4.7e13 | 1.95 | 1.0003e-3 | 1.0058e-3 | 1.0e-6 | 73 | 1.19e-4 |
| 2.2e14 (m_Y at point) | 3.15 | 1.0004e-3 | 1.0082e-3 | 1.6e-6 | 119 | 1.20e-4 |
| 1.12e15 (abs X) | 9.0 | 1.0005e-3 | 1.0107e-3 | 3.0e-6 | 344 | 1.22e-4 |

kappa has no lambda_Y term at one loop, so the stabilizing quartics grow and the Young margin improves: **the
boundedness conclusion survives running.** But lambda_Y reaches 4 pi at 1.43e15 GeV if lambda=1 at M_Y (6.7e16 if at
H_inf): the inflationary point sits near the edge of perturbative control. The stated point has no matching scale.
CW plus RG-improved variants move abs(Y), abs(Sigma) by up to ~50% (Sigma 0.14-0.31H). The light phase stays light:
0.27-0.38H; heavy 2.73-3.20H.

### 3d. Field dictionary [D, C]

With z_i = R_i e^{-i theta_i}, n_i = 2R_i^2 theta_i-dot, V_ang = -A cos psi1 - B cos psi2, A=2kappa R_Y R_X^3 (51.2 H^4),
B = 2 mu R_Sigma R_Y^2 (0.0203 H^4):

  d(a^3 n_X)/dt = -3A a^3 sin psi1,  d(a^3 n_Sigma)/dt = +B a^3 sin psi2,  n_X + 3 n_Y + 6 n_Sigma = 0.

**The charge that becomes X after Y -> 3X decays is N_eff = N_X + 3N_Y = -6 N_Sigma. It is sourced only by
6B sin psi2, exactly at tree level at any amplitude.** In the late EFT, psi1 is slaved and psi2 -> delta = 6thX - thS
= 2psi1 + psi2, recovering the bridge's U6 sin delta with B -> 2|alpha6| R_X^6.

At inflation psi1 is heavy and psi2 is light. The light eigenvector is almost a pure Sigma-phase rotation
(theta_X, theta_Y, theta_Sigma = -0.0019, -0.0037, 2.643 per canonical H), with ddelta = 2.654 and dpsi1 = 0.002.
f_eff = 0.377H ~ sqrt2 abs(Sigma). **So the light inflationary mode is exactly the phase that sources late X charge.**
The gauge direction (mostly theta_X) is eaten (m_A=3.5e-6 H). It carries no gauge-invariant phase.

Consequences:
1. **Alignment is not maintained [C].** Stochastic equilibrium P(psi2) ∝ exp(0.535 cos psi2), so <cos psi2> = 0.26,
   P(|psi2| > pi/2) = 0.34, and relaxation takes ~21 e-folds.
2. **Linear isocurvature [E].** S ≈ 0.036 cot delta0 (N=50; 0.022 for N=60), versus the one-field 0.030 cot phi0.
   That is 4.0e3 (2.5e3) x S_max.
3. **Release instability [C].** The aligned real trajectory gives exactly 0 (Codex CONFIRMED), but it is unstable.
   Angular perturbations grow at 0.263 per conformal time (two tolerances agree to 1e-3), x1e13 between a=35 and 150.
   Onset coincides with kappa abs(X)^2 = M_Y^2 at a≈36 [A: Y-mediated parametric resonance].
   The instability is robust:
   - MD: x2e6 by a=1000;
   - H=1e13 RD: x6e15.

   Any delta0 >= 1e-14 saturates at median |eps_X| = 0.11 (max 0.46). The sign and size are chaotic:
   - 19 sign changes in 40 points across 0.6-1.0 rad;
   - mean eps 5e-4, rms 0.19;
   - a 1e-9-rad shift or a tolerance change decorrelates the runs by a≈98-124.

   Homogeneous N_X does plateau: -104.8±1.2 over a=1000-3000 in one run. A uniform sign would give Y~1e-6,
   about 200 Y_Q, but the sign is not controlled.

**Transfer verdicts.** The random-phase exclusion **transfers and is strengthened**. The light mode is the source
phase, and its map to the charge is chaotic on scales below the 0.036-rad CMB fluctuation. The heavy-angle DRT
escape is absent at this point, and the release instability also threatens it [A, lattice needed]. A tagged
+H^2|Sigma|^2 (not in 7737487) lifts psi2 to 1.0H. It makes the CMB-scale phase fluctuation negligible, but then
the phase aligns. The **remainder numbers fail to transfer.** eps≈4.5 eta does not apply. The charge fraction is
set by the instability (median 0.11, which alone would give 0.9(eps^-4/3 - 1)≈17 if it survived cold). A Y
condensate has E_Y/E_X = 0.43, 0.84 and 2.1 over the windows a=300-600, 600-1000 and 1000-3000, growing ∝ a until
Y->3X [E: Gamma≈3kappa^2 M_Y/256pi^3≈3.8e2 GeV, a≈3.5e5]. That is a symmetric source outside the one-field law.

## 4. Verdicts

| Claim | Codex | Mine | Verdict |
|---|---|---|---|
| M_Y/H=0.0213; kappa X^2/M_Y^2~1100 | as stated | 0.02128; 1104.5 (1254 at min) | CONFIRMED |
| Full tree radial minimum and Higgs 177.5x | 1.1197e15/8.8895e13/1.2548e13 | same to 6 digits | CONFIRMED |
| Phase masses 0.378H, 2.751H; one light | finite-difference | full 6x6 canonical Hessian | CONFIRMED |
| Full potential bounded; truncated runaway is no theorem | Young + quartics | numerical + exact V_eff, lambda_eff=8.81e-4 | CONFIRMED |
| Quantum running "not computed" | — | one-loop: stability survives; lambda_Y non-perturbative at 1.4e15 (mu0=M_Y) | NOT CHECKED by Codex; estimate added |
| One-field numbers need a new dictionary | stated | built; light mode = source phase psi2 | CONFIRMED |
| "Masses do not supply an isocurvature prediction" | — | they do: S≈0.036 cot delta0, then chaos | CORRECTED |
| Aligned phase gives zero charge | exact | exact but unstable (x1e13), and not maintained in inflation | CORRECTED (scope) |
| Our REVIEW_XCHARGE: "isocurvature, wall and remainder constraints apply unchanged" | — | isocurvature: yes (stronger); remainder: no | CORRECTED (self) |

## 5. Limits and open items

The calculation is classical, homogeneous and in flat-space CW. Specific limits:
- Release choices are supplied: instant switch-off, RD or MD, zero velocities.
- No k≠0 modes or lattice: chaotic homogeneous sign plus inhomogeneous phases suggests ±charge domains and a
  suppressed net yield [A].
- No Y decay, thermal or hidden-plasma physics.
- Other light spectators are not examined, in particular D with m≪H.
- **New flag [E/R]:** the vector is light during inflation (m_A/H=3.5e-6). By Graham–Mardon–Rajendran scaling,
  longitudinal production at the vacuum mass 6e4 GeV gives rho/s≈0.3 GeV. Decay to X X̄ would then give ~1e3
  pairs per net charge. This is not computed and is for the leftovers track.

## 6. Credit

Affleck & Dine 1985; Dine, Randall & Thomas 1995/96; Starobinsky & Yokoyama 1994 (stochastic equilibrium);
Coleman & Weinberg 1973; Machacek & Vaughn 1983-85 (general one-loop betas); Greene, Kofman, Linde & Starobinsky
1997 (parametric resonance); Kasuya, Kawasaki & Takahashi 2008 (AD isocurvature); Graham, Mardon & Rajendran 2016
(vector from inflation); Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki (Q-ball formation); Planck 2018 X
[Q]. All [R] unless marked.

## 7. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/gp_matching
python3 code/point_full.py     # 11 s  -> results/point_full.json
python3 code/boundedness.py    # 60 s  -> results/boundedness.json
python3 code/release.py        # ~4 min -> results/release.json
python3 code/instability.py    # ~6 min, 4 cores -> results/instability.json
python3 code/variants.py       # 1 min -> results/variants.json
python3 code/running.py        # 40 s  -> results/running.json
```
