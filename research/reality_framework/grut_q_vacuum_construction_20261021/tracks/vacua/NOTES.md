# Track `vacua`: QV1 with more than one matter vacuum, and the Standard Model relation it can imply

21 October 2026. Stage 2 of the q-theory construction (issue #13). Briefs: `../../BRIEF.md`, `../../BRIEF_STAGE2.md`.
Stage-1 records used: `../equations/NOTES.md` (cited "equations §n"), `../literature/NOTES.md` ("literature §n"),
`../ledger/NOTES.md` ("ledger §n").

**Labels.** [D] derived here; [C] computed by code in `code/` (output in `results/`); [Q] stated by the cited source
(paraphrased, with location and read level); [E] estimate; [A] argued; [R] recalled, not re-read here.
Conventions as in the brief: signature (-,+,+,+), K = M^2 = 1/(8 pi G), M = 2.435e18 GeV, rho_Lambda = 2.52e-47 GeV^4.
`eps*(mu) = sup_q [mu q - eps(q)]` is the Legendre (convex) conjugate of eps.

**Scope.** No new theory, no empirical validation, no confirmatory test is claimed. Every comparison with an existing
measurement is **retrospective**. A comparison is called held out only when the data do not exist yet.

**Reproduction.** `cd code && python3 run_all.py` (about 1 minute; numpy 2.5, scipy 1.18, sympy 1.14, mpmath 1.3).
It writes `results/vacua_q_sector.json` (§1, §3.4), `results/nucleation.json` (§2) and `results/sm_criticality.json`
(§5). `results/sources.json` lists every source with its read level. Downloads (arXiv TeX, PDG files) are in
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/vacua_dl/`. None is in the repository.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdicts in one screen

| # | Result | Label |
|---|---|---|
| 1 | With one global mu and heavy q, matter vacuum i gravitates with Lambda_i = F(V_i - eps*(mu)), F(x) = x + (4/3) x^2/(K m_q^2) + a3 x^3 + ..., F(0) = 0, F' > 0 for x < ~(3/8) K m_q^2. Vacuum-energy **differences gravitate fully**: Lambda_i - Lambda_j = (V_i - V_j)[1 + O(x/(K m_q^2))]; the correction is < 0.5% for vacua down to |V| ~ 1e71 GeV^4, the SM potential at the reduced Planck scale (natural m_q). With K' = 0 it is exact: Lambda_i = V_i - eps*(mu). | [D][C] |
| 2 | mu is continuous across a matter domain wall (the wall carries no 3-form charge). q follows the local quasi-static value when m_q w >> 1 (the physical case: m_q w ~ 1e16 for an electroweak wall); for m_q w << 1, K(q) is continuous with a kink [n.dK] = S/3 = -sigma, the trace of the extrinsic curvature does not jump and neither does the traceless part for a pure-tension wall (the f(R) junction conditions in q-language). | [D][C] |
| 3 | A q-dependent matter vacuum energy needs a coupling of matter to the 3-form, e.g. -int theta(phi) F~ (or a q|H|^2 term, CP-odd since q is a pseudoscalar). Then Lambda_i = V_i - eps*(mu - theta_i): differences no longer equal V_i - V_j. Hierarchy-safe |H|^2 couplings change them by < 1e40 GeV^4, negligible against Planck-deep differences (~1e70). | [D][E] |
| 4 | Nucleation. For any rate Gamma > 0, every worldline eventually lies inside a true-vacuum bubble: P(no bubble) = exp(-Gamma V4), with V4 ~ t^4 (Minkowski, matter era) or ~ (4 pi/3) t/H^3 (de Sitter). In de Sitter the false-vacuum *volume* still grows iff Gamma/H^4 < 9/(4 pi). For the SM rate (~1e-577 Gyr^-1 Gpc^-3) Gamma/H_L^4 ~ 1e-574. An AdS interior crunches within pi sqrt(3K/|Lambda_in|) (3e-41 s for |Lambda_in| = 7.5e70 GeV^4). | [D][C][Q] |
| 5 | Readings of P_q (QV1 says "the final T -> 0 vacuum"): (a) global minimum: Lambda_ours = V_ours - V_true >= 0; (b) our vacuum: self-consistent under QV1 only if Gamma = 0; (c) the vacuum the whole ends in: with the worldline reading the only self-consistent choice is the true vacuum (reduces to (a)); with the volume reading in de Sitter, (b) survives; (d) coexistence at shared mu and zero pressure: **not implied by the action** (mu is shared automatically, equal zero pressures need V_i = V_j, which nothing in the q-sector supplies); it is Froggatt-Nielsen's principle added as an assumption. | [D][A] |
| 6 | Under (a), (c)-worldline, and (b) with QV1 uncalibrated: our vacuum must be the global minimum to within rho_Lambda. With the SM valid to the Planck scale this is absolute stability: **M_t <= 170.94 + 0.495 (M_h - 125.20) + 0.353 (alpha_s - 0.1180)/0.0009 GeV (+-0.11 th)**, the equality being the criticality (multiple-point) line lambda(mu*) = beta_lambda(mu*) = 0, mu* = 4.7e17 GeV. Inverse: M_h,crit = 128.54 + 2.02 (M_t - 172.56) - 0.71 (alpha_s - 0.1180)/0.0009 GeV. | [C] |
| 7 | Validation: lambda(M_Pl) = -0.01437 (Buttazzo et al. -0.0143), g1, g2, g3, y_t at M_Pl to 4 digits, log10 Lambda_V = 9.497 (9.5), Lambda_I/Lambda_V = 12.7 (~13), M_t,crit = 171.38 at Bednyakov et al.'s inputs (their 171.44), 171.00 at PDG 2024-like inputs (Hiller et al. 171.10). | [C][Q] |
| 8 | Retrospective (PDG 2025): M_t(direct) = 172.56 +- 0.31 lies 1.62 GeV above the bound, **z = 3.3** (2.8 with a 0.3 GeV MC-to-pole theory error); with the cross-section pole mass 172.4 +- 0.7, z = 1.8. Under readings (a)/(c), QV1 with the SM to the Planck scale is in retrospective tension at that level. Without the relation, (a) at PDG central values gives Lambda_ours >= 7.5e70 GeV^4 = 3e117 rho_Lambda. | [C] |
| 9 | Held out: FCC-ee threshold m_t (~40 MeV), alpha_s to 0.1%, m_h to 4 MeV give sigma_tot ~ 0.13 GeV (matching-limited): a ~12 sigma decision if central values stay; HL-LHC m_t ~0.2 GeV [R] plus lattice alpha_s halved gives ~5.6 sigma. | [C][Q][R] |
| 10 | Novelty and credit. The SM relation is Froggatt-Nielsen's (criticality) and Degrassi/Buttazzo/Bednyakov's (stability bound), recomputed. "Vacuum energy cancels at the global minimum, so a false vacuum gravitates with its splitting" is AHHKM's premise (astro-ph/0005111) and KV's "true vacuum of global equilibrium" (1102.3152 §6). What this track adds is bookkeeping: the action-level Lambda_i, the junction behaviour of q, and the fixed-point argument that QV1's final-state condition plus nucleation (worldline reading) forces the true-vacuum reading and hence the stability inequality. Nothing is GRUT-native. | [A] |

---------------------------------------------------------------------------------------------------------------------

## 1. Gravitating vacuum energy of several matter vacua (action of BRIEF.md, one global mu)

### 1.1 Set-up [D]

S = int sqrt(-g)[K(q) R/2 - eps(q)] + S_m[phi, g], with S_m containing -int sqrt(-g) V(phi). The matter potential
V(phi) is the full zero-temperature effective potential; its minima phi_i have energies V_i (all loop and condensate
constants included, as in QV1's eps_inf). The matter action does not depend on q (KV's choice; literature §2.6
item 2). The q-equation is (q) eps'(q) - K'(q) R/2 = mu, with one constant mu on the connected spacetime
(equations §1.2). In a static maximally symmetric state of vacuum i (R = 4 Lambda_i/K, H^2 = Lambda_i/(3K)) the
two conditions are

    E1:  eps'(q_i) - 2 K'(q_i) Lambda_i / K(q_i) = mu
    E2:  Lambda_i = V_i + eps(q_i) - mu q_i

### 1.2 K' = 0: exact Legendre form [D]

E1 gives eps'(q_i) = mu, so q_i = q_mu is **the same in every matter vacuum**, and

    Lambda_i = V_i + eps(q_mu) - mu q_mu = V_i - eps*(mu),        Lambda_i - Lambda_j = V_i - V_j   (exactly).

This is the 4-form / Henneaux-Teitelboim result (equations §1.5): the q-sector adds one common constant -eps*(mu) to
every matter vacuum and nothing else. [R] Brown-Teitelboim and Henneaux-Teitelboim for the K' = 0 structure.

### 1.3 K' != 0 [D][C]

Write x_i = V_i - eps*(mu), the vacuum energy measured from the level at which this mu is in Minkowski
equilibrium (c*(mu) of equations §5.2). Solving E1, E2 as series in x with generic Taylor data at q0
(`code/vacua_q_sector.py`, part A; sympy):

    Lambda_i = x_i + a2 x_i^2 + a3 x_i^3 + O(x^4),
    a2 = 2 K'^2/(K^2 eps'') = 4/(3 K m_q^2),   a3 = 8 K'^2 (3 K'' eps'' - K' eps''')/(3 K^3 eps''^3),
    q_i - q0 = 2 K' x_i/(K eps'') + O(x^2),     dK/K (vacuum i vs Minkowski) = (4/3) x_i/(K m_q^2),

with m_q^2 = (2/3) K eps''/K'^2. This reproduces equations §5.2 (same a2, a3) for a constant shift; here it is
applied vacuum by vacuum at a **shared** mu. Exact slope along the static family (implicit differentiation, [D][C]):

    dLambda/dV = K^2 D / (K^2 D - 4 Lambda K'^2),   D = eps'' - 2 Lambda (K K'' - K'^2)/K^2,

equal to 1 at Lambda = 0 and positive for Lambda < K^2 D/(4K'^2) ~ (3/8) K m_q^2. So Lambda_i is a monotonic
function of V_i, and **Lambda_i = Lambda_j if and only if V_i = V_j** (used in §3.4). Differences:

    Lambda_i - Lambda_j = (V_i - V_j) [1 + a2 (x_i + x_j) + O(x^2)].

Numerical check [C] (KV 0806.2805-type eps = (-q^2 + q^4/3)/2, K = 1 + K'(q - 1), units q0 = K0 = eps''(q0) = 1):
for K' = 1 and x = 0.01 the exact static root is Lambda = 0.010190 against 0.0102 from two terms; for K' = 0.3,
0.0100177 against 0.010018; for K' = 0, Lambda = x exactly (`results/vacua_q_sector.json`, A_numeric_KV_eps).

Size at natural scales [E]: m_q ~ 2e18 GeV (equations §7). For a Planck-deep SM region, x ~ 7.5e70 GeV^4
(§5.4): a2 x = (4/3) x/(K m_q^2) ~ 4e-3. For electroweak-scale differences (~1e8 GeV^4): ~6e-66. Differences
gravitate fully, as in GR.

### 1.4 Matter vacuum energy that depends on q [D][E]

The matter sector can feel q only through the 3-form, i.e. S_m must contain F~ = q sqrt(-g):
- a metric-independent axion-type term -int d^4x theta(phi) F~ (theta a pseudoscalar function of the matter fields;
  [R] Dvali hep-th/0507215, Kaloper-Sorbo 0811.1989 for such couplings), or
- a term -int sqrt(-g) xi q |H|^2 (dimension 4, since [q] = 2). Because q is a pseudoscalar density, this one is
  CP-odd ([Q] Kaloper-Westphal 1907.05837, abstract: Higgs-4-form couplings give direct CP violation;
  Giudice-Kehagias-Riotto 1907.05370, abstract: 4-form fluxes correlated with the Higgs mass).

With L_m containing -V(phi, q), mu = eps' + d_q V - K'R/2 is still one global constant (the Gauss law is unchanged).
For the linear case V(phi_i, q) = V_i + theta_i q and K' = 0 [D; sympy check, part B]:

    Lambda_i = V_i - eps*(mu - theta_i),   q_i = (eps*)'(mu - theta_i),
    Lambda_i - Lambda_j = (V_i - V_j) + q_mu (theta_i - theta_j) + O(theta^2).

So the gravitating differences are no longer the matter differences, and q now differs between the phases (this is
what §3.4 needs for any Maxwell construction). Size [E]: a Higgs coupling xi q|H|^2 shifts the Higgs mass parameter
by ~xi q0 ~ xi M^2; keeping that below m_h^2 needs |xi| < ~3e-33, and then the shift of a Planck-deep difference
is |xi| q0 Delta|H|^2 < ~ m_h^2 h*^2 ~ 1e40 GeV^4, i.e. ~1e-30 of V differences ~1e70 GeV^4. Within QV1 (no such
coupling) the case does not arise; it is a different version.

### 1.5 Junction across a matter domain wall [D][C]

Wall of tension sigma and thickness w between vacua i and j; unit normal n from j to i; [X] = X_i - X_j.

1. **mu is continuous.** The A-equation is d_nu(dL/dq) = 0 as a distribution. A wall of the matter field carries no
   3-form charge (no int_W A term), so [mu] = 0 [D]. Contrast: a Brown-Teitelboim membrane with charge e gives
   [mu] = e ([R] Brown-Teitelboim 1988); that is how mu could change, and it is outside the action.
2. **Heavy q (m_q w >> 1).** q tracks its local quasi-static value, K(q(x)) = K_qs(T(x), V(phi(x))); away from the wall
   q -> q_i, q_j with q_i - q_j = 2K'(Lambda_i - Lambda_j)/(K eps'') + ... (zero if K' = 0). Gravity obeys the GR
   Israel conditions with tension sigma, up to O(1/(m_q w)^2). Linear check [C] (part D):
   (d^2/dz^2 - m_q^2) dK = T/3 = -rho_w(z). The deviation from local tracking dK = rho_w/m_q^2 is 1.8% at
   m_q w = 10 and 0.12% at m_q w = 30, scaling as ~2/(m_q w)^2.
3. **Light q, thin wall (m_q w << 1).** The δ-function parts of (g) and (q) give [D]:
   Write Kx_ab for the extrinsic curvature, to keep it apart from the coupling K(q).
   - [K] = 0, i.e. q is continuous;
   - [Kx] = 0 (trace of the extrinsic curvature), because R = 2(eps' - mu)/K' is a function of q and has no δ part;
   - [n.dK] = S/3, where S is the trace of the surface stress; for a pure-tension wall S_ab = -sigma h_ab, so
     [n.dK] = -sigma;
   - -K [Kx_ab] = S_ab - h_ab S/3, which vanishes for a pure-tension wall.

   So the wall's tension is carried by a kink of q, spread over 1/m_q, and the extrinsic curvature does not jump.
   These are the f(R) junction conditions written in q ([R]/abstract: Deruelle-Sasaki-Sendouda 0711.1150; Senovilla
   1303.1408, who also finds the trace condition). The linear check gives a kink of -0.994 sigma at m_q w = 0.01 and
   -0.997 sigma at m_q w = 0.1 (corrected for the e^{-m_q z} fall-off) [C].
4. **Physical regime [E].** An electroweak-scale wall (w ~ 1/m_h) has m_q w ~ 1.6e16, so it is in case 2. An SM
   vacuum-decay bounce of size R ~ 1/(1e17 GeV) has m_q R ~ 20, so the q-mode adjusts adiabatically; the corrections
   are of order (m_q R)^-2 ~ 3e-3. A light q (m_q ~ meV, the Eöt-Wash floor, equations §4.3) would put every particle
   physics wall in case 3.
5. **Consequence [D/A].** With heavy q, bubble nucleation between matter vacua is GR vacuum decay with
   Lambda_i = F(V_i - eps*(mu)) inside and outside, with mu unchanged. The q-sector does not alter the CDL
   exponents beyond O(1/(m_q R)^2), and it does not readjust Lambda inside a bubble.

---------------------------------------------------------------------------------------------------------------------

## 2. Nucleation and "the final vacuum" (`code/nucleation.py`, `results/nucleation.json`)

### 2.1 Past-light-cone 4-volume [D][C]

The 4-volume of the past light cone of a point at time t is V4(t) = (4 pi/3) int a^4 (eta - eta')^3 deta'.

| Background | V4 | Check |
|---|---|---|
| Minkowski, started at t = 0 | pi t^4/3 | analytic |
| matter domination | (16 pi/1485) H^-4 = (81 pi/1485) t^4 = 0.1714 t^4 | [C] to 1e-15 |
| LambdaCDM today (Omega_m = 0.315, H0 = 67.4) | 0.1230/H0^4, with eta0 = 3.18/H0 and t0 = 13.79 Gyr | [C] by two independent integrations; Buttazzo et al. quote ~0.15 [Q 1307.3536 sec. 6.3] |
| de Sitter (late LambdaCDM) | dV4/dt -> (4 pi/3)/H_L^3; computed 7.389896 vs exact 7.389896 (H0 = 1 units) | [C] |

### 2.2 Does every point end inside a bubble? [D]

- **Worldline.** P(no bubble in the past light cone) = exp(-Gamma V4(t)) tends to 0 whenever V4 -> infinity. That
  holds in asymptotically Minkowski space (V4 ~ t^4), in a matter-dominated future (t^4) and in de Sitter space (V4 ~
  (4 pi/3) t/H^3). **For every Gamma > 0, every worldline eventually enters a true-vacuum bubble**, either a bubble's
  interior or its future light cone. Only Gamma = 0 avoids this.
- **Volume (de Sitter only).** The comoving false-vacuum fraction is f = exp[-(4 pi/3)(Gamma/H^4) H t], so the
  physical false-vacuum volume a^3 f grows iff **Gamma/H^4 < 9/(4 pi) = 0.716**. Below this the false vacuum is
  eternal by volume while each worldline still decays ([R] Guth-Weinberg 1983 on non-percolation; [Q abstract]
  Arkani-Hamed-Dubovsky-Senatore-Villadoro 0801.2399 for the window between "not yet decayed" and "not future
  eternal").

### 2.3 Standard Model numbers [Q rates; C]

Rates: Chigusa-Moroi-Shoji 1707.09301 (abstract) give ~1e-577 Gyr^-1 Gpc^-3, with a range of 1e-295 to 1e-1465
from M_t. Andreassen-Frost-Schwartz 1707.08124 (abstract) give a lifetime of 1e139 yr, and > 1e58 yr at 95%.
The two papers differ in method and inputs; I only use the first as a range.

| Rate (Gyr^-1 Gpc^-3) | log10 Gamma (GeV^4) | log10 P(decayed in our past cone) | log10 Gamma/H_L^4 | worldline decay time, dS future | if the future were Minkowski / matter era |
|---|---|---|---|---|---|
| 1e-577 | -741.3 | -574.8 | -573.6 | 1e574 Gyr | 1e145 Gyr |
| 1e-295 | -459.3 | -292.8 | -291.6 | 1e292 Gyr | 1e74 Gyr |
| 1e-1465 | -1629.3 | -1462.8 | -1461.6 | 1e1462 Gyr | 1e367 Gyr |

All rates give Gamma/H_L^4 << 9/(4 pi): with Lambda_ours = rho_Lambda the false vacuum is eternal by volume.

### 2.4 AdS interiors [E]

Inside a bubble with Lambda_in < 0 (open slicing, Lambda-dominated) a(t) = l sin(t/l), with l = sqrt(3K/|Lambda_in|),
so the crunch comes at t_c = pi l. For the following values of |Lambda_in|:

| |Lambda_in| (GeV^4) | Case | t_c |
|---|---|---|
| 7.5e70 | the SM potential at 2.4e18 GeV for PDG central values, §5.4 | 3e-41 s |
| 2e65 | at 1e17 GeV | 2e-38 s |
| rho_Lambda | a vacuum degenerate to within rho_Lambda | 1.7e18 s = 55 Gyr |

Within the SM the interior field keeps rolling toward Planckian values, so these numbers are order-of-magnitude
only. For exact AdS, the t = pi l singularity of the open slicing is a coordinate artifact. With the interior scalar
field evolving, it is a genuine big crunch ([R] Coleman-De Luccia 1980; Abbott-Coleman 1985). **A crunching region has no T -> 0 final state, so P_q has no referent there** [A].

### 2.5 Gravitational quenching near criticality [E]

Coleman-De Luccia ([R] 1980): decay of a Minkowski false vacuum into AdS is forbidden if eps <= (3/4) sigma^2/M^2
(thin wall). On the SM potential just above the critical line, eps exceeds this threshold by 6e3 already 0.001 GeV
above M_t,crit, and by 1e17 at 0.05 GeV (`sm_criticality.json`, CDL_quench_estimate). Gravity therefore moves the
stability boundary by < 1e-3 GeV, which is negligible. In this context "absolutely stable" and "gravitationally stable"
coincide to working precision. The estimate uses a thin wall, an RG-improved tree potential and Planckian field
values, so it is an order of magnitude only ([R] Isidori-Rychkov-Strumia-Tetradis 0712.0242 for the full treatment).

---------------------------------------------------------------------------------------------------------------------

## 3. Readings of P_q with several matter vacua

QV1's P_q (BRIEF_STAGE2): mu = mu0[eps_inf], where eps_inf = eps(q) plus the q-independent vacuum energy of "the
final T -> 0 vacuum". With one matter vacuum the final vacuum is unique. With several, P_q needs an interpretation.
In every reading mu is an integration constant: the dynamics accepts any value (equations §6.2; ledger §1.4). The
question is which value the state condition names, and whether that value is self-consistent with what the whole
then does. With K' corrections negligible (§1.3), write Lambda_i = V_i - V_f, where f labels the vacuum that P_q names.

### 3.1 (a) Final = the global minimum of the matter potential [D/A]

- **Added assumption.** QV1 says "final", not "global". This reading is KV's "true vacuum of global equilibrium"
  (1102.3152 §6, via literature §3.2). It is also AHHKM's premise that an unknown mechanism cancels the vacuum
  energy at the global minimum (astro-ph/0005111, via literature §5.5).
- **Consistency with conserved mu: yes.** It is a value of mu. When Gamma > 0 it agrees with (c), the reading in
  which the universe actually ends there (§3.3).
- **Our Lambda:** Lambda_ours = V_ours - V_true >= 0. The sign is a consequence; that the observed sign is positive
  is retrospective and weak. With QV1-cal's offset, Lambda_ours = V_ours - V_true + rho_Lambda.
- **Size.** Since V_true <= V(h) for every h in the region where the SM holds, at PDG 2025 central values (§5.4)
  Lambda_ours >= |V_eff(2.4e18 GeV)| = 7.5e70 GeV^4 = **3e117 rho_Lambda** (2e65 GeV^4 even if the SM is trusted only
  to 1e17 GeV). This is retrospectively excluded, unless the true parameters lie on the stable side of the critical
  line, i.e. our vacuum is the global minimum to within rho_Lambda. The tolerance is delta lambda ~ 4 rho_Lambda/h*^4
  ~ 1e-118, so the condition is exactly the absolute-stability inequality of §5.
- **Beyond the SM.** Reading (a) also requires that no deeper vacuum exists anywhere in the UV completion. That is
  generically false in landscapes with AdS vacua [A]. The SM inequality is the testable necessary condition.

### 3.2 (b) Final = our vacuum [D]

- **Added assumption** (anthropic or "present" reading).
- **QV1 (uncalibrated).** Lambda_ours = 0, so our vacuum is Minkowski; Lambda_true = V_true - V_ours < 0, so bubbles of
  true vacuum are AdS and crunch (§2.4). In Minkowski V4 ~ t^4, so with Gamma > 0 every worldline ends in a crunch
  (§2.2), and the whole does **not** end in our vacuum. The reading is self-consistent only if Gamma = 0: absolute
  stability. Gravitational quenching gives no relief at the precision used here (§2.5). The SM consequence is the same
  inequality as (a). (QV1's Lambda_ours = 0 already conflicts with rho_Lambda,obs; that is QV1's known residual
  problem, not this track's result.)
- **QV1-cal.** Lambda_ours = rho_Lambda > 0. Each worldline still decays (in ~1e574 Gyr), but by volume the false
  vacuum is eternal (Gamma/H_L^4 ~ 1e-574 << 0.716). Read by volume, (b) is self-consistent and implies **no** SM
  relation: metastability is allowed.

### 3.3 (c) Final = the vacuum the whole ends in, given nucleation [D/A]

This is the literal reading of QV1, but it needs a measure.
- **Worldline measure (added assumption: "every observer's asymptotic vacuum").** Treat P_q as a fixed-point
  condition: mu = mu0 of V_f, and f must be where worldlines end.
  - Try f = ours with Gamma > 0. Then Lambda_true < 0, every worldline enters an AdS crunch (§2.2, §2.4), no T -> 0
    final state exists, and P_q has no referent. **Not self-consistent.**
  - Try f = true. Then Lambda_true = 0 (Minkowski interiors, no crunch) and Lambda_ours = V_ours - V_true > 0 (de
    Sitter). Every worldline ends in the true vacuum. **Self-consistent.**
  - If Gamma = 0, then f = ours and our vacuum is the global minimum: (a) and (b) coincide.

  So with Gamma > 0 the worldline reading **reduces to (a)**, with the same SM consequence.
- **Volume measure (de Sitter).** If Lambda_ours > 0 and Gamma/H^4 < 9/(4 pi), the false vacuum dominates the volume
  forever, f = ours is also self-consistent, and there is no SM relation. Both fixed points exist. Choosing between them
  is the eternal-inflation measure problem, which the action does not settle [A].
- **Consistency with conserved mu.** Yes, both fixed points are values of mu. The selection is a final boundary
  condition (ledger §1.3 d).

### 3.4 (d) Phase coexistence at shared mu and zero external pressure (KV/Volovik thermodynamics; MPP) [D]

**Derivation within the action.**
1. *Shared mu* holds automatically: there is one global integration constant, continuous across walls (§1.5).
2. *Equal pressure.* The vacuum pressure is P_i = -Lambda_i. Two macroscopic phases coexist without net conversion
   only if P_i = P_j. By §1.3 (F monotonic) this holds iff **V_i = V_j**. A planar wall with Lambda_i != Lambda_j is
   pushed toward the higher-Lambda side.
3. *Zero external pressure* in both phases, P_i = P_j = 0, gives V_i = V_j = eps*(mu). That is two conditions on one
   integration constant, so **one relation among matter couplings is needed**: degeneracy, the multiple-point condition.
4. *Does the dynamics force it?* No. V_i and V_j are matter parameters, and the q-sector enters every phase through
   the same eps*(mu). No field equation relates V_i - V_j to anything in the q-sector. **Refuted as a consequence of
   the action**; it is Froggatt-Nielsen's principle added to QV1.

**Fixed extensive quantity in q-language** (FN's justification: a fixed extensive quantity, e.g. int d^4x |phi|^2,
forces a mixture, as for slush at fixed energy [Q hep-ph/9511371 §1]).
- The 4-form's own extensive quantity is Q4 = int F~ d^4x = int sqrt(-g) q d^4x. It is conjugate to mu and fixed by
  Dirichlet data on A ([R] Brown-Teitelboim; Duncan-Jensen).
- Fixing Q4 forces a mixture only if a Maxwell point exists: a mu_c with P_i(mu_c) = P_j(mu_c) and q_i != q_j.
  - With q-independent matter and K' = 0: q_i = q_j = (eps*)'(mu), and P_i - P_j = V_j - V_i for every mu (sympy,
    part C). **No Maxwell point**: Q4 cannot tell Higgs phases apart.
  - With K' != 0: P_i - P_j = -[F(V_i - eps*) - F(V_j - eps*)], which is nonzero for every mu when V_i != V_j.
    Again no Maxwell point.
  - With k 4-forms: dP_i/dmu_a = q_a, identical rows for all phases, so the Jacobian has rank 1 and at most one
    non-degenerate phase can sit at P = 0.

  This is Volovik's counting nu <= k (hep-ph/0309144 §3, full read: "coexisting vacua must have the same chemical
  potentials ... P^(i)(mu_1..mu_k) = 0 ... satisfied simultaneously if nu <= k"). In the 4-form realization the
  counting needs charge densities that differ between phases, i.e. matter coupled to q (§1.4).
- With a coupling theta(phi) F~, a mu_c exists (closed form in `vacua_q_sector.json`, part C). Zero pressure at mu_c
  is then one further relation, involving the new coupling theta, so it is not a relation among SM quantities.
- Volovik's paper reaches the compatible conclusion that the adjustment of the deep vacuum "does not lead to any
  sizable correlation between the parameters of the Standard Model" (§2, ratio E_ew^4/E_Pl^4), and that coexistence
  is a stronger, separate assumption regulated by soft global charges (baryon, lepton, family) [Q].

**If added anyway (reading (d) = MPP + P_q):** all coexisting (or all reachable final) vacua have Lambda = 0, so our
vacuum is degenerate with the deeper one: the **equality** M_t = M_t,crit(M_h, alpha_s) of §5, with mu* free.
Froggatt-Nielsen's original form also fixes mu* ~ M_Pl, giving two relations (§5.4). The "agnostic final" variant
(P_q must hold for whichever vacuum the whole ends in, not known in advance) is equivalent [A].

### 3.5 Summary table

| Reading | Added assumption? | Consistent with conserved mu? | Lambda_ours | SM implication (SM to M_Pl) |
|---|---|---|---|---|
| (a) global minimum | yes (KV 1102.3152 §6; AHHKM) | yes | V_ours - V_true >= 0 | absolute stability: M_t <= M_t,crit |
| (b) ours, QV1 | yes | only if Gamma = 0 | 0 | absolute stability (Gamma = 0) |
| (b) ours, QV1-cal, volume measure | yes | yes | rho_Lambda (calibrated) | none |
| (c) worldline end-point | measure choice | yes; unique fixed point = true vacuum if Gamma > 0 | as (a) | as (a) |
| (c) volume end-point | measure choice | yes (second fixed point) | as (b)-cal | none |
| (d) coexistence at P = 0 / MPP | yes (FN), not derivable | mu shared yes; equal P needs V_i = V_j | 0 (or rho_Lambda with cal) | criticality: M_t = M_t,crit |

---------------------------------------------------------------------------------------------------------------------

## 4. Literature (search, read levels, credit)

**Search.**
- arXiv API:
  - abs:"multiple point principle" (2015-2026);
  - au:Volovik / au:Klinkhamer combined with "degenerate vacua", "multiple point", "false vacuum", "coexist",
    "bubble";
  - abs:"q-theory" combined with MPP / Higgs / false vacuum;
  - vacuum-stability titles 2018-2026;
  - FCC-ee top threshold.
- INSPIRE co-citation:
  - papers citing any core KV q-theory paper (0711.3170, 0806.2805, 0905.1919, 0907.4887, gr-qc/0604062) **and**
    FN 1995: **1 record**, KV 2111.07962. There the MPP is invoked for equal rho_V of the q = 0 and q = q0 vacua, in
    the q-sector, not for the Higgs; read the citation context.
  - core q-theory and (Buttazzo 2013 or Degrassi 2012): 3 records. These are Jegerlehner 1402.3738 and Bass 2508.15289
    plus a Bass book (emergent SM, vacuum stability, Volovik's Gibbs-Duhem nullification). Read at abstract level plus
    citation-context grep.
  - citing Volovik hep-ph/0309144: 19 records; citing KV gr-qc/0503090: 21 records (titles scanned).
- **No paper found derives an SM relation from q-theory's final-state condition.**

**Read levels and what is used.** Full list in `results/sources.json`.

| Source | Read level | Used for |
|---|---|---|
| Froggatt-Nielsen hep-ph/9511371, PLB 368 (1996) 96 | §§1-2 (TeX) | MPP: degenerate minima, second at ~M_Pl, lambda = beta = 0; fixed-extensive-quantity justification; M_t = 173 +- 5, M_h = 135 +- 9 |
| Froggatt-Nielsen-Takanishi hep-ph/0104161; FN(L) hep-ph/0312218, hep-ph/0406110; Bennett-Nielsen hep-th/0503080, hep-ph/9311321; FN-Nevzorov-Thomas 1704.08453; Sidharth et al. 1605.01169, 1801.06979 | abstract | borderline metastability (M_h = 121.8 +- 11); MPP as zero CC in several vacua; non-locality as the origin |
| Volovik hep-ph/0309144, Pisma ZhETF 79 (2004) 131 | full main text (TeX) | coexistence = shared chemical potentials + P = 0; nu <= k; no SM correlation from deep-vacuum adjustment |
| Klinkhamer-Volovik gr-qc/0503090, PLA 347 (2005) 8 | §§1-2 (TeX) | coexisting vacua; a dS bubble stabilised by trapped fermionic charge |
| Klinkhamer 1302.1496, JETP Lett. 97 (2013) 297 | full (TeX) | SM second degenerate minimum ~ E_Planck/sqrt(8 pi xi); M_t ~ 171.4 + (M_h - 126)/2 + (alpha_s - 0.1184)/0.0028 +- 0.5 |
| Klinkhamer-Volovik 2111.07962 | context of the FN citation | MPP for q-vacua |
| KV 1102.3152 §6; AHHKM astro-ph/0005111 | via literature §3.2, §5.5 | reading (a) |
| Degrassi et al. 1205.6497 | abstract | NNLO stability; stability excluded at 98% CL for M_h < 126 |
| Buttazzo et al. 1307.3536, JHEP 12 (2013) 089 | §3, §§4.1-4.3, §§5.1-5.3, §6.3, App. A, B (TeX v4) | matching fits, 3-loop RGEs, lambda_eff, stability formula, lifetime formulas |
| Bednyakov-Kniehl-Pikelner-Veretin 1507.08833, PRL 115 (2015) 201802 | criterion, inputs, results (TeX) | lambda = beta = 0 criterion; M_t,cri = 171.44 +- 0.30 (+0.17 -0.36); 1.3 sigma vs MC mass then |
| Hiller-Höhne-Litim-Steudtner 2401.08811, PRD 110 (2024) 115017 | table and §2 (TeX) | PDG 2024: stability needs M_t < 171.10 (eff.) or alpha_s > 0.1213; MC mass -5.1 sigma one at a time, < 2 sigma (xsec) and 3 sigma (MC) combined; factor 2-3 better errors give 5 sigma |
| Andreassen-Frost-Schwartz 1707.08124; Chigusa-Moroi-Shoji 1707.09301 | abstract | decay rates (§2.3) |
| Arkani-Hamed-Dubovsky-Senatore-Villadoro 0801.2399 | abstract | the decay window and eternal inflation |
| PDG 2025 (pdg-2025-v0.2.3.sqlite, release 2026-05-14; QCD review §9.4.8, rev. Aug 2025) | database values; section | inputs (§5.1) |
| Defranchis et al. 2503.18713; FCC Feasibility Study Report vol. 1, 2505.00272 | abstract; sections (TeX grep) | future precision (§5.6) |
| Senovilla 1303.1408; Deruelle-Sasaki-Sendouda 0711.1150 | abstract | f(R) junction conditions (§1.5) |
| Giudice-Kehagias-Riotto 1907.05370; Kaloper-Westphal 1907.05837 | abstract | q-Higgs couplings (§1.4) |
| Coleman-De Luccia 1980; Guth-Weinberg 1983; Brown-Teitelboim 1988; Henneaux-Teitelboim 1989; Isidori et al. 0712.0242; Elias-Miró et al. 1203.0237; Lebedev 1203.0156 | [R] | as cited |

---------------------------------------------------------------------------------------------------------------------

## 5. The Standard Model relation (`code/sm_rge.py`, `code/sm_criticality.py`, `results/sm_criticality.json`)

### 5.1 Inputs and provenance

| Input | Value | Source (read level) | Role |
|---|---|---|---|
| M_h | 125.20 +- 0.11 GeV (S = 1.36) | PDG 2025 database S126M | measured, retrospective |
| M_t (direct, MC-template) | 172.56 +- 0.31 GeV (S = 1.6) | PDG 2025 Q007TP | measured; MC mass treated as pole mass with optional +-0.3 GeV theory error (Buttazzo's choice) |
| M_t (pole, from cross sections) | 172.4 +- 0.7 GeV | PDG 2025 Q007TP4 | alternative |
| alpha_s^(5)(M_Z) | 0.1180 +- 0.0009 | PDG 2025 QCD review eq. (9.25); FLAG 2024: 0.1183 +- 0.0007 | measured |
| M_W | 80.3692 +- 0.0133 GeV | PDG 2025 S043M | enters g2, gY; effect < 1e-3 GeV |
| y_b(M_t), y_tau(M_t) | 0.0156, 0.0100 | [E] from m_b(m_b), m_tau | effect 1e-4 GeV |
| SM valid to M_Pl = 1.2209e19 GeV (also 2.435e18) | | **supplied [P]** | cut-off |
| Matching theory errors | delta lambda = 3e-4, delta y_t = 5e-4 | [Q] 1307.3536 §3 | theory error |

### 5.2 Method [C]

1. MSbar couplings at mu = M_t from the NNLO fits of Buttazzo et al. 1307.3536 §3 [Q]. Masses are in GeV and
   a3 = alpha_s(M_Z).
   - lambda = 0.12604 + 0.00206 (M_h - 125.15) - 0.00004 (M_t - 173.34)
   - y_t = 0.93690 + 0.00556 (M_t - 173.34) - 0.00042 (a3 - 0.1184)/0.0007
   - g3 = 1.1666 + 0.00314 (a3 - 0.1184)/0.0007 - 0.00046 (M_t - 173.34)
   - g2 = 0.64779 + 0.00004 (M_t - 173.34) + 0.00011 (M_W - 80.384)/0.014
   - gY = 0.35830 + 0.00011 (M_t - 173.34) - 0.00020 (M_W - 80.384)/0.014

   The potential is V = -m^2|H|^2/2 + lambda|H|^4.
2. RG evolution with the **full 3-loop** beta functions of g1, g2, g3, y_t and lambda, transcribed from Buttazzo et al.
   App. A: Mihaila-Salomon-Steinhauser, Chetyrkin-Zoller and Bednyakov et al., credited there. The 4-loop pure-QCD
   term of g3 is included. y_b and y_tau run at two loops. Integration uses DOP853 with rtol 1e-10.
3. Criteria:
   - **'msbar'**: min over mu <= M_Pl of lambda(mu) = 0, i.e. lambda(mu*) = beta_lambda(mu*) = 0 (Bednyakov et al.;
     gauge independent). This is the degenerate-second-minimum condition: **criticality = MPP = boundary of absolute
     stability**.
   - **'eff'**: min over h <= M_Pl of lambda_eff(h) = 0, with the Landau-gauge 1- and 2-loop effective potential of
     Buttazzo et al. App. B. The 2-loop part is in their lambda = 0 approximation.

### 5.3 Validation [C vs Q]

| Quantity | This code | Published |
|---|---|---|
| g1, g2, g3, y_t, lambda at M_Pl (inputs 173.34, 125.15, 0.1184, 80.384) | 0.61537, 0.50546, 0.48731, 0.38248, **-0.01437** | 0.6154, 0.5055, 0.4873, 0.3825, **-0.0143** (1307.3536 §4.1) |
| d lambda(M_Pl)/d M_t; per 0.0007 in alpha_s; d/d M_h | -0.00652; +0.00185; +0.00286 | -0.0066; +0.0018; +0.0029 |
| log10 Lambda_V; Lambda_I/Lambda_V; lambda = 0 crossing | 9.497; 12.7; 5.7e9 GeV | 9.5; ~13; ~2 Lambda_V ~ 6e9 GeV |
| lambda(M_Pl) by RGE order 1 / 2 / 3 loops | -0.0337 / -0.0151 / -0.0144 | — (convergence) |
| M_t,crit at (125.7, 0.1185), 'msbar' | 171.38 | 171.44 (Bednyakov et al.; same criterion, own matching) |
| M_t,crit at (125.20, 0.1180) | 170.94 ('msbar'), 171.00 ('eff') | 171.04 (tree resummed), 171.10 (eff.) (Hiller et al., PDG 2024) |
| M_t,crit at (125.15, 0.1184) | 171.07 / 171.13 | Buttazzo eq. (stability) gives 171.12; their quoted "M_t < 171.53 +- 0.42" does **not** follow from their own linear formula. This is an internal 0.4 GeV inconsistency in the source, recorded and not resolved |

### 5.4 Results [C]

**Stability bound = criticality line** (PDG 2025 point; 'msbar'; 'eff' in brackets):

    M_t,crit = 170.94 (171.00) + 0.495 (M_h - 125.20) + 0.353 (alpha_s - 0.1180)/0.0009   GeV,
    mu* = 4.7e17 GeV (3.0e18 for 'eff'); identical with the cut at the reduced Planck mass.
    Inverse (criticality "prediction" of M_h):
    M_h,crit = 128.54 (128.41) + 2.02 (M_t - 172.56) - 0.71 (alpha_s - 0.1180)/0.0009   GeV.

**Theory error on M_t,crit:**

| Source | Shift |
|---|---|
| delta lambda = 3e-4 | +0.072 GeV |
| delta y_t = 5e-4 | -0.087 GeV |
| **combined matching error (quadrature)** | **0.113 GeV** |
| criterion ('eff' - 'msbar') | +0.064 GeV, half used as an error |
| 2-loop vs 3-loop RGE | -0.109 GeV (recorded; the 4-loop step is expected to be several times smaller [E]) |
| y_b, y_tau | 1e-4 GeV |
| M_W | 4e-4 GeV |

**Froggatt-Nielsen double condition** (lambda = beta_lambda = 0 at a fixed mu*, alpha_s = 0.1180):

| mu* | M_t | M_h |
|---|---|---|
| 1.22e19 GeV | 174.20 | 132.03 |
| 2.435e18 GeV | 172.61 | 128.64 |
| 1e18 GeV | 171.71 | 126.77 |

Notes on this table:
- FN 1995 quote (173 +- 5, 135 +- 9) at two loops, with mu* known only to an order of magnitude.
- The fits are extrapolated in M_h, with an error of ~0.05-0.2 GeV [E].
- With mu* at a Planck mass, M_h misses the measured 125.20 by 3.4-6.8 GeV (retrospective). Only the single-condition
  line with mu* ~ 5e17 GeV passes through the measured M_h.

**SM potential at PDG central values** (M_t = 172.56):

| Scale | lambda | lambda_eff |
|---|---|---|
| zero crossing of lambda | 4.2e10 GeV | |
| 1e17 GeV | -0.0105 | -0.0097 |
| 2.435e18 GeV | -0.0104 | -0.0100 |

V_eff = e^{4 Gamma} lambda_eff h^4/4 is -2.0e65 GeV^4 at 1e17 GeV, -7.5e70 at 2.435e18 GeV and -4.8e73 at 1.22e19 GeV.
These numbers feed §3.1.

### 5.5 Position of the measured values: retrospective tension

z = (M_t,meas - M_t,crit)/sigma_tot. sigma_tot combines in quadrature: sigma(M_t); the propagated sigma(M_h) and
sigma(alpha_s); the matching error 0.113 GeV; optionally the criterion error 0.032 GeV and a 0.3 GeV MC-to-pole error.

| M_t input | Criterion | Delta (GeV) | sigma_tot (GeV) | z | one-sided p |
|---|---|---|---|---|---|
| direct 172.56 +- 0.31 | msbar | 1.62 | 0.49 | **3.3** | 4e-4 |
| direct, + 0.3 GeV MC-pole | msbar | 1.62 | 0.57 | **2.8** | 2.3e-3 |
| direct | eff | 1.56 | 0.49 | 3.2 | 7e-4 |
| xsec pole 172.4 +- 0.7 | msbar | 1.46 | 0.79 | **1.8** | 0.03 |
| xsec pole | eff | 1.40 | 0.79 | 1.8 | 0.04 |

In M_h units, M_h would have to be 128.5 GeV instead of 125.20 +- 0.11. The meaningful statement is the joint one
above, dominated by M_t and alpha_s. Consistent with Hiller et al. 2024 (< 2 sigma with the cross-section mass, ~3
sigma with the MC mass) and with Buttazzo et al. 2013 (2.8 sigma).

**Status of the relation under readings (a), (c)-worldline and (b)-QV1:** the measured values are on the metastable
side at 1.8 to 3.3 sigma, depending on which top mass is used. Per the stage-2 rules this is reported as a
**retrospective tension of that version (QV1 + reading + SM to M_Pl)**. It is not a decisive failure at present
errors, and it is not repaired here.

### 5.6 Held-out tests (data do not yet exist) [Q][R][C]

| Scenario | Inputs | sigma_tot | z if central values stay |
|---|---|---|---|
| FCC-ee | m_t to ~40 MeV (6.8 MeV exp + ~35 MeV N3LO theory, PS mass: 2503.18713; FSR: O(50) MeV theory); alpha_s to 0.1% (+-0.00012, FSR vol. 1); m_h to 4 MeV (FSR vol. 1); matching errors as today | 0.13 GeV | **12** (direct), 11 (xsec) |
| FCC-ee, matching errors halved | as above | 0.09 GeV | 18 |
| HL-LHC + lattice | m_t ~0.2 GeV [R, order of magnitude]; lattice alpha_s halved (+-0.00045, FSR vol. 1: "factor of two in the next decade"); m_h +-0.02 GeV [R] | 0.29 GeV | **5.6** |

These z values illustrate reach. They are not predicted outcomes. The decisive held-out measurement is a theoretically
clean short-distance top mass from the e+e- threshold scan, together with alpha_s from the Z pole or lattice QCD. At
that point the SM matching errors (y_t and lambda at M_t) become the limiting theory input. AFS 1707.08124 (abstract):
excluding stability at 3 sigma needs delta M_t < 250 MeV or delta alpha_s < 0.00025.

---------------------------------------------------------------------------------------------------------------------

## 6. Verdict

1. **Does QV1 imply a relation beyond the calibration target?** Not as stated: with several matter vacua, "the final
   vacuum" in P_q is ambiguous (§3). It does once a reading is fixed:
   - Under **(a) global minimum** or **(c) with the worldline measure** (which reduces to (a) whenever Gamma > 0), and
     under **(b) with uncalibrated QV1**, QV1 implies that our vacuum is the global minimum to within rho_Lambda.
     With the SM valid to the Planck scale (supplied) this is the absolute-stability inequality
     **M_t <= 170.94 + 0.495 (M_h - 125.20) + 0.353 (alpha_s - 0.1180)/0.0009 GeV (+-0.11 th)**, a one-sided relation
     among measurable SM quantities. It does not use rho_Lambda's value, only its smallness [C][D].
   - Under **(d)** (an added assumption: MPP, or "P_q for every reachable final vacuum"), it becomes the equality
     M_t = M_t,crit (M_h,crit = 128.5 GeV at the measured M_t).
   - Under **(b)/(c) with QV1-cal and the volume measure**, no SM relation follows.
2. **Retrospective status.** The measured (M_t, M_h, alpha_s) are on the metastable side at 3.3 sigma (direct M_t),
   2.8 sigma (with the MC-to-pole error) and 1.8 sigma (cross-section pole mass). Readings (a)/(c)-worldline therefore
   carry a retrospective tension at this level, and the equality form (d) carries the same tension. Without the
   relation, reading (a) with the SM to the Planck scale at central values gives Lambda_ours ~ 3e117 rho_Lambda; the
   relation is exactly the condition that avoids this.
3. **Escape.** New physics below the instability scale (lambda = 0 at 4.2e10 GeV in MSbar; Lambda_V ~ 1e9.5-1e10 GeV)
   that makes our vacuum the global minimum: for example a scalar threshold ([R] Elias-Miró et al. 1203.0237, Lebedev
   1203.0156) or the Higgs-portal singlets of Hiller et al. 2401.08811 [Q]. Reading (a) further requires the absence
   of deeper vacua in the UV completion. With the escape, the relation constrains the extension, not the SM, and loses
   predictive force unless the extension is specified. That would be a new version (call it QV2-vacua); QV1's result
   above is kept.
4. **Novelty, and credit.**
   - The SM relation is Froggatt-Nielsen's criticality (hep-ph/9511371) and the stability bound of Degrassi,
     Buttazzo and Bednyakov, recomputed and validated here. The equality form is also in Klinkhamer 1302.1496.
   - "Cancellation at the global minimum, so a false vacuum's Lambda is its splitting" is AHHKM's premise and KV's
     global-equilibrium choice.
   - That q-theory's 4-form does not provide FN's fixed extensive quantity for Higgs phases is consistent with
     Volovik's own counting (hep-ph/0309144).
   - What this track adds:
     - the action-level formula Lambda_i = F(V_i - eps*(mu)) with K' corrections, and the junction behaviour of q
       (standard f(R) junction conditions in q-language);
     - the observation that QV1's own final-state condition, applied with a worldline measure in the presence of
       nucleation, has a unique self-consistent fixed point, the true vacuum. That turns QV1 into the
       stability inequality without adding MPP.

     I found this specific chain in no searched paper (§4). It is bookkeeping on known physics, not a new mechanism.
5. **GRUT's contribution.** None that is physical. GRUT's "the whole has no exterior" does not choose between the
   worldline and volume measures, nor between "global" and "final", so it does not decide the reading. "The vacuum the
   whole ends in" is a final boundary condition, a kind FORMAL_CORE §7 lists as outside its table (ledger §1.3). THM
   2.15 restates the point: conserved mu must be prepared, and preparing it correctly requires knowing the whole's
   decay endpoint. **Nothing is GRUT-native.**
6. **Explicit added assumptions.**
   - Reading (a): global minimum.
   - Reading (b): our vacuum.
   - The measure in (c).
   - MPP in (d).
   - SM validity to M_Pl (supplied).
   - Treating the PDG direct mass as the pole mass (with an optional 0.3 GeV error).

---------------------------------------------------------------------------------------------------------------------

## 7. Limits

- Classical q-sector. Bubble nucleation is treated as GR decay at fixed mu (justified for m_q R >> 1, §1.5). No
  quantum q-sector, no Brown-Teitelboim membranes.
- SM inputs:
  - The matching fits are those of 2013, linearised; I did not redo NNLO matching. The 0.4 GeV internal inconsistency
    of 1307.3536 (171.53 against its own formula) is noted, not resolved.
  - The 'eff' criterion uses the 2-loop lambda_eff in the lambda = 0 approximation, with real parts taken where
    y_t^2/2 < g2^2/4. That region is outside the relevant one.
- The MC-mass interpretation is the dominant non-perturbative issue. Its size (0.3-0.5 GeV) is [R]/[Q], not
  computed.
- HL-LHC precisions are recalled, order of magnitude only.
- Vacuum-decay rates are quoted at abstract level. The two quoted rates differ by many orders of magnitude; only the
  qualitative conclusions of §2 (P_past ~ 0, Gamma/H^4 << 0.716) are used, and they hold over the whole range.
- §2.5 (CDL quench) is a thin-wall, Planckian-field estimate.
- Volume vs worldline measure (§3.3) is not settled by the action. Choosing the worldline measure is an assumption.

## 8. Files

- `code/sm_rge.py`: matching fits, 3-loop RGEs, Higgs anomalous dimension, lambda_eff (1- and 2-loop).
- `code/sm_criticality.py` -> `results/sm_criticality.json`: validation, bound, inverse relation, FN double
  condition, tension, held-out projection, potential depths, CDL estimate.
- `code/vacua_q_sector.py` -> `results/vacua_q_sector.json`: §1.3 series and exact slope, numerical static roots,
  §1.4 coupled case, §3.4 coexistence/Maxwell, §1.5 linear wall.
- `code/nucleation.py` -> `results/nucleation.json`: §2 (V4, rates, crunch times).
- `code/run_all.py`: runs all three scripts.
- `results/sources.json`: sources, read levels, search record.
