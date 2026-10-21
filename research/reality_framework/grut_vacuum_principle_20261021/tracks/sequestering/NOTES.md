# Track: vacuum-energy sequestering against the four-question gate

Packet `grut_vacuum_principle_20261021`, track `sequestering`, 21 October 2026. Written under
`BRIEF.md`. Labels: [D] derived here, [C] computed here (code below), [Q] quoted from a paper at the
stated read level, [E] order-of-magnitude estimate, [A] argued. No new theory, no empirical validation
and no confirmatory test is claimed. Every comparison with Planck, DESI or supernova data is
**retrospective**.

**Bottom line.** Sequestering is the existing proposal closest to "the whole universe fixes the
gravitating vacuum term". It escapes the energy-origin obstruction. It does **not** select the matter
shift c. Instead, it changes the gravitational equations so that gravity is also blind to c. Credit
goes to Kaloper and Padilla, with Stefanyszyn and Zahariade for the local version. The proposal passes
gate items 1 and 3 for matter loops; only the Gauss–Bonnet version also covers graviton vacuum
loops. It fails item 2. In the local and Gauss–Bonnet versions, the residual is a free flux-ratio
constant. In the global version, it is a functional of a whole history whose inputs are supplied,
and those inputs include a logarithm of the crunch cutoff. Item 4 also fails. In its simplest
published global realization, a linear-potential collapse trigger, the one relation the theory
imposes (<R> = 0) cannot be met for any slope I scanned (s = m^3/(M H0^2) from 0.3 to 3) or any
closed curvature Omega_k from -0.95 to 0. Open curvature only makes <R> more negative. This reproduces
Ferreira–Avelino (2015) [C]. The other versions make no independent prediction. Nothing here is
GRUT-native.

---------------------------------------------------------------------------------------------------

## 1. Sources and read levels

Downloads are kept in the scratchpad only (`.../scratchpad/seq_dl/`); no paper text is copied here.

| Ref | Identifier | Read level |
|---|---|---|
| Kaloper & Padilla, "Sequestering the SM vacuum energy", PRL 112, 091304 (2014) | arXiv:1309.6562 | full |
| Kaloper & Padilla, "VES: the framework and its cosmological consequences", **PRD 90, 084023 (2014)** | **arXiv:1406.0711** | §3, §4.1–4.3, §5, §7, §8, §9 read; §2, §6 skimmed |
| Kaloper & Padilla, "'The End'", PRL 114, 101302 (2015) | arXiv:1409.7073 | full |
| Kaloper, Padilla, Stefanyszyn, Zahariade, "Manifestly local VES", PRL 116, 051302 (2016) | arXiv:1505.01492 | full |
| Kaloper & Padilla, "VES and graviton loops" (omnia sequestra, OS), PRL 118, 061303 (2017) | arXiv:1606.04958 | full |
| Avelino, "VES and cosmic dynamics", PRD 90, 103523 (2014) | arXiv:1410.4555 | full |
| Ferreira & Avelino, "The fate of a Universe driven by a linear potential" | arXiv:1508.00631 | full |
| Coltman, Li, Padilla, "Cosmological consequences of Omnia Sequestra" | arXiv:1903.02829 | §2, §3 (residual estimate), conclusions |
| Lawrence & Seery, "UV sensitivity of the cosmological sequester", PRD 103, 064018 (2021) | arXiv:2007.07058 | §II C–§III D, conclusions |
| Lombriser, "Late-time acceleration ... ultimate collapsed structures", JCAP 09 (2019) 065 | arXiv:1805.05918 | §2, §3.4 |
| Lombriser, "On the cosmological constant problem", PLB 797, 134804 (2019) | arXiv:1901.08588 | abstract |
| Kaloper, Padilla, Stefanyszyn, "Sequestering effects on and of vacuum decay", PRD 94, 025022 (2016) | arXiv:1604.04000 | abstract |
| D'Amico, Kaloper, Padilla, Stefanyszyn, "An étude on global vacuum energy sequester" | arXiv:1705.08950 | abstract |
| Bufalo, Klusoň, Oksanen (canonical/path integral, local VES) | arXiv:1606.01282 | abstract |
| Svesko & Zahariade (constraint structure, local VES) | arXiv:1812.11625 | abstract |
| El-Menoufi, Nagy, Niedermann, Padilla (quantum corrections, monodromy), CQG 36, 215014 (2019) | arXiv:1903.07612 | abstract |
| Khoury, Muntz, Padilla, "A lapse in the cosmological constant problem" (2026) | arXiv:2604.08659 | abstract |
| DESI DR2 Results II, PRD 112, 083515 (2025) | arXiv:2503.14738 | Table 5 (section) |
| Andrei, Ijjas, Steinhardt, PNAS 119, e2200539119 (2022) | arXiv:2201.07704 | abstract |
| Luu, Qiu, Tye, "The lifespan of our Universe" (2025) | arXiv:2506.24011 | abstract |
| Van Raamsdonk & Waddell (linear-potential fit to DES-SN5YR + DESI DR2) | arXiv:2505.10662 | abstract |
| Berghaus, Kable, Miranda, PRD 110, 103524 (2024) | arXiv:2404.14341 | abstract |

**Identifier correction.** The task text cites "PRD 89 (2014) 084023, arXiv:1310.0464" for global
sequestering. In fact arXiv:1310.0464 is *The Effective Field Theory of Large Scale Structures at
Two Loops*, an unrelated paper (abstract checked). The sequestering paper is **arXiv:1406.0711, PRD
90, 084023 (2014)**, and it is the one used here.

**Prior GRUT work, not repeated as new.** These points are already in the record:
`grut_closure_fork_20261019/DERIVATION.md:197-215` gives the FRW identity for <R> = 0 and the
conditional need for k > 0, and it flags stiff endpoints as open. `REPORT.md:121-127` says that the
global and local formulations must not be conflated. `CLOSURE_SELECTION_RESPONSE.md:166-168`
withdraws the claim that every sequestering formulation forces positive curvature. Section 5 below
closes the stiff-endpoint item numerically for one realization. It does not reinstate the withdrawn
claim.

## 2. Conventions

The signature is (-,+,+,+). M is the reduced Planck mass, M^2 = 1/(8 pi G). The spacetime average is
<Q> = ∫d^4x sqrt(-g) Q / ∫d^4x sqrt(-g). The trace is tau = tau^a_a = -rho + 3p. The protected-sector
stress tensor is split as T^mu_nu = -V_vac delta^mu_nu + tau^mu_nu, where V_vac is the regularized
vacuum energy at any loop order and tau holds the local excitations.

The residual enters the field equation as M^2 G^mu_nu = tau^mu_nu - Lambda_eff delta^mu_nu. It
therefore enters Friedmann as an energy density: 3M^2(H^2 + k/a^2) = rho + Lambda_eff. A positive
Lambda_eff is de Sitter-like. (The opening paragraph of 1406.0711 §4 writes "-<tau>/4" once. Its eqs.
(19), (44) and (45) fix the sign used here [Q].)

The target is rho_vac,total / M^4. Its observed value is rho_DE0/M^4 = 3 Omega_DE (H0/M)^2 =
7.16e-121 for h = 0.674, Omega_DE = 0.685, M = 2.435e18 GeV [C]. That is rho_DE0 = (2.24 meV)^4.

## 3. The versions

### 3.1 Global sequestering (KP 1309.6562; 1406.0711)

* **Action [Q].** S = ∫d^4x sqrt(-g) [ M^2 R/2 - Lambda - lambda^4 L(lambda^-2 g^{mu nu}, Phi) ] +
  sigma( Lambda / (lambda^4 mu^4) ). The global variables are Lambda (the bare cosmological
  constant) and lambda (the matter-to-Planck scale ratio, m_phys = lambda m). sigma is odd,
  differentiable and not a pure log, and it sits **outside** the integral. mu is a phenomenological
  scale. The protected sector, including its regulator (e.g. Pauli–Villars), must couple to
  g~ = lambda^2 g. Then every vacuum contribution scales as lambda^4.
* **Field equations [Q].** sigma'/(lambda^4 mu^4) = ∫sqrt(-g). 4 Lambda sigma'/(lambda^4 mu^4) =
  ∫sqrt(-g) T. Together these give Lambda = <T>/4, and
  **M^2 G^mu_nu = T^mu_nu - (1/4) delta^mu_nu <T>**. This is all ten equations, trace included; there
  is no hidden integration constant, unlike unimodular gravity. With T = -V_vac + tau:
  M^2 G^mu_nu = tau^mu_nu - (1/4) delta^mu_nu <tau>. The trace and average give **<R> = 0**.
* **Residual [Q].** Lambda_eff = <tau>/4 (1406.0711 eq. 19). For one fluid with fixed w,
  Lambda_eff = -(1-3w) ∫dt a^3 rho / (4 ∫dt a^3) (eq. 45). It is negative if w < 1/3 dominates the
  4-volume-weighted history and positive if w > 1/3 does.
* **Spacetime requirement [Q].** ∫sqrt(-g) must be finite, otherwise lambda → 0 and the protected
  sector becomes massless. That requires compact spatial sections plus a bang and a crunch. Infinite
  universes are solutions, but they are not our world. For a large, old, nearly FRW history the
  authors argue k = +1 (eq. 47 discussion). The endpoint term [a^3 H] is then a cutoff-dependent
  boundary datum.
* **What is sequestered [Q].** Protected-sector (matter) vacuum energy at every loop order, provided
  the regulator couples to g~. Graviton loops are not sequestered (1606.04958 eqs. 5–7). The
  renormalization M^2 → M^2 + lambda^2 M_UV^2 is harmless because <R> = 0 (1406.0711 eqs. 56–60).
  Phase transitions after t* are suppressed by ∫_0^{t*} a^3 dt / ∫a^3 dt (eqs. 49–52).
* **Fixed or free [Q/A].** The residual is not a free integration constant. It is fixed a posteriori
  by the complete history, which KP call a bootstrap or gap equation and a "boundary condition at
  future infinity". That history contains supplied inputs: initial data, the collapse sector and the
  crunch cutoff. In GR terms the theory adds exactly **one relation** to the space of GR histories,
  namely <R> = 0 [D, §4].
* **Relation to cosmic history [Q].** Lambda_eff ≈ O(1) M^2 H_age^2 (eq. 40). It is bounded by the
  critical density if the lifetime exceeds about H0^-1. Lambda_eff cannot itself be the present dark
  energy: if it dominated, collapse would be impossible. So w_DE ≈ -1 must be transient.
* **Consequences claimed [Q].** A future crunch, a spatially closed and compact universe (with
  possible topology signatures), and transient acceleration. No quantitative collapse time is given
  in the general framework.

### 3.2 Global sequestering with a linear trigger ("The End", 1409.7073)

* **Added structure [Q].** A protected scalar with V = m^3 phi. The linear form is protected by the
  shift phi → phi + C, which is absorbed by Lambda. The slope is technically natural, and it is
  **chosen** to be m^3 ≈ M H0^2. Numerically M H0^2 = 5.0e-39 eV^3, or (1.7e-13 eV)^3 [C].
* **Constraint [Q].** In the gauge V = m^3 phi: <phi_dot^2> - 4 m^3 <phi> + <3p - rho> = 0
  (their eq. 9). KP argue that this picks a trans-Planckian phi_in, at least one e-fold of
  acceleration before turnaround, and turnaround at t ≈ H0^-1. They offer this as an answer to
  "why now".
* **Critique [Q].** Avelino (1410.4555): for the linear model the crunch is kinetic-dominated, so
  Lambda diverges unless it is cut off. Ferreira & Avelino (1508.00631): once phi_in is fixed by the
  present Omega_DE, the constraint cannot be met for slopes and curvature allowed by observation.
  Solutions exist only if curvature dominates before the scalar. My §5 calculation reproduces this
  independently.

### 3.3 Local sequestering (KPSZ 1505.01492)

* **Action [Q].** S = ∫sqrt(-g)[ kappa^2(x) R/2 - Lambda(x) - L_m(g, Phi) ] +
  ∫[ sigma(Lambda/mu^4) F + sigma^(kappa^2/M^2) F^ ], with F = dA and F^ = dA^ two 4-forms.
  The functions sigma and sigma^ are smooth and nonlinear. The gauge symmetries of the 4-forms force
  Lambda and kappa^2 to be constant on shell. The action is local and additive, so a Hamiltonian and
  a path integral exist (1606.01282, 1812.11625, abstract level).
* **Equations [Q].**
  - kappa^2 G = (∇∇ - g □) kappa^2 + T - Lambda g.
  - (sigma'/mu^4) F = sqrt(-g) eps.
  - (sigma^'/M^2) F^ = -(1/2) R sqrt(-g) eps.
  - Hence Lambda = <T>/4 + ΔΛ, with ΔΛ = kappa^2 <R>/4 = -(mu^4 kappa^2 sigma^')/(2 M^2 sigma')
    × ∫F^/∫F.
  - Field equation: kappa^2 G^mu_nu = T^mu_nu - (1/4) delta <T> - ΔΛ delta.
  - **Lambda_eff = <tau>/4 + ΔΛ.**
* **Spacetime [Q].** A finite 4-volume is not required. The flux ratio stays finite in an infinite
  universe. With eternal expansion and NEC sources, <tau> → 0, so Lambda_eff = ΔΛ.
* **What is sequestered [Q].** Matter loops only, as in the global version.
* **Fixed or free [Q].** ΔΛ is "a priori arbitrary", a ratio of integration constants. It is
  radiatively stable once measured, and it is to be fitted to observation.
* **Consequences.** None beyond GR + Lambda locally: local solutions coincide with GR. Phase
  transitions in descendant vacua are suppressed by volume ratios (1604.04000, abstract).

### 3.4 Omnia sequestra (KP 1606.04958; Coltman–Li–Padilla 1903.02829)

* **Action [Q].** S = ∫sqrt(-g)[ M^2 R/2 + theta(x) R_GB - Lambda(x) - L_m ] +
  ∫[ sigma(Lambda/mu^4) F + sigma^(theta) F^ ]. The global limit gives:
  - sigma' c/mu^4 = ∫sqrt(-g), and sigma^' c^ = -∫R_GB sqrt(-g);
  - so <R_GB> = -mu^4 sigma^' c^/(sigma' c).
* **Residual [Q].** Field equation: M^2 G = T - <T> g/4 - ΔΛ g. The residual is Lambda_res =
  <tau>/4 + ΔΛ, where
  ΔΛ^2 = (3M^4/8)[ <R_GB> - <W^2> + (4/M^4)<(tau_mu_nu - tau g_mu_nu/4)^2>
  - (1/(6M^4))(<tau^2> - <tau>^2) ].
  ΔΛ has two roots (±), and both are flux-set and free.
* **What is sequestered [Q].** Matter and graviton vacuum loops. The theta shift symmetry is broken
  only by topological terms. Lawrence–Seery find that graviton dressing of the Lambda–volume
  coupling is absorbed here.
* **Costs [Q].** The historic integrals are UV-sensitive unless the universe keeps expanding. For a
  Planckian cutoff that means a_max/a0 ≳ 1e40, about 92 more e-folds. The QCD transition needs about
  23 e-folds. A vanishing flux ratio forces hyperbolic spatial slices. Lawrence–Seery also note that
  the constraint on <R_GB> roughly requires mu ~ H_age if sigma' and sigma^' are O(1).
* **Consequences.** None locally. The required long future runs opposite to the global version's
  imminent crunch. A near-term crunch would make the OS residual cutoff-sensitive (unnatural), not
  strictly forbidden [A].

### 3.5 Other follow-ups (abstract level unless stated)

* Lombriser (1805.05918, §2 and §3.4 read): averaging over "ultimate collapsed structures", with a
  uniform prior on our location in collapse, gives Omega_Lambda = 0.697. The 2019 PLB paper
  (1901.08588) uses a Planck-mass variation and gives 0.704. Both add a halo-patch averaging
  assumption and a location prior, and both are matched **retrospectively** to known data. Neither
  is a confirmatory test.
* D'Amico et al. (1705.08950): global-multiplier proposals that do not treat the counterterm as
  dynamical need retuning order by order.
* Khoury–Muntz–Padilla (2604.08659): a global constraint from a 5d lapse; residual and tests not
  read.

## 4. Relation to the energy-origin obstruction [D]

Under the matter shift L → L + C (equivalently T^mu_nu → T^mu_nu - C delta^mu_nu), <T> shifts to
<T> - 4C. The right-hand side **T - (1/4) delta <T> is exactly invariant**. Gravity therefore becomes
as blind to C as the closure-only quantum data are to H → H + cI. The global variable Lambda
absorbs C (Lambda → Lambda - lambda^4 C), changing the action only by
sigma(Lambda/(lambda^4 mu^4) - C/mu^4) - sigma(Lambda/(lambda^4 mu^4)) ≈ -sigma' C/mu^4. That is a
change in the global constraint, not in the bulk equations.

So sequestering **does not select c**. It makes c unobservable gravitationally. The observable
curvature is then set by a different quantity: <tau>/4 (global) or a flux ratio (local, OS).

* **Q1 answered.** The extra structure is new global or 4-form variables with a non-additive or
  topological coupling. In the global and local versions it also includes a scaling assumption on
  the protected sector and its regulator.
* **The closure feature is neither necessary nor sufficient [D].**
  - *Not sufficient.* In GR on a closed, compact, finite-4-volume spacetime, adding -C ∫sqrt(-g)
    still shifts T by -C g and changes the curvature. A finite whole does not remove the gravitating
    c.
  - *Not necessary.* The local and OS versions cancel V_vac in infinite spacetimes. In the global
    version, finiteness is a consistency requirement (lambda ≠ 0), not the cancelling agent.
* **Single relation [D].** In GR with source tau_GR and no separate constant, M^2 R = -tau_GR. The
  global theory is therefore GR plus the condition <R> = 0. Any constant in a potential is a shift
  of phi.

## 5. Calculation [C]

Code: `code/seq_history.py`. Results: `results/seq_history_results.json`. Runtime is about 1 min.

**Inputs.** These are stated, not predicted: h = 0.674, Omega_m = 0.315 and Omega_r = 4.18e-5/h^2,
which are Planck-2018-like. They lie within 0.3 sigma of the DESI DR2 Table 5 CMB-only LambdaCDM
column, Omega_m = 0.3169 ± 0.0065 and H0 = 67.14 ± 0.47 [Q]. Units are M = H0 = 1. The crunch cutoff
rho_cut is one of M^4 (= 2.87e120 M^2 H0^2), (1 TeV)^4, or rho_s = 1e12 (no tail).

**Part A: LambdaCDM with an instantaneous mirrored collapse.** This is a scale check only, with no
mechanism. Only matter enters tau, and the observed Lambda is identified with Lambda_eff.
- Lambda_eff/rho_Lambda = **-0.408** for collapse at a_c = 1 (Omega_m = 0.315), and -0.396 for
  Omega_m = 0.308.
- |Lambda_eff| = rho_Lambda requires a_c = **0.721** (0.713 for Omega_m = 0.308).
- This reproduces in magnitude the "40 %" and "a = 0.712" quoted by Lombriser (1805.05918 §2) for
  Planck-2015 parameters.
- **Sign:** with KP conventions the residual is negative, because the dust-dominated history has
  w < 1/3. A LambdaCDM history cannot produce the observed positive Lambda this way, and such a
  history would not collapse anyway.

**Part B: the KP linear trigger, solved self-consistently in GR form.**
- *Model.* The protected sector is matter, radiation and phi with V = s phi, where
  s ≡ m^3/(M H0^2). There is no separate constant. phi starts on the radiation-era attractor, and
  phi_in is shot so that rho_phi(a=1) = 3 Omega_DE.
- *Integration.* GR from a = 1e-7, through turnaround, into a kinetic-dominated crunch. The crunch
  is continued analytically: in the w = 1 tail, a^3 H = const = -C, which gives
  ∫a^3 phi_dot^2 dt = C ln(rho_cut/rho_s).
- *Checks.* The identity ∫a^3 R dt = 6[a^3 H] + 6∫(k a - a^3 H^2)dt agrees with -∫a^3 tau dt to
  about 1e-9 relative. The Friedmann drift is 4e-11, and the kinetic fraction at the start of the
  tail is 1 - 1.3e-10.
- *Self-consistency.* Lambda_res ≡ <tau>/4 in this gauge must vanish; equivalently <R> = 0. The
  table below uses Omega_k = 0. CPL fits are taken over 0.5 < a < 1 and are diagnostic only.

| s = m^3/(M H0^2) | w0 | wa (CPL fit) | e-folds to turnaround | turnaround / crunch from now [Gyr] | Lambda_res/rho_DE0: M^4 cut | (TeV)^4 cut | no tail | <R>/H0^2 (M^4 cut) |
|---|---|---|---|---|---|---|---|---|
| 0.3 | -0.997 | -0.005 | 24.6 | 594 / 618 | +4.24 | +2.13 | +0.52 | -34.8 |
| 0.6 | -0.988 | -0.020 | 6.8 | 167 / 184 | +8.45 | +4.23 | +1.01 | -69.5 |
| 1.0 | -0.966 | -0.054 | 2.9 | 71 / 84 | +14.07 | +7.03 | +1.66 | -115.6 |
| 1.5 | -0.924 | -0.122 | 1.6 | 39 / 49 | +21.07 | +10.51 | +2.46 | -173.2 |
| 2.0 | -0.869 | -0.215 | 1.0 | 26 / 35 | +28.07 | +13.99 | +3.26 | -230.7 |
| 3.0 | -0.726 | -0.472 | 0.6 | 15 / 23 | +42.06 | +20.94 | +4.85 | -345.7 |

**Findings [C].**
1. **The residual's scale is the trigger's own scale.** In the late dynamics the only scale is the
   slope: the equations are invariant under t → t/sqrt(s) with rho ∝ s, and earlier epochs are
   diluted. The fit is
   **Lambda_res ≈ m^3 M [0.102 ln(rho_cut/(m^3 M)) + 0.59]**. That is 28.9 m^3 M for a Planck
   cutoff and 14.4 m^3 M for a TeV cutoff. Lambda_res/(s rho_DE0) = 14.1 ± 0.06 across the whole
   grid. The observed vacuum scale comes in through the chosen s ~ O(1). It is not output.
2. **The crunch dominates.** At a Planck cutoff, the w = 1 crunch tail beyond rho = 1e12 M^2 H0^2
   supplies about 88 % of <tau>/4, and kinetic energy overall supplies about 99 %. 1409.7073
   (footnote 4) says the cut-off singular region is small in big universes. Here it is not
   subdominant: it exceeds the rest of the history by a factor of about 7. The stiff tail begins
   just after turnaround, so its weight is the turnaround weight times about
   6 ln(a_*/a_cut) ≈ 250. Halving ln(rho_cut) halves Lambda_res: the dependence on the cutoff is
   logarithmic, not a power.
3. **No consistent solution.** <R> ≈ -115.5 s H0^2 ≠ 0. It is insensitive to curvature: across the
   scan Omega_k ∈ [-0.95, 0] (20 points) for s = 0.6, 1, 1.5 and 2, <R> changes by less than 0.2 %
   and never changes sign. No root exists. Open curvature (k < 0) adds a negative 6k/a^2 to R and
   cannot help [A]; it was spot-checked at Omega_k = +0.0023. Closed curvature (k > 0) can help,
   but its effect is diluted by a_max^-2. The stiff endpoint term [a^3 H]_f = -C < 0 also has the
   wrong sign. <R> → 0 only as s → 0, which is
   infinite 4-volume and lambda → 0, excluded in the global theory. Within the slope scan
   s ∈ [0.3, 3] this is [C]; extrapolation beyond it rests on the scaling argument [A]. This closes,
   for this realization, the stiff-endpoint item left open in `DERIVATION.md:210-211`, and it
   independently confirms Ferreira–Avelino (credited). It is **not** a statement about every
   collapse mechanism.
4. **Phase transitions are suppressed [C/E].** I used the s = 1 history with a sudden transition
   and order-of-magnitude ΔV. The electroweak shift is -2.9e-25 rho_DE0 (ΔV = m_h^2 v^2/8,
   weight 6.2e-80). The QCD shift is -1.4e-20 rho_DE0 (ΔV = (0.2 GeV)^4, weight 2.1e-64). These
   are consistent with KP eq. 52.
5. **Graviton-loop instability estimate [E].** In the global and local versions, the unsequestered
   one-graviton vacuum term ~ M_c^6/M^2 exceeds rho_DE0 unless the sequestered-sector cutoff
   M_c ≲ 23 MeV, or ≲ 0.12 GeV with a (16π^2)^-2 loop factor. At M_c = 1 TeV the ratio is 6.7e27.
   KP 2017 say "thirty orders" [Q]. Lawrence–Seery add that Lambda–volume dressing alone exceeds
   the observed value if the heaviest sequestered scale is ≳ 10 TeV [Q]. OS removes both [Q].

**What a test would have to measure.**
- For the global version with a specified collapse sector:
  1. the present slope of the dark-energy potential, through w(z) and dw/dz;
  2. the curvature sign and magnitude;
  3. a check that the inferred **whole** history, including the future crunch and a stated cutoff,
     satisfies <R> = 0.
  The residual depends on the future and on ln(rho_cut), so no present-epoch measurement fixes it by
  itself. A test is only as sharp as the collapse model and cutoff that are frozen in advance.
- For the linear trigger, this test already fails on exposed data (item 3 above).
- For the local and OS versions, no measurement tests the residual: it is a fitted constant.

## 6. Already-exposed data (retrospective only)

* **Planck 2018 inputs.** Omega_m, h and Omega_r are inputs, so they are not held out.
* **DESI DR2 Table 5 [Q].**
  - LambdaCDM + Omega_K, DESI+CMB: 10^3 Omega_K = 2.3 ± 1.1, slightly open. Closed universes have
    Omega_K < 0.
  - w0wa + Omega_K, DESI+CMB+Pantheon+: 10^3 Omega_K = 1.1 ± 1.3, w0 = -0.853 ± 0.057,
    wa = -0.54 ± 0.22.
  - w0wa, DESI+CMB+DESY5: w0 = -0.752 ± 0.057, wa = -0.86 (+0.23, -0.20).
* **Supernova samples.** Pantheon+, Union3 and DES-Y5 are inside those combinations.
* **Linear-potential fits [Q, abstracts].** Berghaus et al. 2024 and Van Raamsdonk & Waddell 2025
  already fit linear potentials to these data.
* **Linear trigger vs DESI.** The s ≈ 2–3 rows land in the DESI quadrant (w0 > -1, wa < 0),
  retrospectively, but the sequestering constraint fails there. None of these numbers can be called
  held out.
* **Possible future held-out probes.** Later DESI releases, Euclid, Rubin/LSST and Roman. They are
  held out only if the collapse sector, cutoff and predicted (w(z), Omega_k) are frozen before
  unblinding. H0 is not a target: Omega_m, Omega_r and Omega_k are inputs.

## 7. Closure versus the GRUT closed-whole premise [A]

* **Formal analogy.** In GRUT every subsystem's bath is the rest of the closed whole. In global
  sequestering every local Einstein equation contains <tau>, an average over the entire 4D history,
  future included. Lambda is fixed "by the whole universe as the detector" (KP).
* **Difference 1: the kind of closure.** GRUT's closure (FORMAL_CORE DEF 1.1/1.7) is a quantum
  datum on a supplied algebra with a time-independent H. It contains no metric and does not require
  finite 4-volume. KP's closure is classical finiteness of the spacetime (compact space, bang,
  crunch) plus an a-posteriori future boundary datum. Neither implies the other.
* **Difference 2: where the cancellation comes from.** It comes from the global variational variables
  and the λ-scaling of the protected sector, not from finiteness (§4). The local and OS versions are
  just as "closed" in the sense of global constraints, yet they leave the residual free.
* **Consequence for GRUT.** If GRUT invokes "closed whole ⇒ global, non-additive variational
  principle", it adopts KP's assumption and must credit them. It then inherits the crunch
  requirement, and the simplest realization of that requirement fails (§5). I found no
  closure-derived assumption that the literature lacks which would fix ΔΛ (local/OS) or the
  collapse sector (global). Status: **no GRUT-native principle identified**.

## 8. Four-question gate

| Version | Q1 extra structure | Q2 fixes residual? | Q3 stable with heavy threshold/loops? | Q4 independent consequence |
|---|---|---|---|---|
| Global (KP 2014) | Global Lambda, lambda; sigma outside the integral; protected sector and regulator on lambda^2 g | Not free, but a functional of the supplied history; sign set by w near turnaround/crunch; ~O(1) M^2 H_age^2 times a crunch log | Matter loops: yes, every order (heavy-threshold V_vac cancels identically; lambda shifts by O(δV/mu^4)). Graviton loops: no (needs M_c ≲ 0.1 GeV [E]) | Crunch, closed compact space, transient w; no number without a collapse sector |
| Global + linear trigger (KP 2015) | Adds V = m^3 phi with m^3 ≈ M H0^2, chosen | Lambda_res ≈ m^3 M [0.10 ln(rho_cut/m^3 M) + 0.6] [C]: the scale is input through m^3 | As above; the slope is technically natural | Predicts <R> = 0. Violated for s ∈ [0.3, 3] and Omega_k ∈ [-0.95, 0] [C] (open k < 0 makes <R> more negative [A]); also Ferreira–Avelino. Retrospectively fails |
| Local (KPSZ 2016) | Two 4-forms with nonlinear sigma, sigma^; local fields Lambda, kappa^2 | No: ΔΛ is a free flux-ratio constant | Matter loops yes; graviton loops no | None beyond GR + Lambda |
| Omnia sequestra (KP 2017) | theta R_GB plus 4-forms | No: ΔΛ is flux-set and free, with two roots | Matter and graviton vacuum loops yes; Lambda–volume dressing harmless; historic integrals UV-safe only with ≳ 92 more e-folds | Only a naturalness preference for a long future |
| Lombriser averaging (2018/19) | Halo-patch averaging plus a location prior | Gives 0.697 / 0.704 | Inherits global matter-loop sequestering | Retrospective match only |

## 9. Verdict

Sequestering removes the large protected-sector contribution in a radiatively stable way. That is
Q1 and Q3, with OS needed for graviton loops. It does not supply a principle that fixes the total
renormalized rho_vac/M^4 with sign and scale (Q2). The local and OS versions keep an arbitrary
constant. The global version replaces the constant with a functional of a history whose decisive
inputs are supplied: the trigger slope, the cutoff log, and the existence and timing of a crunch. Its
only sharp relation fails in the published realization. Nothing here is new; everything is credited
to the cited authors. The track does not recommend building a GRUT candidate on sequestering unless a
specified, frozen collapse sector and cutoff are written down first, together with the (w(z), Omega_k)
relation they imply.

## 10. Limits

* Homogeneous FRW only, with no perturbations or inhomogeneities. Black-hole interiors were
  estimated as negligible by KP [Q] and were not recomputed.
* One collapse mechanism (the linear trigger). Other triggers, such as a negative constant plus an
  axion or a steeper potential, were not scanned. The no-root result is specific to this family.
* The crunch cutoff is a sharp truncation at rho_cut. Planck-scale physics near the singularity is
  unknown, and Lambda_res depends logarithmically on it.
* Omega_m was held fixed while Omega_k varied. CPL fits are a local diagnostic, not a likelihood
  analysis. Phase-transition ΔV values are order of magnitude. Graviton-loop thresholds use O(1)
  Wilson coefficients.
* Only the global-limit equations of the local and OS versions were used. Their full canonical
  structure and boundary terms were not re-derived.
* Several follow-ups were read at abstract level only (table in §1).
