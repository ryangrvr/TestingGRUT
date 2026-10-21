# Track `unimodular`: unimodular gravity, its quantum versions and everpresent Lambda against the four-question gate

21 October 2026. Brief: `../../BRIEF.md`. Codex's question: issue #13 comment 6324148604; packet
`git show 02b4515:research/reality_framework/grut_scale_feasibility_20261021/NEXT_DIRECTION.md`.

Labels: **[D]** derived here; **[C]** computed by `code/unimodular_checks.py` (record:
`results/unimodular_checks.json`, 39/39 checks pass, seed 20261021); **[Q]** quoted or paraphrased from a
source at the read level in §1; **[E]** order-of-magnitude estimate; **[A]** argued, not a theorem.

Conventions: signature (-,+,+,+), hbar = c = 1. Reduced Planck mass M = (8 pi G)^(-1/2) = 2.435e18 GeV and
l_r = 1/M. Non-reduced M_P = G^(-1/2) = sqrt(8 pi) M and l_P = 1/M_P. **Target quantity: the total
renormalized rho_vac / M^4 (reduced convention).** Everything built from H0, Omega_m, Omega_Lambda, BBN,
CMB or SN data is **retrospective** (the data are already exposed).

Nothing here is a new theory, an empirical validation or a confirmatory test. No candidate is built in
this track. The other tracks of this packet (`obstruction`, `sequestering`) are cited where they already
cover a point; their results are not repeated as new.

---------------------------------------------------------------------------------------------------

## 0. Result

1. **Classical unimodular / trace-free gravity removes the large term and leaves a free constant [D].**
   The trace-free equations do not see T -> T - C g. Bianchi plus conservation then give
   M^2 R + T = 4 Lambda_hat, a constant of integration. The gravitating constant
   rho_eff = rho_vac + Lambda_hat = (M^2 R + tau^mu_mu)/4 is fixed by slice data. In FRW it is
   3M^2(H^2 + k/a^2) - rho_matter. Every value is a solution. Codex Q2 fails: the residual's sign and scale
   are not fixed. This is credited to Einstein (1919, via Weinberg), Weinberg RMP 1989 §7, Ellis et al. 2011
   and Padilla-Saltas 2015.

2. **Quantum unimodular gravity (Henneaux-Teitelboim, Unruh) moves Codex's obstruction one level up [D].**
   In the HT minisuperspace, the compensated shift (rho_vac, Lambda) -> (rho_vac + c, Lambda - c) changes the
   Lagrangian only by -c dT/dt. Here T is unimodular time, which counts 4-volume. The unimodular-time
   Schroedinger generator becomes h -> h - cI, so psi -> e^{icT} psi.
   * This is Codex's fixed-background identity with the external t replaced by unimodular time T.
   * It yields a **Lambda-translation lemma**: the kinematic closure structure is invariant under
     Lambda -> Lambda + c. That structure is L^2(dLambda), the HT inner product, unitarity in T,
     normalizability, uncertainty relations and entropies of |A(Lambda)|^2.
   * Such a principle can therefore fix only translation-invariant features, such as a spread. It cannot
     fix <Lambda>.
   * To fix the mean, a principle needs a Lambda-dependent weight from the gravitational dynamics: a
     Euclidean action, a boundary condition, an on-shell action, or post-selection.
   * Every published weight of that kind either drives Lambda -> 0+ (Baum-Hawking, Ng-van Dam 1990) or is
     ill-defined or supplied (§3.4).
   * The lemma restates Codex's identity in HT variables. Its structure is credited to HT 1989, Unruh 1989,
     Smolin 2009 eq. (24) and Das et al. 2023 footnote 2. It is not new.

3. **Everpresent Lambda is the one proposal here that fixes a scale [Q/D/C], and it fixes only the
   spread.**
   * **Principle.** Causal-set discreteness holds N fixed, so Poisson fluctuations give
     Delta V ~ sqrt(V) l^2. Conjugacy (Delta rho)(Delta V) ~ hbar then gives
     sigma(rho_Lambda) = alpha M^2 / sqrt(V).
   * **Inputs.** V is the 4-volume of the past light cone. <Lambda> = 0 is assumed by all authors. alpha
     is free.
   * **Magnitude, Planck-2018 LCDM [C, retrospective].**
     * V_cone(today) = 0.1229 H0^-4 = 1.01e240 M^-4, so (V M^4)^(-1/2) = 9.9e-121.
     * Observed rho_Lambda / M^4 = 7.15e-121.
     * A 1-sigma match therefore needs alpha_red = 0.72. In the G = 1 convention of Das et al. this is
       alpha = 0.029.
     * In the saturated-Gaussian normalization it corresponds to a discreteness length l = 0.83 l_r =
       4.2 l_P.
     * Other choices of "V" move this by factors of 2 to 23 (§4.2). The agreement is order of magnitude
       only, in the coincidence form rho ~ M^2 H^2.
   * **Statistics per epoch [D/C].**
     * sigma(Omega_Lambda) = 2.73 alpha in the radiation era, 1.81 alpha in the matter era, 2.08 alpha at
       recombination and 0.95 alpha today.
     * The lag correlation is sqrt(V1/V2). In the matter era it is about a^-3 per e-fold. The correlation
       between z = 1 and today is 0.19.
     * The sign is +/- with equal prior probability. Survivor bias gives 51-54% positive.
   * **Monte Carlo of the published Model 1 [C].** The run is converged in step size and started before
     BBN.
     * Realizations terminate at H^2 <= 0 before today in 35% of cases at alpha_red = 0.15 and 95% at
       0.20. At alpha_red >= 0.3, all 20,000 terminate.
     * P(Omega_Lambda0 in [0.64, 0.74]) <= 1e-4 at every alpha scanned.
     * The joint requirement (today window, |Omega(rec)| < 0.05 and |Omega(BBN)| < 0.07) is met by
       0 of 180,000 realizations.
     * An optimistic linear-Gaussian estimate peaks at P = 5e-5.
   * **Gate.**
     * Q1: yes, discreteness plus conjugacy.
     * Q2: partial. The scale is fixed; the sign is random; the mean is assumed.
     * Q3: the spread is insensitive to vacuum loops. The mean is a value statement.
     * Q4: the statistical consequences are well defined: fluctuating w, early dark energy and sign
       changes. Every quantitative test so far used exposed data. The simplest models are strained by the
       CMB (Das et al. 2024: Model 1 chi^2 more than 150 worse than LCDM). With independent horizon
       patches, Barrow's bound applies (my estimate: alpha_red <~ 2e-5, 4 x 10^4 below what is needed).

4. **The GRUT closure argument is suggestive, not sound, and not GRUT-native [A].** The argument is: no
   external clock, so time is internal, so unimodular time is conjugate to Lambda.
   * Step 1 is standard for constrained systems. GRUT's own DEF 1.1 and 1.6, however, supply an external t.
   * Step 2, singling out 4-volume as the clock, is an added postulate in the "clock slot" of FORMAL_CORE §7.
     It also needs the HT fields. Other clocks give inequivalent quantum theories, and Kuchar 1991 criticized
     the step.
   * Step 3, "therefore Lambda is fixed", does not follow, by the lemma of point 2.
   * Unruh 1989, Unruh-Wald 1989, Henneaux-Teitelboim 1989 and Sorkin 1987/1994 made the same argument
     (§5).

5. **Gate summary.** No proposal in this track fixes the total renormalized rho_vac / M^4 with its sign and
   scale and also survives quantum corrections without retuning. Everpresent Lambda comes closest: it fixes
   the fluctuation scale. Its mean is assumed, its amplitude parameter is fitted to exposed data, and its
   simplest realizations make the observed history improbable at the 1e-4 to 1e-5 level.

---------------------------------------------------------------------------------------------------

## 1. Sources and read levels

PDFs and text extracts are in the scratchpad only (`.../scratchpad/um_dl/`; the Weinberg translation is
under `.../scratchpad/lit/vac/`). No paper text is copied into the repository.

| Ref | Identifier | Read level |
|---|---|---|
| Weinberg, "The cosmological constant problem" | Rev. Mod. Phys. 61, 1 (1989), DOI 10.1103/RevModPhys.61.1; authorized Russian translation UFN 158, 639 (1989) | APS returned 403. §7 (unimodular, Teitelboim particle model, ADM path-integral remark, Witten 3-form) read **in the translation**; equations not extractable; statements paraphrased |
| Henneaux & Teitelboim, "The cosmological constant and general covariance" | PLB 222, 195 (1989), DOI 10.1016/0370-2693(89)91251-3 | **abstract** (INSPIRE); action as reproduced by Smolin 2009 §2.2, Padilla-Saltas eq. (3), Afshordi-Magueijo eqs. (1)-(3) |
| Unruh, "A unimodular theory of canonical quantum gravity" | PRD 40, 1048 (1989), DOI 10.1103/PhysRevD.40.1048 | **abstract** (INSPIRE); APS 403 |
| Unruh & Wald, "Time and the interpretation of canonical quantum gravity" | PRD 40, 2598 (1989), DOI 10.1103/PhysRevD.40.2598 | abstract |
| Kuchar, "Does an unspecified cosmological constant solve the problem of time in quantum gravity?" | PRD 43, 3332 (1991), DOI 10.1103/PhysRevD.43.3332 | abstract |
| Smolin, "The quantization of unimodular gravity and the cosmological constant problems" | PRD 80, 084003 (2009), arXiv:0904.4841 | §§1, 2, 4.3, 5, 6 **full**; §§3-4.2 not read |
| Ellis, van Elst, Murugan, Uzan, trace-free Einstein equations | CQG 28, 225007 (2011), arXiv:1008.1196 | abstract, §§1-4 (parts), §6 |
| Padilla & Saltas, "A note on classical and quantum unimodular gravity" | EPJC 75, 561 (2015), arXiv:1409.3573 | **full** |
| Carballo-Rubio, Garay, García-Moreno, status report | CQG 39, 243001 (2022), arXiv:2207.08499 | not re-read here; see `obstruction` track §4.4 |
| Sorkin, "Forks in the road, on the way to quantum gravity" | IJTP 36, 2759 (1997), gr-qc/9706002 (talk 1993) | §VI and §VII(b) (the Lambda "prediction") |
| Sorkin, "On the role of time in the sum-over-histories framework for gravity" | IJTP 33, 523 (1994), DOI 10.1007/BF00670514 | abstract |
| Sorkin, "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?" | AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675 | **full** |
| Ahmed, Dodelson, Greene, Sorkin, "Everpresent Lambda" | PRD 69, 103523 (2004), astro-ph/0209274 | **full** |
| Ahmed & Sorkin, "Everpresent Lambda II: structural stability" | PRD 87, 063515 (2013), arXiv:1210.2589 | abstract, introduction, results on Omega distribution (p. 10), conclusions |
| Barrow, "A strong constraint on ever-present Lambda" | PRD 75, 067301 (2007), **gr-qc/0612128** | **full** |
| Zwane, Afshordi, Sorkin, "Cosmological tests of Everpresent Lambda" | CQG 35, 194002 (2018), arXiv:1703.06265 | §§1-3, 4, start of 5, 8; §§6-7 not read |
| Das, Nasiri, Yazdi, "Aspects of Everpresent Lambda I" | JCAP 10 (2023) 047, arXiv:2304.03819 | abstract, §§2.2-4.3, 5.1, 6 |
| Das, Nasiri, Yazdi, "Aspects of Everpresent Lambda II: cosmological tests" | JCAP 10 (2024) 076, arXiv:2307.13743 | abstract, §§1, 3, 7 |
| Moradi, Yazdi, Zilhão, "Fluctuations and correlations in causal set theory" | CQG 42, 045017 (2025), arXiv:2407.03395 | abstract |
| Nasiri, "Stochastic dark sector" | PhD thesis, Imperial College London (2025), INSPIRE record (no arXiv) | abstract |
| Afshordi & Magueijo, "Lower bound on the cosmological constant from the classicality of the early Universe" | PRD 106, 123518 (2022), arXiv:2209.07914 | abstract, §§I-III (parts) |
| Ng & van Dam, "Possible solution to the cosmological constant problem" | PRL 65, 1972 (1990), DOI 10.1103/PhysRevLett.65.1972 | abstract |
| Ng & van Dam, "A small but nonzero cosmological constant" | IJMPD 10, 49 (2001), hep-th/9911102 | abstract |
| Daughton, Louko, Sorkin, "Initial conditions and unitarity in unimodular quantum cosmology" | gr-qc/9305016 (1993) | abstract |
| Josset, Perez, Sudarsky | PRL 118, 021102 (2017), arXiv:1604.04183 | abstract |
| Perez & Sudarsky, "Dark energy from quantum gravity discreteness" | PRL 122, 221302 (2019), arXiv:1711.05183 | abstract |
| Kaplan & Rajendran, "Redshifting the cosmological constant in unimodular gravity via nonlinear QM" | arXiv:2511.12897 (2025) | abstract |
| Gielen, "Boundary conditions and Hilbert spaces in no-roll quantum cosmology" | CQG 43, 195001 (2026), arXiv:2606.17163 | abstract |
| Gielen & Neves, unimodular de Sitter unitarity | CQG 42, 145001 (2025), arXiv:2412.01907 | abstract |
| Davies & Magueijo, "Teleocosmology and quantum post-selection" | arXiv:2606.02514 (2026) | abstract |
| Gielen & Menéndez-Pidal, "Singularity resolution depends on the clock" | arXiv:2005.05357 | title and listing only |

**Identifier corrections.**
* The task text gives Barrow's paper as "astro-ph/0612128". That identifier is an unrelated paper, *Ionisation-induced star formation I*; I checked the downloaded PDF. The correct identifier is **gr-qc/0612128**.
* The Sorkin 2007 title is "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?".

**Follow-up search for "inhomogeneous everpresent Lambda".** I searched the arXiv API ("everpresent",
"ever-present Lambda", "causal set" plus "cosmological constant") and INSPIRE (titles and the papers
citing ADGS since 2022). I found no published covariant inhomogeneous everpresent-Lambda model up to
October 2026. What exists:
* Das-Nasiri-Yazdi I and II (2023-24): homogeneous Models 1 and 2.
* Moradi-Yazdi-Zilhão (2024-25): tools for correlations of causal-set quantities, presented as paving the
  way for aspects of everpresent Lambda (abstract).
* Nasiri's 2025 thesis (abstract): initial conditions and inflation-like phases, SN and CMB fits.

ZAS 2018 treat inhomogeneities by a CAMB quintessence-fluid stand-in [Q]. The search used the
listed queries only, so coverage is limited.

---------------------------------------------------------------------------------------------------

## 2. Classical unimodular gravity

### 2.1 The vacuum term drops out; Lambda is an integration constant [D; C: A1-A2]

Vary the Einstein-Hilbert plus matter action only at fixed sqrt(-g). Equivalently, add a Lagrange
multiplier, -lambda(x)(sqrt(-g) - eps0) (Padilla-Saltas eq. 1). Only the trace-free part survives:

    M^2 (R_mu nu - (1/4) R g_mu nu) = T_mu nu - (1/4) T g_mu nu .                       (2.1)

* **Blindness.** For any constant C, TF(T - C g) = TF(T) identically. This is check A1, done
  symbolically in C with a random rational Lorentzian metric.
* **Integration constant.** Take the divergence of (2.1) and use nabla^mu G_mu nu = 0, so
  nabla^mu R_mu nu = (1/2) d_nu R, and nabla^mu T_mu nu = 0, which is assumed separately. This gives
  (1/4) M^2 d_nu R = -(1/4) d_nu T. Hence

      M^2 R + T = 4 Lambda_hat  (constant),      M^2 G_mu nu = T_mu nu - Lambda_hat g_mu nu .   (2.2)

* **Effective constant.** Write T = tau - rho_vac g, where tau is everything that is not vacuum. Then

      rho_eff = Lambda_hat + rho_vac = (M^2 R + tau^mu_mu)/4 ,                          (2.3)

  which does not depend on rho_vac at any order (checks A2a-A2d in FRW).
* **FRW form.** The trace-free equation is a''/a - H^2 - k/a^2 = -(rho + p)/(2 M^2). This is Ellis et al.
  eq. (22) with 4 pi G = 1/(2M^2), and rho_vac cancels from rho + p. Conservation makes
  3M^2(H^2 + k/a^2) - rho_tot constant (A2b). The gravitating constant is therefore
  rho_eff = 3M^2(H^2 + k/a^2) - rho_matter on any slice (A2c).

This reproduces Weinberg §7 (eqs. 7.3-7.7 in the translation, paraphrased). There the cosmological
constant appears only as a constant of integration, vacuum fluctuations cancel from the right side of
(7.3), and flat space is a solution. Weinberg states that the remaining problem is why one should pick
the flat solution from among all possible ones [Q, translation]. Ellis et al. (§6) make the same point:
the TFE "does not solve the issue of why the cosmological constant has the value it has today" [Q].

### 2.2 Numerical illustration [C: A3]

I integrated the trace-free second-order ODE for a(t) (dust plus radiation, k = 0) twice: once with
rho_vac = 0 and once with rho_vac = 1e6 in code units. The two histories agree to 0 (machine
identical), because rho_vac never enters. The gravitating constant 3H^2 - rho_matter stays equal to its
initial-slice value (2.1) to 2e-8. **Any value is a solution, and its value is read off the initial data.**
Reading it from (H0, rho_m0) is a retrospective measurement, not a prediction. Codex's note applies: H0 is
not a parameter-free target.

### 2.3 Gate answers for classical unimodular gravity

* **Q1 [D].** The extra structure is a fixed volume form, or a restricted variation, which leaves the
  theory invariant only under transverse diffeomorphisms (TDiff). It removes the trace (constant) mode of
  T from the local equations.
* **Q2 [D/Q].** No. Lambda_hat is an arbitrary integration constant of either sign and any scale.
* **Q3 [D/Q/A].** The local equations are blind to every constant contribution at every loop order (2.1).
  Padilla-Saltas argue that once the restriction is implemented with a multiplier, matter loops shift the
  multiplier's constant boundary value. "The cosmological constant is precisely this (constant) boundary
  value", so it must be retuned exactly as in GR [Q, full read]. The `obstruction` track §4.4 already
  scored this dispute (the two sides hold different quantities fixed). I do not repeat it. In both
  readings there is no predicted value whose stability could be tested.
* **Q4 [Q].** None classically. Padilla-Saltas: classical unimodular gravity is locally gauge-fixed GR
  with identical predictions. Ellis et al. check that stellar junction conditions are unchanged.

---------------------------------------------------------------------------------------------------

## 3. Quantum unimodular gravity: Lambda conjugate to 4-volume

### 3.1 Henneaux-Teitelboim form and the boundary-phase identity [Q action; D/C: B1]

The HT action, as reproduced by Smolin §2.2, Padilla-Saltas eq. (3) and Afshordi-Magueijo eq. (1), is

    S_HT = ∫ [ sqrt(-g) (M^2 R/2 - Lambda) + Lambda d_mu calT^mu ] + S_m .

It has no background structure. Varying calT gives d_mu Lambda = 0. Varying Lambda gives
d_mu calT^mu = sqrt(-g), so T(Sigma) = ∫_Sigma calT^0 counts the 4-volume to the past of Sigma. HT's
abstract describes this global zero mode as a "cosmic time" canonically conjugate to the cosmological
constant [Q, abstract].

FRW minisuperspace (unit comoving volume, code B1):

    L = -3M^2 a a'^2/N + 3M^2 k N a - N a^3 (rho_m(a) + rho_vac + Lambda) + Lambda T' .

The Euler-Lagrange equations are Lambda' = 0, T' = N a^3, and the Friedmann constraint with
rho_vac + Lambda.

* **The vacuum shift is a boundary term [D, B1b].**
  L(rho_vac + c, Lambda - c) - L(rho_vac, Lambda) = -c T' exactly. The action changes by
  -c (T_f - T_i) = -c x (4-volume). The constraint depends only on Lambda_eff = Lambda + rho_vac (B1c).
* **Canonical form.** Take p_T = Lambda. The constraint is a^3 p_T + K(a, p_a) + a^3 rho_vac = 0, where K
  collects the gravitational kinetic and curvature terms and the matter terms. In the quantum theory this
  is a Schroedinger equation in T, with generator h = -K/a^3 - rho_vac. Unruh's abstract: the quantum
  version "has a normal 'Schrödinger' form of time development", and the small present Lambda "becomes a
  genuine question of initial conditions" [Q, abstract].
* **The shift in the quantum theory [D, C: B2].** Under rho_vac -> rho_vac + c the generator becomes
  h_c = h - cI, so psi_c(T) = e^{icT} psi(T) (the sign is convention). Probabilities in unimodular time
  agree to 5e-16, and the relative phase is e^{icT} to 1e-15 (12-dimensional stand-in for h).
  **This is Codex's identity H -> H + cI with the external t replaced by unimodular time T.** The
  `obstruction` track §4.3 made this bridge at the level of e^{-i rho_vac V4}. Here it is shown explicitly
  in the HT canonical form.

### 3.2 The Lambda-translation lemma [D; C: B3]

**Setting.** A wavefunction of the universe in HT form,
psi(q, T) = ∫ dLambda A(Lambda) e^{-i Lambda T} psi_Lambda(q) (Afshordi-Magueijo eq. 4; Smolin eq. 70). It
has the natural inner product <psi1|psi2> = ∫ dLambda A1* A2 (A-M eq. 10) and [Lambda, T] = i (in units of
hbar_eff).

**Lemma.** Let P be any selection rule built only from this kinematic structure: normalizability,
unitarity in T, the inner product, uncertainty relations, entropies or purity of |A|^2, and Gaussianity.
P is invariant under the Weyl translation W_c = e^{icT}: A(Lambda) -> A(Lambda - c), with Lambda -> Lambda + c
and T unchanged. If P selects a state, it selects the whole translation orbit. So P cannot fix <Lambda>,
and it cannot fix <Lambda_eff> = <Lambda> + rho_vac.

**Proof.** The Lebesgue measure dLambda is translation invariant, [Lambda + c, T] = [Lambda, T], and
|psi(T)|^2 changes only by the phase e^{-icT}. **Numerical check (B3).** Under a shift by c = 7.5:
* the norm, sigma_Lambda, sigma_T, the entropy of |A|^2 and the profile |psi(T)|^2 are unchanged to
  <= 1e-10;
* the mean moves by exactly 7.5;
* the unchirped Gaussian saturates sigma_Lambda sigma_T = 1/2 to 1e-13.

**Consequence [A].**
* A closure-type principle in quantum unimodular gravity can at most fix the **spread** of Lambda. The
  causal-set Poisson argument (§4) does exactly this, and only this.
* Fixing the mean requires weights that depend on Lambda non-covariantly under translation. Such weights
  come from the gravitational dynamics psi_Lambda(q), not from closure. Examples:
  * normalizability in a only for Lambda < 0 (obstruction track C7-C8);
  * the Euclidean action;
  * a boundary condition at a = 0;
  * an on-shell action;
  * a final state.

This is the precise form of "unimodular gravity relocates the free constant to state selection"
(`obstruction` §4.3). It is a restatement, not a new result. The shift symmetry T -> T + g C,
phi -> phi - 4 pi G C is Smolin eq. (24). Das et al. I, footnote 2, state that the mean of Lambda is
arbitrary in unimodular gravity, that vacuum energy can be absorbed into the mean of the Fourier variable,
and that "a more complete theory should advise the correct mean" [Q].

### 3.3 What could fix Lambda in the quantum theory? Published state principles [Q/D/C]

| Principle | Extra structure | What it fixes | Read level / status |
|---|---|---|---|
| Lambda eigenstate | none | any value; not normalizable (A-M §II) | [Q] |
| normalizable wave packet A(Lambda) | HT inner product | spread with sigma_Lambda sigma_T >= hbar_eff/2; mean free (lemma) | [D] |
| A-M classicality | requires the unimodular clock to be classical in a patch of comoving volume V_c | relates sigma(Lambda) <= Lambda to the onset of classicality (~7e11 GeV with the observed Lambda as input) | [Q abstract; uses the observed Lambda] |
| no-boundary / tunnelling in unimodular minisuperspace | boundary condition at a = 0, self-adjoint extension | DLS 1993: no-boundary gives non-unitary flux in T. Gielen 2026: with Lambda left free, one-parameter family of boundary conditions; energy eigenstates (Lambda values) remain a state choice | [Q abstracts] |
| Baum-Hawking / Ng-van Dam 1990: integrate Lambda in the Euclidean path integral | Euclidean weight e^{-S_E}, S_E(dS) = -24 pi^2 M^4 / rho_Lambda | Lambda = 0+ "overwhelmingly" (Ng-van Dam abstract). Exponent at the observed value 24 pi^2 / 7.15e-121 = 3.3e122 [C]; weight maximal at 0+, not normalizable. Sign fixed (+), scale zero | [Q abstract; C] |
| Ng-van Dam / Smolin 2009 stationary phase | integrate Lambda; stationarity in total 4-volume; mediocrity | Smolin (84)-(86): Lambda/(4 pi G) ≈ <L_m - T/2>, then Lambda/G ≈ 2 pi rho, "remarkably close for such a speculative argument" [Q]. Ill-defined: see §3.4 | [Q full; D] |
| two-time post-selection (Davies-Magueijo 2026) | final Chern-Simons state, connection variables, unimodular time | acceleration with Lambda = 0 in the forward Hamiltonian; extrapolated w < -1 (phantom) | [Q abstract] |
| nonlinear QM (Kaplan-Rajendran 2025) | state-dependent Hamiltonian terms (average Ricci scalar) | the "shadow" Lambda redshifts away. The abstract describes a late cosmology dominated by matter and radiation, so the observed acceleration would need another source (my reading) | [Q abstract] |
| energy diffusion (Josset-Perez-Sudarsky 2017; Perez-Sudarsky 2019) | violation of nabla T = 0 from Planckian granularity, absorbed by the unimodular Lambda(t) | an observed-order Lambda "from simple dimensional analysis" with a phenomenological parameter | [Q abstracts] |
| causal-set Poisson (Sorkin; §4) | fixed number N of elements; Poisson sprinkling | spread only; mean assumed 0 | [Q full] |

No entry fixes the observed residual's sign and scale from first principles without either the observed
value (A-M), an assumed mean (Sorkin), a supplied final state (D-M), or a contested Euclidean measure
that gives zero (Ng-van Dam).

### 3.4 Why the Ng-van Dam / Smolin relation is not well defined as stated [D; C: C1-C2]

Use S = ∫ sqrt(-g)[M^2 R/2 - rho_Lambda + L_m]. On shell the trace gives M^2 R = 4 rho_Lambda - T, so the
on-shell density is rho_Lambda + L_m - T/2 (C1). Setting the on-shell action per unit 4-volume to zero
(Smolin's (84), up to his normalization) gives rho_Lambda = -<L_m - T/2>. The on-shell value of a
perfect-fluid Lagrangian differs between actions that give the same field equations (L_m = p versus
L_m = -rho):

| fluid | L_m = p | L_m = -rho |
|---|---|---|
| dust | rho_Lambda = -<rho>/2 | +<rho>/2 |
| radiation | -<rho>/3 | +<rho> |

**Both the sign and the magnitude depend on a choice that leaves the dynamics unchanged [D, C2].** Two
further problems:
* Stationarity in Lambda at fixed history gives dS/d rho_Lambda = -V4, which is never zero at finite V4.
  A flat-measure Lambda integral instead imposes delta(V4 - V_fixed), the Fourier conjugacy (C3).
* The relation needs a finite total 4-volume and a measurement theory. Smolin himself calls it "so
  formal" and says it "requires a precise measurement theory for quantum cosmology, which we lack" [Q].

Structurally this is a history average of the matter trace or Lagrangian, the same output class as
global sequestering. Credit goes to Kaloper-Padilla, scored in the `sequestering` track, and the finding
is not repeated here.

### 3.5 Equivalence and the problem of time [Q]

Padilla-Saltas argue that quantum unimodular gravity equals quantum GR if the HT fields are auxiliary,
that is, never on external legs. Any inequivalence is a choice made in the UV. Integrating over the free
boundary value lambda_0 leaves a global integral ∫dlambda_0 f(lambda_0) whose measure f is unspecified
(their eq. 12). That is the state/measure slot again.

Kuchar 1991 (abstract) holds that cosmological time labels only equivalence classes of hypersurfaces
separated by zero 4-volume. Unless complemented by a hypertime, it does not uniquely set the conditions
for measuring geometric variables, classically or quantum mechanically. Padilla-Saltas call the
problem-of-time claim "strongly refuted" by Kuchar [Q].

---------------------------------------------------------------------------------------------------

## 4. Everpresent Lambda

### 4.1 The principle as equations, with every input listed [Q, D]

Sources: Sorkin 1997 §VII(b); Sorkin 2007; ADGS 2004 §II; Das et al. I §§2.2-3.1.

1. **Number fixed.** The sum over causal sets is taken at fixed number N of elements, because "n plays the
   role of time" in sequential growth (Sorkin 2007; ADGS) [Q]. In the continuum this is the unimodular
   path integral at fixed V: Z(V) = ∫ Dg delta(V - ∫sqrt(-g)) e^{iS} (Das I eq. 3.1).
2. **Conjugacy.** Z(lambda) = ∫ dV e^{-i lambda V} Z(V) is the GR path integral with Lambda = lambda, so
   (Delta rho_Lambda)(Delta V) >~ hbar/2 (Das I eqs. 3.2-3.5). Das et al. state that no canonical proof of
   [Lambda, V] exists (their footnote 1) [Q].
3. **Poisson kinematics.** Sprinkling at density 1/l^4 gives V = N l^4 +/- sqrt(N) l^4 [Q].
4. **Spread.** With saturation, sigma(rho_Lambda) = (1/2) / (l^2 sqrt(V)). In reduced units this is
   sigma(rho_Lambda)/M^4 = alpha_red (V M^4)^(-1/2), with alpha_red = (1/2)(l_r/l)^2 [D, B5].
5. **Mean.** <Lambda> = 0 is assumed: "for reasons still to be discovered" (Sorkin 2007), "We must assume
   ... that this 'target value' is zero" (ADGS), and Das I footnote 2 [Q]. This is the slot the lemma (§3.2)
   says closure cannot fill.
6. **Which V.** V is the 4-volume of the past light cone of a point,
   V(tau) = (4 pi/3) ∫_0^tau dtau' a(tau')^4 (tau - tau')^3 (ADGS eq. 3; Das I eq. 4.4). It is a choice
   "somewhat at odds with the meaning that V has in the unimodular context" (ADGS) [Q].
7. **Dynamics (Model 1).** S_{n+1} = S_n + alpha xi_n sqrt(Delta V_n), with xi_n ~ N(0,1), and
   rho_Lambda = S/V. The scale factor follows **only** the Friedmann equation; the second Friedmann
   equation is dropped, because a varying Lambda is incompatible with the Bianchi identity (ADGS §II;
   Ahmed-Sorkin 2013 test other linear combinations) [Q].
8. **Model 2 (ZAS).** Omega_de = tanh(Omega_hat), where Omega_hat is Gaussian in ln a with variance alpha^2
   and correlation exp(-mu (Delta ln a)^2 / 2). It has two parameters (alpha, mu) [Q].
9. **Background.** Spatial homogeneity, k = 0. Inhomogeneity is not modelled; ZAS use a CAMB quintessence
   fluid as a stand-in [Q].

**Every parameter and input.**
* alpha, or l, or Sorkin's l^2 ~ 3 kappa.
* The random seed and realization.
* S_0, V_0 and the start epoch.
* The step recipe, since termination statistics depend on it (§4.4; Das I §5.2 and Fig. 3).
* Omega_m h^2, Omega_r h^2 and flatness.
* For Model 2, mu.
* No principle fixes <Lambda> or alpha.

**alpha conventions [C: D3a].**
* ADGS and Das et al. report alpha ~ 0.01-0.02 and 0.005-0.012 respectively. Das et al. define
  Lambda = 8 pi alpha / sqrt(V) (eq. 3.13), so rho_Lambda / M_P^4 = alpha_G1 (V/l_P^4)^(-1/2).
* The reduced coefficient is alpha_red = 8 pi alpha_G1 (exact by construction; D3a checks it to 4e-15).
* Sorkin 2007's ansatz Lambda = (1/l^2)/sqrt(V) with l^2 ~ 3 kappa gives alpha_red ~ 1/3, that is
  alpha_G1 ~ 0.013.
* With the mapping alpha_red = 8 pi alpha_G1, the published ranges therefore all sit at alpha_red ~ 0.12-0.5.
* My coarse-grid run (§4.4) gives survival fractions of 0.17 and 0.05 at alpha_red = 0.25 and 0.30. Das I
  report 0.32 and 0.054 at alpha_G1 = 0.010 and 0.012. Same order, consistent with the mapping [C/E].
* The ADGS footnote states units with hbar = c = m_p = 1 and m_p the reduced Planck mass. I do not assert
  which convention their code used.

### 4.2 Magnitude today [C, retrospective]

Inputs: Planck 2018 base LCDM with h = 0.6736, Omega_m = 0.3153, T_cmb = 2.7255 K and N_eff = 3.046
(massless approximation), giving Omega_r = 9.2e-5 and Omega_Lambda = 0.6846. The age is t0 = 13.80 Gyr.

**Two checks of the volume integral.** First, nested moment integrals agree with direct quadrature to
<= 1.6e-5 at four epochs (D2). Second, the power-law eras give exact Beta-function values:
* radiation era: H^4 V = 4 pi / 840 = 0.014960 (D1a, D1c);
* matter era: H^4 V = 0.033849 (D1b).

| Quantity | Value |
|---|---|
| observed rho_Lambda / M^4 (reduced) | **7.15e-121** |
| observed rho_Lambda / M_P^4 (G = 1) | 1.13e-123 |
| Lambda l_P^2 (curvature, G = 1) | 2.85e-122 |
| (H0/M)^2 | 3.48e-121 |
| (t_P/t0)^2 (Barrow's "1e-120", t_P = sqrt(G)) | 1.53e-122 |
| H0^4 V_cone(today) | 0.1229 |
| V_cone M^4 | 1.01e240 |
| V_cone / l_P^4 (Sorkin's N if l = l_P) | 6.4e242 |
| (V M^4)^(-1/2) | **9.93e-121** |
| (V/l_P^4)^(-1/2) | 3.95e-122 |
| alpha_red* = (rho_Lambda/M^4)(V M^4)^(1/2) = 3 Omega_Lambda sqrt(H0^4 V) | **0.720** |
| alpha_G1* (Das et al. convention) | 0.0286 |
| discreteness length for 1 sigma (saturated Gaussian) | l = 0.83 l_r = 4.2 l_P |

So "Lambda ~ V^(-1/2) ~ 1e-120" is right to within a factor of order one in the reduced convention:
V^(-1/2) = 9.9e-121 against 7.2e-121. In the G = 1 convention the agreement is within a factor of about
30 (3.95e-122 against 1.1e-123), because the coefficient differs by 8 pi. The ubiquitous "1e-120" spans
1.5e-122 to 7e-121 across conventions.

**Dependence on the choice of V (D4).** Here "ratio" means rho_obs / (V M^4)^(-1/2), the required alpha_red:

| V | ratio |
|---|---|
| past light cone | 0.72 |
| H0^-4 | 2.05 |
| t0^4 | 1.85 |
| Hubble volume x age | 4.10 |
| particle-horizon volume x age | 23.2 |

The order-of-magnitude statement is robust. The O(1) coefficient is not determined by the principle.

**What was ex ante [Q/A].**
* Sorkin's statement that Lambda is "of this magnitude ... 1e-120 in natural units", with whatever
  mechanism drives Lambda to vanish leaving such a residue, was made before the 1998 supernova
  detections: reported from 1987 talks, preprint 1993, arXiv June 1997.
* Sorkin's own reference list cites Krauss-Turner, "The cosmological constant is back" (GRG 1995, year from
  recall), so observational hints of Lambda ~ rho_crit already existed.
* The statement has the coincidence form rho ~ M^2 H^2. It leaves the sign open (+/-) and leaves an O(1)
  coefficient free. Any Lambda ~ H^2 estimate shares this information content (`obstruction` row R7).
* It is a genuine prior order-of-magnitude statement. It is not a confirmatory test of the mechanism.
* alpha was chosen after 1998 (ADGS 2002) so that Omega_Lambda0 ~ 0.5-0.7 is reachable.

### 4.3 Per-epoch statistics [D; C: D2, E1]

Model 1 is a random walk in V. S(V) has Var(S) = alpha^2 V, and rho_Lambda = S/V. Hence:
* std(rho_Lambda) = alpha / sqrt(V) [D];
* Cov(rho(V1), rho(V2)) = alpha^2 / V2 for V1 < V2, so the correlation is sqrt(V1/V2) [D];
* X = S/sqrt(V) is an Ornstein-Uhlenbeck process in ln V with unit stationary variance and correlation
  exp(-Delta ln V / 2) [D]. Its zero set is fractal in the continuum, so "number of sign changes" depends
  on the time resolution. The meaningful statement is the correlation time [A].
* In units of the critical density, sigma(Omega_Lambda) = alpha_red / (3 sqrt(H^4 V)).

| epoch (Planck-2018 LCDM background) | H^4 V | sigma(Omega_Lambda)/alpha_red |
|---|---|---|
| BBN, a = 1e-9 (pure radiation: 4 pi/840) | 0.01496 | 2.725 |
| matter-radiation equality | 0.02063 | 2.321 |
| recombination, z = 1089.8 | 0.02566 | 2.081 |
| z = 2 | 0.0372 | 1.729 |
| z = 1 | 0.0452 | 1.567 |
| z = 0.5 | 0.0609 | 1.351 |
| today | 0.1229 | 0.951 |
| pure matter era (exact) | 0.03385 | 1.812 |

Lag correlations of rho_Lambda (E1):

| pair | correlation |
|---|---|
| z = 0.5 vs 0 | 0.40 |
| z = 1 vs 0 | 0.19 |
| z = 2 vs 0 | 0.060 |
| recombination vs 0 | 8.5e-10 |
| BBN vs 0 | 4e-33 |

Per e-fold of a the correlation is e^-3 = 0.05 in the matter era (t ∝ a^(3/2), V ∝ t^4) and e^-4 in the
radiation era. The model's coherence time is therefore about a third to a half of an e-fold in a, that is
about one Hubble time. The Monte Carlo reproduces sqrt(V1/V2) for z = 1 versus 0 to 0.015 (F2).

"Omega_Lambda ~ alpha-weighted O(1) stochastic" thus has a precise form. With one alpha,
**sigma(Omega)(BBN) / sigma(Omega)(today) = 2.9 is fixed by past-light-cone geometry [D/C].**
* If the observed Omega_Lambda0 = 0.685 is a 1-sigma event, alpha_red ≈ 0.72. Then sigma(Omega) is about
  2.0 at BBN and 1.5 at recombination.
* If instead alpha_red <= 0.0128 (|Omega_BBN| < 0.07 at 2 sigma), Omega_Lambda0 = 0.685 is a 56-sigma
  event in the linear approximation (E2).
* An optimistic linear-Gaussian joint probability peaks at P = 5.0e-5 at alpha_red = 0.42 (E2). The
  conditions are a today window [0.64, 0.74], |Omega_rec| < 0.05 and |Omega_BBN| < 0.07. The estimate
  treats the epochs as independent, which the 1e-9 correlations justify, and ignores termination.
* This is the quantitative form of Barrow's remark that the model "would require ... an effective alpha
  that increased with time" [Q], and of Das II's "suppressing the values of Lambda near the last
  scattering surface" [Q].

### 4.4 Sign and fluctuation statistics: Monte Carlo of Model 1 [C: F]

**Implementation.** Units are rho / (M^2 H_m0^2) and V H_m0^4. In these units the 1e-120 scales cancel
and the model is scale-free. Further choices:
* Omega_m h^2 = 0.1430 and Omega_r h^2 = 4.18e-5 (Planck 2018, retrospective). H_m0 = 37.8 km/s/Mpc is the
  rate today if rho_Lambda = 0.
* Log grid with 80 steps per e-fold.
* V from the nested Cauchy form V = 8 pi J3, with J0 = ∫ a^4 dtau and J_{n+1} = ∫ J_n dtau. This form has
  no cancellation; the moment form lost accuracy and was replaced. The volume error is O(dx^2): 0.58%,
  0.15% and 0.04% at 40, 80 and 160 steps (F0).
* Stationary start S_i ~ N(0, alpha^2 V_i) at a = 1e-10, before BBN.
* Termination when H^2 <= 0 at a grid point; Model 1 has no dynamics beyond that point. The termination
  fraction is converged to 0.7 points between 80 and 160 steps per e-fold (F4).
* 20,000 realizations per alpha.
* At alpha_red = 0.02 the Monte Carlo std of Omega_0 matches the linear prediction alpha/(3 sqrt(H_m0^4 V))
  to 0.3% (F1).

| alpha_red (alpha_G1) | terminated before today | Omega_Lambda0 16/50/84% (survivors) | P(Omega>0 \| survive) | P(today in [0.64,0.74]) | P(\|Omega_rec\|<0.05) | P(\|Omega_BBN\|<0.07) | joint | H0 median (km/s/Mpc) |
|---|---|---|---|---|---|---|---|---|
| 0.02 (0.0008) | 0 | -0.036 / 0.001 / 0.036 | 0.507 | 0 | 0.77 | 0.80 | 0 | 37.8 |
| 0.05 (0.002) | 0 | -0.091 / 0.003 / 0.089 | 0.515 | 0 | 0.37 | 0.40 | 0 | 37.9 |
| 0.10 (0.004) | 0.002 | -0.18 / 0.01 / 0.17 | 0.524 | 0 | 0.19 | 0.21 | 0 | 38.0 |
| 0.15 (0.006) | 0.35 | -0.28 / 0.03 / 0.25 | 0.543 | 0 | 0.086 | 0.13 | 0 | 38.4 |
| 0.20 (0.008) | 0.954 | -0.38 / 0.03 / 0.34 | 0.538 | 1e-4 (2 runs) | 0.005 | 0.069 | 0 | 38.5 |
| 0.25 (0.010) | 0.9999 | (2 survivors) | - | 0 | 0 | 0.023 | 0 | - |
| 0.30, 0.50, 0.72 | 1.000 | - | - | 0 | 0 | <= 0.006 | 0 | - |

* **Termination versus need.** alpha_red* = 0.72 is required for 1 sigma today. Every realization started
  before BBN terminates at alpha_red >= 0.3. The joint condition is met by 0 of 180,000 realizations, so
  the 95% upper limit is about 1.5e-4 per alpha.
* **Survivor bias.** The positive fraction today is 0.507 to 0.543, rising with alpha (F3). This agrees
  with Ahmed-Sorkin 2013, who report positive values at "fifty to sixty percent", and
  Omega_final > 0.7 about 2% of the time in their variant [Q].
* **H0 among survivors.** The median is about 38 km/s/Mpc, the Lambda-free value. This agrees with Das I,
  who find the H0 distribution peaked near 40 and the observed value 2.7-7.9 sigma away [Q].
* **Start-epoch dependence (5,000 realizations each, F_start_epoch_scan).**
  * Starting later suppresses early excursions. With a_start = 1e-3 (after recombination), survival is
    0.40 at alpha_red = 0.3 and P(today window) reaches 0.016. With a_start = 1e-4 it reaches 0.006, and
    with 1e-6 it reaches 0.0012.
  * The "good fits" of ZAS and Das II therefore rely either on starting late, or on selecting seeds with
    already-exposed data (ZAS Method 2 treats the seed as an MCMC parameter), or on suppressing Lambda near
    last scattering.
* **Grid dependence (F_grid_dependence).** Das I's equal-steps-in-a recipe (a_i = i/n, n = 1e4, from
  a = 1e-4) gives survival 0.35, 0.17 and 0.05 at alpha_red = 0.2, 0.25 and 0.3. The converged log grid
  from the same start gives 0.84, 0.40 and 0.10. Das I already report that termination depends on the step
  (§4.3, Fig. 3). Published termination statistics are recipe-dependent at the factor-2 level.

### 4.5 Radiative stability (Codex Q3) [D/Q/A]

* **Spread.** At fixed V the vacuum term contributes a global phase e^{-i Lambda_tilde V / 8 pi G}. It
  cancels from normalized correlators (Das I eq. 2.14) [Q]. In the Fourier picture it shifts the mean of
  lambda and leaves its spread untouched [D]. The spread depends only on l and V. A heavy threshold, of size
  m^4/(64 pi^2) = 4.3e30 rho_Lambda for one electron degree of freedom (H, as in `obstruction` §3.4),
  does not change it.
* **Mean [D, C: H1-H2].** Write Lambda_eff = Lambda_tilde + lambda_bar + delta_lambda.
  * Prescription (i), <lambda> = 0 for the multiplier: the mean is Lambda_tilde. It moves by the full
    threshold shift each time (H1). This is unstable.
  * Prescription (ii), <Lambda_eff> = 0 imposed on the total: the mean is stable (H2) only because the
    total is stated. That is a value statement (obstruction row R1), the "zero" target that ADGS assume.
* **alpha [A].** alpha = (1/2)(l_P/l)^2 involves the renormalized G. Matter loops shift M^2 by
  O(N m^2/16 pi^2), which is negligible for m << M. The model has no mechanism for l itself.
* **Graviton loops.** They are covered by the trace-free blindness for any constant [D]. Whether the HT or
  WTDiff quantum theory preserves this is disputed (Padilla-Saltas versus Carballo-Rubio et al.; see
  `obstruction` §4.4) [Q].

### 4.6 Testable consequences and observational status (Codex Q4)

Consequences, given (alpha, start, recipe) [D/Q]:
* (a) Stochastic dark energy with coherence time of about a Hubble time. corr(rho_Lambda(a1),
  rho_Lambda(a2)) = sqrt(V1/V2), which is about (a1/a2)^3 in the matter era. Generic w(z) far from -1 at
  z < 2, and sign changes (Das I §6: "likely to change sign roughly every Hubble time").
* (b) Early dark energy with |Omega_Lambda| ~ 2-2.7 alpha_red at recombination and BBN.
* (c) Termination or recollapse when H^2 -> 0. The dynamics there is unknown (ADGS §IV; Das I §4.3).
* (d) If Lambda in disjoint past cones is independent, horizon-scale potentials give
  dT/T ~ dOmega/4 at recombination. With the uniform-sphere Sachs-Wolfe estimate (G),
  alpha_red <~ 2e-5 for dT/T <~ 1e-5. That is 3.7e4 below alpha_red* [E]. Barrow: "alpha <~ 1e-5"
  makes Lambda "effectively never-present" [Q].
  * ZAS dispute the independence: everpresent Lambda is meant to be quantum and nonlocal [Q].
  * No covariant inhomogeneous model exists (§1). This consequence is undetermined, not refuted.
* (e) Tracking: |rho_Lambda| scales as a^-4 in the radiation era and a^-3 in the matter era (ADGS Fig. 2;
  reproduced, D1).
* (f) Everpresent behaviour only in 3+1 dimensions (Barrow; Sorkin 2007) [Q]. Not testable.

**Already-exposed data, which cannot be called held out.**
* SNe Ia: 1998 onward, JLA, Pantheon+SH0ES (Das II).
* CMB: COBE, WMAP, Planck 2015 (ZAS) and Planck 2018 (Das II).
* BAO: 6dF, MGS, BOSS DR11 including Ly-alpha (ZAS), and DESI 2024-25, whose evolving-w hints Das II
  mention.
* Local H0 (Riess 2016 in ZAS).
* BBN abundances and Li-7 (ZAS §7).
* High-z quasars (ZAS §6).
* The parameter values in this note (Planck 2018).

Seed selection by MCMC is a fit to these data. Every "good realization" is retrospective.

**Status [Q].**
* ZAS (Model 2, alpha = 0.88, mu = 0.98): chi^2 = 11335.2 against 11334.6 for LCDM on Planck 2015 + BAO,
  for a selected history.
* Das II (Model 1):
  * Supernovae: about 0.14% of seeds fit Pantheon+SH0ES about as well as LCDM, and 0.015% fit better.
  * CMB: the best Model 1 chi^2 is "more than 150.0 greater than ΛCDM". With Lambda suppressed near last
    scattering it is about 20 worse.
  * Favoured variations are "smaller in magnitude than is typical" for everpresent Lambda.
* The ensemble statistic is the model's actual held-out content: P(observed history | alpha). It is
  <= 1e-4 in my converged Model 1 runs and <= 5e-5 in the optimistic linear estimate.

**Possible future (held-out) tests.**
* Dark-energy tomography at z > 2: DESI final, Euclid, Roman, 21-cm intensity mapping (CHIME, HIRAX, SKA;
  listed in Das II §3).
* CMB-S4-class limits on recombination-era dark energy.
* A predicted statistic, for example the two-point function of rho_Lambda(ln a), can be tested only once
  alpha, the start epoch and the recipe are fixed **before** the data.

### 4.7 Gate answers for everpresent Lambda

* **Q1.** Causal-set discreteness (fixed N, Poisson) plus unimodular conjugacy. This is class S/X of the
  `obstruction` escape list. Also a phenomenological stochastic ansatz, the past-light-cone V, and a
  modified Friedmann system that drops the second equation.
* **Q2. Partial.**
  * Scale: fixed as |rho_Lambda| ~ alpha M^2 / sqrt(V) ~ alpha M^2 H^2, with H as input.
  * Sign: random.
  * Mean: assumed zero, the unsolved "first horn" in ADGS's words.
  * Coefficient: alpha is free. The 1-sigma match needs alpha_red = 0.72.
* **Q3.** The spread is stable against heavy thresholds. The mean is stable only as a value statement on
  the total.
* **Q4.** Statistical consequences (a)-(e) are stated. All quantitative comparisons used exposed data.
  The simplest realizations are strained by the CMB and by BBN versus today. They are also strained by
  Barrow's estimate if the patches are independent.
* **Verdict.** The best-developed published proposal that ties Lambda's scale to a principle. It does not
  fix the total renormalized vacuum energy (the mean), and its single-alpha form makes our history
  improbable.

---------------------------------------------------------------------------------------------------

## 5. GRUT closure: "no external clock -> internal time -> unimodular time conjugate to Lambda" [A]

**Who already made it.**
* **Unruh 1989**: unimodular time gives a Schroedinger form, and the small Lambda becomes a question of
  initial conditions [Q abstract].
* **Unruh-Wald 1989**: the conditional-probability interpretation needs a preferred time. They prove that
  for H bounded below no dynamical variable correlates monotonically with Schroedinger time. They then
  seek a nondynamical time parameter and obtain "general relativity with an arbitrary, unspecified"
  cosmological constant. They add that the theory is unlikely to admit sufficiently many observables
  for general spacetimes [Q abstract].
* **Henneaux-Teitelboim 1989**: a cosmic time canonically conjugate to Lambda [Q abstract].
* **Sorkin 1987 and 1994**: the "V = T" rule identifies elapsed 4-volume with the collapse time of the
  sum over histories, which "turns the cosmological constant into an essentially classical constant of
  integration" [Q abstract]. The 1997 and 2007 papers add causal-set number as time [Q].
* **Smolin 2009 §6**: global time [Q].
* **Critique: Kuchar 1991** [Q abstract].

**Step-by-step assessment.**
1. **"A closed whole has no external clock."** For a spatially compact gravitational system this is
   standard: H_total ≈ 0 (`obstruction` §2.1). **But GRUT's FORMAL_CORE DEF 1.1 and 1.6 supply an external t**
   that no subsystem reads. They state that a Page-Wootters clock is an added construction and that
   "nothing here solves the problem of time". Adopting the premise therefore changes GRUT's formalism from
   a unitary closed system to a constrained one. It is not a consequence of closure as GRUT defines it.
2. **"Therefore time is internal (relational)."** Sound within constrained systems. This is standard
   (DeWitt, Page-Wootters). It does not single out a clock.
3. **"The internal time is 4-volume."** Not forced.
   * There are many clocks: matter scalars, volume, York time, Chern-Simons time, unimodular time.
   * In pure GR the 4-volume between slices has no conjugate conserved partner. Making it one requires the
     HT fields or a fixed volume form, which is extra structure (Codex Q1).
   * Different clocks give inequivalent quantum theories. Afshordi-Magueijo note that the canonical
     transformations Lambda -> phi(Lambda) give inequivalent quantizations [Q]. Gielen & Menéndez-Pidal 2020
     make the clock-dependence point (title).
   * Kuchar 1991: cosmological time does not by itself set measurement conditions [Q].
   * This step is a postulate in FORMAL_CORE §7's "clock slot".
4. **"Therefore Lambda is fixed."** Does not follow. With unimodular time, Lambda is the generator of
   T-evolution. The compensated vacuum shift is the phase e^{icT} (§3.1). The kinematic closure structure
   cannot select <Lambda> (lemma, §3.2). The argument reproduces Codex's obstruction one level up. It ends
   in the state-selection slot, where earlier GRUT closure work found no rule (FORMAL_CORE §4 and §7; the
   C1, C2 and CP-12 campaigns). Those negatives are not repeated as new here.

**Verdict.** The chain is **suggestive** as a reason to consider unimodular quantum cosmology. It is
**not sound** as an argument that closure fixes Lambda. Its best-developed published form is Sorkin's, and
that form yields a **spread** (everpresent Lambda), not a value. Nothing in it is GRUT-native.
* A GRUT version would need two explicit postulates: a clock-slot postulate (4-volume time, with HT fields)
  and a state or measure rule with Lambda-dependent gravitational weights.
* Both have published representatives (§3.3). Such a version would have to beat Sorkin, Ng-van Dam and
  Smolin on the gate.

---------------------------------------------------------------------------------------------------

## 6. Scorecard (four-question gate)

| Proposal | Q1 extra structure | Q2 fixes residual (sign, scale)? | Q3 stable under heavy thresholds without retuning? | Q4 independent consequence / exposed data | Verdict |
|---|---|---|---|---|---|
| Classical unimodular / trace-free (Einstein 1919; Weinberg 1989 §7; Ellis et al. 2011) | fixed volume form / TDiff; removes trace mode | **No**: arbitrary integration constant Lambda_hat; rho_eff = (M^2R + tau)/4 from slice data [D] | local equations blind to constants at all orders [D]; multiplier boundary value retuned (Padilla-Saltas) [Q]; nothing predicted to stabilize | none; classically = GR + Lambda [Q] | relocation; does not fix the residual |
| HT / Unruh quantum unimodular (Lambda conjugate to 4-volume time) | global HT fields calT^mu; Lambda a conserved momentum | **No**: A(Lambda) is a free state datum; vacuum shift = phase e^{icT}; Lambda-translation lemma [D, C] | the vacuum shift is a pure phase in unimodular time; distribution of Lambda_eff is a state label | sigma_Lambda sigma_T >= hbar/2; A-M classicality relation uses the observed Lambda as input | obstruction reappears as state selection |
| Lambda-integrated Euclidean measure (Baum 1983, Hawking 1984; Ng-van Dam 1990) | integral over Lambda with weight exp(24 pi^2 M^4/rho_Lambda) | sign + and scale **zero** (not the observed residual); exponent 3.3e122 at the observed value [C]; not normalizable | flux-like label; measure contested (Euclidean conformal factor) [Q/A] | none beyond Lambda = 0+ | fixes the wrong value |
| Ng-van Dam / Smolin 2009 stationary phase | Lambda integral + stationarity in total V4 + mediocrity | formally rho_Lambda ~ history average of matter Lagrangian; **sign and size depend on on-shell L_m choice** [D]; needs finite V4 and measurement theory | same structural class as global sequestering (credit KP; `sequestering` track) | none stated | ill-defined as stated |
| Everpresent Lambda (Sorkin 1987-2007; ADGS 2004; Ahmed-Sorkin 2013; ZAS 2018; Das-Nasiri-Yazdi 2023-24) | causal-set fixed N + Poisson Delta V ~ sqrt(V); unimodular conjugacy; past-cone V; stochastic Model 1/2; Friedmann-only dynamics | **Partial**: spread fixed alpha M^2/sqrt(V) ~ alpha M^2 H^2 (needs alpha_red = 0.72 today [C]); sign random; **mean assumed zero**; alpha free | spread independent of vacuum loops (phase) [D/Q]; mean only by a value statement on the total [D] | fluctuating w with corr sqrt(V1/V2), early DE ~2-2.7 alpha, sign changes, termination; Barrow inhomogeneity bound alpha_red <~ 2e-5 if patches independent [E]; P(our history) <= 1e-4 (MC), 5e-5 (linear) [C]; all quantitative tests used exposed SN/CMB/BAO/BBN/H0 data | best existing scale-fixing proposal; fails on the mean; strained by data |
| Other unimodular state rules (Afshordi-Magueijo 2022; Gielen 2026; Davies-Magueijo 2026; Kaplan-Rajendran 2025; Josset-Perez-Sudarsky 2017, Perez-Sudarsky 2019) | classicality requirement; self-adjoint BCs; final-state post-selection; nonlinear QM; energy-diffusion current | none fixes the observed residual from first principles (abstract level): uses observed Lambda (A-M); leaves Lambda a state choice (Gielen); Lambda = 0 with a chosen final state (D-M); removes Lambda (K-R); dimensional estimate with a free parameter (J-P-S) | not assessed (abstracts) | D-M: phantom w < -1 when extrapolated; K-R: G_N tests; J-P-S: none specific | not evaluated beyond abstracts |
| GRUT closure -> internal time -> unimodular time | none from closure; needs clock-slot postulate (4-volume time, HT fields) + Lambda-dependent state weight | **No** (lemma) | n/a | none | suggestive, not sound; made by Unruh(-Wald) 1989, HT 1989, Sorkin 1987/1994; criticized by Kuchar 1991 |

---------------------------------------------------------------------------------------------------

## 7. Relation to the other tracks

* `obstruction`: its §4 (unimodular is relocation, not refutation), §6.2 (escape classes G, B, S, X,
  R, P) and its R3 and R7 scorecard rows already state the classical result and the coincidence form.
  This track adds four things:
  * the explicit HT canonical form of the identity (vacuum shift = e^{icT});
  * the translation lemma as the reason closure-type rules give a spread and not a mean;
  * the everpresent-Lambda numbers, statistics and Monte Carlo;
  * the assessment of the GRUT time argument.
* `sequestering`: the Ng-van Dam / Smolin relation is a history average of the same type as global
  sequestering. Credit goes to Kaloper-Padilla, and its scoring is in that track.

## 8. Limits

* **Primary 1989 papers.**
  * HT and Unruh are read at abstract level only (Elsevier and APS returned 403). The HT action comes
    from secondary reproductions.
  * Weinberg §7 is read in the authorized Russian translation without equations. Statements are
    paraphrased, and I do not quote the English wording.
  * Unruh-Wald 1989, Kuchar 1991, Sorkin 1994, Ng-van Dam 1990/2001 and DLS 1993: abstracts only.
* **Recent proposals.** Afshordi-Magueijo (partial), Davies-Magueijo, Kaplan-Rajendran, Gielen 2026,
  Josset-Perez-Sudarsky and Perez-Sudarsky are assessed at abstract level. Their scorecard row is
  provisional.
* **ZAS and Das.** ZAS §§6-7 (quasars, Li-7) and Das II §§4-6 (fitting details) were not read. The
  status numbers are quoted from the abstracts and conclusions.
* **Monte Carlo scope.** It reproduces my own implementation of Model 1:
  * stationary start;
  * log grid with 80 steps per e-fold;
  * termination checked at grid points;
  * Planck-2018 Omega_m h^2 and Omega_r h^2.
  It does not reproduce Das et al.'s exact recipe (S_0, V_0, start, steps). Model 2 is not simulated. The
  selection windows are illustrative stand-ins for data likelihoods, not fits.
* **Barrow-type bound.** An order-of-magnitude estimate (uniform horizon sphere, Sachs-Wolfe, dT/T <= 1e-5)
  that is conditional on independent patches.
* **The lemma.** It covers principles built from the kinematic HT structure. It is not a theorem about
  dynamics, measures, or non-Gaussian rules that use Lambda-dependent gravitational input. Those are
  exactly the escape routes listed.
* **Smolin's normalization.** He writes G_ab - Lambda g_ab = 4 pi G T_ab. I did not reconcile it with mine.
  The on-shell-Lagrangian ambiguity (§3.4) does not depend on it.
* **Searches.** Literature searches were arXiv API and INSPIRE queries only. A paper with different
  wording could be missed.
* **Data.** All cosmological numbers are retrospective (Planck 2018). No new data, no blinding, and no
  confirmatory claim.

## 9. Files

* `code/unimodular_checks.py`: sections A-H, numpy, scipy and sympy, seed 20261021, 39 checks, about one
  minute.
* `results/unimodular_checks.json`: every number quoted above, the check list (39/39 pass), and the
  Monte Carlo tables `F_model1_MC`, `F_start_epoch_scan` and `F_grid_dependence`.
* Downloads (not in the repository): scratchpad `um_dl/` and `lit/vac/w.txt` (Weinberg translation).
