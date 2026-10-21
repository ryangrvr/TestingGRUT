# Track `others`: the remaining vacuum-energy proposal classes against the four-question gate

Packet `grut_vacuum_principle_20261021`, track `others`, 21 October 2026. Written under `../../BRIEF.md`.

Labels: [D] derived here, [C] computed by `code/others_checks.py` (record `results/others_checks.json`; 26/26
checks pass), [Q] quoted or paraphrased from a source at the read level in §12, [E] order-of-magnitude
estimate, [A] argued.

Nothing here is a new theory, an empirical validation or a confirmatory test. Every comparison with Planck,
DESI, GW170817 or supernova data is **retrospective**, because those data are already exposed. Prediction
dates are given only to show which published claims were made before a dataset existed. That history does
not make them held-out tests for this program.

Conventions: signature (-,+,+,+); M = (8 pi G)^(-1/2) = 2.435e18 GeV (reduced Planck mass); hbar = c = 1.

The target is the total renormalized rho_vac/M^4. Its observed value, rho_DE0/M^4 = 3 Omega_DE (H0/M)^2,
is 7.16e-121 for h = 0.674 and Omega_DE = 0.685 [C]. Both inputs are retrospective. The value matches
the sibling tracks.

Scope relative to the sibling tracks: the `obstruction` track (§§5-7) has already placed these classes in a
coarse scorecard (rows R4-R7). That scorecard is mostly at abstract or recall level. This track reads the
primary sources more deeply, writes their equations, reproduces their published numbers, and scores each class
on its own. Sequestering and Lombriser's global variation belong to the `sequestering` track. Unimodular gravity
and Henneaux-Teitelboim belong to the `obstruction` track. They appear here only as pointers.

---------------------------------------------------------------------------------------------------

## 0. Bottom line

1. **Every class here escapes the fixed-background energy-origin obstruction.** None is GRUT-native. Each
   escape adds a structure that the shift-invariant matter argument does not contain. The structures are:
   - a rolling scalar that breaks the vacuum's Poincaré invariance (Fab Four, well-tempering);
   - a non-local or acausal filter of the constant mode (ADDG, degravitation);
   - a UV/IR state-counting bound tied to a horizon (CKN, holographic DE);
   - Lambda promoted to a variable, plus a measure (Baum-Hawking, Coleman, no-boundary, landscape);
   - a global variational constraint (Barrow-Shaw, everpresent Lambda, Einstein 1917).
   These are the G, B, S and X items of the obstruction track's minimal escape list.
2. **No class passes Q2 and Q3 together.** No class fixes the sign and scale of the residual from its own
   principle while staying stable under a heavy threshold without retuning.
   - Self-tuning reaches flat space. Well-tempering reaches H_dS = kappa M, where kappa and M are Lagrangian
     inputs [C].
   - ADDG only suppresses: rho_eff = E M^2/(M^2 + Mbar^2). Mbar is fitted to data (4.9e47 GeV for
     E = TeV^4) [C]. At fixed Mbar the residual scales with the vacuum energy E. A 10 TeV scalar threshold
     multiplies it by 17; a 1e16 GeV threshold multiplies it by 1.6e49 [C].
   - Degravitation leaves a residual of order M^2/L^2, with L = 0.70/H0 supplied [C].
   - Holographic DE fixes the scale M^2/R_h^2, but its present Omega_de is an integration constant [Q, C].
   - Euclidean measures fix rho -> 0+. That value is fixed but wrong, and the sign of the on-shell action is
     contested.
   - The landscape gives a distribution, not a value.
3. **Only one class member found here makes an independent numerical prediction beyond Lambda: Barrow-Shaw.**
   Their causal Lambda-variation predicts Omega_k0 = -0.0055 (2 zeta_b) [Q]. I reproduced it as -0.00554 [C].
   With current inputs it becomes -0.0062 to -0.0064 [C]. DESI DR2+CMB gives Omega_K = +0.0023 +/- 0.0011.
   That is 7.7 sigma away for the authors' zeta_b = 1/2, and 2.1 sigma for the sign alone [C, retrospective].
4. **Holographic DE (Li) is a testable w(z) model, not a vacuum-energy solution.**
   - At c = 1, w0 ≈ -0.90 and w rises toward -1/3 in the past. The CPL projection is w_a ≈ +0.44 [C]. The
     DESI DR2+CMB+SN fits have w_a ≈ -0.6 to -1.1 [Q].
   - On DESI DR2 BAO alone, HDE (c = 1) fits as well as LambdaCDM: chi2 10.60 vs 10.55 for 13 data [C]. My
     LambdaCDM BAO fit reproduces DESI's Omega_m = 0.2975 ± 0.0086 as 0.2970 ± 0.0085 [C].
   - The Hubble-scale cutoff gives w = 0 [D]. It is excluded by the BAO shape, Delta chi2 = 1445 [C].
   - Published fits that include the CMB find tension within HDE [Q abstracts].
5. **Two published arguments narrow what any GRUT principle could be.**
   - Bousso-Polchinski: a small present Lambda cannot be set by early dynamics, because it is negligible at
     nucleosynthesis. I compute rho_Lambda/rho_rad = 7e-36 at T = 1 MeV [C].
   - The ADDG/Padilla causality trilemma: a long-distance modification must either cancel late (ruining BBN),
     also cancel short-wavelength sources, or be acausal through a global or future condition [Q].
   - Together these push any candidate toward a whole-history (teleological) functional or a state/measure
     selection [A]. These are exactly the slots where earlier GRUT closure work found no selection rule (§2).
     The "whole universe" premise is therefore not new structure by itself. Each published whole-history
     proposal supplies extra variables and boundary terms, and each fails Q2 or Q4.
6. **Already in the GRUT record (§2); not to be presented as new:**
   - The energy-origin point itself: VACUUM_CLUSTER_MAP.md:19 records D0(ii), "flat-space L -> L + c is an exact
     symmetry; gravity breaks it".
   - Landscape + measure = RELOCATION (rrp/rrp02/RRP_02_VERDICT_01.md:12).
   - Lambda is a datum and a Jacobson-recovery integration constant in GRUT III (GRUT_III/01_PRINCIPLES.md:10, :52).
   - Sequestering k > 0 is conditional only. The all-formulations claim was withdrawn (closure fork).
   - No state or no-boundary selection rule was found (closure fork, s09).
   - Candidate 1's rho_Lambda(mu0) = 0 is a benchmark matching, not a prediction.

## 1. The gate, restated with the minimum each answer must contain

For each proposal I list: the equations; the extra structure (Q1); the residual's sign and scale and whether
they are fixed (Q2); stability when a heavy threshold delta rho = N m^4/(64 pi^2) ln(m^2/mu^2) is added at one
loop (Q3); and the independent consequence together with the exposed data it has already met (Q4). A
reproduction of a published number is marked [C] when the code recomputes it, and [Q] when it is only quoted.

## 2. What earlier GRUT packets already concluded (searched under research/reality_framework and the repository root)

Search terms used (grep -i): sequester, unimodular, cosmological constant, rho_Lambda, vacuum energy, Obadia,
no-boundary, GRUT III, D2. A second pass covered the class names (Fab Four, self-tuning, degravitation,
holographic, Coleman, Baum, wormhole, anthropic, Bousso, landscape, Barrow, Barbour, Mach, Weinberg).
- 105 files under `research/reality_framework` match the first set. Most are Candidate 1 or wall-history
  code, where "vacuum" means the scalar vacuum.
- The root-level hits are VACUUM_CLUSTER_MAP.md, rrp/rrp02, GRUT_III and a few program documents.

| # | Record | What it concluded | Status for this pass |
|---|---|---|---|
| P1 | `VACUUM_CLUSTER_MAP.md` (in initial commit 95ae18f; dated 2026-08-04) | :19 "spine" D0(ii): flat-space L -> L + c is an exact symmetry, and gravity breaks it, so the vacuum term's gravitating premise is load-bearing. :39-45: "120 orders" is convention-dependent (~5 orders, undetermined sign); the regulator-free core is the threshold depths (electron ~1e31, QCD ~1e43-44.5, EW 1e54.675 x rho_obs). :37, :107: Banks relation scored TRADE (-1 parameter, +1 posit). :76: CKN, UV/IR and holographic bounds noted as accepting decoupling but denying the naive bookkeeping of the dimension-zero operator. :88-95: undecided pairs, including (1) is unimodular decoupling stable under quantization, (4) is rho_Lambda typical under a computable non-question-begging measure, (5) w = -1 versus DESI | Codex's obstruction restates D0(ii) in closure language. My electron number (1.7e31 x rho_obs [C]) is consistent with :45. Pair (4) is the anthropic row below. |
| P2 | `rrp/rrp02/RRP_02_VERDICT_01.md` (2026-09-24) | :12 landscape + measure = RELOCATION: the Lambda value is traded for flux integers, a measure, a reference class and a typicality prior. :25 principle II': supplied content is relocated or converted between types, never eliminated. :74 entropic and holographic derivations not yet audited | The landscape verdict stands. This track adds the holographic audit (§5). |
| P3 | `GRUT_III/01_PRINCIPLES.md:10, :52, :84`; ledger s05 | Lambda is a datum. The Jacobson-style GR recovery gives Lambda as an integration constant. GRUT III predictions F1-F6 are shared with LambdaCDM | Same relocation as unimodular gravity (obstruction track §4). |
| P4 | `GRUT_III/10_SECULAR_IR_SURVEY.md`; ledger s05 §3.18; `LITERATURE_PASS_02.md:38` | Graviton secular IR effects in de Sitter (the Tsamis-Woodard screening line): an unreconciled "genuine open window" (reading report :104) versus "not a GRUT one" (author :138). The two-loop IR sign statement is credited to Tsamis-Woodard 1995 | An IR-relaxation route related to class (b). Still unreconciled; not reopened here. |
| P5 | `grut_closure_fork_20261019/` DERIVATION.md:199-215, REPORT.md:85-91, :114-127, CLOSURE_SELECTION_RESPONSE.md:7-12, :28-40, :90-111, :157-177; `grut_closure_control_20261019/STATUS.md` | Global sequestering's <R> = 0 gives k∫a dt > 0 only conditionally. Withdrawn: "all sequestering formulations force k > 0", the sign-only curvature falsifier, and D2 novelty. D2 = the slow-roll coefficient T = H(1-eps)/2pi, which is prior art (Obadia PRD 78, 083532 (2008)); it is not a Lambda claim. C2/P1 restriction selects no state. The no-boundary link to the Euclidean vacuum is a conditional matter-state construction. Maldacena 2024: the no-boundary weight favours too little inflation and too much curvature. The B three-prescription contour control is registered and unexecuted | Class (d) inherits "no state-selection rule". (GRUT III's own "D2" at 02_DERIVATIONS.md:25 is the stress-tensor noise kernel, which is unrelated.) |
| P6 | `grut_candidate_1_quantum_20261018/REPORT.md:7, :98, :134-136`; `grut_issue13_cross_review_20261020/second_agents/derivation_scope/FINDINGS.md:31-35, :73` | rho_Lambda(mu0) = 0 is a benchmark matching. The one-loop scalar vacuum energy is -1.2265e-6 v^4 with no adjustment applied. The H0-envelope exclusion needs rho_vac,tot >= -(2.28 meV)^4, an external calibration tuned to about 1 part in 1e52 | The program's own sector has no vacuum principle. |
| P7 | `grut_closed_whole_theory_20261019/` (s09 ledger; WORKING_THEORY; FORMAL_CORE §4 per the obstruction track) | Closure forced nothing additional. Constants belong to Theta | Baseline. |
| P8 | `GRUT_MINIMUM_CORE.md:19`; `RAI_STRUCTURAL_THEORY_SEARCH.md:156, :183` | rho_Lambda is a measured datum every candidate must reproduce. "The vacuum energy density" is ill-posed without a state choice. Whether vacuum energy gravitates has never been measured | Retained. |

Sibling tracks of this packet: `obstruction` (Claim C escape list G/B/S/X/R/P; scorecard R0-R8) and
`sequestering` (Kaloper-Padilla, Lombriser; Ferreira-Avelino reproduction).

## 3. (a) Self-adjustment

### 3.1 Weinberg's no-go (the baseline theorem)

Primary text: Weinberg RMP 61, 1 (1989). The obstruction track read §§3, 6 and 7 in the authorized Russian
translation; Padilla's lectures (1502.05296) §3.2 were read in full. I cite those readings and did not re-read
the RMP.

The setting is constant fields phi_i and a flat metric, so the residual GL(4) gives L -> det(A) L. Then L =
sqrt(-g) V(phi). The field equations dV/dphi_i = 0 and V = 0 require a tuning. In the degenerate (dilaton) case
L = sqrt(-g) e^{4 phi} V0, and a solution needs V0 = 0 or e^phi -> 0 [Q, via obstruction §5].

This is a theorem that bounds the class, not a proposal. Its loopholes are non-constant adjusting fields (3.2,
3.3), non-local or global structure (§4, §8), and measures (§6, §7).

### 3.2 Fab Four (Charmousis-Copeland-Padilla-Saffin, PRL 108, 051101 (2012), arXiv:1106.2000; read in full)

Equations [Q]: S = ∫sqrt(-g)[V_john G^{mu nu} d_mu phi d_nu phi + V_paul P^{mu nu a b} d_mu phi d_a phi
nabla_nu nabla_b phi + V_george R + V_ringo G_GB] - ∫sqrt(-g) rho_bare + S_m[g; psi]. Matter is minimally
coupled. The self-tuning filter has three parts:
- on the flat (Milne-sliced, kappa < 0) vacuum H^2 + kappa/a^2 = 0, the scalar equation must vanish
  identically;
- the Friedmann equation must depend on phi-dot;
- both must hold through arbitrary jumps of rho_Lambda.
Ringo cannot self-tune alone. George alone is GR or Brans-Dicke with omega = 0, which solar-system tests exclude.

- **Q1.** Extra structure: a time-dependent scalar on the self-tuning vacuum, which breaks Poincaré invariance
  in the scalar sector, plus four free potentials.
- **Q2.** The target is curvature zero, i.e. Minkowski. A non-zero residual needs an extra scale. The proposal
  fixes "0", not the observed residual.
- **Q3.** Classically, any piecewise-constant rho_Lambda is absorbed. Quantum stability is not established. The
  authors say matter loops will probably spoil it; Padilla §7.1 says stability is "by no means guaranteed" [Q].
- **Q4.** Consequences: light-scalar fifth forces that need Vainshtein screening; screening of matter as well
  as vacuum energy (Appleby-Linder 1805.00470 §II item 1, citing Linder 1310.7597 and Copeland-Padilla-Saffin
  1208.3373) [Q]; and tensor speed c_T ≠ 1 from the John, Paul and Ringo-type kinetic mixing.
  - GW170817 bounds |c_T - 1| ≲ 1e-15, which strongly disfavours late-time Horndeski terms that change c_T
    (Ezquiaga-Zumalacárregui PRL 119, 251304 (2017), abstract) [Q].
  - Copeland-Kopp-Padilla-Saffin-Skordis (PRL 122, 061301 (2019); abstract and introduction) find a class
    whose c_T - 1 is proportional to the scalar equation, which rescues it at background level. They then show
    the loophole fails once large-scale inhomogeneities are included [Q].
  - Exposed data: GW170817, solar-system tests, the cosmic expansion history.
- **Verdict:** escapes Weinberg's theorem, not the residual (Q2 fails). Q3 is unproven. Q4 is largely already
  constrained by exposed data.

### 3.3 Well-tempered cosmological constant (Appleby-Linder, JCAP 1807, 034 (2018), arXiv:1805.00470; §§I-III, VII read; §§IV-VI skimmed)

The arXiv id the task left open is **1805.00470**.

Equations [Q]: S = ∫sqrt(-g)[M^2 R/2 + K(phi,X) - G(phi,X) box phi + rho_Lambda], shift-symmetric, with
K = c1 M^3 phi + M^4 A(X) and G = M B(X). Dimensionless variables: psi = phi/M, h = H/M, tau = M t,
mu = M_pl^2/M^2. Degeneracy is imposed so that the scalar and scale-factor equations coincide on h = kappa.
The simplest member is A = const, B = (c1/9kappa^2)(ln psi' + 1/psi').

[C] (sympy, checks B1-B4):
- On h = kappa, the scale-factor equation gives psi'' = 3 kappa psi'. The scalar equation is then satisfied
  identically.
- This agrees with the paper's eq. (28) for F = -3kappa/(1 - psi').
- rho_Lambda enters only the Hamiltonian constraint. There it shifts psi by dpsi/d(M_Lambda^2/M^2) = mu/c1 and
  leaves H unchanged.

Gate:
- **Q1.** Shift-symmetric cubic Horndeski with a tadpole (c1 ≠ 0 is required), plus degeneracy conditions
  (24)-(25) relating A and B. c_T = 1 by construction (alpha_T = 0).
- **Q2.** The de Sitter attractor H_dS = kappa M is set by Lagrangian parameters. The sign is positive by
  construction, but the scale is supplied. The authors leave open why M is such that acceleration happens
  today (§VII) [Q].
- **Q3.** Classically, rho_Lambda jumps, including phase transitions, are absorbed while matter eras are
  preserved [Q]. Shift symmetry is offered as quantum control, but quantum corrections are left for future
  work (§VII) [Q]. Not established.
- **Q4.** A transient w ≠ -1 during a vacuum phase transition, a de Sitter attractor, and possible Vainshtein
  screening. GW170817 was a design input, so it is exposed data and not a test. No independent number is
  fixed.
- **Verdict:** removes a large rho_Lambda classically but keeps a free scale (Q2 fails; Q3 open).

## 4. (b) Degravitation and non-local IR modification

### 4.1 Arkani-Hamed-Dimopoulos-Dvali-Gabadadze (ADDG), hep-th/0209227 (preprint, no journal reference in INSPIRE; §§1-3.2 read, rest skimmed)

Equations [Q]: M^2 (1 + F(L^2 box)) G_mn = T_mn. Here F is a filter with F -> 0 for short wavelengths and
F(0) = Mbar^2/M^2 ≫ 1. In the L -> ∞ limit:

  M^2 G_mn - (1/4) Mbar^2 g_mn Rbar = T_mn,   Rbar = ∫sqrt(g) R / ∫sqrt(g)            (ADDG eq. 8)

Rbar is a spacetime average dominated by the asymptotic future. For T_mn = -E g_mn + matter in an
asymptotically de Sitter history, the trace gives (our sign convention) Rbar = 4E/(M^2 + Mbar^2). So

  rho_eff = E M^2/(M^2 + Mbar^2)                                                        [D from eq. 12]

Equivalently, M^2 G = T - (1/4) g [Mbar^2/(M^2+Mbar^2)] Tbar (eq. 15). This is the unimodular trace-free
equation plus one global equation. Setting the coefficient to 1 gives the Tseytlin-type global subtraction
(their ref. [4], PRL 66, 545 (1991)), the precursor of Kaloper-Padilla global sequestering [Q]. ADDG argue that
SM loops do not renormalize Mbar, because Gamma_SM is local [Q]. They state that the theory is maximally
acausal (§3.2).

- **Q1.** A global, acausal spacetime average of R, and a new IR parameter Mbar.
- **Q2.** Suppression, not cancellation.
  - The residual has the sign of E. Its scale is E (M/Mbar)^2, with Mbar fitted to data.
  - [C] Fitting gives Mbar = 4.85e47 GeV for E = (1 TeV)^4 (ADDG quote ~1e48) and 2.9e78 GeV for E = M^4
    (ADDG quote ~1e80).
- **Q3.** Fails.
  - At fixed Mbar, the residual is linear in E. [C] A single real scalar threshold with |delta E| = m^4/(64
    pi^2) multiplies the fitted residual by 1.0016 for m = 1 TeV, 16.8 for m = 10 TeV, and 1.6e49 for
    m = 1e16 GeV.
  - The value is stable only if nothing heavier than the fitted E exists, or if Mbar is retuned.
- **Q4.** GR is recovered wherever Rbar = 0 (localized or ever-decelerating histories). Beyond that, there is a
  claimed "self-terminating inflation" and no computed separate observable. Exposed data: the expansion
  history (Mbar is fitted to it).
- **Verdict:** a non-local suppression with a fitted parameter. Its x -> 1 limit belongs to the
  sequestering lineage (credit Tseytlin, ADDG, Kaloper-Padilla).

### 4.2 Dvali-Hofmann-Khoury (DHK), PRD 76, 084006 (2007), hep-th/0703027 (§§I-II, III start, XI read; rest skimmed)

Equations [Q]: G_N^{-1}(L^2 box) G_mn = 8 pi T_mn, written as (1 - m^2(box)/box) G_mn = T_mn/(2M^2), with
m^2(box) ∝ L^{2(alpha-1)} box^alpha.
- Unitarity requires alpha >= 0. Degravitation requires alpha < 1.
- A consistent filter is the helicity-2 part of a massive (alpha = 0) or resonance graviton of mass or width
  1/L. It therefore carries extra polarizations, a vDVZ/Vainshtein structure and strong coupling.
- DGP (alpha = 1/2) does not degravitate. In massive gravity, phase transitions can excite a BD-like
  instability.
- §XI: filtering can be trusted only while HL ≫ 1. The most interesting outcome would be extra states that
  stop degravitation and leave a residual of order L^-2 (in M^2 units). That outcome is offered as a
  possibility, not derived.

- **Q1.** An IR graviton mass or width 1/L, with extra polarizations.
- **Q2.** Not fixed: rho_res ~ M^2/L^2 with L supplied. [C] L = 0.70/H0 reproduces rho_obs; the O(1)
  coefficient is not computed in DHK.
- **Q3.** A constant source of any loop order is filtered if the filter itself is stable. The residual, if it
  is M^2/L^2, would not depend on E. But that endpoint is conjectural, and the radiative stability of L and of
  the strong-coupling scale is not shown [Q/A].
- **Q4.** Graviton mass or width ~ H0, modified large-scale gravity and Vainshtein-screened fifth forces. No
  viable complete theory with computed observables was found at the read level achieved. Causality trilemma:
  see §9.
- **Verdict:** relocates the value to a supplied IR scale (Q2 fails; Q3 not shown).

Related GRUT record: the Tsamis-Woodard quantum-gravitational screening line is GRUT III's unreconciled IR
"open window" (P4). It was not reopened.

## 5. (c) Holographic and IR-cutoff dark energy

### 5.1 Cohen-Kaplan-Nelson (CKN), PRL 82, 4971 (1999), hep-th/9803132 (read in full)

Equation [Q]: an EFT with UV cutoff Lambda in a box of size L is trusted only if L^3 Lambda^4 ≲ L M_P^2
(eq. 2). That excludes states inside their own Schwarzschild radius. The bound is stronger than Bekenstein's
L^3 Lambda^3 ≲ L^2 M_P^2.

- [C] For L = 1/H0, Lambda_UV = 1.9 meV (reduced M) or 4.2 meV (non-reduced). CKN quote ~10^-2.5 eV.
- CKN say explicitly that this does not predict the cosmological constant's value, because a constant can
  always be added. They also say it removes the need for fine-tuning [Q].
- Their separate consequence is a minimum QFT uncertainty, e.g. delta(g-2) ~ (alpha/pi)(m_e/M_P)^{2/3}
  ~ 1e-15 alpha/pi [Q].

- **Q1.** A UV/IR state-counting bound: a non-extensive Hilbert space.
- **Q2.** Not fixed (CKN's own words).
- **Q3.** CKN argue that matching between cutoffs never exposes a quartic tuning [Q]. The finite m^4 ln m
  threshold terms of particles that exist are a separate matter.
  - [C] The electron's |delta rho| ≈ 1.7e31 rho_obs. The top quark's is 6.7e53 rho_obs.
  - These are scheme-independent in the sense of P1. The IR bound does not remove them [A].
- **Q4.** The QFT-precision floor (tiny, untested), and the UV/IR mixing pattern.
- **Verdict:** a degree-of-freedom counting argument that explains the natural scale M^2/L^2, not the value.

### 5.2 Hubble cutoff (Hsu, PLB 594, 13 (2004), hep-th/0403052, abstract; Li §1)

With L = 1/H, rho_de = 3c^2 M^2 H^2 and rho_m = 3(1-c^2) M^2 H^2. Their ratio is constant, so w = 0 [D, check E3].
- [C] Against DESI DR2 BAO, the total ∝ a^-3 shape gives Delta chi2 = +1445 relative to LambdaCDM. Excluded
  (retrospective).

### 5.3 Li holographic dark energy, PLB 603, 1 (2004), hep-th/0403127 (read in full)

Equations [Q]:
- rho_de = 3c^2 M^2 R_h^{-2}, with R_h = a ∫_t^∞ dt'/a the future event horizon.
- In ln a: dOmega/dln a = Omega(1-Omega)(1 + 2 sqrt(Omega)/c).
- w = -1/3 - 2 sqrt(Omega)/(3c).

Results:
- [C] Reproduces Li's eq. (29): w = -0.903 + 0.104 z for Omega = 0.73 and c = 1.
- Li argues for c = 1 from the second law and de Sitter symmetry. c < 1 is phantom.
- Li says the asymptotic constant is a boundary condition, or equivalently the initial Omega is an initial
  condition, which he attributes to inflationary e-folds [Q].
- [C] With the same c, different Omega_de0 values solve the same equation (check E5). The present value is an
  input.

Gate:
- **Q1.** CKN saturation plus a global, future-dependent IR length (the event horizon) [A: teleological].
- **Q2.** Scale and positive sign: yes. Value: no. Omega_de0 is an integration constant, and c is a parameter.
- **Q3.** HDE is postulated to be the whole dark energy. The m^4 threshold terms (P1; §5.1) must be removed by
  something else. Not addressed.
- **Q4.** A definite w(z) shape for given (c, Omega_m0). All tests below are retrospective.
  - [C] CPL projection over a ∈ [0.3, 1]: for c = 1, w0 ≈ -0.92 and w_a ≈ +0.44. For c = 0.6, w0 ≈ -1.30 and
    w_a ≈ +0.92.
  - In every case w_a > 0, because w -> -1/3 in the past. DESI DR2 + CMB + {Pantheon+, Union3, DESY5} instead
    give w0 = -0.838, -0.667, -0.752 and w_a = -0.62, -1.09, -0.86 (Table 5) [Q].
  - [C] DESI DR2 BAO alone (Table 4: 13 numbers with per-bin D_M-D_H correlations; H0 r_d profiled; flat):

    | Model | Omega_m | chi2 | Note |
    |---|---|---|---|
    | LambdaCDM | 0.2970 [0.2886, 0.3057] | 10.55 | validates against DESI 0.2975 ± 0.0086 |
    | HDE, c = 1 | 0.272 | 10.60 | |
    | HDE, free c | 0.272 | 10.39 | best c = 0.95 |

    c ≲ 0.6 is disfavoured with profile Delta chi2 ≳ 9. BAO alone does not separate HDE from LambdaCDM.
  - Published joint fits:
    - with CMB, HDE wants c < 1 and phantom crossing beyond 10 sigma (2411.08639) [Q abstract];
    - ACT DR6 versus DESI DR2 + DESY5 shows an early/late tension within HDE, and the evidence favours
      LambdaCDM (2509.02945) [Q abstract].
  - Exposed data: all of the above.
- **Verdict:** a falsifiable phenomenological DE model, not a principle that fixes the total vacuum energy.
  Q2 and Q3 fail. Q4 exists but has met only exposed data, with mixed to adverse results.

## 6. (d) Euclidean, wormhole and no-boundary selection

On-shell Euclidean de Sitter (S^4 of radius r^2 = 3/Lambda) [D, sympy F1]:

  I_E = -(1/16 pi G) ∫sqrt(g)(R - 2Lambda) = -3pi/(G Lambda) = -24 pi^2 M^4/rho_Lambda.

- **Baum (PLB 133, 185 (1983)) and Hawking (PLB 134, 403 (1984)), abstracts** [Q]:
  - Lambda becomes an effective variable, via a scalar minimum (Baum) or a 3-form field strength /
    topological fluctuations (Hawking).
  - The weight exp(-I_E) = exp(+24 pi^2 M^4/rho) is largest as rho -> 0+.
  - [C] At rho_obs the exponent is 3.31e122, which matches the obstruction track's ln d. So
    ln[P(rho_obs/2)/P(rho_obs)] = 3.3e122: the observed value is maximally disfavoured against any smaller
    positive value.
- **Coleman (NPB 310, 643 (1988), abstract; mechanism via Padilla 1502.05296 §5, read in full)** [Q]:
  - Wormholes shift couplings by alpha. Then P(alpha) ∝ exp(-(1/2) D alpha alpha + Z(lambda + alpha)), with
    Z ≈ exp(3/(8 G^2 Lambda)) = exp(24 pi^2 M^4/rho).
  - The result is a double exponential at rho -> 0+, claimed exact for all lower-energy interactions.
  - Problems: P is also peaked at |xi| -> ∞ for the curvature-squared couplings, and the giant-wormhole
    catastrophe (Fischler-Susskind, Kaplunovsky) [Q, Padilla §5].
- **Duff, PLB 226, 36 (1989), abstract** [Q]: for a 3-form potential, the coefficient of sqrt(-g) in the action
  and the Lambda in the field equations differ by a sign. This invalidates Hawking's argument. Bousso-Polchinski
  §2.1 describes the needed boundary term [Q, read]. The sign of the measure is therefore contested.
- **No-boundary and tunnelling wavefunctions** [A]: |Psi|^2 ~ exp(±24 pi^2 M^4/V). The no-boundary sign favours
  the smallest V; the tunnelling sign favours the largest. Neither gives 7e-121. Maldacena 2024 (2403.10510,
  abstract; already in the GRUT closure fork P5): the no-boundary weight predicts too little inflation and too
  much curvature. The GRUT closure work found no state-selection rule, and the B contour control is unexecuted
  (P5).

Gate:
- **Q1.** Lambda made a variable (3-form flux, wormhole alpha, potential minimum), plus a Euclidean measure or
  state.
- **Q2.** Fixes sign and scale as rho = 0+. That is the wrong value.
- **Q3.** The selection acts on the total effective Lambda, so loop shifts are absorbed [A]. The measure's sign
  and the conformal-factor problem are contested.
- **Q4.** Lambda = 0, contradicted by the measured acceleration. Coleman in principle also predicts other
  couplings, but they are not computable.
- **Verdict:** fixes a value, and the value is excluded (retrospective since 1998). The measure is contested.

## 7. (e) Anthropic and landscape baseline

- **Weinberg, PRL 59, 2607 (1987), abstract** [Q]: for non-recollapsing universes, the anthropic bound
  is that Lambda must not prevent gravitationally bound structures.
  - [D, sympy G1] Spherical collapse in matter + Lambda: turnaround requires Lambda/3 <= K^3/(27 (GM)^2). With a
    growing-mode shell (K = (5/3) delta H^2 r^2), this gives rho_V <= (500/729) delta^3 rho_m. That is the
    standard spherical-collapse coefficient of the Weinberg line; I derived it here and did not check it
    against the Weinberg or MSW text.
  - [E] For a fluctuation collapsing at z_c = 4 (EdS extrapolation delta_0 = 1.686 (1+z_c)), rho_V,max ≈ 411
    rho_m0 ≈ 189 rho_obs. Bousso-Polchinski quote the window -1e-120 < lambda < 1e-118 (upper bound about 100
    times observed) [Q].
  - A crude flat prior gives P(rho <= rho_obs) ~ 0.5% [E]. Martel-Shapiro-Weinberg (ApJ 492, 29 (1998),
    astro-ph/9701099, abstract) weight by the collapsed fraction and find rho_V up to ~3 rho_m "reasonably
    likely" [Q].
- **Bousso-Polchinski, JHEP 0006, 006 (2000), hep-th/0004134** (§§1, 2.1, 2.4, 3.3, 4 read) [Q]:
  - lambda = lambda_bare + (1/2) Σ_i n_i^2 q_i^2 over J four-form fluxes.
  - The spacing near zero is Delta lambda_min = D Π q_i / (omega_{J-1} |2 lambda_bare|^{J/2-1}).
  - [C] For J = 100, lambda_bare ~ 1 and D = 1: q = 10^-1.43 (BP: ~10^-1.6; same order).
  - Membrane nucleation populates the values. Selection is anthropic in the weak sense.
  - BP's own conclusion: a small present value cannot come from dynamics alone, because at nucleosynthesis the
    present Lambda is negligible. [C] rho_Lambda/rho_rad = 7.1e-36 at T = 1 MeV and 7.1e-32 at T = 0.1 MeV
    (BP's "1e-30" is looser).
- **Causal entropic principle** (Bousso-Harnik-Kribs-Perez, PRD 76, 043513 (2007), hep-th/0702115, abstract)
  [Q]: claims excellent agreement with the observed value under an entropy-production weighting. The
  observers-per-baryon weighting prefers larger values. The result is measure-dependent.

Gate:
- **Q1.** A discretuum (fluxes and compactification), a population mechanism, a measure, a reference class and
  a prior.
- **Q2.** No dynamical value. It gives a distribution, with the sign weakly favoured positive by the asymmetric
  window and the scale ~ the window edge.
- **Q3.** Stable. The flux quanta are not locally renormalized. A loop shift of lambda_bare only moves the
  shell, and the density of vacua near zero is unchanged (BP: a larger |lambda_bare| makes cancellation easier
  for fixed charges) [Q/A].
- **Q4.** Correlated, measure-dependent predictions.
  - Open curvature from CDL tunnelling (Freivogel et al., JHEP 0603, 039 (2006), hep-th/0505232, abstract) [Q].
    DESI DR2+CMB Omega_K = +0.0023 ± 0.0011 (exposed) is consistent in sign, but no measure-free amplitude
    exists.
  - Bubble collisions: not assessed.
- **Prior GRUT record:** RELOCATION (P2); undecided pair (rho_Lambda, typicality) (P1).
- **Verdict:** the baseline that explains the order of magnitude statistically but does not fix the value
  dynamically. Not new.

## 8. (f) Global, Machian and whole-universe constraints

### 8.1 Barrow-Shaw (PRL 106, 101302 (2011), arXiv:1007.3086; GRG 43, 2555 (2011), arXiv:1105.3105; both read in full)

Principle [Q]: promote the bare lambda to a field, restrict the action to the observer's past light cone
(boundary terms: GHY on the light cone; initial-surface data independent of Lambda), and impose
dS_cl/dLambda = 0. This is equivalent to V_M = (1/2)∮ ... (a 4-volume = boundary term relation), which gives
Lambda ~ t_U^-2. In FRW with baryon on-shell Lagrangian L_b = -zeta_b rho_b (zeta_b ≈ 1/2 in the chiral bag
model):

  k = kappa ∫_0^{tau0} (tau0-tau)^3 a^4 zeta_b rho_b A dtau / ∫_0^{tau0} a^2 (tau0-tau)^2 (4(tau0-tau) + 6 tau) A dtau,
  A(tau) = (a H/6) ∫_0^tau dtau*/H^2.                                           (BS eq. 3)

Here k > 0, and Omega_k0 = -k/(a0^2 H0^2).

[C] check H1: with their inputs (Omega_Lambda = 0.73, Omega_b = 0.0423, T_CMB = 2.725 K, zeta_b = 1/2;
h = 0.70-0.72), Omega_k0 = -0.00554. They print -0.0055.

[C] Current inputs, evaluated as written:

| Inputs | Omega_k0 |
|---|---|
| Planck-2018-like (Omega_Lambda = 0.6847, Omega_b = 0.0493, h = 0.6736) | -0.00637 |
| DESI DR2+CMB LambdaCDM (Omega_m = 0.3027, omega_b = 0.02218, h = 0.6817) | -0.00619 |

DESI DR2+CMB in LambdaCDM+Omega_K gives Omega_K = +0.0023 ± 0.0011 (Table 5) [Q]. The Barrow-Shaw history is
GR with a constant Lambda, so that is the appropriate fit.
- Tension: 7.7 sigma for zeta_b = 1/2, and 7.1 sigma for their printed -0.0055.
- For the sign alone (zeta_b -> 0+, which makes k -> 0+): 2.1 sigma [C].
- This is retrospective for us. The prediction itself was published before Planck 2013 and DESI.

Gate:
- **Q1.** A global, causal-past variational constraint on a Lambda field, plus boundary-term choices.
- **Q2.** The scale Lambda ~ t_U^-2 comes only as an order of magnitude (positive). The precise statement is a
  relation between Lambda, k and tau0: Lambda is traded for curvature, not predicted absolutely.
- **Q3.** The constraint selects the total effective Lambda, so loops in lambda_bare are absorbed [A].
  However, the relation depends on the hadronic input zeta_b and on the boundary-term prescription [Q]. No
  heavy-threshold stability analysis of zeta_b or of the boundary terms was found.
- **Q4.** Yes: Omega_k0 = -0.0055 (zeta_b/0.5), a genuinely separate observable. Against exposed data it is
  disfavoured at 7.7 sigma, and the sign at 2.1 sigma.
- **Verdict:** the only class member found that makes an independent numerical prediction. That prediction
  is disfavoured by data now exposed. (I did not survey published critiques; INSPIRE citation queries failed,
  see §11.)

### 8.2 Everpresent Lambda (Sorkin; Ahmed-Dodelson-Greene-Sorkin PRD 69, 103523 (2004), astro-ph/0209274; follow-ups)

Read level: abstracts of the original (astro-ph/0209274), Zwane-Afshordi-Sorkin (CQG 35 (2018), 1703.06265),
Das-Nasiri-Yazdi Part I (JCAP 10 (2023) 047, 2304.03819) and Part II (JCAP 10 (2024) 076, 2307.13743).

Equations [Q]: unimodular gravity makes Lambda conjugate to 4-volume V. Causal-set Poisson fluctuations give
Delta V ~ sqrt(V), hence Lambda ~ ±V^{-1/2} (Planck units), fluctuating in sign with a vanishing mean.
- [E] With V ~ H0^-4, rho ~ M^2 H0^2 ≈ 0.49 rho_obs.

- **Q1.** Unimodular gravity plus fundamental discreteness.
- **Q2.** Scale: yes (≈ ambient density at every epoch). Sign: fluctuating. The **mean is postulated to
  vanish**, so the large-contribution problem is moved into that postulate [A].
- **Q3.** Rests on the zero-mean postulate. Not shown.
- **Q4.** Stochastic, time-varying Lambda tracking the dominant density.
  - Fits SN Ia for a small fraction of seeds.
  - Struggles with the CMB unless Lambda is suppressed near last scattering.
  - Observed H0 and Omega_Lambda are not typical outcomes [Q abstracts].
  - Exposed data: SN, CMB.
- **Verdict:** scale-linking with a supplied mean (credit Sorkin et al.). Q2 is partial; Q4 is under pressure.

### 8.3 Einstein 1917 static universe (recall; not re-read)

Lambda = 4 pi G rho_m = 1/R^2 fixes Lambda from a global matter density in a static closed universe. The
solution is unstable (Eddington) and is contradicted by the expansion. This is the historical precedent: a
"Machian" global constraint that did fix Lambda, in a falsified cosmology.

### 8.4 Zero total energy / Hamiltonian constraint (Tryon 1973, recall)

A closed universe has total energy zero because H ≈ 0. The obstruction track (§2, checks C1-C3) shows the
constraint fixes the generator's value but leaves Lambda a parameter related to state constants. Nothing to
add. Credit there.

### 8.5 Barbour-type scale-invariant (conformal) gravity (Anderson-Barbour-Foster-O Murchadha, CQG 20, 1571 (2003), gr-qc/0211022; gr-qc/9911071; abstracts only)

Best matching with respect to scale removes the volume as a dynamical variable [Q]. At the read level
achieved, I found no determination of rho_vac. Not scored beyond "no value located".

### 8.6 Pointers credited elsewhere

- Lombriser global M_Pl variation (Omega_Lambda = 0.704 with a location prior): `sequestering` track.
- CosMIn (Padmanabhan-Padmanabhan, IJMPD 22, 1342001 (2013), 1302.3226, abstract): Lambda L_P^2 ~ 3.4e-122 from
  the postulate N ≈ 4 pi plus inflation and equality inputs [Q]; `obstruction` track row R7.
- Henneaux-Teitelboim and unimodular gravity: `obstruction` track §4.

## 9. Two cross-class constraints and what they imply for a GRUT proposal [A]

1. **Bousso-Polchinski dynamics argument** [Q + C]. A mechanism acting at or before nucleosynthesis cannot
   know the present Lambda, because rho_Lambda/rho_rad ~ 7e-36 at 1 MeV. A dynamical selection must therefore
   act late, or be global in time, or the value must be selected (by measure or state) rather than evolved.
2. **ADDG/Padilla causality trilemma** [Q, ADDG §3.2; Padilla §7.2]. A causal IR modification cannot tell
   "constant" from "long wavelength" until t ≳ L. So it either:
   - (i) cancels the vacuum late (incompatible with BBN), or
   - (ii) also cancels short-wavelength sources (screens matter; compare the Fab Four's matter screening), or
   - (iii) is acausal through a global or future condition (ADDG's Rbar, Li's event horizon, Barrow-Shaw's
     light-cone variation, global sequestering's spacetime averages).
3. **Implication.** A GRUT "whole closed system" principle that escapes the energy-origin obstruction must be
   of type (iii), a whole-history functional, or a state/measure selection. The closure work already recorded
   that "nothing outside" does not by itself supply either (P5, P7; obstruction track §6.3). Each published
   member of those types has been scored above or in the sibling tracks. Each needs extra variables, boundary
   terms or a measure, and fails Q2 or Q4.

   A future proposal would count as GRUT-native only if it uses an assumption none of these papers use. It
   would also have to state a separate observable before its numbers are computed. Barrow-Shaw's curvature
   prediction shows the form such an observable takes, and also how it can fail.

## 10. Scorecard

| Class / proposal | Q1 extra structure | Q2 fixes residual (sign, scale)? | Q3 stable under heavy threshold, no retuning? | Q4 independent consequence; exposed data | Verdict |
|---|---|---|---|---|---|
| Weinberg no-go (baseline) | none (theorem on local constant-field adjusters) | n/a | n/a | n/a | bounds the class; loopholes = non-constant fields, non-local/global, measures |
| Fab Four self-tuning | rolling scalar breaking vacuum Poincaré invariance; 4 Horndeski potentials | tunes to 0 (Minkowski); no residual | classical yes; quantum doubtful (authors) | fifth forces, matter screening, c_T; GW170817 exposed | escapes theorem, not residual |
| Well-tempering | shift-symmetric cubic Horndeski + tadpole + degeneracy | H_dS = kappa M supplied (sign +) | classical yes [C]; quantum open | transient w ≠ -1, dS attractor; c_T = 1 by design (exposed) | free scale; Q2 fails |
| ADDG non-local Rbar | acausal spacetime average; Mbar | suppresses: rho_eff = E M^2/(M^2+Mbar^2), Mbar fitted | no: residual ∝ E (x16.8 for a 10 TeV scalar) [C] | none computed; expansion history exposed | suppression with fit; sequestering precursor |
| DHK degravitation | massive/resonance graviton 1/L, extra polarizations | rho_res ~ M^2/L^2, L = 0.70/H0 supplied [C] | filter acts on all constants; endpoint and stability not shown | IR-modified gravity at L ~ 1/H0; not computed | relocates to L |
| CKN UV/IR bound | non-extensive state counting | no (authors) | quartic piece only; m^4 thresholds (e: 1.7e31 rho_obs) untouched | QFT precision floor ~1e-15 alpha/pi (untested) | explains natural scale, not value |
| Hubble-cutoff HDE | L = 1/H | w = 0 [D] | — | excluded by DESI BAO shape, Delta chi2 = 1445 [C] | excluded |
| Li HDE (event horizon) | CKN saturation + future event horizon (teleological) | scale and sign yes; Omega_de0 integration constant, c free | assumes other vacuum terms removed | w(z): c = 1 gives w_a ≈ +0.44 vs DESI w_a < 0; BAO-only no discrimination [C]; CMB-joint fits adverse [Q]; all exposed | testable DE model, not a vacuum principle |
| Baum-Hawking Euclidean | Lambda variable (3-form / scalar min) + Euclidean measure | yes: rho -> 0+ (wrong) | selection on total, loops absorbed; measure sign contested (Duff) | Lambda = 0 contradicted (exposed since 1998) | fixed value excluded |
| Coleman wormholes | wormhole alpha-parameters + double-exponential measure | rho -> 0+ | claimed exact below wormhole scale; giant-wormhole catastrophe | other couplings fixed in principle, not computable; Lambda = 0 excluded | fixed value excluded |
| No-boundary / tunnelling | wavefunction prescription + contour | 0+ or largest V; no 7e-121 | prescription-dependent | too little inflation, too much curvature (Maldacena 2024) | no value; GRUT closure found no selection rule |
| Anthropic / landscape (Weinberg; BP) | discretuum + population + measure + reference class | distribution; window ~ -1e-120 to 1e-118, upper edge ~100-190 rho_obs [C/E] | yes (flux quanta, dense spectrum) | open curvature, bubble collisions; measure-dependent; Omega_K exposed | baseline; RELOCATION (GRUT P2) |
| Barrow-Shaw causal Lambda variation | Lambda field + past-light-cone action + boundary terms | order t_U^-2; precise value only as Lambda(k, tau0) relation | total Lambda selected; depends on zeta_b and boundary terms | Omega_k0 = -0.0055 (2 zeta_b); updated -0.0062; DESI+CMB +0.0023 ± 0.0011: 7.7 sigma (sign 2.1 sigma) [C] | only independent prediction found; disfavoured (retrospective) |
| Everpresent Lambda (Sorkin) | unimodular + causal-set discreteness | scale ~ H^2 [E]; sign random; mean postulated zero | rests on zero-mean postulate | stochastic w; SN fits for few seeds, CMB difficult [Q] | scale-linking with supplied mean |
| Einstein 1917 static | global static closure | Lambda = 4 pi G rho_m (+) | matter-defined | static universe falsified | historical precedent only |
| Zero total energy / H ≈ 0 | reparametrization constraint | no (obstruction track) | n/a | none | no value |

## 11. Limits

- **Read levels** are in §12. These were not re-read:
  - Weinberg 1989, taken from the obstruction track and from Padilla.
  - Hawking 1984, Baum 1983, Coleman 1988, Duff 1989, Weinberg 1987, MSW 1998: abstracts only.
  - Einstein 1917 and Tryon 1973: recall.
  - Later critiques of Barrow-Shaw: not surveyed. INSPIRE "refersto" queries returned unfiltered results.
- **DESI BAO fit**: BAO only, diagonal blocks with per-bin D_M-D_H correlations, flat geometry. H0 r_d is
  profiled. Li's HDE equation is used without radiation (negligible at z < 2.4; the LambdaCDM fit changes by
  3e-4 in Omega_m). No SN, no CMB, no MCMC. This is not a likelihood analysis and not a new result: published
  HDE-DESI analyses exist (§5.3).
- **CPL projection** of HDE is a least-squares fit of w(a) over a ∈ [0.3, 1]. It is not a posterior.
- **Barrow-Shaw**: I implemented their printed eq. (3) on a flat background, which is first order in k. It
  reproduces their number, but I did not re-derive the equation. The radiation and neutrino treatment
  (N_eff = 3.046, massless) is mine. zeta_b is their hadronic-model input. The DESI Omega_K tension uses the
  LambdaCDM+Omega_K Gaussian summary only.
- **Thresholds**: |delta rho| = N m^4/(64 pi^2) with the log factor set to 1. These are order-of-magnitude
  magnitudes and scheme-dependent in sign.
- **Anthropic window**: EdS extrapolation at z_c = 4 [E]. It is not the MSW collapsed-fraction distribution.
- **Not assessed**: the radiative stability of Horndeski self-tuning functions beyond the authors' statements;
  the strong-coupling scale of degravitating theories; the bubble-collision phenomenology of the landscape.
- No new theory, no empirical validation, no confirmatory test. All data comparisons are retrospective.

## 12. Sources and read levels

| Source | Identifier | Read level | Used for |
|---|---|---|---|
| Charmousis, Copeland, Padilla, Saffin | PRL 108, 051101 (2012), arXiv:1106.2000 | full | Fab Four equations, self-tuning filter, loop caveat |
| Appleby, Linder | JCAP 1807, 034 (2018), arXiv:1805.00470 | §§I-III, VII; §§IV-VI skimmed | well-tempering equations; free H_dS; open quantum corrections |
| Copeland, Kopp, Padilla, Saffin, Skordis | PRL 122, 061301 (2019), arXiv:1810.08239 | abstract + introduction | c_T loophole fails with inhomogeneities |
| Ezquiaga, Zumalacárregui | PRL 119, 251304 (2017), arXiv:1710.05901 | abstract | GW170817 constraint |
| Padilla, lectures | arXiv:1502.05296 | §5, §7.1-7.2 (full sections); §3 via obstruction track | Coleman mechanism, catastrophe; Fab Four stability; causality trilemma |
| Weinberg | RMP 61, 1 (1989) | via obstruction track (§§3, 6, 7 translation) | no-go |
| Arkani-Hamed, Dimopoulos, Dvali, Gabadadze | hep-th/0209227 (preprint) | §§1-3.2 full; rest skimmed | non-local Rbar equation; Mbar numbers; acausality |
| Dvali, Hofmann, Khoury | PRD 76, 084006 (2007), hep-th/0703027 | §§I, II, III start, XI | filter, extra polarizations, residual ~ L^-2 |
| Cohen, Kaplan, Nelson | PRL 82, 4971 (1999), hep-th/9803132 | full | UV/IR bound; does not predict the value; g-2 floor |
| Hsu | PLB 594, 13 (2004), hep-th/0403052 | abstract (+ Li's restatement) | w = 0 for the Hubble cutoff |
| Li | PLB 603, 1 (2004), hep-th/0403127 | full | HDE equations, w(z), integration-constant statement |
| Li T.-N., Li Y.-H., Du, Wu, Feng, Zhang (HDE after DESI 2024) | EPJC 85, 608 (2025), arXiv:2411.08639 | abstract | HDE + CMB: c < 1, crossing |
| Wu, Li, Du, Zhang (HDE/RDE vs ACT DR6 + DESI DR2) | Chin. Phys. C 50, 045105 (2026), arXiv:2509.02945 | abstract | early/late tension; evidence favours LambdaCDM |
| DESI DR2 Results II | PRD 112, 083515 (2025), arXiv:2503.14738 | Tables 4, 5, 6 and eqs. (26)-(28) (section; local HTML copy in the scratchpad) | BAO data; w0wa; Omega_K |
| Baum | PLB 133, 185 (1983), DOI 10.1016/0370-2693(83)90556-7 | abstract (INSPIRE) | minimum Euclidean action |
| Hawking | PLB 134, 403 (1984), DOI 10.1016/0370-2693(84)91370-4 | abstract (INSPIRE) | Lambda probably zero |
| Coleman | NPB 310, 643 (1988), DOI 10.1016/0550-3213(88)90097-1 | abstract (INSPIRE) + Padilla §5 | wormholes |
| Duff | PLB 226, 36 (1989), DOI 10.1016/0370-2693(89)90284-0 | abstract (INSPIRE) | sign objection |
| Maldacena | arXiv:2403.10510 | abstract | no-boundary curvature problem |
| Weinberg | PRL 59, 2607 (1987), DOI 10.1103/PhysRevLett.59.2607 | abstract (INSPIRE) | anthropic bound |
| Martel, Shapiro, Weinberg | ApJ 492, 29 (1998), astro-ph/9701099 | abstract | likely values |
| Bousso, Polchinski | JHEP 0006, 006 (2000), hep-th/0004134 | §§1, 2.1, 2.4, 3.3, 4 | discretuum, window, dynamics-at-BBN argument |
| Bousso, Harnik, Kribs, Perez | PRD 76, 043513 (2007), hep-th/0702115 | abstract | causal entropic principle |
| Freivogel, Kleban, Rodríguez Martínez, Susskind | JHEP 0603, 039 (2006), hep-th/0505232 | abstract | open-curvature landscape signature |
| Barrow, Shaw | PRL 106, 101302 (2011), arXiv:1007.3086 | full | eq. (3), Omega_k0 prediction |
| Barrow, Shaw | GRG 43, 2555 (2011), arXiv:1105.3105 | full | causal Lambda-variation; Lambda ~ t_U^-2 |
| Ahmed, Dodelson, Greene, Sorkin | PRD 69, 103523 (2004), astro-ph/0209274 | abstract | everpresent Lambda |
| Zwane, Afshordi, Sorkin | CQG 35 (2018), arXiv:1703.06265 | abstract | tests; vanishing mean |
| Das, Nasiri, Yazdi | JCAP 10 (2023) 047, arXiv:2304.03819; JCAP 10 (2024) 076, arXiv:2307.13743 | abstracts | typicality; SN/CMB tests |
| Anderson, Barbour, Foster, Ó Murchadha | CQG 20, 1571 (2003), gr-qc/0211022; gr-qc/9911071 | abstracts | scale-invariant gravity |
| Padmanabhan, Padmanabhan | IJMPD 22, 1342001 (2013), arXiv:1302.3226 | abstract | CosMIn pointer |
| Tseytlin | PRL 66, 545 (1991) | citation only (via ADDG ref. [4]) | lineage |
| Einstein 1917; Tryon, Nature 246, 396 (1973) | — | recall | historical rows |

Downloads are kept in the scratchpad (`.../scratchpad/others_dl/`). No paper text or PDF is in the repository.

## 13. Files

- `code/others_checks.py`: sections A-H; numpy/scipy/sympy; deterministic, about 5 s. Run with
  `python3 -I code/others_checks.py results/others_checks.json`.
- `results/others_checks.json`: every [C] number above, the 26 checks (all pass), versions, and scope.
