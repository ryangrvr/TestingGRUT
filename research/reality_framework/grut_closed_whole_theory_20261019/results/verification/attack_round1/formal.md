# Round 1, lens: equation and theorem checker

Packet: `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019`
Files read in full: FORMAL_CORE.md, WORKING_THEORY.md (all, section 2 line by line), CALCULATIONS.md,
code/closed_whole_demo.py, code/marginal_compatibility.py, OPEN_CLAIMS.md, SOURCES.md, README.md,
the handoff (NEXT_CHAT_HANDOFF.md, START_NEXT_CHAT.txt).
Primary texts fetched and read (pdftotext, in `dl/`): Fewster-Verch arXiv:1106.4785 (AHP 13, 1613);
Masanes-Galley-Muller arXiv:1811.11060 (Nat. Commun. 10, 1361); Witten arXiv:1803.04993 (RMP 90,
045003); Brunetti-Fredenhagen-Verch math-ph/0112041; abstracts of Pollock et al. 1512.00589 and
1801.09811, Cotler-Penington-Ranard 1702.06142, Linden-Popescu-Short-Winter 0812.2385.
Scripts: `checks.py` (output `checks_output.txt`), `nz_check2.py` (output `nz_check2_output.txt`),
both run with `python3 -I` against the packet's exact 8-qubit datum (re-implemented, not imported).

Nothing under /home/user/TestingGRUT was modified.

---

## Findings (most severe first)

### F1. THM 2.3(ii): the "iff" for the inhomogeneous term is false -- MAJOR
- File / location. FORMAL_CORE.md:107-108: "the inhomogeneous term vanishes iff Q rho(0) = 0, i.e. iff
  rho_Omega(0) = rho_A(0) (x) rho_ref". Repeated in WORKING_THEORY.md:103: "the inhomogeneous term
  vanishes only for product preparations".
- Problem. Only the "if" direction holds. I(t) = P L e^{QLt} Q rho(0) vanishes identically whenever
  Q rho(0) != 0 lies in a QL-invariant subspace annihilated by PL; the simplest case is a correlation
  chi with Tr_Abar chi = 0 and [H, chi] = 0. The packet's own datum has such an operator: the global
  Z-parity P = Z^{(x)8}, which it states is conserved.
- Deciding check. `nz_check2.py`: rho(0) = (1/2)(x)(1/128) + eps P with eps = 1/512 (the packet's own S4
  value; ||eps P||_1 = 0.5, the same size as its chi = Y0X1 case): min eig = 1/512 > 0;
  ||Q rho(0)||_F = 0.03125 != 0; ||L Q rho(0)||_F = 0.0 exactly, so e^{QLt} Q rho(0) = Q rho(0) for all t
  and I(t) = P L Q rho(0) = 0 for all t; max_t D(rho_A^(1)(t), rho_A^(2)(t)) = 3e-16 over [0, 40].
  Contrast in the same run: the packet's chi = Y0X1 gives ||Tr_Abar L(eps chi)||_F = 0.35355
  (reproduces CALCULATIONS S4). The WORKING_THEORY wording is wrong in the other direction too: a
  product preparation rho_A (x) sigma with sigma != rho_ref has Q rho(0) = rho_A (x) (sigma - rho_ref)
  != 0 and generically I != 0.
- Fix. FORMAL_CORE: "(ii) I == 0 if Q rho(0) = 0, i.e. if rho_Omega(0) = rho_A(0) (x) rho_ref; the
  converse fails (e.g. rho_Omega(0) = rho_A (x) rho_ref + eps chi with Tr_Abar chi = 0 and [H, chi] = 0,
  such as the conserved Z-parity of calc 1, has Q rho(0) != 0 and I == 0)". WORKING_THEORY: "vanishes
  when rho_Omega(0) = rho_A(0) (x) rho_ref".

### F2. COND 5.3 misstates the hypotheses of Fewster-Verch's natural-state theorem -- MAJOR
- File / location. FORMAL_CORE.md:216-221: "A covariant rule assigning a preferred state to every
  spacetime is excluded for locally covariant theories satisfying dynamical locality with non-trivial
  relative Cauchy evolution (Fewster-Verch, AHP 13, 1613 (2012), extending Brunetti-Fredenhagen-Verch
  2003)."
- Problem. The theorem has three further hypotheses the packet omits, and its conclusion is weaker
  than "excluded". FV Definition 6.12 / Theorem 6.13 (`dl/fv2012.txt` lines 2238-2247): A dynamically
  local in LCT, with a natural state (natural under ALL Loc morphisms); IF there is a spacetime M with
  NONCOMPACT Cauchy surfaces such that omega_M induces a FAITHFUL GNS representation with the
  REEH-SCHLIEDER property, then the relative Cauchy evolution is trivial in M; the theory is trivial
  only if in addition it obeys extended locality. The noncompact-Cauchy-surface hypothesis is used in
  the proof ("here we use the noncompactness of the Cauchy surfaces", line 2262) and is salient for a
  program whose campaigns concern closed universes: as proved, the theorem does not address a category
  restricted to spatially compact spacetimes. The intro (lines 231-236) states the same three extra
  conditions. "Extending BFV 2003" overstates: FV say it "brings to sharper form an argument sketched
  in BFV and [32] for the free scalar field" (line 2230); BFV's own text (`dl/bfv2003.txt` 836-875)
  only sketches a free-field heuristic and attributes the non-existence result to their ref. [24].
  SOURCES.md line "Primary texts actually inspected ... [[LITERATURE]]" is unfilled, consistent with
  the handoff's warning (NEXT_CHAT_HANDOFF.md: "Only search/abstract metadata was checked ... inspect
  the relevant primary text before applying its theorems") not having been acted on.
- Deciding check. `dl/fv2012.txt` lines 2226-2280; abstract lines 36-38 ("together with typical
  assumptions").
- Fix. Restate: "FV 2012 Thm 6.13: a dynamically local theory with a natural state (covariant under
  every Loc embedding) that, in some spacetime with noncompact Cauchy surfaces, is faithful with the
  Reeh-Schlieder property, has trivial relative Cauchy evolution there; with extended locality the
  theory is trivial. It excludes embedding-covariant state rules under these hypotheses; it says
  nothing about rules covariant only under isometries, about finite models, or (as proved) about
  categories of spatially compact spacetimes." Fill [[LITERATURE]] with the texts actually read.

### F3. COND 5.3 and section 7 apply the natural-state obstruction where it does not apply; the real finite-class obstruction is unstated -- MAJOR
- File / location. FORMAL_CORE.md section 7 table row "state selection: rho_Omega(0) = F(H, access) ...
  must survive COND 5.3 and the marginal problem"; WORKING_THEORY.md section 7 "a positive, normalizable
  state functional rho_Omega(0) = F(H, access) that survives the natural-state obstruction (FORMAL_CORE
  COND 5.3) and the marginal problem".
- Problem. (a) FV's theorem is about functorial covariance over Loc in locally covariant QFT. In the
  finite class the packet works in, rules such as rho_Omega(0) = e^{-beta H}/Z or the ground projector
  of H exist and violate nothing in COND 5.3. (b) The genuine finite-class obstruction is different and
  provable: a rule F depending on H alone and equivariant under unitary conjugation,
  F(U H U^+) = U F(H) U^+, satisfies [F(H), H] = 0 (take U = e^{iHs}), so the selected global state is
  stationary and rho_A(t) = rho_A(0) for every subsystem -- the finite analogue of "trivial relative
  Cauchy evolution". A rule must therefore use the access structure (and so break unitary
  equivariance) to produce any local dynamics. (c) "Survive the marginal problem" is a category error
  for a global rule: any global state has consistent marginals automatically (THM 2.6-2.7); the
  marginal problem obstructs bottom-up assembly of local families, not a top-down F.
- Deciding check. Algebraic argument in (b); `checks.py` section E: Gibbs(beta=1) and ground projector
  of the packet's H give ||[F, H]||_max ~ 1e-15 and max_t D(rho_A(t), rho_A(0)) ~ 1e-15.
- Fix. Replace the state-selection cost with: "in the finite class, a unitarily equivariant rule F(H)
  is stationary (no local dynamics), so a rule with content must depend on the access structure; in
  locally covariant QFT an embedding-covariant rule is excluded under the hypotheses of FV Thm 6.13;
  bottom-up local families must additionally pass the marginal problem."

### F4. Section 7 "exhaustive by DEF 1.7" is false on its own terms -- MAJOR
- File / location. FORMAL_CORE.md:245-246 "The possibilities are exhaustive by DEF 1.7"; table rows:
  state selection, dynamical constraint, access constraint, modified global law, clock constraint,
  measurement constraint. WORKING_THEORY.md:87 "lists exhaustively".
- Problem. DEF 1.7 has six slots: (1) algebra or H_Omega, (2) H, (3) rho_Omega(0), (4) {M_A},
  (5) instruments and records, (6) clock. The table has rows for slots 2, 3, 4, 6; "modified global
  law" and "measurement constraint" are not slots of Theta (they modify DEF 1.1 and DEF 1.8); slots 1
  and 5 have no row. Missing law types: (i) kinematic constraints on the algebra (finite dimension /
  holographic bound, type, superselection and symmetry structure, spacetime dimension);
  (ii) constraints on physical instruments and records (conservation-law limits of Wigner-Araki-Yanase
  type, einselected record bases); (iii) measure or typicality postulates over Theta, which force
  statements with high probability rather than by selection (Popescu-Short-Winter, Nat. Phys. 2, 754
  (2006); Goldstein-Lebowitz-Tumulka-Zanghi, PRL 96, 050403 (2006)); (iv) relational constraints
  linking slots, e.g. locality of H selecting the access structure. Item (iv) is a known partial result
  directly relevant to the access row's cost "must not presuppose the factorization it derives" and to
  DEF 1.2's "no rule below derives it": Cotler-Penington-Ranard, "Locality from the spectrum", Commun.
  Math. Phys. 368, 1267 (2019), arXiv:1702.06142, argue that for generic H the spectrum determines a
  local tensor product structure when one exists, up to stated dualities (their claims are arguments
  with stated hypotheses, not full proofs; the packet should cite and scope them rather than omit them).
- Deciding check. Count of DEF 1.7 slots (FORMAL_CORE.md:69-72) against table rows
  (FORMAL_CORE.md:248-255); abstract of arXiv:1702.06142 (fetched).
- Fix. "The slots of Theta are exhaustive by construction; a law may constrain one slot, relate two
  slots, or impose a measure over Theta." Add rows for the algebra slot, the instrument/record slot,
  typicality/measure postulates, and cross-slot relational constraints (with CPR 2019 cited and scoped
  in the access row).

### F5. Section 3 "Exactly the theorems of section 2" omits what closure plus a generic condition does force locally -- MAJOR
- File / location. FORMAL_CORE.md section 3 "Exactly the theorems of section 2"; section 4 row
  "universal or monotone relaxation NOT FORCED"; WORKING_THEORY.md section 0 "Closure constrains
  nothing local".
- Problem. Closure with one generic, checkable hypothesis on H (non-degenerate energy gaps) forces
  time-averaged equilibration of every subsystem for every initial state:
  <D(rho_A(t), omega_A)>_t <= 1/2 sqrt(d_A^2 / d_eff), omega the dephased initial state,
  d_eff = 1/Tr omega^2 (Linden-Popescu-Short-Winter, PRE 79, 061103 (2009); Reimann, PRL 101, 190403
  (2008); Short, NJP 13, 053009 (2011) for mixed initial states). This is a local consequence of
  closure not in section 2, standard, and conditional in exactly the packet's COND sense. "Monotone
  relaxation not forced" is correct; "nothing local" and "exactly section 2" are not the right boundary.
- Deciding check. `checks.py` section D on the packet's S1 datum: non-degenerate gaps verified
  (0 coincident ordered gaps among 256*255 at 1e-9; min gap 5.9e-4); d_eff(omega) = 4.613; LPSW bound
  0.4656; grid time-average of D(rho_A(t), omega_A) over [0, 40] = 0.3717 (max 0.509), consistent with
  the theorem. The bound is loose here because the supplied preparation has small d_eff, which is
  itself the right lesson: the size of forced equilibration is set by the supplied preparation.
- Fix. Add "COND 2.10 (equilibration on average)" with the hypothesis and bound, and add canonical
  typicality as a COND statement; change section 3 to "the theorems of section 2, plus the conditional
  equilibration and typicality theorems"; change the section 4 row to "pointwise or monotone
  relaxation: NOT FORCED; time-averaged equilibration: FORCED under non-degenerate gaps, with magnitude
  set by d_eff".

### F6. DEF 1.8: "operational closure" is not Masanes-Galley-Muller's assumption -- MINOR
- File / location. FORMAL_CORE.md:81-82: "Masanes-Galley-Muller, Nat. Commun. 10, 1361 (2019), whose
  'operational closure' is a composition property of measurements".
- Problem. The phrase does not occur in the paper. MGM derive the measurement postulates from the three
  "unitary QM" postulates (states, transformations, composite systems), the "possibility of state
  estimation" assumption, and structural requirements on outcome-probability-function sets: closed
  under mixtures, under unitaries, and "closedness under system composition" (Definition 9;
  associativity of the star-product, eq. 14).
- Deciding check. `dl/mgm2019.txt` lines 195-212 and 640-660 (Definition 9), line 231 (state
  estimation); grep for "closure" finds only "closure of F_d" (completeness) and "closedness under
  system composition"; the arXiv abstract page has no "closure".
- Fix. "... whose 'possibility of state estimation' and 'closedness under system composition' (Def. 9)
  are operational assumptions about measurements, not the dynamical closure of DEF 1.1".

### F7. C5 / THM 2.7: "the obstruction is odd frustration, not cyclicity" overgeneralizes -- MINOR
- File / location. CALCULATIONS.md:183 and code note in marginal_compatibility.py ("The obstruction is
  an odd number of anticorrelated edges around the cycle, not cyclicity as such"); FORMAL_CORE.md:150
  "(odd frustration)".
- Problem. On a cycle, pairwise-consistent families fail to extend exactly when an n-cycle
  (Bell/CHSH-type) inequality is violated (Araujo-Quintino-Budroni-Terra Cunha-Cabello, PRA 88, 022118
  (2013)); perfect odd anticorrelation is the deterministic extreme, not the characterization.
- Deciding check. `checks.py` section C: 4-cycle with <AB> = <BC> = <CD> = 1/sqrt2, <DA> = -1/sqrt2,
  uniform singles; every pair marginal is a valid distribution, overlaps agree, no edge is perfectly
  (anti)correlated; LP for a joint distribution: infeasible (CHSH = 2.828).
- Fix. "a cycle permits failure (violation of an n-cycle inequality, e.g. an odd number of perfect
  anticorrelations) but does not force it".

### F8. THM 2.3: TCL "equivalently" needs an invertibility hypothesis; P needs Tr rho_ref = 1 -- MINOR
- File / location. FORMAL_CORE.md:103-104 "or equivalently rho_A-dot = K_TCL(t) rho_A + I_TCL(t)".
- Problem. The TCL generator exists only where the Q-resolved propagator (1 - Sigma(t)) -- equivalently,
  for product preparations, Phi_t -- is invertible; this can fail at finite times (Breuer-Petruccione,
  The Theory of Open Quantum Systems, section 9.2; Shibata-Takahashi-Hashitsume 1977). The packet's S3
  uses Phi-dot_t Phi_t^{-1} and reports condition numbers, i.e. the hypothesis is used tacitly.
  P^2 = P requires rho_ref to be a unit-trace operator on H_Abar, not stated.
- Deciding check. Definition of K_TCL = Phi-dot Phi^{-1} in code/closed_whole_demo.py (`tcl_generator`);
  standard derivation.
- Fix. "equivalently, wherever 1 - Sigma(t) is invertible, ..."; "rho_ref a density operator on H_Abar".

### F9. COND 5.1: "are type III_1" cited to Witten -- MINOR
- File / location. FORMAL_CORE.md:201-203 "(Witten, RMP 90, 045003 (2018), sections 6-7)".
- Problem. Witten section 6.5 says the local algebras "are of Type III ... In fact, they are believed
  to be of Type III_1" and gives a heuristic (`dl/witten.txt` 3027-3029, 3095). The rigorous sources are
  Fredenhagen, Commun. Math. Phys. 97, 79 (1985) and Buchholz-D'Antoni-Fredenhagen, Commun. Math. Phys.
  111, 123 (1987), under a scaling-limit hypothesis. Section numbers are otherwise right: no
  factorization is section 6.1 (lines 2796-2798), split/factorized states section 7.
- Fix. Cite BDF 1987 for III_1 with its hypothesis; keep Witten as the review.

### F10. COND 5.2: KMS sign convention unstated -- MINOR
- File / location. FORMAL_CORE.md:208 "KMS at beta = 1 with respect to its own modular flow".
- Problem. With sigma_t(A) = Delta^{it} A Delta^{-it}, a faithful normal state is KMS at beta = -1
  (Takesaki; Bratteli-Robinson II, Thm 5.3.10), i.e. at beta = 1 w.r.t. sigma_{-t}; Connes-Rovelli
  choose the opposite sign so that beta = 1.
- Fix. State the convention.

### F11. THM 2.5 COND: the influence-functional formula is the factorized-preparation form -- MINOR
- File / location. FORMAL_CORE.md:124-131: F[x, x'] = Tr_Abar[U_Abar[x] rho_Abar U_Abar[x']^+] with the
  COND "prepared in a Gaussian state (factorized or Gaussian-correlated)".
- Problem. For correlated preparations the functional carries a preparation function and is not of the
  quoted form (Grabert-Schramm-Ingold, Phys. Rep. 168, 115 (1988); Hakim-Ambegaokar, PRA 32, 423 (1985)).
- Fix. Restrict the formula to factorized preparations or cite GSI for the correlated case.

### F12. DEF 1.1 vs DEF 1.6: clock wording contradictory -- MINOR
- File / location. FORMAL_CORE.md:29 "no reservoir, drive, measurement device or clock outside Omega";
  FORMAL_CORE.md:65 "the clock that reads it is not part of Omega".
- Problem. If no clock is outside Omega and none inside, there is no clock; t is a bare parameter.
- Fix. "t is a parameter of U(t) read by no subsystem of the model".

### F13. Section 4 rows that are definitional are labelled as countermodels -- MINOR
- File / location. FORMAL_CORE.md section 4 rows "a pure global state" (countermodel: "DEF 1.1 admits
  mixed rho_Omega"), "the Born rule" ("DEF 1.8"), "the observed constants" ("elements of Theta"); row
  "persistence of records: stable iff [H, Z_R] = 0".
- Problem. Section 0 defines an answer as a derivation or a countermodel in class M; these rows settle
  the question by how DEF 1.1/1.8/1.7 were written. For purity the substantive argument (no exterior to
  be correlated with, so mixedness is epistemic) is excluded by definition, not refuted. "Stable iff
  [H, Z_R] = 0" holds for all states and times; for the particular recorded state stability can hold
  with [H, Z_R] != 0 (the packet's kappa X_R case keeps optimal readout 1.000).
- Fix. Label such rows "BY DEFINITION (not a countermodel)"; write "stable for all states iff".

### F14. THM 2.7 and 2.4 citation/wording nits -- MINOR
- THM 2.7: the Bell-chain argument (pure marginal => product) does not use Coffman-Kundu-Wootters; cite
  CKW only for the quantitative tangle inequality. Liu 2006's QMA-completeness is under Turing
  reductions; Karp-reduction completeness is Broadbent-Grilo, SIAM J. Comput. 51, 1400 (2022).
- THM 2.4: "Phi_{t2} != V_{t2,t1} Phi_{t1} with V completely positive" should read "there need not exist
  a completely positive V with Phi_{t2} = V Phi_{t1}"; BLP 2009 defines the information-backflow
  measure, RHP 2010 the divisibility criterion.
- THM 2.3 / S4: the "t = 0 value of the inhomogeneity" quoted (0.3536) is ||Tr_Abar L(eps chi)||_F,
  i.e. the A-factor of I(0) without the (x) rho_ref factor; say so.

---

## Confirmed (checked and found correct as stated)
- THM 2.1, 2.2 (cyclicity of the partial trace), 2.6: correct.
- THM 2.3 NZ equation: signs, kernel K(tau) = P L e^{QL tau} Q L P and I(t) = P L e^{QLt} Q rho(0) are the
  standard Schrodinger-picture forms for time-independent H; references Nakajima PTP 20, 948 (1958),
  Zwanzig JCP 33, 1338 (1960), Shibata-Takahashi-Hashitsume JSP 17, 171 (1977) correct.
- S4 hand computation: [g X0X1, Y0X1] = 2ig Z0, -i Tr_Abar[H, eps chi] = (1/4) Z0, Frobenius
  sqrt2/4 = 0.35355 (reproduced to 1e-16).
- THM 2.7 Bell-chain monogamy argument: rho_BC = (1/2) (x) rho_C, fidelity 1/4; "chain" is the right
  contrast (the chain hypergraph is acyclic, so Vorob'ev's theorem makes every consistent classical
  family extend there); Vorob'ev reference (Theory Probab. Appl. 7, 147 (1962)) and the acyclicity
  characterization correct.
- Calc 2 numbers re-derived independently (cvxpy/CLARABEL, own code): max eig of
  Phi+_AB (x) 1 + 1 (x) Phi+_BC = 1.5; isotropic chain p* = 2/3; isotropic triangle p* = 1/3; Werner
  chain 2/3 and Werner triangle 1/3 (all within 1e-8).
- Pollock et al.: PRA 97, 012127 (2018) "Non-Markovian quantum processes: complete framework and
  efficient characterisation"; PRL 120, 040405 (2018) "Operational Markov condition for quantum
  processes" -- titles, journals, content as described (process tensor; operational Markov condition).
- Fewster-Verch AHP 13, 1613 (2012) and MGM Nat. Commun. 10, 1361 (2019) journal references correct.
- Witten RMP 90, 045003 (2018) reference and sections 6-7 as the relevant sections (no factorization
  6.1, type III 6.4-6.5, factorized states 7).
- Bocchieri-Loinger Phys. Rev. 107, 337 (1957) "Quantum recurrence theorem": correct; THM 2.9's
  almost-periodicity statement correct in finite dimension; Lieb-Ruskai 1973, Araki-Lieb inequality
  correct.
- Pechukas PRL 73, 1060 (1994); Stelmachovic-Buzek PRA 64, 062106 (2001); Rivas-Huelga-Plenio PRL 105,
  050403 (2010); Breuer-Laine-Piilo PRL 103, 210401 (2009); Caldeira-Leggett Physica A 121, 587 (1983)
  kernels eta ~ int J sin, nu ~ int J coth cos; Page-Wootters PRD 27, 2885 (1983); Connes-Rovelli CQG
  11, 2899 (1994); Gleason 1957; Busch PRL 91, 120403 (2003): all correct.
- COND 5.1's "GRUT III's P1 is exactly this identification": matches GRUT_III/01_PRINCIPLES.md:23-30
  (bath = commutant A(O)', = A(O') where Haag duality holds; type III_1; no reduced density matrix).
- Equilibration hypothesis holds in the packet's model (non-degenerate gaps), so the omitted theorem
  (F5) applies to calc 1 as supplied.
