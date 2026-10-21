# Track `grut_native`: the minimal GRUT-native vacuum principle, attempted

Packet `grut_vacuum_principle_20261021`, track `grut_native`, 21 October 2026. Written under `../../BRIEF.md`.

## P. Pre-registration (written 2026-10-21T20:00Z, before any code was written or run)

This section was written and hashed (`results/prereg.json`) before any numerical or symbolic
calculation of this track. The route assessment behind it was done on paper (R1-R3 below), and no
number below was computed yet. Later sections may disagree with it; if they do, this section is not
edited, and the disagreement is reported in the results section.

**Most promising route, chosen in advance: R1 + R1b, "clock closure with a counting clock".**

- R1 (clock closure). GRUT's FORMAL_CORE DEF 1.6 leaves the time parameter external and supplied.
  Removing that external clock, and choosing as internal clock the variable conjugate to the energy
  origin, gives the Henneaux-Teitelboim / Unruh unimodular structure: the total effective vacuum
  term Lambda_eff becomes the eigenvalue label of the whole's state with respect to the 4-volume
  clock T, and every constant matter vacuum shift becomes a phase exp(-i c T).
- R1b (counting clock). GRUT III's own map (GRUT_III/13_SYNTHESIS.md §7, layer 1) defers the
  microstructure to "causal order plus counting" (number ~ 4-volume, statistically). With Poisson
  counting, Delta V ~ sqrt(V l^4), and the conjugate relation Delta Lambda_eff Delta V >~ 1 gives
  rho_Lambda(t) = xi(t) alpha M^2 V(t)^(-1/2), with xi a zero-mean, unit-variance stochastic process
  and alpha an O(1) amplitude that absorbs the discreteness scale l relative to 1/M.
- Credit, fixed in advance: this combination is Sorkin's everpresent Lambda (Sorkin; Ahmed, Dodelson,
  Greene, Sorkin 2004) with Henneaux-Teitelboim and Unruh for the clock. It is a rederivation, not a
  GRUT-native principle. Its zero mean is a postulate, not a derivation. No assumption that the cited
  literature lacks was identified at pre-registration time.

**Independently testable consequence, written before any calculation.**

1. The dark-energy density is not a constant. It fluctuates with r.m.s. amplitude of order
   alpha M^2 V(z)^(-1/2), a roughly fixed fraction of the critical density at every epoch, with sign
   changes allowed and a correlation length of about one Hubble time (Delta ln V of order 1).
2. Its mean is zero. In any epoch, a negative value is about as likely as a positive one.
3. The test is held out only in this form. alpha is calibrated once on exposed data (a one-parameter
   fit), then frozen. The held-out statistic is the bin-to-bin structure of rho_DE(z) reconstructed in
   independent redshift bins over 0 < z < 2.5 from data not released when alpha was frozen: DESI
   releases after DR2, Euclid cosmology releases, Rubin LSST supernovae, Roman supernovae. If those
   data show rho_DE(z) constant across bins at a precision much smaller than the frozen
   alpha rho_crit(z), with a probability below 1% over the model's own realizations, the route fails.
4. Already exposed, therefore not held out: Planck 2018 CMB, ACT DR6, SPT-3G, DESI DR1 and DR2 BAO
   (including its w0 > -1, wa < 0 preference), Pantheon+, Union3 and DES-Y5 supernovae, BBN
   abundances, distance-ladder H0, and the current Omega_K constraints. Any agreement with these is
   retrospective. Rho_obs itself is exposed, so reproducing its order of magnitude is not a test.
5. Route-level statement independent of the stochastic part: under R1 alone, Lambda_eff is a
   whole-state label decoupled from constant matter vacuum terms at every loop order, and it has no
   classical consequence beyond GR + Lambda. R1 alone predicts nothing new.

**Checks fixed in advance (consistency checks on exposed numbers, not tests).**

- K1. The exponent p in rho/M^4 = (V M^4)^(-p) that reproduces rho_obs, with V the 4-volume of our
  causal past. R1b's scale claim requires p within [0.45, 0.55] for every convention tried (reduced
  or non-reduced Planck length; causal-past or Hubble 4-volume). Outside that window, R1b fails K1.
- K2. Whether the Gibbs (KMS) state of the whole with respect to the internal 4-volume clock is
  normalizable, and what it selects.
- K3. Whether a fixed finite total 4-volume (the minimal one-variable global constraint) is compatible
  with a constant positive Lambda_eff and the present acceleration.
- K4. The sign that an "energy budget of the whole is zero" rule gives the vacuum term.
- K5. Which state rules f(H)/Tr f(H) are invariant under H -> H + cI.
- K6. A Page-Wootters toy with an internal clock: is the matter shift a pure clock phase, and is the
  Lambda_eff distribution shift-invariant?

## 0. Result

Labels: [D] derived here; [C] computed by `code/grut_native_checks.py` (record
`results/grut_native_checks.json`, 29/29 checks pass, seed 20261021, about 4 s); [Q] quoted or
paraphrased from a source at the read level in §9; [E] order-of-magnitude estimate; [A] argued, not a
theorem. Nothing here is a new theory, an empirical validation or a confirmatory test. Every number
built from H0, Omega_Lambda, Omega_m, BBN, CMB, BAO or SN data is **retrospective**.

**No GRUT-native principle was found. The result of this track is the obstruction, stated precisely
in §5.** In short:

1. **R1, clock closure, escapes Codex's obstruction for the gravitating term, at the cost of an added
   clock postulate [D, C: A1-A4, B1-B5].** Closing GRUT's clock slot (no external t, FORMAL_CORE DEF
   1.6) with the clock conjugate to the energy origin gives the Henneaux-Teitelboim (HT) / Unruh
   structure. Every constant vacuum term then enters as h -> h + cI in the unimodular-time generator, a
   pure clock phase. The gravitating Lambda_eff becomes the spectral label of a c-independent operator
   h0 that the whole's state occupies: a quantum number of the state, not a law parameter. Q1 and Q3
   pass for constant terms; Q2 and Q4 fail. This is (ii), a rederivation (HT 1989; Unruh 1989;
   Unruh-Wald 1989; Sorkin), and it agrees with the `unimodular` track's §3, which I do not repeat.
2. **R1a, state principles, fix nothing, fix the wrong value, or move the value into a supplied
   number [D, C: C1-C4, D1-D4].** The only equivariant state rules f(H)/Tr f(H) invariant under
   H -> H + cI are Gibbs states. In R1 the Gibbs (KMS) state with respect to the 4-volume clock is not
   normalizable with the Gibbons-Hawking density of states: all its weight sits at the cutoff as
   Lambda_eff -> 0+, which is Baum-Hawking. With a flat density, it gives <Lambda_eff> = 1/|beta|, so
   the value moves into the clock temperature beta. A self-consistent KMS condition built from M and H
   alone is either Planckian or an identity. Rules anchored on the zero of h0 give Lambda_eff = 0.
   This is (ii) or (iii).
3. **R1b, fluctuation, needs discreteness that GRUT does not supply [A, C: E2-E5].** GRUT's finite
   Hilbert space is a regulator. Its CP-12 graph runs in external time. GRUT III defers its
   microstructure layer to causal-set theory ("GRUT has no dynamics of its own" there). Adopting that
   deferral imports Sorkin's Poisson kinematics and gives everpresent Lambda (credited; scored in the
   `unimodular` track). Pre-registered check K1 passes (p = 0.499-0.506) but is weak: it only confirms
   the coincidence form rho ~ M^2 H^2.
4. **R2, global self-consistency, is a rederivation, and its minimal form is not Kaloper-Padilla
   (KP) [D, C: F1-F3, G1-G2, H1-H3].**
   - With one global variable, Lambda decouples from rho_vac only if sigma is linear. That is a fixed
     total 4-volume V0, i.e. HT with Dirichlet data on T. V0 is then a supplied substitute for
     Lambda_eff.
   - A constant Lambda_eff > 0 with the present acceleration gives infinite 4-volume, so this minimal
     form is retrospectively incompatible with the expansion unless a collapse sector is added.
   - The "bath = rest" reading (measure the vacuum relative to the whole) has two inequivalent covariant
     realizations. They agree on pure vacuum and differ on matter: KP's (1/4)<T>, which needs a second
     global variable, and Tseytlin's S/V, which gives <tau/2 - L_tau> with an on-shell-Lagrangian
     ambiguity. GRUT selects neither.
5. **R3, closure of the energy budget, is circular or the wrong sign [D, C: I1-I4].**
   - With gravity included, the whole's energy is the Hamiltonian constraint, which vanishes for every
     Lambda (iii).
   - With gravity excluded, the rule gives rho_Lambda = -<rho_matter> <= 0, and it requires a preferred
     frame. On a k = +1 slice the non-gravitational energy can never vanish.
   - For flat dust plus radiation with Lambda < 0, the history-averaged condition has no solution at
     any scale: <rho>/|rho_Lambda| = 2.00 to 2.92, never 1.
6. **Most promising route, pre-registered (§P): R1 + R1b. It is credited to Sorkin et al. and is
   already strained by exposed data** (the `unimodular` track's Monte Carlo puts P(our history | alpha)
   <= 1e-4). Its consequence, a stochastic sign-changing dark-energy density tracking rho_crit, can be
   held out only if alpha, the start epoch and the update recipe are frozen before data released after
   DESI DR2.

## 1. Setup

**Conventions.** Signature (-,+,+,+), hbar = c = 1, M = (8 pi G)^(-1/2) = 2.435e18 GeV (reduced),
M_P = sqrt(8 pi) M. Target: the total renormalized rho_vac/M^4. Observed, retrospective:
rho_obs/M^4 = 3 Omega_Lambda (H0/M)^2 = 7.15e-121 (Planck 2018: h = 0.6736, Omega_m = 0.3153,
Omega_r = 9.2e-5) [C: E0]. The siblings quote 7.16e-121 for h = 0.674 and Omega_Lambda = 0.685.

**GRUT datum.** Theta = (algebra, dynamics, preparation, access structure, instruments, clock)
(FORMAL_CORE DEF 1.7). DEF 1.1 and 1.6 make t external and supplied; COND 5.4 notes that with gravity
there is no preferred t and no conserved total energy. "Every subsystem's bath is the rest" is DEF 1.3
(bath = commutant), a definition with no further property. GRUT III treats Lambda as a datum, recovered
as a Jacobson integration constant (01_PRINCIPLES.md:10, :52). Its pre-registered N2 ("does any
horizon principle fix a renormalisation constant, Lambda in particular, without a fit?") is open with
no candidate (13_SYNTHESIS.md:366).

**What counts as GRUT-native (BRIEF).** A route is GRUT-native only if it uses an assumption that the
cited literature does not. Each route below records which slot of Theta it changes (FORMAL_CORE §7).

**Not repeated as new.** The following are already established and are only cited here:
- The fixed-background identity and its local-covariant form (`obstruction` §§1, 3).
- The minimal escape classes G/B/S/X/R/P (`obstruction` §6.2).
- The HT canonical identity and the Lambda-translation lemma (`unimodular` §3).
- Everpresent-Lambda numbers and Monte Carlo (`unimodular` §4).
- Sequestering versions and the linear-trigger failure (`sequestering`).
- Barrow-Shaw, Euclidean measures and the landscape (`others`).
- From the GRUT record: "nothing outside" does not imply finite 4-volume (FORMAL_CORE §4); no state
  selection rule survived (C1, C2, s09); the withdrawn all-sequestering k > 0 claim (closure fork).

## 2. Route R1: "no external clock"

### 2.1 From the GRUT clock slot to the HT action [D/A]

On a fixed background, the energy origin c is conjugate to the external t: e^{-i(H+c)t} =
e^{-ict} e^{-iHt}, and c is invisible because no subsystem reads t (DEF 1.6). The clock-closure premise
says that a closed whole has no external clock, so t must be an observable of the whole. If the clock
is internal, the generator of its translations is a dynamical variable. The shift H -> H + c then
becomes a shift of the clock's conjugate momentum, and the origin becomes physical relative to it. In
Page-Wootters form, (p_T + H_rest) Psi = 0 and H_rest -> H_rest + c is p_T -> p_T - c.

With a dynamical metric, the constant couples through sqrt(-g) (`obstruction` §1.5). The variable
conjugate to it is therefore the 4-volume, d_mu calT^mu = sqrt(-g). The minimal action is HT's:

    S = ∫ d^4x [ sqrt(-g) (M^2 R/2 - Lambda) + Lambda d_mu calT^mu ] + S_m[g, psi].        (2.1)

**Added assumptions, explicitly:**
- (R1-i) The time parameter of the law is an observable of the whole (clock closure). This removes
  DEF 1.6's supplied external t.
- (R1-ii) The clock is the one conjugate to the energy origin, i.e. 4-volume time. Clock closure alone
  does not force this (`unimodular` §5 step 3: other clocks give inequivalent theories; Kuchar 1991).
- (R1-iii) The global fields calT^mu and Lambda(x) are in the algebra slot.

GRUT's own COND 2.20(ii) gives a further GRUT-internal reason for (R1-ii): no internal clock runs
monotonically if its generator is bounded below. In (2.1) the generator is unbounded, because
Lambda_eff ranges over the reals. This is the Unruh-Wald 1989 motivation [Q abstract via `unimodular`
§5], and it is credited to them.

### 2.2 FRW minisuperspace and the generator [D, C: A1-A4]

With unit comoving volume:

    L = -3M^2 a a'^2/N + 3M^2 k N a - N a^3 (rho_m(a) + rho_vac + Lambda) + Lambda T'.

- **Equations [C].** Varying T gives Lambda' = 0 (A2). Varying Lambda gives T' = N a^3 (A3). The lapse
  equation depends on rho_vac and Lambda only through Lambda_eff = Lambda + rho_vac (A4).
- **The compensated shift is a boundary term [C: A1].** (rho_vac, Lambda) -> (rho_vac + c, Lambda - c)
  changes L by exactly -c T'. This agrees with the `unimodular` track's B1.
- **Unimodular-time generator [D].** With p_T = Lambda and lapse N = a^-3, the constraint is
  p_T + h0 + rho_vac ≈ 0, where h0 = (gravitational and matter Hamiltonian)/a^3 does not contain
  rho_vac. Quantum mechanically this reads

      i dPsi/dT = (h0 + rho_vac) Psi,   Psi(T) = ∫ dLambda_eff phi(Lambda_eff) e^{-i(rho_vac - Lambda_eff)T} psi_{Lambda_eff},
      (h0 + Lambda_eff) psi_{Lambda_eff} = 0.                                                (2.2)

- **The energy origin is now a quantum number of the state [D].** rho_vac shifts the generator by cI,
  which is Codex's identity with t replaced by the internal clock T. The gravitating Lambda_eff is the
  spectral label of the c-independent operator -h0. It is not a coupling in the law. It is the support
  of phi, which is the state.

### 2.3 Toy check of (2.2), the access structure and which rules are invariant [C: B1-B5]

Setup: h0 is a random 8 x 8 Hermitian matrix, phi a random state, and c = 7310, a stand-in for a heavy
threshold.
- **B1-B2.** Conditional (Page-Wootters) expectations of a random observable agree to 8e-11 for
  T in [0, 30]. The shifted state equals e^{-icT} times the unshifted one to 3e-11.
- **B3, access structure.** No local subsystem has T in its algebra. HT's T(Sigma) = ∫_Sigma calT^0
  integrates over a whole slice. The state seen through the access structure is therefore the
  T-average. That average equals the Lambda_eff mixture sum |phi_k|^2 |k><k| to 1.2e-4 (finite
  averaging window) and is c-independent. For local observers Lambda_eff is a classical,
  superselected random variable with weights |phi|^2 [D/A].
- **B4-B5, which state rules are invariant.**
  - A rule stated on the bare multiplier label ("select the sector with Lambda = Lambda_eff - rho_vac
    nearest 0") moves Lambda_eff with c: 0.152, 1.233, -1.506 and 5.494 for rho_vac = 0, 1, -2 and 7310.
  - A rule anchored on h0's own zero ("Lambda_eff nearest 0") is invariant: 0.152 in every case.
  - This is the computed form of the Padilla-Saltas versus Carballo-Rubio dispute (`obstruction`
    §4.4): retuning is needed exactly when the state is specified by the bare label.

### 2.4 Why GRUT's closure theorems cannot select Lambda_eff [A, with D components]

- **Universal quantification.** Every FORMAL_CORE statement in §§2-5 has the form "for every Theta
  in class M (with hypothesis h), property P holds". After B3, the clock-closed whole decomposes for
  local observers into Lambda_eff sectors, and each sector is an ordinary Theta with definite
  Lambda_eff. Each theorem therefore holds sector by sector. A theorem that is universally quantified
  over a parameter cannot select a value of that parameter.
- **Hypothesis-carrying statements.** The COND statements hold or fail by their hypotheses, which do
  not depend on Lambda_eff:
  - non-degenerate gaps: generic;
  - KMS of the static patch at T_GH = H/2 pi: true for every Lambda_eff > 0 (COND 5.2, imported);
  - relaxation for an infinite complement: no Lambda-dependence in its premises.
- **Kinematic statements.** Statements about the T-clock kinematics are covered by the `unimodular`
  track's Lambda-translation lemma: they can fix a spread, never a mean.

### 2.5 Gate answers for R1

- **Parameters and inputs.** M; the HT fields; the state phi(Lambda_eff); the operator ordering of h0;
  the boundary condition at a = 0; the matter content. No cosmological input is needed to state the
  route. Reading Lambda_eff off (H0, Omega_m) is a retrospective measurement.
- **Q1 [D].** The extra structure is the internal 4-volume clock: the HT global fields, which make
  sqrt(-g)'s constant mode conjugate to a dynamical variable.
- **Q2 [D].** No. Lambda_eff is a free state label of either sign and any scale.
- **Q3 [D/Q].**
  - Every constant term c ∫ sqrt(-g), of any loop order or origin in the effective action, is a clock
    phase (B2). No retuning is needed, provided the state is specified by its h0 support (B5) and not
    by the bare label (B4).
  - Whether the full quantum HT theory beyond minisuperspace keeps this is disputed (`obstruction`
    §4.4) [Q].
- **Q4 [A].** None. In each Lambda_eff sector local physics is GR + Lambda. Stated in FORMAL_CORE COND
  2.18 terms: no local observable has a different response, state or higher cumulant. This agrees with
  Padilla-Saltas's "classical unimodular gravity = GR" [Q via `unimodular`].
- **Novelty.** (ii). HT 1989, Unruh 1989, Unruh-Wald 1989, Sorkin 1987/1994 (V = T), Smolin 2009/2010.
  Unimodular quantum cosmology with the clock conjugate to Lambda is also in Chiou-Geiller 2010
  [Q abstract]. GRUT adds only the reading that DEF 1.6's external t is the slot being closed.

## 3. Route R1a: a state principle for the whole

### 3.1 Which state rules survive the shift? [D, C: C1-C4]

Take an equivariant rule rho = f(H)/Tr f(H) with f a fixed function. Invariance under H -> H + cI for
all c requires f(x + c) = g(c) f(x). Differentiating in c at 0 gives f' = -beta f, so f = A e^{-beta x}
(C4).
- **Numerically (C1-C2).** Gibbs deviates by 7e-16. Square, logistic and Lorentzian rules deviate by
  0.52-0.81.
- **Anchored rules (C3).** A rule that uses a spectral anchor, f(H - E0(H) I) with E0 the ground energy
  (for example the ground-state projector), is also invariant, to 1e-15. It needs the anchor to exist.
- **In R1.** The generator h0 + rho_vac has no ground state, because Lambda_eff is unbounded above, so
  it has no bottom anchor. The only anchors available come from gravitational dynamics: h0's zero
  (Lambda_eff = 0), or an AdS-side extreme for fixed matter content (Lambda_eff < 0). So closure-native
  rules are Gibbs rules. Anchored rules are dynamics-native and select 0 or a negative value.

### 3.2 Gibbs/KMS state of the whole with respect to the 4-volume clock [D, C: D1-D4]

The generator eigenvalue is rho_vac - Lambda_eff, so e^{-beta h} is proportional to e^{+beta Lambda_eff}.
Write lam = Lambda_eff/M^4. Two densities of states:

- **(a) Gibbons-Hawking weight, g = exp(S_dS), S_dS = 24 pi^2/lam.** This is imported. It is the
  Baum-Hawking / no-boundary weight.
  - Z(beta) = ∫ exp(24 pi^2/lam + beta lam) dlam diverges at lam -> 0+ for either sign of beta.
  - With a cutoff lam_min: ln Z = 24 pi^2/lam_min to 1% for lam_min = 0.1 to 1e-4 and beta = ±1 (D1).
  - All the weight lies within a factor 2 of the cutoff, P = 1.0 (D2).
  - There is no normalizable KMS state of the whole. Regulated, it selects Lambda_eff -> 0+, which is
    Baum-Hawking (credited; value excluded since 1998, retrospective).
  - At rho_obs the exponent is 3.31e122 (D0, as in the siblings).
- **(b) Flat density on lam > 0.** It is normalizable only for beta < 0, and then <Lambda_eff> = 1/|beta|
  (D3). The value moves into the clock temperature beta, which has dimension (4-volume), a supplied
  number.
- **(c) Self-consistent KMS.** Fix beta from the state's own scales: beta = kappa M^-a H^(a-4), and
  require 1/beta = 3M^2 H^2.
  - a = 0, 1, 3, 4: H = sqrt(3 kappa) M, 3 kappa M, M/(3 kappa) and M/sqrt(3 kappa) respectively, all
    Planckian for kappa = O(1).
  - a = 2: an identity when kappa = 1/3, with H free (D4).
  - To land on the observed H with a = 0 needs kappa = 1.16e-121. The large number reappears in the
    coefficient.
- **(d) Maximum entropy of the whole.** With S_dS as the entropy, this is (a) at beta = 0, so
  Lambda_eff -> 0+. Bounding the entropy by a finite ln d (DEF 1.1 finite case) relocates the value to d
  (Banks/Bousso N-bound; `obstruction` §6.3: ln d = 3.3e122 needed, retrospective).
- **(e) Eigenstate or normalizable packet.** Any value, or the mean is free (`unimodular` §3.2-3.3).
- **(f) No-boundary or tunnelling.** 0+ or the largest value (`others` §6; `unimodular` §3.3). GRUT's
  C1/C2 work found no selection rule (s09). Not repeated.

### 3.3 Gate answers for R1a

- **Parameters.** beta and the density of states (Gibbs); the cutoff; d; the prescription.
- **Q1.** R1's clock, plus a state rule in the preparation slot.
- **Q2.** It fixes 0+ (wrong; excluded, retrospective), or a value set by a supplied beta or d, or
  nothing.
- **Q3.** Rules act on Lambda_eff, so they are stable (§2.3). The supplied beta or d is not
  renormalized by matter loops in this setting [A].
- **Q4.** Lambda_eff = 0+ is contradicted by the observed acceleration (exposed). Otherwise none.
- **Novelty.** (ii) Baum 1983, Hawking 1984, Ng-van Dam 1990, Banks 2001, Bousso 2000. Or (iii): a
  supplied beta or d is a value statement. The Gibbs-uniqueness step (§3.1) is elementary (Cauchy's
  functional equation) and is not claimed as new.

## 4. Route R1b: a fluctuation principle

### 4.1 Does GRUT contain discreteness that gives Delta V ~ sqrt(V)? [A; record read]

Inventory of what GRUT has:
- **FORMAL_CORE DEF 1.1 finite H_Omega.** A regulator (§6 item 1: "with its regulator if
  continuum"). Hilbert-space dimension carries no number-volume correspondence.
- **A finite-dimensional clock.** In the finite class the clock is finite-dimensional, so h is bounded
  and COND 2.20(ii) makes the T-clock imperfect. A clock with d_T distinguishable readings over the
  history resolves Delta T ~ V/d_T. Matching the Poisson spread Delta T ~ sqrt(V l^4) needs
  d_T ~ sqrt(V M^4) = 1.0e120 [C: E5], which is a supplied number.
- **CP-12 graph Hilbert space (ledger s06).** Bounded-degree, evolved in external time. GRUT III §3 shows
  that an external clock picks a frame, and Bombelli-Henson-Sorkin's no-preferred-direction result
  concerns Poisson sprinklings, which CP-12 is not. CP-12 gives no Poisson number-volume relation.
- **GRUT III 13_SYNTHESIS §7 layer 1.** "causal order plus counting, not adjacency", with "number ~
  volume" holding "only statistically". Its conclusion: "It adds nothing yet to layer 1, where causal-set
  theory has the right kinematics and GRUT has no dynamics of its own." N3 defers further work to
  Mathur-Singh-Surya and Carlip-Carlip-Surya.
- **COND 2.13 (Wigner-Araki-Yanase, Ozawa) applied to the T-clock [D/A].** Lambda is conserved, and
  T does not commute with it. A record of elapsed 4-volume made inside the whole with error eps_V
  therefore needs a Lambda_eff spread >~ 1/(2 eps_V) in the whole's state. The observed age is known to
  ~0.2%, which gives a spread >~ 250/V ~ 3e-238 M^4 [E, using E1]. That is about 117 orders below
  rho_obs. A GRUT-flavoured reading, empty for the scale.

**Answer.** GRUT supplies no discreteness or finite-information assumption that gives
Delta V ~ sqrt(V). Adopting GRUT III's own layer-1 deferral imports Sorkin's Poisson kinematics. With
R1 and a zero-mean postulate, the result is everpresent Lambda: rho_Lambda = xi alpha M^2/sqrt(V).

### 4.2 Pre-registered check K1 and what it shows [C: E1-E4]

- **V_cone.** V_cone H0^4 = 0.1229, i.e. V_cone M^4 = 1.01e240 (reproduces `unimodular` §4.2). The age
  is t0 H0 = 0.950.
- **Exponent p in rho_obs/M^4 = (V M^4)^(-p):**

  | V | p (reduced) | p (non-reduced) | alpha_red for 1 sigma |
  |---|---|---|---|
  | past light cone | 0.5006 | 0.5064 | 0.720 |
  | H0^-4 | 0.4987 | 0.5045 | 2.05 |
  | t0^4 | 0.4989 | 0.5046 | 1.85 |

- **K1 passes:** p is in [0.45, 0.55] for every convention.
- **K1 is weak, and I say so.** A window of ±0.05 in p spans an amplitude factor of 1e12 (E3). K1
  confirms only rho ~ M^2 H^2, which every V^(-1/2) law shares. Its one sharp content is that a sharp
  Fourier clock, Delta V ~ V, would give rho ~ 1/V = 1e-240 M^4 (E4). Among power laws, only the
  central-limit exponent 1/2 lands near rho_obs. That fact is Sorkin's 1987-1997 order-of-magnitude
  statement. It predates 1998 but is retrospective for us.

### 4.3 Gate answers for R1b

These are as scored by the `unimodular` track §4.7, which I do not recompute.
- **Q1.** Poisson discreteness plus conjugacy.
- **Q2.** The scale is fixed, but the coefficient (alpha_red = 0.72 for 1 sigma) is not. The sign is
  random. The mean is postulated to be zero. In R1 the zero is at least shift-invariant (the h0
  anchor, §3.1), so it needs no retuning. It is still a value statement.
- **Q3.** The spread is stable. The mean is stable only as a statement about the total.
- **Q4.** Stochastic dark energy (§8). The simplest realizations are strained by exposed data.
- **Novelty.** (ii): Sorkin; Ahmed-Dodelson-Greene-Sorkin 2004; Zwane-Afshordi-Sorkin 2018;
  Das-Nasiri-Yazdi 2023-24. A symmetry argument for the zero mean (energy -> -energy) also exists in
  the literature (Kaplan-Sundrum 2006 [Q abstract]). GRUT supplies none.

### 4.4 R1c, a variant: the access-restricted clock [A]

GRUT's DEF 1.2 lets a prediction use only observables in the access structure. HT's T(Sigma) is a
whole-slice integral (§2.3, B3). The clock an observer can read is therefore the 4-volume of their own
causal past.
- That is the V choice ADGS made ad hoc, which they themselves called "somewhat at odds with the
  meaning that V has in the unimodular context" [Q via `unimodular` §4.1].
- It is also the causal restriction of Barrow-Shaw's variation (`others` §8.1). Their Omega_k0 =
  -0.0055 (2 zeta_b) prediction is disfavoured at 7.7 sigma by DESI DR2+CMB (retrospective).
- So GRUT's access rule supplies a reason for an existing choice. It is not a new principle (ii).

## 5. Route R2: global self-consistency

### 5.1 The minimal one-variable global constraint [D, C: F1-F3]

Add one global variable Lambda (no derivatives) and one global term:

    S = ∫ sqrt(-g) (M^2 R/2 - Lambda + L_m) + sigma(Lambda/mu^4),
    =>  M^2 G = tau - Lambda_eff g,   Lambda_eff = Lambda + rho_vac,   sigma'(Lambda/mu^4) = mu^4 V4.   (5.1)

Hold V4 as a function of Lambda_eff, that is, of the history. Then dLambda_eff/drho_vac =
sigma''/(sigma'' - mu^8 V4'(Lambda_eff)) (F1).
- **Nonlinear sigma.** This is nonzero (F3), so rho_vac leaks into Lambda_eff and the obstruction is
  not escaped.
- **Linear sigma, sigma = v Lambda/mu^4 (F2).** The derivative is exactly 0. The Lambda-equation becomes
  **V4 = V0 = v/mu^4: the total 4-volume of the history is fixed.** The action changes by
  c(V0 - V4) = 0 on shell, and Lambda_eff is whatever value makes the history's 4-volume equal V0.
  This is HT's (2.1) with Dirichlet data on T at both ends, since T_f - T_i = V4.

**Conjugate-pair relocation [D].** (Lambda_eff, V4) is the HT canonical pair. R1 supplies
phi(Lambda_eff). R2a supplies V0, the Fourier-dual boundary datum. R1b supplies a relation between
their spreads plus a mean. In each case the energy-origin obstruction moves into the supplied item.

### 5.2 Is a fixed finite V0 compatible with the observed expansion? (K3) [D, C: G1-G2]

Finite V4 needs compact spatial sections and a crunch.
- **Constant Lambda_eff > 0 [C: G1].** With Omega_Lambda in [0.60, 0.75] and Omega_k in
  [-0.05, 0.05], E^2 a^3 = Omega_m + Omega_k a + Omega_Lambda a^3 has no zero for a > 1. The cubic is 1
  at a = 1 and increasing whenever 3 Omega_Lambda > -Omega_k. So V4 = infinity: 0 of 1271 grid points
  recollapse.
- **Consequence.** R2a with its own constant Lambda_eff is retrospectively incompatible with the
  present acceleration. It survives only with an added collapse sector (dark-energy dynamics that turn
  the expansion around), which is the same requirement as global sequestering.
- **Closed dust, to show the map [C: G2].** Units: 3M^2 = 1, A = rho_m a^3 = 1, lam = Lambda_eff/(3M^2),
  lam_E = 4/27.

  | lam/lam_E | -10 | -1 | 0 | 0.5 | 0.9 | 0.99 | 0.999 |
  |---|---|---|---|---|---|---|---|
  | V4 | 3.32 | 18.9 | 33.9 | 56.7 | 144 | 333 | 555 |

  - The lam = 0 value is the closed-dust 35 pi^3/32 = 33.91.
  - V4 increases monotonically and diverges at lam_E.
  - So V0 <-> Lambda_eff is one-to-one at fixed matter content: V0 is Lambda_eff in other units, and
    the value is supplied.

### 5.3 "Bath = rest" read as a relational vacuum: measure the vacuum relative to the whole [D, C: H1-H3]

The reading: a subsystem's vacuum energy is physical only relative to its bath (Casimir structure,
`obstruction` §1.5). The whole has no bath, so its vacuum level should not gravitate. Covariantly,
subtract the whole's average.
- **The trace subtraction (KP).** The only linear covariant scalar that sees -C g is the trace, so the
  minimal form subtracts (1/4)<T>:

      M^2 G = T - (1/4) <T> g   =>   <R> = 0,   Lambda_eff = <tau>/4.                    (5.2)

  This is exactly KP's global sequestering bulk equation (sequestering §3.1) [D].
- **It is not the only realization.** Tseytlin's S/V action (PRL 66, 545 (1991) [Q abstract]) gives
  M^2 G = T - <M^2 R/2 + L_m> g. With the trace, M^2<R> = 4 Lambda_eff - <tau>, so
  **Lambda_eff = <tau/2 - L_tau>** (H1). Both subtractions absorb any c exactly (H2). Per unit
  <rho> (H3):

  | fluid | KP (1/4)<tau> | Tseytlin, L = p | Tseytlin, L = -rho |
  |---|---|---|---|
  | dust | -1/4 | -1/2 | +1/2 |
  | radiation | 0 | -1/3 | +1 |

  The two agree on pure vacuum (T = 4L there) and differ on matter. Tseytlin's carries the
  on-shell-Lagrangian ambiguity found for Smolin's relation (`unimodular` §3.4), and is the same object.
- **From the equation to an action.** Equation (5.2) needs two global variables: one coupled to
  ∫sqrt(-g), one to ∫sqrt(-g) T. That is KP's (Lambda, lambda) pair with the protected sector on
  lambda^2 g. The alternative is a nonlocal ratio (Tseytlin). One variable gives only §5.1.
- **Equivalence statement.** The minimal GRUT-motivated global constraint, one variable, is **not**
  equivalent to Kaloper-Padilla. It is HT with a fixed V0. The relational reading reproduces KP's bulk
  equation only after choosing the trace subtraction and adding KP's second variable and scaling
  assumption. GRUT closure motivates neither choice, and the closure fork's "nothing outside" does not
  imply finite 4-volume (FORMAL_CORE §4, SCOPE row).

### 5.4 Gate answers for R2

- **R2a.**
  - Parameters: v (that is, V0), mu, the matter content, a collapse sector.
  - Q1: a global multiplier for ∫sqrt(-g).
  - Q2: not fixed. Lambda_eff = F(V0, history).
  - Q3: matter-loop constants are absorbed (F2); graviton loops likewise, at the effective-action level
    [A].
  - Q4: the universe must recollapse. With a constant Lambda_eff > 0 that is retrospectively excluded
    (G1).
  - Novelty: (ii) HT 1989 with fixed boundary data; a global unimodular constraint.
- **R2b-KP.**
  - Parameters: sigma, mu, the protected-sector scaling, the collapse sector and its cutoff.
  - Q1: two global variables.
  - Q2: the residual is a history functional; the linear trigger fails <R> = 0 (sequestering §5).
  - Q3: matter loops yes; graviton loops no.
  - Q4: crunch and transient acceleration.
  - Novelty: (ii), KP.
- **R2b-Tseytlin.**
  - Q2: sign and size depend on the on-shell L choice.
  - Q4: none stated.
  - Novelty: (ii), Tseytlin 1991 and ADDG.

## 6. Route R3: closure of the energy budget

- **(a) With gravity included [D, C: I1].** On a compact slice the whole's energy is the Hamiltonian
  constraint, which vanishes on every solution. Symbolically, 3M^2 H^2 - rho_m - Lambda = 0 holds
  identically on the flat dust + Lambda solution for every Lambda. "The vacuum is the zero of the
  whole's energy" is then satisfied for every Lambda: circular, (iii). Credit Tryon 1973 [recall] and
  `obstruction` §2.1.
- **(b) With gravity excluded, on a slice [D].** Require ∫_Sigma (rho_m + rho_Lambda) sqrt(h) = 0. The
  constraint gives rho_tot = 3M^2(H^2 + k/a^2).
  - k = +1: this is positive on every slice, so the condition can never hold.
  - k = 0: it holds only at a turnaround, with rho_Lambda = -rho_m(t*) < 0.
  - k = -1: it holds where adot = 1, also with rho_Lambda < 0.
  - It also needs a preferred slice (class B in the `obstruction` list).
- **(c) With gravity excluded, averaged over the history [D, C: I3-I4].** rho_Lambda = -<rho_matter>_V4
  <= 0 for matter satisfying the weak energy condition. That is the wrong sign against the observed
  acceleration (retrospective). It is also not self-consistent:
  - Flat dust with Lambda < 0 gives <rho_m>/|rho_Lambda| = 2.000 at every scale (I3).
  - Adding radiation (ratio r at turnaround) gives 2.00, 2.01, 2.12, 2.54 and 2.92 for r = 0, 0.01,
    0.1, 1 and 10 (I4). The required value 1 is never reached, and no scale enters.
  - The weighting uses rho = T_mu nu u^mu u^nu, so a frame must be supplied.
- **(d) "Relative to what?"**
  - Relative to the bath: §5.3, i.e. KP or Tseytlin.
  - Relative to an asymptotic vacuum: a closed whole has no asymptotic region, and asymptotic flatness
    fixes 0 by fiat (`obstruction` §2.1).
  - Relative to the whole's ground state: the constraint has none.
  - Relative to h0's zero in R1: Lambda_eff = 0.
- **Gate.**
  - Q1: none in (a); a frame plus a global constraint in (b)-(c).
  - Q2: nothing in (a); the wrong sign in (b)-(c).
  - Q3: (c) as a global constraint would absorb loops, but it is moot.
  - Q4: rho_Lambda <= 0, excluded (retrospective).
  - **Verdict:** (iii) circular with gravity; sign-excluded without it.

## 7. The result: the lifted energy-origin obstruction for GRUT [D/A]

**Hypotheses.**
- (H-a) A GRUT model instance Theta whose clock slot is closed in the minimal way that makes the energy
  origin dynamical: an internal clock conjugate to the coefficient of ∫sqrt(-g) (HT fields), or the
  equivalent one-variable global constraint.
- (H-b) Matter is a perturbative local QFT, and gravity is a local diffeomorphism-invariant effective
  action, apart from the global HT variables.
- (H-c) Every additional rule is built from GRUT closure data: the FORMAL_CORE theorems, restriction,
  the access structure, the HT kinematics, and Gibbs/KMS or entropy conditions of the whole.

**Statements.**
1. **[D, C: A1, B2]** Every constant vacuum contribution enters as h -> h + cI on the internal-clock
   generator, i.e. a clock phase. The gravitating Lambda_eff is the h0-spectral label carried by the
   state. Codex's law-level obstruction is escaped, which answers Q1 and Q3 for constant terms.
2. **[A/D, §2.4; `unimodular` lemma]** No rule under (H-c) selects Lambda_eff's mean.
   - Theorems quantified over Theta hold in every Lambda_eff sector.
   - Kinematic rules are translation-invariant and fix at most a spread.
   - The unique closure-native equivariant shift-invariant rule is Gibbs. It is non-normalizable with
     the GH density of states and concentrates at 0+ (D1-D2), or its value is a supplied beta (D3).
3. **[D, C: D4, E3-E4; credit `obstruction` item R]** Any selection with O(1) coefficients built from M
   and the state's own Lambda_eff (or H) gives 0, a Planckian value, or an identity. Reaching
   7.15e-121 requires one supplied element: a dimensionless input of size 1e120-1e122
   (sqrt(V M^4) ~ 1e120, ln d ~ 3e122, |beta| M^4 ~ 1e120, d_T ~ 1e120), an exponential suppression
   e^-276, or a fluctuation law with exponent exactly 1/2 together with a stated mean.
4. **[D, §5.1]** For the conjugate pair (Lambda_eff, V4), closure can turn one of them, or a relation
   between their spreads, into a datum. It cannot fix both. Whichever is supplied is where the
   obstruction now sits.
5. **[record]** GRUT's Theta contains none of the items in statement 3:
   - the state slot is empty (C1, C2, s09);
   - finite d is a regulator;
   - layer 1 is deferred to causal sets (GRUT III §7);
   - finite 4-volume is not implied (FORMAL_CORE §4).

**Therefore** a GRUT vacuum principle that fixes the total renormalized rho_vac/M^4 must add, in some
slot of Theta, a Lambda_eff-dependent weight or one of the supplied elements of statement 3. Every such
addition located here or in the sibling tracks already has a published owner:
- HT, Unruh;
- Baum-Hawking, Ng-van Dam;
- Banks-Bousso;
- Sorkin, ADGS;
- Tseytlin, KP;
- Barrow-Shaw.

No assumption that the cited literature lacks was found. This is the scoped result. It is not a no-go
theorem for quantum gravity, for nonperturbative or nonlocal UV physics, or for principles that lie
outside (H-b) and (H-c).

## 8. Pre-registered route: outcome against §P

- **Checks K1-K6.** All ran.
  - K1 passed, but weakly (§4.2).
  - K2: no normalizable Gibbs state with the GH density; with a flat density, the value moves into beta.
  - K3: incompatible with a constant Lambda_eff > 0 and the present acceleration.
  - K4: circular with gravity; rho_Lambda <= 0 without it.
  - K5: only exponential rules are invariant, plus anchored ones.
  - K6: the matter shift is a pure clock phase; only rules anchored on h0 are invariant.
  - None of these changes the pre-registered choice of route.
- **Status of the pre-registered consequence.** The stochastic, sign-changing rho_DE tracking
  rho_crit is Sorkin's everpresent Lambda. The `unimodular` track's converged Model 1 Monte Carlo
  (20,000 realizations per alpha) finds:
  - P(Omega_Lambda0 in [0.64, 0.74]) <= 1e-4 at every alpha scanned;
  - 0 of 180,000 realizations meeting the joint today/recombination/BBN window.
  Das et al. 2024 report a CMB chi^2 more than 150 worse than LCDM for Model 1 [Q via `unimodular`].
  All of this is retrospective. The route is therefore already strained by exposed data in its
  simplest form.
- **What a held-out test would need.** Freeze, before any data released after DESI DR2 are used, all
  of the following:
  - alpha;
  - the start epoch;
  - the update recipe;
  - any suppression near last scattering;
  - the predicted two-point function of rho_DE(ln a) (correlation sqrt(V1/V2)).

  Then compare with redshift-binned rho_DE(z) from DESI after DR2, Euclid, Rubin LSST and Roman, plus
  z > 2 tomography (21-cm intensity mapping). Each later modification of the recipe is a supplied
  choice, and a "good seed" chosen with exposed data is a fit.
- **Exposed, never held out.**
  - Planck 2018, ACT DR6 and SPT-3G CMB.
  - DESI DR1 and DR2 BAO, including the w0 > -1, wa < 0 preference.
  - Pantheon+, Union3 and DES-Y5 supernovae.
  - BBN abundances, the distance-ladder H0, and the DESI+CMB Omega_K = +0.0023 ± 0.0011.

## 9. Scorecard

| Route | Slot of Theta changed | Equations | Parameters / inputs | Q1 extra structure | Q2 fixes residual? | Q3 stable at stated order? | Q4 independent consequence; exposed data | Novelty | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| R1 clock closure | clock (+ algebra: HT fields) | (2.1), (2.2) | M, HT fields, state phi(Lambda_eff), ordering, a = 0 BC | internal 4-volume clock conjugate to the origin | no: free state label | yes for constant terms (clock phase, B2), if the state is given by its h0 support (B4-B5); full quantum HT disputed | none: GR + Lambda per sector (COND 2.18 reading) | (ii) HT, Unruh, Unruh-Wald, Sorkin, Smolin, Chiou-Geiller | escapes the law-level obstruction; relocates it to the state |
| R1a Gibbs/KMS w.r.t. the T-clock | preparation | §3.2 | beta, density of states, cutoff | R1 + Gibbs rule (the unique shift-invariant equivariant rule, C4) | GH density: 0+ (non-normalizable, D1-D2); flat: 1/abs(beta) (supplied) | yes (acts on Lambda_eff) | 0+ excluded (retrospective) | (ii) Baum-Hawking / (iii) | fails Q2 |
| R1a max-entropy / finite d | preparation (+ algebra: d) | S_dS = 24 pi^2 M^4/rho | d | entropy bound | 0+, or relocated to ln d = 3.3e122 | yes | none | (ii) Banks, Bousso | relocation |
| R1a self-consistent KMS | preparation | 1/beta(M,H) = 3M^2H^2 | kappa, a | dimensional closure | Planckian or identity (D4) | n/a | none | (iii) | empty |
| R1a anchored on h0 = 0 | preparation | h0 phi = 0 | none | gravitational anchor | Lambda_eff = 0 (wrong) | yes | Lambda = 0 excluded | (ii) | fails |
| R1b counting clock | algebra (layer 1) + preparation (mean) | rho = xi alpha M^2/sqrt(V) | alpha, l, mean, start, recipe, V choice | Poisson discreteness (not in GRUT; GRUT III defers it) | scale yes (K1, coincidence form); coefficient free (0.72); sign random; mean postulated | spread yes; mean only as a value statement | stochastic DE; P(our history) <= 1e-4 [`unimodular` MC, retrospective] | (ii) Sorkin, ADGS, ZAS, Das et al. | most promising; strained |
| R1c access-restricted clock | clock + access | V = past-cone 4-volume | as R1b, or Barrow-Shaw's | access rule picks causal-past V | as R1b / Barrow-Shaw relation | as R1b | Barrow-Shaw Omega_k = -0.0055, 7.7 sigma off (retrospective) | (ii) ADGS's choice, Barrow-Shaw | reason for an existing choice |
| R2a one global variable | dynamics (global term) | (5.1), sigma linear | V0 (= v/mu^4), matter, collapse sector | multiplier for ∫sqrt(-g) | no: Lambda_eff = F(V0); nonlinear sigma leaks rho_vac (F3) | yes (F2) | must recollapse; with constant Lambda_eff > 0 retrospectively excluded (G1) | (ii) HT with fixed boundary data | relocation to V0 |
| R2b relational vacuum, trace form | dynamics (two global variables) | (5.2) | sigma, mu, lambda-scaling, collapse sector, cutoff | KP global variables | history functional <tau>/4; linear trigger fails <R> = 0 (sequestering §5) | matter loops yes; graviton loops no | crunch, transient w | (ii) KP (Tseytlin/ADDG lineage) | not GRUT-forced |
| R2b relational vacuum, action form | dynamics (nonlocal) | S/V | on-shell L choice | nonlocal ratio action | <tau/2 - L_tau>: sign depends on L choice (H3) | constants yes | none stated | (ii) Tseytlin 1991; = Smolin / Ng-van Dam form | ill-defined residual |
| R3 energy budget with gravity | none | H ≈ 0 | none | none | no (identically satisfied, I1) | n/a | none | (iii) Tryon | circular |
| R3 energy budget without gravity | dynamics + frame | rho_Lambda = -<rho> | frame, averaging domain | global constraint + preferred frame | sign <= 0; no self-consistent flat dust/radiation solution (I3-I4) | moot | rho_Lambda <= 0 excluded (retrospective) | (iii)/(ii) | excluded |

## 10. Limits

- **Minisuperspace.** R1 is analysed in homogeneous minisuperspace and in a finite-dimensional toy
  (B: d = 8, random h0). The toy is not a Wheeler-DeWitt operator. It checks the algebra of (2.2), not
  gravitational dynamics. Operator ordering, the a = 0 boundary condition and the self-adjoint
  extension are left unspecified.
- **The universal-quantification argument (§2.4).** Structural [A]. It covers FORMAL_CORE's theorems
  as stated. It does not cover a future GRUT rule that refers to Lambda_eff-dependent gravitational
  weights. Such a rule would be a new postulate (statement 5 of §7).
- **Gibbs/KMS (§3.2).** It uses the semiclassical GH density of states, imported, and a sharp cutoff.
  No other density was tried beyond the flat one. Beyond minisuperspace and the semiclassical weight,
  the state space of the whole is not specified, so "non-normalizable" is a statement about this
  ansatz.
- **The scale lemma (§7 item 3).** Dimensional analysis with O(1) coefficients. Exponentially
  suppressed mechanisms are named but not analysed.
- **K3 (G1).** Uses constant Lambda_eff, dust + Lambda, and Omega_k in [-0.05, 0.05]. Collapse
  sectors are not scanned here; the `sequestering` track scanned the linear trigger.
- **G2.** Closed dust only, in units A = 1. Illustrative of the V0 <-> Lambda_eff map.
- **R3 (I3-I4).** Flat dust and radiation with Lambda < 0. k = +1, stiff matter and cutoffs are not
  included.
- **R1b.** No Monte Carlo was run here. The numbers are the `unimodular` track's (retrospective), and
  everpresent-Lambda gate answers are cited, not recomputed.
- **Literature.** Read at abstract level here: Tseytlin 1991, Chiou-Geiller 2010, Smolin 2010,
  Tiwari 2006 and Kaplan-Sundrum 2006. Everything else is cited through the sibling tracks' read levels.
  Felten-Isaacman-type recollapse criteria were not consulted; G1 is a direct computation. Searches were
  one web search plus the arXiv API and INSPIRE, so a differently worded paper may have been missed.
- **Novelty assessments.** Relative to the literature located here and in the sibling tracks.
- **Claims.** No new theory, no empirical validation, no confirmatory test. All data comparisons are
  retrospective.

## 11. Sources and read levels

| Source | Identifier | Read level | Used for |
|---|---|---|---|
| GRUT FORMAL_CORE, WORKING_THEORY, CAMPAIGN_LEDGER | `research/reality_framework/grut_closed_whole_theory_20261019/` | full (FORMAL_CORE, WORKING_THEORY); ledger by search | datum, DEF 1.1-1.8, COND 2.13, 2.18, 2.20, §4, §7 |
| Closure fork | `grut_closure_fork_20261019/DERIVATION.md:197-215`, `CLOSURE_SELECTION_RESPONSE.md:140-177` | sections | prior sequestering and k > 0 record |
| GRUT III synthesis | `GRUT_III/13_SYNTHESIS.md` §§3, 7, 9 (lines 100-135, 240-300, 350-372); `01_PRINCIPLES.md` | sections | layer-1 deferral to causal sets; N2, N3; Lambda as datum |
| Codex direction | `git show 02b4515:.../NEXT_DIRECTION.md` | full | obstruction and gate |
| Sibling tracks | `../obstruction`, `../sequestering`, `../others`, `../unimodular` NOTES.md | full | cited results; their read levels apply |
| Tseytlin | PRL 66, 545 (1991), DOI 10.1103/PhysRevLett.66.545 | abstract (INSPIRE/APS) | S/V action; equations in §5.3 derived here |
| Chiou, Geiller, "Unimodular loop quantum cosmology" | PRD 82, 064012 (2010), arXiv:1007.0735 | abstract | unimodular clock in quantum cosmology; semiclassicality needs small <Lambda> |
| Smolin, "Unimodular loop quantum gravity and the problems of time" | arXiv:1008.1759 | abstract | Unruh-Wald-Sorkin time proposal |
| Tiwari, "Thermodynamics of spacetime and unimodular relativity" | gr-qc/0612099 | abstract | Lambda thermodynamic fluctuations suggested |
| Kaplan, Sundrum, "A symmetry for the cosmological constant" | JHEP 0607, 042 (2006), hep-th/0505265 | abstract | energy -> -energy symmetry (zero-mean precedent) |
| Henneaux-Teitelboim; Unruh; Unruh-Wald; Kuchar; Sorkin; ADGS; ZAS; Das et al.; Smolin 2009; Padilla-Saltas | see `unimodular` §1 | as recorded there | R1, R1b |
| Kaloper-Padilla (global, local, OS); Avelino; Ferreira-Avelino | see `sequestering` §1 | as recorded there | R2b |
| Baum; Hawking; Coleman; Bousso; Banks; Barrow-Shaw; ADDG | see `others` §12 and `obstruction` §9 | as recorded there | R1a, R1c, R2b |
| Page-Wootters; Wigner-Araki-Yanase; Ozawa; Unruh-Wald | as cited in FORMAL_CORE DEF 1.6, COND 2.13, 2.20 | via FORMAL_CORE | clock closure, record limits |
| Tryon | Nature 246, 396 (1973) | recall | R3 |

Downloads are kept in the scratchpad (`.../scratchpad/gn_dl/`). No paper text is in the repository.

## 12. Files

- `code/grut_native_checks.py`: sections A-I; numpy, scipy and sympy; seed 20261021; 29 checks.
  Run `python3 -I code/grut_native_checks.py results/grut_native_checks.json`.
- `results/grut_native_checks.json`: every [C] number above, the check list (29/29 pass) and versions.
- `results/prereg.json`: SHA-256 of §P as written before any code. To verify, hash the text from
  "## P." up to, but not including, the newline that precedes "## 0. Result".
