# Track `review_scope`: skeptical referee report on novelty, scope and testability

Packet `grut_vacuum_principle_20261021`, 21 October 2026. Written under `../../BRIEF.md`.

**Read:**
- `BRIEF.md`;
- all five track notes in full: `obstruction`, `sequestering`, `unimodular` (751 lines, read in two parts),
  `others`, `grut_native`;
- `grut_native/results/prereg.json`;
- Codex's packet (`git show 02b4515:research/reality_framework/grut_scale_feasibility_20261021/NEXT_DIRECTION.md`).

**Inspected** for load-bearing checks: the WDW code in `obstruction/code/obstruction_checks.py` (lines 295-313),
the track result JSONs, file timestamps and `git log` of this packet.

Labels: [D] derived here; [C] computed by `code/review_checks.py` (record `results/review_checks.json`, 17/17
checks pass, deterministic, no Monte Carlo, about 20 s); [Q] quoted at the read level in §7; [E] estimate;
[A] argued. All data comparisons are **retrospective**.

This track builds no candidate and claims no new theory, empirical validation or confirmatory test. No wall,
population or Candidate 1 run was launched.

## 0. Verdict in brief

1. **Novelty (question 1).** No track claims a GRUT-native route, and I found none hidden.
   - Every GRUT-motivated route relabels a published one:
     - R1 is Henneaux-Teitelboim / Unruh / Unruh-Wald / Sorkin;
     - R1a is Baum-Hawking / Banks-Bousso;
     - R1b is Sorkin's everpresent Lambda;
     - R1c is ADGS / Barrow-Shaw;
     - R2a is HT with fixed total 4-volume;
     - R2b is Kaloper-Padilla / Tseytlin / ADDG;
     - R3 is Tryon / Einstein 1917.
   - **The credit lists miss the closest published closure argument.** Klinkhamer-Volovik q-theory: a
     self-sustained vacuum with no external environment has zero gravitating vacuum energy in equilibrium.
     Volovik also explicitly contests the "overall energy shift" premise of the obstruction. Any GRUT
     "closed whole ⇒ the vacuum does not gravitate" argument must be credited to them (§1).
2. **Testability (question 2).** No proposed consequence is both independent of exposed data and written
   before the construction was adjusted.
   - The one pre-registration (`grut_native` §P) was written after the sibling everpresent-Lambda Monte Carlo
     was on disk. It freezes no number.
   - Its statistic (binned rho_DE(z), 0 < z < 2.5) is already measured by exposed DESI DR2 and SN data, on the
     same single realization.
   - My linear estimate puts the route at or beyond its own pre-registered failure level on exposed data [E]
     (§2).
3. **Overall conclusion (question 3).** The scoped obstruction and its "relocation, not refutation" reading are
   sound. Several statements overreach or are wrong:
   - "For Lambda > 0 there is no normalizable spectrum" (WDW) is false. The endpoint is limit-circle [C] (§3.1).
   - "R1 escapes the obstruction; Q1 and Q3 pass" overstates. Q2 fails, so Q3 has no value to protect.
   - The minimal-escape statements of `obstruction` and `grut_native` disagree on the R item, and their class
     lists omit published routes.
   - Four tracks name four different "best" routes.
   - The Barrow-Shaw "7.7 sigma" is specific to one dataset.
4. **Best route.** None. No surveyed route fixes the total renormalized rho_vac/M^4, stays stable without
   retuning, and gives a separate test. The packet's honest outcome is **an identified, scoped obstruction plus
   a credited map of escape classes**, with per-criterion leaders named only together with their failures
   (§4).

## 1. Question 1: is anything claimed as GRUT-native actually new?

### 1.1 Route-by-route (grut_native's own routes, checked against the cited literature) [A]

| Route | Assumption added | Already in the literature | GRUT content left |
|---|---|---|---|
| R1 clock closure | 4-volume clock conjugate to Lambda; HT fields | HT 1989; Unruh 1989; Unruh-Wald 1989; Sorkin 1987/1994 ("V = T"); Smolin 2009; Chiou-Geiller 2010 | labelling FORMAL_CORE DEF 1.6 as "the slot being closed" |
| R1a Gibbs/KMS, max entropy, anchored state | state rule on Lambda_eff | Baum 1983; Hawking 1984; Ng-van Dam 1990; Banks 2001; Bousso 2000 | the Cauchy-equation uniqueness of Gibbs rules (elementary; not claimed) |
| R1b counting clock | Poisson discreteness, zero mean | Sorkin; ADGS 2004; ZAS 2018; Das et al. 2023-24; zero mean also Kaplan-Sundrum 2006 | adopting GRUT III's own deferral to causal sets |
| R1c access-restricted clock | V = 4-volume of the observer's causal past | ADGS's choice; Barrow-Shaw's past-cone variation | see §3.6: GRUT's access structure is metric-free, so this is an added identification, not a "reason" |
| R2a one global variable, linear sigma | fixed total 4-volume V0 | HT with Dirichlet data; global unimodular constraint | none |
| R2b "bath = rest", relational vacuum | subtract the whole's average, or Tseytlin's S/V | KP 2014; Tseytlin 1991; ADDG 2002; **also Klinkhamer-Volovik 2008 (not credited)** | none |
| R3 energy budget of the whole | the whole's energy is zero | Tryon 1973; Einstein 1917; **the equilibrium version is q-theory (not scored)** | none |

The `obstruction` and `unimodular` tracks add no GRUT-native claim:
- The Lambda-translation lemma is credited to HT, Unruh, Smolin eq. 24 and Das et al.
- The escape taxonomy G/B/S/X/R/P is offered as [A].

The taxonomy restates the structure of Weinberg 1989 §§6-8 and of the Padilla (1502.05296) and Burgess
(1309.4133) reviews. It should be credited as a restatement, not presented as the packet's own classification.

### 1.2 The missing precedent: q-theory is the published "closed whole" argument [Q abstracts; A]

**Volovik, gr-qc/0604062 (IJMPD 15, 1987 (2006)), abstract.**
- Lists as "myths" the "dependence of the vacuum energy on the overall shift of energy" and "the absolute value
  of energy only has significance for gravity", plus the Casimir-plates argument.
- Argues that they fail once the thermodynamics of a self-sustained quantum many-body ground state is used.
- This is the leading published counter-position to Codex's framing. It also bears directly on `obstruction`
  §1.5 (the Casimir lemma).

**Klinkhamer-Volovik, PRD 77, 085015 (2008), arXiv:0711.3170 [abstract].**
- A conserved vacuum variable q characterizes a "self-sustained" Lorentz-invariant vacuum.
- Its self-tuning "nullifies the effective energy density which enters the low-energy gravitational field
  equations".
- The observed value is attributed to an imperfect (perturbed) vacuum.

**PRD 78, 063528 (2008), arXiv:0806.2805 [abstract].**
- A four-form realization.
- The vacuum energy oscillates and decays.
- Its present average is claimed to agree in order of magnitude with the observed density, which is the
  coincidence form and retrospective.

Where this sits, and why it matters [A]:
- In the `obstruction` taxonomy it is **S** (Lambda dynamical through the conserved q, a 4-form) plus an
  **equilibrium rule**. "No external environment, so zero external pressure" acts like a boundary or value
  condition and selects 0. A non-equilibrium residual is added on top.
- It does not refute the scoped obstruction. The volume or q is a dynamical variable coupled to the constant,
  as in `obstruction` §1.5.
- It is, however, almost verbatim the reading `grut_native` §5.3 calls "bath = rest": "The whole has no bath,
  so its vacuum level should not gravitate."
- It is a third covariant realization besides KP and Tseytlin: local, through a 4-form and Gibbs-Duhem, not a
  spacetime average. So §5.3's "two inequivalent covariant realizations" is incomplete.
- `grut_native` R3(d) says "relative to the whole's ground state: the constraint has none". q-theory answers
  "relative to the equilibrium of the self-sustained vacuum", which gives 0. R3's verdict ("circular or wrong
  sign") therefore misses the strongest published version of the energy-budget route.
- That version gives the wrong value (0) in equilibrium, as Baum-Hawking does. Its residual is
  dynamics-dependent.

### 1.3 Other unscored escape routes (affect the completeness of the class lists) [Q abstracts; C]

- **Wang-Zhu-Unruh, PRD 95, 103504 (2017), arXiv:1703.00543.**
  - Large, inhomogeneous vacuum fluctuations give H ∝ Lambda_UV exp(-beta sqrt(G) Lambda_UV).
  - With Lambda_UV = M_P, matching H0 needs an exponent ln(M_P/H0) = 140 [C], an O(100) parameter.
- **Carlip, PRL 123, 131302 (2019), arXiv:1809.08277.**
  - A large Lambda is hidden in Planck-scale spacetime foam and averages away macroscopically.
- **Classification.** Both violate homogeneity of the metric at the Planck scale, i.e. W3 and the semiclassical
  <T> coupling in H2. They belong in class B and are not credited in any track.
- **Abbott, PLB 150, 427 (1985) [recall], and Graham-Kaplan-Rajendran, PRD 100, 015048 (2019),
  arXiv:1902.06793 [abstract].**
  - Dynamical relaxation by a slowly rolling field: a large positive CC is driven to a small negative value,
    the universe crunches, and a bounce is needed.
  - These are class B and relevant to the Bousso-Polchinski BBN argument (`others` §9). Not scored.
- **Arkani-Hamed-Hall-Kolda-Murayama, PRL 85, 4434 (2000), astro-ph/0005111 [abstract].**
  - A matter-sector scale rho^(1/4) ~ M_EW^2/M_Pl.
  - (1 TeV)^2/M = 0.41 meV reduced (0.082 meV non-reduced) against 2.24 meV [C].
  - This shows that `obstruction`'s "every available scale is cosmological" holds only because item R is
    *defined* to be independent of matter couplings (§3.4).

## 2. Question 2: are the "independently testable consequences" independent of exposed data and pre-registered?

### 2.1 The only pre-registration (`grut_native` §P, `results/prereg.json`) [record]

- **Hash.** The text of §P hashes to the stored SHA-256 (21ab80b3...e1fed, checked here). It has not been
  edited since hashing. The hash is self-attested: it sits in the same folder, was not committed to git before
  the code ran, and rests on file mtimes only.
- **Timing.** §P was written at 20:00Z and hashed at 20:01Z. The sibling `unimodular` everpresent Monte Carlo
  results are dated 19:53 and its notes 19:57. The `obstruction` (19:27) and `others` (19:48) records also
  existed. §1 of `grut_native` cites them.
  - So §P is "before this track's code", not "before the construction was confronted with data".
  - Pre-registered check K1 (p in [0.45, 0.55]) was computable from the sibling numbers already on disk
    (V_cone M^4 = 1.01e240, alpha_red* = 0.72). It is not blind. I recompute p = 0.5006 [C].
- **Nothing numeric is frozen.**
  - §P item 3 says alpha "is calibrated once on exposed data, then frozen", with no value.
  - The start epoch, step recipe, V choice, Model 1 versus Model 2, suppression near last scattering, binning
    and likelihood are not fixed. They appear as a to-do list only in §8, after the code.
  - The sibling track shows each of these moves the outcome:
    - the start epoch moves P(today window) from <= 1e-4 to 0.016;
    - the V choice moves alpha* by factors of 2-23;
    - the step recipe changes survival by a factor of about 2.
  - These are open forking paths.
- **The construction was already adjusted to exposed data in the literature.** alpha was chosen after 1998
  (ADGS 2002), seeds were selected by MCMC (ZAS), and late starts or suppression near last scattering were
  added (Das II). The pre-registration inherits that.

### 2.2 The pre-registered statistic is already decided by exposed data [C/E, retrospective]

**Setup.**
- Model 1: rho = S/V, with S Brownian in the past-cone V. Planck-2018 LCDM background; V recomputed here.
  - H0^4 V(today) = 0.1229.
  - Corr(z = 1, today) = 0.189.
  - alpha* = 0.720, the 1-sigma-today calibration.
  - These reproduce `unimodular` §4.2-4.3.
- Conditional structure [D]: given rho at larger V, rho at smaller V has mean equal to it and variance
  sigma^2(V_small) - sigma^2(V_large) (Brownian-bridge Markov property).
- At alpha*, the conditional bin-to-bin standard deviations are about 2.3, 4.7 and 16 rho_DE0 for 0 -> 0.5,
  0.5 -> 1 and 1 -> 2.
- Exposed DESI DR2 w0wa central curves (DESI+CMB+{Pantheon+, Union3, DESY5}, Table 5 via `others`) give
  rho_DE(z)/rho_DE0 of:
  - 1.06-1.18 at z = 0.5;
  - 0.98-1.06 at z = 1;
  - 0.73-0.76 at z = 2 (extrapolated).

**Probability that every bin-to-bin change stays within delta (units of rho_DE0)** [C, `results/review_checks.json`]:

| delta | P(bins z = 0.5, 1) at alpha* | P(bins z = 0.5, 1, 2) at alpha* | max over alpha of P(today window) x P(bins z <= 1) |
|---|---|---|---|
| 0.2 | 2.4e-3 | 2.4e-5 | 1.6e-4 |
| 0.3 | 5.4e-3 | 8.1e-5 | 3.5e-4 |
| 0.5 | 1.5e-2 | 3.7e-4 | 9.7e-4 |

- **Tolerances.** delta = 0.5 is generous for z <= 1. The z = 2 bin is weakly constrained, because dark
  energy is about 6% of rho_tot there. Loosening its tolerance to 1.0 rho_DE0 gives 2.7e-4 at delta = 0.3.
- **Reading [E/A].**
  - §P's failure condition is "near-constant rho_DE(z) with probability below 1% over the model's
    realizations". It is at or beyond its threshold on **already-exposed** data at the 1-sigma calibration.
  - Calibrating alpha to the binned data instead makes today's value improbable. The joint probability is
    below 1e-3 for every alpha.
  - Later DESI, Euclid, Rubin and Roman data re-measure the same realization over the same redshifts. Their
    outcome is predictable from DR2 and the SN samples, so they are **not an independent test**.
- **Caveats.** Linear-Gaussian. No termination, BBN or CMB (each would only lower these probabilities). CPL
  central values stand in for a likelihood.
  - **Owed:** a likelihood-level check on binned rho_DE(z) posteriors.

### 2.3 The other tracks' consequences

| Consequence | Written before exposure? | Independent of exposed data? | Assessment |
|---|---|---|---|
| Barrow-Shaw Omega_k0 = -0.0055 (2 zeta_b) (`others` §8.1) | yes (2011, before Planck 2013 and DESI) | separate observable, but now exposed | the only genuinely ex-ante numerical consequence found; retrospectively disfavoured at 3.3-7.9 sigma depending on dataset [C] (§3.5) |
| Sorkin's 1987-97 "Lambda ~ 1e-120" (`unimodular` §4.2) | yes (pre-1998) | coincidence form only; sign open | low-information; correctly labelled |
| KP linear trigger: w0/wa, crunch (`sequestering` §5) | slope chosen as m^3 ~ M H0^2 | DESI quadrant exposed; crunch untestable | the track labels it retrospective; correct |
| KP "future collapse, future-testable in principle" (`obstruction` R4 row) | — | the collapse is 15-600 Gyr ahead (`sequestering` table); only precursors (w(z), Omega_k) are testable, and they are exposed | overstated (§3.7) |
| Everpresent "z > 2 tomography, CMB-S4" (`unimodular` §4.6) | — | DESI Lyα BAO (z = 2.33) and Planck/ACT recombination constraints are exposed and already adverse | improvements of exposed, already-failing comparisons |
| Li HDE w(z); Hubble-cutoff HDE; well-tempering (`others`) | c and the design use exposed data (GW170817 as a design input) | exposed | correctly labelled retrospective |

**Answer.** Barrow-Shaw is the only consequence written before the data that now test it, and it is
disfavoured. Nothing in the packet is held out.

## 3. Question 3: overclaims and errors in the stated conclusions

### 3.1 WDW, Lambda > 0: "no normalizable spectrum" is false [C, D]

`obstruction` §2.2 (lines 196-198), §0 item 2 and finding F5 state that for lam > 0 "there is no normalizable
spectrum and no quantization at all". Its check C8 passes if the lowest Dirichlet-box level moves by more than
1 between boxes 12 and 16. That cannot tell a continuum from a limit-circle endpoint.

**Tested here (A1-A2)** for -psi''/2 + (x^2/2 - lam x^4/2) psi = E psi with lam = 0.05:
- Both independent solutions (psi(0) = 0 and psi'(0) = 0) are square-integrable at large x for E = 1, 3 and
  -20.
- The norm increments follow the x^-2 tail: increment ratio 0.660 against 0.667 predicted.
- The extrapolated N(∞)/N(30) is 1.03-1.18.
- For contrast, lam = -0.1 and lam = 0 at a non-eigenvalue E give N(8)/N(4) = 2.6e44 and 6.8e18.

**Consequence [D/Q].**
- The endpoint is Weyl limit-circle. The operator -d²/dx² - x^4 is the textbook example of a Hamiltonian that
  is not essentially self-adjoint.
- Every self-adjoint extension, i.e. a supplied boundary condition at large a, has a discrete spectrum that is
  unbounded below. The Dirichlet box is one such extension, which is why its levels move with the box.
- So for Lambda > 0, normalizability *does* quantize E_r = E_n(Lambda; theta), but only after an extension
  parameter theta is supplied (compare Gielen 2026, cited at abstract level in `unimodular`).
- The track's conclusion (Lambda not fixed by closure) survives, and is sharpened: one more supplied datum.
  The stated mathematical claim must be corrected.

### 3.2 "R1 escapes the obstruction; Q1 and Q3 pass" [A]

`grut_native` §0 item 1, §2.5, §7 statement 1 and the §9 row R1 say R1 "escapes Codex's obstruction for the
gravitating term" and that "Q1 and Q3 pass".
- R1 decouples constant matter terms (a clock phase) and leaves the gravitating Lambda_eff a free state label.
  Q2 fails, so Q3 ("is the claimed value stable") has no value to test.
- "Stable if the state is specified by its h0 support" (B4-B5) is the Carballo-Rubio side of the dispute that
  `obstruction` §4.4 records as contested.
- The siblings say "relocation, not refutation". This should read "decouples constant terms; relocates the
  value to the state; Q2 fails; Q3 not applicable".
- The packet should define "escape" once, as Q1's wording ("constrains the gravitational vacuum term"), and not
  use it for "decouples" or for "evades Weinberg's theorem".

### 3.3 Two different "minimal escape structure" statements [A]

`obstruction` Claim C (G/B/S/X plus R and P) and `grut_native` §7 statement 3 do not match:
- Claim C's R item asserts that "every available scale is cosmological".
- `grut_native` instead lists four supplied elements:
  - an O(1e120) dimensionless input;
  - an exponential suppression e^-276;
  - an exponent-1/2 fluctuation law with a stated mean;
  - (implicitly) any of these.

Claim C's "complete relative to H1-H4 by construction" stands only if:
- B is read to include a Planck-scale inhomogeneous metric (WZU, Carlip) and rolling relaxation (Abbott, GKR);
- R is read as "some small dimensionless input".

That input can be a cosmological scale (the coincidence form), an O(100) exponent (WZU: 140 [C]), or a
matter-sector ratio (AHKM). The last two need their own stability argument. The synthesis should state one
harmonized list, with credits.

### 3.4 Smaller labelling and wording issues

- **`grut_native` §7 statement 4** ("closure ... cannot fix both", labelled [D]) is not derived. Classically V4
  and Lambda_eff are one-to-one at fixed matter content (its own G2), so fixing one fixes the other. Relabel
  [A] and reword: "closure supplies at most one as a datum; neither is selected by closure".
- **`grut_native` §0 item 6** quotes "P(our history | alpha) <= 1e-4" without the start-epoch dependence
  (0.016 for a start at a = 1e-3, from `unimodular` F_start_epoch_scan).
- **`others` §0 item 2**, "No class passes Q2 and Q3 together", conflicts with its own scoring of
  Baum-Hawking (fixes 0+, stable, wrong). It should read "no class fixes a nonzero residual compatible with
  observation while staying stable".
- **`others` §0 item 3**, "only Barrow-Shaw makes an independent numerical prediction", conflicts with its Li
  HDE row (c = 1 gives a definite w(z)). It should read "the only ex-ante independent prediction".
- **`obstruction` §1.5** (Casimir lemma): standard. Credit, for example, Jaffe, PRD 72, 021301 (2005),
  hep-th/0503158 [title verified], and Volovik 2006, who argues the opposite reading.

### 3.5 Barrow-Shaw tension depends on the dataset [C]

Omega_k0 prediction (-0.0055 printed; -0.0062 to -0.0064 with current inputs) against exposed Omega_K
summaries:

| Exposed summary | Tension |
|---|---|
| DESI DR2+CMB, LambdaCDM+Omega_K (0.0023 ± 0.0011) | 7.1-7.9 sigma |
| DESI DR2+CMB+Pantheon+, w0wa+Omega_K (0.0011 ± 0.0013) | 5.1-5.7 sigma |
| Planck 2018 TT,TE,EE+lowE+lensing+BAO (0.0007 ± 0.0019, recall) | 3.3-3.7 sigma |

`others` §8.1 and `grut_native` §4.4 quote only "7.7 sigma". Quote the range, with datasets.

### 3.6 R1c is "compatible with", not "a reason for" [A]

GRUT's access structure (FORMAL_CORE DEF 1.2) is defined on a supplied algebra with no metric. Identifying the
accessible clock with the 4-volume of the observer's past light cone is an added identification. The
causal-past choice is also Sorkin's own causal-set heuristic.

### 3.7 KP collapse "future-testable in principle" (`obstruction` R4 row)

Turnaround lies 15-600 Gyr ahead (`sequestering` §5 table). The residual depends on an unfrozen collapse
sector and on ln(rho_cut). Only exposed precursors are testable now. Write "not testable without a frozen
collapse sector; precursors exposed".

### 3.8 Consistency with the user's public statement

The user wrote "No further wall or population simulations were launched". The `unimodular` everpresent-Lambda
runs (180,000 realizations, plus start-epoch and grid scans) are stochastic ensembles of a published
dark-energy toy model, not Candidate 1 wall or population runs.
- The synthesis should describe them in those words, so that the public statement stays accurate.
- It should also restate the user's correction: Candidate 1 is not independently test-ready, and universal
  equivalence to cold dark matter is unproved. No track contradicts it.

### 3.9 What is sound

- `obstruction` §1 (fixed-background invariance) and §3 (local covariant form; Wald axiom iv fixes only the
  split).
- `obstruction` §4 and `unimodular` §§2-3 (relocation to state selection).
- `sequestering` §4 (blindness of T - <T>g/4) and §5, with its own limits.
- `unimodular` §3.4 (on-shell-Lagrangian ambiguity) and §4 numbers. I reproduce 0.1229, 0.19 and 0.72.
- The repeated statement that no route is GRUT-native.

These are correctly scoped and credited.

## 4. Verdict on "best route"

There is no best route in the gate's sense. Every surveyed route fails at least one of the three requirements:
fixing the total, stable without retuning, and an independent test.

The four tracks' "best" labels answer different questions:
- `obstruction`: sequestering is the "strongest existing class" (Q1, Q3).
- `grut_native`: R1 + R1b is "most promising" (chosen in advance; already strained).
- `unimodular`: everpresent Lambda is "best-developed" (Q2 scale only).
- `others`: Barrow-Shaw is "the only independent prediction" (Q4).

The synthesis should report per-criterion leaders together with their failures:
- **Q1/Q3: Omnia Sequestra.** Matter and graviton loops are sequestered; the residual is free and two-rooted.
- **Q2, scale only: everpresent Lambda.** The mean is assumed and alpha is fitted. On exposed data it is at or
  beyond its own pre-registered failure level (§2.2).
- **Q4, ex-ante specificity: Barrow-Shaw.** Disfavoured at 3.3-7.9 sigma.

The packet's conclusion should be **an identified obstruction, stated with its scope**: the one Codex proposed
and the user restated. It needs the harmonized escape list (§3.3) and the added credits (§1.2-1.3). It should
not name a GRUT route or an existing route as the best candidate. Before any later claim that "closure fixes
the vacuum energy", q-theory must be scored on the gate.

## 5. Missing items (owed before a synthesis)

1. A q-theory row (Klinkhamer-Volovik 2008a, 2008b; Volovik 2006) on the four questions, read beyond the
   abstract.
2. Rows for WZU 2017, Carlip 2019, Abbott 1985 and GKR 2019.
3. Credits for the escape taxonomy (Weinberg 1989; Padilla 2015; Burgess 2013).
4. Before calling any statistic held out:
   - a numeric freeze file (alpha, V choice, start epoch, recipe, model, last-scattering suppression, binning,
     likelihood);
   - committed to git **before** computation;
   - for a statistic not already constrained by DESI DR2 and the SN samples.
5. A likelihood-level version of §2.2 using published binned rho_DE(z) reconstructions.
6. A corrected WDW statement with an explicit self-adjoint-extension parameter, and a replacement for check C8.
7. One synthesis document that defines "escape", harmonizes the minimal-escape list and states the
   per-criterion ranking.
8. Primary texts still unread for load-bearing points: HT 1989, Unruh 1989 (abstracts only), Wald 1977
   (secondary only) and Weinberg 1989 equations (translation only).

## 6. Checks (`code/review_checks.py`, 17/17 pass)

- **A1-A2:** WDW limit-circle at Lambda > 0; limit-point contrast at Lambda <= 0.
- **B0-B4:** past-cone V and correlation reproduce the `unimodular` numbers. The pre-registered bin
  statistic is evaluated against exposed central values.
- **C1-C3:**
  - rho_Lambda^(1/4) = 2.240 meV;
  - K1 p = 0.5006;
  - the electron one-degree-of-freedom threshold is 4.29e30 rho_Lambda.
  The Barrow-Shaw tension range, the WZU exponent (140) and the AHKM scale (0.41 meV) are recorded without
  pass/fail.

Run: `python3 -I code/review_checks.py results/review_checks.json`.

## 7. Sources and read levels (new in this track; others as recorded in the sibling tracks)

| Source | Identifier | Read level |
|---|---|---|
| Volovik, "Vacuum energy: myths and reality" | IJMPD 15, 1987 (2006), gr-qc/0604062 | abstract (arXiv API) |
| Klinkhamer, Volovik, "Self-tuning vacuum variable and cosmological constant" | PRD 77, 085015 (2008), arXiv:0711.3170 | abstract |
| Klinkhamer, Volovik, "Dynamic vacuum variable and equilibrium approach in cosmology" | PRD 78, 063528 (2008), arXiv:0806.2805 | abstract |
| Wang, Zhu, Unruh | PRD 95, 103504 (2017), arXiv:1703.00543 | abstract |
| Carlip, "Hiding the cosmological constant" | PRL 123, 131302 (2019), arXiv:1809.08277 | abstract |
| Graham, Kaplan, Rajendran, "Relaxation of the cosmological constant" | PRD 100, 015048 (2019), arXiv:1902.06793 | abstract |
| Abbott, "A mechanism for reducing the value of the cosmological constant" | PLB 150, 427 (1985) | recall |
| Arkani-Hamed, Hall, Kolda, Murayama | PRL 85, 4434 (2000), astro-ph/0005111 | abstract |
| Burgess, "The cosmological constant problem: why it's hard to get dark energy from micro-physics" | arXiv:1309.4133 | title/listing |
| Jaffe, "The Casimir effect and the quantum vacuum" | PRD 72, 021301 (2005), hep-th/0503158 | title/listing |
| Planck 2018 VI, Omega_K with lensing + BAO | A&A 641, A6 (2020) | recall (0.0007 ± 0.0019) |
| Reed-Simon II, -d²/dx² - x^4 not essentially self-adjoint | textbook | recall; checked numerically (A1) |

Abstract downloads are kept in the scratchpad (`review_scope_dl/`). No paper text is in the repository.
