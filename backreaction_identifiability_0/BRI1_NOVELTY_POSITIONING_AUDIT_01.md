# BRI1 NOVELTY / POSITIONING AUDIT 01 — V2
## Executed under `BRI1_PUBLICATION_VERIFICATION_CHARTER_01.md` (frozen `8d765ea`), V2 lane only
## Hostile prior-art search; the goal is to kill the novelty claim, not rescue it

**Frozen candidate contribution (not broadened):**

> For a finite reciprocal anharmonic environment, interventional reduced-force
> laws across a shared clamp family can escape every single shared causal
> signed-affine exogenous forcing representation through a protocol-dependent
> standardized response-law shape invariant; for frozen X1 the witness is a
> third-cumulant/skewness effect of order 1/N_B, and it disappears in the
> reservoir limit.

**Presumed known (per charter §1):** anharmonic baths produce non-Gaussian
noise; finite baths have non-Gaussian corrections; higher cumulants occur;
GLEs arise from oscillator baths; Gaussian reservoir limits occur; nonlinear
response can change noise statistics.

---

## 2. Search executed (threat order A → B → C)

Approximately 30 distinct web searches with heavy synonym variation, plus
direct inspection of the closest candidate papers' full text / abstracts via
arXiv. Search strings varied along every axis of the charter's synonym
firewall (§3): bath/reservoir/environment/thermostat; back-reaction/feedback/
reciprocal coupling; clamp/intervention/forcing protocol/controlled
trajectory; exogenous/random/latent/stochastic drive; affine modulation/
location-scale/multiplicative noise; skewness/third cumulant/higher-order
statistic/non-Gaussian shape; identifiability/distinguishability/observational
equivalence.

---

## 4–5. Nearest prior-art table (papers actually inspected)

### Paper 1 (the strongest physical-domain threat found)

| Field | Content |
|---|---|
| 1. Citation | *Analysis and Simulation of Generalized Langevin Equations with Non-Gaussian Orthogonal Forces* (arXiv:2505.15665) |
| 2. Model class | Butane dihedral-angle dynamics in explicit solvent; classical MD; GLE coarse-graining via Mori and Zwanzig projections |
| 3. Coupling | Reciprocal in the underlying MD (real solvent), but this is not studied as a variable |
| 4. Bath size | Effectively infinite (explicit solvent); no finite-bath axis is studied |
| 5. Intervention family | **No interventional protocol family on a fixed system.** The comparison axis is *projection formalism* (Mori vs. Zwanzig vs. generalized) applied to the same system; and a second physical setup (butane with two frozen carbon atoms) which is a different molecule/simulation, not a protocol applied to one fixed environment |
| 6. Observable | Orthogonal-force distribution (non-Gaussian shape), two-point correlations, mean first-passage times |
| 7. Competitor model class | Markovian embedding with Gaussian noise (as the thing that fails numerically); **no shared affine/location-scale exogenous competitor class is formulated** |
| 8. Identifiability theorem? | **No.** Zero occurrences of "shared", "exogenous", "affine" in the full text; no theorem about whether the orthogonal force can be represented by a common modulated noise across multiple interventions |
| 9. Explicit finite-size rate? | No (no 1/N axis; "thermodynamic limit" count = 0) |
| 10. Overlap with BRI1-X1 | Overlaps on: anharmonic environment → non-Gaussian orthogonal force (ingredient B/C); that Gaussian approximations can fail for kinetic observables; explicit-solvent classically-exact framework |
| 11. Exact difference | BRI1-X1 is a **theorem** about a *shared* exogenous competitor across an *interventional protocol family* on a *fixed* environment, with an explicit *finite-bath escape rate* and *reservoir-limit disappearance*. Paper 1 has: no protocol family (projection-formalism axis, not intervention axis), no shared-exogenous competitor class, no finite-bath axis, no theorem, no rate. **Not an equivalent result; the closest physical-domain paper found.** |

### Paper 2 (quantum non-Gaussian bath review-type)

| Field | Content |
|---|---|
| 1. Citation | Papers studying non-Gaussian baths in open quantum systems (e.g. quantum Brownian motion with anharmonic baths, higher-order influence functionals; located via search but not individually deep-read) |
| 2–7 | Quantum system-bath models; reciprocal coupling typical; baths often finite or truncated; protocols usually single; observable: reduced dynamics/correlations |
| 8. Identifiability theorem? | Not found in inspected abstracts; focus is simulation accuracy / correlation corrections, not shared-exogenous distinguishability |
| 9. Explicit finite-size rate? | Not found in inspected abstracts |
| 10–11. Overlap/difference | Ingredient overlap only (non-Gaussian corrections exist, finite baths matter). No shared-affine-competitor theorem found. **Weaker than Paper 1 as a threat.** |

### Paper 3 (statistical/causal domain — closest conceptual analogue)

| Field | Content |
|---|---|
| 1. Citation | *Towards identifiability of interventional stochastic differential equations* (arXiv:2505.15987) and adjacent nonlinear-ICA / interventional-SDE literature |
| 2. Model class | Abstract SDEs / structural models; not physical system-bath |
| 3–5. Coupling/bath/protocols | N/A in physics sense; "interventions" = distribution shifts |
| 6. Observable | Drift/diffusion parameters, latent variables |
| 7. Competitor model class | Observational equivalence of SDE parameterizations under interventions |
| 8. Identifiability theorem? | Yes (in the abstract statistical sense: what interventions render identifiable) |
| 9. Explicit finite-size rate? | No |
| 10. Overlap with BRI1-X1 | Conceptual overlap at the level of "interventions can distinguish components of a stochastic model" — but no shared-affine-quotient class, no physical back-reaction structure, no shape-invariant witness, no finite-size escape rate |
| 11. Exact difference | Statistical identifiability results concern parameter identifiability of a *specified* model family; BRI1-X1 is a *non-representability theorem* against a *defined competitor class* (all shared causal affine modulations), which is a different mathematical object. **Not equivalent; positioning-relevant.** |

Also inspected (quickly): non-Gaussian orthogonal-force papers via Mori
projection (Zwanzig/Mori literature); active-matter/bath papers measuring
higher-order environmental statistics; microrheology papers distinguishing
bath types from probe statistics; quantum process-tensor identifiability
papers. None found formulating a shared-affine-quotient escape theorem.

---

## 6. Ingredient provenance matrix

| Ingredient | Status | Notes |
|---|---|---|
| A. finite reciprocal anharmonic bath | KNOWN (as a physical model) | Studied in MD / classical GLE literature |
| B. non-Gaussian environmental force | KNOWN | Presumed known per charter |
| C. third/higher cumulants | KNOWN | Presumed known |
| D. protocol-dependent response | KNOWN-IN-NEARBY-FORM | Response depends on drive; nonlinear response theory |
| E. multi-protocol/interventional comparison | KNOWN-IN-NEARBY-FORM | Standard in causal/statistical literature; in physics, less formalized |
| F. one shared exogenous competitor across protocols | KNOWN-IN-NEARBY-FORM | Implicit in "Markovian embedding"/"exogenous noise" framings; not formalized as a shared-quotient class in the physical literature found |
| G. affine/location-scale/sign quotient | NOT-FOUND (as a defined physical competitor class) | This quotient formalization is what makes BRI1-X1 a theorem |
| H. shape invariant surviving affine quotient | NOT-FOUND (as a formalized witness class) | "Standardized force law is invariant" exists implicitly in the projection-formalism comparison (Paper 1 compares orthogonal-force *shapes* across formalisms) but not as a defined quotient-invariant witness |
| I. theorem-level failure of the shared representation | NOT-FOUND | No paper found proving non-representability against a shared affine class |
| J. explicit O(1/N_B) escape | NOT-FOUND (as an escape rate) | 1/N corrections to non-Gaussianity are structurally expected (CLT); but its role as the *escape rate from a defined competitor class* is not found |
| K. Gaussian/effective-harmonic reservoir recovery | KNOWN | Presumed known per charter |

**The candidate novelty lives in G+H+I+J as a combined statement.**

## 7. Makri / effective-harmonic comparison

The effective-harmonic-bath / Gaussian influence-functional literature
(Makri and others) says: a bath with weak/nonlinear coupling can be
*approximated* by an effective harmonic bath at the level of two-time
correlation functions, under weak-coupling or thermodynamic-limit
assumptions. What it does **not** do (as far as inspected):

- define an intervention/protocol family and ask whether the *same*
  effective harmonic representation works across all of them;
- quantify the finite-bath departure from the effective-harmonic
  description as an escape rate from a *defined* competitor class;
- formulate the competitor as a shared affine/location-scale exogenous
  process.

So: the **reservoir-limit side** of BRI1-X1 (K: Gaussian recovery) is
consistent with known physics. The **finite-bath identifiability theorem**
(I+J against G) is not found there.

## 8. Non-Gaussian-bath comparison

Papers where nonlinear/anharmonic baths produce non-Gaussian noise (Paper 1
is the strongest found): they **derive** nonzero higher cumulants and
**demonstrate** that Gaussian approximations fail numerically for kinetic
observables. None found proves the statement "no single shared affine
exogenous process can reproduce the response laws across multiple
interventions." Paper 1 in particular does not formulate that question.

## 9. Identifiability comparison

Statistical literature (Paper 3, nonlinear ICA, interventional SDEs) proves
identifiability of *parameters* or *latent variables* under interventions.
BRI1-X1 is a different mathematical object: a *non-representability theorem*
against a *defined competitor class* (shared causal affine modulations), not
a parameter-identifiability result. No abstract theorem directly subsuming
BRI1-X1 was found.

## 10. Finite-size-rate comparison

The 1/N_B scaling of non-Gaussian cumulants is **structurally expected**
(central-limit behavior for sums of i.i.d. bath contributions; standard
finite-size corrections). It is **not claimed as novel**. Its novelty (if
any) is its role as the *quantitative escape rate from the defined shared
affine quotient* — a framing not found in the literature.

---

## 11. Earliest-priority rule

No older paper was found proving an essentially equivalent result.

## 13. One-paragraph delta

> Relative to the nearest prior work, BRI1-X1 contributes a **theorem**,
> not a simulation finding: for a finite reciprocal anharmonic environment,
> it proves that the *standardized shape* of the environmental force — a
> quantity invariant under any shared deterministic rescaling and sign
> modulation of a common random process — differs between two intervention
> protocols applied to the same environment, so that no single shared
> affine-modulated random process can reproduce the measured force laws
> across the protocol family. It further proves that this distinguishing
> signature is of order 1/N_B in the number of bath degrees of freedom and
> vanishes in the thermodynamic (large-bath) limit, where the standard
> Gaussian linear-response description is recovered. Prior work has
> established that anharmonic environments produce non-Gaussian forces and
> that Gaussian approximations can fail numerically; it has not, to the
> best of this audit's search, formulated the shared-affine competitor
> class as a defined object, proven its failure across a protocol family,
> or given the finite-bath escape rate.

## 14. Adversarial referee section

**Referee objection 1: "This is already known because anharmonic baths
produce non-Gaussian noise, and it is well known that Gaussian
approximations fail for coarse-grained dynamics."**

Response: Correct as stated — and BRI1-X1 does not claim otherwise
(ingredient B/C/K are marked KNOWN). The theorem's claim is different in
kind: it defines a *competitor class* (all shared causal affine modulations
of one exogenous process — which includes non-Gaussian ξ, arbitrary memory
M, and arbitrary sign/scale modulation G), and proves that *no member of
that class* can reproduce the force laws across a protocol family, using a
witness (standardized shape) that is invariant under the entire class. This
is a distinguishability statement about a *quotient*, not a statement about
Gaussianity. The nearest paper (arXiv:2505.15665) compares *formalisms* on
the same system, not *protocols* against a *shared competitor*.

**Referee objection 2: "This is merely a change of language, not a new
result — 'standardized force law differs across protocols' is just saying
the orthogonal force is non-Gaussian in a protocol-dependent way."**

Response (partial concession, partial rebuttal): The *observation* that
non-Gaussianity is protocol-dependent is indeed implicit in nonlinear
response physics. What is not implicit, and is the theorem's content, is:
(a) the formalization of the competitor as a *shared* class across
protocols (so that a *per-protocol* fit is ruled out by construction);
(b) the proof that a single *shape invariant* — not a fitted parameter —
decides membership; (c) the explicit 1/N_B rate at which the escape
disappears, connecting the identifiability statement to the known
Gaussian-reservoir limit. If the referee still holds that (a)–(c) are
routine, the correct disposition would be KNOWN-RESULT-NEW-FRAMING; this
audit's judgment is that (a) is a genuine formalization step, so
NOVEL-THEOREM-SHAPE is defensible but the margin is not large.

## 15. V2 disposition

**NOVEL-THEOREM-SHAPE** — with the qualifier that the margin over
"KNOWN-RESULT-NEW-FRAMING" rests on the formalization of the shared affine
quotient as a defined competitor class and the theorem-level escape proof,
not on any ingredient being new.

## 16. Search sufficiency

Literature access was sufficient for the strongest physical-domain
candidate (full text inspected) and adequate for statistical-domain
candidates (abstracts/related work). A deeper search of the quantum
open-system literature (higher-order influence functionals, process
tensors) was lighter; no shared-affine-quotient theorem was found there
either, but the margin is thinner.
