# GRUT CONJECTURE MODE CHARTER 01
## with Annex A (refiled historical cards) and Annex B (Card #1 rerun requirements, NOT EXECUTABLE)

**Status:** **ACCEPTED IN SUBSTANCE** by owner/reviewer ruling on the amended draft. Three editorial repairs have been applied:
- ER-1: the Prior-Exposure Disclosure is made permanent and general;
- ER-2: owner-lock authority is clarified;
- ER-3: the expected-verdict sentence is removed from the Card #1 requirements.

The governance state is **awaiting owner review / freeze**.

**Home:** `grut-program-governance-1`, branched from `grut-selector-screen-1-frozen @ a1194546c08d81c110a6ad8cc7583fddeedbebd1`. This file is **program law**.

**Authority split:**
- Part I (including §7) is authoritative here.
- Annex A records the owner-declared historical card statuses, as intake for the future execution branch.
- Annex B states the minimum requirements binding on any future Card #1 spec.
- Card ledgers, registries, individual card specs, data provenance, code and results belong to a **separate descendant branch** (e.g. `grut-conjecture-mode-1`). That branch is **not** created yet and opens only after this governance state is reviewed and frozen.

**Execution status:** **CARD #1 IS NOT EXECUTABLE.**
- No DESI product may be accessed, and no network preflight may be run.
- No threshold is frozen in this document. Every numerical kill boundary is an **OWNER-LOCK slot**, left empty on purpose.
- The exact Card #1 response law is an empty C0 field.

**Amendment 01 summary** (all accepted by the owner):

| # | Amendment | Where |
|---|---|---|
| A1 | New status **CONJECTURE-GENERAL** | §2, §2A |
| A2 | New C0 field **C0.2 Origin / GRUT rationale**; later fields renumbered C0.3–C0.7 | §4 C0, §6, C-F9 |
| A3 | Card #1 primary claim = **dissipative branch only (ε < 0)**; ε > 0 is a labeled mathematical control | Annex A Card #1, Annex B §2.2, §2.6 |
| A4 | **Owner-lock gate:** C0 spec commit → OWNER THRESHOLD RULING → data access | §4A, Annex B §2.7, §5, §10 |
| A5 | **D3 primary**; CPL is a diagnostic only; explicit downgrade if D3 cannot be done | Annex B §4 |
| A6 | **Network-access preflight** → CARD-01-ACCESS-BLOCKED, with no improvised substitute | §4B, Annex B §3.0 |
| A7 | **Result states** for nested-null cards (EXPLANATORY-SURVIVES / EXPLANATORY-KILLED / MODEL-KILLED / …) | §4C, Annex B §5, §8 |
| A8 | **Prior-Exposure Disclosure**: permanent; applies to all human and AI participants (ER-1) | §4D |
| ER-2 | Owner-lock **authority roles**: executor proposes, reviewer audits/recommends, scientific owner approves | §4A, C-F10 |
| ER-3 | An executable card spec contains no forecast of its own verdict | §4D, Annex B §5 |

---

# PART I — GRUT CONJECTURE MODE CHARTER 01

## 1. Purpose

GRUT currently has a strong audit and derivation architecture. It has no formal channel for introducing a new physical law that is deliberately **postulated rather than derived**.

This charter creates that channel.

The core rule:

> **A postulate may be unearned. It may never be hidden.**

Conjecture Mode does not relax the scientific firewall. It changes the order of operations:

1. Introduce a specific postulate openly, **originating in GRUT's own structure** (C0.2).
2. Price its new information.
3. Test its immediate consistency.
4. **Have the owner lock its kill thresholds before any data is opened.**
5. Force it into data contact.
6. Apply the full GRUT audit machinery only if it survives.

This separates hypothesis generation from hypothesis certification.

---

## 2. Epistemic statuses

Add the statuses **POSTULATED** and **CONJECTURE-GENERAL** to the existing scientific vocabulary.

The relevant statuses are:

- **SUPPLIED**: required as an input to the framework.
- **POSTULATED**: deliberately introduced as a new law or value for testing; not derived.
- **CONJECTURE-GENERAL**: a bold generalization across cases, used to steer what gets tested next (§2A). It is not a law and not a result.
- **DERIVED**: follows from declared premises.
- **CONSTRUCTED**: an explicit mathematical realization exists.
- **VERIFIED / REPRODUCED**: independently checked at the declared verification grade.
- **PREDICTED**: a preregistered observable consequence follows from a frozen postulate and projection.
- **VALIDATED**: survives external empirical or independent scientific review.

A POSTULATED claim must never be described as derived, selected, emergent, or uniquely required.

## 2A. CONJECTURE-GENERAL

Bold generalizations are **allowed and encouraged**. Example: *"long memory tends to smooth rapid transitions."*

A CONJECTURE-GENERAL entry is valid only if it carries all four fields:

| field | content |
|---|---|
| G1 Supporting cases | the specific cards, theorems or runs it generalizes from, cited by record ID |
| G2 Bold generalization | the general statement, in one sentence |
| G3 Current scope | the class in which the supporting cases actually live (models, regimes, approximations) |
| G4 Breaking test | a specific counterexample or test that would refute it |

Permissions and limits:

- **May** guide which card, campaign or test is proposed next. It may be cited in a card's C0.2 Origin field (as long as C-F9 is met).
- **May not** itself earn PREDICTED, DERIVED or VALIDATED status.
- **May not** support a class-level exclusion ("the class cannot …").
- **May not** be quoted outside the program without its G3 scope attached.
- If its G4 test is run and fails, it is banked as **CONJECTURE-GENERAL-REFUTED** (C-F8 applies).

CONJECTURE-GENERAL entries live in `conjecture_mode/CONJECTURE_GENERAL_REGISTER.md`.

---

## 3. Two legal doors for future scientific campaigns

### Door D — Derivation / Selector Mode

A derivation campaign may open only if it:

1. names an exact entry in the frozen residual ledger that it intends to eliminate or render downstream;
2. names a mechanism capable of distinguishing at least two currently admissible residual members;
3. plausibly passes the selector-screen S1 information-elimination test;
4. plausibly passes S2 witness discrimination before any calculation.

If these conditions fail: **NO DERIVATION CAMPAIGN.**

The frozen selector-screen gate applies to Door D.

### Door C — Conjecture Mode

A conjecture campaign does **not** need to eliminate the residual entry. Instead it must:

1. name the exact residual entry it instantiates;
2. **show that the postulate originates in GRUT's own structure (C0.2), not in the shape of a dataset anomaly;**
3. introduce one explicit physical postulate;
4. state its information price;
5. freeze the map from the postulate to an observable;
6. name a real dataset or experiment;
7. preregister a kill condition, **with its numerical boundaries fixed by owner ruling before data access (§4A).**

Conjecture Mode is permitted precisely because the new law is openly labeled POSTULATED.

---

## 4. Conjecture pipeline

### Stage C0 — Conjecture Card

A card is invalid unless all **seven** fields are frozen before any computation against target data.

#### C0.1 Residual target
The exact residual-ledger entry being instantiated.

#### C0.2 Origin / GRUT rationale *(new, Amendment A2)*
State which piece of GRUT's own structure the postulate comes from. It must be at least one of:

- (a) a **SUPPLIED primitive**;
- (b) an **unresolved or declared fork** in the frozen record, cited by ID;
- (c) a **residual coupling**: a frozen-ledger freedom that the postulate fixes;
- (d) the **core responsiveness / memory principle**, together with the specific way the postulate instantiates it.

Rules:

- "The data anomaly has this shape" is **not** a valid origin and cannot create a GRUT card.
- External physics may be used to **instantiate or test** the postulate (a standard fluid description, a known likelihood). It cannot serve as the origin.
- **Preference rule:** when choosing between cards, prefer targets where GRUT's distinctive structure is load-bearing. That means a target where removing the GRUT ingredient changes the prediction. A card that any standard dark-energy parameterization reproduces with the same information price is low priority.
- The field must also say **which features of the postulate are not dictated by the origin**. Those features count toward C0.4 (information price).

#### C0.3 Postulate
One explicit law, equation, value, kernel, boundary rule or parameter relation. It must be precise enough to calculate with.

#### C0.4 Information price
Count every new:

- scalar parameter;
- dimensional scale;
- dimensionless constant;
- function;
- kernel shape;
- boundary or initial datum;
- branch or sign choice;
- projection convention that is physically load-bearing.

A fitted parameter is not free merely because it is conventional.

A branch choice that the origin (C0.2) fixes on stated grounds is recorded as **"fixed by origin"**, with the grounds given. If those grounds are later withdrawn, the choice becomes a priced bit.

#### C0.5 Frozen projection
Freeze every convention needed to map the postulate to the observable, where relevant:

- background cosmology;
- redshift or time interval;
- observable definition;
- fitting basis;
- weighting measure;
- likelihood;
- nuisance treatment;
- normalization;
- parameter priors;
- interpolation;
- numerical tolerance.

A postulate does not predict a number until its projection is frozen.

#### C0.6 Observable + dataset
Name the observable and the exact public dataset or likelihood product.

#### C0.7 Kill condition
State the empirical or mathematical result that kills the card, using the result-state vocabulary of §4C. The **form** of the kill condition belongs in C0. Its **numerical boundaries** are fixed only by the owner ruling of §4A.

---

### §4A. Owner-lock gate *(Amendment A4)*

Each data-contact card follows a strict three-step order:

```
(1) C0 SPEC COMMIT        — all seven C0 fields, kill-condition FORM, threshold SLOTS empty
        ↓
(2) OWNER THRESHOLD RULING — owner fills every numerical slot; committed as its own file
        ↓
(3) NETWORK PREFLIGHT → DATA ACCESS
```

**Authority roles (ER-2):**

| role | who | may |
|---|---|---|
| **Executor** | Claude Code (or any agent running the card) | compute, tabulate and propose candidate thresholds from **pre-data considerations only** (the model's own manifold, numerical tolerance, generic statistical conventions, the information price). Each proposal is labeled `CLAUDE-PROPOSAL (non-binding)`. The executor **never** grades itself and cannot authorize any slot. |
| **Reviewer** | the owner's designated scientific reviewer | audit the proposals and recommend values. A recommendation is non-binding. |
| **Scientific owner** | the program owner / user | **explicitly approve** every final numerical value (S1–S8 for Card #1; the corresponding slots for any card). Only the owner's approval authorizes a threshold. |

- **The executor may not** choose the final value of any of the following:
  - `epsilon_null`;
  - any probability or credible-region level;
  - any Δχ² / likelihood-ratio cut;
  - any "non-negligible manifold" threshold;
  - any other kill boundary.
- A threshold is valid **only if** it appears in a committed file `…/CARD_NN_OWNER_THRESHOLD_RULING.md`, and that commit is an ancestor of every commit touching data products.
- The ruling file records:
  - the owner's approved values, with an explicit statement of owner approval;
  - for each value, whether it adopts, modifies or rejects the executor's proposal and the reviewer's recommendation;
  - the commit hash of the C0 spec it rules on.
- The executor may transcribe and commit the ruling file. It may not author, alter or fill any value that lacks the owner's explicit approval.
- A data-access step whose commit history does not contain a valid owner ruling is **void**. Any result it produced is banked as **PROCEDURE-VOID**, not as evidence.

### §4B. Network-access preflight *(Amendment A6)*

Before any data work, run and log a reachability check against the **official** hosts and the exact products named in C0.6. The check uses HEAD / listing requests only and does not download likelihood or chain content.

- **All reachable** → proceed.
- **Any required official product unreachable** → the card's result is **ACCESS-BLOCKED**. Stop.
- **Forbidden substitutes** after a blocked preflight:
  - mirrors that are not official;
  - values read off plots;
  - numbers quoted in papers or from memory;
  - Gaussianized summaries;
  - reconstructed ellipses.
  Using any of them requires a separate, explicit owner ruling that downgrades the card to C2-R **before** it is used.

### §4C. Result states for cards that nest a null model *(Amendment A7)*

When the card's family contains the standard null model (for example, ε = 0 ⇒ ΛCDM), a non-rejection **is not support**. The data can always retreat to the null. The allowed result states are:

| state | meaning |
|---|---|
| **EXPLANATORY-SURVIVES** | A **non-null** region of the card's allowed (primary-branch) family remains viable at the owner-locked threshold, **and** it accounts for the targeted signal to the owner-locked degree. |
| **EXPLANATORY-KILLED** | Only the owner-locked null neighborhood (e.g. `|ε| < ε_null`) remains viable, **or** the allowed deformation cannot reach the targeted signal region. The card is dead **as an explanation of the targeted signal**. |
| **MODEL-KILLED** | Even the nested null and the whole allowed manifold conflict with the data at the owner-locked threshold. This is reserved for that much stronger case. |
| **INCONCLUSIVE** | Data cannot discriminate at the locked thresholds (e.g. the posterior is too broad). |
| **ACCESS-BLOCKED** | Required official products are unreachable (§4B). |
| **RECONNAISSANCE** | Only C2-R-grade contact was possible. Wording is limited to "NOT FOUND WITHIN THE TESTED APPROXIMATION". |
| **PROCEDURE-VOID** | The owner-lock order (§4A) was violated. |

Mandatory wording: an EXPLANATORY-SURVIVES result must be reported together with the null comparison (C-F7). It must never be summarized as "the data support Card NN" without stating Δ relative to the null and the information price (C-F4).

### §4D. Prior-Exposure Disclosure *(Amendment A8; permanent)*

A card that tests against an already-published dataset is **not a blind experiment**. Participants generally know the dataset's published headline results.

The Conjecture Mode firewall is therefore **not** blindness to the published headline. It is this:

> **the postulate, origin, branch, projection, thresholds and analysis procedure are all frozen and committed before the likelihood products are accessed and before the test is executed.**

Requirements:

- **All participants, human and AI** (owner, reviewer, executor and any sub-agent), disclose in the card spec what they knew about the target dataset's published results when the spec was written. This includes preferred regions, quoted significances and any prior reconnaissance contact.
- The disclosure is a **record of exposure**. It is not a forecast. An executable card spec must not forecast its own verdict (ER-3). Model-derived structural facts (e.g. which side of w = −1 a branch lies on) may be stated. The likelihood test produces the verdict.
- Threshold proposals and recommendations must not be tuned to the disclosed headline. The owner-lock gate (§4A) is the safeguard.
- A comparison against existing published data remains a retrospective test (Annex B §6). Only a preregistered future holdout can earn PREDICTED.

---

### Stage C1 — Minimal Firewall

Run before any fit to the target dataset.

Required checks:

1. Dimensions.
2. Conservation laws.
3. Well-posed evolution.
4. Initial and boundary data specified.
5. **Equivalence registry:** prove the card is not a reparameterization of an already-tested card at the declared scope.
6. **Gross consistency constraints:** e.g. early-universe behavior, positivity, singularities, causal structure, **energy conditions where a fluid reading is used**.
7. **Parameter identifiability:** the data must constrain what the card introduces.
8. **No post-selection:** confirm the functional form and kill condition were not chosen after seeing this card's residuals.
9. **Known-physics comparison:** identify the nearest standard model or class, and what, if anything, is distinctive.
10. Information-price audit.
11. **Origin audit (C0.2):** confirm the origin is admissible under (a)–(d) and that C-F9 holds.

Outcomes:

- **C1-PASS**
- **C1-REPAIR**
- **C1-EQUIVALENT**: bank the equivalence and do not run a duplicate card.
- **C1-KILL**

No target-data fit is allowed before C1 passes **and** the owner threshold ruling is committed.

---

### Stage C2 — Data Contact

Use the real likelihood or posterior products where publicly available.

#### C2-A — Full data contact
Uses the actual published likelihood, chain or data products appropriate to the model.

#### C2-R — Reconnaissance
Uses approximations such as:

- Gaussianized posterior ellipses;
- assumed correlations;
- proxy CMB anchors;
- coarse parameter grids;
- background-only substitutes for a perturbation-complete model.

Reconnaissance may:

- kill a card;
- expose a flaw;
- motivate a repair;
- rank candidates for a real run.

Reconnaissance may **not** bank a quantitative prediction or a class-level exclusion.

Required wording: **NOT FOUND WITHIN THE TESTED APPROXIMATION**, never **THE CLASS CANNOT…**.

---

### Stage C3 — Full Audit

Only a card that reaches **EXPLANATORY-SURVIVES at C2-A** triggers the existing full GRUT audit machinery.

C3 may include:

- hostile controls;
- theorem audit;
- derivation-vs-postulate accounting;
- comparison with standard model classes;
- an independent code path;
- external primary-source review;
- reproducibility checks;
- a preregistered follow-up prediction.

A card cannot become a GRUT prediction merely because it fits data.

---

## 5. Anti-curve-fitting rules

**C-F1 — Card budget.** At most **3 genuinely distinct cards per target** before a mandatory zoom-out. An equivalence discovered before fitting does not use up a new-model slot, but it is banked in the equivalence registry.

**C-F2 — Structural motivation.** A new card may not be introduced merely because the previous card missed a residual in parameter space. Its motivation must be one of:

- an independent physical principle;
- a structural failure of the previous model;
- a previously declared unresolved fork;
- a different exact residual entry.

**C-F3 — Equivalence registry.** Every proven collapse or reparameterization is permanently banked. A mathematically equivalent model may not return as a new card under different notation.

**C-F4 — Information-price dominance.** If two models perform comparably in statistical terms, the one with the lower information price is preferred. No claim of explanatory progress may ignore the information price.

**C-F5 — No same-data invention and validation.** A functional form invented after inspecting a dataset cannot use that same dataset as confirmatory evidence. It remains exploratory until tested on new or withheld data.

**C-F6 — No moving projection.** Once C0 is frozen, the projection cannot be changed after seeing the result without creating a new card or version.

**C-F7 — Null-model comparison.** Every empirical card must be compared against the appropriate null or standard model under the same likelihood treatment.

**C-F8 — Negative results are banked.** A killed postulate, a refuted CONJECTURE-GENERAL entry, and the reasons they died all remain part of the scientific record.

**C-F9 — No target-chasing *(new, A2)*.** A postulate whose functional form, sign or branch was chosen **to move a model toward a dataset's preferred region** fails C0.2, whatever origin is claimed. Any branch or sign choice must be justified from GRUT structure (C0.2), not from which side the data prefer.
- If the GRUT-justified branch predicts the side the data **disfavor**, that branch is still the one tested.
- The disfavored-side outcome is a legitimate, sharper kill and is banked as such.

**C-F10 — Owner-locked thresholds *(new, A4 / ER-2)*.** A kill boundary is valid only with the scientific owner's explicit approval in a committed ruling that predates data access. Values proposed by the executor or recommended by the reviewer are not valid on their own (§4A).

---

## 6. Conjecture card record format

| field | required content |
|---|---|
| Card ID | stable identifier |
| Target residual entry | exact frozen-ledger entry (C0.1) |
| **Origin / GRUT rationale** | C0.2 category (a)–(d), with citation |
| Status | POSTULATED |
| Physical law | exact equation (C0.3) |
| **Primary branch / controls** | the preregistered branch; any labeled mathematical controls |
| Information price | parameters, functions and choices (C0.4) |
| C0 projection | frozen conventions (C0.5) |
| Dataset | exact public release or product (C0.6) |
| Kill condition | preregistered form (C0.7) |
| **Owner threshold ruling** | file + commit hash; must predate data access |
| **Prior exposure** | headline knowledge at spec time (§4D) |
| **Network preflight** | log + verdict |
| C1 verdict | PASS / REPAIR / EQUIVALENT / KILL |
| C2 grade | FULL or RECONNAISSANCE |
| **Result state** | §4C vocabulary |
| Distinctive? | yes/no, and the comparator |
| Audit status | not opened / C3 |
| Reproduction status | explicit |
| Notes | scope only |

---

## 7. Governance interaction

This charter modifies the proposed program-wide campaign gate as follows:

> **No new DERIVATION/SELECTOR campaign may open unless it names the frozen residual entry it intends to eliminate or render downstream and a named mechanism plausibly capable of distinguishing at least two frozen-admissible members.**

> **No new CONJECTURE campaign may open unless it names the frozen residual entry it instantiates, an admissible GRUT origin (C0.2), an explicit postulate, its information price, a frozen observable projection, a dataset, and a preregistered kill condition whose numerical boundaries are fixed by a committed owner ruling before data access.**

The selector-screen terminal governs the first door, not the second. This is not permission for unrestricted parameter fitting.

*Placement:* this charter's §7 text is incorporated into `program_governance/PROGRAM_CAMPAIGN_GATE_01.md` on `grut-program-governance-1`, which is the authoritative home for program-wide campaign governance. It is binding once the owner reviews and freezes the governance state.

---

# ANNEX A — EXISTING DARK-ENERGY CARDS: REFILED STATUS

These are historical or reconnaissance cards. They are not retroactively upgraded. This annex is the owner-declared **intake** for the future execution branch's `CARD_LEDGER.md`. Once that ledger exists, it becomes the authoritative card record.

**Origin audit note (A2):** Cards #2–#4 were written before C0.2 existed. They are entered at their current grades with **C0.2: NOT AUDITED**. The owner has observed that Cards #2–#4 drifted toward DESI's preferred region. Any future rerun of #2–#4 must first pass the C0.2 origin audit and C-F9 as a C1 item.

## Card #1 — Horizon relaxor

**Postulate:** the inserted second relaxation scale is tied to the present Hubble scale,

`tau_2 = H_0^{-1}`,

with the previously used single-mode toy response.

**Primary branch (A3):** **dissipative branch, ε < 0 only**, fixed by the record's second-law interpretation of the relaxor (C0.2 / C0.4 "fixed by origin"). The ε > 0 branch is **not** part of the preregistered GRUT claim. It is retained only as a labeled **mathematical control**.

**Current grade:** C2-R — RECONNAISSANCE.

**Current safe statement:**

> Under the tested projection, the model produced a positive-slope one-parameter manifold through LambdaCDM and did not enter the current DESI+CMB+SN preferred `w0 > -1, wa < 0` quadrant except near LambdaCDM.

**Re-reading under A3** (bookkeeping only; nothing new computed): the reconnaissance manifold contained both branches. Restricted to the dissipative branch, the historical result concerns the branch on the w ≤ −1 side. The reconnaissance verdict is not upgraded or downgraded by this restriction.

**Not yet banked:** `wa/(1+w0) = 1.37`, because the ratio depends on the projection convention.

## Card #2 — Matter-responsive vacuum

**Current grade:** C1-REPAIR. **C0.2: NOT AUDITED.**

The prior CMB anchor was insufficient for an interacting model. Required repairs before a real C2 run:

- separate baryons from interacting CDM;
- a consistent background evolution from recombination;
- an explicit interaction four-vector / perturbation prescription `Q^mu`;
- a physical early-universe parameterization;
- a full comparison using proper cosmological likelihood machinery.

Do not bank the previous approximate chi-square values as class-level results.

## Card #3 — Linear geometric target

**Current grade:** C1-EQUIVALENT at narrow scope. **C0.2: NOT AUDITED.**

Banked equivalence:

> In flat homogeneous FRW with pressureless matter and a `w=-1` vacuum, and at the declared background level, targets linear in `rho_m`, `rho_v`, `H^2` and `dot H` collapse to the same linear-response family after coefficient redefinition.

Do not extend this equivalence to:

- radiation;
- curvature;
- perturbations;
- different interaction laws;
- nonlinear target functions.

## Card #4 — Long-memory variants of Card #2

**Current grade:** C2-R — RECONNAISSANCE, inheriting Card #2's C1 defect. **C0.2: NOT AUDITED.**

Safe statement:

> The particular tested exponential and power-law kernels did not improve the approximate fit relative to the tested single-pole implementation.

Do not generalize to arbitrary memory kernels.

*Possible CONJECTURE-GENERAL seed (not registered; for owner decision):* "long memory tends to smooth rapid transitions". To be registered it would need G1–G4, and G1 cannot rest on Card #4 alone, given the inherited C1 defect.

---

# ANNEX B — CARD #1 REAL-LIKELIHOOD RERUN REQUIREMENTS

> **NOT EXECUTABLE.** This annex defines the minimum requirements and the empty slots for a future Card #1 spec. It does not authorize data access. The card spec itself will live on the execution branch. Execution requires the following, in order:
> 1. owner review and freeze of this governance state;
> 2. creation of the descendant execution branch (e.g. `grut-conjecture-mode-1`) by separate authorization;
> 3. a committed C0 spec with every non-threshold field filled;
> 4. a committed OWNER THRESHOLD RULING;
> 5. the network preflight.

## 1. Objective

Upgrade Card #1 from reconnaissance to a reproducible test of its **dissipative branch** against the actual public DESI DR2 cosmology products.

This task must **not** search for a better Card #1 projection or branch. The projection and the branch are frozen first.

## 2. Hard stop before data acquisition: freeze the spec

Create `conjecture_mode/cards/CARD_01_HORIZON_RELAXOR_SPEC.md` and commit it before the owner threshold ruling. Both must be committed before **any** DESI product is accessed.

### 2.1 Background model
State explicitly:

- flat FRW or otherwise;
- radiation treatment;
- neutrino treatment;
- whether `Omega_m` is fixed or fitted;
- whether `H0` is fixed or fitted;
- every other background parameter used to generate `H(z)`.

No value may be changed after data contact without a new card version.

### 2.2 Exact response law
Write the precise toy law being postulated. No shorthand such as "the record's relaxor form" is acceptable. If the law is `w(a) = -1 + ε X(a)`, define X(a) analytically and define exactly how `tau_2 = H_0^{-1}` enters.

State:

- the domain;
- **the sign convention for ε, such that the dissipative branch is ε < 0 ⇔ w ≤ −1** (or state the record's equivalent condition exactly);
- whether the `H` in the response is the model's self-consistent H or a background H;
- the normalization.

**Required analytic check (pre-data, in the spec):** prove from the written law that ε < 0 implies `w(z) ≤ −1` on the projection domain. If the implication fails anywhere on the domain, the card goes to CARD-01-REPAIR before the owner ruling.

*Slot status:* **EMPTY.** The exact law must be transcribed from the frozen record. No participant has reconstructed it in this document; reconstruction requires separate authorization.

### 2.3 Parameter status
Identify:

- postulated quantities;
- fitted quantities;
- externally fixed quantities.

If ε is fitted on ε ≤ 0, Card #1 is a **one-parameter, one-sided shape family**, not a parameter-free model. The potentially parameter-free object is the **shape relation after projection**, conditional on the frozen projection and the branch.

Information price entry for the branch: **"fixed by origin (second-law interpretation)"**, per C0.4.

### 2.4 Direct model evaluation (primary route, see §4 D3)
Specify how the card's own `H(z)` (or `w(z)`) is passed to the likelihood code:

- which code;
- which likelihood modules;
- the interface (tabulated H(z), or a w(z) table into a Boltzmann/background code);
- how the nuisance and early-universe parameters are handled.

These are frozen here.

### 2.5 CPL projection (diagnostic and visualization only)
Freeze **one** CPL projection definition, e.g. a least-squares fit of `w_model(z)` to `w0 + wa z/(1+z)` on `0 ≤ z ≤ 2` with norm `J = ∫ [w_model − w_CPL]² W(z) dz` and a declared W.

- Pick W before any data contact.
- Do not switch between local-derivative and global-fit definitions.
- The CPL projection **is not the likelihood** and cannot by itself decide a kill under the primary route.

**Required pre-data check:** w ≤ −1 pointwise does **not** automatically give a fitted `w0 < −1`. A weighted least-squares CPL fit can return w0 > −1 if the deviation is concentrated at high z. The spec must therefore compute the frozen-projection image of the ε < 0 branch and report whether it lies in `w0 < −1`. This is a model-only computation and needs no data.

### 2.6 Branch structure
- **Primary (preregistered GRUT claim):** ε ∈ [ε_min, 0], where ε_min is an **OWNER-LOCK slot** (or "unbounded below", with the prior range declared).
- **Control (labeled, not part of the claim):** ε > 0. It may be evaluated and reported as `CONTROL — NOT A GRUT CLAIM`. A favorable control result cannot rescue the card, cannot change the result state, and cannot open a new card under C-F9 unless an independent C0.2 origin is supplied first.

### 2.7 Owner-lock slots (all EMPTY; filled only by the owner ruling)

| slot | meaning | executor may propose? |
|---|---|---|
| S1 `epsilon_null` | half-width of the ΛCDM neighborhood on the dissipative side | yes (pre-data only) |
| S2 primary-route statistic | e.g. profile Δχ² of the best ε<0 point vs ε=0, or a posterior on ε | yes |
| S3 EXPLANATORY-SURVIVES threshold | level at which a non-null ε<0 region counts as viable *and* "meaningfully accounts" for the signal | yes |
| S4 MODEL-KILLED threshold | level at which the nested null / whole manifold conflicts with data | yes |
| S5 INCONCLUSIVE band | when the data cannot discriminate | yes |
| S6 dataset combination(s) | which official DR2 BAO + CMB + SN combination(s) count as primary | yes |
| S7 ε prior range / ε_min | the primary-branch domain | yes |
| S8 CPL diagnostic credible level | used only for the D1/D2 diagnostics | yes |

Every value requires the scientific owner's explicit approval (§4A; the reviewer may audit and recommend). Executor proposals carry the label `CLAUDE-PROPOSAL (non-binding)` and may be computed only from pre-data considerations (§4A).

## 3. Dataset

### 3.0 Network preflight (A6)
Before any data step, log reachability of the official DESI public-release hosts and the exact products named in S6. Use HEAD / listing requests only.

- If any required product is unreachable: **CARD-01-ACCESS-BLOCKED**. Stop. No substitutes (§4B).
- *Known environment fact (disclosed):* the current execution proxy denies arxiv.org, inspirehep.net, doi.org and api.crossref.org. DESI hosts have **not** been tested. ACCESS-BLOCKED is a live possibility.

### 3.1 Products
Use the official public DESI DR2 cosmology likelihood / chain products. Record the following in `CARD_01_DATA_PROVENANCE.md`:

- exact product names;
- release identifiers;
- URLs;
- SHA-256 hashes.

Do not use any of:

- plots read by eye;
- manually assumed covariance ellipses;
- guessed correlations;
- numbers quoted from memory or from papers.

## 4. Data-analysis hierarchy (A5: D3 primary)

1. **D3 — Direct model likelihood (PRIMARY).** Evaluate the card's actual `H(z)/w(z)`, dissipative branch, directly against the public DESI DR2 BAO + CMB + SN likelihoods, using the interface frozen in §2.4. The result state is decided **here**, using S1–S5.
2. **D1/D2 — CPL diagnostics (SECONDARY, non-decisive).** Report the posterior geometry in (w0, wa) and overlay the frozen CPL image of the ε ≤ 0 branch (plus the labeled control). Use direct chain / KDE density evaluation. Use Mahalanobis distance only if Gaussianity is demonstrated. These outputs are labeled `DIAGNOSTIC — NOT THE KILL TEST`.
3. **Downgrade rule.** If D3 cannot be done with the public products, state this explicitly. The card is then graded **CARD-01-C2R (RECONNAISSANCE)**, and its result is worded "NOT FOUND WITHIN THE TESTED APPROXIMATION". The CPL chain overlay may **not** be promoted to an exact or C2-A test. Any CPL-only kill boundary must have been locked separately in the owner ruling as a reconnaissance criterion.

## 5. Preregistered Card #1 kill condition (form only; numbers are owner-locked)

Because ε = 0 is ΛCDM, Card #1 cannot be killed as a whole nested model while ΛCDM remains viable. Only the explanatory claim can die.

> **CARD-01 EXPLANATORY-KILLED:** under the frozen spec and the D3 primary route, either
> - only the owner-locked neighborhood `−ε_null < ε ≤ 0` remains viable at threshold S3, **or**
> - no viable ε < 0 point accounts for the targeted dynamical-dark-energy signal to the S3 standard.
>
> Card #1 is then rejected **as an explanation of the dynamical-DE preference**.

> **CARD-01 EXPLANATORY-SURVIVES:** a viable ε < −ε_null region exists at S3 **and** improves on ε = 0 by the S2/S3 standard, reported with its information price (C-F4).

> **CARD-01 MODEL-KILLED:** even ε = 0 and the whole ε ≤ 0 manifold conflict with the data at S4.

> **INCONCLUSIVE / ACCESS-BLOCKED / RECONNAISSANCE / PROCEDURE-VOID:** per §4C.

*Model-derived structural facts (no data used):*
- On the dissipative branch, w(z) ≤ −1. This is subject to the §2.2 proof from the exact law.
- The CPL image of the branch is a separate, model-only computation (§2.5).
- Because ε = 0 is ΛCDM, MODEL-KILLED requires the nested null itself to fail at S4.

Prior exposure to published DR2 results is recorded only in the card's §4D disclosure. Per ER-3, the spec contains no forecast of its own verdict.

## 6. Future prediction / holdout

Do not call the DR2 comparison a prediction: DR2 already exists.

If Card #1 reaches EXPLANATORY-SURVIVES on DR2 sufficiently far from ε = 0, preregister a future holdout against the DESI complete-five-year dark-energy cosmology release. Use the release name current at preregistration time. The holdout must be registered before those results are public.

A dissipative-branch holdout statement must also be owner-locked before registration.

## 7. Required outputs (on the execution branch, when execution is eventually authorized)

The charter itself stays on the governance branch as `program_governance/CONJECTURE_MODE_CHARTER_01.md`. The execution branch references it and does not copy it as an authority.

- `conjecture_mode/EQUIVALENCE_REGISTRY.md`
- `conjecture_mode/CARD_LEDGER.md`
- `conjecture_mode/CONJECTURE_GENERAL_REGISTER.md` *(new)*
- `conjecture_mode/cards/CARD_01_HORIZON_RELAXOR_SPEC.md`
- `conjecture_mode/cards/CARD_01_OWNER_THRESHOLD_RULING.md` *(new; values explicitly approved by the scientific owner, §4A)*
- `conjecture_mode/cards/CARD_01_NETWORK_PREFLIGHT.log` *(new)*
- `conjecture_mode/cards/CARD_01_DATA_PROVENANCE.md`
- `conjecture_mode/cards/CARD_01_RESULT.md`
- `conjecture_mode/cards/code/` (fresh code and logs)

Cards #2–#4 enter the ledger at their current grades only, with C0.2: NOT AUDITED. Do not rerun them.

## 8. Card #1 grading

| grade | meaning |
|---|---|
| CARD-01-C1-PASS | Spec fully specified, branch check (§2.2) passed, no equivalence or firewall defect |
| CARD-01-C2A-EXPLANATORY-SURVIVES | D3 on real public likelihoods; §5 survives |
| CARD-01-C2A-EXPLANATORY-KILLED | D3 on real public likelihoods; §5 explanatory kill met |
| CARD-01-C2A-MODEL-KILLED | D3; the nested null and the whole manifold conflict at S4 |
| CARD-01-C2A-INCONCLUSIVE | D3; inside the S5 band |
| CARD-01-C2R | D3 infeasible; CPL / reconnaissance only |
| CARD-01-ACCESS-BLOCKED | preflight failed |
| CARD-01-REPAIR | spec defect before a valid C2 test (incl. failure of the ε<0 ⇒ w≤−1 check) |
| CARD-01-PROCEDURE-VOID | owner-lock order violated |

## 9. Scientific firewall

Whatever the result, preserve all of the following:

- Card #1 is POSTULATED, not derived.
- `tau_2 = H_0^{-1}` is POSTULATED.
- The toy response law is POSTULATED unless derived independently elsewhere.
- The dissipative branch is fixed by the record's second-law **interpretation**. That is an origin argument, not a derivation.
- The ε > 0 control is not a GRUT claim in any outcome.
- EXPLANATORY-SURVIVES is not support for GRUT, because ε = 0 is ΛCDM. It must be reported with the null comparison and the information price.
- A good fit does not derive the postulate.
- A bad fit kills only the declared card and branch.
- No result from Card #1 alone establishes primitive vacuum ontology, or GRUT as a fundamental theory.
- A comparison against existing DR2 data is a retrospective test, not a prediction. Only a preregistered future holdout may receive PREDICTED status.

## 10. Execution order (each step needs separate owner authorization)

1. **Done:** charter accepted in substance (A8 kept; ER-1–ER-3 applied). Committed on `grut-program-governance-1` together with `PROGRAM_CAMPAIGN_GATE_01.md`.
2. Owner review and **freeze** of the governance state.
3. Create the descendant execution branch (e.g. `grut-conjecture-mode-1`) from the frozen governance snapshot. Add the equivalence registry, the card ledger and the CONJECTURE-GENERAL register, and enter Cards #1–#4 at the Annex A statuses.
4. Write and commit the Card #1 C0 spec:
   - all fields;
   - the exact law transcribed from the record;
   - the ε<0 ⇒ w≤−1 proof;
   - the frozen-projection image of the ε≤0 branch (model only);
   - the §2.4 interface;
   - **threshold slots empty**;
   - the §4D Prior-Exposure Disclosure from all participants;
   - the executor's non-binding proposals in a separate section, if requested.
5. **STOP → OWNER THRESHOLD RULING.** The reviewer audits/recommends, and the scientific owner explicitly approves S1–S8. The ruling is committed as `CARD_01_OWNER_THRESHOLD_RULING.md`.
6. Network preflight. On failure: ACCESS-BLOCKED, stop.
7. Acquire the official DESI DR2 products; record provenance and hashes.
8. Run D3 once (plus the D1/D2 diagnostics and the labeled control).
9. Write the result in §4C vocabulary.
10. Do not change the projection, branch or thresholds in the same card.
11. Stop for owner review.

The following are **not** part of this task:

- any Card #2 repair;
- a new cosmology card;
- a TT / conformal-infinity campaign;
- a selector campaign;
- any canonical GRUT modification;
- any PR or merge.

---

## 11. Terminal principle

Conjecture Mode exists to permit

> **specific wrong ideas that can die cleanly, originating in GRUT's own structure,**

without allowing

> **hidden assumptions to masquerade as derivations, or data shapes to masquerade as GRUT origins.**

The standard is not that a postulate must be earned before it is tested. The standard is that its **origin, information price, branch, projection and owner-locked kill thresholds** are visible before the result is known.
