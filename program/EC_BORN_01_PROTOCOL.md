# E&C-BORN-01 — BORN-RULE PROVENANCE AUDIT PROTOCOL (FROZEN)

> Freeze-only artifact. Created before any case work, source adjudication,
> scoring, or classification. The audit object is the *provenance* of Born-rule
> derivations in the existing literature — not GRUT, not a new interpretive
> framework, and not a GRUT-specific Born rule.

---

## 0. ROLE AND BOUNDARY

This protocol governs a future execution of E&C-BORN-01. The freeze itself
performs **no literature adjudication, no scoring, no classification, and no
case evidence collection**. It fixes: the question, the route families, the
required fields, the classification vocabularies, the firewalls, and the
admissibility rules.

The research target is REALITY. GRUT is the umbrella research program, not a
privileged theory. **No GRUT-designed Born rule is authorized by this
protocol** — during freeze or execution.

Never: bare `git push`. Push target (only if separately authorized):
`testinggrut master`. Canonical GRUT-RAI: UNMODIFIED.

---

## 1. PRIMARY QUESTION

> For each serious Born-rule derivation route in the literature, what exactly
> does the derivation assume, and what exactly does it produce?

Not: "Can we explain quantum mechanics?" Not: "Which interpretation do we
prefer?"

## 2. GOVERNING RULE (verbatim)

$$
\boxed{\text{A derivation of } p_i=|\alpha_i|^2 \text{ is insufficient unless
the assumptions needed to obtain it are themselves classified.}}
$$

A mathematical argument reaching $|\alpha_i|^2$ does not count as derivation
of the Born rule if the argument assumes the very probability structure it
claims to explain. Equally, an audited argument whose premises are
independently classified as weaker than Born content **earns a positive
result**.

## 3. ADVERSARIAL HYPOTHESIS (tested, not assumed)

    H_BORN:
        Some existing Born-rule derivations import probability structure
        equivalent to what they claim to derive.

H_BORN is a hypothesis to be tested per route. It may fail for some routes.
Do NOT preregister "all derivations are circular" as the expected result. The
two-directional attack requirement (§6) exists precisely to prevent both
reflexive circularity-calling and reflexive credulity.

---

## 4. FROZEN TARGET SEPARATION

Every route must be assessed independently for each of:

| Target | Meaning |
|---|---|
| `WEIGHT_REPRESENTATION` | A theorem uniquely fixes a functional form for a probability/measure assignment |
| `OBJECTIVE_CHANCE` | The framework explains physical single-run chances or frequencies |
| `RATIONAL_CREDENCE` | The framework constrains what a rational agent should expect |
| `BRANCH_MEASURE` | The framework defines a measure over Everett branches |
| `OUTCOME_SELECTION` | The framework explains why one outcome occurs (if it claims one does) |

**No target automatically transfers credit to another.** Frozen non-implications:

    WEIGHT_REPRESENTATION ⇏ OBJECTIVE_CHANCE
    RATIONAL_CREDENCE     ⇏ OBJECTIVE_CHANCE / physical probability
    BRANCH_MEASURE        ⇏ OUTCOME_SELECTION / one realized outcome
    OUTCOME_SELECTION     ⇏ derived Born weights

A Deutsch–Wallace result may at most target RATIONAL_CREDENCE and/or
BRANCH_MEASURE; it must not receive credit for single-outcome selection. A
collapse model may select outcomes but receives no automatic credit for
deriving its stochastic measure.

---

## 5. FROZEN ROUTE FAMILIES (broad-family level; no substitutions after
case work begins)

    R1. Gleason / Busch / noncontextual-measure routes
    R2. Envariance / symmetry routes
    R3. Deutsch–Wallace Everett decision-theoretic routes
    R4. Typicality / frequency / self-locating uncertainty routes
    R5. Quantum reconstruction / operational axiomatization routes
    R6. Bohmian quantum-equilibrium / typicality routes
    R7. Bohmian nonequilibrium and counterexample class
    R8. Collapse stochastic-law routes (DP admitted ONLY as audited
        contextual prior work; DP-R0-01 conclusions are context, not
        evidence for R8's classification)

Newly discovered serious routes during execution are recorded as
`OUT_OF_SCOPE_CANDIDATES` for a future protocol version — never added
mid-run. No GRUT Born-rule route exists in this registry.

---

## 6. TWO-DIRECTIONAL ATTACK REQUIREMENT

For every route, execution must run BOTH attacks:

    Attack 1: attempt to show the derivation imports Born-equivalent
              probability structure.
    Attack 2: attempt to show the imports are genuinely weaker than the
              Born rule for the route's claimed target scope.

A route surviving both attacks with independent premises earns
`DERIVED_WITH_INDEPENDENT_PREMISES` (or `CONDITIONAL_NONCIRCULAR`) and the
protocol must award it without embarrassment. Skepticism bias is a protocol
failure, as is credulity bias.

---

## 7. MANDATORY PER-ROUTE SCHEMA FIELDS

    route_id
    target_claims
    formal_result_type
    mathematical_domain
    exact_assumptions
    assumption_type
        MATHEMATICAL | PROBABILITY | ONTOLOGICAL | DECISION_THEORETIC |
        PHYSICAL_DYNAMICAL | TYPICALITY | SEMANTIC
    assumption_status
    assumption_independence_status
    derived_conclusion
    what_is_not_derived            (MANDATORY — most information-dense field)
    circularity_status
    relation_to_outcome_selection
    relation_to_Cbp_boundary
    empirical_leverage
    criticisms_and_replies         (primary-source objections + replies;
                                    preserve disputes, do not paraphrase
                                    into strawmen)
    criticism_admission_rule
    reply_admission_rule
    final_disposition

Every assumption must be CLASSIFIED, not merely listed. Illustrative
classifications (freeze-time illustrations, not judgments):

- "noncontextual measure" (R1): substantive PROBABILITY assumption —
  whether it is circular depends on the claimed target;
- "rational preferences over quantum games" (R3): DECISION_THEORETIC
  assumption, not a physical collapse mechanism;
- "quantum equilibrium" (R6): TYPICALITY/measure assumption, distinct from
  the guidance law;
- "typicality" (R4): may yield a conditional measure claim without
  explaining unique realized outcomes;
- a collapse model's stochastic measure (R8): occurs in an evolution
  equation, but occurrence is not derivation.

---

## 8. FROZEN CIRCULARITY_STATUS VOCABULARY

    NOT_CIRCULAR
    CONDITIONAL_NONCIRCULAR
    PROBABILITY_STRUCTURE_IMPORTED
    TARGET_REDEFINED
    CIRCULAR_OR_TARGET_ASSUMED
    DISPUTED
    UNRESOLVED

"Uses probability somewhere" does NOT automatically imply circular.
"Contains a probability assumption" does NOT automatically disappear into a
vague mathematical premise. Each route's circularity status must be reached
through the premise classification of §7.

## 9. FROZEN FINAL-DISPOSITION VOCABULARY

    DERIVED_WITH_INDEPENDENT_PREMISES
    CONDITIONAL_DERIVATION
    WEIGHT_REPRESENTATION_ONLY
    RATIONAL_CREDENCE_ONLY
    BRANCH_MEASURE_ONLY
    OUTCOME_SELECTION_ONLY
    POSTULATED_MEASURE
    REFRAMED_TARGET
    ASSUMPTION_DEPENDENT_AND_CONTESTED
    INSUFFICIENTLY_SPECIFIED

Different routes may successfully establish different things. Do NOT force
one universal final verdict; per-route dispositions are the primary output.
A cross-route synthesis (form/weights/frequencies/outcomes coverage matrix)
is produced only AFTER per-route classification (execution Phase 5).

---

## 10. Cbp BOUNDARY FIREWALL (verbatim)

    Cbp boundary:
      scope information only.
      It is not disconfirmation evidence against a route unless the route
      demonstrably lies within the theorem's premises and claims the same
      target.

The Experiment-P / hostile-replication boundary result (commits 10681b0,
8ebac8b, 793143e) shows that internal dynamics of the branch-preserving
unitary class preserve inherited branch amplitudes. It does NOT establish a
universal no-go. Gleason, envariance, decision-theoretic, typicality,
reconstruction, Bohmian, and collapse routes are NOT rejected by it in
advance. It may define a comparison class, nothing more.

---

## 11. HARD FIREWALLS (all frozen; all independently testable)

    decoherence ≠ outcome selection
    outcome selection ≠ probability law
    Born-weight reproduction ≠ Born-weight explanation
    diagonal reduced density matrix ≠ outcome
    branch existence ≠ outcome selection
    Envariance's factorization choice is a premise, not a result
    rationality axioms are premises, not physics
    typicality measures are premises unless independently justified
    Bohm |psi|^2 equilibrium is a hypothesis; the quantum-nonequilibrium
      literature (R7) is preserved as the standing counterexample class
    collapse models supply weights BY POSTULATE unless a route derives
      the stochastic measure
    reproduction of observed statistics is not explanation of the weights
    a derivation inside framework F does not establish that reality
      requires F (WITHIN_FRAMEWORK discipline, inherited from E&C-01 RA1)
    derived ≠ assumed ≠ effective ≠ empirical_input ≠ contested
    mathematical uniqueness of a measure ≠ physical origin of objective
      probability

---

## 12. FROZEN EXECUTION ORDER (never reversed)

    PHASE 1  secure primary sources per route
    PHASE 2  premise ledger construction (full §7 fields)
    PHASE 3  two-directional attack
    PHASE 4  per-route classification (§8, §9 vocabularies)
    PHASE 5  cross-route matrix: form / weights / frequencies / outcomes
             coverage
    PHASE 6  hostile audit of the execution
    PHASE 7  integration with escape-door map + compression matrix
             (read-only references; no score modification)

Inherited committed artifacts (Experiment P, boundary map, escape-door map,
DP-R0-01) are CONTEXT ONLY. None of their conclusions may serve as
evidence for a route's classification.

## 13. SOURCE-ADMISSION AND CRITIQUE-ADMISSION RULES

Admissible: primary research papers; authoritative reviews; rigorous
mathematical literature; canonical books by the derivation's principal
authors. For every disputed route, at least one source representing the
principal objection must be admitted alongside the derivation's own source.

Not admissible as sole evidence: AI-generated overviews, blog summaries,
secondary paraphrases. No substantive source adjudication occurs during
freeze — admission rules only.

---

## 14. OUT OF SCOPE (explicit)

    - constructing a GRUT Born rule;
    - proposing new interpretive frameworks;
    - scoring/ranking interpretations for preference;
    - altering the escape-door map or compression matrix;
    - expanding E&C-01 chains;
    - cosmology, gravity, constants, observer theory;
    - repairing the quarantined cosmetic-rewrite set or the five malformed
      comparative-foundations JSONs (separate process).

## 15. FUTURE EXECUTION ARTIFACT NAMING

    program/EC_BORN_01_SOURCES.json/.md
    program/EC_BORN_01_ROUTE_R1.json ... EC_BORN_01_ROUTE_R8.json
    (one per frozen family, .md companions)
    program/EC_BORN_01_MATRIX.json/.md
    program/EC_BORN_01_REPORT.json/.md
    program/EC_BORN_01_HOSTILE_AUDIT.json/.md

Execution runs only after separate authorization following inspection of
this frozen protocol.

---

*Canonical GRUT-RAI: UNMODIFIED. This and all companion artifacts:
TestingGRUT laboratory record. Freeze-only run; nothing staged, committed,
or pushed.*
