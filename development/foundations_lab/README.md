# Foundations lab — issue #3

Start with [the decision report](REPORT.md). Recommendation:
**CONTINUE GRUT BUT BROADEN/REFRAME**, preserving the operational audit while
requiring an extra measurable physical restriction before fundamental promotion.

**DRAFT — external scientific review pending.** This is an isolated
foundations/control lane, based on engineering head
`b5d0493a23e5b741397f4dd3867bd8a09da08858`. It neither repairs nor clears the
engineering blockers and does not implement CR-5. No candidate K, mechanism
score, L0 price, IP-12 draft, card, canonical ruling or bank entry is created.
Card 1 remains closed; budget consumed **zero**. Scripts replace notebooks
to keep reproduction dependency-free.

| File | Purpose |
|---|---|
| [REPORT.md](REPORT.md) | Required eight-part synthesis, three different routes, budget-type recommendations and decision |
| [TARGETS.md](TARGETS.md) | Ten distinct missing-physics questions and eight-target ranking |
| [QIH_AUDIT.md](QIH_AUDIT.md) | Primary PDF audit, central equations, revisions, limits and hostile conclusions |
| [FRAMEWORKS.md](FRAMEWORKS.md) | 25 framework families, selection mechanisms and strongest GRUT arguments |
| [EXPERIMENTS.md](EXPERIMENTS.md) | Dated primary-source opportunity sweep, evidence categories and limitations |
| [OBSERVERS.md](OBSERVERS.md) | Detector, records, agents, functional and phenomenal consciousness |
| [PROOFS.md](PROOFS.md) | Exact standard proofs, counterexamples, source references and scope |
| [RESULTS.json](RESULTS.json) | Machine-readable rational/numerical witnesses and code identities |
| [SOURCES.json](SOURCES.json) | 86-source ledger with URLs, dates, inspection levels and source hashes |
| [VALIDATION.json](VALIDATION.json) | Local reproduction and protected-input comparison evidence |
| [controls.py](controls.py), [test_controls.py](test_controls.py), [run.py](run.py) | Dependency-free reproducible control package |

## Reproduce

From repository root with **Python 3.12.14**:

```sh
python development/foundations_lab/run.py --check
```

This executes 13 meaningful standard/control checks and verifies that the
checked-in JSON matches current code and the pinned Python version. It fails
on missing/stale/different output. To regenerate after a disclosed code change:

```sh
python development/foundations_lab/run.py
```

Fractions serialize as `{"rational": [numerator, denominator]}`. Transcendental
evaluations are labelled FLOAT; their exact formulas/proofs are separate.
Finite samples do not certify an infinite theorem. Same-author calculation
and literature review are not independent reproduction or scientific approval.

## Controls and how they matter

| ID | Witness | What it prevents |
|---|---|---|
| C01 | Same Z probabilities, different X probabilities; QFI/FS convention | Calling one qubit marginal a generated physical readout theory |
| C02 | Reciprocal cosine >1 and SI unit mismatches | Accepting physical identifications without units or valid geometry |
| C03 | Positive plus-sign finite determinant, incompatible with real xi zeros | Calling an internally inconsistent spectral target a selector |
| C04 | Chiral trace=0 and divergent Gibbs trace | Inferring thermal existence/vacuum stress from self-adjointness alone |
| C05 | 3-to-2 dimensional encoding cannot be an isometry | Hiding a code-subspace/capacity premise |
| C06 | SLD versus BKM: 4 versus 4 ln 3 | Claiming one quantum metric follows from monotonicity alone |
| C07 | Marginal purity changes under a different tensor factorization | Demanding a unique subsystem from a bare Hilbert space |
| C08 | PR box: normalized, positive, no-signalling, CHSH=4 | Claiming no-signalling alone reconstructs quantum theory |
| C09 | Real-QM states agree locally but differ globally | Losing the composition assumption in reconstruction claims |
| C10 | Reversible models: same statics/unforced generator, different finite-drive transition | Treating missing microscopic input as demonstrated missing fundamental physics |
| C11 | Conditional anomaly constraints fix ratios but leave a normalization | Confusing genuine conditional selection with generation of all premises |
| C12 | Prime arithmetic and X-dephasing counterexample | Inferring dark mass or universal protection from an integer label |
| C13 | Unit-circle speed=1 versus claimed omega-dependent speed | Accepting inconsistent simultaneous definitions even at c=1 |

## Scope and review

The full higher-jet/three-state kinetic manuscript and fixture ZIP were not
retrieved here; their claims are identified as supplied-record claims. Five
QIH PDFs were downloaded and hashed; their text was used for targeted audit,
not copied into this repository. PDF equation pages for the spectral target
were visually inspected. Source abstracts/reviews support the wider
comparison; paywalled full-text access is not implied.

Issue #2's expected-red manifest remains byte-identical. The new workflow
only tests standard/control reproduction. Existing engineering checks retain
their separate BLOCKED / REVIEW_REQUIRED statuses; a green foundations job
does not weaken, hide or adjudicate those findings.

All files live under the development namespace (plus this workflow and the
coordination log). The canonical GRUT repository was read but never edited.
The draft PR is a noncanonical review surface; it is not a scientific merge.
