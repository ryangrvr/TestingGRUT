# WAVE-1 PROBE CHARTERS (pre-registered before any probe runs)

Every probe reports:
- *structure inserted* vs *structure forced*;
- the measure used (MEASURE_ORIGIN_LEDGER);
- the access / intervention set used;
- the information-accounting row.

## S2-1 — System individuation (RUN FIRST)

- **Question.** Starting from (Hilbert space, H) with **no declared tensor-product structure (TPS)**, is a subsystem
  decomposition forced?
- **Tests:**
  - (a) **No-go baseline (#46):** without a criterion, every TPS of a given total dimension is unitarily equivalent.
    So H alone (its spectrum) fixes nothing beyond `dim` and the spectrum.
  - (b) **Locality selector (#1, Cotler–Penington–Ranard-type):** the criterion is "H is 2-local w.r.t. a qubit TPS".
    Compute the Jacobian rank of the map (2-local couplings) → (spectrum), modulo local unitaries (dim 3n), for
    n = 3…9 qubits.
    - The 2-local TPS is **locally unique** iff rank = P − 3n, where P is the number of couplings.
    - It is **non-unique** iff the fibres have positive dimension.
  - (c) **Explicit hostile (#2):** for small n, construct two **inequivalent** 2-local Hamiltonians with the same
    spectrum. Inequivalent means different local-unitary invariants. Their TPSs differ.
- **Prices to record:**
  - the locality criterion (k);
  - the local factor dimension (qubits);
  - the n-threshold;
  - that "most spectra admit no local TPS" (a measure statement: which measure?).
- **Outcomes:** SYSTEM INDIVIDUATION SELECTOR (scoped) / SUBSYSTEM NON-UNIQUENESS / both by regime.

## S2-2 — Composition rule

- **Question.** Given two individuated systems (state spaces as convex sets), which composite is forced?
- **Tests:**
  - compare the min tensor, max tensor, quantum tensor, Cartesian, graded and direct-sum composites by:
    (i) independent preparability; (ii) existence of reversible interacting dynamics; (iii) no-cloning / copying;
  - squares (gbits) and discs.
- **Firewall:** local tomography is **not** an input.
- **Outcomes:** COMPOSITION FORCED BY [primitive] / NONUNIQUE.

## S2-3 — Probability from deterministic microdynamics

- **Question.** Does deterministic (chaotic) dynamics + coarse-grained access force accessible statistics to be an
  affine probability rule with an **earned** measure?
- **Firewall:** an invariant measure existing ≠ a probability rule. A reason experiments *must* use it is required.
- **Tests:**
  - mixing maps (Arnold cat, baker, logistic r = 4) vs non-mixing (rotation);
  - frequencies of coarse cells from many initial conditions **drawn how?**;
  - sensitivity to the initial-condition measure (Lebesgue vs singular).
- **Outcomes:** MEASURE-PRICED / TRUE DERIVATION (scoped) / NO-GO.

## S2-4 — Convexity / mixtures

- **Question.** Is `state = preparation equivalence class` + randomized preparation ⇒ convexity? Can endogenous
  deterministic randomness (a chaotic coin) replace supplied randomness?
- **Hostile:** a coin built from a deterministic map is affine only relative to a measure on its seed. Locate it.
- **Outcomes:** CONVEXITY DERIVED (scoped) / MEASURE-PRICED / DEFINITIONAL.

## S2-5 — Local-tomography emergence

- **Question.** What primitive property removes the hidden global degrees of freedom that break LT in real QM,
  fermionic QT and superselected theories?
- **Tests:**
  - real QM + a global rebit reference reproduces complex QM statistics (ABW-type encoding): the non-local parameter
    is a **missing shared reference frame**;
  - superselection lifted by a reference system;
  - a fermionic ancilla.
- **Target:** a primitive reason (e.g. "no unobservable global frame / universal access / factor completeness"), not
  another operational axiom.
- **Outcomes:** LT FROM [primitive] / ACCESS-PRICED / CONTENT REPACKAGING.

## S2-6 — Basin / measure selection

- **Question.** Can generic dynamics select preparation / state information (a physical measure) without a supplied
  measure?
- **Tests:**
  - SRB / physical measures for chaotic maps: uniqueness of the measure attracting Lebesgue-a.e. initial data;
  - hostile: ergodic-decomposition non-uniqueness (several invariant measures);
  - the selector "attracts almost every initial condition" is **Lebesgue-priced**.
- **Outcomes:** BASIN SELECTED (MEASURE-PRICED) / NONUNIQUE.

## S2-7 — Consistent-histories set selection (observer / access)

- **Question.** Do consistency conditions select a unique quasi-classical history set (thus systems / records)?
- **Hostile:** Dowker–Kent-type many inconsistent alternative consistent sets.
- **Outcomes:** SELECTOR / NONUNIQUE.

## S2-8 — Darwinism / records

- **Question.** Does redundancy of environmental records select pointer observables **and** a system/environment
  split, or does it presuppose the split?
- **Outcomes:** SELECTOR (given split) / SPLIT-PRICED.
