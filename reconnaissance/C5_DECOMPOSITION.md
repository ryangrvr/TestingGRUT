# C5 DECOMPOSITION

C5 is not one object. Each component below is stated as **what would count as deriving it**, together with its known
alternatives and the hostile baseline.

| ID | Component | What would count as deriving it | Known inequivalent alternatives | Hostile baseline (known failure / non-uniqueness) |
|---|---|---|---|---|
| **C5-A** | individuation / system identity | from primitive dynamics or relations, a **unique** partition into subsystems (plus environment) that is fixed by a stated, non-arbitrary criterion | any unitary / bijective relabeling of a finite state space gives a different factorization (Zanardi: TPS are observable-induced) | for a bare Hilbert space, *every* factorization `d = d₁d₂` is unitarily equivalent: no individuation without extra structure (an algebra of accessible observables, or H) |
| **C5-B** | composition rule | given individuated parts, the composite state / effect space is **forced** | Cartesian (classical), tensor (QM), graded tensor (fermions), min vs max GPT tensor products (boxworld), Jordan composites, direct sum (superselection) | GPT composites range between the min and max tensor products; the choice is free unless further principles are added |
| **C5-C** | locality / independence | operational independence of A and B (separate specifiability) **forced**; the six notions kept distinct: spatial, operational, statistical, causal, commutativity, factorization | commuting-subalgebra (Tsirelson / field-theory) vs tensor factorization; the split property | the Tsirelson problem: commuting-operator ≠ tensor-product correlations in infinite dimensions (MIP* = RE) |
| **C5-D** | state / effect structure | operational states and effects arise as equivalence classes of preparations / tests of an underlying dynamics, with **no** declared probabilistic framework | ontological models (Spekkens), Koopman–von Neumann, epistemic restrictions (Spekkens toy) | "state = equivalence class" needs an equivalence criterion, i.e. access (circularity risk) |
| **C5-E** | probability (affine rule) | `p(mixture) = mixture of p` forced for accessible statistics, with the measure **earned** | frequentist (needs a typicality measure), Bayesian (agent priors), invariant measure (needs initial data) | the measure problem: typicality presupposes a measure |
| **C5-F** | local tomography | global states determined by local statistics, **from** individuation + composition + access | real QM, fermionic QT, superselection, min-tensor GPTs | real QM ≅ complex QM + a global rebit reference (a hidden global phase): LT fails because of a missing shared reference frame |
| **C5-G** | dimensionality | graph / spectral / spacetime / state-space / Hilbert / capacity dimension selected | any dimension consistent with the axioms; causal-set dimension estimators; spectral dimension running | most dimension notions are inputs (the SCOUT-1 W1-S continuum) |
| **C5-H** | basin / preparation / measure | initial-condition and basin data selected without relocating the assumption | maximum entropy, typicality, SRB / physical measures, attractors, cosmological initial states (Past Hypothesis) | "almost all" is measure-relative; SRB measures need Lebesgue-a.c. initial data (a measure price) |

**Dependency sketch (target chain):**

    A individuation → B composition → C independence → D states/effects → E affine probability → F tomography
                                                                        ↑ H measure / basin enters D and E
    G dimension enters A (capacity) and B (K_AB rule)

**Circularity hazards** (flagged for every probe):
- A defined by C (access);
- D defined by access;
- E defined by a supplied measure (H);
- F derived from purification, which already contains B.
