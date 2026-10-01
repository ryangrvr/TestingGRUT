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

---

# REPAIR 01 PRE-REGISTRATIONS (owner audit; written before any run)

Frame: T2-1′ — C5 → D ⊕ H ⊕ A. D = primitive-dynamics structure; H = state / basin / measure / preparation;
A = access / intervention / readout / reference-sharing. Each probe below attacks one wall. Every verdict states which
of D / H / A it moved, and which it priced.

## S2-3b — Zero-entropy uniquely-ergodic mixing (horocycle-type); attacks H

**Firewall: three notions, pre-registered as distinct. None may be identified with another.**

| Notion | Definition | Not the same as |
|---|---|---|
| **A. UNIQUE MEASURE** | the dynamics admits exactly one invariant Borel probability measure | B (B can hold off a null/exceptional set while A fails) |
| **B. TIME-AVERAGE UNIVERSALITY** | for every (or every explicitly characterized) initial point, Birkhoff averages of continuous observables converge to the same value | C (B is a single-trajectory statement; it says nothing about ensembles) |
| **C. PREPARATION FORGETTING** | pushforwards of distinct preparations become indistinguishable to the *observables an agent can access* | B, and mixing as such |

Hard rule: for an **invertible measure-preserving** flow, fine-grained densities remain exactly distinguishable
(TV / L¹ distance conserved). Any "forgetting" is therefore a statement about an **observable class (A)** and a
**preparation class (H)**, never about fine-grained information.

**System.** Horocycle flow h_t = [[1,t],[0,1]] on X₂ = SL(2,ℝ)/SL(2,ℤ) (unimodular lattices in ℝ²), simulated by
Gauss reduction. X₂ is **non-compact**: this is the hostile variant. The compact-quotient case (Furstenberg: uniquely
ergodic; Marcus/Ratner: mixing with rates) is KNOWN RESULT IMPORT, not simulated.

**Measurements (separately reported):**
1. **A:** exhibit or exclude multiple invariant measures (periodic horocycles: lattices with a horizontal vector;
   Dani's classification — SECONDARY).
2. **B:** Birkhoff averages of `1[|v₁|² < s]` vs the Haar value `3s/π` (s ≤ 1, analytic) from several starting lattices:
   Haar-random, irrationally rotated ℤ², nearly horizontal, exactly periodic.
3. **C-i decay of correlations:** Haar-sampled ⟨f · g∘h_t⟩ − ⟨f⟩⟨g⟩ for t up to O(10²).
4. **C-ii weak convergence:** ⟨f⟩ under pushforwards of (a) two disjoint a.c. blob preparations, (b) a Dirac
   preparation, (c) a periodic-orbit preparation.
5. **C-iii fine-grained conservation:** the pulled-back indicator `1_{B₁} ∘ h_{−t}` distinguishes the pushed blobs with
   contrast 1 at every t; forward/back round-trip error; growth with t of the number of coarse cells needed to
   describe h_t(B₁) (zero entropy → expected polynomial, contrasted with the doubling map's exponential growth).

**Outcomes (pre-registered):**
- **H-MEASURE SELECTED / A-COARSE FORGETTING:** A (compact case) earned from D; C holds only for accessible
  observables and a.c. preparations.
- **H UNTOUCHED:** forgetting requires a separately supplied reference measure.
- **NO-GO (scoped):** uniqueness and forgetting cannot co-occur without special (homogeneous, rigid) D.
- **D-PRICED:** the selection rests on compactness / algebraic homogeneity that is itself inserted.

## S2-1b — Can locality be selected without supplying graph distance, local dimension or k? Attacks D

**Question.** Given only an abstract dynamics (a Hamiltonian as a spectrum, or a matrix in an arbitrary basis on
`ℂ^N`), do competing locality criteria select the **same** tensor-product structure?

**Criteria, all run on the same abstract dynamics:**
1. minimal k (smallest k such that the Hamiltonian is k-local in some TPS);
2. sparsest interaction graph;
3. Lieb–Robinson velocity / light-cone sharpness;
4. MDL (shortest description of H as a sum of local terms);
5. stability (TPS robust to perturbation of H);
6. locality-maximizing factorization (Zanardi / CPR-type);
7. predictive autonomy (subsystems whose reduced dynamics is most nearly autonomous).

**Hostile.** Construct dynamics where two criteria choose **different** TPSs (including different factorizations
`N = d₁·d₂…`). Then ask: **what selects the objective function?** If nothing in the dynamics does, the criterion is
an inserted D-item (or an A-item, if it is defined by what an agent can control).

**Outcomes:** LOCALITY SELECTED (criteria agree, no inputs) / CRITERION-PRICED (D) / ACCESS-PRICED (A) / NONUNIQUE.

## S2-7 — Consistent histories: set selection. Attacks A

- **Question.** Does consistency (decoherence functional) select a unique quasi-classical set of histories, and
  hence records and subsystems?
- **Hostile:** Dowker–Kent: many mutually incompatible consistent sets; consistent sets that are not quasi-classical.
- **Firewall:** if the selection requires a fixed coarse-graining, a fixed system/environment split or a fixed time
  sequence of projectors supplied from outside → **A-PRICED**.
- **Outcomes:** SELECTOR / NONUNIQUE / A-PRICED.

## S2-8 — Quantum Darwinism / records. Attacks A

- **Question.** Does the redundancy of environmental records select pointer observables **and** the
  system/environment split, or only the pointer observable given the split?
- **Firewall:** if redundancy is defined only after the split (S | E₁ ⊗ … ⊗ E_m) is specified → **A-PRICED**
  (and possibly D-priced via S2-1b). A test that the split itself can be chosen by redundancy maximization over all
  TPSs is required before any claim of selection.
- **Outcomes:** SELECTOR (split and pointer) / SELECTOR (pointer, given split) / A-PRICED.

## S2-G — Dimension selection (C5-G). Kept active

Four notions, tested separately and never identified:
1. **local Hilbert dimension** d (the factor size of the TPS; S2-1 / S2-1b);
2. **graph / spectral dimension** of the interaction graph;
3. **spacetime dimension** (from Lieb–Robinson cones or correlation decay);
4. **information capacity** (log dimension per unit region / per record).

**Outcomes per notion:** SELECTED / PRICED (D, H or A) / NONUNIQUE.
