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

---

# WAVE 2 PRE-REGISTRATIONS (owner review of `59d77f6`; written before any run)

## S2-Σ (+ S2-G1) — factorization AND local dimension selection without hidden weights

**Hostile question.** Can (H, ψ) select a factorization *and* a local dimension in an objective-independent way?
Or does every proposed selector import a utility function, a scale or a factor dimension?

### Σ-0 Firewall

- A unique maximum of a hand-chosen weighted score F = w_L L + w_Q Q + … does **NOT** count as a selected TPS. Such
  an F is **SELECTOR-PRICED** unless its weights and normalizations are derived.
- Criteria are computed independently. Each is oriented so that larger is better:

| criterion | definition |
|---|---|
| **L** locality | −k_eff, the ‖H‖²-weighted mean factor-support size of H |
| **Q** quasiclassical stability (Carroll–Singh-type) | −(mean over factors of the entanglement-entropy growth S_f(τ) − S_f(0)) under H for the state ψ |
| **R** record redundancy | max over factors f of #{g ≠ f : I(f:g)(τ) ≥ (1−δ) S_f(τ)}, counted only when S_f(τ) > 0.1 bit |
| **P** predictive autonomy | the ‖H‖² fraction on single-factor supports |
| **M** description length | −(dimension of the support-closed operator space containing H) / N² |

- All criteria are LU-invariant (exact-factor-support decomposition; factor entropies).
- **Admissible candidates** have ≥ 2 factors. This is a supplied constraint: the trivial factorization wins L, P and M
  degenerately (S2-1b).
- **Candidate set** (a supplied finite family; recorded in Σ-9): every set partition of n = 6 qubit slots into ≥ 2
  blocks (202 groupings, so local dimensions 2…32 and mixed), crossed with frames:
  - the identity;
  - random Clifford circuits;
  - one Haar-random global unitary;
  - the commutant frame W = e^{−iHs}.

**Outcomes (more than one may hold):**
- **A. ROBUST TPS SELECTOR:** one candidate weakly dominates all others.
- **B. PARETO NONUNIQUENESS.**
- **C. DIMENSION-PRICED.**
- **D. SCALE-PRICED.**
- **DEGENERACY-PROTECTED NONUNIQUENESS.**

### Σ-1 Pareto front first

- Compute V(Σ) = (L, Q, R, P, M) for every candidate, with no scalarization.
- Report:
  - the front size;
  - whether a single candidate weakly dominates all the others;
  - the block-size types (local dimensions) and frames present on the front.

### Σ-2 Monotone-transformation hostile

- Pareto dominance is invariant under every strictly increasing componentwise transform. That invariance is stated
  and checked.
- Scalar winners are recomputed under transforms (log1p of shifted values, square, sqrt) and normalizations (min-max,
  z-score, rank).
- Report the number of distinct winners.

### Σ-3 Weight simplex (only after Σ-1)

- Sample Dirichlet(1) weights on the 4-simplex. This sampling measure is supplied and named.
- Report the win share per candidate under each normalization, and whether any candidate wins > 90%.
- Calling a weight vector "natural" is forbidden.

### Σ-4 Local dimension as a variable

- No fixed "qubits". Block-size multisets range over {1⁶}, {2,2,2}, {3,3}, mixed, and so on.
- **TPS choice** and **local-dimension choice** are recorded separately.
- If the objective family cannot jointly select d → **LOCAL DIMENSION REMAINS A Σ PRIMITIVE**.

### Σ-5 Scale hostile

Vary:
- the time horizon τ;
- the redundancy threshold δ (record resolution);
- the fragment definition (single factors vs pairs);
- the perturbation strength ε (GUE added to H; measure named).

A selector is strong only if its front / winner is stable over a nontrivial window. Otherwise → **SCALE-PRICED**.

### Σ-6 State dependence

- The same H is evaluated in several states:
  - a product state;
  - a random-local product state;
  - the ground state;
  - a mid-spectrum eigenstate;
  - a Haar-random state;
  - a record-forming state;
  - a Gibbs state (β = 1).
- If ψ changes the winner, that cost moves explicitly into **H** (the state chooses the subsystems).

### Σ-7 Symmetry / degeneracy hostile

- **(a)** A translation-symmetric ring: groupings related by translation.
- **(b)** The **commutant frame** W = e^{−iHs}. H-only criteria (L, P, M) are exactly invariant under it, while
  ψ-criteria transform as a time shift.
- Identical vectors for inequivalent TPSs → **DEGENERACY-PROTECTED NONUNIQUENESS**.

### Σ-8 Literature conflict, decided by toy models rather than authority

| work | claim | assumptions to identify from the toy models |
|---|---|---|
| Cotler–Penington–Ranard (arXiv:1702.06142) | the spectrum determines a local TPS when one exists | which assumptions this needs |
| Carroll–Singh (PRA 103, 022213) | quasiclassicality defines a preferred factorization | which assumptions this needs |
| Stoica (arXiv:2103.15104) | structures from (H, ψ) alone cannot be both physically relevant and unique | which assumptions this needs |

Each literature statement is classified by scope.

### Σ-9 Information accounting

Separate every entry:
- H;
- ψ;
- total dimension N;
- local dimension;
- number of factors;
- objective family;
- weights;
- normalization;
- time scale;
- coarse scale / δ;
- fragment definition;
- the measure over candidate TPSs (the candidate family).

**Σ DERIVED** only if every material entry except H and ψ is eliminated or shown to be gauge.

## S2-G — dimension, split into four questions (kept separate in the ledger)

| Question | Content | Run with |
|---|---|---|
| **G1** local Hilbert factor dimension | jointly with Σ (above) | S2-Σ |
| **G2** graph / spectral dimension | from N(r) ~ r^d, spectral density, return probability. The graph, metric or diffusion operator used is priced | after the Σ zoom-out |
| **G3** spacetime dimension | not inferred from graph dimension. Does causal / dynamical propagation fix an effective spacetime dimension uniquely? | after the Σ zoom-out |
| **G4** operational / information capacity | the maximum distinguishable-state count: derived from D + H + Σ, or supplied? | after the Σ zoom-out |

## S2-H2 — preview (do NOT run before the Σ/G1 zoom-out)

- **Target:** a unique global physical attractor / state from every admissible initial condition, with no supplied
  measure.
- **Compare:**
  - contractive dissipative dynamics;
  - primitive Markov semigroups;
  - gradient systems with a unique minimum;
  - invertible Hamiltonian / unitary controls.
- **Firewall:** a closed invertible microscopic theory cannot erase fine-grained information. If a unique basin needs
  dissipation, openness or coarse-graining, that price goes into D or A.

---

# S2-H2 PRE-REGISTRATION (owner review of `c063723`; written before any run)

## S2-H2 — state / basin selection by dynamics (attacks H)

**Question.** Can primitive dynamics uniquely determine the physical state or basin from **every** admissible
initial condition, with no supplied measure or preparation class? "Almost every" counts as **MEASURE-PRICED**.

**H2-0 — three notions, kept separate:**

| Notion | Question |
|---|---|
| **H2-A ATTRACTOR UNIQUENESS** | does one asymptotic state / orbit / gauge class exist? |
| **H2-B GLOBAL REACHABILITY** | does every admissible initial state converge to it? |
| **H2-C INFORMATION ERASURE** | are distinct initial states unrecoverable from the *exact* final microstate, including all degrees of freedom of the model (environment / bath where present)? |

A and B do not imply C.

**H2-5 — admissible state spaces, fixed now. No exclusions after the run.** Any later exclusion is priced as
H / Σ / A.

| Model | Dynamics | Admissible states |
|---|---|---|
| H2-1 contraction | F(x) = c·R(x) + b on ℝ², c < 1, R a rotation | all of ℝ² |
| H2-1′ dilation | the same contraction realized unitarily: system + fresh ancilla register (swap-type collision model, qubit) | all system density matrices; ancillas fixed in \|0⟩ |
| H2-2a gradient, unique minimum | V(x) = x⁴/4 + x²/2 − a·x | all of ℝ |
| H2-2b gradient, double well | V = x⁴/4 − x²/2 + εx | all of ℝ, **including** the unstable stationary point |
| H2-2c gradient, symmetric well | V = x⁴/4 − x²/2 | as in H2-2b; outcome tested modulo the gauge x ↔ −x |
| H2-3 Markov | primitive, doubly stochastic 5-state chain; also a generic primitive chain | all probability vectors |
| H2-3′ quantum channel | primitive amplitude-damping + dephasing qubit channel | all density matrices |
| H2-4 unitary hostile | random Hamiltonian on 6 qubits | all pure states |
| H2-7 alignment / frame field | (i) dissipative Kuramoto-type alignment of N planar frames on complete and ring graphs; (ii) Ising frame field under zero-temperature Glauber dynamics | all phase / spin configurations, including twisted and domain configurations |

**Required in every model:**
- H2-6 robustness: small generic perturbations; whether uniqueness survives (otherwise **TUNING-PRICED**);
- H2-8 accounting row;
- H2-9 classification:
  - **H → D RELOCATION** (the attractor value appears as an explicit parameter of the law); or
  - **TRUE H COMPRESSION INTO D** (the attractor is fixed by the symmetry or class of D, with no state-valued
    parameter).

**H2-10 strong targets:**
- consensus / alignment, unique modulo a global gauge;
- the uniform stationary state of a doubly stochastic primitive process;
- the ordered shared reference, unique modulo J ↔ −J.

**Verdict vocabulary:**
- H-SELECTED-FROM-D;
- H-SELECTED-MOD-GAUGE;
- H-MEASURE-PRICED;
- H-TUNING-PRICED;
- H-ACCESS-PRICED;
- NO FINE-GRAINED H SELECTION;
- H → D RELOCATION;
- TRUE H COMPRESSION INTO D.

---

# S2-D-ARROW PRE-REGISTRATION (owner-approved; theorem-first; written before any numerics)

**Primary question.** Can a closed unitary universe produce a nontrivial effective subsystem arrow from **EVERY**
admissible global state, without a low-entropy, low-correlation or otherwise special preparation?

**Secondary question.** Is Σ needed before "the arrow" can even be defined?

Each supplied item is priced if it is needed:
- a low-entanglement / product start;
- low subsystem entropy;
- special correlations;
- a measure / typicality distribution;
- a TPS;
- a coarse-graining.

| Step | Content |
|---|---|
| **D0** | A theorem: no nonconstant continuous arrow functional is monotone along any orbit of finite-dimensional unitary dynamics (recurrence), together with a window version via time reversal. Scope: finite-dimensional closed unitary only |
| **D1** | Every-state hostile for S(ρ_S(t)), using: product; Haar; the time-reversed post-maximum state; an energy eigenstate; a recurrence-near state |
| **D2** | Typicality ≠ arrow: Haar / energy-shell typical states are already near-maximal at t = 0. Measure the available increase and the direction asymmetry |
| **D3 / D6** | Same marginals, different correlations: product vs correlated (time-reversed evolved) states with an equal or near-maximally-mixed bath marginal |
| **D4** | Collision models with pure / thermal / maximally mixed / classically correlated / GHZ-entangled ancillas. Five separate questions: attractor; its value; Markovianity (trace-distance non-increase); entropy arrow; records |
| **D5** | Freshness firewall: fresh vs finite reused bath (M = 1, 2, 4) vs recycled ancilla |
| **D7** | Record arrow (S2-8 diagnostics) vs the antiunitary time-reverse of a record-forming final state |
| **D8** | The same global state and H in inequivalent TPSs (Clifford / Haar frames; **not** commutant frames, which are H-relative gauge after REPAIR 02) |
| **D9** | Access: full state vs local marginal vs coarse macro-variable vs records |
| **D10** | Janus: entropy on both sides of a special middle state; classify BOUNDARY CONDITION ELIMINATED vs TIME ORIENTATION NOT SELECTED; whether a law-fixed special state can carry dynamics |
| **D11** | STRONG (every state) and WEAK (typical: measure, macrostate, TPS, coarse-graining, preparation class) verdicts kept separate |
| **D12** | ARROW_ORIGIN_LEDGER |
| **D13** | Outcomes A–E as given by the owner |
| **D14** | Theorem target, scoped |

**Admissible state space:** all pure (and, for D3 / D4 / D6, mixed) states of the stated finite-dimensional system. No
exclusions.

**Literature** (owner-cited; primary verification attempted, see the literature ledger):
- Popescu–Short–Winter (Nat. Phys. 2006);
- Goldstein–Lebowitz–Tumulka–Zanghì (canonical typicality, PRL 2006);
- Goldstein–Tumulka–Zanghì (PRD 94, 023520);
- Bocchieri–Loinger (quantum recurrence);
- Barbour–Koslowski–Mercati (Janus point).

---

# S2-G2 / G3 / G4 PRE-REGISTRATION — DIMENSION ORIGIN CAMPAIGN (owner-approved; written before any run)

G1 (local Hilbert factor dimension) is already adjudicated: **Σ PRIMITIVE**. "Dimension" stays split. Emergence of one
notion **never** counts as emergence of another.

**Central question:** which dimensions are derived from lower-level dynamics / relations, which are observer- or
probe-dependent, and which remain supplied?

## G2 — graph / spectral dimension

**G2-0 estimators, kept separate:**

| estimator | definition |
|---|---|
| growth d_g | N(r) ~ r^{d_g} (BFS balls) |
| spectral d_s | P_ret(t) = (1/N) Tr e^{−tΔ} ~ t^{−d_s/2}; d_s(t) = −2 d ln P / d ln t |
| walk d_w | ⟨r²(t)⟩ ~ t^{2/d_w} (heat-kernel second moment in graph distance) |
| Hausdorff | = growth dimension for graphs (noted, not separate) |

**Steps:**
- **G2-1** positive controls: chain, square and cubic tori.
- **G2-2** same order of N, different graphs:
  - chain;
  - square;
  - random 3-regular;
  - small-world (Watts–Strogatz);
  - binary tree;
  - Sierpinski gasket (known d_f = 1.585, d_s = 1.365, d_w = 2.322);
  - comb.
- **G2-3** scale hostile:
  - anisotropic 2D (bundled chains, transverse weight ε);
  - layered 3D;
  - small-world perturbation of a ring.
- **G2-4** diffusion-operator hostile, on the same graph:
  - standard vs random-weighted vs anisotropic Laplacian;
  - fractional generator Δ^{α/2} (same graph, different dynamics);
  - long-range Lévy-type rates on a ring.
- **G2-5** structural:
  - d_s is a Laplacian-spectral invariant (blind to every isospectral rearrangement);
  - a same-growth / different-diffusion pair (comb vs square).

**Verdicts:** G2-DERIVED-FROM-D / G2-SCALE-PRICED / G2-PROBE-DYNAMICS-PRICED / G2-NONUNIQUE.

## G3 — spacetime / causal dimension

**G3-0 firewall:** a spacetime claim needs a causal order or propagation cones. A diffusion exponent is never
"spacetime dimension".

**Steps:**
- **G3-1** sprinkled causal intervals in 1+1, 2+1 and 3+1 Minkowski; Myrheim–Meyer ordering-fraction estimator.
- **G3-2** estimator comparison:
  - Myrheim–Meyer;
  - midpoint scaling;
  - longest-chain scaling (N^{1/d});
  - sensitivity to region shape (interval vs slab) and to N.
- **G3-3** same spatial graph dimension, different causal structure: lattice spacetimes over the same 2D spatial lattice
  with Manhattan (L1) vs Chebyshev (L∞) cones, and different propagation speeds; MM estimates compared.
- **G3-4** propagation cones from local quadratic hopping Hamiltonians:
  - chain;
  - square lattice;
  - triangular lattice;
  - chain with next-nearest hopping.

  Measure cone speed per direction and cone volume growth.
- **G3-5** Lorentz hostile: anisotropic cones and non-relativistic dispersion at the same dimension. DIMENSION SELECTED
  and LORENTZ STRUCTURE SELECTED are scored separately.

**Verdicts:** CAUSAL DIMENSION DERIVED FROM ORDER / G3-D-PRICED / ESTIMATOR-PRICED / GRAPH DIMENSION ≠ SPACETIME
DIMENSION / LORENTZ STRUCTURE NOT DERIVED.

## G4 — operational / information capacity

**Definitions, kept separate:**

| symbol | meaning |
|---|---|
| d | Hilbert dimension (where defined) |
| K | affine dimension of the state space |
| N | the maximum number of perfectly distinguishable states (no-restriction effects unless stated) |
| log₂ N | classical bit capacity |

**Steps:**
- **G4-1** controls with matched N and different state spaces: classical simplex, rebit, qubit, gbit (square),
  polygons, spin-factor balls. N is computed by LP over vertex subsets for the polytopes.
- **G4-2** capacity from dynamics:
  - closed classes of a Markov chain;
  - the noiseless-subsystem structure of collective SU(2) noise on 3 qubits (commutant algebra).
- **G4-3** access hostile: the same state space with restricted effect sets (collective-only readout; unsharp effects).
- **G4-4** noise: exact capacity vs ε-distinguishability capacity vs Holevo / channel capacity for a depolarized qubit.
- **G4-5** composition:
  - N_AB vs N_A N_B for classical, real and complex QM;
  - two gbits under min-tensor (local polytope) and max-tensor (boxworld), by LP;
  - K_AB vs K_A K_B (cross-reference SCOUT-1 D3 / S2-2).

**Verdicts:** CAPACITY DERIVED FROM D / CAPACITY A-PRICED / COMPOSITION-PRICED / CAPACITY DOES NOT FIX
REPRESENTATION.

## Critical cross-hostiles (actively sought)

| pattern | test case |
|---|---|
| G2 same, G3 different | L1 vs L∞ cones over the same spatial lattice |
| G3 same, G4 different | the same causal set carrying different local systems |
| G4 same, G1 / representation different | capacity 2 for bit / rebit / qubit / gbit / polygons / spin factors; capacity 2 for an access-restricted qutrit |

The ledger is `ledgers/DIMENSION_ORIGIN_LEDGER.md`. **ZOOM_OUT_07** comes afterwards, with the owner's eight questions.

---

# S2-ΣH PRE-REGISTRATION — JOINT FACTORIZATION / BOUNDARY-CONDITION SELECTOR (owner-approved; first saturation candidate)

**Question.** Can one principle on the closed pair (H, ψ) jointly select **both** a physical factorization Σ and the
independence boundary condition, with no hidden weights, scales, epochs or access assumptions?

| Step | Content |
|---|---|
| **ΣH-0** | Theorem-first: a productness no-go, plus a joint-compatibility proposition (see the result file) |
| **ΣH-1** | Pareto front over (L_H, −C_ψ), not combined. L_H = −k_eff (the S2-Σ locality). C_ψ(Σ, t) = Σ_f S(ρ_f(t)) − S(ρ(t)) (total correlation; for pure ψ, Σ_f S_f) |
| **Candidates** | the 202 groupings of 6 qubit slots (local dimensions 2 … 32, mixed) × frames: identity; 3 random Clifford circuits; a Haar unitary; the commutant frame e^{−iH·0.7}; **the ψ-product frame W_ψ** (a Householder map sending ψ(0) to \|0…0⟩). That is 1414 candidates |
| **ΣH-2** | local dimension variable: the grouping types on fronts and among winners |
| **ΣH-3** | lexicographic rules A, B, C, D as specified by the owner. C and D use k_max, the MDL of H's support closure, and autonomy P |
| **ΣH-4** | ε-constraint: min C s.t. L ≥ L_min, and max L s.t. C ≤ ε, swept over quantiles |
| **ΣH-5** | epoch variable on a symmetric grid t ∈ [−10, 10], step 0.25: C₀ (t = 0, a supplied epoch); C_min (inf over the grid: picks an epoch); C_avg (uniform time average on the grid: a time measure); C_typ (median over the grid: a typicality measure). Kept separate |
| **ΣH-6** | time translation ψ → e^{−iHs}ψ, s = 2.3 |
| **ΣH-7** | time reversal Θ = K for real H; ψ vs Θψ |
| **ΣH-8** | arrow quality R only after ΣH-1 … 7: the fraction of steps with non-decreasing C over a forward horizon h ∈ {2, 5} from the selected epoch |
| **ΣH-9** | records, with one fragment rule for all TPSs: every factor in turn is the system, every other single factor is a fragment; R = the max redundancy count. Comparability across factor numbers is checked |
| **ΣH-10** | MDL = #(H Pauli terms > δ in Σ) + c_ψ · #(ψ amplitudes > δ in a declared product basis of Σ). Varied: the basis language (computational vs local Hadamard), c_ψ ∈ {¼, 1, 4}, δ ∈ {10⁻³, 10⁻⁶, 10⁻⁹} |
| **ΣH-11** | symmetry: quotient candidates by Sym(H, ψ) (translation on the ring model); bare-TPS degeneracy and physical equivalence reported separately |
| **ΣH-12** | positive control: H strongly local (mixed-field Ising chain), ψ product in that same frame |
| **ΣH-13** | hard hostiles: (1) H local in Σ_A, ψ product in Σ_B (a Clifford frame); (2) the ψ-product-many-TPS case (W_ψ); (3) a translation-symmetric ring; (4) a record-forming state; (5) an eigenstate and a Gibbs state; (6) the Janus mean-field state; (7) Haar |
| **ΣH-14** | the ledger `ledgers/JOINT_BOUNDARY_LEDGER.md` |
| **ΣH-15** | success standard A – E as specified by the owner |
| **ΣH-16** | theorem target, scoped |

ZOOM_OUT_08 follows as a **FORMAL SATURATION GATE**, brought to the owner. No GRUT bridge and no automatic freeze.
