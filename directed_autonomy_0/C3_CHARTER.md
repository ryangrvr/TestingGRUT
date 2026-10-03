# DA0 · C3 CHARTER — REFLEXIVE ARCHITECTURE GENERATION AND INTERNAL-MODEL CLOSURE

> **C3 CHARTER REPAIR 01** (C3R-01 … 07; `DA0_CORRECTION_LEDGER.md`).
> - **Owner review:** CONDITIONAL PASS. **Primary model ACCEPTED** (local adaptive-coupling Ising). **DO NOT RUN YET.**
> - **Pre-repair charter boundary:** `f266be2f6432e4c8de6c73a183ea5135fe816a8e`.
> - Where wording differs, the repaired text and the ledger take precedence.

**Status:** CHARTER ONLY — **NOT RUN; AWAITING OWNER REVIEW.** No C3 code exists and no C3 experiment has been run.

**Base.** `grut-directed-autonomy-0 @ 6587169`.
- Reviewed C2 boundary: `99428ff7462645fb788840cbab979ae737459410`.
- RA0-frozen: `ab4fd860ebbc72a6a3c0c5dd1ff2cef52616b7f5`.

**Not modified:** C1 and C2 scientific claims.

**Prohibitions:**
- no consciousness interpretation;
- no quantum, collapse or Born-rule work;
- no PR or merge.

**Labels (as in DA0):** KNOWN / REDERIVED, PROVED HERE (internally), NUMERICAL ILLUSTRATION, CONJECTURE, INTERPRETATION.
**Rule:** hypotheses may be introduced aggressively; their information price may not be hidden.

## 1. The surviving obstruction and why C3 joins two questions

**From C2.** Blind recovery tolerates K_N → ∞ (Theorem C2-A1a), but **every growing architecture recovered was supplied by
its family**. C2 found no fixed simple rule that **generates** growing collective architecture.

**Reserved in the DA0 charter (§2.4).** C3 was reserved for genuinely reflexive internal modelling, together with the
hostile question: *can finite recursive closure exist without an externally supplied hierarchy or stopping depth?*

**C3 joins them.** Can feedback of a system onto its own effective dynamics **generate** the architecture C2 could only
**recover**? And can the result support internal-model closure **without** a supplied model register, factorisation or
meta-depth?

## 2. Central hypothesis and the fixed-law form

**Allowed form.** An augmented physical state Z = (X, W) with one fixed autonomous law:

  dX/dt = F(X, W),  dW/dt = G(X, W)   (or the stochastic / jump-process analogue).

- X is the ordinary physical state.
- W is a **local** interaction / memory / bond / material / coupling degree of freedom.
- The effective X-dynamics evolves because W evolves. **The complete Z-law is fixed.** This is not described as "changing
  the laws of physics".

**Priced working-theory primitives** (assumed, not derived):
- the existence and local state space of W;
- the fixed update law G;
- the microscopic locality / substrate;
- the fixed microscopic constants;
- any external nonequilibrium drive or reservoir.

**Compression criterion.** C3 earns information-accounting compression only if this **bounded** description generates
architecture whose description **grows with N**, without inserting that architecture into the initial conditions, labels,
coupling schedule or parameter scaling. A new primitive is not a defeat. A hidden one is.

**C3-F0 (substrate split is supplied).** The split Z = (X, W) is part of the supplied substrate. **Identifying W as
"the model" and X as "the world" is relabelling (C2-F6), not an internal-model result** (see §5).

## 3. C3-A — endogenous architecture generation

**Question.** Can a homogeneous local adaptive law, with every microscopic parameter fixed in N, generate growing
**collective** metastable or dynamically autonomous architecture?

### 3.1 Mandatory entry requirements (primary model)

1. a fixed local microscopic rule;
2. a bounded description, independent of N;
3. parameters fixed in N;
4. no tree, level, basin, module, block or hierarchy labels;
5. no N-dependent coupling schedule;
6. no supplied macro-units or recovery grains;
7. no clustering, crispness or partition optimisation;
8. (X, W) is a well-defined physical stochastic or deterministic dynamical system;
9. if X alone is non-Markovian, Z = (X, W) is Markov, with an explicit state space;
10. the **thermodynamic / energetic cost** of maintaining and adapting W is stated, not ignored.

### 3.2 Collective-architecture firewall (C3-F1; strengthens C2-F6)

**Not a positive.** None of the following yields a C3-A positive:
- one independent memory bit per site;
- one frozen bond label per edge;
- isolated defects;
- microscopic-site relabellings;
- independent product factors whose number merely scales with N.

**[C3R-04] The pairwise criterion below is WITHDRAWN** (kept as a record). It fails its own independent-bit control.
- **The counterexample:** two random macrostates of N independent bits differ on ~N/2 bits, so PR ~ N/2 grows linearly.
  The median pairwise PR would therefore certify the product-bit architecture this firewall forbids.
- **The replacement is §3.2a**, which measures collectivity on elementary transitions.

**Original (withdrawn) collective criterion.** For each pair of recovered macrostates a, b, let δ(a, b) be the profile of
conditional mean differences of the local variables over the microscopic sites. Its **participation number** is

  PR(a, b) = (Σ_i δ_i²)² / Σ_i δ_i⁴,

which is basis-free with respect to relabelling sites. The criterion has two parts:
- **Collective:** the median PR over recovered pairs must **grow with N**. The growth is judged by the sign of a fitted
  scaling exponent across sizes. This is a diagnostic: no threshold is chosen post hoc, and none enters any definition.
- **Bookkeeping:** structures whose PR stays bounded in N are microscopic bookkeeping, graded **C3-A3**.

**Growth requirement.** A positive needs **both**:
- the number of recovered collective structures grows with N, **or** a filtration of collective structures deepens;
- per-structure participation grows with N.

### 3.2a Collective criterion on elementary transitions [C3R-04] — replaces the pairwise criterion

**Edges.** Elementary transitions a ↔ b between symmetry-inequivalent candidate macrostates, defined canonically in §3.5′
(lowest-barrier saddle / minimum-action connection). No arbitrary pair is used.

**Two conditions on each relevant edge,** both using **gauge-invariant** local observables O_i (C3-F2), with
δ_i^{ab} = ⟨O_i⟩_a − ⟨O_i⟩_b:

1. **Participation.**
   PR_edge(a, b) = (Σ_i (δ_i^{ab})²)² / Σ_i (δ_i^{ab})⁴ must **grow with N**.
   Independent product bits flip one bit per elementary edge, so PR_edge = O(1) and they **correctly fail**.
2. **Non-vanishing contrast.** With
   S_ab = Σ_i (δ_i^{ab})² / Σ_i [Var_a(O_i) + Var_b(O_i)],
   require **liminf_{N→∞} S_ab > 0**, or an equivalent asymptotic signal-to-noise condition. This is a limit, not a
   fitted threshold. It stops O(N) infinitesimal changes from faking collectivity.

### 3.2b C3-F2 — exact-symmetry / gauge-orbit firewall [C3R-03]

**The gauge symmetry.** The local law has the exact Ising gauge symmetry σ_i → η_iσ_i, J_ij → η_iη_jJ_ij (η_i = ±1).
Because J is dynamical, a single adaptive minimum can have **exponentially many gauge copies**.

**Counting rule.** **Architecture is counted modulo every exact symmetry of the microscopic law:** local gauge,
translations, lattice automorphisms and global flip.
- **Gauge-invariant observables** are used for the primary audit, e.g. J_ijσ_iσ_j, J_ij², plaquette products
  Π_{b∈∂p} J_b, and loop / frustration quantities.
- **No arbitrary gauge fixing may generate K-growth.**
- The question is whether the adaptive law generates **gauge-invariant frustration architecture**, not random-looking
  gauge copies of bond signs. Gauge-invariant frustration, not gauge copies, carries spin-glass structure.

### 3.3 Candidate audit (before model selection; metadata and abstracts only)

For each class: is growing architecture **generated**, or encoded in microscopic types, topology, stored templates,
level-dependent rules, initial conditions, external drive or a supplied objective?

| class | representative | verdict |
|---|---|---|
| 1. adaptive / coevolutionary networks | coevolving voter model (Holme–Newman 2006): three phases, including an absorbing **fragmented** phase | **encoded / not growing.** The number of opinions is G = N/γ, i.e. the **state alphabet grows with N** (supplied types). Rewiring to random nodes is non-local. Binary-opinion variants fragment into O(1) components. The final state is absorbing (a frozen topology), not metastable architecture |
| 1′. adaptive networks: SOC by local rewiring | Bornholdt–Rohlf 2000: quiet nodes gain links, active nodes lose them, converging to critical connectivity | **partially generated:** a global critical state from a local rule. But the output is a scalar statistic (mean connectivity → K_c), not a recoverable collective macro-partition. The network is asymmetric and non-spatial (topology supplied as a random graph). Weak fit |
| 1″. **adaptive couplings with time-scale separation** | Coolen–Penney–Sherrington (PRB 48, 1993): Ising spins (fast) with Hebbian-type slow interactions at a **different temperature**; adiabatic elimination gives a replicated equilibrium with replica number **n = T_fast/T_slow**, set physically | **genuinely generative candidate.** One fixed law generates a coupling field whose statistics are not supplied. [C3R-01 corrected] At **n = 1** the joint dynamics satisfies detailed balance. Integrating out J leaves the σ-marginal uniform, **but the joint and J-marginal laws retain bond / loop correlations**: J_b | σ ~ N(s_b/μ, T/μ), and ⟨Π_{∂p} J_b⟩ = 1/μ⁴ on a plaquette. The **CPS mean-field model has an ordered q > 0 phase for n ≤ 2, including n = 1.** So n = 1 is an **equilibrium comparator, not a no-architecture null**. [C3R-07] The known theory is **mean-field / infinite-range**; the local lattice version is a **new construction, not validated by CPS** |
| 2. local physical memory / self-interacting dynamics | reinforced random walks; self-interacting diffusions (aging, long memory) | **not growing collective architecture:** localisation on a few sites (site-type structure, C3-F1), or aging without a recoverable partition |
| 3. SOC / driven dissipative | sandpiles; driven self-assembly | **drive-encoded, or transient.** Avalanches are scale-free but form no persistent macrostates. Architecture in driven assembly is typically templated by the drive |
| 4. reaction–diffusion / active matter | Turing patterns; motility-induced phase separation | **K grows, participation does not:** the domain count ∝ N/ℓ*^d, with ℓ* fixed by parameters, so per-domain participation stays bounded. **Fails C3-F1.** Coarsening gives the opposite trend |
| 5. algorithmic / hierarchical self-assembly | tile assembly; programmed hierarchies | **encoded** in tile types / level rules (C2-F3′). Supplied-architecture comparator only |

### 3.4 Proposed primary C3-A construction (structural grounds; subject to owner review)

**PROPOSED: LOCAL TWO-TEMPERATURE ADAPTIVE-COUPLING ISING MODEL** (a CPS-type law on a lattice).

**Rule:**
- spins σ_i = ±1 on a d-dimensional hypercubic lattice with periodic boundaries;
- **one real coupling J_ij per nearest-neighbour bond** (W);
- **fast:** heat-bath Glauber for σ at temperature T, with H(σ | J) = −Σ J_ij σ_i σ_j;
- **slow:** Langevin for each J_ij at temperature T′:
  dJ_ij = ε(σ_iσ_j − μ J_ij) dt + √(2εT′) dB_ij, with ε → small being the time-scale separation;
- **fixed constants:** (d, T, T′, μ, ε), independent of N;
- **homogeneous initial ensemble** (e.g. J ≡ 0 or i.i.d. small); **no labels**.

**Why it is chosen** (structural, before any computation):
1. It is the only audited class where the generator of disorder or architecture is **itself the dynamics** of a local
   field, with no types, templates, level rules, opinion alphabet or drive pattern.
2. Z = (σ, J) is an explicit Markov process (requirement 9).
3. **Thermodynamics is explicit (requirement 10)** [C3R-01 corrected].
   - n = T/T′ = 1 is the **equal-temperature equilibrium comparator** (A1′). It is **not** a no-architecture null.
   - For T ≠ T′, **two reservoirs are supplied (priced)**. At finite ε the entropy production / heat flow is
     **calculated or measured**, not assumed nonzero by wording.
   - The strict adiabatic reduced J-process has an effective equilibrium description (§3.5′), even though the underlying
     two-temperature construction is physically different.
4. **It connects to the C2 boundary.** It asks whether the quenched disorder that C2-E2 would have needed can be
   **generated** by one fixed adaptive law rather than supplied. It also tests whether that disorder carries growing
   collective metastable architecture.

**Registered risks (expected hostile outcomes):**
- The local model may only produce domain structure with a parameter-fixed size (**C3-A3 / A2**).
- It may reduce to an effective EA-like problem whose architecture is the open droplet-versus-RSB question
  (**C3 INDETERMINATE**).
- In the adiabatic limit, J may simply track ⟨σ_iσ_j⟩ (bond bookkeeping, C3-F1).

**The primary model is a proposal.** It is fixed only on owner approval of this charter.

### 3.5′ Primary analytic object and resolution of C3-M1 [C3R-02, C3R-05] — supersedes §3.5 where they differ

**Adiabatic landscape** (exact in the ε → 0 limit; to be derived in the preflight).

  **𝓕_T(J) = (μ/2) Σ_b J_b² − T log Z_T(J)**, where Z_T(J) = Σ_σ e^{(1/T) Σ_b J_b s_b} and s_b = σ_iσ_j.

- Since ∂_{J_b} log Z_T = ⟨s_b⟩_J / T, the adiabatic drift is ε(⟨s_b⟩_J − μJ_b). So dJ/dt = −ε∇_J𝓕_T(J) + slow noise
  √(2εT′).
- The stationary law of the reduced process is ∝ e^{−𝓕_T/T′} = Z_T^{T/T′} e^{−(μ/2T′)ΣJ²} (CPS-type).
- **The stationary points of the deterministic adaptive flow depend on T and μ only, not on T′.** T′ changes the
  slow-noise weighting and the transitions.

**Hessian (exact):** H_{bb′} = μδ_{bb′} − (1/T)·Cov_J(s_b, s_{b′}). Stability is lost when the covariance operator
overwhelms μI.

**Homogeneous state.** At J = 0 the spins are independent, so ⟨s_b⟩ ≈ J_b/T and dJ_b/dt ≈ ε(1/T − μ)J_b. This gives
the **exact linear instability T_lin = 1/μ** of J = 0.

**C3-M1 is resolved, model-specifically.** The primary detector is the **canonical adaptive-landscape / communication
structure**:

| element | definition |
|---|---|
| nodes | symmetry-inequivalent (C3-F2) **stable minima of 𝓕_T** |
| elementary relations | canonical **minimum-barrier / minimum-action saddles** |
| metastability | the barrier ΔF_ab (or the action / exit time) **diverges** in the appropriate thermodynamic scaling |
| growth | the number of symmetry-inequivalent **collective** minima grows with N, **or** a collective communication filtration deepens |

**Forbidden as the primary certificate:** generic clustering, PCCA, a chosen K, chosen lag-time partitions, or a fitted
macrostate count. Simulation may **approximate** these preregistered objects; it may not **define** them by clustering.

**C3-A1 then requires all three:**
1. a growing number or depth of symmetry-inequivalent structures;
2. collective elementary transitions (§3.2a);
3. growing dynamical stability / barriers.

### 3.5″ Analytic preflight before any simulation [C3R-06]

The first C3-A calculation, **once separately approved**, is analytic only:
1. derive the local adiabatic landscape 𝓕_T;
2. derive its exact symmetries and the gauge quotient;
3. derive the stability of the homogeneous state (T_lin = 1/μ and beyond);
4. carry the high-temperature / loop expansion of log Z_T(J) far enough to find the **first coupling between adaptive
   bonds**. Bond-local terms come first; genuinely collective, gauge-invariant terms enter through closed loops, starting
   at plaquettes on the hypercubic lattice;
5. decide whether the first generated structure is **bond-local bookkeeping** or **genuinely loop / collective**.

**No large-lattice simulation until this preflight is reviewed.**

### 3.5 Method firewall (C3-M) — scalable, no 2^N diagonalisation as the primary method [superseded in part by §3.5′]

**Proposed methods:**
1. **Analytic adiabatic limit (ε → 0):** derive the effective stationary law of J. CPS-type: weight ∝ Z_T(J)^n·e^{−β′·(μ/2)ΣJ²}.
   State exactly what is KNOWN (mean-field) and what is REDERIVED or new on the lattice.
2. **Large-scale stochastic simulation** of Z on lattices up to the feasible L, at fixed constants, with trajectory-based
   metastability diagnostics (residence times, return statistics). **Every finite threshold is declared diagnostic-only.**
3. **Architecture detection without clustering optimisation.** **OPEN METHOD ITEM (C3-M1), to be resolved and
   preregistered before any run:** a canonical, scalable macro-architecture detector compatible with C0 / C2-F rules.
   - Candidates include empirical transfer operators on a basis-free observable family, and the RA0 / C2 recovery map
     applied to a certified slow space.
   - Two-replica overlap distributions may be used **as a priced diagnostic only**: the overlap is a chosen observable,
     not an endogenous partition.

**No K, ε (cut), tolerance, depth or module count enters any physical definition.**

### 3.6 Required hostile controls

| ID | control | purpose |
|---|---|---|
| **A0** | **frozen interaction:** J held at its initial value | is the adaptive feedback responsible for any architecture? |
| **A1** | **decoupled memory:** J given the same marginal noise / relaxation but no σ-dependent drift (or σ given no J-dependence) | separates genuine reflexive coupling from merely adding slow variables |
| **A1′** | **EQUILIBRIUM ADAPTIVE-COUPLING COMPARATOR:** n = T/T′ = 1 [renamed per C3R-01] | distinguishes structure that needs reservoir mismatch from structure generated by adaptive coupling itself. It is **not** a null. A0 and A1 remain the actual null controls |
| **A2** | **supplied-architecture positive control:** explicit modular / hierarchical couplings (P1-type) | validates the detector. Graded **ARCHITECTURE SUPPLIED**, never a C3 success |
| **K1 – K4** | permanent DA0 comparators | as chartered |

Controls are not multiplied without a stated purpose.

## 4. C3-B — genuinely reflexive internal modelling (only after C3-A1)

**Gate.** C3-B can earn a positive **only if C3-A has produced endogenous collective architecture**. No variable is named
"self", "model", "observer", "agent" or "world".

**Question.** Does the endogenous architecture contain a **dynamically identifiable** structure that does all three of
the following?
1. **Representation / prediction:** it carries information about the future of other endogenous structure, not
   reducible to a K1 memory bit.
2. **Update:** it is modified by the dynamics whose future it represents.
3. **Causal efficacy:** intervening on it, with the relevant remaining state fixed, changes the subsequent endogenous
   dynamics.

**Firewalls:**
- None of 1 – 3 may rest on a **supplied** decomposition into model and modelled variables.
- **C3-F0 applies:** the supplied X / W split may **not** serve as the model / world split. CPS-type couplings trivially
  "learn" spin correlations; that is the substrate, not a result.
- If evaluating 1 – 3 requires a factorisation, subsystem boundary, observable channel or model register, the terminal
  is **REFLEXIVE CLOSURE REQUIRES SUPPLIED FACTORISATION / MODEL REGISTER**.
- No optimisation over candidate subjects or subsets (C0-P2).
- **K2:** rediscovering predictive equivalence classes is reported as K2, not as a new result.

## 5. Recursive-closure hostile test (retained from DA0 §2.4)

**Test.** Can a finite internal model close on itself **without**:
- an infinite model-of-model hierarchy;
- an externally specified recursion depth;
- a chosen stopping rule;
- a supplied meta-level?

**Possible positive forms:** a genuine fixed point, a finite cycle, algebraic closure, or another invariant **generated
by the physical dynamics**.

**If every construction needs a chosen meta-depth,** that is recorded as the result. Nothing is engineered around it.
G1's lesson stands: exact finite closure under a predictor-of-the-future map collapses to lumpability / bisimulation.

## 6. Preregistered terminals

| ID | terminal | meaning |
|---|---|---|
| **C3-A1** | **ENDOGENOUS GROWING COLLECTIVE ARCHITECTURE — CONSTRUCTED** | a fixed bounded local adaptive rule generates growing collective architecture at fixed parameters (C3-F1 passed). A construction result, not a GRUT prediction |
| **C3-A2** | **ADAPTIVE DYNAMICS PRESENT — NO GROWING COLLECTIVE ARCHITECTURE FOUND** | |
| **C3-A3** | **APPARENT GROWTH IS MICROSCOPIC BOOKKEEPING** | growth only through site / bond / local-memory labels (PR bounded) |
| **C3-B1** | **FINITE REFLEXIVE INTERNAL-MODEL CLOSURE — CONSTRUCTED** | needs C3-A1, closure beyond K1 – K4, and no supplied self / world factorisation |
| **C3-B2** | **ENDOGENOUS ARCHITECTURE WITHOUT ENDOGENOUS MODEL / BOUNDARY** | |
| **C3-B3** | **REFLEXIVE CLOSURE REQUIRES SUPPLIED FACTORISATION / MODEL REGISTER / META-DEPTH** | |
| | **C3 INDETERMINATE** | only for a clearly stated technical obstruction (e.g. C3-M1 unresolved, or the droplet / RSB question undecidable at accessible sizes) |

## 7. Information accounting (every positive)

**Derived:** whatever architecture, hierarchy, model relation or closure actually follows from the fixed law.

**Supplied:**
- the substrate and lattice;
- the existence of W, and the X / W split;
- the update law G;
- the constants (d, T, T′, μ, ε);
- the reservoirs (if T ≠ T′). The heat current / entropy production is calculated, not assumed (C3R-01);
- the initial ensemble;
- any regularity assumptions.

**Not earned, unless actually derived:**
- generic subsystem factorisation;
- consciousness, awareness or selfhood;
- quantum outcome selection or Born weights;
- a universal architecture selector;
- TRUE COMPRESSION relative to frozen GRUT.

**Compression arithmetic.** If a bounded adaptive law replaces architecture whose description grows with N, C3 computes
the information-accounting improvement explicitly. It is **not** called TRUE COMPRESSION unless an actual frozen GRUT
residual is eliminated. A successful small adaptive law is retained as a **candidate working-theory constitutive law**,
for later audit.

## 8. Hard stop (charter stage)

**Charter only.** Before any C3 run, the owner reviews:
- the proposed primary model (§3.4);
- the open method item C3-M1;
- the collective criterion (§3.2).

**Sources consulted** (metadata / abstracts; full texts not re-read):
- Coolen, Penney & Sherrington, *Phys. Rev. B* 48, 16116 (1993); also the NeurIPS 1993 companion, *Coupled dynamics
  of fast neurons and slow interactions*. [C3R-07] **CPS motivates the mechanism; it does not validate the local lattice
  theory.** The analysed CPS model is infinite-range;
- Holme & Newman, *Phys. Rev. E* 74, 056108 (2006);
- Bornholdt & Rohlf, *Phys. Rev. Lett.* 84, 6114 (2000);
- Gross & Blasius, adaptive coevolutionary networks review (arXiv:0709.1858).
