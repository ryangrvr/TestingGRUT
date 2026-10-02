# Families E (MaxEnt), F (information-theoretic coarse-graining), G (algorithmic simplicity)

**Frozen witnesses:**
- **W-Hm:** Gibbs at T_b vs a non-Gibbs stationary GGE (BR3-01).
- **W-A2:** coarse-grained sets O₁ / O₂ (E-G7).
- **W-Σ1 / W-Σ2** (family A).
- **W-dim:** modular / geometric patterns of different types (GI-07, GI-09).

## E1 — MaxEnt / minimum cross-entropy (Jaynes; Shore–Johnson)

| C | content |
|---|---|
| C0 | the state is the maximum-entropy (or minimum-cross-entropy) distribution consistent with given expectation constraints |
| C1 | H_marginals |
| C2 | W-Hm |
| C3 | the constraint set; the prior / reference measure (continuous case: an invariant measure m(x)) |
| C4 | Gibbs if the constraint is energy only; GGE if more charges are imposed |
| C5 | Shore–Johnson (abstract): uniqueness **given constraints and a prior**. Which charges are constrained is exactly the Gibbs-vs-GGE choice |
| C7 | twin: World A constrains ⟨H⟩ only, World B constrains all conserved charges. Both are admissible (frozen stationary non-Gibbs family) |
| C8 | vary the constraint set → the state changes (Gibbs ↔ GGE; Ilievski et al. memory-grade: local-only charges incomplete) |
| C12 | **KILLED — S5** (constraint set + prior) |

## E2 — Generalized Gibbs ensemble (Rigol et al.)

| C | content |
|---|---|
| C0 | relaxed state = MaxEnt with all (relevant) conserved charges |
| C3 | the set of charges and their initial values |
| C5 | the GGE "carries more memory of the initial conditions" (abstract), i.e. state-priced |
| C12 | **KILLED — S5** (charge set + initial state) |

## E3 — Maximum caliber / maximum entropy production

| C | content |
|---|---|
| C0 | path distributions by MaxCal; steady states selected by MEP |
| C1 | H_marginals / A_time / orientation |
| C3 | path constraints; the reference path measure; for MEP, the flux constraints |
| C5 | Dewar 2009 (abstract): MEP "is not a physical principle" but MaxEnt inference translating assumptions. Grinstein–Linsker: derivation error. Bruers: min or max depending on constraints |
| C9 | orientation: entropy *production* presupposes the direction of the time parameter in the path measure (the supplied orientation) |
| C12 | **KILLED — S5** (constraints / measure). MEP grade: contested, not a principle |

## F1 — Information bottleneck

| C | content |
|---|---|
| C0 | compress X into T keeping maximal information about a relevance variable Y, at trade-off β |
| C1 | A_resolution |
| C2 | W-A2 |
| C3 | p(x, y); the relevance variable Y; β; the cardinality of T |
| C5 | Y and β are the selecting inputs. The owner rule says a tunable trade-off parameter is not selection |
| C6 | a family of solutions along β; degenerate (Kolchinsky et al.: trivial solutions at every point when Y = f(X)) |
| C12 | **KILLED — S5** (relevance variable + β) |

## F2 — Causal states / ε-machines; Shalizi–Moore macrostates — *the strongest A_resolution candidate*

| C | content |
|---|---|
| C0 | the coarse-graining of histories into **causal states**: the unique minimal maximally predictive partition. Shalizi–Moore: macrostates = "the unique maximal partition" consistent with the observations and Markovian |
| C1 | A_resolution / coarse-graining |
| C2 | W-A2 |
| C3 | a stationary process **over a given observable / alphabet**; the choice to predict the future of that same observable |
| C4 | the unique causal-state partition, given the process |
| C5 | uniqueness is real and theorem-grade *given* the observable channel. That channel is A_interface (A_readout) + A_partition, the very data the frozen record lists as supplied. Shalizi–Moore: "start with a given set of observables" |
| C7 | twin: the same microdynamics observed through two different readouts gives two different causal-state partitions. Both are admissible (frozen W-A1: 112 / 276 seed pairs) |
| C8 | vary the measurement partition → the ε-machine changes |
| C12 | **KILLED — S5** (observable channel = A_interface / A_partition). A **class C reconstruction** (unique given inputs); analogous to the frozen A_closure = f(D, A_seed; R_closure) |

## F3 — IB / RSMI-optimal RG coarse-graining (Koch-Janusz–Ringel; Gordon et al.; Lenggenhager et al.)

| C | content |
|---|---|
| C0 | choose the RG coarse-graining that maximizes the mutual information between block variables and the distant environment |
| C3 | Monte Carlo samples (state); block / buffer / environment geometry; the number of coarse variables; the compression level |
| C5 | geometry and compression level are A_partition / A_resolution inputs. "Optimal" is relative to the RSMI objective |
| C12 | **KILLED — S5** (block / buffer geometry + compression level) |

## F4 — Predictive information (Bialek–Nemenman–Tishby)

| C | content |
|---|---|
| C0 | the divergent part of I_pred is the unique complexity measure (by required properties) |
| C1 | A_resolution |
| C3 | a stationary process and an observation variable |
| C5 | a measure of complexity, not a selector of resolution. Uniqueness comes from axioms on the measure |
| C12 | **KILLED — S2** (does not choose between O₁ / O₂) |

## G1 — Solomonoff / Kolmogorov / MDL simplicity

| C | content |
|---|---|
| C0 | prefer the description (Σ, dimension, law, state) of minimal algorithmic complexity / description length |
| C1 | Σ / dimension / A_resolution |
| C2 | W-Σ1 (both nets have finite descriptions) |
| C3 | the universal machine / description language / coding convention |
| C5 | invariance holds only up to an additive constant (or multiplicative constants for priors, Wood–Sunehag–Hutter). For finite witness pairs that constant can reverse the ranking. Müller 2010 (abstract): the attempt to remove machine dependence via a stationary distribution fails |
| C9 | **S6:** change the universal machine (an equivalent representation) → the preferred member of a finite pair can flip |
| C12 | **KILLED — S6** (convention-priced: machine / language) |

## G2 — "Law without law" (Müller 2020)

| C | content |
|---|---|
| C0 | an algorithmic prior over observer states yields an emergent external world |
| C3 | a universal monotone machine; the encoding of observer states; the current observer state |
| C5 | machine choice (G1) + observer-state encoding (= state / partition) |
| C12 | **KILLED — S6** (machine / encoding convention). Its asymptotic machine-independence claims, if any, were not retrieved: not verdict-bearing |
