# Architecture, regress, selection, and supplied structure

**DRAFT — pending Claude Code / external review.** Issue #5, sections A–F.
This is a bounded audit, not an exhaustive impossibility or novelty search.
Primary-source inspection levels are explicit in [SOURCES.json](SOURCES.json).
The exact mathematical claims are in [PROOFS.md](PROOFS.md).

## A. Ten versions of changing laws

| Architecture | Complete mathematical meaning | Fixed representation | What could still be physically selective? |
|---|---|---|---|
| Fixed state evolution | x′=f(x,a), or supplied kernel | Already fixed | Restrictions in f/kernel |
| Time-dependent parameters | x′=f(x,cₙ,a), cₙ supplied or c′=g(c,t) | Include c and any required clock | A specific g and coupling to x; arbitrary schedules may need unbounded state |
| Stochastic couplings | Joint kernel for (x,c) | Controlled Markov kernel on joint state | Specified noise spectrum, drift, correlations and conservation constraints |
| Branching effective theories | Regime label b plus transition kernel and conditional dynamics | State (x,b), or disjoint union of regime state spaces | Which transitions are allowed, their probabilities and operational observables |
| Law-space evolution | L′=F(L,H), physical update under L | State containing L and sufficient history H | A particular F with physical consequences; “law space” alone imposes none |
| Self-reference/fixed points | Consistency constraints or update on rule descriptions | Fixed constraint problem or interpreter dynamics | Existence, uniqueness, stability and observational selection must be proved |
| Selection/replicators | Population distribution, fitness and mutation/inheritance | Fixed population dynamics | Restrictions from actual fitness/mutation mechanisms; population measure supplied |
| Path-dependent transitions | Next kernel depends on history | Finite memory lift or whole-history state | Independently bounded memory, locality or resources may exclude a restricted comparator |
| RG/coarse-graining | Effective couplings depend on resolution scale | Fixed microscopic theory and coarse-graining transformation | Universality/exponents can be forced conditionally; scale flow is not chronological drift |
| Law/state unification | One joint state with no preferred split into fast state and slow rule | Fixed joint map on that state; no product ontology needed | A compact joint dynamics can be substantive, but absence of a split does not remove the update law |

T1 proves intervention-preserving representation for the specified causal
kernel class. T3 handles unbounded histories and the finite-state caveat.
Neither assumes a physically privileged law/state partition. An emergent,
initial-condition-dependent distinction between slow and fast variables can
be valuable effective physics while the joint dynamics remains fixed.
Smolin's actual matrix recurrence gives an explicit example, not a merely
verbal objection [M01]. The coordinate or state representation may be much
more expensive than the adaptive description; representation equivalence
does not guarantee equal simplicity, efficiency or mechanistic explanation.

## B. The meta-law regress: exact boundary of the result

Once F and its domain are fixed, state augmentation terminates the
**mathematical hierarchy** at one update. If F itself changes according to
G, include its current description and supervisor state. Every specified
finite hierarchy becomes a single joint update. A computable indefinite
hierarchy can be a program/data process under a fixed interpreter. This
does not answer **why that joint update**, interpreter, probability rule or
initial state is physically realized.

No theorem here establishes that a unique law must explain itself, or that
no such explanation can exist. The narrower result is that moving the word
“law” to an evolving register does not supply that explanation. Leaving
novel-event probabilities undefined avoids a complete kernel, but also
leaves predictive likelihoods undefined. A unique fixed point can reduce
freedom only after its constraints, existence and uniqueness are supplied
or proved. Gödel/Turing limits are not blanket prohibitions on deriving
constants or self-consistent physical theories [M06,M07].

| Examined route | Supplied commitments | What it actually offers | Regress/exclusion assessment |
|---|---|---|---|
| Smolin/Unger dilemma and law/state unification [M01] | Matrix domain, recurrence, initial pair, effective timescale interpretation | Explicit joint evolution with approximate state/law split | Pair-state lift is exact; explanatory choice of joint update remains |
| Principle of precedence [M02] | Definition of same preparation/measurement, ensemble of precedents, novelty rule | Proposed historical dependence of quantum experiments | New-event and small-precedent prescription is a substantive obligation; urn analogy is not its derivation |
| Cosmological natural selection [M08,M10] | Universe reproduction, inheritance/mutation, population measure, fitness proxy | Conditional selection, possibly including replication fidelity | Ordinary branching/selection dynamics if fully specified; black-hole fecundity mechanism is not established by analogy |
| Varying-constant fields [M14,M15] | Field, action, coupling functions and initial cosmology | Testable drift/attraction within fixed extended physics | Genuine predictions may exclude frozen-parameter models; not all fixed models |
| Causal-set sequential growth [M11] | Order, growth, covariance/causality axioms and nonnegative coupling sequence | Restricted transition family | Conditional reduction is real; remaining tₙ are not all selected by those axioms |
| Energetic causal sets [M12] | Causal events, momenta, conservation, amplitudes, probability relation | Proposed emergent spacetime from processes | Larger fixed process/constraint system, with substantial physical postulates |
| Histories/path-integral formulations | History domain, action/amplitude weight, composition and boundary conditions | Global correlations and interference | Reformulation alone changes no predictions; nonlocal/global constraints can be substantive but are still specified laws |
| Universal induction/MDL [M06,M07] | Reference language/machine, prior, loss or description criterion | Conditional predictive convergence and compression criteria | Not an observed physical selector; uncomputability and language-dependent finite constants matter |
| RG and universality [M20] | Microscopic degrees, action, coarse-graining and scale | Robust effective restrictions insensitive to some microscopic details | Real conditional freedom reduction; relevant directions and matching data remain; not law change in time |
| Replicator/mutation selection [M10] | Accessible types, fitness and transmission | Convergence or concentration under specified hypotheses | C10 shows absence/tie obstructions; C11 shows rate-dependent stationary concentration |
| Self-organized criticality/attractors [M21] | Driving, dissipation, boundaries, update | Robust organized regimes in particular model classes | Does not select arbitrary measured constants or the dynamics that creates criticality |
| Adaptive networks [M09] | Nodes, edges, rewiring and node rules | Topology/state coevolution and emergent organization | Joint state includes topology; fixed representation exact |
| Bayesian/variational learning | Hypothesis space, likelihood, prior, updating objective | Data-dependent parameter inference | Learning a constant is not nature deriving it; treating posterior as a physical state requires additional dynamics |
| Fixed-point/self-consistency | Supplied constraint/map and domain | Unique solution only if an existence/uniqueness theorem applies | C08 versus identity/translation maps separates real conditional selection from an empty slogan |

These assessments are inference from the stated constructions and inspected
sources. Abstract-only sources establish what authors propose; their full
theorems or empirical claims are not independently certified by this audit.

## C. Have any requested physical constants been selected?

Within this bounded pass, **no adaptive-selector mechanism is established
that uniquely computes the measured joint set** of constants and structures
below without further physical inputs. This is NOT a proof that no such
mechanism exists. Genuine conditional attractors and reconstruction
theorems must be retained at their proper scope.

| Requested quantity | Strong examined mechanism | Category earned | Remaining supplied/fitted information |
|---|---|---|---|
| Fine-structure constant | Dilaton/least-coupling attraction [M15] | Dynamically attracted under specified mass/coupling functions and basin assumptions | Extremum position and coupling value are not uniquely derived; the rational C08 target 1/137 is deliberately supplied and is not the measured value |
| Mass ratios | Standard Model Yukawa sector; varying scalar-dependent masses [M14,M27] | Fitted parameters or conditional field-dependent relations | Matter content, Yukawa data, strong scale and matching conditions |
| Cosmological constant | Flux discretuum [M16]; sequestering [M22] | Scanned/environmentally conditioned, or dynamically constrained under a specified modification | Charges/flux sectors, measure, global constraints and residual cosmological data; no exact measured Λ derived here |
| Dimensionality | Causal triangulation effective geometry [M26] | Emergent effective behavior within a defined ensemble | Microscopic simplices, allowed topology/action and ensemble are inputs; not a unique structure-free dimension selector |
| Gauge group/matter | Consistency/anomaly constraints on supplied representations; reviewed in PR #4 | Conditional restriction | Candidate representation set and charge normalization; consistency alone is not evidence for a unique measured group/content |
| Symmetry-breaking scales | Higgs potential or dimensional transmutation [M27] | Dynamics yields a scale conditional on parameters and matching | Potential parameters or reference coupling/scale; eliminating one scale in favor of a coupling is not eliminating all input information |
| Born structure | Operational reconstruction [M17]; pilot-wave relaxation [M18] | Derived conditional on axioms; numerical/coarse-grained attraction in a supplied model | Preparation/composition axioms or wave/guidance dynamics, initial distribution and coarse graining; no universally proved selection from zero structure |
| Arrow-of-time boundary | Low-entropy boundary, typicality/ensemble, irreversible effective descriptions | Boundary/measure conditioned in the examined descriptions | Time orientation or special preparation remains an obligation; memory and irreversible dynamics cannot be silently inserted |

The selection controls give exact positive and negative results. C08 has
initial-condition-independent convergence, with the target encoded in its
map. C09 has multiple basins. C10 cannot create an absent fitter type and
does not select under equal fitness. C11's stationary distribution depends
on specified mutation probabilities. T6 shows a standard historical
reinforcement process whose asymptotic parameter is distributed rather
than a unique number. None is promoted into a physical candidate.

## D. Replace “best path” with priced physical content

“Observer” below means an external agent is needed to *define* a criterion,
not merely that experimenters are needed to measure a prediction. All
criteria need a supplied mathematical domain and operational calibration.

| Criterion | Local/global; temporal character | External observer required? | Selects one branch/constant? | Measurable consequence and hidden premise |
|---|---|---|---|---|
| Stationary action | Global trajectory constraint, often time symmetric | No | Usually many stationary paths; boundary data needed | Trajectory/phase relations; action and boundaries supplied; stationary is not always minimum |
| Entropy production | Local rates or integrated trajectory; asymmetric in effective ensembles | No | No universal maximization rule established here | Dissipation/fluxes; reservoir, entropy convention, driving and allowed variations supplied |
| Thermodynamic viability/free energy | Ensemble/local or global; relaxation can be asymmetric | No for physical free energy; yes for an arbitrary design score | Depends on phase, constraints, temperature and domain | Work/heat/stationary distributions; thermostat/ensemble cannot be unexplained inputs |
| Variational inference/free energy | Global statistical objective with sequential updates | Agent/model selection normally supplied | An optimum need not be unique or true physics | Predictive losses; generative model, likelihood and priors supplied |
| Stability/persistence | Local linear or global basin property; dynamical time oriented | No | Often many attractors, including trivial ones | Decay/growth rates and basin tests; dynamics, perturbation norm and persistence threshold supplied |
| Compression | Global coding/model criterion; no inherent time arrow | Code/user choice supplied | Many tied codes/models; shortest programs need not be computable | Held-out code length; fixed language/precision and penalization supplied |
| Predictive sufficiency | Conditional distributions for future records; time oriented | Interface/protocol definition supplied | Identifies predictive equivalence classes, not unique microphysics | Held-out predictions; accessible readouts, interventions and memory scope supplied |
| Causal consistency | Local or global; order presupposed | No | Excludes inconsistent processes but admits many | Forbidden signalling/cycles; causal domain and compositional axioms supplied |
| Universe reproduction | Population/global generations; asymmetric | No if physically defined | Conditional, with mutation/fitness tradeoffs | Reproduction statistics would be needed; fecundity, inheritance and measure supplied |
| Robustness | Local perturbation or global ensemble; either temporal symmetry | No | Depends on perturbation repertoire and objective | Sensitivity bounds; which perturbations/counting measure supplied |
| Computational depth/complexity | Global program/history; time orientation by computation | Machine/encoding supplied | No unique physical constants follow | Resource predictions require a physical computer model; universal description does not identify nature's program |
| Constraint/anomaly satisfaction | Global algebraic consistency; no inherent arrow | No | Restricts supplied sectors, generally not unique | Forbidden amplitudes/charges; representations, gauge content and admissible baseline supplied |
| Decoherence/branch stability | Open-system/dynamical; effectively asymmetric | No consciousness premise needed | Suppresses interference conditionally; does not uniquely select an outcome or constants by itself | Interference/record correlations; system/environment split, interaction and initial state supplied |

C14 demonstrates that changing the supplied objective changes the selected
branch on an unchanged feasible set. A physically derived objective can be
explanatory; naming an optimum without its physical derivation cannot be.

## E. What “0→1” already needs

| Claim | Minimal structure required to state/test it |
|---|---|
| No distinction | A domain and a relation that identifies its elements, or a designated empty object |
| First distinction | At least two distinguishable alternatives, or a grammar for creating distinct objects and an equality rule |
| Relation | Objects and a relation/arrow with incidence conditions |
| Ordering | An order relation and its axioms |
| Causality | An operational precedence/influence rule, not order vocabulary alone |
| Probability | An event algebra and normalized measure or kernel |
| Composition | A specified operation, typing and associativity/compatibility conditions |
| Persistence | Successive states/relations, an identity criterion, an interval and dynamics |

Boolean/distinction calculi supply syntax and operations. Category/process
frameworks supply objects/arrows/composition; constructor theory supplies
tasks and possible/impossible transformations [M13]. Causal sets supply
order and growth rules [M11]; energetic versions add momenta/amplitudes
[M12]. Graph rewriting and cellular automata supply grammars/neighborhoods
and updates. Quantum automata supply local quantum algebras and unitary
causal updates. Spin/tensor networks supply graph, labels/tensors and
contraction; GPTs supply state/effect spaces and composition. Algorithmic
generation supplies a machine and input. None of these descriptions
literally has no supplied structure. This is an input audit, not a claim
to refute every framework's independent physical results.

T7 gives the narrow exact obstruction: literal ∅ cannot carry a normalized
initial probability. An empty graph, vacuum state or tensor unit is a
well-defined state of an already structured theory. A supplied birth rule
is ordinary dynamics; it does not derive probability, ordering or the
meaning of a first distinction.

## F. Historical memory and changing effective laws

| Apparent law change | Conventional fixed-law mechanism | What observation alone establishes |
|---|---|---|
| Hysteresis after equal visible reset | Retained internal memory, domains, defects or slow variables | Visible state is insufficient; full-state reset not certified |
| Branch-dependent constants | Scalar background or different vacuum/sector | Different effective parameters; not changing fundamental update |
| Irreversible symmetry choice | Symmetry-breaking dynamics and boundary/initial data | A realized basin/phase, not a selected universal constant |
| Topological/metastable memory | Conserved sector, barrier or long relaxation time | Restricted accessible transitions given state and timescale |
| Nonergodic parameter frequencies | Quenched disorder, mixtures, conserved sectors | Ergodic assumption failed; no universal law evolution inferred |
| Repetition-dependent response | Reinforcement, apparatus learning, environment aging | A mechanism may be testable through intervention; T6 limits passive identification |
| Running coupling | RG resolution dependence | Effective-theory scale dependence; compare at a common physical resolution before inferring chronological drift |

This architecture can describe history without teleology. It has not earned
a distinct fundamental physical prediction merely by containing history.
