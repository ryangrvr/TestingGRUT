# Comparator classes and standard exclusions

**DRAFT — external review pending.** This is a bounded primary-source audit
and exact standard-control package. Inspection depths and failed full-text
access are recorded in [SOURCES.json](SOURCES.json). It is not an exhaustive
QFT/gravity consistency proof or a new physical selector.

## A. The hierarchy is a partial order, with independent constraint axes

Inclusions below mean operational realization after declared resources and
coarse graining, not equality of primitive theories. Fixed elapsed time,
readout, intervention location and resource budgets cannot be changed to
manufacture an inclusion. "Finite" always needs a quantitative bound.

| Class | Additional premises | What the fixed lift establishes / does not establish |
|---|---|---|
| C0 arbitrary causal stochastic kernels | Normalization, nonnegativity, past-dependent settings; arbitrary measurable retained state | PR #6's architecture reduction is exact here; no physical realization follows |
| C1 finite retained classical state | At most n jointly retained states; fixed operational instruments | Eight-state and four-state controls have finite lifts. The count-state urn and unbounded clock do not have uniform finite-resource lifts established for all horizons |
| C2 computably specified kernels | Algorithm for amplitudes/probabilities to declared precision and effective operations | Orthogonal to locality/thermality. A finite table containing arbitrary noncomputable reals is not a finite algorithm. Computability is a separately stated premise, not a theorem here about all nature |
| C3 local finite-range classical updates | Frozen lattice/graph, finite radius, fresh settings independent of initial world; finite step count | Exact causal cones; reversible finite-horizon control with scratch and schedule works (T3,T8). A global update table is not automatically local at the same latency |
| C4 finite positive CTMC | Q_ij≥0, Q1=0; bounded exit rates, continuous physical time; fixed readout | Exact finite-time NOT is excluded even with finite hidden memory under T6's no-jump premises. This is narrower than all classical or quantum dynamics |
| C5 stationary classical detailed balance | Stationary μ and μP=P^Tμ, passive time-even common readout, no drive | Oriented pair current excluded for any hidden dimension (T5). Not equivalent to Gibbs preservation, equilibrium with odd time-reversal variables, or all thermal operations |
| C6 closed classical Hamiltonian | Specified symplectic phase space, Hamiltonian flow, initial ensemble; Liouville measure conservation | Injective dynamics cannot discard distinguishable information without another degree of freedom. Finite reversible control is a logic-level construction, not a proved smooth autonomous Hamiltonian compilation with fixed resources |
| C7 finite-dimensional quantum instruments | Linear CP outcome maps, sum TP, supplied quantum state and setting/record semantics | Every finite C1 kernel has an exact diagonal realization (T1). Quantum dimension is not interchangeable with classical hidden-state count; T2 bounds perfect decoding |
| C8 closed finite-dimensional unitary | System plus all environment/controller, unitary evolution, initialized ancillas, energy/time constraints if invoked | Dilates C7 for a finite stroke with enough ancilla. Exact irreversible resets accumulate distinguishable information; finite faithful thermal resources cannot produce exact pure reset (T3) |
| C9 local quantum lattice dynamics | Tensor geometry, local interaction terms/norms, finite propagation envelope or circuit depth | Strict cones for circuits; LR tails for continuous Hamiltonians (T8). A generic Hermitian logarithm of a dilation is generally not local |
| C10 relativistic QFT/AQFT | Local algebras, spacelike commutativity, covariance, spectrum/positivity, admissible states; consistency of field content | Fresh spacelike signalling excluded. Neither a finite CPTP table nor a lattice H proves a relativistic completion. No universal embedding theorem is claimed |
| C11 controlled EFT with a specified UV class | Cutoff/domain/error, locality/symmetries, unitary amplitudes; analyticity/growth/crossing when used | Conditional positivity excludes some Wilson coefficients (T10). EFT locality alone is weaker than the full UV hypotheses; massless gravity needs separate treatment |
| C12 equilibrium KMS/passivity class | Specified algebra/time evolution; KMS temperature or cyclic work-extraction conditions, with complete passivity when invoked | Restricts equilibrium states/operations. Does not determine every kinetic mode/projector. A KMS state is not an arbitrary DB condition on invasive observed records |
| C13 generally covariant/gravitational model | Specific action/constraints, gauge-invariant relational observables and causal structure; quantum/UV assumptions stated separately | Coordinate redundancy must not become a physical setting or privileged clock. No universal finite-state lift into gravity or quantum gravity is proved here |

C1⊂C0, and C1 has a diagonal operational embedding into C7 with sufficient
dimension. C8 produces C7 after specifying ancillary state and discarding
it. A finite local classical circuit embeds in a finite local quantum
circuit on orthogonal states. C4 and C5 are independent restrictions on
appropriate classical models. C6 is not contained in C5; C9 is not C10;
C10 does not imply an exact finite-dimensional C7 representation; C11 is
a controlled domain description, not an exact superset of every model.
C2, capacity, geometry, energy, clock and thermality are intersecting axes.
There is no honest single "more physical" nesting containing all rows.

## B. What happens to the actual PR #6 controls

| Control | Exact conventional realization or obstruction | Remaining supplied resources / boundary |
|---|---|---|
| Eight-state adaptive kernel | Fixed C1 table and fixed C7 instruments; all feedback records preserved. Reversible XOR computes its deterministic part (T1,T3) | The prescribed rule, preparation, labelled readout, clock, setting selection, random bits, scratch/ancilla and any erasure. No identity or law selection emerges |
| Two-bit shift register | Fixed four-state classical register or four orthogonal quantum states reproduces I4. d<4 quantum storage cannot perfectly reproduce the common probe (T2) | Preparation dependence in external side channels invalidates a capacity claim. Four-state larger fixed comparator restores equality |
| Pólya α=β=1 | Fixed count-state urn reproduces both passive records and the same physical forced-success/count-increment operation. Static latent-p comparator does not | The count state is unbounded across horizons; fixed finite-horizon implementation needs count registers up to that horizon. Reinforcement/update readout is supplied, not an unpriced universal field |
| Power-of-two deterministic record clock | Fixed unbounded integer counter produces the record; no finite autonomous deterministic-record HMM does (PR #6) | Clock/history memory and elapsed step information are resources. External timestamps cannot be silently handed to a purported autonomous finite comparator |
| Smolin-style integer matrix recurrence | Fixed update on a pair of preceding matrices | Entries may grow without bound; finite precision/energy/capacity is not certified by formal computability. No Lorentz, KMS or conserved-Hamiltonian completion established |

T1 does not give a free finite-temperature realization of every one of these
controls. T3 gives finite logic-level reversibility; T5 and T6 exhibit actual
subclass obstructions. When a standard comparator outside the rejected
subclass survives, the failed subclass is not evidence for fundamental
law evolution. Across the tested controls there is no observation outside
all conventional physical model classes.

## C. Restrictions already supplied by established physics

Every row states premises, excluded possibilities, invariant content and
operational consequences. Additional premises are not generated by GRUT.
T1–T10 are proved in [PROOFS.md](PROOFS.md); other rows are scoped source
audits, not newly proved universal theorems. All source summaries are
paraphrases; no paper's prose is reproduced.

| Known restriction | Premises and excluded possibilities | Representation / measurable content and limit |
|---|---|---|
| Complete positivity [P01] | An operation extends consistently to arbitrary independent entangled ancillas. Positive transpose fails | Choi spectrum/ancilla output positivity invariant under unitary basis changes; Bell partial-transpose witness −1/2. Initially correlated restricted domains need separate treatment (T7) |
| Quantum vs no-signalling [P05] | Separated quantum tensor factors, compatible setting independence and no communication | CHSH≤2sqrt2 excludes PR=4 in every quantum dimension. No-signalling alone is insufficient; ordinary quantum physics already narrows the set (T7) |
| Finite memory dimension [P01,P13] | Same probe, retained resource boundary independent of preparation, no preparation side channel | Equiprobable decode success≤d/M; perfect I4 requires d≥4 (T2). Relabelling/coherent basis changes preserve it; changing resource boundary does not |
| CTMC embedding [P04] | Finite positive bounded-rate generator in physical time | detP>0; no-jump flip error≥exp(−ΛT). Quantum/reversible deterministic gates can evade by leaving this class (T6) |
| Causal locality / LR [P06] | Local algebra/gate geometry or bounded interaction norm/range, fresh independent setting | Strict circuit/CA cones; bounded Hamiltonian commutators with tails; AQFT spacelike operations cannot signal. Arbitrary graph relabelling does not preserve a frozen physical distance (T8) |
| Unitarity / reversible information conservation [P01] | Closed unitary including environments and controlled inputs | Orthogonal input histories must stay distinguishable; exact d-level N-reset environment dim≥d^N under T3. Marginal irreversibility needs retained/discarded information, not impossible physics |
| Thermal erasure / finite bath [P02] | Initial independent Gibbs bath, global unitary, heat and all resources counted | βQ=ΔS+I+D; finite faithful bath cannot exact pure-reset full-rank memory. Initial side information/nonthermal ancilla changes premises, not a universal counterexample (T3,T4) |
| DB / Onsager / KMS / passivity [P03,P07,P12,P26] | Respective stationary, microscopic time-reversal, equilibrium analyticity or cyclic-work hypotheses | DB pair current zero; Onsager–Casimir coefficients transform with odd variables/field reversal; KMS links positive/negative frequency spectra; complete passivity is stronger than one-copy passivity. These are not interchangeable. Rates and measured equilibrium dynamics may remain free (T5,T9) |
| Fluctuation relation [P20] | Forward/reverse protocols, microscopic reversibility, conjugate work convention and equilibrium preparations | Crooks p_F(W)/p_R(−W)=exp[β(W−ΔF)], where path support permits; it excludes arbitrary work distributions under those premises. No prediction for a supplied bath becomes a new GRUT law |
| Noether / gauge consistency [P09,P21] | Specified symmetry, variational dynamics/boundaries; chiral gauge theory field content and quantum Ward identities | Conserved charges and gauge constraints restrict dynamics; anomalous chiral U(1) content needs cancellation. Σq=Σq³=0 excludes (1,1) but both (1,−1) and (2,−2) pass: no charge-scale selector. Added anomaly-cancelling sectors change the premise |
| Topological quantization [P22] | Appropriate compact bundle or gapped band/projector over parameter space | Integral Chern class gives quantized transport in the applicable setting. Gauge frame changes preserve the integer; gap closure or changed topology can change it. It does not fix all Hamiltonian parameters |
| Reflection positivity [P10] | Euclidean correlators plus the full applicable OS regularity/covariance/symmetry hypotheses | Negative reflected Gram form obstructs a positive Hilbert-space reconstruction. RP alone is not the full existence theorem; this audit does not independently reprove the OS construction |
| EFT analyticity/unitarity/crossing [P08] | Mass gap, dispersion/growth, optical theorem, pole-subtracted forward amplitude | Nonnegative a2 in T10 excludes negative-sign coefficients in that UV class. On-shell data are invariant under admissible field redefinitions; massless gravity forward poles block naïve transfer |
| Bootstrap islands [P11] | Specified symmetry, operator gaps, crossing and positivity, numerical truncation/certification | Bounds constrain spectra/OPE data in that stated class. A numerical island with assumed gaps is not a proof of unique microscopic reality; no new bootstrap was run |
| Quantum error correction / holography [P14] | Supplied code subspace/operator algebra and a specific duality/background regime | Reconstruction/error-correction constrains accessible operator relations in that model. Code choice and boundary theory remain inputs; no universal metric or field-content generator follows from the audit |

Unitary basis changes, physical relabellings and gauge transformations leave
appropriate observable restrictions invariant. They do not authorize
changing energy, spacetime separation, resource counts, intervention
semantics or a readout's time-reversal parity. Interface covariance is an
operational obligation, not a way to discard the excluded observable.

## E. Consequences for GRUT and the minimum scientific comparator

R1 remains conditional on an independently certified common carrier.
Π, readout class and T cannot be freely changed between interventions;
quotient separation is not a measure of all physical influence. Supplying
a CPTP detector or declaring a resource bound does not derive that carrier.

The kinetic g/eigenprojector results demonstrate coarse operational
underdetermination, not missing microscopic law. Positive reversible
conductances still span E=N(N−1)/2 parameters for fixed π. The GKLS control
keeps Gibbs and thermal ratios fixed while ordinary relaxation/dephasing
couplings vary. If complete labelled short-time kernels or full microscopic
dynamics are supplied, these freedoms need not remain unknown. No additional
kinetic theorem family is added here.

TC1-style scaling comparisons need a frozen physical size, locality,
resource/energy bound, elapsed-time normalization and measurement error.
A model can match a scaling curve by supplied microscopic parameters or
resource growth. This audit neither certifies TC1 nor sets its coefficients.

Joint identity/response generativity remains an ambition. A new model must
not attach response rates to a selected slow-mode partition and call both
generated. The current evidence does not establish that an extra universal
joint generator is required by nature.

The **narrowest useful operational comparator** is a fixed finite-horizon
instrument family with the actual intervention/readout contract and an
independently certified total retained dimension, clock/latency, energy,
reset and geometry budget. In a quantum experiment use quantum instruments
or process tensors with all environment memory counted; add DB/KMS only
if the preparation/protocol earns it, and relativistic/QFT premises only
if relevant. This is a conditional definition, not a certified available
Stage-3 lock. All ordinarily adjustable instrument/bath parameters within
the declared independent calibrations must remain allowed.

For a future law to reduce physical freedom beyond standard consistency,
there must exist a conventional comparator satisfying **all** independently
earned standard conditions and the same operations, but violating a
preregistered relation. Capacity/heat/LR/positivity relations already forced
by those conditions cannot serve as new-law evidence. No such extra
relation has been obtained or candidate-tested here. This conclusion is
an evidence boundary, not a universal impossibility theorem about GRUT.
