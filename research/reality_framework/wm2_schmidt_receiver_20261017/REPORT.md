# WM2-S8 — coherent Schmidt receiver and charge-rank allocation

**Completed campaign:** both original S8 recovery components stop before
the target. The separately frozen S8-G1 global allocation completes one
component, while the other stops on the Schmidt allowance. Neither method
passes the two-component recovery gate; the six-event target is not run.
These are numerical-method outcomes, not evidence against physical formation.

S8 implements a charge-conserving MPS for the complete history register
and receiver chain. It retains emitted photons and their correlations,
instead of projecting a receiver orbital onto vacuum. Every discarded
Schmidt component is charged positively against full physical-state error,
with explicit gate, SVD reconstruction, canonical-environment, arithmetic,
field-factorization and receiver-boundary accounts. This is a numerical
representation of the unchanged S4 law, not a new microscopic Hamiltonian.

## Frozen protocols and unchanged physical target

Original charter `0e42a5d` predates its controls and trajectories. Only
this folder changes; S4-S7, CP-12 and SCOUT-1 stay frozen. All ordered
detuning-type histories, reverse transitions, original root weights and
exact photon charge remain. Boundary gates act coherently on the complete
history register. A history is never normalized individually.

Physical L_N=64N^2 and the original neutral/receiver-vacuum preparation
remain. The computational prefix R=min(L,ceil(12T/5)+32) is a **paid
approximation** to that physical receiver, embedded in its remaining
vacuum modes. A full-H weighted-charge/Duhamel argument, including every
possible return, bounds the discrepancy by sqrt(N-1)/2^32. It is not an
unexplained small-bath substitution. The receiver's supplied chain is
not claimed as dynamically selected physical-space geometry.

The N=5 target remains T=100, physical L=1600, g=.9, common detuning
zero, offsets +/-1/4, original root weights (3/5,2/5), and m>=3. Full
assembly m=5 is a separate event. A method must complete both components
and obtain a paid root-weighted lower bound >=1/2 before its six-event
test opens. Unknown endpoints do not count as T=100 observations.

## Original S8 stops

| Remaining root counts | Completed whole-step time / required 100 | Wall seconds | Peak MiB | Stop |
|---|---:|---:|---:|---|
| (2,2) | 30.0654 | 53.022 | 56.305 | Schmidt cap / truncation budget |
| (3,1) | 57.3144 | 185.399 | 58.746 | Schmidt cap / truncation budget |

The 16-per-charge cap bound while total retained bond rank was only 55
and 53 respectively, below the total ceiling 96. Both consumed almost
all .35 of allowed Schmidt truncation. Neither exhausted the 900-second
worker or 2048-MiB memory budget. Both stopped within a sweep, so their
JSON explicitly distinguishes completed whole-step time from a complete
stored endpoint state. No target-time event interval is provided.

These are failures of the frozen *numerical recovery protocol*, not
negative physical formation observations. The stopped runs are preserved
and are never promoted to passes by a later rank-allocation rule.

## Independent physical and arithmetic controls

Checks include exact rational positive norms; checked Gram enclosures;
complete field-operator factorization comparisons; exact first-three-mode
commutator algebra with every spectator-photon block; and six full small
physical wavefunctions at T=.5,2,8. Their errors lie inside the accounts.

The physical reference compares **all 333,505 amplitudes** of the original
N=3, L=576, T=36 matrix model. The tensor's 119-site computational prefix
is embedded into the full physical vector, including zero approximate
amplitudes outside the prefix. Its wavefunction error is **0.00234663**,
below the evaluated account **0.14468294**. This tests the full state,
not a selected population or a cancellation-prone missing-norm subtraction.

The arithmetic claim is conditional on the declared IEEE/BLAS model in
ARITHMETIC.md. SVDs and eigensolvers are checked by direct residuals and
Gram defects; library output is not treated as exact. Controls do not
formally verify the software stack or independently validate the physics.

## Separately frozen S8-G1

G1 charter `b15fcf0` was published before G1 controls or trajectories.
It changes one allocation rule: choose the globally largest singular
values across exact charge blocks under the **same total ceiling 96**.
The per-charge ceiling becomes 96 instead of 16. Error, memory, time,
operator audits, physical law, targets and root weights stay fixed.
This is a different numerical protocol, not a pass under original S8.

At N=5 the original five charge blocks also imply an effective total
maximum 5*16=80, whereas G1 can use 96. Thus G1 changes both the admissible
distribution and usable total rank at this size; it is not an equal-rank
comparison. PROTOCOL_CLARIFICATIONS.md records this without changing a gate.

Its source is checked to differ from original S8 only at this allocation
rule/default and its import. G1 reruns the full controls, including the
576-mode physical reference. It also checks an exact rank-capacity witness
in RANK_ALLOCATION.md: 18 or 33 orthogonal histories can require that many
joint Schmidt vectors in one charge while each conditional history needs
one. Those packets are allowed states, not demonstrated neutral-seed
outputs or a no-go for a physical phase.

The first G1 component stops at completed-step time **88.34422658** rather
than 100, after **839.938 seconds**, with **90.102 MiB** peak RSS. It reaches
global rank 96, a charge rank of 50, and 996,160 retained tensor entries.
Its Schmidt ledger is **0.3499998723** of .35; the next truncation cannot
fit. Runtime and memory caps do not cause this stop. Its accumulated ledger
totals about .464346, but this mid-sweep number is not a certificate for
physical evolution at the last whole-step time. No target-time interval is
reported. Failure of this component keeps G1's six-event gate closed,
regardless of its other component.

The second G1 component, remaining counts (3,1), completes T=100 in
**778.287 seconds**, at **72.445 MiB** peak RSS. Its maximum bond rank is
88, maximum charge rank 35, and maximum retained tensor count 503,999.
Its state-error account is .33527211, including Schmidt .20583146,
Strang .11219854 and field factorization .01701847. The paid probability
lower bounds, rounded downward, are **88.7514%** for m>=3 and **82.9695%**
for m=5 under the declared arithmetic model.

This is a completed root component of weight 2/5, not the original
root-weighted family. The unresolved component has weight 3/5. Its missing
target-time result prevents G1 recovery and keeps the six-event gate closed.
No favorable roots from different numerical methods are combined to pass
either charter. No extra numerical allocation/rank/budget variant is run.

RESULT_TABLES.md and the raw JSON give every component's status, costs,
error accounts and probability intervals. The fresh immutable S7 workers
provide the completed reference family. The methods differ in integration,
rank constraints, error allocation and time ceilings: this is not an
equal-error or equal-rank speed comparison. A stopped run is not a cheaper
solution of the target. G1 improves the reached time in its harder component
but does not solve the receiver-retention problem within its frozen limits.

| Recorded two-component protocol | Complete components | Sequential wall s | Maximum worker peak MiB | Weighted macro/full lower bounds |
|---|---:|---:|---:|---|
| Original S8 | 0/2 | 238.421 | 58.746 | No target-time family bound |
| S8-G1 | 1/2 | 1618.225 | 90.102 | No target-time family bound |
| Fresh S7 reference | 2/2 | 220.482 | 86.496 | 69.8713% / 57.9462% |

Displayed probability lower bounds are rounded downward. Wall sums include
each worker's imports, setup, operator audits and propagation, but not the
separate controls. Peak RSS is the maximum of sequential workers, not their
sum. Timers stop before final result serialization/printing. The reference
agrees in scope with S7's earlier result; small differences in evaluated
endpoints are retained rather than replacing the fresh observations with
archived favorable numbers.

## Separate analytical construction for projection accounting

SWEEP_BOUND.md states a conditional error lemma for a single one-way
sweep: when kept prefix ranges are nested and later gates commute with
earlier prefix projectors, transported orthogonal losses are pairwise
orthogonal. Their norm can be square-summed; numerical in-range defects
must still be added linearly. The document specifies a direct prefix-Gram
residual for that missing numerical defect.

Dense controls check the nested-projector identity and a nonorthogonal
prefix metric. A two-step rotation/projection counterexample gives actual
error about .198919 versus an invalid all-time square-sum .140834: returns
can add coherently. Consequently the tighter account applies only to an
individually certified monotone sweep, not across reversals or all steps.
This adapts a standard tensor-SVD orthogonality principle; it is not a
new quantum law, a production numerical certificate or a physical prediction.

**Neither S8 nor G1 uses this lemma.** No stopped ledger is lowered, no
extra primary family is run and no target completion is inferred from it.
The next constructive task is to bound the in-range projection defect
in production arithmetic before a new recovery-first campaign.

## Scope

This campaign tests whether receiver information can be retained at useful
cost. It does not derive dimension, original CP-12 matter universality,
a real-energy consistency gap, gravity or a distinctive empirical prediction.
Even a successful finite recovery or six-event result would leave general
independent-disorder formation and all-size resource scaling unresolved.
