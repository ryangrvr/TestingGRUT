# WM2-S8 — coherent Schmidt receiver and charge-rank allocation

**Interim checkpoint:** original S8 recovery failed its numerical limits;
the separately frozen S8-G1 allocation variant is under evaluation. No
six-event formation result is asserted by this checkpoint.

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

Its source is checked to differ from original S8 only at this allocation
rule/default and its import. G1 reruns the full controls, including the
576-mode physical reference. It also checks an exact rank-capacity witness
in RANK_ALLOCATION.md: 18 or 33 orthogonal histories can require that many
joint Schmidt vectors in one charge while each conditional history needs
one. Those packets are allowed states, not demonstrated neutral-seed
outputs or a no-go for a physical phase.

G1 target-time outcomes, fresh S7 whole-worker costs and any gated six-event
result will be evaluated with their actual completed/stopped status.
No additional numerical variant is permitted by this packet's charter.

## Scope

This campaign tests whether receiver information can be retained at useful
cost. It does not derive dimension, original CP-12 matter universality,
a real-energy consistency gap, gravity or a distinctive empirical prediction.
Even a successful finite recovery or six-event result would leave general
independent-disorder formation and all-size resource scaling unresolved.
