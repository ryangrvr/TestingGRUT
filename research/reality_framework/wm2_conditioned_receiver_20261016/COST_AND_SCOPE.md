# Computational cost and physical scope

S6 changes a numerical representation of the unchanged S4 formation
law. A receiver frame can depend on the complete type word. It contains
all L_N physical frequency coordinates, and its Fock coefficients have
exactly the photon number fixed by that history. Original labels are
reduced only by S4's exact repeated-detuning class isometry.

## Exact finite storage counts

For frames of ranks r_u and k=|u|, coefficient count is
sum_u binomial(r_u+k-1,k). Physical frequency-coordinate storage is
L_N sum_u r_u. The vacuum word's frame is computationally redundant
and is retained by the prototype for uniform bookkeeping.

At rank ten the balanced N=4,5,6 root components with non-root counts
(1,2), (2,2), (2,3) have coefficient ceilings 846, 5,851 and 28,951.
These are exact combinatorial counts, not a proof that ten orbitals
suffice at the frozen error or time. S5's common rank 23 at the last
component requires 974,051 coefficients. Comparing these numbers shows
a possible reduction in coefficient storage, not a runtime speedup or
an equal-error comparison. S5 and S6 use different numerical allocations.

The complete prototype also stores per-history physical frames, dense
cross-Fock maps below a fixed threshold, matrix-free SVD/Givens plans,
and reusable occupation/rotation indices. Coefficient count alone is
not total RAM. The 8,000,000-entry cap applies to cached dense Fock
overlap matrices. Rotation indices and SVD plans are counted by the
formulas below, rather than represented as that dense-cache metric.

## Work introduced by distinct history frames

An emission/removal application needs one parent/child overlap per edge.
The full residual additionally needs siblings and grandparent/grandchild
overlaps, because amplitudes from those sources can enter the same target
history. Dropping that work drops physical interference.

Constructing one-particle overlaps costs O(L_N r_u r_v). A dense
Gamma_k(G) stores D(r_v,k)D(r_u,k) entries and can be built recursively
using occupation creation maps. This cost can exceed the coefficient
cost. Matrix-free evaluation factors G=USV^dagger and lifts the two
unitaries into normalized Fock space. Each unitary uses at most O(r^2)
two-mode rotations. Each rotation groups D(r,k) coefficients into blocks
of length at most k+1. Application work is at most

    O((k+1)[r_u^2 D(r_u,k)+r_v^2 D(r_v,k)]).

Reusable pair-rotation indices cost O(r^2 D(r,k)) per occupied (r,k)
space. SVD/Givens plan storage includes O(r_u^2+r_v^2) one-particle
parameters, phases and occupation-diagonal data. These plans are exact
algebraic evaluations of contractions; no singular value is silently
discarded. Cost is paid again if an adaptive frame changes.

S5 can apply its sparse emission operator in a shared frame. S6 trades
that sparse shared representation for many cross-frame contractions.
Thus fewer coefficients can coexist with slower propagation. Measure
target-time completion rather than assuming compression is a speedup.

## Size dependence that remains

With all N-1 non-root detunings distinct, the complete history count is

    sum_(k=0)^(N-1) (N-1)!/(N-1-k)!.

It is asymptotic to e(N-1)!. Binary repeated-detuning families still have
exponentially many complete order histories. Neither can be replaced by
an attached subset: different orders have different physical histories.

The physical-H Strang step count for a fixed error scales as T^(3/2);
at T=4N^2 this is O(N^3). Every step includes all represented histories,
their transitions and three exact full-action residual coefficients.
The cumulative residual and mode-projection errors depend on the actual
state. No size-uniform bound on their ranks, computational work or
formation probability is proved here.

## Error and scientific interpretation

Within exact arithmetic, the residual theorem includes every omitted
future return and all coherent target-history interference. A physical
formation statement additionally needs completion at T, the original
root weights and an informative event lower bound. Numerical policy
failure, rank exhaustion and a wall-time stop are solver failures, not
theorems of physical non-formation.

S6's local residual quota is stricter than the cumulative theorem.
S6-G1 is separately frozen and pays the same .20 global residual budget
without that local quota. Preserve both outcomes. G1 cannot be relabeled
an original-policy success, a smaller full-state error, or a new law.

Ordinary double precision, residual norm cancellation and finite-frame
Gram defects are disclosed. The analytic formulas are not a complete
roundoff enclosure. Full-matrix comparisons are by the same Builder;
external independent validation is not claimed.

Any finite formation result joins only S4's actual conditional theorem:
record every accepted word, retain its actual receiver, and pay the
acceptance probability and recording/probe work. This numerical method
does not select dimension, make original charged CP-12 matter obey the
path, restore a real-energy gap, or derive a distinctive prediction.
