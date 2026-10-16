# WM2-S6: exact cross-history receiver calculus; formation frontier unchanged

S6 derives and tests the missing transition and error calculus for
history-conditioned receiver frames. The prototype is correct in the
tested complete models and uses substantially fewer field coefficients.
It does not add a controlled primary formation result beyond S5. The
rank-ten heuristic encounters residual, projection/rank and compute
limits. Six-event formation at T=144 remains unresolved.

The physical S4 law, original neutral preparation, physical receiver,
root weights and observation protocol are unchanged. These solver
limits are not a physical no-formation theorem. Indeed, S5 already
established finite formation at the same four- and five-event centers.
All-size formation under independent storage disorder remains OPEN.

## What remains exact and what is compressed

Each ordered detuning-type word has its own orthonormal receiver frame.
All allowed words, reverse transitions and charge-allowed occupations
are retained. Every frame column has all L_N=64N^2 physical frequency
coordinates. Only numerical receiver subspaces are compressed; the
model has no substituted small bath, additional photon cutoff, forgotten
attachment order, selected root or irreversible emission law.

Different frames cannot be treated as equal. The complete append block
from parent u with k photons to child v is

\[
 -\lambda_v a^\dagger(W_v^\dagger c)
       \Gamma_k(W_v^\dagger W_u),
 \qquad \lambda_v={g\over N}\sqrt{n_a-k_a(u)}.
\]

Gamma is the full normalized rectangular bosonic Fock contraction.
Removal is its exact adjoint. Dense and matrix-free SVD/Givens
implementations retain every phase and bosonic multiplicity.

Although source histories are orthogonal, their Hamiltonian images can
enter the same target history. S6 computes every sibling and
parent/grandchild interference term in the full action norm. Direct
physical-matrix controls change the squared action norm by up to about
.0132 when those cross terms are removed. They cannot be silently
discarded in a history-conditioned solver.

## A physical-H error bound, not a guessed residual

With isometry V, P=VV^dagger, C=V^dagger B V and R=(1-P)BV, Duhamel gives

\[
 \|e^{-iBh}Vx-Ve^{-iCh}x\|
 \le\int_0^h\|Re^{-iCs}x\|\,ds.
\]

The prototype evaluates the first three full residual coefficients and
bounds the entire remaining tail using ||B||<=g. It adds the original
physical-H Strang error, all receiver projections, finite-exponential
tails, frame changes and global normalization costs. Physical unitarity
pays every omitted component's later return and interference.

No history is separately normalized. A discarded orbital in a rare
history has its original global weight, and the full state is normalized
only globally. Unknown target-time components keep zero trivial lower
bound at their original root weights, not a renormalized success sample.

## Two openly separated numerical protocols

Original S6 was frozen at e627bb5f8287efb47d2ab4c79efbad6fd2984132 before
its controls/trajectories. Its state-error envelope is .60: physical-H
Strang <=.15, projections <=.20, complete projected-emission residual
<=.20, other contributions <=.05. It assigns a local share of the
remaining residual budget to each future step.

That local quota stops the smaller disordered runs while most of the
cumulative residual allowance is unused. S6-G1 was separately frozen
at 23de49378a64c16c405db87c7b1a703952db9813 before its reruns. It uses
the same .20 cumulative residual bound directly, without the extra
per-step quota. Every original stop remains in RESULTS.json. A G1
outcome is not an original-quota outcome or a smaller full-state error.

Both protocols retain rank ten per history, 2,000,000 coefficients,
8,000,000 cached dense overlap entries, degree-eight exponentials,
three residual coefficients and a 180-second step-boundary compute
limit. Both test balanced +/-1/4 storage offsets at N=4,5,6, common
detuning 0 or 1/2, g=.9, and the original T_N=4N^2. Every unique
primary root component is attempted with the original root weights.

## The physical-receiver controls work; the primary solver remains limited

The complete N=3 reference uses all 576 physical modes, dimension
333,505, and the full T=36. The history-conditioned approximation's
full-wavefunction error is .0281479931, below the evaluated analytic
account .3079119582. The comparison uses the entire two-photon tensor.
It is not merely a population, norm or small substituted-receiver test.

The physical N=4 uniform control reaches T=64, with full-state error
about .00010648 against the known exact quadratic solution. It uses
only ten coefficients at peak. This is a positive implementation
control, not independent-disorder formation evidence.

For the primary disordered families the original policy reaches local
quota limits at N=4,5 and compute limits at N=6. G1 advances the smaller
cases further, then reaches rank/projection or cumulative residual
limits. Its six-event target remains unresolved under the compute cap.
The original and G1 ledgers give all root-component endpoints and error
costs. No intermediate population is promoted to an observation at T_N.

The endpoint ranges below cover all attempted root components and both
common detunings at each size; they are not uncertainty intervals.

| N | Required T_N | Original S6 completed-time range | G1 completed-time range |
|---:|---:|---:|---:|
| 4 | 64 | 15.42-17.02 | 29.80-31.97 |
| 5 | 100 | 13.93-14.84 | 21.09-25.90 |
| 6 | 144 | 6.34-6.45 | 14.44-14.61 |

Every primary family remains UNRESOLVED at its required observation.
The zero trivial family lower bounds used for unknown roots are not
claims that the true formation probability is zero.

No new primary preparation-recording-probe lower bound is claimed from
these stopped trajectories. S5's evaluated finite formation bounds and
their S4-based joined protocol remain unchanged and remain the strongest
solved disordered instances in this construction family.

## Why fewer coefficients did not settle scaling

For balanced non-root counts (1,2), (2,2), (2,3), ranks at most ten
give exact coefficient ceilings 846, 5,851 and 28,951. A shared rank-23
six-event frame would instead use 974,051 coefficients. These are
storage counts, not an equal-error speed comparison.

Separate history frames introduce complete cross-Fock contractions.
Residuals also require sibling and grandparent/grandchild contractions.
At six events the matrix-free contractions are computationally costly.
The full method additionally stores physical frame coordinates,
occupation rotations and overlap plans; coefficient count is not total
RAM. COST_AND_SCOPE.md derives these costs explicitly.

For all distinct non-root detunings the history count remains
sum_k (N-1)!/(N-1-k)!, asymptotic to e(N-1)!. The physical-H splitting
needs O(N^3) steps at T_N=4N^2 for fixed error. No uniformly sufficient
rank, polynomial runtime or all-pattern formation theorem is obtained.

## A further exact construction, not substituted into the runs

The recorded prototype chooses borrowed directions from separate
neighbor densities. ACTION_ENRICHMENT.md instead derives the one-body
density of the complete coherent Hamiltonian action, including the
cross-history phases. It can be evaluated through smaller, one-fewer-
photon Fock contractions and a low-rank factorization, without an
L_N by L_N physical density matrix. One added direction captures the
whole missing current action in a one-photon target.

More generally, if enlarged history frames contain B^l Phi for l=0..q,
the first q complete residual coefficients vanish. A sufficient explicit
construction unions receiver frames within q primitive word moves and
contact. For q=3 at T=4N^2 and O(N^3) Strang steps, the summed remaining
residual tail is O(1/N). This is a conditional error-scaling theorem.
Its enclosure ranks and Fock costs may still grow prohibitively.

The action-density, cancellation, one-photon capture and Krylov-support
identities are checked against complete physical matrices. This module
was derived after the heuristic tests began and was not used to alter
either primary protocol. No new formation number, efficient rank bound
or original-gate pass follows from its algebraic controls.

## Verification, execution and boundaries

The archive contains unequal-rank complex Fock checks, independent
physical append/pop and residual comparisons, complete interference
controls, six small propagated wavefunctions, the full 576-mode reference,
the physical uniform control, larger matrix-free contractions and the
exact common-frame S5 limit.

The final matrix-free evaluator precomputes fixed two-mode rotation
matrices; original S6 computed those same matrices on demand. Cached
and uncached arithmetic is bit-identical in the saved controls. G1 uses
the cached evaluator. Wall-time comparisons therefore cannot isolate
only the residual-allocation policy. Both policies' original arithmetic
results and stops remain preserved.

Arithmetic is ordinary double precision. In particular residual norm
subtraction can suffer cancellation; tiny negative squared residuals are
recorded and clipped, rather than enclosed by interval arithmetic.
NUMERICAL_LIMITS.md makes this debt explicit. Analytic reduction bounds
and complete-matrix consistency do not constitute an outward-rounded
proof or independent mathematical validation.

This is a tested numerical/theoretical tool for the unchanged candidate
model. It does not restore a real-energy gap, establish size-uniform
formation, select dimension, make original charged CP-12 matter share
the geometry, or predict a distinctive observation about nature.
