# S5-P1: certify formation without resolving every field detail

The physical target is a geometry-event probability. The post-recording
S4 theorem already holds for every admissible conditional receiver
state. Consequently the formation calculation need not reproduce each
receiver observable accurately to establish the joined protocol.
It must instead bound the geometry-event probability under the original
Hamiltonian and seed, including every omitted field contribution.

## A second elementary bound

Let normalized states psi (physical) and phi (computed) obey
||psi-phi||<=e. Let G be the same macroscopic event projector frozen
in S5, and p_hat=||G phi||^2. S5 used

\[
 p_{\rm physical}\ge(\sqrt{p_{\rm hat}}-e)_+^2.
\]

Applying exactly the same triangle inequality to 1-G instead yields

\[
 \|(1-G)\psi\|\le\sqrt{1-p_{\rm hat}}+e,
 \qquad
 \boxed{p_{\rm physical}\ge
 \max\{0,1-(\sqrt{1-p_{\rm hat}}+e)^2\}.}
\]

Use the maximum of the two valid lower bounds and zero. For p_hat=1,
the complementary-event bound is 1-e^2, whereas the first bound is
(1-e)^2. Thus an e=.55 state-error account can still establish a
formation probability near .6975 when the computed event is nearly
certain. The error account itself is neither reduced nor discarded.

This theorem requires normalized states, the identical event, and a
valid complete state-error account. It does not infer an uncomputed
target-time probability from a large intermediate population. It
does not assume any field state is favorable or inaccessible.

## Separate numerical charters, one physical test

S5's original envelope (state error <=.20, Strang <=.04, projections
<=.15) and every resource stop remain in RESULTS.json. A P1 success
is not a success under that original envelope.

PROBABILITY_CHARTER.md, published before P1 trajectories, instead
allows state error <=.60, Strang <=.20 and projections <=.35. The rank,
coefficient and step-boundary compute ceilings are retained. The
Hamiltonian, original root weights, initial state, disorder families,
T_N=4N^2, geometry event and required probability >=1/2 are unchanged.
Only the sufficiency certificate and numerical allocation change.

Completed S5 states are reused. All unresolved primary components are
attempted under P1, without selecting them by their intermediate
formation population. Unknown components keep zero lower bound at
their original root weight. Continuous finite-N comparison boxes use
the same full-H norm comparison described in MATHEMATICS.md.

For those P1 boxes, the available per-root error allowance is

\[
 a(p_{\rm hat})=\max\{\sqrt{p_{\rm hat}}-1/\sqrt2,
                         1/\sqrt2-\sqrt{1-p_{\rm hat}}\}.
\]

Use the smallest positive margin \(s=\min_{\rm roots}(a-e)\), then
\(h=s/[4T(N-1)]\) and \(q=s/(4T)\). The added Duhamel error is
at most \(s/2\), so every original root component still has formation
probability at least one half. This is the complementary-event version
of the original S5 comparison, not the smaller reverse-triangle margin.

The projection, product-formula and finite-exponential errors still
cover every omitted return path. No irreversible bath, missing word,
unpaid virtual state, altered observation time or postselected root is
introduced. Floating-point evaluations remain non-outward-rounded.

## Scope of the joined conclusion

If the weighted finite-instance lower bound is p_*, then the exact
preparation-recording-probe protocol has formation probability >=p_*,
with S4's conditional response/protection bounds on every accepted
record and its actual field. The nearby signal per original trial at
t=.005 is at least .00102794836247... times p_*.

Any such result is still restricted to the solved finite families
and their explicit finite-N comparison neighborhoods. It establishes
neither a fixed-width all-pattern large-N formation theorem nor a
spacetime phase. This numerical refinement does not restore the lost
real-energy consistency gap or supply a distinctive prediction.
