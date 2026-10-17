# Separate analytical construction: orthogonality within one-way sweeps

This analysis was developed during the frozen G1 campaign. **It is not used
by either S8 or G1 trajectories, truncation choices, ledgers or gates.** No
extra primary variant is launched here. It identifies a tighter projection
account to implement and verify in a future campaign.

The orthogonality principle is standard in tensor-train SVD error analysis
(Oseledets, 2011, Theorem 2.2). The lemma below states the conditions needed
when local unitary gates are interleaved with the projections. It is a
numerical-analysis construction, not a new physical law or GRUT prediction.

## Exact nested-projector lemma

Let U_1,...,U_m be unitaries in a one-way sweep. At cut j choose an
orthogonal kept-prefix projector P_j, with P_0=I. Suppose

* P_j <= P_(j-1), after embedding all prefix spaces in the full Hilbert space;
* every later U_k, k>j, commutes with P_j;
* the stored state x_j belongs to range(P_j).

Define the orthogonal projection loss and the in-range implementation error

\[
 e_j=(I-P_j)U_jx_{j-1},\qquad
 f_j=x_j-P_jU_jx_{j-1}.
\]

Then

\[
\boxed{\|U_m\cdots U_1x_0-x_m\|
 \le\sqrt{\sum_j\|e_j\|^2}+\sum_j\|f_j\|.}
\]

Proof: telescope the local differences e_j-f_j through the remaining
unitaries. For j<k, the propagated e_j remains orthogonal to range(P_j),
because all later gates commute with P_j. The propagated e_k stays in
range(P_(k-1)) subset range(P_j). Thus all propagated e_j are pairwise
orthogonal. Apply Pythagoras to their sum and the triangle inequality to
the propagated f_j. No assumption about irreversible amplitude loss enters.

For an exact left-to-right tensor sweep, a new kept prefix is embedded
through the preceding kept prefix, so its range is nested. Later gates
act only to its right. The right-to-left version follows by exchanging
prefix and suffix. The lemma applies to one such monotone sweep; it does
not combine arbitrary forward/backward sweeps into one orthogonal batch.

## Why one cannot square-sum losses over all time

Take x_0=(1,0), U a two-dimensional real rotation by theta and P=|0><0|.
Apply U then P twice. The exact state is U^2 x_0; the stored state is
cos(theta)^2 x_0. For theta=.1, the actual norm error is about .1989,
whereas the root-sum-square of the two local discarded norms is about
.1408. The subsequent rotation revisits the discarded coordinate and does
not commute with P. Returning errors can add coherently.

Therefore a future S8-style implementation may combine losses *within each
certified monotone sweep*, but must add the resulting sweep errors across
reversals, boundary operations and physical steps. The full physical-H
reference and its coherent returns are retained.

## The numerical obligation still missing from S8/G1

Floating SVD factors are not exact projectors. For a future implementation:

1. Define P_j from the exact span of the stored numerical prefix map.
   These ranges are structurally nested, even if the columns are not
   perfectly orthonormal.
2. Enclose the actual orthogonal loss e_j and the in-range defect f_j using
   checked prefix Gram matrices, SVD reconstruction, gate errors and all
   arithmetic. A small raw singular tail alone does not enclose f_j.
3. Accumulate the positive squared e_j bounds inside a monotone sweep;
   add all f_j bounds linearly. Reset the squared batch on every reversal.
4. Freeze this implementation and its cost/error gates before repeating
   recovery and any six-event target.

A concrete way to enclose f_j avoids constructing a physical prefix vector.
Let L be the stored prefix map before the gate, H=L†L its Gram matrix,
K the kept local left-factor matrix, and X the untruncated two-site tensor
in those coordinates. Set H_new=K†(H tensor I)K. The exact orthogonal
projection coefficient is

\[
Y=H_{\rm new}^{-1}K^\dagger(H\otimes I)X.
\]

If Z is the stored retained coefficient factor, its in-range defect can
be bounded positively through

\[
w=H_{\rm new}Z-K^\dagger(H\otimes I)X,\qquad
\|(L\otimes I) K(Z-Y)\|\le
 {\sqrt{\|H_{\rm new}\|}\over\lambda_{\min}(H_{\rm new})}\|w\|.
\]

with the right-environment factor and all construction/gate errors added.
Checked near-isometry bounds control the inverse; no inaccessible global
spectral gap is assumed. Computing a direct residual w plus its arithmetic
enclosure is different from subtracting two large squared norms. This is
a specification for a future implementation, not an implemented S8/G1 fee.

The current ledgers conservatively sum all local tail bounds linearly.
They do not compute the required orthogonal/in-range decomposition, so this
lemma cannot retroactively lower those ledgers or complete a stopped run.

`sweep_control.py` checks a full 16-dimensional monotone-gate construction,
including explicit nested projectors and transported error inner products,
an in-range error example, and the coherent-return counterexample. Its
arithmetic is ordinary double precision; the proof is the lemma above.
