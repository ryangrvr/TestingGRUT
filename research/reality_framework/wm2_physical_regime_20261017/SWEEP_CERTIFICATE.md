# Production certificate: one-way projection errors

This is the pre-implementation specification. Implementation and executed
control status will be recorded separately; no success is presumed here.

## Geometry of the numerical ranges

Inside one left-to-right sweep, each stored prefix is the previous stored
prefix times a kept local factor. Exact spans of these floating matrices
are nested. Every subsequent ideal gate is to their right and commutes
with the earlier prefix projector. Reverse the roles for right-to-left.
The history site remains a coherent register, not a selected history.

Let y_j be the physical tensor represented by the computed pre-truncation
center, P_j the exact orthogonal projector onto its new stored prefix,
x_j the retained tensor. Split U_j x_(j-1)-x_j into

    e_j=(I-P_j)y_j, f_j=P_j y_j-x_j,
    a_j=U_j x_(j-1)-y_j.

The e_j remain in the preceding nested range and are orthogonal after
all later ideal gates in the batch. Therefore one batch costs

    sqrt(sum_j ||e_j||^2) + sum_j (||f_j||+||a_j||).

Add batch costs linearly across every reversal and return. Skipped exact
vacuum gates are identity on the stored state and do not widen a batch.

## Direct in-range defect, without physical prefix vectors

In each charge block let X be the computed two-site center, K the kept
factor, Z the actual stored coefficient factor, and R=X-KZ. For a left
sweep use the adjoint-oriented X,K,Z. Let H=L^dagger L be the unknown
physical prefix Gram matrix. S8's checked local Gram bounds imply
||H-I||<=h=s/(1-s), s=sum eta_j<1e-6. Let eta enclose ||K^dagger K-I||.

Compute positively, with multiplication/subtraction fees,

    R=X-KZ, w_0=K^dagger R.

Then

    ||w|| <= ||w_0|| + h sqrt(1+eta)||R||,
    w=K^dagger(H tensor I)R.

The true orthogonal projection coefficient Y satisfies
H_new(Y-Z)=w, H_new=K^dagger(H tensor I)K. Consequently

    ||(L tensor I)K(Y-Z)||
      <= sqrt((1+h)(1+eta))/[(1-h)(1-eta)] * ||w||.

Multiply by the checked right-environment factor. This is a direct
residual account with a rigorous unknown-Gram correction, not subtraction
of two squared norms or reliance on SVD convergence. Charge-block defects
combine by their positive Euclidean norm. The orthogonal loss is at most
the physical norm of R. Gate/input arithmetic is linear and distinct.
The Gram inverse is numerical conditioning, not a physical energy gap.

Exact spans of stored factors define the projectors. SVD reconstruction,
residual evaluation, new-core scaling, isometry and product rounding still
require their declared positive fees. Production code must check these
formulas against dense nonorthogonal-prefix controls before any primary.
