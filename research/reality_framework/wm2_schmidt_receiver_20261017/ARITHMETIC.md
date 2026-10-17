# S8 arithmetic and physical-error contract

This is a roundoff-aware account under declared IEEE/reduction assumptions,
not formally verified software, a hardware certificate or external review.
Numerical checks test implementations; they do not establish these assumptions.

## Scalars and positive quantities

Use binary64 round-to-nearest with unit roundoff 2^-53 and correctly rounded
square root. The frozen S7 strict C++ kernel accumulates positive squared
components in `long double`; its runtime significand width is read, and
conversion to binary64 is charged. Positive reductions use gamma_n=nu/(1-nu)
plus an operation-count times the smallest normal binary64 underflow pad.
Scalar endpoints are rounded outward. Arithmetic controls compare these
enclosures against exact rational sums and independently computed Gram norms.

BLAS and LAPACK interfaces must return the declared complex128 buffers;
their matrix products are assumed to obey standard absolute reduction bounds,
including under cancellation. No accuracy claim about an SVD or eigensolver
is accepted solely from that library's convergence status.

## Local commutator certificate

Every allowed spectator-photon block is included. Hermitian matrix norms
use all computed eigenvalues, explicit `K*V-V*diag(lambda)` residuals and
outward Gram bounds for V. Matrix-construction error, diagonal symmetrization,
products and reductions are covered by
gamma(4096*d+4096)*(1+||G||_infinity+||F||_infinity)^3, d the block dimension,
and a separate eigenvector-product bound. These deliberately overcount
operations in the finite nested commutators. Exact primitive rational/square-
root coupling errors are dominated by this allowance. No favorable state or
low-energy spectral weight enters the certificate or timestep choice.

## Complete gates and field factorization

Boundary G gates use a degree-32 Hermitian Taylor polynomial for each complete
local charge block. Its integral remainder is (||G||*|s|)^33/33!. The matrix
polynomial rounding account propagates positive infinity-norm term bounds,
with gamma(512*d+128) per multiplication and separately paid sum reductions.
Primitive coupling errors are covered by gamma(1024)*(1+||G||_infinity).
Both append and its exact reverse are retained in the real symmetric input.

Receiver gates use the frozen S5 normalized finite Fock polynomials. S7's
fast absolute polynomial bound pays every coefficient before cancellation.
Independent interval sine/cosine enclosures bound the one-particle gate
phase error. Thus the local reference gates are true unitaries, even though
their floating coefficients are approximate. The complete ideal product is
also audited against the finite one-particle free exponential.

That exponential is evaluated by a degree-32 tridiagonal Taylor recurrence,
with a Hermitian remainder and positive row-norm rounding propagation. For
the full sweep product, positive row sums and row-error bounds are propagated
through every plane rotation. The direct difference matrix is enclosed by
sqrt(||D||_1*||D||_infinity), with subtraction, reference and product arithmetic
added. A separate phase sum bounds the ideal-unitary versus floating product.
The audited operator error is lifted to k<=M photons by the telescoping bound
k*delta*(1+delta)^(k-1). No sampled-state comparison substitutes for it.

## Schmidt truncation, SVD checks and canonical environments

Each two-site tensor contains all allowed physical charges at its outer
bond charges. Its SVD is performed separately by the intermediate conserved
charge, using LAPACK `zgesdd`. The positive norm of discarded singular values
is bounded with the full U/V Gram defects. The *direct* full factorization
residual `X-(U*s)*Vh` is separately measured positively. It is not inferred
by subtracting total and retained squared norms. Both full factor matrices'
isometry defects are checked; a retained subset inherits their bound.

Left/right canonical sites have stored operator Gram bounds eta_j. The
outside environments have amplification at most 1/(1-sum eta_j), provided
sum eta_j<1e-6. This is checked, not assumed from an ideal SVD. Integer
charge-sector copies remove no physical occupation. Cache reuse stores only
index layouts, never state amplitudes or projected physical dynamics.

The per-gate arithmetic allowance
gamma(8192*(largest_factor_dimension+d1*d2+96+1))
times environment_amplification*(1+center_norm) covers the explicit complex
center contraction, gate multiplication, residual evaluation, new-core
scaling and positive reductions. Absolute-product Frobenius bounds control
cancellation. This overcounts the small charge-block products; the actual
factorization residual and gate-operator errors are still added separately.

The only physical information discarded is the paid Schmidt tail. The
history register is never dephased or normalized component by component.
Any skipped gate acts on an exactly vacuum stored suffix: the incoming
bond charge M forbids photons there. It is the identity on that stored
state, not a declaration that omitted amplitude can never return.

## Whole-state ledger and final probabilities

Exact unitary reference gates telescope the local state errors, including
all later returns. Physical Strang, finite-field factorization, Schmidt
truncation and arithmetic components are distinct. Intermediate physical
norms must remain <=1.01. One global normalization is applied per physical
step, with its actual state change and division rounding charged.

At the final boundary center, if s=sum eta_j, squared norms are enclosed
by center squared norms times [1-s,1/(1-s)]. Formation event intervals use
these environment factors for both numerator and total norm. Final norm
defect is paid before the normalized reverse/complementary-event probability
inequalities. The exact time-grid mismatch is charged with
||H||<=(N-1)*(2+max|delta|)+g. The analytic physical-receiver cut cost is
included, not hidden in state compression.

Complete error must be <=.60, with the separate frozen ceilings in CHARTER.
An incomplete target earns no target-time probability claim even when its
intermediate populations or compressed representation look favorable.
