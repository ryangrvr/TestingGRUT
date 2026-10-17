# Arithmetic contract and error ledger

S7 is a roundoff-aware implementation under a declared floating-point model.
It is not a formally verified compiler/BLAS stack or an independently validated
numerical theorem. Small controls test the implementation; they do not prove
every hardware or library assumption.

## Model

Binary64 scalar arithmetic uses round-to-nearest, unit roundoff `2^-53`.
Native positive sums use C++ `long double` with its reported significand width
(64 bits on this run); conversion back to binary64 is charged separately.
`native.py` checks the rounding mode. Native compilation disables fast-math and
contraction. Subnormal/underflow contributions are padded by the smallest
*normal* binary64 number times an operation count, a deliberately loose bound.
No numerically observed negative residual is clipped to obtain a certificate.

For an n-operation reduction use

\[
\gamma_n=\frac{nu}{1-nu},\qquad
|\operatorname{fl}(\sum a_k)-\sum a_k|
 \le\gamma_n\sum|a_k|.
\]

The BLAS products used for contact coordinates and Gram matrices are assumed
to obey this standard absolute reduction bound (reordered or fused operations
are allowed when they obey it). NumPy's complex128 buffers must have two
binary64 components. The ctypes native struct has a matching size assertion;
matrix/adjoint controls check the actual buffer interpretation. A violation of
these runtime/library assumptions is outside the enclosure claim.

## Constants, phases, positive norms

`Fraction` supplies exact rational endpoints and endpoint tests for scalar
square roots. Machin's identity with rational alternating-series remainders
encloses pi. Interval Horner evaluation of 16 sine/cosine terms, with an
outward order-32 Taylor remainder, encloses the physical receiver contact,
all finite-chain frequencies and every applied phase. No library sine/cosine
result is treated as an exact physical constant.

`positive_norm` sums the *squares of the input components*, using extended
precision. `positive_raw_bounds` encloses this with a four-operations-per-entry
bound plus conversion and underflow terms. Vacancy probabilities use the
same positive accumulation on distinct occupied modes, not number-weighted
occupation. A computed normalized event interval divides its positive event
norm enclosure by the total norm enclosure. Formation and full assembly use
different coefficient slices.

For contact leakage, directly form `c_float-W*z` and add bounds for the
complex dot products, subtraction and contact-constant uncertainty to its
positive norm. This encloses distance to the *span* of W regardless of its
orthogonality. For coordinates in the exact polar frame, separately pay
the Gram/polar discrepancy and contact dot-product error. These are different
costs; neither is inferred by subtracting two almost equal squared norms.

## Frame and coefficient transformations

An outward row-sum bound for `W*W-I`, including complex reduction errors,
encloses its operator norm eta. Require eta<=1e-9. The polar-distance bound
from MATHEMATICS.md is conservatively replaced by eta. Before/after frame
changes pay the Fock telescoping bound, using `2*(N-1)` as a small-eta safety
factor. Free evolution also pays the exact-phase enclosure and rounded
complex products. A numerical enrichment vector need not equal an exact
Gram-Schmidt vector: keeping the old columns preserves their span, while
the new Gram and contact checks determine the actual accepted representation.

Two-mode Fock matrices use S5's finite normalized polynomial, not an
unverified matrix exponential. There are at most six input occupation states
per block in this campaign. `fast_matrices` bounds the absolute sum of its
polynomial terms *before cancellation*, then uses gamma(4096) per entry and
row/column bounds. This overcounts the explicit complex operations, integer
binomials, powers and square-root factors for n<=5. `rotation.matrices`
independently evaluates their full rational/complex interval polynomials in
the controls; its smaller bounds must lie below the fast bound. The accepted
bound is not based on their agreement. Nonunitarity of the proposed two-mode
matrix, native block products and the numerical physical-frame update are
charged in `basis`.

Eigenvectors and eigenvalues propose a mode to remove. Acceptance always
uses the *actual positive discarded-probability enclosure* after the proposed
rotations. Integer embedding, restriction and mode permutations copy entries
exactly. The retained vector is normalized globally; input norm deviations
and division rounding are charged separately. No history is normalized alone.

## Emission polynomial and remaining components

The ideal coefficient operator using the declared numerical coefficients is
Hermitian. Rational/square-root endpoint comparisons bound coefficient errors
relative to the physical g/N law. Its norm is bounded using the complete
absolute-entry emission bound, not a sampled vector norm. The Hermitian
cubic tail `(G*h)^4/24` belongs in `polynomial`.

The native sparse action has at most `(number of types)*rank+N-1` incident
contributions per coefficient. Its explicit complex products and sums, three
polynomial actions, scalar multiplications and additions are dominated by
gamma(`4096*(rank*(types+1)+N-1+10)`) times
`norm(input)*(1+G*abs(h))^4`. Coordinate and primitive-coefficient perturbations
are separately paid by the unitary Duhamel bound in `arithmetic`.

`normalization` pays the actual vector change from scaling plus its arithmetic
enclosure. `strang` uses the unchanged physical commutator bound. `time_grid`
pays the exact rational discrepancy between Nsteps times binary64 h and the
declared T using `||H|| <= (N-1)*(max|delta|+2)+g` in the shifted sector.
`final_normalization` pays the final norm defect before applying normalized
event inequalities. All error additions and scalar interval endpoints are
rounded outward; the whole ledger is summed with an additional gamma bound.

Only a completed target time with every gate satisfied earns a physical
probability interval. Stopped endpoints and their error ledgers are retained,
but earn no target-time formation claim. A large evaluated enclosure remains
a large enclosure even when a reference comparison has a much smaller error.
