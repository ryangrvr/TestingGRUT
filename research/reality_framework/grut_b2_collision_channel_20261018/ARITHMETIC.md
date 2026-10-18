# Evaluated propagation-error account

This account is conditional on the stated arithmetic model. It is not an
outward-rounded interval certificate or validation of the numerical stack.
The code records every computed residual/orthogonality term separately.

For the exactly specified rational physical matrix H, let computed real
spectral data be V,D. Define R=HV-VD, delta=||V^*V-I|| and

    y(t)=V exp(-itD)V^*|source>.

Exact Duhamel, without assuming V unitary, gives

    ||U_H(t)|source>-y(t)||
       <=delta+|t| ||R|| sqrt(1+delta).

Here ||VV^*-I||=||V^*V-I|| because V is square. Computed spectral accuracy
is therefore assessed after the eigensolve; no stability theorem for that
eigensolver is assumed. Frobenius norms upper-bound the needed operator
norms. The source is an exact basis vector, so V^*|source> is one stored
matrix row, not a separate cancellation-prone inner product.

Arithmetic contract: binary64, unit roundoff u=2^-53, no overflow/underflow,
dot products and matrix operations satisfying the conventional gamma_n=
nu/(1-nu) absolute-error model, correctly represented input integers,
and complex exponential evaluation within the conservative phase allowance
32u(1+max|D|*|t|). This includes multiplication of D by t and assumes an
appropriate library range-reduction error bound. The contract is stated,
not empirically proved for every hardware/software component.

Positive conservative additions used by the code are

* delta <= computed_F(V^*V-I)
  +2 gamma_(4n+16) (||V||_F^2+sqrt(n));
* ||R|| <= computed_F(H_float V-VD)
  +2 gamma_(4s+16) (||H_float||_F+max|D|)||V||_F,
  where s is maximum row nonzero count;
* ||H-H_float|| <=128u(1+max row absolute sum of H_float), for this
  bounded-entry rational builder. Integer/quarter field, mass and density
  entries are exact; the nonbinary fuel/counterterm coefficients and their
  few additions are covered by this deliberately loose row bound;
* phase and final matrix-vector contributions are added positively, including
  2 gamma_(4n+16)||V||_F sqrt(1+delta).

The factors of two cover the smaller positive norm-evaluation roundoff at
these frozen dimensions (n<=773). No residual is computed by subtracting
two squared norms, and no tiny negative square is clipped. The displayed
evaluated error is used in projection-amplitude enclosures

    max(0,sqrt(p_hat)-e)^2 <= p <= min(1,sqrt(p_hat)+e)^2.

These do not enclose arbitrary hidden implementation errors. Gates use
the enclosures, not bare midpoint populations. If the arithmetic contract
cannot be established by a future implementation, report numerical
agreement and an unvalidated error account rather than a certified theorem.

The two-qubit channel quality is derived algebraically, not estimated from
one favorable input's fidelity. Its sole numerically estimated parameter
on the endpoint path is the joint outgoing probability p_R.
