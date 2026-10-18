# B2-BB: actual bulk symbol and small boundary kernel; operation unresolved

**The actual bulk-plus-boundary decomposition is derived. The frozen
operating gate remains unresolved.** All new mathematical claims are
author-derived with independent external proof review pending.

The user's derivative and monotonicity bounds check:

    Sigma(E)<=0, Sigma'(E)<=0,
    ||Sigma(E)||<=25/6432,
    ||Sigma'(E)||<=25/4489, E<=2.5.

The proposed 36/133010089 source-tail bound is also correct. It is weaker
than B2-SO's 36/27613^2, which can be rederived using the closed-sector
spectral separation. The stronger original bound stays in the joint
probability account; no parameter or error budget is relaxed.

## Actual spatial reduction

The whole-line extension of the SAME two-charge law gives the exact
translation-invariant bulk symbol

    v(k)=jp^2+jm^2+2jp jm cos k,
    sigma(z,k)=2v(k)/[z-13/4-v(k)/[z-4-v(k)/[z-5-...]]].

Below 2.5 there is one bulk low branch, determined by E-2-sigma(E,k)=0.
The infinitely many transverse states remain in the continued fraction.
Bulk propagation does not establish the original fuel's launch weight.

The ACTUAL quadrant self-energy is exactly

    Sigma_+(z)=T_+(sigma(z))+W(z),
    |W_nm(z)| <=(75/20636)(5/72)^(2max(n,m)), Re z<=2.5.

This is derived by retaining all virtual words and isolating those that
leave the quadrant. It is not an assumed finite-rank source defect.
The root counterterm gives U=W-A|0><0|. Its leading word cancels at
the reference energy 2; the energy mismatch and all higher words remain.
The explicit two-hop/four-hop coefficients and complete remainder prove

    ||U(E)||_1 <= (1506254603/840983040000)J
               <J/500, E in [2-5J,2].

This interval contains the known first bulk band. It does not contain
every unresolved source component by assumption. The norm is trace norm
of the co-pair boundary correction; the order-J fuel coupling is still
retained separately and exactly.

## Source-to-receiver inverse with phases retained

For fixed Im z>0, the ACTUAL scalar symbol kappa=z-2-sigma admits an
explicit Fourier-log factorization. Its triangular Toeplitz factors give
a spatial kernel g_T(n,m), with no replacement by 1/Im z. Then exactly

    S=[I-g_T U]^-1 g_T,
    F=[z-Omega-J^2<0,S0>]^-1,
    F_n=-J F<n,S0>.

[SOURCE_KERNEL.md](SOURCE_KERNEL.md) provides the factor coefficients
and their exact convolution formula. These equations retain the actual
boundary, fuel, physical co-pair detector and propagation phases.

Their uniform real-axis limit and coherent Fourier inversion are still
unproved. In particular ||U||_1<J/500 does not imply an inverse bound
on the propagating band. Source poles, exceptional singular measures,
edge behavior and the fuel denominator still have to be included or
bounded. No new arrival or late-retention inequality has been established.

## Verification and decision

All 66 frozen exact-rational checks passed; the one replay was byte
identical. Actual primitive virtual paths of lengths two and four agree
with the continued-fraction coefficients and boundary matrices at the
declared energies. Rational arithmetic verifies the derivative/tail and
trace-norm margins. The analytical path-length estimate supplies the
infinite-size remainder, not a finite grid of host lengths or times.

The child workers used 0.061856 seconds combined and 10.125 MiB peak RSS.
No dynamics, spectral quadrature, physical operating probabilities,
external reviewer or independent implementation ran.

Stop at the exact remaining theorem in [REVIEW_PROBLEM.md](REVIEW_PROBLEM.md):
real-axis factor/scattering control and fuel denominator estimates strong
enough for the coherent receiver event on n~L,t~L/J. The existing finite-host
comparison still pays reflections and returns; the correct 34% one-source
to >11.458225% joint implication remains conditional. B2's failed finite
endpoint, the unchanged law and the original gates are preserved.
