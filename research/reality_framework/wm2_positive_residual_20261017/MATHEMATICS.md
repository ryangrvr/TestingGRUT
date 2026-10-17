# Shared-frame action enclosure and a positive full residual bound

This construction changes the solver, not S4's physical Hamiltonian.
Use S5's exact ordered-type-word representation in probe vacuum, with
H=A+B and k=|u| receiver photons in word u. Root components are averaged
with their original weights because root is invariant and the formation
event is block diagonal; this does not dephase or select the preparation.

## Complete B enclosure with one common frame

Let W be an orthonormal common frame, V=direct_sum_u F_|u|(W), P=VV^dagger.
If physical contact c lies in span W, creation and annihilation at c
preserve this complete Fock space. Hence

\[
 (1-P)BV=0,\qquad P B^l V=B^l V\quad(l\ge0).
\]

This is the shared-group special case of S6's exact Krylov enclosure.
It includes all word histories and interference; common coordinates
make their full Fock overlaps identity. It inserts no physical metric.

For imperfect capture set c_perp=(1-WW^dagger)c. Annihilation never
leaves the span. Append creates one outside photon, with no outside
bosonic enhancement because that mode begins in vacuum. Every child
has a unique parent. Summing exact type multiplicities gives

\[
 \|(1-P)BV\|\le {g\sqrt{N-1}\over N}\|c_\perp\|.
\]

Indeed the squared outgoing norm from a k-photon parent block is
(g/N)^2 (N-1-k)||c_perp||^2 times its squared norm. The maximum is k=0.
Orthogonal target words pay all children without a word-count factor.
This proof includes both directions of the unchanged Hamiltonian.

For any full-column-rank numerical frame W_f and any coefficient z,

\[
 \operatorname{dist}(c,\operatorname{span}W_f)
 \le\|c-W_fz\|.
\]

Thus a directly enclosed positive vector norm bounds leakage, without
assuming W_f is exactly orthogonal and without squared-norm subtraction.
Let W_bar be its exact polar isometry. With C=V_bar^dagger B V_bar,

\[
 \|e^{-iBh}V_{\rm bar}x-V_{\rm bar}e^{-iCh}x\|
 \le |h|{g\sqrt{N-1}\over N}\sigma_c\|x\|,
 \qquad \sigma_c\ge\|c-W_fz\|.
\]

The physical unitary Duhamel identity pays every future return. If a
rank cap prevents enrichment, sigma_c remains in the budget. Nothing
allows the exact-enclosure/tail-only conclusion when this term is nonzero.

## Finite Hermitian polynomial and splitting

For a Hermitian coefficient operator C_num with ||C_num||<=G, the
unitary integral Taylor remainder gives

\[
 \|e^{-ihC_{\rm num}}-\sum_{j=0}^3{(-ihC_{\rm num})^j\over j!}\|
 \le (G|h|)^4/24.
\]

There is no exponential amplification factor: the integral contains a
unitary commuting with C_num^4. Code coefficients/contact approximating
C are paid separately by |h| ||C-C_num||. Arithmetic error in polynomial
evaluation and subsequent normalization is also separate.

Keep the physical-H Strang constant from S5,

\[
 C_S={g(d^2+2d+\sqrt2)\over24}+{2g^2(d+1)\over12},\quad d\le3/4.
\]

The global splitting bound is C_S T^3/n^2. All actual time-grid and
phase rounding is accounted separately. Only representation changes,
not physical terms, are compressed.

## Projection and polar-frame bookkeeping

For normalized abstract coefficients, an orbital-vacuum projection with
discarded global probability q changes the globally normalized state by
sqrt(2q/[1+sqrt(1-q)]). Positive discarded and retained coefficient norms
can instead bound the actual projection/renormalization directly. No
history is normalized independently. Every eigenvector is a numerical
proposal; accept removal only using the enclosed discarded norm.

If eta>=||W_f^dagger W_f-I||<1, the polar factor differs by at most
eta/(1+sqrt(1-eta)) in one-particle norm. Fock lifts on k<=N-1 particles
are Lipschitz with the telescoping bound
k delta max(||W_1||,||W_2||)^(k-1). Bound this conservatively using eta
before/after enrichment, rotation, projection or reorthogonalization.
Do not silently replace the numerical frame by its polar factor.

Exact integer coefficient permutations introduce no rounding. For a
two-mode proposed rotation T, enclose T^dagger T-I and its finite Fock
lifts, and charge coefficient/frame arithmetic and nonunitarity. Frame
Gram checks decide whether the declared envelope is satisfied.

## Arithmetic scope and final probabilities

The numerical implementation must specify directed scalar enclosures,
positive reduction bounds, exact-constant/phase enclosures, matrix-vector
rounding bounds and basis transformation costs before a primary claim.
That finite arithmetic model is distinct from formal verification of
the compiler/BLAS stack and external mathematical review.

Given normalized computed probability p_hat in [p_lo,p_hi] and complete
state-norm error e, use monotonic outward evaluation of

\[
 p_{\rm phys}\ge\max\{0,(\sqrt{p_{\rm lo}}-e)_+^2,
               1-(\sqrt{1-p_{\rm lo}}+e)^2\}.
\]

The corresponding upper bound is the minimum of 1,
(sqrt(p_hi)+e)^2 and 1-(sqrt(1-p_hi)-e)_+^2. Formation and full assembly
are distinct projectors. Incomplete target times have no such new claim.

This calculus does not prove affordable rank growth or formation. It
makes the action enclosure and its arithmetic testable before numbers
are promoted to physical results.
