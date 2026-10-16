# Exact action-directed enrichment; separate from the primary prototype

The primary S6/G1 trajectories use the neighbor-density heuristic frozen
in their charters. This note derives a more informed numerical basis
construction after observing that heuristic's limitations. It is checked
only as an algebraic module here. It is not substituted into those
trajectories, and supplies no new formation number or gate pass.

## One-body density of a coherent sum of receiver packets

Inside a fixed target history, write the complete Hamiltonian action as

\[
 y=\sum_a F_k(W_a)z_a.
\]

For a removal contribution, z_a is the annihilated child field. For an
append contribution, temporarily enlarge its parent frame by contact
and create a photon there. This is an exact finite representation,
not a projection. Include every actual coupling and phase in z_a.

Let D_a have columns a_l z_a in the (k-1)-photon occupation basis of
W_a. With Gamma the complete rectangular Fock contraction, define

\[
 M_{ab}=\left[D_b^\dagger
       \Gamma_{k-1}(W_b^\dagger W_a)D_a\right]^T.
\]

Then the exact physical one-body density, with
gamma_rs=<a_s^dagger a_r>, is

\[
 \gamma(y)=\sum_{a,b} W_a M_{ab} W_b^\dagger=Z M Z^\dagger,
 \qquad Z=[W_1\ \cdots\ W_s].
\]

This follows by annihilating one photon from each side of the coherent
sum before taking their complete cross-frame inner product. In
particular tr gamma=k||y||^2 and gamma is positive. The off-diagonal
blocks are essential. Separate incoming densities omit their physical
interference and can select a different enrichment direction.

No L_N by L_N matrix is needed. QR of Z reduces diagonalization to at
most sum_a r_a dimensions. For the leading direction outside a current
frame W, use Z_out=(1-WW^dagger)Z and diagonalize its reduced covariance.
The cross-Fock contraction here has k-1 photons, one fewer than a
direct coherent field overlap. Exact Fock interference is retained.

For a one-photon target, the outside covariance has rank at most one.
Adding its nonzero eigenvector captures the entire current Hamiltonian
action in that target. At general k, adding every positive-density
orbital also captures y exactly: zero expected occupation outside the
chosen support implies vacuum outside it. There is no uniform bound
on the number of such orbitals.

## A systematic Krylov enclosure

There is a stronger sufficient construction. At the start of a B step,
let Phi be the actual represented state. If enlarged history frames give
projector P' with

\[
 P'B^l\Phi=B^l\Phi,\qquad l=0,\ldots,q,
\]

then the first q complete residual coefficients vanish. Indeed, writing
C'=V'^dagger B V', induction gives V'C'^l x'=B^l Phi for l<=q.
Consequently (1-P')BV'C'^l x'=0 for l<q. The physical B-step error is
then bounded by the full uncomputed tail

\[
 \rho_h\le e^{gh}(gh)^{q+1}/(q+1)!.
\]

This includes every longer virtual word and return. It is not a
truncated physical Hamiltonian. One can build the enclosure from the
supports of the exact action densities of Phi,B Phi,...,B^q Phi.
Alternatively, an explicitly sufficient but potentially expensive
construction is

\[
 W'_u=\operatorname{span}\{c,W_v:
       \operatorname{dist}_{\rm word}(u,v)\le q\}.
\]

Each primitive B transition changes the word by one append/pop and
adds/removes only a contact photon. Therefore every contribution to
(B^l Phi)_u for l<=q lies in the displayed span. This word-neighborhood
is solely a numerical stencil in configuration space; it supplies no
new physical spatial metric.

For q=3, T=4N^2 and n=O(N^3) physical-H Strang steps, h=O(1/N).
The summed residual tail is at most

\[
 n\rho_h\le e^{gh}g^4 T h^3/24=O(N^{-1}).
\]

This is a conditional error-scaling theorem for exact Krylov enclosure.
It does not prove a scalable solver. If the maximum word degree is
Delta and initial ranks are at most r, a crude enclosure rank is
min(L_N,1+r sum_(l=0)^q Delta^l). Fock dimensions and the full history
count can still be enormous. The tested rank-ten heuristic does not
satisfy this exact enclosure hypothesis; its recorded nonzero residuals
must not be replaced by this tail-only estimate.

The new module checks coherent action densities, complete cancellation,
one-photon action capture, and the Krylov support identities against
independently lifted finite physical matrices. These are algebraic
checks, not a later change to the frozen formation protocols, not a
proof of sufficient small ranks, and not an all-size formation result.
