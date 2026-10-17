# Pre-implementation construction

S8 changes the receiver representation only. Standard MPS/Schmidt methods
are numerical tools, not new fundamental laws or a derivation of quantum
mechanics. Full derivations and arithmetic details accompany the final packet.

## Charge and receiver cut

In a word u of length k, M-k free events plus k receiver photons have
conserved total M=N-1. The boundary assisted transition trades one free
event for a photon at receiver site zero. Define the positive weighted charge

K_mu=(M-k)+sum_r exp(2*mu*r)n_r.

It commutes with detuning and the full boundary emission/removal term.
Receiver hopping obeys |d<K_mu>/dt|<=4*sinh(mu)<K_mu>, hence from the
neutral vacuum seed <K_mu(t)><=M*exp(4*sinh(mu)t). Boundary Duhamel for
a receiver cut after R sites, including every return in the original law,
gives

||U_L(T)psi0-U_R(T)psi0||
 <=sqrt(M)*(exp(c*T)-1)/c*exp[-mu*(R-1)], c=2*sinh(mu).

Set mu=1, use 2<c<12/5 and R>=12*T/5+32. Since e>2, the bound is
less than sqrt(M)/2^32. If no physical cut is made it is zero. This paid
embedding retains the original physical L_N receiver; it is not an exact
small-bath substitution or a photon cutoff.

## Local commutators are not truncated dynamics

F is the full free receiver hopping and G is word detuning plus the
unchanged contact emission/removal. [G,F] acts only on receiver modes
zero and one. [G,[G,F]] uses those same modes; [F,[F,G]] additionally
uses mode two. Every more distant mode is a spectator in these operators.
Decompose by spectator photon number r=0,...,M and retain *all* local
occupations with n_local=|u|-r. Each finite matrix is an exact block of the
full commutator algebra, not the Hamiltonian used for propagation. Taking
the maximum certified block norm bounds the complete physical operator.

For Hermitian K and a computed eigenvector matrix V/eigenvalues lambda,
eta>=||V*V-I||<1 and e>=||KV-V diag(lambda)|| imply

||K|| <= (max|lambda|*sqrt(1+eta)+e)/sqrt(1-eta).

Pay matrix-construction and multiplication errors separately. This uses all
eigenvalues and a residual, not a hand-picked low-energy sector.

## Tensor truncation and returns

The history site has charge M-|u|, receiver site r has charge n_r and total
charge is M. In a sweep, all sites outside the two-site center have checked
left/right isometry defects eta_j. Their environments have norm bounded by
products of (1+eta_j); a conservative factor is 1/(1-sum eta_j).

At a center SVD, the positive discarded singular-value norm, direct full
factorization residual and checked Gram bounds upper-bound the physical
discarded state norm after that factor. Multiplication/gate errors are paid
as well. Exact charge-sector copies and assignments remove no occupation.
Only one global normalization is used. Sum local physical-state errors by
telescoping the exact unitary reference gates; every possible future return
is paid, rather than declaring a discarded photon permanently gone.

Receiver factors approximate exp(-iFh) in the finite prefix. If the audited
one-particle error is delta and its norm is <=1+delta, its k-photon Fock
lift differs by at most k*delta*(1+delta)^(k-1), k<=M. Include numerical
Fock-gate coefficient errors separately. No sampled vector norm substitutes
for this complete operator account.
