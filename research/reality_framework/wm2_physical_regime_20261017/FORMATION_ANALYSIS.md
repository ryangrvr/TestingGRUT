# Formation analysis — exact history route, before checks

This records a preliminary analytic route at freeze, not a proved claim.
The observation protocol/class are fixed in PHYSICAL_TARGET.md. Final
derivation, constants, checks and remaining gaps will be recorded here.

Fix a root, M=N-1. Let S_v append unused event v to a word, and let
S=sum_v S_v. On the complete history Hilbert space,

    S_v^dagger S_w=delta_vw P_(v free), S^dagger S=M-K.

No subset/radial reduction is used. The relative formation law is

    H=D+dGamma(h_L)-lambda(S a_0^dagger+S^dagger a_0),
    lambda=g/N, D=sum_v delta_v P_(v attached).

h_L has off-diagonal -1 and zero diagonal on the physical L=64N^2 chain.
Its physical frequencies are 4+sigma(h_L); the relative law subtracts a
conserved charge. The seed is the empty word and receiver vacuum.

Investigate the exact Heisenberg receiver equation

    a_0(t)psi_0=i lambda integral_0^t C_L(u) S(t-u)psi_0 du,
    C_L(u)=<0|exp(-i h_L u)|0>.

The half-line endpoint has semicircle density
rho(x)=sqrt(4-x^2)/(2pi), |x|<2. Direct integration by parts with endpoint
cutoffs appears to give |C_infinity(u)|<=min(1,2u^(-3/2)). Finite-chain
closed-walk tails may control C_L-C_infinity throughout t<=8N^2.
If so, contact occupation is O(1/N) uniformly over history and disorder.

Next remove only the explicit detuning phase from the retarded S. The
interaction commutator includes **both** [S,A]a_0^dagger and
[S^dagger,A]a_0. The first is generally nonzero for ordered histories;
it cannot be deleted by analogy with independent emitters. Bound both
using the actual contact occupation and ||A||<=sqrt(M).

The proposed comparison drift is

    d<K>/dt=(g/N)^2 sum_v sqrt(4-delta_v^2)<P_(v free)> + r_N(t).

It is not an assumed irreversible law. Derive and bound r_N from the
full closed Hamiltonian, including initial-slip, finite-chain and coherent
return terms. An O(N^(-5/4)) current remainder plus an integrable
initial-slip term would suffice for an O(N^(-1/4)) macroscopic error on
the frozen N^2 time scale. This scaling is a target to check, not a
conclusion at freeze.

Attempted counterexamples must retain the actual receiver and detuning
class. Band-edge detunings outside |delta|<=3/4, an empty fuel sector,
or a closed small receiver are different models. A failed proof estimate
does not refute physical formation.
