# B2-RA fixed-path operating theorem — self-contained statement

**Status:** author-derived at 38fa6a3, independent mathematical review
pending. This statement consolidates that result; it does not strengthen
its validation status or change its law. Units: hbar=1.

## Complete model

Fix L in 24N, L>=2016. The geometry is one PREPARED, FIXED path with
events x=0,...,L-1 and only consecutive occupied bonds. These ranks are
coordinates for a given record, not a theorem of dimensional selection.
Geometry does not evolve in this theorem.

At each event there is a hard-core plus slot with states vac,0,1 and a
hard-core minus slot with states vac,occupied. A plus carries one qubit;
a minus has no internal qubit. At each endpoint there is an additional
neutral hard-core fuel slot with states vac,0,1. Different events and
species commute. No fermionic exchange signs are supplied. Electric links
have fixed capacity E_x in {-1,0,1}, positive energy E_x^2, and Gauss law.

An explicit scalar basis is (P,M,F), where P,M are subsets of the events,
F is a subset of {left,right}, and

    k=|P|=|M|=2-|F|,
    E_x=sum_(y=0)^x [1_(y in P)-1_(y in M)], x=0,...,L-2,
    |E_x|<=1.

All such configurations are retained, including separated constituents.
Opposite species may co-locate. Set exactly

    jp=1/64, jm=1/96, kappa=1/4,
    A=(jp^2+jm^2)/(5/4)=13/46080,
    J=lambda=2jp jm/(5/4)=1/3840, Omega=2-2A.

The scalar Hamiltonian H_L has diagonal

    2k+Omega|F|+sum_x E_x^2
      +kappa sum_(x=0)^(L-2)
          [1_(x in P)1_(x+1 in M)+1_(x in M)1_(x+1 in P)]
      -A[1_(0 in P intersect M)+1_(L-1 in P intersect M)].

Its only off-diagonal entries are:

1. Move one plus (minus) by one path edge to an empty slot of the same
   species, with matrix element -jp (-jm), provided the resulting field
   still satisfies the displayed capacity and Gauss constraints.
2. At either enabled endpoint, convert an occupied fuel into a co-located
   plus and minus, when both matter slots there are empty, with matrix
   element -J; retain the exact reverse conversion, which requires an
   empty fuel. A neutral co-pair changes no electric field.

There are no other entries. The matrix is finite and Hermitian. Its
colored lift transports the plus color with a hop and transfers the color
between fuel and plus during conversion; all coefficients are color-blind.
No virtual state, return transition or collision is deleted.

The conserved resources are N_++F=N_-+F=2, F the number of excited fuels.
The two color carriers occupy the extended line
fuel_left,event0,...,event_(L-1),fuel_right. Its allowed moves cannot pass
one another. Sorting them left to right therefore proves, rather than
postulates, H_colored=H_L tensor I_4 on this sector.

## Preparation, instrument and observables

At t=0, P=M=empty, E=0 and both fuels are excited. Their colors encode
an arbitrary density matrix rho_ABZ, including any external reference Z.
The spatial preparation is independent of that input. The geometry record,
fuel loading and final readout are supplied. Initial added mean energy is
2Omega and variance is 2J^2; subsequent evolution is time independent and
unitary. There is no later switch, reset, favorable measurement or sink.
Any decoupled spectator/environment can be retained and purified.

Define

    I_L=[L/(8J),L/(6J)],
    X_L={L/8,...,3L/8-1}, Y_L={5L/8,...,7L/8-1}.

The success projection R_L requires empty fuels and P=M={x,y}, with
x in X_L,y in Y_L. Sum over all allowed unmeasured final configurations;
read plus colors in left/right order and discard the spatial registers.
The output instrument flags every other outcome as erasure. On success,

    E_success,t(rho_ABZ)=p_out(t) rho_ABZ.

The failure output may be represented as (1-p_out)|erase><erase| tensor
rho_Z after discarding the unread carriers. Success probability and
conditional channel quality are separate. The latter is exactly zero
trace-distance error, with trace distance defined as half the trace norm.

Retention is the event k=2 (equivalently F=0), including ALL admissible
separated-charge states. It is not restricted to co-located composites.

## Uniform conclusion

For the exact model and original preparation above,

    p_out(t)>514089/4000000=0.12852225, every t in I_L,
    p_ret(t)>434281/1000000=0.434281,
                              every 2/J<=t<=L/(6J).

The statements are uniform one-time measurement probabilities, not a
claim about an unmeasured classical trajectory or repeated measurements.
They meet the frozen 10% joint-output, 25% retention and 1% conditional
information-error task. No parameter interval or disorder class beyond
the exact nominal values is asserted.

## Complete comparison account and proof interface

The analytical source is h on C f + ell2(N_0^2): h_ff=Omega,
h_f,00=-J, diagonal 2+|i-j|+(1/4)1_(|i-j|=1)-A1_(i=j=0), and constituent
hops -jp,-jm in their respective coordinate directions, with adjoints.
This half-line operator is a reference; the physical finite matrix is H_L.

| Step | Exact identity or paid bound | Domain |
|---|---|---|
| Actual source spectral kernel | Closed-arm divided-difference factor; weighted boundary inverse; fuel denominator differs from -J exp(-ik) by <J/100 | Complete first bulk band, thresholds included |
| Source spectral accounting | k-density between 24/25 and 26/25 times (2/pi)sin^2 k; all unclassified mass <=1/25 | Original fuel seed; no physical filtering |
| Source return | abs(G(t))<=2/25+1/(J abs(t)), hence <=29/50 for J abs(t)>=2 | All such half-line times |
| Source capture | Orthogonal ac-central weight >4/5; complete missed squared norm <3/20, including fuel, separated states, both phases and global local remainder <1/200; coherent complement paid | Every t in I_L; gives A_L(t)>9/25 |
| High-source spectral tail | omega<=36/27613^2, sigma=sqrt(omega)<=6/27613 | Conserved in exact source evolution |
| Source truncation | M=L/2-2; tau_L<1.01 exp(.3-L/105)<1/1000 | 0<=t<=L/(6J) |
| Finite two-source state comparison | delta_L=12/27613+1500 exp(-L/105)<1/1000 | Same physical Hilbert space; all collisions, reflections and returns retained |

Let C_M retain fuel and all |i,j> with max(i,j)<=M, and let R_X retain
the co-pair sites |x,x> with x in X_L. For
u=exp(-ith)1_(-infinity,2.5)(h)f, set a=||R_X u||^2, g=<f,u> and
m=||P_matter C_Mu||^2. Exact reconstruction gives

    a>=A_L-2sigma, |g-G|<=omega,
    m=1-omega-|g|^2-tau_L^2,
    p_out>=[a-delta_L]_+^2,
    p_ret>=[m-delta_L]_+^2.

The joint reference output is a^2, NOT a. The source-filter price 2sigma
and the initial mismatch already inside delta_L are both conservatively
paid. There is no additional unlisted fixed-path approximation error.
These inequalities yield a-delta_L>717/2000 and m-delta_L>659/1000,
whose squares are the stated guarantees. The exact successful color
channel is a separate no-passing identity for this SAME output projection.

For a review of the proofs behind these bounds, use REVIEW_DOSSIER.md.
The 110 archived exact-arithmetic checks verify constants/finite identities;
they are not independent review of this analytic chain or a physical test.

## Explicit boundary

No moving-geometry error has been paid here: geometry is fixed by the
theorem's domain. Forming the host, recording it, internal fuel sites,
extensive density, autonomous preparation, a continuum limit, gravity and
empirical discrimination are outside this theorem. The original eight-event
B2 failure and earlier intermediate limitations remain unchanged.
