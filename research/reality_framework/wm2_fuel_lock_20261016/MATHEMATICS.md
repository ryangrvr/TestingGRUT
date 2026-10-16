# S2: resonant formation, complete protection and the remaining failures

These are statements about the declared auxiliary S2 invariant component,
at most one probe, and exact resonance unless explicitly stated. No
uniform finite-density, original-charged-matter or all-WM2 result follows.

## 1. Exact component and energy-defined band

The finite-register append/pop rules leave a rooted word w with distinct
labels, just as in S1. With a probe all 2^(N-1) non-root storage patterns
are retained for every word and probe position. This is not a storage
cutoff. In probe vacuum the unsynchronized interaction is identically
zero, so the initially consistent component stays exactly consistent.

Completeness with one probe also has a constructive proof. At the root,
assisted append followed by bare pop flips a fresh storage bit to spent;
bare append followed by assisted pop flips it back. Keep the probe at
the root for these two-move cycles to set every non-root storage bit
arbitrarily. Then build any target word using bare appends while walking
the probe with the head, preserving those bits. Walk it back to any
attached target position. For a free target position, append that event,
hop the probe onto it and bare-pop it again. These reversible sequences
reach every listed word/storage/probe configuration; root charge cannot
change. The finite graph checks do not carry the all-N proof alone.

At delta=0 subtract UN. With C_i=f_i+q_i-1, the diagonal is

`D = sum_i [Lambda C_i^2 + U C_i]`.

P0=1_(all C_i=0) includes every assembly stage. A bound/fresh mismatch
costs Lambda-U=48 and a free/spent mismatch costs Lambda+U=80. Hence
D_Q>=Lambda-U=48. The protected sector is not the set K=0 and contains
the bondless seed, enabling formation within the same energy band.

Let A+ be assisted append. Every child has a unique parent including its
storage state, and each parent has at most N possible children. Its
coefficient g/sqrt(N) gives ||A+||<=g. The probe-assisted bare append
similarly has norm <=epsilon, since its one-probe occupancy factor is
at most one. Current-word matter hopping has norm <=2J. Therefore

`||T|| <= t0 = 2g+2epsilon+2J = 3`,
`||P0 H P0|| <= aP = 2g+2J = 2.5`.

No factorial configuration factor is omitted: these identities sum the
whole branching in Hilbert-space norm. B=Q0 H P0 comes only from the
unsynchronized interaction. Each one-mismatch destination identifies a
unique consistent source, so B^dagger B is diagonal and ||B||<=epsilon.

Thus QHQ>=45, while PHP has maximum <=2.5. Interlacing derives a low
band Pi of dimension dim(P0), with separating gap at least **42.5**.
For completeness let C=QHQ, L=H restricted to Pi, X=Q restricted to Pi
and Y=P restricted to Pi. The invariant-subspace equation is
CX-XL=-BY. Here C>=45 and L<=2.5. Its integral solution has norm
bounded by ||B|| integral_0^infinity exp(-42.5 t) dt. Equal subspace
dimensions then give ||Pi-P0||=||X||<=epsilon/42.5. Initial consistent
states have high-band weight at most its square. Splitting their unitary
evolution into Pi and its complement gives ||Q exp(-itH)P||<=2||Pi-P||,
so mismatch probability is at most four times that square at every time.
This is a newly derived auxiliary S2 band; the S1 gap has not been copied.

## 2. Operational distance with real assembly moves included

For a recorded word s of any length m, define d_op^s(i,j) between its
attached events by the minimum protected primitive configuration moves
from (s,i) to (s,j). Assembly and storage moves are allowed along the
route; geometry is not frozen in this definition.

Let k=m-min_t|w_t| be the largest removal depth from the initial word
along a path (zero if no initial event is removed). The initial prefix
of length m-k is immutable. If both endpoints remain in that prefix,
the number h of matter hops is at least their original separation d.
If exactly one endpoint lies in the removed suffix, d<=h+k-1. If both
do, d<=k-1. At least k geometry operations occur. Therefore every
configuration path, even with arbitrary final word and storage pattern,
has total length >=d. The static route uses exactly d matter hops.
Consequently **the operational metric equals separation along s**.

This argument proves the metric property rather than assigning distance
from a selected amplitude. Its event degree is at most two and balls in
the interior grow linearly; that dimensional class still follows from
the supplied single-end grammar.

## 3. Virtual excursions retain the stronger d+2 bound

An excursion through Q0 starts and ends consistent but contains a
mismatch. It must contain at least one pop and at least one append.
If only appends occurred, the first bound/fresh mismatch could not be
removed: its storage bit can change only when its event reattaches.
The analogous free/spent mismatch cannot be removed using only pops.

With removal depth k>=1, the number G of geometry moves is therefore
at least k+1. If k=0, G>=2. Applying the three prefix cases gives:

* both endpoints in the prefix: length>=h+G>=d+2;
* one endpoint removed: length>=h+k+1>=d+2;
* both removed: length>=k+1>=d+2.

Hence for **every** consistent initial/final geometry, with both initial
events in the source word,

\[
k_{\rm virtual}\ge d_{\rm op}+2.
\]

All storage states and all possible intermediate assembly sizes are
included. The original S1 diagonal activation barrier U d does **not**
survive: assisted pops can refund energy into free events resonantly.
S2 protection relies on this word barrier and the bounded weighted sum,
not on retaining that obsolete S1 energy claim.

## 4. Full connected Schur bound, all receiver records

For every real E<=2.5, D_E=D_Q-E>=48-E>=45.5. The exact complement is
QHQ-E=D_E-T_Q. The convergent inverse includes all repeated paths:

\[
\Sigma(E)=-B^\dagger D_E^{-1}
\sum_{r\ge0}(T_QD_E^{-1})^rB.
\]

Let P_(i,d;j) select all consistent source configurations with probe at
i, i,j attached and d_op^s(i,j)>=d. Let R_j^P include **all** consistent
receiver configurations, of every size, with probe at j. The d+2
lemma makes the first d terms vanish as an entire operator block.
Using the norm bounds for every remaining term gives

\[
\boxed{\|R_j^P\Sigma(E)P_{i,d;j}\|_2
\le\frac{\epsilon^2}{48-E-3}
\left(\frac3{48-E}\right)^d
\le\frac1{680}\left(\frac6{91}\right)^d.}
\]

This is uniform in N, all low-band energies, arbitrary coherent source
geometry packets and summed receiver records. It pays for all storage
mismatches and geometry paths. A global norm alone would not provide
the first d vanishing terms. All off-diagonal coefficients have one
sign, but no cancellation assumption is used in the operator bound.

For d>1 the direct PHP receiver block vanishes, so this is the complete
H_eff bound. For adjacent events the direct matter map is -J times an
isometry on the same geometry records. Arbitrary neighboring source
packets therefore retain receiver action >=J-1/680>0.2485. Individual
same-word nearest elements are reinforced by the same-sign correction.

PHP also contains real resonant assembly; H_eff(E) is not an
energy-independent propagator. No exponentiation of H_eff(0), nonlocal
eigenbasis definition of events, or discarded mismatch is used.

## 5. Full-H matter propagation and bulk geometry persistence

The unrestricted configuration word bound length>=d gives a full-H
Dyson estimate. In the interaction picture of D, the r-th term has
norm <=(3|t|)^r/r!. For R_j including **all** final geometries and
storage states,

\[
\|R_j e^{-itH}P_{i,d;j}\|_2
\le\sum_{r\ge d}(3|t|)^r/r!
\le e^{3e^\mu|t|-\mu d}.
\]

Thus real resonant pops have not secretly opened an instantaneous
matter shortcut. This gives a source-resolved statement too. For any
recorded word of length at least two, inject
(|s,i>+|s,j>)/sqrt(2) at adjacent attached events, and pulse exp(-i theta n_i).
The source variance is 1/4. Write A_v=R_k exp(-itH)|s,v>. The exact
receiver response is chi_k=Im <A_j,A_i>. Consequently

`|chi_k(t)| <= tail_(d(i,k))(3|t|) tail_(d(j,k))(3|t|)`,

where tail_d(z)=min(1,sum_(r>=d) z^r/r!) and tail_0=1. All final
geometry and storage records are included, including detached receivers.
In particular this is a controlled distance-dependent response front,
not just a bound on an effective matrix element.

Nearby response remains nonzero uniformly in N. The interaction-picture
Heisenberg expansion gives first term chi_j=Jt and remainder at most
exp(6t)-1-6t: each nested commutator costs at most 2||T||, and commuting
the Hermitian remainder with the source projector does not increase
its norm. The diagonal commutes with the number projectors. Thus

`chi_j(t) >= Jt-[exp(6t)-1-6t]`.

At t=0.005 this lower bound is greater than 0.000795. It holds for
every recorded word of length at least two, so the macroscopic formation
protocol below does not need a selected geometry or favorable spectral
weight for nearby response. Injection and recording remain declared
operations, not spontaneous creation of matter. Small full-H controls
include macroscopic partial records and verify the +/-pi/2 pulse identity.

There is a sharper **bulk metric persistence** bound. Fix an initial
word s. Decompose the entire Hamiltonian by longest common prefix
length ell(w,s). All changes of ell are across the unique pair of
geometry words s[:ell] and s[:ell+1]. The operator on storage/probe
registers across each such block has norm <=(g+epsilon)/sqrt(N).
The total off-block norm is therefore
tau_N<=2(g+epsilon)/sqrt(N); every other term, including large diagonal
energies and ordinary probe motion, belongs to unitary diagonal blocks.

Destroying a prefix r levels below the original end needs at least r
off-block moves. Its probability is bounded by

\[
\left[\sum_{a\ge r}(\tau_N|t|)^a/a!\right]^2,
\]

capped at one. A prefix occupying a fixed fraction of a formed
macroscopic component has r=O(N). The bound tends to zero even for
probe times proportional to N, since tau_N t=O(sqrt(N)). This is
conditional stability of bulk neighborhoods; it is not irreversible
locking of the entire word. It depends on the declared normalization.

## 6. Exact resonant formation, not a favorable selected time

In probe vacuum the mismatch interaction vanishes and M=0 remains
exact. For root fixed, the uniform state of each assembly length spans
an exact N-dimensional invariant radial component. Let m=1,...,N.
After the constant energy shift, its matrix has diagonal zero and
adjacent-length coefficients -g sqrt((N-m)/N).

In reversed basis k=N-m this is the finite oscillator
H_rad=-g(a+a^dagger)/sqrt(N), with truncated ladder k=0,...,N-1 and

`[a,a^dagger]=I-N|N-1><N-1|`.

The seed is |N-1>. In any eigenstate, the expectation of
[H_rad,a-a^dagger] vanishes. The displayed commutator forces every
eigenstate to have seed weight exactly 1/N. The irreducible Jacobi
matrix has simple spectrum. Its time-dephased seed density is therefore
I/N in the radial component. In particular,

\[
\overline p(m\ge L)=\frac{N-L+1}{N},\qquad
\overline p(m=N)=\frac1N.
\]

This proves macroscopic formation, not nonvanishing full-assembly
probability or a stationary single geometry. Full connectivity of all
N events remains unproved; unused events are isolated in that record.

There is also a finite observation window, with no favorable instant
chosen. The scaled eigenvalues are roots of the probabilists' Hermite
polynomial. For y=exp(-x^2/4)He_N(x),
y''+(N+1/2-x^2/4)y=0. Between consecutive zeros, integration and
Poincare's inequality give spacing at least pi/sqrt(N+1/2).
Hence consecutive radial energies have separation at least

\[
\Delta_N=\frac{\pi g}{\sqrt{N(N+1/2)}}.
\]

Let H_(N-1) be the harmonic number and choose
T_N=32 H_(N-1)/Delta_N=O(N log N/g). A time average on [0,T_N]
has off-diagonal spectral terms bounded by 2/(T_N|E_a-E_b|).
Since all seed eigenweights are 1/N, any projector probability differs
from its infinite-time average by at most

`4 H_(N-1)/(T_N Delta_N) = 1/8`.

For L=max(2,ceil(N/2)), the probability of finding a macroscopic path
at one uniformly sampled time in that window is at least **3/8**.
Moreover at least 1/6 of the window has instantaneous success
probability >=1/4 (use probability<=1 to bound the measure of that set).
This is an unconditioned observation-time protocol, not selection of
a rare label ordering or first-passage/continuous-measurement claim.
Formation duration grows with N; no constant-time assembly is claimed.

## 7. Complete energy and recording account

At consistency q-count=m and K=N-m. At resonance every configuration
has binding plus stored energy UN. A net attachment releases U from
the free-event counter and deposits U in its event's storage bit. The
seed's U(N-1) binding energy becomes stored charge; no energy is
created. In the infinite-time ensemble mean m=(N+1)/2, so on average
U(N-1)/2 of that binding energy has transferred to storage.

At detuning, relative to the initial state,

`U(N-1)=U(N-<m>)+(U+delta)(<m>-1)+<H_kin>`.

The finite controls check the complete ledger, including kinetic energy.
A complete geometry record includes every length, every word and root.
For radial populations p_m its outcome entropy is
H(p)+sum p_m log2[(N-1)!/(N-m)!]+log2(N). Its mean energy change is
-<H_kin>, zero at resonance for this seed. No statement that a real
device has zero cost or heat follows. All successful relational records
are retained. The injected probe and its preparation are additional
operations; number conservation does not create it from vacuum.

## 8. Detuning: the resonance is not a robust phase

The perturbation delta sum q_i is local, label symmetric and respects
the declared head/probe conservation rules. It is not forbidden by any
stated S2 principle. In the consistent radial component it adds
delta(m-1). With seed energy zero, residual width
sigma=g sqrt((N-1)/N), and boundary b_L=g sqrt((N-L+1)/N), projection
of the conserved residual onto m>=L yields

\[
p(m\ge L,t)\le
\left[\frac{\sigma+b_L}{|\delta|(L-1)-2g}\right]^2
\]

whenever the denominator is positive. For every fixed nonzero detuning
and L=ceil(N/2), this is an all-time O(1/N^2) upper bound. Exact
resonance, or a shrinking detuning window, is essential to the proved
formation mechanism. A symmetry making this equality unavoidable has
not been derived. The zero-detuning spectral/Schur theorem cannot be
automatically asserted for the extensive detuned band.

Spent charge also does not produce an irreversible lock: assisted pop
has nonzero matrix element g/sqrt(N) and converts the charge back into
free-event energy. The bulk theorem supplies conditional persistence
instead of pretending that reverse paths vanished.

## 9. Naive branching fails the operational dimension test

The separate two-child tree grammar permits leaf deletion and attaching
that free leaf to any bound vertex with spare local child capacity.
Take a tree with two arms of length L. Its two leaves i,j have raw tree
distance 2L. The five moves

`delete j -> attach j to i -> probe i to j -> delete j -> restore j`

restore the **entire** original tree and all charged storage bits while
transferring the probe. Every step remains M=0. Moving i instead gives
a second path. Independent exact full-H evaluation gives unit fifth
coefficient -2 and zero lower coefficients for 2L>5. With amplitudes
g/sqrt(N) and J, that coefficient is -2g^4 J/N^2.

This is primarily an exact word/metric obstruction; the Taylor
coefficient alone is not a uniform finite-time large-N theorem. The
true protected configuration distance between leaves is <=5, and cannot
be identified with their growing raw tree distance. In a complete
binary-tree record with O(N) leaves, a radius-five operational ball at
any leaf contains O(N) events. A uniform finite-dimensional ball-growth
claim therefore fails on this allowed branching class.

No supplied lattice was used in the counterexample. It does not exclude
all multi-defect architectures or prove these trees are typically formed.
It rejects this naive route to dimension and warns that the metric must
be recomputed after changing primitive support.

## Status

At exact resonance S2 supplies explicit energy transfer, appreciable
macroscopic assembly, a re-derived gap and complete conditional virtual
protection with bulk persistence. It does **not** supply a robust
untuned spatial phase, a dynamically selected dimension, universal
CP-12 matter coupling or a distinctive empirical prediction. These
remaining failures are part of the result, not deferred qualifications.
