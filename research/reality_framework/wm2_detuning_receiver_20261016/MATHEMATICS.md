# S3 derivation: fixed-interval formation, full transport and the lost gap

These are project-specific statements in ordinary quantum mechanics.
They are not an independent validation, an empirical prediction or an
assertion that standard bosonic emission and resolvent methods are new.
The declared perturbation class is delta in [-1/2,1/2], g in [9/10,11/10].

## 1. Exact component and how the S2 inequality is defeated

Geometry transitions still produce rooted words of distinct event labels.
The charge sum q - n_R=1 is exact, so all occupations have at most N-1
quanta. With a probe the component includes all root-preserving partial
words, all non-root storage bits, every receiver occupation with
n_R=sum q-1, and every probe position. The independent register checks
test the complete finite graphs.

Completeness has a direct construction. With the probe at the root,
assisted append followed by bare pop flips a fresh q to spent and adds
one receiver quantum. The inverse cycle removes it. Receiver motion can
bring a quantum to contact. Build arbitrary target storage patterns at
the root, then use bare appends while walking the probe with the head
to build a target word without changing those patterns. End at an
attached probe location by walking; for a free target, append it, hop
onto it and bare-pop it. Every listed occupation distribution is connected
by receiver mode moves. Root storage cannot flip.

In probe vacuum bare moves vanish identically. Thus q equals word
membership and n_R=m-1 exactly. This is a physical invariant component,
not a numerical omission of reachable mismatches.

S2's detuning argument required the macroscopic compression to have a
diagonal separation approximately |delta|m while its kinetic norm was
O(1). Here at attachment count k the receiver relative spectrum spans
[-2k,2k]. The combined delta k plus receiver energy can stay on shell
for every |delta|<=1/2. Its extensive energy range defeats that inverse
bound. The initial residual width actually shrinks as
g sqrt(N-1)/N; no large unexplained initial energy variance is used.

## 2. Exact quadratic reduction

Fix a root, put M=N-1, and let G_k be the normalized uniform superposition
of all ordered k-label extensions of that root, with consistent storage.
Their number is P_k=M!/(M-k)!. The unnormalized append operator A_+
satisfies

\[
A_+G_k=\sqrt{M-k}\,G_{k+1}.
\]

Every child has a unique parent; the adjoint relation holds too. Replace
G_k by a bookkeeping pump occupation |M-k> and retain all receiver Fock
states of total k. On the fixed-M bosonic sector this is an exact isometry,
not an added physical pump or a favorable initial state.

The reduced Hamiltonian is the second quantization of the one-particle
Jacobi operator K:

\[
K_{00}=0,\quad K_{rr}=\delta\ (r\ge1),\quad
K_{01}=K_{10}=-\lambda,\quad
K_{r,r+1}=-1,\quad \lambda=g/N.
\]

The seed is M bosons in the bookkeeping pump. Linear evolution gives a
single-particle amplitude vector z(t)=exp(-itK)|0>. If s(t)=|z_0(t)|^2,
the exact net-attachment distribution is

\[
k\sim {\rm Binomial}(M,1-s(t)).
\]

This includes the entire finite receiver and every stimulated-emission
and reabsorption path. Bosonic enhancement is paid for by the declared
g/N normalization. There is no master equation, continuum-rate fit,
Markov approximation or independent-decay assumption.

The finite small full wavefunction is checked against its multinomial
Fock amplitudes, not just against one averaged observable.

## 3. Uniform continuum comparison without a bound state

First extend the receiver's mode chain to a half-line, only as a solvable
comparison. Its relative one-quantum band is [delta-2,delta+2].
For the convention (z-K)^(-1), the contact Weyl function is

\[
m(z)=\frac{z-\delta-\sqrt{(z-\delta)^2-4}}2,\qquad
G_{00}(z)=\frac1{z-\lambda^2m(z)}.
\]

The square-root branch has m(z)~1/z. Outside the band the denominator
is strictly increasing. At its upper edge it is delta+2-lambda^2>0,
and at its lower edge delta-2+lambda^2<0 throughout our class, since
lambda^2<=0.3025<1/3. There is no isolated pump pole. Inside the band
its imaginary part is nonzero; neither edge is a threshold pole.
The normalized pump spectral density is therefore

\[
\rho(E)=\frac{\lambda^2}{2\pi}
\frac{\sqrt{4-(E-\delta)^2}}
{(1-\lambda^2)E^2+\lambda^2\delta E+\lambda^4}.
\]

The denominator D obeys

\[
D=(1-\lambda^2)
\left(E+\frac{\lambda^2\delta}{2(1-\lambda^2)}\right)^2
+\lambda^4\left[1-\frac{\delta^2}{4(1-\lambda^2)}\right]
\ge\frac9{10}\lambda^4.
\]

Write y=sqrt(4-(E-delta)^2). It vanishes at the edges, is at most two,
and has total variation four. D is convex, so the variation of 1/D is
at most 2/min D. Hence

\[
\|\rho'\|_1
\le \frac{\lambda^2}{2\pi}
\left(\frac4{\min D}+2\frac2{\min D}\right)
\le \frac{40}{9\pi\lambda^2}.
\]

Integrating the Fourier transform by parts gives the exact, uniform
comparison bound

\[
|z_0^\infty(t)|\le
\min\left(1,\frac{40}{9\pi\lambda^2|t|}\right).
\]

For every t>=4N^2, lambda^2 t>=81/25, so

\[
|z_0^\infty(t)|\le A=\frac{1000}{729\pi}<0.437.
\]

This result is uniform over the continuous parameter intervals, not a
claim inferred from the numerical parameter grid. The density calculation
is standard Jacobi/Friedrichs analysis specialized to the displayed law.

## 4. Finite closed receiver and a nonshrinking formation interval

Return to the actual finite L_N=64N^2 model. Both half-line and finite
one-particle operators have geometry-mode hopping norm at most two.
The diagonal delta term can be resummed exactly. A pump-to-pump path
that notices the cut beyond mode L must use at least 2L+2 hops.
Consequently, for 0<=t<=8N^2,

\[
|z_0^L(t)-z_0^\infty(t)|
\le2\sum_{r\ge2L+2}\frac{(2t)^r}{r!}
\le2e^{-120N^2}.
\]

For the last inequality use
sum_(r>=n) x^r/r! <= exp(x)(ex/n)^n,
x<=16N^2 and n>=128N^2; the exponent is smaller than -120N^2.
This is a return-amplitude comparison, not a silent infinite bath.

For every N>=2, every parameter in the frozen class and every
t in [4N^2,8N^2], the finite pump survival is therefore less than 1/4.
Let L_event=max(2,ceil(N/2)). Failure of macroscopic assembly requires
at least N-L_event+1 bookkeeping pump quanta to remain. Markov's
inequality applied to their exact binomial distribution gives

\[
\Pr[m<L_{\rm event}]
\le\frac{(N-1)s(t)}{N-L_{\rm event}+1}
\le2s(t)<\frac12.
\]

Thus the primary finite-model theorem is

\[
\boxed{\Pr[m\ge\max(2,\lceil N/2\rceil)]\ge\frac12}
\]

throughout the entire declared window and parameter class. The ideal
comparison gives the sharper 1-2A^2=0.6186926... lower bound; the finite
error leaves ample margin for the conservative one-half statement.

No favorable observation time, receiver state or label ordering is
selected. All successful macroscopic records are accepted. Full-N
connectivity at this window, an all-time finite-model guarantee and
an irreversible stationary phase do not follow.

## 5. Complete energy, receiver and recording account

The physical receiver spectrum is fixed at [2,6], independent of delta.
Only the combined bookkeeping operator has band [delta-2,delta+2].
The emitted physical frequency Delta-delta lies inside that fixed
band. The law is never retuned to a measured detuning.

Mean net attachments are (N-1)(1-s). Let k denote that expectation.
The exact conserved-energy identity is

\[
16(N-1)=16(N-1-k)+(12+\delta)k
+\langle H_R\rangle+\langle H_{\rm int}\rangle.
\]

Equivalently the relative receiver kinetic energy plus delta k plus
interaction energy is zero. Initial energy, finite capacity N-1 quanta,
positive frequency band, reverse processes and all receiver occupations
are explicit. There is no energy source in the decoupled CP-12 factor.

A geometry record has probability p_k/(N P_k) for each root and word.
Its full outcome entropy is

\[
H_{\rm record}=H(p_k)+\sum_k p_k\log_2 P_k+\log_2 N.
\]

This cost grows with system size. The word measurement commutes with
diagonal and receiver energies and deletes assisted off-word coherence,
so its mean projective work is -<H_int>, bounded in magnitude by g.
Zero or bounded mean projective work is not zero device cost.
Conditional on any recorded word, the receiver is the k-boson Fock
state in the actual evolved emitted mode. No reset is needed for the
probe theorem, which holds for arbitrary initial receiver states anyway.

## 6. Metric and all-state virtual-word barrier

Let s be a recorded word of length m and let i,j be its attached events.
Along any configuration path let k=m-min length be its maximum removal
depth (zero if no initial event is removed). The original prefix of
length m-k is immutable. With h matter hops:

* both endpoints in that prefix require h>=d;
* exactly one removed endpoint gives d<=h+k-1;
* two removed endpoints give d<=k-1,

where d is their original path separation. At least k geometry operations
occur, so the total number of geometry/probe operations is at least d.
The fixed-word route uses exactly d hops. Hence the protected operational
metric equals path separation even when real assembly and receiver
motion are allowed. Resumming receiver-only operations at zero cost
cannot make it smaller.

A P0-to-Q0-to-P0 excursion must contain a pop and an append. An append-only
path cannot repair the first bound/fresh mismatch; a pop-only path
cannot repair the first free/spent mismatch. Receiver motion does not
change those flags. If k>=1, geometry count is at least k+1; if k=0 it
is at least two. The three prefix cases then prove

\[
k_{\rm virtual}\ge d+2
\]

for all final consistent words, of all sizes, with every allowed receiver
state. This is an all-state configuration theorem, not an energy-gap claim.

## 7. Uniform complete full-H propagation

Put the reservoir Hamiltonian and all geometry/probe-preserving
diagonals into D, and the other interactions into V. D need not have a
small global norm; it preserves each word/storage/probe block.

For assisted append, a parent with s stored bits has at most N-s fresh
destinations and at most s creation enhancement. Unique parents give

\[
\|W_+\|^2\le(g/N)^2s(N-s)\le g^2/4.
\]

Thus assisted append plus adjoint has norm <=g. The bare append pair
has norm <=2epsilon, and path matter hopping <=2J. Uniformly in the class,

\[
\|V\|\le\nu=1.1+0.5+0.5=2.1.
\]

All occupations and all possible geometry words are paid for. The
extensive receiver norm is resummed inside D; it is not discarded.

For sources P_(i,d;j) over arbitrary coherent consistent geometry and
receiver packets, and R_j summing all final storage/geometry/receiver
records, the word barrier and interaction-picture Dyson expansion give

\[
\|R_j e^{-itH}P_{i,d;j}\|
\le {\rm tail}_d(\nu|t|)
\le \exp(\nu e^\mu|t|-\mu d),\quad \mu>0.
\]

Here tail_d(x)=min(1,sum_(r>=d)x^r/r!), tail_0=1. Constants are independent
of N, L, delta and g in the frozen class. This is an actual full-unitary
response front; no effective matrix at one energy is exponentiated.

For an injected adjacent probe (|s,i>+|s,j>)/sqrt(2), the source number
variance is 1/4. A pulse exp(-i theta n_i) has response

\[
\chi_\ell(t)={\rm Im}\langle A_j,A_i\rangle,\qquad
A_v=R_\ell U(t)|s,v;\Phi_R\rangle.
\]

Thus |chi_l| is bounded by the product of the two corresponding tails.
The receiver includes all storage and receiver configurations. The
Heisenberg series gives nearby first term Jt and a remainder bounded by
exp(2nu t)-1-2nu t. At t=0.005,

\[
\chi_j(t)\ge0.0010279483...
\]

uniformly for every recorded word of length at least two and arbitrary
receiver packet. Therefore successful macroscopic formation does not
require a favorable conditional reservoir state to retain nearby response.

## 8. Finite probe leakage and conditional bulk persistence

Only the bare interaction couples P0 to Q0; its one-mismatch destination
identifies a unique consistent source. Hence ||QHP||<=epsilon.
Without claiming a gap, Duhamel gives

\[
\|Qe^{-itH}P\|\le\min(1,\epsilon|t|).
\]

This is a finite-time leakage estimate, not S2's all-time mismatch floor.

For a fixed initial word, decompose by longest common prefix. Its
ell-changing edges are the unique pairs of successive initial prefixes.
Across each pair, the receiver creation norm is <=sqrt(N), so the
assisted coefficient is <=g/sqrt(N). Adding the bare term gives a total
prefix-crossing norm at most 2(g_max+epsilon)/sqrt(N)=2.7/sqrt(N).

Destroying a prefix r levels below the original end has probability at
most tail_r(2.7|t|/sqrt(N))^2. For a macroscopic interior prefix and
probe times proportional to N, this tends to zero. The conditional
geometry need not be a fixed entire word and is not irreversibly locked.
The result depends on the declared normalization and recording protocol.

## 9. Complete virtual bound at fixed energy resolution

No small extensive reservoir norm or copied S2 gap is used. For
z=E+i eta, eta>0, take the Q configuration-block distance from the
support of B P_source, where B=QHP. The virtual-word theorem puts
B R_receiver at distance at least d. Clip the distance at d.

D commutes with this weight. In the occupation basis all off-block
coefficients have one sign, and the norm of their absolute kernel is
bounded by nu. Conjugating by exp(mu distance) changes QVQ by at most
nu(exp(mu)-1). Choose mu=log(1+eta/(2nu)); the weighted resolvent norm
is at most 2/eta by a Neumann perturbation of the self-adjoint resolvent.
Using ||B||<=epsilon on both ends yields the complete result

\[
\boxed{
\|R_j^P\Sigma(E+i\eta)P_{i,d;j}\|
\le\frac1{8\eta}
\left(1+\frac{5\eta}{21}\right)^{-d},\quad \eta>0.
}
\]

This includes all repeated reservoir, geometry and storage paths and
sums all consistent receiver records. It is a distance-weighted bound,
not an assertion that a global extensive Hamiltonian norm is small.

Its correlation length diverges as eta tends to zero. Fixed positive
energy resolution cannot be silently turned into S2's real-energy theorem.

## 10. Explicit failure of a uniform real-energy complement gap

This limitation is constructive even for the finite L_N model. Set all
q_i=1, use a word of length N-1, and put the probe on the unique free
event v. There is one free/spent mismatch. With M=N-1 receiver quanta,
the relative diagonal is

\[
80+\delta M+d\Gamma(h_{\rm rel}).
\]

Matter hopping vanishes because v is free. The only bare append goes
back to P and is removed in QHQ; bare pop vanishes. Assisted append is
blocked because all q=1. Assisted pop annihilates a contact quantum.

Put all M bosons into a normalized modulated sine packet phi supported
on the last R=L_N/2 modes. Its contact amplitude is exactly zero.
For N>=64 choose carrier energy e=-delta-80/M, which lies in (-2,2).
The one-particle packet obeys

\[
\|(h_{\rm rel}-e)\phi\|\le2\pi/(R+1).
\]

The normalized Q product state consequently satisfies

\[
\|QHQ\Psi_Q\|\le
\frac{2\pi(N-1)}{32N^2+1}\longrightarrow0.
\]

By the spectral theorem, dist(0,sigma(QHQ)) is at most this value.
An isolated complement gap uniformly around the seed energy is therefore
impossible for this S3 family. The explicit lower inverse-norm bound
grows with N when the inverse exists. A direct million-mode residual
check tests the packet formula without constructing a factorial Q matrix.

These are also **coupled** real-energy poles, not merely dark Q levels.
For this packet B^dagger Psi_Q has norm epsilon/sqrt(N): its unique bare
append returns to the corresponding full word. Let r_N be the displayed
residual upper bound and F the Q spectral projection onto
|E|<=2sqrt(N)r_N. The spectral theorem gives
||(1-F)Psi_Q||<=1/(2sqrt(N)). Since ||B||<=epsilon,

\[
\|B^\dagger F\Psi_Q\|\ge\epsilon/(2\sqrt N)>0.
\]

Thus some Q eigenvalue in that shrinking interval has a nonzero
sandwiched-resolvent residue. Sigma(z) has an actual real pole there.
For N>=64 this interval is already inside [-2.5,2.5]. A uniform
bounded real-energy Schur operator throughout that window is false for
the declared P0. Redefining P to absorb these levels would require a
new operational and spectral analysis. Other energy-defined sectors,
including possible ground-state sectors, are not ruled out by this test.

This is a spectral obstruction. It does NOT prove an off-metric matter
shortcut or that every particular off-distance block diverges; the
full-H prefix response bound still holds.

## Status

S3 earns a restricted, nonshrinking interval for finite-model macroscopic
formation, full-H operational propagation, nearby response and conditional
bulk persistence, with explicit scalable receiver resources. It does
not retain S2's isolated real-energy low-band protection. Generic local
perturbations, selected dimension, original-matter universality,
equilibrium spacetime and empirical predictions remain unearned.
