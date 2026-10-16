# Complete virtual sums, exact dressing and the geometry obstruction

All statements below concern the explicitly defined finite Hamiltonians.
Uniform constants mean independent of vertex number in the declared
bounded-degree family, not uniform over arbitrary graph ensembles or
over the full CP-12 gauge Hilbert space.

## 1. Bare and exact projectors are different

Let P0 project onto the microscopic L register and Q0 onto H for the
static module. P0 is the reference subspace of a detuned bare band; it is
not the exact spectral projector Pi. It has P0 H Q0=-gI, so its virtual
term cannot be discarded. With D=Delta I-vA,

\[
H_F(E)=P_0HP_0+P_0HQ_0(E-Q_0HQ_0)^{-1}Q_0HP_0
=-g^2[(\Delta-E)I-vA]^{-1}.
\]

If one instead chooses Pi, Pi H(1-Pi)=0 identically. That true identity
does not bound transport in the original localized microscopic basis.
An explicit localized dressing is needed, as constructed below.

## 2. All walks and their multiplicity

For an undirected simple graph of maximum degree z, v>0 and
a=Delta-E>zv, the absolutely convergent identity is

\[
\langle j|H_F(E)|i\rangle
=-\frac{g^2}{a}\sum_{m=0}^{\infty}
\left(\frac va\right)^m(A^m)_{ji}.
\]

(A^m)_ji counts **all** length-m walks, including returns and repeats.
For d=d_G(i,j), it is zero for m<d and at most z^m. Therefore

\[
|H_F(E)_{ji}|\le\frac{g^2}{a-zv}
\left(\frac{zv}{a}\right)^d.
\]

The prefactor already pays for all path multiplicity. No order of
perturbation theory is dropped. For adjacent i,j there is also the
strict lower bound |H_F(E)_ji|>=g^2 v/a^2, since every walk term has the
same sign. A gap has not been bought by eliminating all transport.

## 3. Derive the actual protected spectrum

Diagonalize A. For each adjacency eigenvalue alpha the two energies are

\[
f(\alpha)=\frac{\Delta-v\alpha-
\sqrt{(\Delta-v\alpha)^2+4g^2}}2<0,
\quad
h(\alpha)=\frac{\Delta-v\alpha+
\sqrt{(\Delta-v\alpha)^2+4g^2}}2>0.
\]

Write D0=Delta-zv>0. There are exactly N low energies in
[-g^2/D0,0) and N high energies at least D0. The separating spectral
gap is at least D0, independently of N. These assertions follow from
the exact roots, not a supplied projection rule. The same-sign nearest
bound above is uniform throughout the low-energy interval by replacing
a with Delta+g^2/D0.

For the spectral subspace Pi, the maximal principal angle obeys
||Pi-P0||<=g/D0. An initially bare L state has high-register probability
at most 4g^2/D0^2 for every time: in each adjacency mode its transition
amplitude is at most 2g/sqrt(D_alpha^2+4g^2). This is a one-probe result;
no finite-density many-body gap is inferred.

## 4. An energy-independent effective Hamiltonian with local dressing

The exact low-band dynamics are not obtained by exponentiating H_F(0).
For the canonical isometry from the microscopic L space into Pi, set

\[
F=f(A),\quad X=-F/g,\quad
W=\begin{pmatrix}I\\X\end{pmatrix}(I+X^2)^{-1/2}.
\]

Then W^dagger W=I, H_G W=WF, and F is an exact energy-independent
Hamiltonian in an orthonormal dressed basis tied to the original events.

To control its locality, use the weighted row norm

\[
\|M\|_\mu=\sup_i\sum_j|M_{ij}|e^{\mu d_G(i,j)}.
\]

It is submultiplicative by the triangle inequality. Thus
||A||_mu<=z exp(mu), and ||D^-1||_mu<=1/D_mu, where
D_mu=Delta-zv exp(mu). If D_mu>2g, the Catalan expansion of the exact
negative root converges absolutely in this norm:

\[
F=\sum_{k\ge1}(-1)^k C_{k-1}g^{2k}D^{-(2k-1)},
\qquad
\|F\|_\mu\le B_\mu
=\frac{D_\mu-\sqrt{D_\mu^2-4g^2}}2.
\]

Consequently |F_ji|<=B_mu exp(-mu d_G(i,j)) at every N. Let
x=B_mu/g<1. The binomial series gives the separate dressing bounds

\[
\|W_L-I\|_\mu\le(1-x^2)^{-1/2}-1,\qquad
\|W_H\|_\mu\le x(1-x^2)^{-1/2}.
\]

So the basis change has controlled exponential tails; it does not assign
arbitrary global eigenvectors to events. Retaining the first term and
bounding all higher Catalan terms also gives, for every adjacent pair,

\[
|F_{ji}|\ge g^2v/\Delta^2-
\left[B_0-g^2/D_0\right].
\]

The chosen parameters make this lower bound strictly positive.

For Delta=8,v=1,g=1/4,z=4 and mu=log(3/2):

* spectral gap >=4;
* |F_ji|<=0.031754164 exp[-d/2.466304];
* every nearest coupling >=0.00091504579;
* bare high-register probability <=0.015625 at every time.

These constants are analytic bounds evaluated in ordinary floating
point for display. They are not measurements of nature.

## 5. Full microscopic propagation and a physical source

For the static module, distance is fixed independently of the response.
Conjugate the full two-band propagator by exp(mu d_G(source,.)), using
the same weight on L and H. The detuning and conversion commute with
this weight. The Hermitian part of the weighted generator has row sum
at most zv sinh(mu). Gronwall's inequality gives

\[
|\langle j,b'|e^{-itH_G}|i,b\rangle|
\le\exp[zv\sinh(\mu)|t|-\mu d_G(i,j)].
\]

This controls the entire microscopic propagator, including high-band
excursions, and can be combined with its trivial bound 1. It derives
finite propagation from the **supplied** bounded-degree graph; it does
not derive that graph's locality or a relativistic speed.

The numerical source is a bare local L-number phase, with a supplied
initial dressed packet W(|0>+|1>)/sqrt(2). For x=U psi and
y=U n_L0 psi, chi_j=2 Im(x_Lj^* y_Lj). The exact +/-pi/2 finite-pulse
identity checks it independently. The source variance is positive.
This packet is expressly not claimed to arise from CP-12 neutral
formation. Signed chi is a susceptibility, not an arrival probability.

## 6. A large gap still permits a virtual geometry shortcut

For GC1 take P0=(q=0,L) and Q0 its **entire** complement, including
(q=1,L). Leaving those states out would omit an essential virtual path.
For a cycle, both instantaneous graphs have degree at most three.
Let

\[
m=\min(\Delta,\delta)-(3v+g+\lambda)>0,\quad
b^2=g^2+\lambda^2.
\]

The Q0 block is at least mI by its diagonal and absolute row sums.
The full Hamiltonian has N low eigenvalues in [-b^2/m,0] and the rest
at least m. The lower endpoint follows by completing the square in
m||q||^2-2b||p||||q||; the high-band bound and low-band count follow
from interlacing and the zero P0 block. Thus a uniform protected band
really exists in this control.

For any initially bare P0 state its conserved mean energy is zero.
The same quadratic inequality gives, at **every** time,

\[
\Pr(Q_0)\le\frac{4b^2}{m^2+4b^2}.
\]

At the primary parameters m=4.5 this is 0.024096386; in particular
the original q=0 geometry survives with probability at least
97.5903614%. Small geometry excursions do not suffice for locality.

For any E in the entire possible low-energy interval, Q0 H Q0-E is
diagonal positive minus a nonnegative off-diagonal matrix. Its inverse
has a convergent, entrywise positive walk expansion. All Schur paths
therefore have one sign. Four five-edge walks reach the antipodal
receiver and restore q=0, for example

`(0,L,i) -> (0,H,i) -> (1,H,i) -> (1,H,j) -> (0,H,j) -> (0,L,j)`.

Keeping just those four terms gives the all-orders lower bound

\[
|H_F(E)_{ji}|\ge
\frac{g^2\lambda^2v}{(\Delta+\delta-E)^2}
\left[\frac1{\Delta-E}+\frac1{\delta-E}\right]^2.
\]

All omitted terms increase its magnitude. Substituting E=-b^2/m
gives a strictly positive uniform floor for every low energy and every
N with an antipodal chord. At the primary parameters that floor is
**9.4380596377e-7**, while d_cycle(i,j)=floor(N/2) diverges. Hence no
finite size-independent C,xi can satisfy the requested exponential
bound in this original-cycle basis. The numerical Schur solve retains
every matrix element and confirms the floor through N=512. Exact
canonical dressed-band controls through N=64 also retain the distant
coupling; those dressed numbers are numerical controls, not a separate
all-N lower-bound theorem.

This is a full all-orders counterexample to **gap + bounded instantaneous
degree + persistent initial graph implies metric protection**. It is
not an all-orders no-go theorem for CP-12 or for every Path B law.

## 7. Normalizing dense hopping also does not create a ring metric

With A the complete graph and high-band hopping v/(N-1), its row budget
and gap Delta-v remain constant. The exact off-diagonal Schur element is

\[
|H_F(0)_{ji}|=
\frac{g^2v}{(N-1)(\Delta-v)(\Delta+v/(N-1))}.
\]

Its 1/N tail is not exponential in a separately fixed ring distance
of order N. It is local only in the complete-graph metric, whose
diameter is one. A shrinking aggregate or pair amplitude alone does
not provide an extended spatial metric.

## 8. Exact activation in the actual WM1 extension of CP-12

Use the allowed charged ring-plus-chord states from the preceding
packet: chord fields (1,0), a particle at each endpoint, and the
compensating outgoing ring hole at j. One original a_j -> b_chord hop
opens the through-field gate; its reverse closes it. For N>=8 the
antipodal ring distance is at least four. Between the same initial and
final geometry state, the leading L_i -> L_j matrix element is

\[
\langle s,L_j|H_{\rm WM1}^5|s,L_i\rangle=-4g^2v.
\]

There must be two conversions, one mediator chord hop and two original
geometry hops. The four valid orderings each have coefficient -g^2v.
Other five-step words cannot reach the receiver and restore geometry;
diagonal insertions would leave too few steps. The bare-band Feshbach
word with all intermediate states in H is C-G-T-G-C, with coefficient
-g^2v. Exact rational computation checks all 30 permutations at
N=8,12,16,24 and flux 1,3,8 against the original unregulated transition
law, with no electric cutoff.

An independent meet-in-the-middle evaluation computes the complete
`<H^2 target|H|H^2 source>` in exact integer/rational arithmetic. It
retains every original diagonal, rewiring, conversion and high-register
hop and confirms the same fifth coefficient and zero lower coefficients
through order four in all 12 controls. Thus the word selection is also
checked against the full operator, not only against a hand-picked path.

At each fixed N the transfer amplitude begins i g^2 v t^5/30, and its
probability begins g^4 v^2 t^10/900. These are exact Taylor coefficients
on supplied charged states. Delta affects later terms and finite-time
interference. This statement alone supplies neither a uniform
finite-time large-N remainder bound nor an all-orders WM1 resolvent
lower bound. The separate GC1 theorem supplies the stronger all-orders
failure within its own explicitly declared class.

## 9. What a positive coupled construction must bound

For a general proposed Q sector write QHQ-E=D_E-T with D_E diagonal
and D_E>=aI>0. Fix the receiver metric before evaluating transport and
extend it to virtual states by their event coordinate. If the weighted
row norms satisfy ||T||_mu<=t_mu<a, and the two rectangular coupling
norms ||PHQ||_mu<=b_L, ||QHP||_mu<=b_R, the same walk expansion gives

\[
\|PHQ(E-QHQ)^{-1}QHP\|_\mu
\le\frac{b_Lb_R}{a-t_\mu}.
\]

The norm includes **every accessible virtual geometry**, its path
multiplicity and its change of event coordinate. In GC1 the opened
chord contributes v exp(mu d_cycle), so no uniform weighted budget is
available even though unweighted row budgets and gaps are uniform.
Deriving this budget, or a comparably strong cancellation/code theorem,
is the next construction obligation. It cannot be replaced by a claim
that the initial graph is sparse.

This row-norm condition is sufficient, not necessary. In an extensive
many-body system disconnected vacuum processes can make a global norm
grow even when connected transport remains local; a controlled linked
cluster expansion or resummation is another acceptable route. Failure
to satisfy this particular sufficient condition is not a no-go theorem.

Finally an exact spectral weight is conserved: for [Pi,H]=0,
||Pi exp(-itH) psi0||^2=||Pi psi0||^2. A large gap alone cannot increase
the initial preparation's low-energy weight. A changing geometry
observable can still form within that weight, and a subsystem can cool
by exchanging energy with an explicitly included reservoir. Neither
process follows from declaring a gap. No such new formation mechanism
is claimed here.
