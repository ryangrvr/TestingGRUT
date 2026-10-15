# Exact charged formation with graph choice

These are finite-order identities of the unchanged full CP-12 Hamiltonian,
not a proposal to replace it with a free-fermion law. They calculate the
original neutral ray's actual transition amplitudes. The target is the
label-invariant connected four-regular, fully circulating sector with at
least ceil(N/4) vertex fermions, as defined in 1801afe. For N=5,6,7 that
threshold is two fermions.

## 1. Why the first coefficient is a free-fermion calculation

Let Q=2N and let electric L1 be the sum of absolute half-leg fields. Every
fully circulating configuration with two vertex fermions has L1 at least
2Q+2. Only matter hopping changes this quantity, by one per event. Therefore

\[
P_G H^j\psi_0=0\quad(j<d),\qquad
P_G H^d\psi_0=P_G K^d\psi_0,\quad d=2Q+2.
\]

At this order no diagonal electric/mass/rewiring term can occur. The graph
cannot change during those d events. For a target on a specified graph g,
the only contributing neutral parent is that same g, with amplitude
1/sqrt(D), where D=binom(M,Q). This is a calculation in the specified initial
ray, not an unrestricted spectral norm or a different moment-matched state.

At minimum L1 every occupied edge has fields (1,-1), up to reversal, and
each hole has (2,-1), up to reversal. Since degree four is even, the Gauss
parity argument in 1801afe forces the two holes to have distinct positive
endpoints, and the unit baseline orientation has two incoming/two outgoing
edges at every vertex. Let S be the two positive endpoints. Each vertex of
S has two choices for which outgoing edge carries its hole: four choices
in total. Every connected Euler orientation has every edge on a directed
circuit. The minimum electric-square resource is 2Q+6.

## 2. Determinant identity before Fourier-coefficient extraction

Fourier transform the electric rotors, retaining their quantum amplitudes.
At fixed angles the hopping matrix is

\[
h=\begin{pmatrix}0&B\\B^\dagger&0\end{pmatrix},\qquad C=BB^\dagger.
\]

Starting with all Q active edge fermions filled and all vertex modes empty,
the fermion evolution is a Slater determinant. Its edge-to-edge block is
cos(T sqrt(B^dagger B)); its edge-to-vertex block is
-i sin(T sqrt C) C^(-1/2) B, interpreted by its entire power series.
Replacing k filled edge rows H by vertex rows S gives, up to the known
Fock row-order sign,

\[
\mathcal A_{S,H}(T;\theta)=(-i)^k
\det\cos(T\sqrt C)\,
\det\left[\frac{\tan(T\sqrt C)}{\sqrt C}B\right]_{S,H}.
\]

This follows by factoring the edge-to-edge block and using the determinant
of a row-replaced matrix. Its inverse exists near T=0; the identity is
used as a formal finite-order series, so no global inverse or convergence
claim is needed. Quantum electric-field amplitudes are the Fourier
coefficients of this expression, not averages of its modulus square.

For a prescribed minimal electric target, each hole's unmatched B factor
must be on its positive endpoint. Choosing its negative endpoint would
require two additional hops to reverse that increment. The remaining
edge variables each occur exactly once around the oriented graph.
In C=degree-A, using a diagonal degree term would spend two hops without
creating a required edge variable. Thus at the minimum total order we
replace C by -A. Here A has the directed edge variables, and its transpose
convention does not change a principal minor. Apart from a unit-modulus
Fock/incidence sign, the Taylor amplitude is the coefficient of the
product of all Q variables in

\[
F_S(A)=\det\cosh\sqrt A\,
\det\left[\frac{\tanh\sqrt A}{\sqrt A}\right]_{S,S}.
\]

The two-hole choices at fixed S therefore have identical magnitudes.
All coefficient arithmetic below is rational. Multiplication by d!
recovers an integer H^d matrix-element magnitude.

## 3. Signed cycle/path covers

Define rational weights by

\[
\log\cosh\sqrt x=\sum_{l\ge1}w_l x^l,\quad
\frac{\tanh\sqrt x}{\sqrt x}=\sum_{r\ge0}u_r x^r,
\quad u_r=2(r+1)w_{r+1}.
\]

For example w1=1/2,w2=-1/12,w3=1/45, and
u0=1,u1=-1/3,u2=2/15. The determinant factor is
exp(sum_l w_l Tr A^l). A closed directed trail of l distinct edges has l
cyclic starting points and consequently carries l*w_l. The principal
minor expands into directed paths between the marked vertices, carrying
u_r for a path of r edges and the sign of its endpoint permutation.

A bijection with a small finite enumeration makes the extraction exact:

* At every unmarked vertex pair its two incoming edges with its two outgoing
  edges (2! choices).
* Add a source and sink marker at each marked vertex. Pair the three incoming
  ports with the three outgoing ports (3! choices).
* The resulting permutation of edge/marker ports decomposes into closed
  trails and marked paths, using each ordinary edge once.
* Multiply closed-trail and path weights and the marked-endpoint permutation
  sign; sum over pairings.

The total pairing count for two marks is 6^2*2^(N-2), rather than a full
reachable-state Hilbert space. Repeated visits to a vertex are permitted;
repeated uses of an edge are not. Expanding the exponential supplies exactly
the unordered-cycle combinatorics; the marker permutation supplies the
determinant sign. This is an amplitude sum with cancellations, not a count
of permitted paths.

The implementation checks independent monotone integer half-leg recurrences
on extremal and representative configurations at N=5,6,7, plus degree-two
cycles. The largest integer-DP layer has 22,295 states. All 16 independent
checks agree. The N5 aggregate reproduces the earlier all-960-state result.

## 4. Restore the original neutral normalization and sum all graph choices

Let c_g be the sum over Euler orientations and marked vertex pairs of four
times |F_S|^2 for a fixed labelled graph. A graph-isomorphism class has
N!/|Aut(g)| labelled copies. Its original-ray joint coefficient is

\[
\gamma_{[g]}=\frac{N!/|\mathrm{Aut}(g)|}{\binom M Q}c_g.
\]

Automorphisms reduce the orientation/marked-pair calculation. They do not
substitute an unrestricted quantum symmetry sector: the normalization and
amplitudes above come directly from the actual neutral parent on each graph.
Independent integer checks include graphs with absent edges and the full
fermion sea, so inactive-edge ordering signs are tested in the magnitudes.

For N6 every four-regular simple graph is K6 minus a perfect matching.
There are 15 labelled octahedra among 455 neutral graph parents. The exact
full-target result is

\[
\boxed{p_{G_6}(T)=
\frac{28937809141133489}{173110482737718477580671905802392578125}
T^{52}+O(T^{54}).}
\]

This proves nonzero original-cyclic-subspace access when graph choice exists.
It remains a short-time, extremely small coefficient, not appreciable weight
or a long-time formation result. Reality makes the probability even in T.

For N7 the complement of every four-regular simple graph is two-regular.
Its cycle components must be C7 or C3 disjoint C4; these exhaust the target
host classes. There are respectively 360 and 105 labelled hosts among
116,280 neutral graph parents. Both have 60 Euler orientations, 5,040
minimum-order configurations per host, Q14, cycle rank8, n_v2 and minimum
electric-square resource34. Their per-host formation coefficients differ:

\[
\frac{c_{\overline{C_7}}}{c_{\overline{C_3\sqcup C_4}}}
=\frac{23516510026538557850597839}{39328225756852078421560506}
\simeq0.5979550202.
\]

The rare charged-readout fraction of the first class tends as T->0 to
0.6721452827, versus its host prior 24/31=0.7741935484. The total target
probability is 1.3866817939621972e-26 T^60+O(T^62).

This is an exact topology-sensitive **joint formation kernel** from an
unfiltered neutral input. It is not selection of a dominant unconditioned
geometry, and no spatial dimension can be inferred from seven vertices.
At this first order inter-graph coherence is not used: an incoherent mixture
of the same neutral parents has the same coefficient. The interference
within each fermion transition amplitude is essential.

## Scope

The finite-order identities require no electric cutoff or large-mass
approximation. Initial finite-support states are in every finite power's
domain, which is sufficient for these Taylor orders without asserting an
infinite analytic series. No coefficient is extrapolated to long times.
All earlier all-time vacuum/low-matter bounds remain unchanged.

The determinant method is standard free-fermion mathematics; its CP12
Fourier-coefficient specialization and graph-choice sums are the calculation
here. A primary methodological reference is I. Klich, [Full Counting
Statistics: An elementary derivation of Levitov's formula](https://arxiv.org/abs/cond-mat/0209642).
