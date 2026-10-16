# S5: a history-exact formation computation with explicit omitted-state cost

This is a numerical construction for the unchanged S4 law. The exact
error statements below are in ordinary finite-dimensional quantum
mechanics. Their floating-point evaluation is not outward-rounded
certification or independent mathematical validation.

## 1. The exact formation component

Use the original bondless, uniform-root seed and receiver/probe vacuum.
In probe vacuum the bare interaction is exactly zero. Storage remains
consistent, and k attachments imply k receiver bosons. Thus no reachable
storage inconsistency is omitted in this formation calculation. The
decoupled frozen CP-12 factor provides no energy. Fixing one root changes
only an orthogonal invariant component; average its original weight
afterward, never optimize the root.

For equal-detuning classes a, of non-root sizes n_a and values delta_a,
retain every type word u=(a_1,...,a_k), k_a(u)<=n_a. S4's exact isometry
has word multiplicity product_a (n_a)_(k_a). The relative Hamiltonian is

\[
H=A+B,\qquad
A=\sum_u\Big(\sum_a\delta_a k_a(u)\Big)|u\rangle\langle u|
  +d\Gamma(h),
\]

\[
B=-\frac gN\sum_{u,a}\sqrt{n_a-k_a(u)}
\left(|ua\rangle\langle u|a_0^\dagger+\mathrm{h.c.}\right).
\]

Here h is the actual L_N=64N^2 receiver chain with off-diagonal -1 and
zero relative diagonal. Its physical frequencies are 4+sigma(h), in
[2,6]. The displayed relative diagonal uses the exact conserved charge.
All detunings distinct simply gives the original factorial word space.
No radial, subset-only, Markov, independent-decay or bath-reset reduction
is made.

## 2. A complete moving receiver subspace

Let W:C^r -> C^(L_N) be an orthonormal one-particle frame. Its Fock
lift is an isometry in each charge-allowed sector k<=N-1. The numerical
state retains every word u and every r-mode occupation of total k=|u|.
Its coefficient dimension is

\[
 D_r=\sum_{k=0}^{N-1}W_k\binom{r+k-1}{k},
\]

where W_k counts complete type words at length k. No extra excitation
cutoff appears. The physical receiver can occupy all L_N modes through
the columns of W; r is a numerical active dimension, not a new reservoir.

Evolution by A is exact in this representation: phase each word by its
own energy, and evolve every column W by exp(-ith). The finite chain's
discrete sine eigenbasis gives eigenvalues -2cos(pi j/(L_N+1)) and
contact components sqrt(2/(L_N+1))sin(pi j/(L_N+1)). Every physical
frequency remains in these vectors. This costs no Fock-space truncation.

Before a B step append the orthogonal part of the contact vector,
(1-WW^dagger)|0>. The enlarged span contains contact exactly, and B
preserves its complete Fock space. The implementation may skip a
residual of norm <=1e-13, charging at most g T 1e-13 in total state
norm by Duhamel. Both emission and
absorption are retained. Its exponential can therefore be evaluated
without the physical receiver's enormous occupation enumeration.

## 3. Size-independent commutator error

For an arbitrary one-particle vector v and bounded attachment weights f_a,
let B_f(v) be the corresponding Hermitian creation/annihilation pair.
Unique parents and exact charge give

\[
\|B_f(v)\|\le g\,\max_a|f_a|\,\|v\|.
\]

Indeed an append from k bosons has squared norm at most
(g/N)^2(N-1-k)(k+1)||v||^2; each directed half has norm at most
g||v||/2. This pays for every child and bosonic enhancement.

With d=max|delta_i|<=3/4, ||h|0>||=1 and ||h^2|0>||<=sqrt(2),

\[
\begin{aligned}
\|B\|&\le g,\\
\|[A,B]\|&\le g(d+1),\\
\|[A,[A,B]]\|&\le a_2=g(d^2+2d+\sqrt2),\\
\|[B,[B,A]]\|&\le b_2=2g^2(d+1).
\end{aligned}
\]

The first double commutator follows by applying the word increment and
receiver h twice: the three terms carry delta_a^2, 2delta_a h|0> and
h^2|0>. The second uses ||[B,C]||<=2||B||||C||. Extensive diagonal or
receiver norms do not enter these constants.

For S_2(dt)=exp(-iA dt/2)exp(-iB dt)exp(-iA dt/2), the standard unitary
Strang commutator bound gives

\[
\|S_2(dt)-e^{-iHdt}\|
 \le C(g)\,|dt|^3,\qquad C(g)=a_2/24+b_2/12.
\]

See Childs et al., arXiv:1912.08854v3, Section 5.1, Eq. (150).
For n steps at T=n dt the accumulated bound is C(g)T^3/n^2.
For g=.9 and d=.75, C=.36662675858899... . Choose n before evolution
so this part is <=.04. The step count is O(N^3) at T=4N^2 for fixed
error. This bound is for the physical Hamiltonian, not merely a split
active-space model.

## 4. Receiver projection is charged in state norm

The one-body receiver density has entries gamma_ij=<b_j^dagger b_i>.
Choose a normalized eigenvector v with a small eigenvalue lambda. A
unitary change of active orbitals puts v last, using exact bosonic
two-mode transformations. Let q be the total state probability with
at least one boson in that last orbital. Then q<=lambda, and projecting
onto its vacuum and normalizing changes the complete state by exactly

\[
 e_{\rm drop}=\sqrt{2\big(1-\sqrt{1-q}\big)}
 =\sqrt{\frac{2q}{1+\sqrt{1-q}}}\le\sqrt{2\lambda}.
\]

The code computes q from all coefficient words and occupations, rather
than treating a small eigenvalue as zero. No geometry word is removed.
Selecting a numerical basis from the propagated state is not selecting
a favorable physical word, receiver state or observation time.

At each step the unspent .15 projection budget is apportioned over the
remaining steps. At the fixed rank ceiling a projection may instead
spend any still-unspent part of that same budget. This change uses no
formation outcome. If there is insufficient rank, coefficient capacity,
error budget or wall time, the instance is UNRESOLVED.

## 5. Why errors cannot hide return paths

For normalized numerical states phi_l at successive physical steps,
insert the exact physical unitary on the preceding state. Its norm is
one. The triangle inequality telescopes over any adaptive choice of
receiver frames and projections:

\[
 \|\psi(T)-\widetilde\psi(T)\|
 \le C(g)T^3/n^2+\sum_l e_{{\rm drop},l}
       + e_{\rm exp}+ e_{\rm other}.
\]

No contraction or linearity of the adaptive projection is assumed.
An omitted component's entire future unitary contribution, including
arbitrary geometry/receiver returns, is bounded by its discarded state
norm. This is why compression is legitimate as an approximation but
cannot itself count as physical suppression.

For the finite B exponential, a degree-p Taylor step has norm remainder
at most exp(g|dt|)(g|dt|)^(p+1)/(p+1)!. Normalizing its result costs
the measured |norm-1| as an additional state change. All such terms are
summed. Polar reorthogonalization of a frame with Gram defect eta has
state-change upper bound M eta exp(M eta), M=N-1. It too is charged.

These exact-arithmetic bounds do not certify all floating-point
operations. Report their double-precision evaluation, norm/Gram defects
and independent wavefunction comparisons separately. A successful
saved verifier is not an outward-rounded numerical proof.

The earliest saved trajectories did not print the contact-threshold
term. They are kept verbatim. The final audited ledger adds its
conservative g T 1e-13 upper bound before interpreting any physical
probability. This is at most 1.296e-11 for the primary instances and
changes no gate verdict. New executions print the term directly.

## 6. A formation lower bound and the root average

For the macroscopic event projector G, set p_hat=||G phi||^2. If the
state-error account is e, the reverse triangle inequality gives

\[
 p_{\rm physical}\ge\max\{0,\sqrt{p_{\rm hat}}-e\}^2.
\]

The same inequality applies separately to full assembly or any other
predeclared event. It is stronger than subtracting 2e from a probability.
The complete original root ensemble has the original root-type weights
n_a^global/N; average the lower bounds with these weights. Never
renormalize over the easier roots. An unresolved component contributes
only the trivial lower bound zero, with its missing weight displayed.

Repeated-detuning families use label-permutation covariance exactly.
There is also a useful sign identity. In the receiver chain basis let
P multiply a configuration by (-1)^(k+sum_r r n_r). It changes the sign
of receiver hopping and assisted emission, and preserves diagonal word
energies. Thus P H_delta P=-H_(-delta). The Hamiltonian and seed are
real in that basis. Complex conjugation then proves identical formation
populations under global detuning sign reversal, including finite L_N.
This permits explicitly labeled reuse of mirror root components; it
does not make arbitrary different disorder patterns equivalent.

A finite complete instance also gives a small continuous comparison
region, without pretending its detunings still repeat. After subtracting
the conserved root shift, perturbing each storage energy by at most h
and common g by at most q changes the physical Hamiltonian norm by at
most (N-1)h+q. Duhamel adds T[(N-1)h+q] to the state-error account.
If every root component has margin

\[
 s=\min_{\rm roots}\{\sqrt{p_{\rm hat}}-e-1/\sqrt2\}>0,
\]

then h=s/[4T(N-1)] and q=s/(4T) leave at least half the margin and
guarantee formation probability >=1/2 in that finite neighborhood,
intersected with the frozen perturbation class. It includes genuinely
distinct site detunings. This is a finite-N, small-neighborhood corollary,
not coverage of all width-1/4 patterns or a nonshrinking all-size result.

## 7. Joining a finite formation result to S4

Accept every macroscopic word G, record it and inject S4's adjacent
root-pair neutral probe with the actual conditional field. Do not reset
that field. If p_G>=p_* is earned, then for any Borel ensemble of
accepted records the conditional S4 bounds hold with
r_min=max(1,ceil(N/2)-2). In particular,

\[
 \Pr(G\ \mathrm{and\ storage\ consistent\ at}\ t)
 \ge p_*\,[1-\min\{1,D_{r_{\min}}(t)\}^2].
\]

The source phase is applied after recording, so the formation weights
do not depend on its phase theta. Give rejected trials zero probe
response. The success-weighted relational nearby response per original
trial satisfies

\[
 \chi_{\rm trial}(t)=\sum_{w\in G}p_w\chi_{w_1}^{(w)}(t)
 \ge p_*\,\big[Jt-(e^{2\nu t}-1-2\nu t)\big]
\]

whenever the bracket is positive, with J=1/4 and nu=2.1. At t=.005
the bracket is .00102794836247... . This pays the formation acceptance
probability and accepts every qualifying geometry, not a favorable
ordering or field.

A finite-instance lower bound supplies only a finite-instance joined
protocol. It cannot produce an all-N independent-disorder theorem.
Recording entropy/work and the external neutral probe preparation
remain supplied operations. Original charged matter universality and
dimension selection are not added.

## Research boundary

The method addresses the immediate missing computational ingredient:
exact attachment history with an explicit full-receiver error cost.
It has not proved bounded active rank, polynomial complexity for all
independent detunings, a uniform formation probability at large N, a
real-energy consistency gap, an equilibrium spatial phase or a
distinctive prediction about nature.
