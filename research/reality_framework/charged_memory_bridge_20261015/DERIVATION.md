# A sourced finite-time charged-sector memory bridge

Parent: `3538c9a44afa34a7d8e5ce94d620b00328c3b2df`. The frozen CP-12
Hamiltonian, its parameters and original neutral preparation are unchanged.
This is a projection of unitary dynamics, not a new fundamental memory law.

## 1. The source and return kernel

Let P be a gauge-invariant charged-network projector, diagonal in the physical
graph/occupation/electric basis. Its definition need not specify degree,
dimension, graph labels or a spatial lattice. Write Z=I-P, and set

\[
A=PHP,\qquad C=PHZ,\qquad D=ZHZ,\qquad
x(T)=P\psi(T),\quad z(T)=Z\psi(T).
\]

Here Z is a projector, not the fixed edge count Q. The electric, mass and
rewiring-diagonal operators commute with P. The finite-N hopping and rewiring
terms are bounded, so C is bounded. Removing the off-diagonal blocks is a
bounded self-adjoint perturbation of H; the diagonal blocks A,D therefore
define self-adjoint generators on their respective domains. This justifies
the unitary propagators below despite the unbounded integer electric fields.

For the original neutral ray, x(0)=0 and z(0)=psi_0. Solving the complement
equation exactly and inserting it into the charged equation gives

\[
\boxed{\dot x(T)=-iAx(T)-if(T)
-\int_0^T K(T-s)x(s)\,ds,}
\]
\[
\boxed{f(T)=Ce^{-iDT}\psi_0,\qquad
K(T)=Ce^{-iDT}C^\dagger.}
\]

The inhomogeneous source is essential: deleting it with x(0)=0 would predict
no charged formation at any time. The kernel retains coherent excursions out
of, and back into, the charged sector. It need not be positive, short-lived,
time-independent in a rate representation, or interchangeable with a thermal
transition matrix. The target probability is ||x(T)||^2; these equations do
not close on that probability alone.

This is the standard time-domain Feshbach/projection construction, specialized
here to the CP-12 charged event. Related established methods include
[Chruscinski and Kossakowski (2013)](https://arxiv.org/abs/1302.6218) and
[Kelly et al. (2016)](https://arxiv.org/abs/1603.01905). No novelty is claimed
for projection methods themselves. Computational savings must be demonstrated,
rather than inferred from writing a formal memory equation.

## 2. A low-rank interface for the Volterra history

Variation of constants gives a second-kind Volterra equation

\[
x(T)=g(T)-\int_0^T L(T-s)x(s)\,ds,
\]
\[
g(T)=-i\int_0^T e^{-iA(T-s)}Ce^{-iDs}\psi_0\,ds,
\quad
L(T)=\int_0^T e^{-iA(T-s)}Ce^{-iDs}C^\dagger\,ds.
\]

If C=U Sigma V^dagger has rank r, with U,V isometries, define

\[
J(T)=\int_0^T e^{-iA(T-s)}Ce^{-iDs}V\Sigma\,ds,
\quad y(T)=U^\dagger x(T).
\]

Then L(T)=J(T)U^dagger, and only r history variables are needed:

\[
\boxed{y(T)=U^\dagger g(T)
-\int_0^T U^\dagger J(T-s)y(s)\,ds.}
\]

Recover the full charged amplitude from
x(T)=g(T)-integral J(T-s)y(s) ds. This preserves the source and returns,
without fitting a population rate or dephasing the state. In finite spectral
coordinates the required integral is evaluated stably as

\[
\int_0^T e^{-i\alpha(T-s)}e^{-i\delta s}\,ds
=T e^{-i(\alpha+\delta)T/2}
\operatorname{sinc}[(\delta-\alpha)T/2],
\]

where sinc(u)=sin(u)/u. `memory.py` handles degenerate eigenvalues with this
identity, and uses composite trapezoid on the causal integral. The integrated
kernel vanishes at T=0, as does y(0), so its grid recurrence is explicit.

For the N4,Q4 invariant benchmark with electric-square budget 256, the
86-by-491 coupling matrix has numerical rank seven, with factorization
spectral residual below 7e-16. At budget 128 its rank is six. Seven is thus
a property of this regulated control, not a universal number of physical
degrees of freedom or a claim about the unbounded-rank limit.

Full complement spectral preprocessing is still used in this validation.
The history dimension is reduced; a large-N complexity improvement has not
been proved. The reservoir construction below supplies an alternative to
retaining the whole complement during time evolution.

## 3. A general first-return sign and an exact CP-12 check

Suppose x(T)=a T^d+O(T^(d+1)) with a nonzero. The first-transfer amplitude
g(T) has the same leading term. Since L(T)=CC^dagger T+O(T^2), the exact
Volterra identity gives

\[
x(T)-g(T)=-\frac{CC^\dagger a}{(d+1)(d+2)}T^{d+2}
+O(T^{d+3}).
\]

Consequently the first return correction to the charged probability is

\[
\boxed{\|x(T)\|^2-\|g(T)\|^2
=-\frac{2\|C^\dagger a\|^2}{(d+1)(d+2)}T^{2d+2}
+O(T^{2d+3}).}
\]

Its first coefficient is nonpositive. This is a local-in-time consequence
of the exact memory equation, not a monotonicity theorem or a large-system
formation bound. For real H, real input and real P the probabilities are even,
so the next possible power in their difference is 2d+4.

The validation event is all four active edges on a directed spanning electric
cycle, with n_v>=1. At N4,Q4 the original ray contains all 15 graph parents;
none is selected initially. Electric L1>=2*4+1 implies a minimum of nine
hops. The exact full unregulated integer recurrence finds 24 nonzero minimum
targets, verifies first-transfer/full equality through amplitude order ten,
and first detects their difference at order eleven. It gives

\[
\boxed{p(T)=\frac{1682}{40186125}T^{18}
-\frac{10368151}{128595600000}T^{20}+O(T^{22}),}
\]
\[
\boxed{p(T)-\|g(T)\|^2
=-\frac{3364}{736745625}T^{20}+O(T^{22}).}
\]

Independently applying C^dagger to the exact leading vector gives
||C^dagger a||^2/||a||^2=6, reproducing the general return coefficient.
The integer recurrence has no electric cutoff and does not depend on the
floating-point invariant-sector reduction. These series are short-time
boundary data, not estimates at T=5 or T=10.

## 4. Reservoir reduction with a full-law acceptance criterion

The complement Krylov space is seeded by the original input and every
coupling port: span{psi_0, C^dagger, D psi_0, D C^dagger, ...}. It preserves
the physical source and does not discard the charged sector. Let E be the
resulting isometric embedding of the entire retained charged space plus the
compressed complement. Set h=E^dagger H E and F=HE-Eh. Duhamel gives

\[
\boxed{\|e^{-iHT}\psi_0-Ee^{-ihT}v_0\|
\le \|\psi_0-Ev_0\|+
\int_0^T\|F e^{-ihs}v_0\|\,ds.}
\]

The electric-budget error must also be included. If H_B is the complete
retained physical-budget Hamiltonian and B_out holds **all** outgoing
unregulated transitions, the residual is represented by the orthogonal stack

\[
F=\begin{bmatrix}H_BE-Eh\\B_{\rm out}E\end{bmatrix}.
\]

No rotor state is silently declared nonexistent. The inherited `Evolution`
routine evaluates an analytic Duhamel upper sum, expanding the leakage vector
on midpoint intervals through fifth order with a sixth-order remainder. Its
evaluation and the isometry checks use ordinary floating point. A separate
5e-10 numerical reporting allowance is retained; these are not outward-rounded
interval certificates. An exact-isometry statement underlies the analytical
bound; observed numerical isometry defects are about 1e-14.

If the resulting norm bound is epsilon, contraction of any physical event
projector gives

\[
\max(0,\sqrt{p_{\rm red}}-\epsilon)^2
\le p_{\rm full}\le
\min(1,(\sqrt{p_{\rm red}}+\epsilon)^2).
\]

The same unitary-dilation idea also bounds a proposed approximate memory
solution: reconstruct its complement using the exact complement propagator.
If its charged equation residual is r(T), the full-state error is at most
the integral of ||r(T)||, plus initial error. Replacing f or K requires
including their errors in that residual. A memory cutoff, Gibbs weight or
Markov rate is acceptable only after this error test, not by analogy.
