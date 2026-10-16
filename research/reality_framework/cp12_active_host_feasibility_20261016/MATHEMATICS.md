# Full active hosts: forced geometry, recording and persistence

All statements concern the frozen finite-N CP-12 law, with
M=N(N-1)/2 possible edges and fixed Q active edges. They do not select
these inputs or constitute distinctive predictions about nature.

## 1. A dense-host obstruction at the currently accessible primary sizes

If two vertices have distance at least three, they are nonadjacent and
have no common neighbor. Their degrees sum to at most N-2. The remaining
N-2 vertices contribute at most binomial(N-2,2) edges. Hence

`Q <= binomial(N-2,2)+(N-2) = binomial(N-1,2)`.

A disconnected graph also has at most binomial(N-1,2) edges: put all
possible edges into a component of N-1 vertices. Consequently, if
`Q>binomial(N-1,2)`, every host is connected and has diameter at most two.

The primary N=5,Q=10 graph is K5, diameter one. At N=6,Q=12 all 455
labelled hosts are connected, diameter two. Every such host has all-root
mean balls `[1,5,6]`, and mean pair distance `2-Q/M=1.2`. These statements
hold at every time under every physical state in those sectors; they are
fixed-Q facts, not dynamical selection. They do not exclude larger-N
phases, where that edge-count condition no longer holds.

## 2. Complete host record: no circulation filter

Let w be the evolved original coherent neutral ray and let P_g be each
labelled active-host projector. Record **all** outcomes. Hopping and the
physical diagonal preserve g; only raw off-diagonal rewiring R_off
connects different hosts. Therefore

`Delta E_record = sum_g <P_g w,H P_g w> - <w,Hw> = -<w,R_off w>`.

Writing P=(-1)^nv, hopping flips P while rewiring and all diagonals
preserve it. Thus `K=(H-PHP)/2`. If D_phys is the non-rewiring diagonal
and D_R the full escape diagonal, then

`Delta E_record = <D_phys+D_R> - (<H>+<PHP>)/2`.

This evaluates raw rewiring even when a signed permutation quotient
turns a physical graph change into an orbit diagonal. For an isomorphism
class with probability p_a and L_a labelled members, the entropy is
`-sum_a p_a log2(p_a/L_a)`. Class, labelled-record and fixed-source
peripheral probabilities are reported separately.

On the initial neutral ray the work is exactly
`Q(M-Q)/(10M)` and the entropy is `log2 binomial(M,Q)`. The record has a
nonzero energy effect; its entropy is not a claim about dissipated heat
or a derived measuring device. At Q=alpha N and fixed alpha the initial
work is asymptotic to alpha N/10, while the record entropy grows as
alpha N log2 N+O(N). A classical host pointer is not earned by using it
as a diagnostic intervention.

## 3. Exact cross-host norm identity

For one specified host let `B_g=(1-P_g) H P_g`. A removable active edge
has zero fields on both half-legs and a filled edge mode. Rewiring
replaces it by an absent neutral edge with amplitude `-1/(10M)`.
For a fixed initial g, each removed/added pair has a different destination
host. Rewiring preserves the complete fermion and electric configuration;
different source basis configurations therefore also have orthogonal
destinations. If Pi_e is the neutral/removable projector,

\[
B_g^*B_g=\frac{M-Q}{100M^2}\sum_{e\in g}\Pi_e,
\qquad
\|B_g\|=c_{N,Q}=\frac{\sqrt{Q(M-Q)}}{10M}.
\]

The norm is attained by a neutral configuration, which satisfies the
Gauss constraints. For every initial state entirely in that host,
including any charged state obtained by a host-only record,

\[
\|(1-P_g)e^{-iHt}\phi_g\|\le c_{N,Q}|t|,
\qquad
\langle P_g\rangle_t\ge \max(0,1-c_{N,Q}^2t^2).
\]

Use Duhamel's formula, or integrate the two-block probability-flow
bound, to obtain the inequality. It includes all matter dynamics and
all repeated graph returns. It is not a stochastic rewiring rate.

At N=6,Q=7, lag 1/2, the universal lower retention bound is 0.99937778.
This holds for **any** supplied host. At Q=alpha N,
`c^2 ~ alpha/(50N)`, so all recorded hosts become nearly static at fixed
time as N grows. That is a direct consequence of the supplied rewiring
normalization `1/(10M)`, not evidence that spatial neighborhoods have
been dynamically selected or protected.

## 4. Conditional transport relative to a frozen host

Let H_g=P_g H P_g, with its original escape diagonal retained. The
embedded frozen-host trajectory is a comparison used only for a bound:

`||exp(-iHt)phi_g-exp(-iH_g t)phi_g|| <= c |t|`.

The same applies to n0 phi_g, whose norm is at most one. Consequently,
for the density-phase response,

`|chi_full-chi_frozen_host| <= 4 c |t|+2(c t)^2`.

Actual source runs always evolve full H. No deleted-rewiring comparison
counts as a geometry-formation success. Within a fixed bounded-degree
host, microscopic hopping is local on its vertex/edge incidence graph;
this locality is an input of H. The graph-rewiring correction above is
uniform in the prepared state and vanishes at fixed time for Q=O(N).
It does not provide a useful sharp front at the small sizes here, select
a low-dimensional host or validate a relativistic continuum.

These single-host bounds must **not** be applied to the unrecorded
coherent superposition of many hosts. Off-diagonal coherence changes
both its dynamics and the energy cost of recording. Original-neutral
host marginals are computed independently under the full H.

## 5. Uniform fixed-Q hosts do not supply a spatial phase

For the neutral uniform labelled ensemble define I as the number of
isolated vertices. Exact moments are

\[
p_1=\frac{\binom{M-(N-1)}Q}{\binom M Q},\quad
p_2=\frac{\binom{M-(2N-3)}Q}{\binom M Q},
\]
\[
\mathbb E I=Np_1,\quad
\operatorname{Var}I=Np_1+N(N-1)p_2-N^2p_1^2.
\]

For Q=alpha N, p1 tends to exp(-2alpha), and p2/p1^2=1+O(1/N).
Thus Var(I)=O(N), while E(I) is proportional to N. Chebyshev gives
`P(connected) <= P(I=0) <= Var(I)/E(I)^2 = O(1/N)`.
The exact finite-size moment bound and seeded graph samples are reported.
This is a fact about the neutral **combinatorial baseline**, not a
theorem that full CP-12 dynamics cannot bias its host distribution.

Average degree four is not dimension. On the same N=64, a ring with
next-nearest neighbors and a periodic square graph are both connected
four-regular Q=2N hosts. Their radius-1,2,3 balls are `[5,9,13]` and
`[5,13,25]`. These are supplied counterexamples, not formed geometries.

## 6. Error and claim boundaries

Neutral trajectories use complete electric-boundary and Krylov relation
residuals inherited unchanged from the published parent. Source branches
use signed rooted lifts, with no matter/electric filter. Numerical tail
trimming contributes its discarded norm to the conditional state error;
all unretained hopping and rewiring rows remain in the residual action.
Receiver responses use the bilinear norm intervals and exact finite-phase
identity from the parent. These calculations use ordinary floating point
and stated allowances, not outward-rounded formal certificates.

One final numerical allocation closes every remaining physical row reached
from the initial numerical source support, within the original orbit cap.
It gives exactly zero initial outgoing residual on both approximate source
vectors. The additional Galerkin evolution error then starts at second
order in time for those finite-support vectors. The conditional preparation
error remains, as do later outgoing residuals. This is numerical basis
closure, not physical suppression of rewiring or an exact full trajectory.

None of these identities fixes a coupling or quotient value from GRUT's
earned structure alone. They are consequences of the explicitly supplied
CP-12 Hamiltonian and preparation. The frozen SCOUT-1 conclusions are
preserved, not reopened through a change of terminology.
