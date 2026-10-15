# Source responses and complete relational records

All equations use the frozen CP-12 law at parent `4e7d6c5`. Numerical
projections below approximate this law; they introduce no transport rule.

## 1. Preparation and a core pointer

Let `w=exp(-i H t_w) psi_neutral`, where the neutral ray is the original
equal coherent sum over labelled fixed-Q active hosts. In a physical
occupation/electric basis state s, let C(s) be the full edge mask of
directed circulation, using the inherited corrected through-edge test
`E_i E_j<0`, including one-leg hole edges such as `(2,-1)`.

For each spanning simple undirected cycle C, P_C projects onto states
with **C(s)=C**, not merely onto states containing C. It filters neither
outside electric states nor matter occupations. The complete instrument
is `{P_C for all C, 1-sum_C P_C}`. A successful branch is

\[
\phi_C=P_Cw/\sqrt{p_C},\qquad p_C=\|P_Cw\|^2.
\]

This is a diagnostic record followed by conditional evolution. It does
not derive autonomous observers or classical geometry. A ring is an
analyst-specified candidate event; its shape is not a discovered spatial
phase. The pointer C is held fixed when defining receiver distances.

In the auxiliary Q=N+1 sequence, a successful host has one extra active
chord f. Refine the instrument to `{P_(C,g), P_fail}` if the initial host
is also recorded. This includes **every** successful host, with no
favourable electric or matter-state preparation. The endpoint control
uses a class with chord distance r and source 0 at one endpoint.

Full permutation symmetry makes all records in a class equiprobable.
There are `(N-1)!/2` undirected spanning cycles. The core-host record
counts are

\[
M_r=\frac{(N-1)!}{2}
\begin{cases}N,&2r<N,\\N/2,&2r=N.\end{cases}
\]

One labelled outcome has probability p_r/M_r, while the probability of
the class with a **fixed** source vertex incident to its chord is
`2 p_r/N`. Relabelling a computed branch is a use of exact symmetry,
not free access to a rare specified outcome. All probabilities are
reported separately.

## 2. A nontrivial gauge-invariant source

Use A=n_0, a gauge-invariant projector, and pulse `exp(-i theta A)`.
For a branch phi, define

\[
x(t)=e^{-iHt}\phi,\quad y(t)=e^{-iHt}A\phi,\quad
\chi_O(t)=2\operatorname{Im}\langle x(t),O y(t)\rangle.
\]

Then `chi_O` is exactly the derivative at theta=0 of the expectation of
O. The source variance is `mu(1-mu)`, `mu=||A phi||^2`. A positive lower
variance bound excludes the number-eigenstate/global-phase failure.
For a projector source the finite-strength check is also exact:

\[
\chi_O(t)=\frac{\langle O\rangle_{+\pi/2}
                         -\langle O\rangle_{-\pi/2}}2.
\]

Each receiver is O=n_j. A reported shell is its per-vertex mean over
fixed initial distances d_C(0,j). Each individual receiver is evaluated
directly, or averaged only over its exact stabilizer orbit: vertices
in that orbit have identical expectations under H, source and initial
branch. Shells need not be stabilizer orbits after recording a host.

For Q=N+1, recording a chord with endpoints 0 and floor(N/2) separates
an actual bypass from dilution over the possible source locations. At
N=5, receivers 2 and 3 both have core distance 2, but only 2 has a
direct initial source edge. At N=6, receiver 3 has core distance 3 and
initial host distance 1, while receivers 2 and 4 have core distance 2
and initial host distance 2. Every subsequent hopping and rewiring
amplitude is retained.

## 3. Survival without a hidden final selection

S_C is the diagonal projector that the same original C remains a
directed ring; it permits additional circulating edges. Report its
probability, the raw density response, and the joint response of S_C O.
If a final conditional density is also requested, its derivative is

\[
\chi_{O\mid S}=
\frac{\chi_{SO}}{p_S}
-\frac{\langle SO\rangle\chi_S}{p_S^2}.
\]

The second term is essential. S_C is a global observable, so responses
conditioned on S_C can include selection effects. Raw local receiver
responses, not joint-survival signals alone, are the locality test.
Survival at sampled lags is not a continuously monitored first-passage
survival probability; no intermediate measurement or Zeno evolution
is inserted.

## 4. Complete record work and entropy

For a complete instrument `{P_a}`, the unconditional record changes
the energy by

\[
\Delta E=\sum_a\langle P_aw,H P_aw\rangle-\langle w,Hw\rangle.
\]

Different exact spanning-cycle masks have no one-step H coupling.
Rewiring preserves all electric fields and hence C. A hop can preserve
a directed arc, add an arc, or remove an arc; it cannot reverse an
integer through arc in one step. Adding an arc cannot destroy an
existing cycle, and removing one cannot create a different cycle.
Thus one hop cannot replace one exact spanning cycle by another.
Writing P=sum_C P_C gives

`Delta E_core = -2 Re <P w,H(1-P)w>`.

Refining by g additionally removes rewiring coherence between hosts
inside successful C records. Hopping preserves g. With R_off denoting
the raw off-diagonal rewiring transitions,

`Delta E_(core,host) = Delta E_core - <P w,R_off P w>`.

Raw transitions must be used here: a permutation quotient can turn a
physical rewiring into an orbit diagonal. Inspecting only changes of
the canonical host label would incorrectly lose that contribution.
A full labelled N=4 calculation independently checks both identities.

For core-only outcomes the entropy is `H2(p)+p log2(M_C)`. For refined
records it is

\[
S_{\rm record}=-(1-p)\log_2(1-p)-\sum_r p_r\log_2(p_r/M_r).
\]

These are the energy change and Shannon entropy of the specified
instrument. They are not a microscopic measurement-device model, a
heat-dissipation theorem or a free-preparation claim. The full labelled
control also measures how the nonselective record changes a later
fixed-receiver response, retaining all failure outcomes.

## 5. Numerical errors

For a retained basis Q with `A=Q^* H Q`, include **all** rows of
`B=(1-QQ^*) H Q`. Physical escape diagonals stay unchanged. The exact
Duhamel bound is

\[
\|e^{-iHt}v-Qe^{-iAt}v\|
\leq\int_0^t\|Be^{-iAs}v\|\,ds.
\]

The inherited midpoint-Taylor integral upper sums keep terms through
order five, a sixth-derivative remainder and a separate lower-order
bound. The smaller of the two valid upper bounds is retained.
Adaptive extension selects numerical outgoing rows, never physical
histories; every unselected and newly outgoing row remains in B.

For neutral six-site preparation, a real Krylov basis V reduces the
complete finite electric-budget matrix H_b. It uses two-pass full
reorthogonalization and a Hermitian `T=V^T H_b V`. Explicitly evaluate

`H_b V - V T = d e_last^T + E`.

The original full-law residual is bounded by the sum of
`||B V c(t)||`, `||d|| |c_last(t)|`, and `||E||_F ||c(t)||`.
Thin QR of B V replaces it by R with `||B V c||=||R c||`.
Integrate the first two actions with the same analytic upper sums and
add `t ||E||_F`. Hence a small Krylov matrix does not ignore the charged
electric boundary or assume irreversible memory loss. Smaller-system
full-basis comparisons test this implementation, including deliberately
incomplete Krylov bases whose measured errors are not negligible.

If the original neutral approximation has error epsilon, exact
permutation symmetry and disjoint equivalent records imply an error
at most `epsilon/sqrt(M)` in each unnormalized recorded branch. Add
the norm of discarded initial numerical amplitudes, which is not a
physical filter. If the approximate branch norm is a and rho=e/a<1,
normalization with a free common global phase obeys

\[
\delta=\sqrt{2-2\sqrt{1-\rho^2}}
=\frac{\sqrt2\rho}{\sqrt{1+\sqrt{1-\rho^2}}}.
\]

Both phi and A phi share this phase. Evolving their approximations
with norm bounds e_x,e_y gives the receiver interval radius

\[
2\big[e_x\|O\hat y\|+e_y\|O\hat x\|+e_xe_y\|O\|\big].
\]

Positive diagonal probabilities use the square-root norm bounds;
ratios and final conditional derivatives use interval arithmetic.
Work uses `||K|| <= sqrt(2Q min(N,Q))` and
`||R_off|| <= Q(M-Q)/(10M)`. Entropy bounds optimize the scalar terms
over their probability intervals. All calculations use ordinary
floating point with stated reporting allowances, not outward-rounded
interval arithmetic or a machine-checked proof.

## 6. Distances and what these tests can establish

The intrinsic ball observable is `B_C(i,r)=#{j:d_C(i,j)<=r}`. A ring
has `min(N,2r+1)` vertices in such balls; this is a property of the
chosen ring event, not a derived dimension. The separate unfiltered
scan reports ball growth of the **entire** circulating mask in the
original neutral evolution, conditioned explicitly on spanning
connectivity. No coordination number is equated with dimension.

A decreasing shell mean or outside/core action ratio does not bound
any specified distant receiver. Conversely, nonmonotone finite-system
susceptibilities alone are not a general no-go theorem for locality:
interference and initial correlations can create nodes in local
systems. The endpoint controls test a concrete candidate metric with
an active bypass, not all possible collective metrics or phases of
CP-12. Only growing-distance early-time envelopes and robust formation
scaling could support the proposed spatial-phase milestone.
