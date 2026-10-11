# CP-04 follow-up: geometry, symmetry and scale-response separation

2026-10-11. Baseline CP-04 f4c4de4912453d1e5e47ec2310aeef597d445312. This is follow-up interpretation of the frozen construction, not a replacement action, rerun or CP-05. The user supplied successive reviews, including a reported Claude check. That reported check is retained as user-provided review evidence, not authenticated independent scientific clearance or evaluation-kit approval.

## Verified source passages

Full primary source: Kaloper–Padilla, [arXiv:1406.0711v3](https://arxiv.org/pdf/1406.0711), 39 PDF pages. Printed p.23, PDF index 23 (24th PDF page), §4.3, Eq.(47), contains the FRW integration-by-parts equation. Exact short quotation:

> a large old approximately FRW universe does require k = 1

The surrounding paragraph distinguishes singular endpoints from a curvature cutoff, discusses sensitivity of small universes to omitted boundary contributions, and limits the large-old-cosmology inference accordingly. This is verified prior literature, not a GRUT novelty claim.

Printed p.27, PDF index 27 (28th PDF page), §6, following Eq.(55), says:

> Having an exponential form of σ reduces this sensitivity to logarithms

Its context is cosmological-history sensitivity and hierarchy accommodation, not a claim that every fixed-input derivative vanishes. Printed pp.26–27 discuss affine sigma and the external mu scale. These citations were inspected directly; the user's reviewer reported inability to access them, not evidence against their existence.

## Exact FRW statement and its boundary

For a>0 on an interval, finite nonzero volume, integrable curvature and well-defined endpoint limits,

\[
\int a^3R\,dt=6\left([a^3H]-\int a^3H^2dt+k\int a\,dt\right).
\]

If the leading CP-04 constraint gives average R=0 and [a^3H]<=0, then a genuinely nonstatic history has

\[
k\int a\,dt=\int a^3H^2dt-[a^3H]>0,
\]

so k>0. Compact flat T^3 changes the spatial volume factor, not this sign. The endpoint sign must be part of the theorem. **Finite four-volume alone does not prove bang-to-crunch endpoints**: CP-04's own finite flat slab is a counterexample to that implication. Nor does integrability automatically survive a singular kinetic-dominated crunch. For a~s^(1/3) near a crunch, s=t_crunch-t, a^3R~ -2/(3s), so the integral diverges logarithmically; its cutoff/global matching cannot be silently discarded. No universal claim about all UV completions follows.

Dust control (M=1, positive density, negative effective cosmological constant):

\[
a=\sin^{2/3}t,\quad H=\tfrac23\cot t,\quad
\rho=\tfrac43\csc^2t,\quad\Lambda_{eff}=-\tfrac43,
\qquad 0<t<\pi.
\]

Here integral a^3 dt=pi/2, integral a^3 rho dt=4pi/3, and integral a^3 H^2 dt=2pi/9. Consequently average rho=8/3, average H^2=4/9, average R=-8/3. This is an exact recollapsing GR solution, **not** a full CP-04 solution: its Lambda_eff=-average rho/2 differs from the global dust requirement -average rho/4. No numerical cosmological solver was run for this follow-up.

## Full-action shift and what sigma actually selects

Under the specified fixed-lambda transformation v->v+delta v, Lambda->Lambda-lambda^4 delta v, the bulk action is invariant while z->z-d, d=delta v/mu^4. The frozen cubic function gives

\[
\sigma(z-d)-\sigma(z)=-d(1+z^2)+zd^2-d^3/3.
\]

For a differentiable function on the full real line, requiring this difference to be independent of z for every d implies sigma'(z-d)=sigma'(z) for every z,d; sigma' is constant, so sigma is affine. C2 regularity then permits writing sigma''=0. This is only a theorem for the stated transformation/domain, not all possible sequestering symmetries.

On the flat-vacuum branch z=-v/mu^4,

\[
\lambda^4=\frac{\sigma'(z)}{\mu^4\Omega},\qquad
\left.\frac{\partial\ln\lambda}{\partial\ln\Omega}\right|_{v,\mu}=-\frac14,
\qquad S_v=\left.\frac{\partial\ln\lambda}{\partial\ln|v|}\right|_{\Omega,\mu}
=\frac{z\sigma''}{4\sigma'}.
\]

Both affine and nonlinear choices have the same explicit fixed-v volume sensitivity. Affine sigma has zero v-sensitivity; the frozen cubic has v^2/[2(mu^8+v^2)]. The added v-sensitivity belongs to the selected function, not to every sequestration model. This does not imply in-universe mass drift from ordinary renormalization bookkeeping.

## Hierarchy ledger, with no free z knob

Along a sign-fixed nonzero z interval where sigma'>0,

\[
\int S_v\,d\ln|z|=\frac14\ln\frac{\sigma'(z_1)}{\sigma'(z_0)}.
\]

If sigma' grows from 1 to 10^240, the integral is 60 ln10=138.15510557964276. This is an integrated budget, not a theorem requiring large endpoint sensitivity. A plateau or supplied large normalization can move the cost into the function. Sigma is not selected by this identity.

Illustrative flat-branch bookkeeping for lambda~1, mu~M_P, mu^4 Omega~10^240: affine slope 1 requires moving the hierarchy into mu~10^-60 M_P; the frozen cubic requires |v|/mu^4~10^120. For the **review's analytic comparison**, sigma'=cosh z requires |z|=553.313569499131, S_v~138.32839237478274, and a 1% increase in |v| at fixed volume gives a mass factor ~3.9879763578032428. These log estimates use ln cosh z=|z|-ln2 with exponentially negligible error at this magnitude. The comparison is not a replacement CP-04 law, freeze, solver run or model selection. No alternative sigma is adopted. z is fixed by v/mu^4, not a free adjustment.

## The necessary total-response condition

For a differentiable self-consistent family at fixed mu and with sigma'>0,

\[
\frac{d\ln\lambda}{d\ln\Omega}=-\frac14+
\frac14\frac{\sigma''}{\sigma'}\frac{dz}{d\ln\Omega}.
\]

Exact scale stability requires (sigma''/sigma') dz/dlnOmega=1. For the review's sinh comparison, this is tanh(z) dz/dlnOmega=1, hence dz/dlnOmega=coth z: approximately +1 on the large positive branch and -1 on the large negative branch. For v>0, the flat branch has z<0. This sign must not be silently lost.

More generally on a flat-vacuum solution family, fixed mu but variable Omega(v),

\[
\frac{d\ln\lambda}{d\ln|v|}=S_v-\frac14\frac{d\ln\Omega}{d\ln|v|}.
\]

Scale stability therefore requires the derived volume response dlnOmega/dln|v|=z sigma''/sigma'. Imposing that response is not deriving it. In the frozen cubic case it requires 2v^2/(mu^8+v^2), approaching 2; the finite-slab fixed-volume control does not satisfy it. No cosmological family satisfying it has been calculated or independently calibrated here. At fixed v,mu on the flat branch dz/dlnOmega=0, so no exponential form alone changes -1/4.

## Disposition

Retain CP-04's conditional verdict and frozen packet unchanged. The sharpening is a three-way distinction: geometrical-history restriction, full-action symmetry, and observable scale selection. A protective cancellation must not be credited with selecting the physical function, history or matter calibration. The literature result, algebraic identities and conditional response requirement are not a new GRUT law. Stage 3 inactive, Card 3 unopened, no CP-05 or canonical integration. Fundamental construction remains active.
