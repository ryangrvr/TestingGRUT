# CP-05: a working covariant response model and its actual boundary

**Verdict: WORKING LEADING SCALAR–TENSOR MODEL; NO NEW GRUT LAW OR EMPIRICAL VALIDATION.** One new, separately accounted action. Seven author check groups pass. No claim to explain all reality, originate spacetime, or complete quantum gravity. Prior packets untouched.

This branch prioritizes an explicit physical response over another proposed delay kernel. The extra response is produced by a local degree of freedom. Its propagation, static range, matter coupling and memory are calculated together. The same calculation establishes that this is conventional scalar–tensor EFT: calling it a responsive substrate adds no new physical exclusion.

## 1. THE LAW

The frozen action [LAW.md](LAW.md), CP-05/D006, is

\[
S=\int d^4x\sqrt{-g}\left[\frac12F(\phi)R
-\frac12g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi
-\frac12m^2\phi^2\right]+S_m[g,\psi]+S_{\rm GHY},
\quad F=M^2+2g_c\phi.
\]

M,m>0, g_c real, F>0. The Lorentzian metric, Einstein kinetic term, matter inventory and universal metric coupling, causal preparation, quantization/probabilities and three matching constants are POSTULATED. They are not generated. F>0 is the positive-tensor-kinetic domain, not a global dynamical guarantee. Weak fields about phi=0, asymptotically flat retarded boundary data and energies below an EFT cutoff define the calculation's domain.

All applicable covariant EFT operators and their finite matching coefficients must be admitted in quantum calculations. The displayed finite action specifies the leading tree approximation; it is not a UV-complete three-parameter quantum theory. A constant/cosmological matching term must be renormalized to support the selected asymptotic vacuum. The scalar mass is not symmetry-protected by the specified action. There is no formal Stage-3 L0 compression score for this separate program; no credited information reduction from hiding the standard physics in the word substrate.

## 2. THE DERIVATION: local equations and positivity

Exact classical variations give

\[
F G_{\mu\nu}=T_{\mu\nu}+\partial_\mu\phi\partial_\nu\phi
-g_{\mu\nu}\left[\tfrac12(\partial\phi)^2+\tfrac12m^2\phi^2\right]
+\nabla_\mu\nabla_\nu F-g_{\mu\nu}\Box F,
\qquad (\Box-m^2)\phi=-g_c R.
\]

Matter is covariantly conserved on its equations because it is minimally coupled in this Jordan frame. The metric Noether identity and scalar equation jointly enforce total consistency; a bare phenomenological dissipative stress is not inserted into Einstein's equation.

Define g^E_munu=(F/M^2)g_munu and A(phi)=M/sqrt(F). The scalar kinetic coefficient in the Einstein frame is

\[
Z_E(\phi)=\frac{M^2}{F}+\frac{3M^2}{2}\left(\frac{F'}F\right)^2>0.
\]

Potential U_E=M^4 m^2 phi^2/(2F^2) is nonnegative on this branch. Around phi=0,

\[
Z=1+\frac{6g_c^2}{M^2},\quad \chi=\sqrt Z\,\phi,
\quad \mu^2=\frac{m^2}{Z},\quad
\beta=M\frac{d\ln A}{d\chi}\bigg|_0=-\frac{g_c}{M\sqrt Z}.
\]

Two massless tensor helicities plus one positive-kinetic massive scalar propagate in the weak-field theory. The two tensor helicities are supplied by Einstein dynamics; the scalar is also an input. Regular GR constraints, including lapse/shift elimination, leave these three local degrees of freedom. This is correspondence, not emergence of constraints.

An independent Jordan-frame route checks the physical mass. To first order the metric trace is -M^2 R=T-6g_c Box phi. Substituting into the scalar equation gives

\[
\boxed{Z\Box\phi-m^2\phi=\frac{g_c}{M^2}T},
\qquad \boxed{(\Box-\mu^2)\phi=\frac{g_c}{M^2 Z}T.}
\]

The original mass m and the mixed physical mass mu cannot be interchanged.

## The covariant memory is derived, not chosen

Let G_m^ret solve (m^2-Box_x)G_m^ret(x,y)=delta_g(x,y), with zero support outside the causal future of y. The exact scalar equation on a given metric gives

\[
\phi(x)=\phi_{hom}(x)+g_c\int d^4y\sqrt{-g(y)}G_m^{ret}(x,y)R(y).
\]

The homogeneous solution records independently specified initial data. Setting it to zero is a preparation, not a theorem erasing information. After accounting for linear metric mixing, the physical source-response kernel instead uses G_mu^ret and phi=-g_c G_mu^ret*T/(M^2 Z).

In flat spacetime,

\[
G_\mu^{ret}(\omega,\mathbf k)=
\frac1{\mathbf k^2+\mu^2-(\omega+i0)^2},\qquad
G_\mu^{ret}(t,\mathbf k)=\theta(t)\frac{\sin(\omega_k t)}{\omega_k},
\quad \omega_k=\sqrt{\mathbf k^2+\mu^2}.
\]

With s=t^2-r^2 its four-dimensional position-space form is

\[
G_\mu^{ret}(t,r)=\theta(t)\left[
\frac{\delta(t^2-r^2)}{2\pi}
-\theta(t^2-r^2)\frac{\mu J_1(\mu\sqrt{t^2-r^2})}{4\pi\sqrt{t^2-r^2}}
\right].
\]

This proves causal support, including the massive tail inside the light cone. Lorentz invariance is compatible with retarded response and time orientation. No preferred fluid velocity is needed. At fixed wave number memory oscillates; it is not a positive decaying exponential with an integrable first moment. At late timelike separation the Bessel tail is oscillatory and power-law. Retarded preparation does not make the underlying local theory fundamentally irreversible.

Formally eliminating the Gaussian scalar contributes g_c^2 R(m^2-Box)^-1R/2 to a symmetric or in-out effective action. **A retarded Green function cannot simply be inserted into an ordinary single-history action and claimed to yield that causal equation:** variation symmetrizes the kernel. The physical response above comes from the local initial-value problem, or an appropriately specified in-in calculation. No path-integral determinant is omitted while claiming an exact quantum effective action; that determinant would add curvature counterterms/nonlocal loop terms.

## 3. THE RESTRICTION: same coefficients across observables

Write weak Jordan metric ds^2=-(1+2Psi)dt^2+(1-2Phi)d x^2. For a weak, unscreened point source of mass M_s,

\[
\Psi=-\frac{G_T M_s}{r}(1+\eta e^{-\mu r}),\quad
\Phi=-\frac{G_T M_s}{r}(1-\eta e^{-\mu r}),
\quad G_T=\frac1{8\pi M^2},\quad
\eta=2\beta^2=\frac{2g_c^2}{M^2+6g_c^2}.
\]

Jordan conformal transformation of the Einstein-frame tensor potentials adds opposite-sign scalar pieces; the static scalar equation gives phi=g_c M_s e^-mu r/(4pi M^2 Z r). This independently gives the coefficients above. The scalar-mediated force is attractive regardless of the sign of g_c.

Thus, at leading order,

\[
0\le\eta<\tfrac13,\quad
\frac{\Phi}{\Psi}=\frac{1-\eta e^{-\mu r}}{1+\eta e^{-\mu r}},\quad
\frac{\Psi+\Phi}{2}=-\frac{G_T M_s}{r},\quad
\frac{|a_r|}{G_T M_s/r^2}=1+\eta(1+\mu r)e^{-\mu r}.
\]

The dispersion of propagating scalar radiation is omega^2=k^2+mu^2. Static inverse range and wave threshold therefore share the same independently calibrated mu. A scalar radiation signal still needs a radiating source; a universal weak coupling can cancel leading dipole radiation for equal charge/mass ratios, so no unjustified dipole forecast is made. Compact-body sensitivities/nonlinear scalarization lie outside the point-source formula.

For eta>0, a static measurement of eta and mu plus tensor G_T determines

\[
g_c^2=\frac{M^2\eta}{2(1-3\eta)},\qquad
m^2=\frac{\mu^2}{1-3\eta}.
\]

The sign of g_c is a phi sign convention here. At eta=0 the scalar decouples and its mass is not inferred from gravitational data. No scalar-wave measurement may then be claimed as a parameter-free prediction.

These are restrictions of the specified family, not definitions. But the canonical Jordan kinetic term is load-bearing: other healthy scalar-tensor kinetic/coupling functions can broaden the range. The bound eta<1/3 is not a theorem for all gravity or a novel GRUT law.

## 4. THE COUNTERMODEL and novelty audit

GR with the same tensor G_T but no coupled scalar has eta=0, no Yukawa slip and no scalar radiation. It is a legitimate empirical alternative, not the entire conventional baseline. Full conventional scalar-tensor EFT already produces every displayed CP-05 observable, because CP-05 is a member of that class. It is therefore not possible to exhibit a same-assumption conventional model excluded by an extra GRUT principle here: there is no such extra principle.

CP-05's specific predictions can exclude other phenomenological models if calibrated and tested, but this does not establish originality or explain why nature chose its scalar and couplings. Renaming a local scalar as vacuum memory would hide this identification; the report explicitly refuses it. A horizon thermality derivation, including GRUT III, does not establish this extra scalar exists either.

Prior literature: Salgado, [gr-qc/0202082](https://arxiv.org/abs/gr-qc/0202082), analyzes covariant/Jordan linear scalar-tensor gravity and weak-field deviations. The primary abstract was inspected; exact normalizations here were derived from the displayed action rather than attributed to an unread equation. No originality claimed.

## 5. THE TEST and attempts to break it

Run `python research/constructive_physics/attempt05/check.py`; standard-library only. Seven groups cover exact mixing/inverse calibration, positivity/strength boundary, independently integrated retarded mode versus analytic step response, static slip/lensing, positive pole diagnostic, common static/dynamical mass and prior/freeze integrity.

The RK4 solver with frozen dt=0.001 and five sample times agrees with q(t)=(1-cos omega_k t)/omega_k^2 to 1.33e-13 maximum absolute error. This independently checks the time-domain response, not a fitted spectral reconstruction. It is neither a quantum loop calculation nor a physical measurement.

Frozen software control M=1,g_c=1/2,m=1 gives Z=5/2, mu^2=2/5, eta=1/5. At r=1/mu the potential scalar fraction is eta/e, the slip is about 0.863 and the force enhancement about 14.7%. Units and numbers are invented controls, not an experimental anomaly or allowed physical parameter region.

An operational test would independently match tensor G_T, source masses and the static (eta,mu), then predict slip/lensing and the scalar-wave threshold without refitting. Calibration must distinguish tensor G_T from a short-range effective Newton constant. Extended sources require integrating the Yukawa Green function; backgrounds and apparatus sensitivities must be quantified. None is supplied here. Existing fifth-force/GW/cosmological constraints must be applied to a physical point before calling it viable in nature; no unverified claim of experimental permission is made.

Failure cases: F<=0 loses the stated healthy tensor domain; m^2<0 would be tachyonic and is excluded by the frozen law; F crossing, strong sources, finite background phi, self-interaction and EFT operators can invalidate the linear forecast. Integrating out the scalar leaves a nonlocal reduced description, but reinstating it gives local second-order evolution: memory alone is not a new architecture. Quantum corrections generate additional interactions and curvature terms; this work does not reproduce CP-01's mistake of claiming exact finite-family closure.

## 6. THE VERDICT and scope earned

The leading theory is mathematically explicit and has positive local weak-field modes, covariant causal response, conservation and falsifiable cross-observable relations conditional on independent matching. It is a **working effective physical model**, not an empirically established model of our universe and not a fundamental theory of reality.

It makes an actual advance over the phenomenological delay proposal: memory has a specified degree of freedom, causal light-cone support, an explicit retained-state interpretation and a common static/dynamic pole. However, established scalar-tensor physics supplies the whole result. The omitted origin question remains genuinely omitted: no Lorentzian signature, Einstein field, quantum rules, observer partition, fundamental arrow or coupling values are generated.

Author evidence only; independent review pending. One execution, no fixture failures or post-freeze law changes. No automatic CP-06 or Stage-3 Card 3. Latest old research and GRUT III preserved; no merging or banking. The broad scientific objective remains active, with no false claim that conventional correspondence has solved it.
