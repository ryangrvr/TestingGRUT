# Candidate 1: nonlinear scalar evolution and conditional escape/fission bounds

18 October 2026. Derivation for the existing Candidate 1 action at its frozen couplings. This file supplies equations and interpretation limits for the new stability calculations. It does not alter the candidate action or the original candidate reports.

“Exact” below means the nonlinear flat-space classical equations of the already specified leading scalar action, with reporter \(D=0\), ordinary matter in its prepared vacuum, and gravity omitted at this order. It does not mean exact quantum dynamics, nonlinear Einstein evolution, or ultraviolet completion. The existing weak-field metric has potentials of order \(10^{-9}\); its feedback is a separate correction, not part of the following exact flat-space conservation statements.

## 1. Variables and normalization

Let \(g=g_X\), and use

\[
\tau=gv t,\qquad \mathbf y=gv\mathbf x,\qquad
A(\tau,\mathbf y)=\frac{\sqrt2 X(t,\mathbf x)}v,
\qquad S(\tau,\mathbf y)=\frac{\phi(t,\mathbf x)}v.
\]

Set \(a=\lambda_\phi/g^2=1\), \(b=\lambda_X/g^2=0.1\). A dot means \(\partial_\tau\), and gradients are with respect to \(\mathbf y\). The complex field \(A=u+iw\) has two canonically normalized real components in the dimensionless density

\[
\mathcal L_*=
\frac12|\dot A|^2-\frac12|\nabla A|^2
+\frac12\dot S^2-\frac12|\nabla S|^2-V(S,A),
\]

\[
V=\frac a4(S^2-1)^2+\frac12S^2|A|^2+\frac b4|A|^4.
\]

The physical scalar action is \(g^{-2}\int d\tau\,d^3y\,\mathcal L_*\). The overall factor does not change the classical equations. Physical energy and charge are

\[
E=\frac vg\bar E,\qquad Q=\frac{\bar Q}{g^2}.
\]

Thus the free unit-charge carrier has dimensionless rest energy one: \(m_X=gv\) and \(E_{\rm free}=m_XQ\) correspond to \(\bar E_{\rm free}=\bar Q\). There is no additional \(4\pi\) in the three-dimensional integrals below; it appears only after reducing a radial integral.

## 2. Exact nonlinear initial-value equations

Varying the displayed density gives

\[
\boxed{\ddot A-\nabla^2 A+S^2A+b|A|^2A=0,}
\]

\[
\boxed{\ddot S-\nabla^2S+a(S^2-1)S+S|A|^2=0.}
\]

For an implementation using real fields,

\[
\ddot u=\nabla^2u-[S^2+b(u^2+w^2)]u,
\quad
\ddot w=\nabla^2w-[S^2+b(u^2+w^2)]w,
\]

\[
\ddot S=\nabla^2S-[a(S^2-1)+u^2+w^2]S.
\]

Initial data are \((A_0,\Pi_0,S_0,P_0)\), with \(\Pi=\dot A\) and \(P=\dot S\). The finite-energy vacuum is \(A,\Pi,P\to0\), \(S\to1\), and \(\nabla S\to0\). The real reporter stays exactly zero in these classical equations if both of its initial data are zero. No new controller, forcing term, dissipation law, or evolution rule is introduced.

The existing stationary profile is

\[
A_*(\tau,\mathbf y)=e^{-i\Omega\tau}F(r),\qquad
S_*(\tau,\mathbf y)=S_*(r),\qquad r=|\mathbf y|,
\]

with \(\Omega\simeq0.8000000019763\). Substitution recovers exactly the already solved radial equations

\[
F''+\frac2rF'=(S_*^2+bF^2-\Omega^2)F,
\qquad
S_*''+\frac2rS_*'=[a(S_*^2-1)+F^2]S_*.
\]

## 3. Conserved energy, charge, and fluxes

The dimensionless energy and positive-charge convention are

\[
\boxed{\bar E=\int_{\mathbb R^3}d^3y\left[
\frac12|\Pi|^2+\frac12|\nabla A|^2+
\frac12P^2+\frac12|\nabla S|^2+V\right],}
\]

\[
\boxed{\bar Q=-\int_{\mathbb R^3}d^3y\,\operatorname{Im}(A^*\Pi).}
\]

For \(A_*=F e^{-i\Omega\tau}\), the charge density is \(\Omega F^2>0\). The exact local continuity laws are

\[
\partial_\tau q+\nabla\cdot\mathbf J_Q=0,
\quad q=-\operatorname{Im}(A^*\Pi),
\quad \mathbf J_Q=\operatorname{Im}(A^*\nabla A),
\]

\[
\partial_\tau e+\nabla\cdot\mathbf J_E=0,
\quad \mathbf J_E=-\operatorname{Re}(\Pi^*\nabla A)-P\nabla S.
\]

For a fixed finite region \(B\),

\[
\frac{d\bar Q_B}{d\tau}=-\int_{\partial B}\mathbf J_Q\cdot\mathbf n\,d\Sigma,
\qquad
\frac{d\bar E_B}{d\tau}=-\int_{\partial B}\mathbf J_E\cdot\mathbf n\,d\Sigma.
\]

Charge/energy leaving a core is not a violation of the global conservation laws. An absorbing layer deliberately changes the closed equations in that layer; its removed energy and charge must be included in a balance ledger. A reflecting/periodic box conserves a box quantity while preventing genuine escape to spatial infinity.

Total momentum is \(\mathbf P_{\rm tot}=\int\mathbf J_E\,d^3y\). A perturbation may introduce translation, boost, phase drift, or angular momentum without destroying the soliton. Shape comparisons should remove the appropriate translation and global phase rather than interpreting these exact symmetries as internal instability.

The reference numbers from the existing `grut_candidate_1_20261018/code/results.json` are

\[
\bar Q_*=399.4384624223012,\qquad
\bar E_*=363.9229985516867,
\]

\[
I_*=\int F^2d^3y=499.29807679441456,\qquad
T_{\rm grad,*}=66.558341736314.
\]

They are numerical values, not interval-certified constants. \(\bar Q_*=\Omega I_*\) supplies a useful charge-normalization check.

## 4. Preparing genuinely fixed-charge perturbations

### General exact charge projection of the momentum

First choose the perturbed fields \(A_\epsilon,S_\epsilon\) and provisional velocities \(\Pi_{\rm raw},P_\epsilon\). Define

\[
N_\epsilon=\int|A_\epsilon|^2d^3y>0,\qquad
Q_{\rm raw}=-\int\operatorname{Im}(A_\epsilon^*\Pi_{\rm raw})d^3y.
\]

Then

\[
\boxed{\Pi_\epsilon=\Pi_{\rm raw}
-i\frac{\bar Q_{\rm target}-Q_{\rm raw}}{N_\epsilon}A_\epsilon}
\]

has exactly \(\bar Q=\bar Q_{\rm target}\). This is the orthogonal projection of the provisional velocity onto the charge constraint in the real \(L^2\) velocity inner product. Apply it on the **actual evolution grid with the actual quadrature**, not only to a radial continuum profile. Record the achieved charge and resulting energy after interpolation, truncation and projection.

For fixed fields and fixed charge, the least kinetic energy occurs for

\[
\boxed{\Pi_\epsilon=-i\Omega_\epsilon A_\epsilon,\qquad
\Omega_\epsilon=\bar Q_{\rm target}/N_\epsilon,\qquad P_\epsilon=0.}
\]

This follows from Cauchy–Schwarz. More generally decompose

\[
\Pi=-i\frac{\bar Q}{N}A+\Pi_\perp,
\qquad -\operatorname{Im}\int A^*\Pi_\perp\,d^3y=0.
\]

The energy is exactly

\[
\bar E=\frac{\bar Q^2}{2N}
+\int[\tfrac12|\nabla A|^2+\tfrac12|\nabla S|^2+V]d^3y
+\tfrac12\|\Pi_\perp\|_2^2+\tfrac12\|P\|_2^2.
\]

This identifies the fixed-charge energy functional used in the earlier Hessian calculation. It also shows that arbitrary extra velocity perturbations add nonnegative energy once their charge component has been separated.

If instead one insists on retaining the original \(\Omega\), rescale the chosen shape by

\[
A_\epsilon\mapsto
\left[\frac{\bar Q_{\rm target}}{\Omega N_\epsilon}\right]^{1/2}A_\epsilon,
\qquad \Pi_\epsilon=-i\Omega A_\epsilon.
\]

These are different physical perturbation families. Changing the amplitude alone while retaining \(\Pi=-i\Omega A\) without normalization changes the charge. No general perturbation simultaneously preserves charge, energy, momentum and a prescribed frequency automatically. A real field deformation with purely rotating velocity and \(P=0\) has zero net momentum; a spatially varying complex phase can impart momentum and must be checked separately.

### An analytic nonradial perturbation with a known energy cost

An especially useful continuum family uses a volume-preserving linear deformation of **both** profiles:

\[
M_\epsilon=\operatorname{diag}(e^\epsilon,e^{-\epsilon},1),\quad
A_\epsilon(\mathbf y)=F(|M_\epsilon\mathbf y|),\quad
S_\epsilon(\mathbf y)=S_*(|M_\epsilon\mathbf y|),
\]

\[
\Pi_\epsilon=-i\Omega A_\epsilon,\qquad P_\epsilon=0.
\]

Since \(\det M_\epsilon=1\), the continuum \(N,\bar Q\), potential energy and rotational kinetic energy are unchanged. Spherical symmetry of the original gradient tensor gives

\[
\boxed{\bar E_\epsilon-\bar E_*
=\left[\frac{\operatorname{tr}(M_\epsilon M_\epsilon^T)}3-1\right]T_{\rm grad,*}
=\frac43T_{\rm grad,*}\sinh^2\epsilon.}
\]

For \(\epsilon=0.05,0.1,0.2,0.3\), the continuum energy increments are approximately \(0.2220461,0.8904067,3.597362,8.229505\). This supplies a nonradial fixed-charge perturbation with an independently known initial energy cost. Its first-order deformation is proportional to \(rF'(r)(\hat y_x^2-\hat y_y^2)\), an \(\ell=2\) perturbation; finite \(\epsilon\) includes higher multipoles. Finite-grid interpolation and clipping alter the exact identities; use on-grid charge projection and the actual discretized energy rather than substituting the analytic values into simulation diagnostics. The integrating worker selected \(\epsilon=0.05,0.2\) in a protocol amendment before any physical evolution trajectory was run; no scalar coupling changed.

## 5. Linearized equations and the meaning of the fixed-charge Hessian

For comparison with the earlier spectral work write

\[
A=e^{-i\Omega\tau}(F+u+iw),\qquad S=S_*+z,
\]

with real small \(u,w,z\). The exact first-order linearization is

\[
\ddot u+2\Omega\dot w+
[-\nabla^2+S_*^2+3bF^2-\Omega^2]u+2FS_*z=0,
\]

\[
\ddot w-2\Omega\dot u+
[-\nabla^2+S_*^2+bF^2-\Omega^2]w=0,
\]

\[
\ddot z+[-\nabla^2+a(3S_*^2-1)+F^2]z+2FS_*u=0.
\]

The first-order charge perturbation is

\[
\delta\bar Q=2\Omega\int Fu\,d^3y-\int F\dot w\,d^3y.
\]

Thus setting \(\dot w=0\) with an arbitrary amplitude perturbation is generally **not** the fixed-charge linear problem. Minimizing over the charge-compatible velocities reproduces the positive amplitude rank-one term

\[
\frac{4\Omega^2}{I_*}|(F,0)\rangle\langle(F,0)|
\]

in the \((u,z)\) Hessian, in addition to its local operator. Only the spherical amplitude component couples to this global rank-one term. The phase operator has the exact zero mode \(w\propto F\); translations give the usual \(\ell=1\) modes. Static constrained Hessian eigenvalues are not automatically the squares of all dynamical frequencies because the rotating-frame equations include the displayed velocity couplings.

If a continuum constrained energy has a strict local minimum modulo phase/translations and controls a specified distance to the soliton orbit, conserved energy and charge can support a Lyapunov/orbital-stability argument. A finite-grid Hessian and finitely many evolutions do not establish those coercivity, compactness or all-time continuum hypotheses. No such theorem is asserted here.

## 6. A sharp universal lower bound and a free-escape threshold

For the present binding class \(ab<1\), define

\[
\Omega_0=(ab)^{1/4}=0.5623413251903491.
\]

The following exact square identity holds for every real \(S\) and complex \(A\):

\[
\boxed{V-\frac12\Omega_0^2|A|^2=
\frac14[\sqrt a(S^2-1)+\sqrt b|A|^2]^2
+\frac12(1-\sqrt{ab})S^2|A|^2\ge0.}
\]

For positive total charge, completing the velocity square gives

\[
\bar E-\Omega_0\bar Q=\int d^3y\left[
\frac12|\Pi+i\Omega_0A|^2+\frac12P^2
+\frac12|\nabla A|^2+\frac12|\nabla S|^2
+V-\frac12\Omega_0^2|A|^2\right]\ge0.
\]

Changing the sign of the square's rotation treats negative charge. Consequently every finite-energy state of this leading scalar model obeys

\[
\boxed{\bar E\ge\Omega_0|\bar Q|.}
\]

This is a rigorous bound for the displayed model, not a lower bound fitted to the found branch. It is much weaker than the target's actual energy. For the target, \(\Omega_0\bar Q_*=224.6207543\), while \(\bar E_*=363.9229986\).

In contrast, a free complex Klein–Gordon field about \(S=1\) has dimensionless mass one and obeys \(\bar E_{\rm free}\ge|\bar Q_{\rm free}|\), by the same velocity square with frequency one. Neutral scalar radiation costs additional nonnegative energy. Therefore **complete asymptotic scattering into free vacuum excitations** is energetically impossible when

\[
\boxed{\bar E_{\rm initial}<|\bar Q_{\rm target}|.}
\]

The target's available positive perturbation budget for this conclusion is

\[
\boxed{\Delta\bar E<\bar Q_*-\bar E_*=35.5154638706145.}
\]

For the analytic ellipsoidal family above, this remains satisfied for \(|\epsilon|<0.5965885\), using the numerical continuum profile. This is an energy exclusion, not a theorem that every such strongly deformed state returns near the original soliton.

If an eventual separated state consists of free radiation plus other bound structures, let their net charges be \(Q_f,Q_b\), with \(Q_f+Q_b=\bar Q_*\). Additivity at asymptotic separation and the two energy bounds give

\[
\bar E_{\rm initial}\ge|Q_f|+\Omega_0|Q_b|,
\qquad
|Q_b|\ge\frac{\bar Q_*-\bar E_{\rm initial}}{1-\Omega_0}
\quad(\bar E_{\rm initial}<\bar Q_*).
\]

At the unperturbed energy this lower bound is \(|Q_b|\ge81.1487717\), about \(20.3157\%\) of the target charge. “Bound structures” here means all charge not converted to the free radiation component; it can include several outgoing solitons far from the original core. This conditional decomposition is not a proven scattering description of arbitrary solutions. In particular the bound does not guarantee central-core retention, prevent fission, fix a final radius, select the original profile, or supply a lifetime.

The proposed analytic perturbations have the following **continuum reference** budgets:

| \(\epsilon\) | \(\Delta\bar E\) | \(\bar E_{\rm initial}\) | \(\bar Q_*-\bar E_{\rm initial}\) | Conditional minimum \(|Q_b|\) | Fraction of target charge |
|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 363.9229986 | 35.5154639 | 81.1487717 | 0.2031571 |
| 0.05 | 0.2220461 | 364.1450446 | 35.2934178 | 80.6414218 | 0.2018870 |
| 0.2 | 3.5973618 | 367.5203603 | 31.9181021 | 72.9292116 | 0.1825793 |

Replace these reference energies by the measured discrete initial energy for the actual run. Their many printed digits do not certify continuum interpolation, mesh, tail or time-step errors.

## 7. Conditional fission and subadditivity tests

Let \(E_{\rm br}(q)\) denote a numerical stationary branch at the **same** couplings and matching prescription. For two same-sign branch fragments, the rest-energy fission difference is

\[
\Delta_2(q)=E_{\rm br}(q)+E_{\rm br}(\bar Q_*-q)-E_{\rm br}(\bar Q_*),
\quad 0<q<\bar Q_*.
\]

If the tested initial state has energy \(\bar E_{\rm initial}\), that listed channel is energetically excluded only when

\[
E_{\rm br}(q)+E_{\rm br}(\bar Q_*-q)>\bar E_{\rm initial}.
\]

The criterion must use the perturbed energy, not automatically the stationary energy. Relative-motion kinetic energy and radiation increase the final threshold, so the rest sum is a conservative channel lower bound **if** the fragment energies being used are the true lowest accessible fragment energies. A known stationary branch supplies available competitor states, not automatically lower bounds on every possible fragment. A negative splitting difference supplies a lower-energy separated competitor, but does not establish a barrier-free dynamical path or its timescale. The present scan uses continuous classical charge; exact quantum charge sectors have physical integer \(Q\), corresponding to increments \(g^2=0.01\) in \(\bar Q\), and require quantum corrections before a precise quantum fission claim.

For several branch fragments and free charge, the corresponding proposed final rest energy is

\[
E_{\rm channel}=\sum_{j=1}^n E_{\rm br}(q_j)+|q_f|,
\qquad \sum_jq_j+q_f=\bar Q_*.
\]

Checking only equal halves or a handful of partitions does not test all fission channels. Charge sectors below the known branch's minimum charge may be represented by free carriers **as a comparator**, but absence of a numerical solution there is not a proof that no lower-energy object exists. Opposite-sign fragments, disconnected solution branches, multi-core bound configurations and quantum channels require their own treatment if included in a global claim.

For a differentiable stationary branch at fixed couplings, variation of the action gives the first law

\[
\frac{dE_{\rm br}}{dq}=\Omega(q),\qquad
\frac{d}{dq}\left(\frac{E_{\rm br}}q\right)
=\frac{\Omega q-E_{\rm br}}{q^2}.
\]

At the target, the second expression is approximately \(-2.7810671\times10^{-4}\). This is a local slope fact, not a global subadditivity proof. If the scan establishes that \(E_{\rm br}(q)/q\) is strictly decreasing over **all** lower positive fragment charges on the relevant branch, then every finite positive-charge partition within that branch has

\[
\sum_jE_{\rm br}(q_j)>E_{\rm br}(\bar Q_*).
\]

Free fragments with ratio one also lie above the target ratio \(0.91108652\). These conclusions are conditional on branch coverage and interpolation/error control; they do not exclude an uncomputed lower branch.

For emission of a free positive charge \(q\) leaving a remnant on the same continuous branch,

\[
E_{\rm br}(\bar Q_*-q)+q-E_{\rm br}(\bar Q_*)
=\int_{\bar Q_*-q}^{\bar Q_*}[1-\Omega(q')]\,dq'.
\]

This is positive wherever the entire remnant interval has \(\Omega<1\); infinitesimal emission at the target costs about \((1-\Omega_*)q\simeq0.2q\). A perturbed state can pay this cost from its excess energy, and fission can involve different remnants, so this observation does not prove indefinite charge retention.

For a genuinely global variational statement define

\[
\mathcal E(q)=\inf\{\bar E:\bar Q=q\}
\]

over the full continuum finite-energy state class. Widely separated trial configurations imply the non-strict inequality

\[
\mathcal E(q_1+q_2)\le\mathcal E(q_1)+\mathcal E(q_2).
\]

Strict binding \(\mathcal E(Q)<\mathcal E(q)+\mathcal E(Q-q)\) against every nontrivial split is a relevant ingredient of a concentration-compactness/orbital-stability argument. Existence, compactness, coercivity and control of all admitted competitors must still be shown. A finite stationary branch scan approximates part of that variational problem; it is not a proof that \(E_{\rm br}=\mathcal E\).

## 8. Finite-grid and boundary controls

### Fourier collocation used by the integrating worker

For a periodic cubic grid with spacing \(h\), \(n\) total sites, and unnormalized forward DFT \(\widehat A_k=\sum_i A_i e^{-i\mathbf k\cdot\mathbf y_i}\), the Hamiltonian consistent with a spectral Laplacian is

\[
\bar E_h=h^3\sum_i[\tfrac12|\Pi_i|^2+\tfrac12P_i^2+V_i]
+\frac{h^3}{2n}\sum_k |\mathbf k|^2
(|\widehat A_k|^2+|\widehat{S-1}_k|^2),
\]

\[
\bar Q_h=-h^3\sum_i\operatorname{Im}(A_i^*\Pi_i),
\qquad \mathbf k=2\pi\mathbf m/L.
\]

An equivalent Parseval-normalized implementation is fine; real-FFT storage requires the appropriate multiplicity weights. Use this spectral gradient energy, rather than the finite-difference expression below, for the actual Fourier-collocation runs. The real local nonlinear force on \(A\) is a coefficient times \(A\), so its exact velocity kick preserves charge pointwise. The exact linear Klein–Gordon Fourier flow also preserves global charge. Their symmetric split composition therefore preserves \(\bar Q_h\) up to floating-point error; it does not preserve the full \(\bar E_h\) exactly. The energy/time-step comparison is an independent check. Cubic aliasing can alter continuum accuracy even if this discrete charge conservation is excellent; inspect high-frequency spectral power and refinement, and document any filtering/dealiasing rule.

The spectral derivative has global grid support, so a continuum causal boundary estimate is not a strict theorem of the collocation scheme. Domain and refinement comparisons are still required.

### Alternative nearest-neighbor check

For a standard nearest-neighbor seven-point Laplacian on a cubic mesh, use the associated edge/forward-difference gradient energy,

\[
\bar E_h=h^3\sum_i\left[\tfrac12|\Pi_i|^2+\tfrac12P_i^2+V_i\right]
+\frac h2\sum_{i,j}|A_{i+e_j}-A_i|^2
+\frac h2\sum_{i,j}(S_{i+e_j}-S_i)^2,
\]

with each undirected edge counted once and boundary handling stated. Centered first-derivative squares give a different lattice Hamiltonian and can report spurious energy drift when paired with that Laplacian. The semi-discrete Hamiltonian system conserves its own energy and charge; a time integrator introduces additional, separately measured errors. Record whether velocities are integer-time or half-step variables when evaluating diagnostics.

### Boundary return and interpretation

For nearest boundary distance \(B\) from the center and measured core radius \(R_c\), let \(d=B-R_c\). In the continuum, a signal emitted from the core and reflected from the boundary, or wrapped through a periodic box, cannot revisit the core earlier than approximately \(2d\). An initial tail/boundary mismatch or exterior forcing already present at \(\tau=0\) can reach the core after \(d\). Quantify the initial tail residual separately; the earlier \(2d\) estimate is not a blanket guarantee that boundary errors are absent. These are continuum causal estimates, not exact claims about the dispersive finite-grid scheme.

An absorbing buffer can reduce returns but introduces non-Hamiltonian terms there. Report cumulative sink/flux, test buffer/domain dependence, and keep the diagnostic core separated from the buffer. A closed periodic/reflected run that lasts beyond the return time can reabsorb radiation and appear artificially persistent. It cannot establish continuum escape or damping. Conversely, loss from a finite diagnostic sphere alone is not evidence that bound charge has disintegrated.

The useful reported quantities are actual initial \((\bar Q_h,\bar E_h)\), their numerical conservation/balance errors, core and total charge/energy, outgoing flux, translation/phase-adjusted shape measures, domain/boundary horizon, mesh/time-step comparisons, perturbation amplitude and duration. The reference rotation period is \(2\pi/\Omega\simeq7.85398\); counting periods is descriptive and does not turn a finite observation into an all-time result.

The defensible numerical conclusion is persistence or failure of the stated packets under the stated perturbations, resolution, boundary conditions and finite observation interval. Continuum nonlinear orbital stability, all-partition fission exclusion, quantum lifetime and empirical validity remain separate claims unless further proofs/calculations establish them.

## Provenance

The evolution, normalization, square identities, ellipsoidal energy increment and conditional inequalities above are direct derivations from the frozen Candidate 1 action. Numerical reference values are copied from the existing `grut_candidate_1_20261018/code/results.json`; they are not new observations. The established conserved-charge soliton framework is identified by S. Coleman's original paper, *Q-balls*, https://doi.org/10.1016/0550-3213(85)90286-X . No general nonlinear-stability theorem or historical FLS numerical bound is imported into this particular branch calculation.
