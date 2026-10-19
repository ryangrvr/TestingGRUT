# Next discriminating calculation: the same finite-charge quantum state

This is a specification of work **not executed in this packet**. The completed scalar-loop calculation is reported in [REPORT.md](REPORT.md). The purpose of the next calculation is to decide whether the demonstrated classical object has a controlled quantum counterpart at the same action/matching. Its result may reject that proposed regime.

## Fixed inputs and quadratic starting point

Retain the original four scalar fields, benchmark couplings and matched operators. Use the hashed classical profile at \(\Omega=\omega/(g_Xv)\simeq0.8\), with \(Q\simeq39943.85\), as a starting saddle. Do not replace it by a constant interior. Write

\[
X=\frac{e^{-i\omega t}}{\sqrt2}\,[f(r)+u(t,\mathbf x)+i w(t,\mathbf x)],
\qquad \phi=\varphi(r)+\eta(t,\mathbf x),\quad D=\zeta(t,\mathbf x).
\]

With the polynomial coefficients defined in [RG_DERIVATION.md](RG_DERIVATION.md), the tree linearized equations derived from the same action are

\[
\ddot u+2\omega\dot w+L_u u+M\eta=0,\qquad
\ddot w-2\omega\dot u+L_w w=0,
\]

\[
\ddot\eta+L_\eta\eta+M u=0,\qquad
\ddot\zeta+L_D\zeta=0,
\]

where

\[
L_u=-\nabla^2+m_X^2+c\varphi^2+3\lambda_X f^2-\omega^2,
\quad L_w=-\nabla^2+m_X^2+c\varphi^2+\lambda_X f^2-\omega^2,
\]

\[
L_\eta=-\nabla^2+m_\phi^2+3\lambda_\phi\varphi^2+c f^2,
\quad L_D=-\nabla^2+m_D^2+d\varphi^2+\alpha f^2,
\quad M=2c\varphi f.
\]

These equations display the essential difference from a static local Hessian: gradients, \(-\omega^2\), and the gyroscopic \(2\omega\partial_t\) mixing all enter. For modes \(e^{-i\nu t}\), the \((u,w,\eta)\) operator pencil is

\[
\begin{pmatrix}
L_u-\nu^2&-2i\omega\nu&M\\
2i\omega\nu&L_w-\nu^2&0\\
M&0&L_\eta-\nu^2
\end{pmatrix}.
\]

This is a derivation of the starting operator, not a computed spectrum or determinant. In the asymptotic vacuum the two charged rotating-frame branches have frequencies \(\sqrt{k^2+m_{X,\mathrm{vac,tree}}^2}\pm\omega\); they must not both be mistaken for a single neutral static mass.

## Charge, zero modes and the energy definition

With the charge convention used in Candidate 1,

\[
\delta Q=\int d^3x\,[2\omega f u-f\dot w]
\]

at linear order. A fixed-charge projection has to include the momentum/phase contribution; removing only a radial amplitude direction does not specify the full quantum constraint.

The exact background equations give a global-phase zero mode \(w=f\) at \(\nu=0\). The three translation zero modes have \(u=\partial_i f\), \(\eta=\partial_i\varphi\), \(w=0\). A discretization must recover these modes with its stated boundary/volume errors. They require collective-coordinate treatment, not a finite arbitrary eigenvalue inside a logarithm.

One useful formal organization is a charge-constrained semiclassical Routhian \(K(\omega)\), with \(Q=-\partial_\omega K\) and \(E(Q)=K+\omega Q\). If its renormalized one-loop construction and saddle expansion are valid, then at fixed Q the first-order correction can be evaluated from \(K_1\) at the classical conjugate frequency: the first-order displacement of that frequency cancels in the Legendre transform. This identity is conditional on constructing the correct constrained quantum object. It does not authorize replacing \(K_1\) by an unconstrained static determinant or by a sum of unweighted rotating-frame frequencies. A grand-canonical negative direction can occur even for a charge-constrained classical minimum; the integration contour, charge projection and any negative modes must be resolved explicitly.

## Bounded sequence and decision rules

1. **Fluctuation operator.** Derive the quadratic action and charge projection directly, implement radial partial waves and continuum boundary conditions, and recover the phase/translation modes. Verify vacuum rotating-branch thresholds. Check for growing classical modes in the full coupled pencil, with radial/volume/angular convergence. A converged growth mode rejects the proposed classical saddle as the quantum starting state.
2. **Physical thresholds.** Compute the vacuum tadpole and momentum-dependent two-point function in the same matching scheme, obtaining the charged particle pole threshold. Keep total gravitational vacuum matching explicit. The zero-momentum curvatures in this packet cannot supply that threshold.
3. **Renormalized finite-charge correction.** Use a fluctuation determinant or spectral-density/phase-shift method with a specified fixed-Q prescription. Subtract the same vacuum, handle the charge branches and collective modes, and match ultraviolet Born subtractions and counterterms to the already derived scalar beta functions. Demonstrate regulator, partial-wave, radial and volume convergence. Repeat scale/order checks on the finite-charge quantity rather than importing the homogeneous-potential scale residual as its error bar.
4. **Binding and survival decision.** Compare the resulting quantum energy to physical free-particle emission and resolved charged-fragment channels with uncertainties. Report an inconclusive interval if corrections/errors overlap a threshold. A computed energy above an allowed channel rejects binding in that regime. A controlled bound below the tested channels would support a quantum bound-state candidate; it would not by itself determine every tunneling/decay rate or a global lifetime.

Select a finite numerical workload before executing these steps. If the determinant is not controlled, publish that limit with the actual spectrum/threshold evidence; do not add fields, move the couplings, ignore negative eigenvalues, or relabel the local CW integral as the missing result. All empirical and apparatus gaps listed in the candidate remain distinct from this calculation.
