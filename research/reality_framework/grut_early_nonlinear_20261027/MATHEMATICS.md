# Preparation, reduced dynamics and measured quantities

Natural units hbar=c=1 are supplied. All quantities here use H_i units, rather than the late cold pilot's m_X units. Define chi_phys=(sqrt2 Re(aX),sqrt2 Im(aX),a phi), code fields chi=chi_phys/H_i and dimensionless conformal time s=H_i eta. The code momentum is pi=dchi/ds=chi_phys'/H_i², where the prime in the last expression differentiates dimensional conformal time.

## Finite-band initial data

At tau=0, diagonalize the six-real-field conformal mass matrix M0 from the audited full X/Y/Sigma model. On the retained band1.8<=k/H_i<=8, every k²+eig(M0) is positive. The specified instantaneous conformal Gaussian has

\[
C_{qq}=U\,\mathrm{diag}(1/(2\omega_j))U^T,\qquad
C_{pp}=U\,\mathrm{diag}(\omega_j/2)U^T,\qquad C_{qp}=0,
\quad\omega_j^2=k^2+\lambda_j(M_0).
\]

The independent phi mode starts with the same instantaneous positive-frequency convention, with its actual g_X²|chi_X|²−a″/a curvature. Propagate all these covariances to tau20 with the full scalar transfer and the phi transfer. e=0 defines the leading scalar limit used for matching and evolution. The positive spectral square roots, all12 charged columns and both phi columns are preserved in MATCH_L12.npz and MATCH_L16.npz. The nonlinear stage uses the light marginal; heavy contributions to that marginal are retained, while later heavy dynamics are not.

For each retained Fourier pair±n, draw independent complex standard Gaussians for all14 initial coordinates/momenta, map them through those square roots and transfers, retain X/phi, and impose complex conjugacy for the negative wavevector. Forward Fourier coefficients carry1/sqrt(V). Thus each initially specified mode has the intended continuum covariance. A stable seed keyed to(seed,L,nx,ny,nz) makes the continuous finite Fourier realization identical across spatial grids. Continuum momentum determines preparation; the FD4 symbol determines subsequent lattice dynamics.

No noise energy or sample charge is retuned. A given realization's global light current differs slightly from the homogeneous reference and is reported. This classical finite-band ensemble is not a selected full quantum state: its missing IR/UV completion, Gauss projection, quantum counterterms and later quantum evolution are absent. In particular, sampling a Gaussian with vacuum-like half quanta and evolving it classically does not itself establish a renormalized quantum particle population.

## Leading early radiation Hamiltonian

After matching, use

\[
\mathcal H=\tfrac12\sum_A\pi_A^2+\tfrac12\sum_A|\nabla\chi_A|^2
+\tfrac{\lambda_X}{4}(\chi_1^2+\chi_2^2)^2
+\tfrac{g_X^2}{2}\chi_\phi^2(\chi_1^2+\chi_2^2)
+\tfrac{\lambda_\phi}{4}\chi_\phi^4,
\]

with lambda_X=.001, g_X²=.01, lambda_phi=.03. The neutral vacuum scale is negligible in H_i units at these times. This massless early model has no late false-vacuum energy reservoir added to it. It omits the heavy-mediated corrections, gauge evolution, residual curvature, the late phi bias and actual decay/collision sources. Their absence is part of the model statement, not a full-EFT certification.

The full/reduced comparison on tau20–50 measures a0.547% orbit-scaled background difference, about0.98–1.03% preparation-phase tangent difference and1.48–2.72% neutral-transfer difference. The charged force differs by0.206% at matching. These are linear/homogeneous diagnostics; they do not bound all nonlinear heavy-field errors. A frozen-Sigma, static quadratic-Y estimate is also reported:

\[
\left|\Delta F_X/F_{\lambda_X}\right|
\lesssim {3\kappa^2\,\max|\chi_X|^2\over
4\lambda_X (M_Y^2-\sqrt2\mu h)\,a^2}.
\]

This follows by bounding the eliminated quadratic Y response with its smallest eigenmass. It assumes Sigma fixed and does not include later Sigma/gauge dynamics. Its measured maximum is about0.45% in the pilot. It is an approximation diagnostic, not an integrated error theorem.

The lattice uses the same symmetric fourth-order Laplacian and matched energy as8159556. Yoshida's fourth-order kick/drift composition preserves the active global X current exactly up to roundoff. The radiation clock is

\[
a(s)=a_m+Cs,\qquad \tau(s)=20+a_ms+\tfrac12 Cs^2,
\]

with s the elapsed conformal time, and a_m,C matched at20. Each integration step stops at the next physical output time, so the actual bulk records satisfy Delta tau<=.05; the half-step control uses.025. No sparse-grid interpolation generates those bulk measurements.

## Active energy and exchange

Book the interaction in E_X. Then

\[
E_\phi=K_\phi+G_\phi+U_\phi,\quad
E_X=E-E_\phi,\quad {dE_X\over ds}=\int g_X^2\chi_\phi\pi_\phi|\chi_X|^2d^3x.
\]

The last equation is reversible Hamiltonian work. Its opposite appears in dE_phi/ds. For comparison with a neutral linear-oscillator control we separately report E_phi+E_mix; this moves the complete interaction to that diagnostic. The two conventions are not added together. All are raw classical field energies, with no quantum vacuum subtraction.

The spatial mean is the exact Fourier zero mode. Its classical charge is Q0=V(bar chi1 bar pi2−bar chi2 bar pi1); the nonzero-mode charge is Q−Q0. Q0/Q can exceed1 when the fluctuation sector carries negative net current. It is not a retained-object fraction. Evaluate the same Hamiltonian on the spatial means to define E0. E−E0 includes nonlinear interaction/correlation energy and is not automatically a free-particle bath.

At fixed total current density the bare quartic model has a homogeneous circular lower bound E_min/V=(3/4)lambda_X^(1/3)|n_Q|^(4/3). This gives an energy budget, not a prediction that it is reached. On infinite space the fixed-charge scaling chi_R(x)=chi(x/R)/R, pi_R(x)=pi(x/R)/R² keeps Q fixed and sends E→E/R. Thus this leading massless vacuum model has no absolutely stable localized positive-energy minimum at fixed Q. Transient structures and later finite-density/binding physics require their own analysis; this does not prove an absence theorem for the full model.

## Export and particle limit

SPECTRAL_EXPORT.npz stores binned real symmetric classical covariances of(chi1,chi2,pi1,pi2,chi_phi,pi_phi), separately weighted with k_eff², with k0 removed and the means retained separately. Kinetic energy, exact discrete gradient energy and current reconstruct from these matrices. The full potential stays outside momentum bins.

For reference only, a common chosen omega_k gives

\[
n_\pm(k)=\tfrac12\left[
{ |\pi_1|^2+|\pi_2|^2+\omega_k^2(|\chi_1|^2+|\chi_2|^2)\over2\omega_k}
\pm\operatorname{Re}(\chi_1^*\pi_2-\chi_2^*\pi_1)\right].
\]

We compare omega²=k_eff² and omega²=k_eff²+2lambda_X<|chi_X|²>+g_X²<chi_phi²>. Both reproduce the same nonzero-mode charge; their sums differ substantially. Neither contains a measured retarded pole, an assured quasiparticle basis, a quantum zero-point subtraction or a complete kinetic distribution. The physical nbar, f_sym and threshold remain null.
