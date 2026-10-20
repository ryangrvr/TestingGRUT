# Finite objects, symmetries and an empirical exclusion

This is an additive analysis of the unchanged Candidate-1 leading scalar action,
not a new action chosen to match a target. Natural units and a=1,b=.1 are used.
The simulation evolves u,w,S with all three spatial coordinates; D=0, flat local
gravity. Quantum and strong-gravity claims require further work.

## What unrestricted perturbations mean here

The Cartesian runs impose no azimuthal or reflection symmetry. The generic seed
uses an obliquely oriented, volume-preserving ellipsoidal deformation, a spatially
varying phase, and independent radial and S velocity kicks. Pair seeds also have
opposite transverse offsets and a velocity-direction perturbation. Both signs of
the seed are retained. The frozen protocol defines their exact amplitudes and
the finite observation interval. No finite campaign can certify survival after
arbitrary energy injection, every initial condition, or for infinite time.

With V=(S^2-1)^2/4+S^2 F^2/2+b F^4/4 and F^2=u^2+w^2,

\[
\ddot u=\nabla^2u-(S^2+bF^2)u,\quad
\ddot w=\nabla^2w-(S^2+bF^2)w,\quad
\ddot S=\nabla^2S-(S^2-1+F^2)S.
\]

The code uses the Hamiltonian associated with its nearest-neighbour differences.
All positive-direction bonds are counted once, including every fixed-vacuum
face bond. Kick-drift-kick conserves the discrete U(1) charge to roundoff and the
Hamiltonian to its time integration error. Boundary reflection is real in this
conservative finite-box problem, so conservation alone cannot exclude it.

## Symmetry tangents are moving solutions

For a solution B=(u,w,S,pu,pw,pS), apply generators to **all six canonical
components**. Exact continuum tangents include

\[
T_x B=-\partial_xB,\qquad R_yB=(z\partial_x-x\partial_z)B.
\]

In an axisymmetric background their real m=1 coefficients are -partial_rho B
and z partial_rho B-rho partial_z B. Reversing a generator's overall sign does
not change its interpretation. The full rotation rotates positions, field
patterns and velocities. Opposite source displacements with velocities held
fixed change the encounter's impact geometry and are not that symmetry.

Their canonical positive norm is a diagnostic, not a conserved Noether charge:

\[
N_\delta^2=\int d^3x\sum_A[(\delta q_A)^2+(\delta p_A)^2].
\]

A time-dependent encounter can change N_delta even for a pure symmetry tangent.
The deciding control is the evolving tangent's residual against the generator
of the contemporaneous background, plus spatial/time/boundary convergence.
Cartesian differences and fixed boundaries break continuous spatial symmetries.
Their residuals must be measured; a nonvanishing late residual is not evidence
of a physical instability. The separate large-domain script audits the original
linear m=1 discretization, including local core displacement fits. It does not
replace the unrestricted nonlinear evolutions.

## The isolated object's fixed-charge energy test

For a stationary nodeless object X proportional to F(r) exp(-i Omega tau), set
I=integral F^2 d^3x and Qbar=Omega I. At fixed charge, eliminating canonical
momenta by Cauchy-Schwarz gives the shape energy

\[
E_Q[F,S]=\frac{Q^2}{2I}+\int d^3x\,[|\nabla F|^2/2+|\nabla S|^2/2+V].
\]

For a spherical harmonic with unit angular norm, transform the radial
perturbation to v=r delta F, s=r delta S. The ordinary imaginary/phase block is

\[
H^-_\ell=-\partial_r^2+\ell(\ell+1)/r^2+S^2+bF^2-\Omega^2.
\]

The real amplitude/S block is

\[
H^+_\ell=\begin{pmatrix}
-\partial_r^2+\ell(\ell+1)/r^2+S^2+3bF^2-\Omega^2&2SF\\
2SF&-\partial_r^2+\ell(\ell+1)/r^2+3S^2-1+F^2
\end{pmatrix}.
\]

For ell=0, fixed Q adds the positive rank-one term

\[
H^+_{0,Q}=H^+_0+
\frac{4\Omega^2}{\int r^2F^2dr}\,|rF,0\rangle\langle rF,0|.
\]

The radial L2 inner product contains dr; the discrete rank-one matrix includes
that quadrature weight. Omitting the fixed-charge term would test the wrong
constraint: the unconstrained ell=0 block does have one negative direction.

The phase zero vector is rF. Translation gives the ell=1 zero vector
(rF',rS'). For a strictly decreasing F and increasing S, flipping the sign of
the amplitude component makes this vector positive and the off-diagonal
coupling nonpositive. The ground-state/nodal argument identifies the translation
as the bottom of that coupled sector. Furthermore

\[
H^+_\ell=H^+_1+[\ell(\ell+1)-2]/r^2\quad(\ell\ge2),\qquad
H^-_\ell=H^-_0+\ell(\ell+1)/r^2.
\]

Thus higher angular harmonics cannot insert a lower negative direction. The
remaining numerical question is the fixed-charge radial block and the correct
symmetry limits. `spectral.py` re-solves the profile, uses dr=.1,.05,.025 and an
enlarged outer radius, and records eigenpair residuals. Small negative sampled
phase/translation eigenvalues must approach zero as dr^2; they must not be
reported as a physical tachyon. The continuum vacuum threshold in the X block
is 1-Omega^2=.36. Finite-box scattering levels slightly above .36 move toward
it as the box grows; they are not new discrete internal resonances.

These computations support a local constrained-energy minimum of the isolated
reference branch, modulo phase and translation. They are not rigorous interval
spectral enclosures, a global-minimizer proof, a claim about every Q, or a
certificate for a time-dependent colliding pair. At D=0 the additional D
quadratic block is positive with mass squared proportional to 1+.01 S^2 at the
stated benchmark; finite D excitation and metric evolution were not simulated.

## A conserved-charge chain that can be excluded without a mass function

For a,b>0, direct minimization gives

\[
\inf_{F>0,S}\frac{2V}{F^2}=\min(1,\sqrt{ab}),\quad
\mu_*=(ab)^{1/4}=0.5623413252\quad(ab<1).
\]

Completing the kinetic squares then proves the leading classical bound
Ebar >= mu_star |Qbar|. Physical energy per charge has one extra m_X factor:
E_phys >= mu_star m_X |Q_phys|. This also holds on an instantaneous spatially
flat FRW slice: scale-factor factors in the gradient energy remain positive.

One explicit additional preparation uses g_X=.1,v=1000 GeV,m_X=100 GeV and a
primordial circular nonthermal condensate F=.5,S=sqrt(.75),Omega=sqrt(.775),
with fully 3D perturbations specified in the production protocol. Its charge
yield is **postulated** to equal the late baryon yield inferred independently
from D+3nu BBN. No hidden charge generation from an initially uncharged vacuum,
no baryogenesis mechanism and no reason for the equality are derived.

Exact adopted U(1) and conserved visible entropy imply Y_Q=n_Q/s is constant.
Specific composites may fission or evaporate into charged X; neither process
erases net Q. Consequently their unknown production efficiency and individual
quantum lifetime cannot remove the energy floor

\[
\omega_{\rm hidden}\equiv\Omega_{\rm hidden}h^2
\ge \frac{\mu_*m_XY_Qs_0}{\rho_c/h^2}
=\frac{\mu_*m_X}{\bar m_b}\,\omega_b
\qquad(Y_Q=Y_B).
\]

BBN omega_b=.02233+-.00036 and the stated .938 GeV mean-baryon-mass convention
give a floor 1.3387. A 1% mass-conversion allowance is retained. Initial
H/m_X is about 4.33e-10; omitting expansion during tau240 is consistent at
H delta t about 1.04e-7. Cold hidden initial conditions are compatible with a
much hotter visible plasma because this benchmark has no leading direct portal.

If the late relic is cold collisionless matter, this exceeds the held-out
Planck omega_c=.1200+-.0012 by at least a factor about 11. That comparison is
conditional on the late matter interpretation; a full CMB likelihood is absent.

A less restrictive observable follows directly from the flat-FRW constraint
with nonnegative other densities (rho_Lambda=0 is the most conservative point):

\[
H_0\ge100\sqrt{\omega_b+\omega_{\rm hidden,minimum}}
\ \mathrm{km\,s^{-1}\,Mpc^{-1}}.
\]

No late hidden equation of state enters this instantaneous constraint. It gives
H0>=116.66; using a three-sigma-low BBN density and the mass-conversion allowance
still gives H0>=113.25. Independent geometric megamaser measurements give
73.9+-3.0, with three-sigma upper edge 82.9. This excludes this flat,
nonnegative-energy, conserved-yield preparation at the adopted leading EFT
order. The published values were held out from parameter fitting, not hidden
from the analyst, and no new raw maser likelihood has been evaluated.

Changing the primordial yield, scale or cosmology is a new hypothesis. The
action alone leaves them free; this exclusion is not evidence against every
member of the scalar family and is not a successful observation of composites.

Primary calibration/comparison sources and exact adopted values are retained
in the three data JSON files. The BBN estimate has no CMB baryon prior. The
megamaser observable uses ordinary geometric distances and velocities, not the
hypothetical D reporter or a claimed direct hidden-scalar detector.
