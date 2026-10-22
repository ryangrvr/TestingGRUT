# From a parameter derivative to a visible scattering amplitude

Natural units. E means vacuum-relative classical energy of a stationary
fixed-Q Candidate 1 object. v means v_phi, not v_H. The parent action and
its numerical profiles are unchanged. The portal below is proposed only.

## 1. Historical credit, with the operator specified

The established background-field/low-energy-theorem method relates a soft
insertion to a derivative of the mass or amplitude with respect to the
parameter that insertion actually changes. This is the appropriate context
for the parent's single-scale identity dE/dv=E/v. The Higgs/nucleon trace-
anomaly calculation is an important related example, not a GRUT discovery.
See the Shifman--Vainshtein--Zakharov and subsequent sources in SOURCES.md.

However, **varying the Lagrangian parameter v is not automatically the same
operation as adding a physical radial-field fluctuation**. Candidate 1 has
an explicit potential parameter v; its radial field is not a dilaton that
rescales every operator. In particular, write s=phi-v and m_s^2=2 lambda v^2.
The exact classical static equation is

\[
(-\nabla^2+m_s^2)s=-j_s,\qquad
j_s=g_X^2\phi f^2+\lambda(3vs^2+s^3).
\]

Define the zero-spatial-momentum static source J_s(0)=integral j_s. Integrate
the equation with the localized boundary conditions and use the parent's
mass/depletion identity:

\[
J_s(0)=-m_s^2\int s,\qquad
\frac Ev=-m_s^2\int s-\lambda v\int s^2.
\]

Thus

\[
\boxed{J_s(0)=\frac Ev+\lambda v\int s^2.}
\]

This is a precise countercheck on identifying this static source with E/v.
For the three parent solutions J_s(0)/(E/v) is 1.16318, 1.28186, 1.37837.
The additional term is positive and charge-dependent. Nor is J_s(0) by itself
the exterior Yukawa-tail coefficient of an extended object: that also
involves the appropriate spatial form factor. No new quantum vertex theorem
or on-shell radial-mode emission amplitude is asserted here.

The original E/depletion and D first-Born identities are unaffected. Credit
is warranted for the background-derivative method; a universal experimentally
accessible scalar force does not follow without specifying its operator.

## 2. A portal that probes exactly the depletion operator

As an explicitly separate extension, add to the potential

\[
\Delta U=\frac{\kappa}{2}(\phi^2-v^2)
\left(H^\dagger H-\frac{v_H^2}{2}\right),\qquad
H=\frac1{\sqrt2}\binom{0}{v_H+h}.
\]

This preserves the stated scalar sign symmetry and the X charge. The
subtractions are a specified tree-level vacuum/matching convention, not a
solution to vacuum energy or permission to change Candidate 1. Counterterms,
other radiatively generated portals and running remain necessary.

The h-linear perturbation is kappa v_H h(phi^2-v^2)/2. At first order in
kappa, stationarity at fixed Q removes the implicit profile variations from
the energy derivative. The worldline mass perturbation is delta E=g_Qh(0) h,
with

\[
\boxed{g_{Qh}(0)=\frac{\kappa v_H}{2}\int(\phi^2-v^2)
=-\frac{\kappa v_H E}{2\lambda v^2}.}
\]

An alternative route completes the square for a uniform J_H=H^dagger H-v_H^2/2:
v_eff^2=v^2-(kappa/lambda)J_H. The position-independent completion term cancels
in object-minus-vacuum energy. Differentiate E(v_eff)=E(v)*v_eff/v at J_H=0;
the same expression follows. This route varies the correct operator, rather
than assuming a physical s field is interchangeable with v.

An isolated object sources h at O(kappa); its reaction back on phi starts
at O(kappa^2). The displayed source is evaluated on the unperturbed saddle.
Loops, renormalization-group generated H-X/H-D interactions, finite portal
backreaction, gravity and inelastic processes are not included in this
tree-level elastic calculation.

For positive kappa the object has a negative h mass derivative, while an
ordinary nucleon's derivative g_Nh=f_N m_N/v_H is positive. Their scalar
potential is repulsive at this order. The cross section alone loses that
sign. This differs from blindly identifying every scalar charge with mass.

At tree level the h-s mass-mixing entry is kappa v v_H. Away from degeneracy,
theta approximately kappa v v_H/(m_h^2-m_s^2), under the convention that the
SM-like mass eigenstate is h+theta*s at leading order. With SM potential
lambda_H(H^dagger H-v_H^2/2)^2, the sufficient positive-quartic condition
for the two vacuum displacements is kappa^2<4 lambda lambda_H. At the tiny
illustrative coupling all threshold and mixing shifts are negligible; a
finite-coupling completion must rematch measured masses explicitly.

## 3. Finite momentum: structure supplies a shape, not a new fit parameter

At three-momentum transfer q, the leading elastic source is

\[
g_{Qh}(q)=g_{Qh}(0)F_Q(q),\qquad
F_Q(q)=\frac{\int d^3x\,(v^2-\phi_Q^2)e^{i\mathbf q\cdot\mathbf x}}
{\int d^3x\,(v^2-\phi_Q^2)}.
\]

For a radial solution, rho=g_X v r and S=phi/v,

\[
F_Q(q)=\frac{\int d\rho\,\rho^2(1-S^2)j_0(q\rho/(g_Xv))}
{\int d\rho\,\rho^2(1-S^2)},\quad
F_Q(0)=1,\quad F_Q(q)=1-\frac{q^2\langle r^2\rangle_{\rm dep}}6+O(q^4).
\]

The profile is not refitted to detector data. For given couplings, v and Q,
the mass, zero-momentum source and recoil-shape factor are linked. The
overall amplitude kappa and the actual object's preparation/population
remain free. This is falsifiable conditional structure, not a parameter-
free prediction or a signature unique among all composite dark-matter models.

## 4. Nucleon and nuclear recoil normalization

With g_Nh=f_N m_N/v_H, the static potential in Fourier space is
V(q)=-g_Qh(q)g_Nh/(q^2+m_h^2). The nonrelativistic Born amplitude is
f(q)=-mu V(q)/(2 pi). Therefore

\[
\boxed{\sigma_{QN}(0)=\frac{\mu_{QN}^2}{\pi}
\left[\frac{\kappa f_N m_N E}{2\lambda v^2 m_h^2}\right]^2.}
\]

For three different charges at common v,
sigma_QN(0)/(mu_QN^2 E_Q^2) is the same at this order. The laboratory
kinematics, not the relativistic D Klein--Gordon normalization, determine
the reduced-mass factors here.

Using an approximately isoscalar Higgs/nucleus response,

\[
\sigma_{QA}(0)=A^2(\mu_{QA}/\mu_{QN})^2\sigma_{QN}(0),\qquad
\frac{d\sigma}{dE_R}=\frac{m_A\sigma_{QA}(0)}{2\mu_{QA}^2\beta^2}
F_A(q)^2 F_Q(q)^2\left(\frac{m_h^2}{m_h^2+q^2}\right)^2,
\quad q=\sqrt{2m_AE_R}.
\]

F_A is an illustrative Helm nuclear form factor normalized to one. A real
likelihood must use the experiment's nuclear/isotope and response conventions.
The cited f_N=.308(.018) includes the reference treatment of approximately
coherent two-nucleon effects; the full momentum-dependent chiral response
is not computed here. Its uncertainty alone is about 12% in normalization
at linear error order, before halo, detector and model uncertainties.

A sufficient control on repeated Born scattering uses the positive source.
Concentrating all its magnitude at a point gives the upper bound

\[
I_{QA}\le\frac{\mu_{QA}A|g_{Qh}(0)g_{Nh}|}{2\pi m_h}.
\]

This bounds the absolute zero-energy Green-kernel norm; spreading either
positive density cannot enlarge it. For kappa=1e-10 it is at most 2.0e-4
over the listed examples. Scaling to the *1 and 10 GeV* formal one-event
couplings also keeps it below .002. This controls elastic multiple
scattering, not cosmological viability or full quantum corrections.

## 5. Ideal count calculation and its limits

The code integrates a specified truncated Maxwellian halo with v_0=220,
v_E=232 and v_esc=544 km/s. Its lab speed density normalizes to one; analytic
mean inverse speed is checked against direct quadrature. The count is

\[
\frac{dR}{dE_R}=N_A^{\rm target}\frac{f_{\rm halo}\rho_{\rm local}}{E_Q}
\int d^3v\,f_{\rm lab}(\mathbf v)|\mathbf v|
\frac{d\sigma}{dE_R}.
\]

Here N_A^target is nuclei per tonne, not Avogadro's constant, and beta=v/c
in the cross section. Unit conversion, recoil integration of the contact
cross section, raw-spectrum factorization and recoil-grid decimation are
checked separately. The illustration takes f_halo=1, rho_local=.3 GeV/cm^3,
xenon-131, perfect efficiency, no background, and a 5-40 keV recoil window.
These are **not accepted detector events or an LZ sensitivity forecast**.

Counts scale as f_halo*kappa^2 at fixed other inputs. Shape ratios remove
that normalization but retain dependence on v, Q, the halo and nuclear
response. No new observed density or population distribution is derived.
An already-published point-WIMP limit cannot simply be applied when F_Q
changes the recoil spectrum. The detector likelihood has not been run.
