# What the soft receiver must distinguish

All quantities below use dimensionless distance rho=m_X r, time tau=m_X t
and momentum k=p/m_X. This packet does not alter the scalar action.

## The existing numerical layer

For a vacuum linear charged field, U=rho X removes the spherical 2/rho
derivative. The nodal production solver alternates its Hamiltonian evolution
with a half-step velocity multiplier exp[-dt eta(rho)/2] at each end of the
step. Its continuum limit in the layer is

\[
U_{tt}=U_{xx}-U-\eta(x)U_t,\qquad
\eta(x)=\eta_{\max}(x/D)^3,\quad 0<x<D.
\]

The outer nodal endpoint is Dirichlet. Here x measures distance from the
layer entrance, and D=L-rho_abs. Original N runs have D=90 and eta_max=.03;
the old larger box has D=120. This benchmark is the linear continuum limit,
not a full nonlinear or split-time finite-grid scattering calculation.

For U=e^{-i omega t} f(x), omega=sqrt(1+k^2),

\[
f''+[k^2+i\omega\eta(x)]f=0,\qquad f(D)=0.
\]

In the physical region f=e^{ikx}+R e^{-ikx}. Backward integration of the
homogeneous layer solution with f(D)=0, f'(D)=1 gives

\[
A=\frac12\left[f(0)+\frac{f'(0)}{ik}\right],\quad
B=\frac12\left[f(0)-\frac{f'(0)}{ik}\right],\quad R=B/A.
\]

Incident and reflected waves share omega and k, so their energy ratio is
|R|^2. This includes both entrance reflection and the attenuated wall return;
these amplitudes interfere and cannot generally be added as probabilities.
No unknown collapse population is weighted by these test-wave results.

The current j=Im(f* f') obeys j'=-omega eta |f|^2. With incident amplitude
one and the reflecting wall,

\[
1-|R|^2=\frac{\omega}{k}\int_0^D\eta(x)|f(x)/A|^2 dx.
\]

The positive loss integral independently checks normalization and the sign
of absorption. Zero damping gives R=-exp(2ikD). A constant-damping slab has
q=sqrt(k^2+i omega eta), Z=f'(0)/f(0)=-q cot(qD) and
R=(ik-Z)/(ik+Z), providing another closed-form control.

A separate second-order finite-difference boundary-value solve uses
f'(0)+ik f(0)=2ik and f(D)=0. Six spacings check convergence to the ODE
result; it is not a rerun of the same initial-value solver. Neither scheme
is presented as the full production time-stepper's exact discrete boundary.

## The trace audit

All six traces are newly loaded and matched to their parent's float64
hashes. No new collapse runs were executed. The nodal probe is exactly
round(rho_probe/h) h, unlike the earlier cell-centered solver's radius.

The measured energy and charge currents at that surface are

\[
J_E=-4\pi\rho^2\operatorname{Re}(\dot X^*\partial_\rho X),\quad
J_Q=4\pi\rho^2\operatorname{Im}(X^*\partial_\rho X).
\]

Their window-squared time integrals are compared with the discrete Fourier
Parseval sums made from the independently recorded field, radial derivative
and momentum. These are signed measured currents, including returning waves.

Conditional outgoing/incoming amplitudes assume a vacuum harmonic spherical
wave, X=(a_out exp(ik rho)+a_in exp(-ik rho))/rho. With signed k matching
the signed frequency, they are proportional to

\[
a_{out,in}=\frac12\left[F\pm\frac{D+F/\rho}{ik}\right].
\]

The existing local control requires k^2>10 max|S^2-1+.1|X|^2|, as well
as k^2>1e-4. Passing that cut alone does not remove finite-time windows,
returning radiation, coherent nonlinear evolution, or missing exterior fields.
For a 12000-unit record, the ordinary Fourier frequency spacing is about
5.24e-4. At k=.01 the vacuum frequency excess sqrt(1+k^2)-1 is only 5.00e-5.
The present FFT method therefore lacks comparable near-threshold frequency
resolution; its characteristic required duration is about 1.26e5, before
other source/window errors. Merely refining the spatial cell size cannot
fix that time-window limitation.
All-band pair numbers remain conditional diagnostics. The difference
between energy computed from measured momentum and the vacuum harmonic
out-minus-in construction is retained, not normalized away.

N0125_box and N0125_far have bitwise identical final fields and physical
histories; only their extraction radius changes. Their matched traces
provide a radius control without changing the underlying realization.

## A matched-layer design reference, not an implemented fix

Complex coordinate stretching is established numerical analysis, not a new
GRUT interaction. In the homogeneous linear exterior, take

\[
\partial_x\longrightarrow s^{-1}\partial_x,\quad
s=1+i\sigma(x)/\omega,\quad \sigma=\sigma_{max}(x/D)^3.
\]

The continuum outgoing solution is exp[ik integral s dx], with no entrance
reflection. A Dirichlet wall returns energy

\[
|R_{wall}|^2=\exp\left[-4\frac{k}{\omega}\int_0^D\sigma(x)dx\right]
=\exp[-(k/\omega)\sigma_{max}D].
\]

For D=90 and sigma_max=40 this is below 2.4e-16 at k=.01. It establishes
a possible continuum design target. A real implementation still needs a
causal time-domain realization, initial-memory convention, discrete stability
and reflection controls, plus a conserving field-to-source interface. No new
time-domain PML is installed or claimed here. Nonlinear exterior forces also
require direct control before extending a vacuum match to formation fields.

The unresolved soft source stays unresolved. No reflection deconvolution,
particle completeness, cosmological production epoch or new abundance is
inferred from the boundary benchmark.
