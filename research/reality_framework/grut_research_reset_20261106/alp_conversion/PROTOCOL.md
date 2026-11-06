# Frozen ALP charge conversion test

Retain the [paused H55 feedback decision](../feedback_source/REPORT.md).
Use the fixed 110 micro-eV retrospective comparator from that packet.
Keep the haloscope likelihood out of all source, state, and numerical
choices. Its exact-frequency retrospective comparison is separate
and follows the conversion decision.

## EFT and physical preparation

Use one real angular field with kinetic coefficient $f_a^2$ and
$V=f_a^2m_a^2(1-\cos\theta)$, constant $m_a$ and $f_a$.
At $\mu=T_{\rm EW}=130$ GeV set $c_W=-1$, $c_Y=c_g=c_f=0$
in the Co–Hall–Harigaya EFT convention. The leading values
are $c_B=12/79$ and $c_\gamma=-1$. Use $Y_B=8.7\times10^{-11}$
and $d_a=0.44$ eV for the all-ALP abundance comparator. Estimate
the size of radiative matching effects separately without changing
the anomaly texture. Freeze the late scalar dynamics to the cosine
self-interaction and a prescribed radiation-dominated background.
Neglect further charge exchange with the SM, gravitational
perturbations in the evolution, and changes in radiation comoving
entropy. The mass torque itself changes the rotation charge.

Use $a=1$ at the adiabatic homogeneous separatrix, with the early
dimensionless comoving charge $q_\theta=4/\pi$. This is an
action-based normalization of the early rotating state. It fixes
$s_*=(4/\pi)c_B T_{\rm EW}^2m_a/Y_B$. Use idealized $g_*=g_{*s}=61.75$
near this epoch and 106.75 at electroweak freeze-out. The physical
$H_*/m_a$ is a derived scale. In the specified classical preparation
the field equations and the early angular speed are independent of
$f_a$. Convert the measured $C$ to $f_a$ and $Y_\theta$ only at
the end. If a fluctuation or interaction changes with $f_a$,
the matched calculation must solve for that dependence jointly.

The physical linear diagnostic starts at $a=0.25$ on the rotational
adiabatic trajectory. Compare two specified subhorizon seeds, a
Gaussian phase-mixed state as defined below and an adiabatic vacuum
with canonical normalization $f_a$. Use $f_a$ corresponding to
$C=0.5,2,3$ for vacuum onset sensitivity only. These are inputs to a
diagnostic, with no claim that they bound physical $C$. Record where
an unbackreacted, slowly swept Floquet calculation first puts one
percent of homogeneous energy into fluctuations. Stop the linear
inference at that point. A transfer from primordial curvature to
these modes requires a charge-generation and perturbation model.

## Numerical ladder and particle definition

First check homogeneous trapping with the pendulum action and
independent radiation-era integration. Require the slow homogeneous
control at $H_*/m=0.001$ to recover $C=2$ to 3%, with small late
amplitude. Check the rotating and librating Lamé bands with monodromy
integration and require an inside-band Floquet rate to agree with
a small-amplitude 3D plane-wave control to 2%. Check stable modes.

In the nonlinear ladder use periodic 3D boxes, initial $a_i=0.6$,
phase zero at the bottom of the potential, rotation set by the
homogeneous action, and $a_f=4$. The Gaussian ensemble has zero
perturbation mean and independent position and velocity quadratures
with free equipartition at the initial phase. The per-logarithmic
velocity variance is $A_s\dot\theta_i^2 W(k)$, where $A_s=2.1\times10^{-9}$,
$W(k)=\exp[-(k/4)^4-(0.30/k)^4]$, and comoving $k$ is in units of $m_a$.
It is a controlled phase-mixed classical state with a fractional
velocity spectrum normalized to the size of primordial curvature.
It is a prescribed transfer at $a_i$, not the actual high-hierarchy
cosmological transfer. The lower seed scales the RMS by 0.1.
Use shared Gaussian mode labels across box resolutions and report
changes across independent seeds.

| Run | $H_*/m$ | Sites per side | Box length in $m^{-1}$ | Conformal step in $m^{-1}$ | Seed | RMS multiplier |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Base | .01 | 48 | 32 | .05 | 314 | 1 |
| Spacing | .01 | 64 | 32 | .05 | 314 | 1 |
| Volume | .01 | 72 | 48 | .05 | 314 | 1 |
| Time | .01 | 48 | 32 | .025 | 314 | 1 |
| Seed | .01 | 48 | 32 | .05 | 2718 | 1 |
| Lower seed | .01 | 48 | 32 | .05 | 314 | .1 |
| Faster expansion | .02 | 48 | 32 | .05 | 314 | 1 |
| Slower expansion | .005 | 48 | 32 | .05 | 314 | 1 |

Use the radiation conformal variable $\chi=a\theta$, a projected
sine force and a time-reversible second-order integrator.
Track comoving energy and pressure work. Only assign particles
on a single global vacuum branch, with a negligible anharmonic
energy fraction and adiabatic free frequencies. Report the comoving
occupation spectrum, momentum moments, fraction at zero momentum,
ultraviolet tail and the uncertainty across grids and seeds.
The primary occupation is the adiabatic free-field action. Compare
against the instantaneous free-energy estimator.

For a conditional conversion on a run require the entire $3.5\le a\le4$
window to have less than 2% action drift and less than 2% disagreement
between dictionaries. Require less than 2% anharmonic energy, fewer
than $10^{-4}$ of sites outside the single branch, less than 1% number
in the outer 20% of the retained component wavenumber range, and
comoving energy residual less than 1% of late energy. Use a 2%
spacing, volume, and timestep tolerance on the base conversion.
Treat seed, lower seed and hierarchy as physical sensitivities.
If these criteria miss, retain the diagnostics and the miss.
Do not extrapolate a reduced hierarchy to the physical point
without a matched source and controlled hierarchy error.
