# Flat-space feedback source gate

Start from the strong-drive comparison pinned at
[9ed95fd](https://github.com/ryangrvr/TestingGRUT/commit/9ed95fd285890d922bb2f5f76bc03ed42e765e2a).
The source target is the force density conjugate to the single real
amplitude of \(\theta(t,z)=A(t)\cos(kz)\) for one Majorana field in flat
space. Work density is \(\dot A F\). Use the in-vacuum before a smooth
drive. Gravity, expansion, a bath, and the dynamics of other scalar
Fourier modes are outside this calculation.

Restore the full \(\frac12{\rm Tr}(C\,\partial_A H)\) before renormalizing
the source. Subtract the initial free energy constant in the energy
ledger. The old pointwise free-covariance subtraction is a diagnostic
comparator. It integrates to zero work for the closed pulse and can
change the instantaneous force.

Use the continuum retarded response with \(K(0)=K'(0)=0\), where \(s\)
is the invariant \((\omega+i0)^2-k^2\). This is a finite matching convention
for the scalar kinetic normalization. It uses the same one-Majorana
and spatial-average factors as the preceding weak benchmark. The
primary finite-cutoff method matches the continuum static value
and temporal slope at \(\omega=0\) for this spatial mode, by two
quadratic scalar terms. The independent diagnostic replaces the
whole finite-cutoff linear response with the continuum retarded
one, keeping the nonlinear residual from the same split solver.

## Qualification

- Check the analytic retarded kernel against the independent Feynman
  parameter quadrature on both sides of threshold. Require relative
  agreement \(10^{-9}\) or better.
- Compute linear work on the inherited compact pulse from the time
  response and independently from the two-body spectral integral.
  Require \(0.2\%\) agreement and pre-pulse response below \(0.2\%\)
  of peak. Vary time padding and sampling.
- At amplitude 0.002 on a small grid, compare the complete finite
  force with the tangent of the same split map, including the
  contact term. Require a whole-curve difference below \(0.2\%\)
  of peak after division by amplitude. Check constant-phase
  zero static susceptibility in a separate algebraic control.

## Prescribed source

Use the inherited amplitude-1.2, half-width-12, frequency-2.4,
wavenumber-1.6 compact pulse with \(M=1\).

| Case | Bands to each side | Bloch nodes | Radial nodes per panel | Radial cutoff | Timestep | Boundary |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Base | 5 | 8 | 4 | 5 | 0.04 | Open |
| Time | 5 | 8 | 4 | 5 | 0.02 | Open |
| Spacing | 5 | 12 | 5 | 5 | 0.04 | Open |
| Longitudinal | 8 | 8 | 4 | 5 | 0.04 | Open |
| Transverse | 5 | 8 | 4 | 8 | 0.04 | Open |
| Fine | 8 | 12 | 5 | 8 | 0.02 | Open |
| Tapered | 8 | 12 | 5 | 8 | 0.02 | Tapered |

The tapered regulator changes only the last two links at each end
of the band matrix. Multiply the outermost link by \(1/4\) and the
next link by \(3/4\). All low-band links and the kinetic term agree
with the open regulator. Each regulator has its own vacuum matching.
The grid family retains a finite sharp radial cutoff.

Measure differences on the common time interval, at common
midpoints or by interpolation with the smooth end values.
Require that all one-axis shifts, the base-to-fine shift, the
fine boundary shift, and the fine difference between the two
subtraction methods are below \(3\%\) of the fine primary
force peak in supremum norm. Require primary net work to
agree with free outgoing energy to \(0.2\%\) at the fine grid.
Record raw and matched force histories including both signs.

## Conditional live feedback test

Run only if the prescribed gate passes. Keep the matching convention
fixed, change the spatial wavenumber to 1.3 and use the new forcing
\(J=Z(\ddot A_{\rm ref}+\Omega^2 A_{\rm ref})\) with
\(Z=1\), \(\Omega^2=k^2+1.1^2\), and a compact \(A_{\rm ref}\) of
amplitude 0.85, half-width 8 and carrier frequency 2.7.
Start in the free vacuum with zero scalar amplitude and speed.
The zero-fermion reference is \(A_{\rm ref}\). Evolve to time 12,
including four unforced time units, with the base, time, fine and
tapered grids at the new wavenumber. Use the scalar and fermion
Hamiltonians with the same matching terms in the equations and energy.
Require the refined maximum total-energy minus injected-work defect
below \(0.2\%\) of the peak total energy, and require the effect
on scalar amplitude to exceed three times the sum of the base–fine
and fine–taper trajectory differences. This certifies only the
one-mode flat-space pilot at the stated matching.

If the gate misses, preserve the miss and identify its dominant
cause before proposing a larger simulation. No abundance enters
the matching or the decision.
