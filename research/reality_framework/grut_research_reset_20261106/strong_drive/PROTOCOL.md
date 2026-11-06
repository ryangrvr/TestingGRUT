# Frozen strong-drive propagation comparison

The pulse and decision rule are inherited from the [research reset](../AUDIT.md)
at `338d56459279294b79cba2cdb3133087f247c1ee`. Use a prescribed mass phase
`A(t) cos(kz)` with `A(t) = 1.2 cos^4(pi t / 24) cos(2.4 t)` on
`-12 < t < 12`, zero outside, and `k = 1.6`. Units set the fermion mass to
one. Both arms use the free in-vacuum and free outgoing number and energy.
Transverse phase space is continuous and integrated radially. The local
arm treats each phase of the spatial wave as an independent homogeneous
problem, with the same momentum integration and spatial average. The external
drive is common. This experiment excludes gravity, bath, expansion and
backreaction by design.

The propagating arm is periodic in space with a continuous Bloch momentum
and an odd number of longitudinal bands. Use the full four-component spinor.
At finite band cutoff take the exponential of the projected cosine in the
mass term. Use a cyclic momentum regulator at the combined refined grid to
measure sensitivity to the boundary prescription. Both finite mass matrices
are unitary before subtraction of the free mass. Time evolution uses exact
free and exact residual substeps in a symmetric split. Compute work from
the explicit time derivative of the Hamiltonian against the evolving vacuum,
with the free covariance subtracted. This constant subtraction integrates to
zero for the closed pulse. Use half the doubled spinor trace for one Majorana
field. Record the energy both from this work integral and from the free
outgoing particles and holes. Number and work are the decision metrics. The
outgoing free-energy moments are diagnostic.

## Gate before the reserved run

- The tangent to the propagation integrator at zero drive must reproduce all
  four archived weak Gaussian energy densities and their cell-local values
  to 0.2% relative at the finer grid. Report the coarse-to-fine differences
  and use numerical midpoint Fourier sums at the finer grid for the
  propagating response. Require a static vacuum with zero resolvable
  production and agreement for a spatially constant drive within roundoff.
- For a separate, weak compact pulse at amplitude 0.002, compare the
  normalized finite-amplitude production with the solver's tangent on
  the same finite momentum grid. Require 0.5% agreement on work-energy
  density. A smaller amplitude may be used if explicitly shown to improve
  this tangent check.
- Independently compare one propagating sector with an unsplit DOP853
  integration. Use a separate compact pulse of amplitude 0.7 with
  half-width 5, frequency 2.3, wavenumber 1.4 and two sidebands on each
  side. Compare number and energy to 0.2% at timestep 0.01. Use exact
  spinor evolution in the homogeneous and vacuum controls.

## Registered strong grids

`L` is the number of bands to each side of the central band.
Use Gauss–Legendre quadrature on the first Brillouin zone and on transverse
panels with edges 0, 0.5, 1, 2, 3, 4, 5, 6 and 8. Local spatial phases are
uniform midpoint samples in [0, pi].

| Case | L | Bloch nodes | Nodes per transverse panel | Transverse cutoff | Spatial phases | Timestep | Momentum boundary |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Base | 4 | 6 | 3 | 4 | 16 | 0.04 | Open |
| Time | 4 | 6 | 3 | 4 | 16 | 0.02 | Open |
| Spacing | 4 | 12 | 5 | 4 | 32 | 0.04 | Open |
| Longitudinal | 7 | 6 | 3 | 4 | 16 | 0.04 | Open |
| Transverse | 4 | 6 | 3 | 6 | 16 | 0.04 | Open |
| Fine | 7 | 12 | 5 | 6 | 32 | 0.02 | Open |
| Regulator | 7 | 12 | 5 | 6 | 32 | 0.02 | Cyclic |

For each observable, estimate a conservative measured numerical allowance by
adding the absolute one-axis shifts from the base in both arms, the absolute
crossed-grid remainder in both arms, and the change under the regulator
switch. For work add the fine independent energy-budget residuals in both
arms. This is a measured refinement envelope. Check work against outgoing
energy in each fine arm and require a relative residual below 0.2%. The
fine work/number comparison is decisive for a large difference if its
absolute gap exceeds 1% of the larger value and three times the measured
allowance for both observables. A small difference requires the gap plus
three times the allowance below 1% for both. Other cases stay unresolved.
Any further runs get separate labels and preserve the original comparison.
