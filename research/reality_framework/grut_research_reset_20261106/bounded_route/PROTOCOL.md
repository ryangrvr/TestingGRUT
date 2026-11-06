# Bounded decision on the physical ALP conversion route

Use the [frozen ALP EFT](../alp_conversion/PROTOCOL.md) and leave H55
feedback paused. Keep the physical radiation hierarchy from the
110 micro-eV comparator and reserve all prospective observations.
The goal is a controlled route to particle conversion.

## Conditional physical preparation and match

Use a constant-mass cosine and the same angular action as before.
At $a=0.49$, set the homogeneous phase to zero and use its
adiabatic rotational speed. Use a Gaussian positive-frequency
vacuum for the angular perturbations at this time. Use
$f_a=2.42927\times10^8$ GeV from the homogeneous normalization
only as an input to the vacuum seed. Earlier linear amplification
must be negligible on the first handoff's energy scale. The
test applies to this vacuum preparation with that phase and
normalization. Primordial and charge-generation perturbations
would require different data.

Propagate the two real canonical conformal quadratures and
the homogeneous phase at the physical $H_*/m_a$. Retain
$2\leq k/m_a\leq3$ in the matched covariance. Use a slowly
swept Lamé calculation on a wider radial grid to check earlier
growth and the power omitted by this band. Use the Gaussian
draws and shared Fourier labels from the prior lattice.
Project the highly amplified band into the real classical
field. The common random quadratures must set both the
field and its canonical momentum. Zero-point modes outside
the band are left outside this classical field. Record
the retained energy, mode participation and a covariance
convergence check with a tighter ODE tolerance.

The two handoffs are when the reference $48^3$, length-24
linear quadrature reaches $10^{-4}$ and $10^{-3}$ in
quadratic fluctuation energy relative to the homogeneous
energy. Before the lattice test require the earlier
vacuum amplification below $10^{-3}$ of the first threshold
and the omitted band below 1% of the estimated amplified
energy at both seams. Require a change of less than 2%
in the matched covariance in a norm weighted by free energy
when tightening the ODE tolerance at the same times.

## Fixed short nonlinear window

The common endpoint is 1,500 conformal mass units after
the first seam. Allow five cases with the listed maximum
steps and a six-minute wall cap per case. A case that hits
the cap has a missing endpoint and fails qualification.

| Case | Handoff | Sites per side | Length | Maximum conformal step |
| --- | --- | ---: | ---: | ---: |
| Base | First | 48 | 24 | .05 |
| Spacing | First | 64 | 24 | .05 |
| Volume | First | 64 | 32 | .05 |
| Time | First | 48 | 24 | .025 |
| Later seam | Second | 48 | 24 | .05 |

Retain the physical expansion during the window. Track
comoving energy and pressure work, coherent kinetic energy
and the variance and gradient energy of nonzero modes.
Use the last 100 conformal units to test a handoff to
free-particle transport. Require a fixed global vacuum
branch through that window, less than 2% anharmonic
energy, fewer than $10^{-4}$ of sites outside that
branch, less than 1% of the quadratic count in the
cutoff tail, and less than 2% drift and dictionary
disagreement. Require a comoving energy-budget residual
below 1% and terminal counts that agree within 5%
on the seam, timestep, spacing and volume comparisons.

If these conditions hold, the short window establishes
a route for a separate redshift and dilute-transfer test.
A failure preserves the matched states, measures the
blocker and ends this ALP lane under the resource boundary.
No further extension or detector comparison is part of
this test.

## Conditional replacement

If the ALP route misses, open one bounded independent
comparator. Use equilibrium thermal freeze-out of a
stable real scalar $S$ with $m_S=40$ GeV and
$V\supset\lambda S^2H^\dagger H/2$ as its only portal.
Assume radiation domination and all dark matter in $S$.
Use the thermal abundance to calibrate the portal or
test a conservative overabundance bound at the collider
ceiling. Use published Higgs invisible data as the
independent observable, with no halo assumption.
Select the mass before calculating either side. Preserve
an exclusion of this fixed point if the two requirements
cannot be satisfied together.
