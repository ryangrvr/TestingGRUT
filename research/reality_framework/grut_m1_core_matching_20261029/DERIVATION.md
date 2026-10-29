# Derivations and interpretation boundaries

All times, radii, energies and momenta below use `m_X=1` unless stated.

## Energy accounting for the newly observed source

For the inherited homogeneous kinetic prescription, the initial residual is

    R = E_pulse - M_stationary(Q_initial) - E_initial_gas.

The later radial source changes both `Q_initial` and `E_initial_gas`.
Therefore

    R_new - R_old = -[M_stationary(Q_new)-M_stationary(Q_old)]
                    -[E_gas,new - E_gas,old].

This is why subtracting the late emitted energy alone from the old residual
is incorrect. At the final scale factor, the reported modeled total is

    [M_stationary(Q) + N_+ + N_- + K_gas + R] / E_ref - 1.

The fine 96-node components are recorded in `SCENARIO_AUDIT.json` and add
to 0.8944364%. The tiny signed ODE residual for `N_-` is a solver error,
not a physical negative population. The excitation is not dynamically
coupled to the gas in this model.

More generally, for a closed sector in an FRW background,
`d(a^3 rho)/d ln a = -3 a^3 P`. If the entire actual sector energy is
included and the integrated pressure work is nonnegative, comoving energy
is nonincreasing. If the cold rest/soliton energy is bounded above by that
total, that total supplies a conditional bound on later cold mass. Internal
scattering and catalytic core conversion rearrange the energy budget rather
than create new energy. The completeness, stress/pressure and dynamical
conditions are not demonstrated by the current source prescription.

## Radial shell conservation

The parent finite-volume Hamiltonian assigns each cell's gradient energy to
its outer face. Applying its own volume, potential and gradient conventions
to `40<=r<80` gives the change of stored shell energy and charge. The
separately recorded probe currents are evaluated with centered spatial
derivatives at the nearest cell centers. Thus the shell ledger and probe
flux should agree up to matching and discretization errors, not identically.
`MATCHING_SURFACES.json` records both sides; their agreement tests net
conservation but does not bound energy-conserving reactions.

## Continuous-energy Poisson cloud

For a source flux `F(t,E)` per energy emitted at `r0`, radial free streaming
gives

    n(r,T) = [1/(4 pi r^2)] integral dE F(T-(r-r0)/v,E) / v.

To avoid representing a finite energy bin as one velocity, change variables
from `E` to emission time `u=T-(r-r0)/v(E)`. With `v=(r-r0)/(T-u)` and
`gamma=(1-v^2)^(-1/2)`,

    dE/du / v = gamma^3 v/(T-u),
    n(r,T) = [1/(4 pi r^2)] integral du F(u,E(v)) gamma^3 v/(T-u).

The code uses piecewise-linear interpolation in energy and flat half-bin
endcaps; those are assumptions. For independent Poisson centers of density
`n_c`, Campbell's formula gives

    mean(n) = n_c integral d^3r n_1(r),
    Var(n) = n_c integral d^3r n_1(r)^2.

At large `T` with a fixed smooth velocity spectrum and neglected `r0/T`,
write `N(E)=integral du F(u,E)`. Then

    CV^2 ~ [integral dv N(E(v))^2 gamma^6]
            / [4 pi n_c T^3 (integral dE N(E))^2].

This gives the `T^(-3/2)` coefficient-of-variation scaling of the smooth
model. It is not uniform in the limit of spectral lines. For a fixed single
speed `v` and constant source of duration `D`, only the shell from
`a=r0+v(T-D)` to `b=r0+vT` contributes. Its exact Poisson variation is

    CV = [4 pi n_c (v D) a b]^(-1/2),

which scales as `T^-1`. In particular, counting the entire ball would use
the wrong support for this counterexample. The committed knot-sweep code
integrates the parent discrete-speed shells exactly for the reconstructed
piecewise-linear time functions and checks the original radial-grid result.

These formulae give pointwise ensemble density variation, not a scattering
rate. Independent-source collision terms and same-source correlations are
different objects; in particular angular correlations are not determined by
the scalar density variance.

## Spherical projections

For a complex charged field or real neutral field on a finite-radius sphere,
project `q_lm(t,r)=integral dOmega Y_lm^* q(t,r,Omega)`. Orthogonality
implies the angular-gradient identity

    integral dOmega |grad_Omega q|^2
       = sum_lm l(l+1) |q_lm|^2.

`multipoles.py` uses this to report a finite-time field-gradient ratio and
uses a finite difference in radius for the radial derivative. A ratio of
these two classical field norms is not, by itself, a probability
distribution for particle momenta. A hypothetical plane-wave energy `k`
can convert `l(l+1)/r^2` to a ratio of squared gradients, but it cannot
create a phase-space distribution or determine two-species correlations.

For the perturbation difference field `delta q`, `angular_flux.py` reports
the signed time-integrated `-Re[(partial_t delta q_lm)^* partial_r delta q_lm]`
by `l`. This is not the full difference of perturbed and unperturbed fluxes:
the latter also contains interference with the background field. The
angular multipoles, phases, time-frequency resolution and incoming/outgoing
character all matter for any later particle matching.
