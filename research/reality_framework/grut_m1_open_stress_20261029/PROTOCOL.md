# M1 pressure, open multipoles, and coherent return controls

Parent GitHub commit `6dc95562635287f3aca47de078395d3ef0fbec97`, tree
`4e9df14fee6fe0912a19eac71832ecd793998dbc`. The supplied review did not run
repository code. Test its proposals with explicit numerical controls; do not
accept any conclusion about future pressure or the 1% gate by assumption.

1. Instrument the inherited spherical finite-volume evolution to record the
   integrated spatial stress trace and the virial and boundary-flux terms at
   radii 20 and 60. Recreate the same 10% pulse and compare against the
   published energy/charge histories. Measure time-resolved and averaged
   stress, include grid/time controls, and compare a same-lattice stationary
   solution. The 0–20000 data must be regenerated because the earlier raw
   fields were temporary; record any interval or resolution not obtained.
2. Derive and test the exact finite-region conservation identity
   `d/dt integral x_i T^0i = integral T^i_i - surface x_i T^ji n_j` in the
   no-absorber region. A finite-radius pressure trace, an isolated
   Minkowski evolution, or a short unweighted mean does not by itself give a
   coupled FRW stress or a bound on all later expansion work.
3. Construct an even-Legendre, axisymmetric radial-angular solver using the
   same local nonlinear action and an open outer absorber. Verify its pure
   monopole limit against the original spherical solver; test radial and
   multipole/angle resolution. The old 3D disturbance is proportional to
   `x^2-y^2` and is not a rotated axisymmetric `Y20`; any new axisymmetric
   pulse is an explicitly different control, not a reproduction of it.
4. Measure monopole/quadrupole charge and energy flux at `r=60`, including
   negative-signed components and outgoing/incoming ambiguity; don't infer
   particle pair angular correlations from field multipole intensity.
5. Run a small reflecting spherical-cell control to test coherent returning
   radiation and compare against an open run. Its density, coherent phase,
   boundary and arrival times are not a calibrated cosmological gas.
6. Audit the proposed 33.7% energy-conversion extrapolation. It is a
   two-trajectory finite difference, not a demonstrated rule for remaining
   energy emitted at a different time or with a different spectrum.
