# Core matching, late emission and mixing audit

29 October 2026. Parent published commit `8f67a83be4177f6baaa475c94565aab3b7c0cd94`,
tree `851ea4c4de2bbfce4cdae4ba88c3190f6ae814cd`. The external analytical
review is a prompt for tests, not independent numerical reproduction. Keep the
M1 action, initial perturbations and density fixed; preserve discrepancies and
negative results.

1. Extend the existing `h=.025,dt=.003125` and `h=.05,dt=.00625` radial
   field checkpoints from `m_X t=10000` to `20000` using the same absorber
   and equations. Check restart continuity and conserved charge and energy
   including the new absorber ledger. Measure late signed/positive-projected
   flux and compare the core at equal charge with same-grid stationary energy.
   No extrapolation beyond the observation interval is a proof of no future
   emission. The parent held remainder is a ledger residual, not a uniquely
   identified core oscillator.
2. Compare identical-run `h=.05` spectra and low-energy bins at probes 40,
   60 and 80; examine time-window and matching sensitivities. A boundary flux
   difference is not, by itself, a bound on interactions inside 40 or on the
   final energy or particle budget.
3. Recreate the supplied 3D 5% quadrupole-disturbed object and undisturbed
   controls at 64 and 96 spatial nodes. The periodic box has side 64, so
   fully contained spherical probes require radius below 32. Measure
   near-core field multipoles and outgoing/incoming diagnostics with controls.
   The 5% 3D excitation is not the 10% radial source or a derived formation
   history. Do not equate single-field multipole weights with the two-particle
   collision-angle distribution.
4. Re-evaluate the proposed Poisson mixing extrapolation with an explicit
   continuous-speed source as well as the previous discrete-energy model.
   Compare pointwise covariance and model clocks; distinguish ensemble
   smoothing from a guarantee for every point, isotropy, or actual mixing.
5. Separate the energy-ledger theorem (conditional on complete sector energy
   and nonnegative integrated pressure work) from source/formation and core
   dynamics. Catalytic conversion cannot create energy in a closed system,
   but the completeness and pressure assumptions remain open.
