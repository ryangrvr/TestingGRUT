# Claude Code review: energetic conversion products

This is the continuation of your provisional review of bb211d1. The
independent review has not been represented as a rerun or certification.

The new calculation keeps the imposed source and frozen Q=200000 objects.
It provides first-conversion histories, reverse reactions, the time/momentum
source at exit from the reactive band, analytic post-escape horizons and
rest-mass accounting. The key correction is omega/m_X=.861195 rather than
the average E/(Qm_X)=.923316. The kinematic concern survives.

Please independently check:

1. First-event survival along an unconverted trajectory, separate from net
   anti-X abundance after inverse conversion. At z_source=1e6, the initial
   rate is about 27 H and half first-convert by z approximately 974500.
2. The delta=K_anti/m_X variable and the absorbing positive-particle
   boundary at E_X/m_X=1+2Omega. Verify number/charge/energy and the exact
   comoving momentum assigned at each crossing.
3. The escaped source horizons, number weighting, causal prefix bound,
   and the distinction between traveled path and spatial transfer function.
   The active-grid path-moment estimator failed its no-reaction control
   by 8.1% and is explicitly excluded. Do not bank those diagnostic fields.
4. The mass fraction: initial anti-X are 3.175% of the old calibrated rest
   budget; complete conversion yields converted X at 3.359% of remaining
   rest matter, and all free X at 6.718%. The old abundance calibration
   loses about 5.469% in this fixed-profile completion. No yield is retuned.

The next decisive calculation is a linear mixed-matter perturbation
calculation using the supplied time/momentum source and appropriate
elastic/conversion angular kernels, followed by a named data likelihood.
It must retain the imposed-source/frozen-profile scope. The actual soft
formation source and general hot-core leakage remain open, as does vacuum
density. Earlier packets and canonical GRUT-RAI are unchanged.
