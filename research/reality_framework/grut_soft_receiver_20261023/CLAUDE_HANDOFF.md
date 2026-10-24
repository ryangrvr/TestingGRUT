# Claude review: the slow-field receiver has a measured boundary problem

Continuation of 59a1a71. This calculation tests the numerical boundary
independently of the collapse integrator, using existing author traces. It is not independent-author
replication or a claim that Claude has executed this packet.

Please check the exact damper and endpoint in the old radial_nodal.cpp against
the free continuum equation derived in MATHEMATICS.md. The original cubic
layer returns about 57.1%, 31.7%, 16.7%, 3.41% and 14.8% of test-wave energy
at p/m_X=.01,.02,.03,.05,1 respectively. These are unit-wave reflection
coefficients, not fractions of the collapse population. A positive absorber
energy/charge ledger does not certify a reflection-free escaping source.

Challenge the initial-value result using the independently implemented
boundary-value solver, the constant-slab closed form, zero damping and the
flux-loss integral. Numerical controls are separated from the failed source
gate. Wider ordinary cubic dampers improve part of the band but the tested
candidates all fail the one-percent target at k=.01.

The N0125_box and N0125_far physical realizations are bitwise identical.
Use that pair to audit radius/window dependence. At the original probe the
conditional outgoing linear-mode energy passes the existing local cut only
at about 25%; this is not a recovered 75% particle population. Conditional
all-band antiparticle estimates also change under absorber-strength controls.
Do not convert those conditional numbers into a physical pair abundance.

The concrete next build is a causal massive-wave exterior receiver with
explicit initial-memory conventions, discrete stability/reflection controls,
and charge/energy plus time/momentum output. The ideal complex-stretch
formula is only a frequency-domain design reference, not an implemented
time-domain fix. Nonlinear exterior forces and any remaining bound component
must be controlled separately. Existing finite-time core-charge convergence
and the separately conditional conversion-history packet are not replaced.
