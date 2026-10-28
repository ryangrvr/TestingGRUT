# Isolate the object from periodic returns

Added while the hot spatial/temporal controls run, before any radial result. The
periodic three-dimensional object test is retained. Add an open spherical control
to distinguish emitted radiation from box recurrences, using the previously
reviewed finite-volume radial engine with the M1 cubic force/energy.

Use the exact M1 reference profile at Q200000. Radial disturbance is the l=0 part
of the three-dimensional five-percent disturbance: F multiplied by
1+.05 exp[-(r/12)^2], S increased by .05(1-S) times that Gaussian. Adjust the
rotating momentum to keep Q exactly200000. T=1000/m_X. Primary dr=.2, dt=.025,
outer radius200, absorber begins100, damping coefficient .2. Include unperturbed,
dr=.1/dt=.0125, dt=.0125 at dr=.2, and outer radius300/absorber begins200 controls.

No damping is applied to the central region. Every removed unit of energy and
charge goes to an explicit receiver ledger. Check Q+Q_abs and E+E_abs, report
the actual finite-volume convergence and compare core quantities at r=12,20,30.
Do not identify this radial pulse with the nonspherical three-dimensional pulse
or with actual cosmological formation. It tests nonlinear scalar survival and
outgoing radiation, not the final cosmological antiparticle population.

The radial engine uses q_scaled=(F_x,F_y,S), energy/charge multiplied by g_X²
relative to the canonical m_X-unit three-dimensional engine, and the opposite
phase-rotation convention. Its reported charge is y px-x py, so initialize py<0.
Convert Q_phys=Q_scaled/g_X² and E/(m_X)=E_scaled/g_X² once. These conventions
are disclosed rather than equating the two arrays or their raw amplitudes.
