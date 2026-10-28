# A single-vacuum binding sector: additional construction test

Written after the first scale-matching runs, before solving the changed profiles.
The owner authorizes developing a version that works. Address the exact-phi-parity
wall problem with one allowed scalar-potential deformation, not a new particle.

Use

V_phi = lambda_phi/4 (phi^2-v^2)^2 - epsilon v (phi-v)^3.

At phi=v, V=V'=0 and V''=2lambda_phi v^2: the existing v and vacuum masses remain
the renormalized targets. In u=phi/v, the derivative factors into (u-1) times
lambda_phi*u*(u+1)-3epsilon*(u-1). The latter quadratic has no real zero when
epsilon/lambda_phi is between1-2sqrt(2)/3 and1+2sqrt(2)/3. The positive potential
factor also has no second zero for0<epsilon<2lambda_phi. Thus epsilon=.002 at
lambda_phi=.03 gives one stationary vacuum, rather than two degenerate phases.

Freeze the test grid epsilon=0,.001,.0018,.002,.0025,.003 at fixed Q200000, gX=.1,
lambdaX=.001, lambda_phi=.03, v=.924089637375362GeV. The primary changed point is
epsilon=.002, chosen by the analytic single-stationary-point criterion before its
profile is computed. Use Q1000000 only as a labelled existence/charge sensitivity.
Do not select a final Q distribution. Continue from the archived refined profiles.

Required construction checks: bounded positive scalar potential; unique stationary
vacuum at the primary point; a converged nodeless stationary object with E/Q<mX;
neighboring-charge chemical slope; virial and first-law checks; spatial/domain
refinement. A full nonradial fluctuation spectrum and cosmological formation remain
separate tests. The absence of two vacua removes this specific domain-wall mechanism,
not every possible defect or nonequilibrium interface.

If the primary point passes, recompute its energy and form factor, then repeat the
physical-reheating abundance calibration using that energy. Do not keep the old
26.7keV profile node or old17.064TeV mass. Record epsilon as a new supplied EFT input,
not a consequence of closure. The very early force is suppressed by(v/H)^3; verify
the size before transporting the early source approximation. Vacuum density still
needs its explicit renormalized matching input.
