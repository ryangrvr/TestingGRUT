# M1 expansion-clock diagnostic

28 October 2026. Parent: c4815306753ece00dad2b03412f0ed9fe89b8f92.

This independent branch investigates the prompt's evolving cosmology while a
separate worktree studies late radial emission and the stationary family. No M1
parameter or yield is retuned. The initial source is the original supplied 10%
pulse in the parent packet. The epoch T=20 MeV, Y_Q and one Q=200000 core per
reference charge are hypotheses, not a formation calculation.

First derive the relation between dimensionless collision time and the FRW scale
factor, and implement the parent's frozen-profile, prompt-h reduced collision
operator with radiation expansion. Use the distribution per reference core,
natural units and a full energy/work ledger. Work in x=ln(a/a0), approximate
constant g_star=10.75 and radiation domination, and evolve the redshift term
with a conservative first-order upwind energy-space generator. Its diffusion
is a known physical limitation; compare against exact free-streaming
characteristics and vary the energy grid. Include independent early-time
comparisons to the parent static solution, and check positivity, charge,
energy/work balance, support losses, ODE tolerance and angular controls as
needed.

The main run holds the Q=200000 scattering strengths and frequency fixed as in
the parent; Q evolves only in the charge ledger. A cold-mass estimate may use
the stationary BVP at the computed Q, but this does not turn the history into
a coupled moving-core calculation. The residual core excitation is retained
as unprocessed energy, with no invented emission rate. Late radial radiation is
initially excluded, then included only as a clearly labelled sensitivity if a
validated source can be obtained. Do not silently clip support or erase energy.

No claimed actual formation, final abundance, full freeze-out, independent
author review, or successful 1% gate may follow from this restricted model.
If redshift or collision-grid convergence fails, preserve and report the
failure. Numerical control is not a physical error bound.
