# Lower-density control, fixed before running it

The primary F=1 result at T=300 has large, merging connected charge regions rather
than clean isolated spheres. Keep that result. It is not legitimate to replace its
mask fraction with a formed-object abundance or to call an embedded low-E/Q patch
an isolated object.

Before computing this auxiliary control, choose F=.5, L=64, T=600, dt=.04, N=64;
N=48 spatial control and seed2. All potential parameters and seeds' construction
are unchanged. The longer time follows from the already-published linear growth
rate .02922 m_X (versus .11428 at F=1). The larger box samples the longer wavelength.
The lower density reduces the expected volume occupancy. These choices test the
mechanism under another supplied cold preparation; they do not constitute a new
cosmological onset, a successful replacement for a failed primary, or a retuning
of M1's couplings or mass criterion.

The elliptic N64 control ends with 20% of kinetic-plus-gradient energy in the
axial high-momentum shell. This is an unresolved numerical concern. Preserve its
outcome but do not claim a continuum hot-population verdict. Add N96 at the same
L48 and T300, dt=.04; also N64 with dt=.02 to separate time integration from the
spatial cascade. These controls respond to the measured numerical concern, not
to an attempt to recover a preferred physical abundance.

Also compute the zero-pressure infinite-volume core by minimizing
Omega_bulk²=S²+2 sqrt[a U_phi(S)], with a=lambda_X/gX²,
U_phi=lambda_phi/(4gX²)(S²-1)²-epsilon/gX²(S-1)³. Its density and energy per charge
are thermodynamic diagnostics, not the finite-object size distribution.
