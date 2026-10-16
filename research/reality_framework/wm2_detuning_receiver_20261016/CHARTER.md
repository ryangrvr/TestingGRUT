# WM2-S3: finite receiver and a fixed detuning interval

Base: WM2-S2 9cbfef8b70f90c875df7b87b59b945c7bc99d02e.
Preserve S2, S1, WM1, CP-12, SCOUT-1 and canonical GRUT documents.

Freeze this charter before the S3 numerical controls. Keep the single-end
prefix grammar. Test mismatch accommodation, not an asserted symmetry
that silently equates the storage energies.

## Declared law and resources

Add a shared bosonic energy receiver with L modes, initial vacuum.
Its one-particle energies lie in [2,6], represented by a Jacobi mode
basis with diagonal Delta=4 and neighboring mode coefficient -1.
Mode index is internal receiver structure; it is not an event coordinate.
There is no map from event labels to receiver modes. The same contact
mode couples to attachment at the movable end for all ordered event pairs.

Storage energy is U-Delta+delta, U=16. Emission-assisted append is
-(g/N) sum X_uv q_v^+ a_1^dagger plus its adjoint. The normalization
changes from S2 and must be declared. Retain S2's probe-assisted bare
append with epsilon=1/4 and matter hopping J=1/4, and mismatch penalty
Lambda=64. All reverse terms and all reservoir/storage configurations
in the invariant component are included.

The new exact conserved charge is sum q - n_receiver = 1, so the
receiver has at most N-1 quanta. This is a derived excitation cap,
not a numerical cutoff. A finite implementation uses L_N=64 N^2 modes.
No bath reset, irreversible rate, drive, Markov limit or work schedule
is introduced.

## Perturbation class frozen in advance

Uniform on-event storage detuning delta in [-1/2,1/2] and common
emission coefficient g in [9/10,11/10]. The interval is independent of N.
Other coefficients remain as above. This is a nonshrinking class of
permutation-symmetric perturbations. It is NOT all independent site
disorder, arbitrary bath changes or arbitrary new Hamiltonian terms.

The required primary target is p(m>=max(2,ceil(N/2))) >= 1/2 throughout
[4 N^2,8 N^2] for this finite closed model, uniformly over that class.
The original CP-12 neutral ray is tensored with the declared uniform-root
bondless seed, root q=1, other q=0, and receiver/probe vacuum.

## Two separate protection gates

A. Operational gate: complete full-H source-response and propagation
bounds, with constants independent of N and the declared perturbations,
all receiver/storage states retained, and nonzero nearby response.
Real reservoir motion may be resummed inside configuration blocks.

B. Spectral gate: do not inherit S2's isolated low band or real-energy
Schur bound. Test for Q near-resonances. Any positive-imaginary-energy
resolvent bound must explicitly keep its energy-resolution parameter;
it cannot be passed off as the eta->0 result.

A positive A result does not establish B or a low-energy spatial phase.

## Budget and honesty rules

Derive the exact quadratic emission reduction and finite-receiver error
before increasing trajectories. Cross-check it against unreduced small
register matrices. Include exact charge, energy and recording accounts.
Test Q excursions and fixed-energy-resolution Schur blocks with every
allowed storage pattern. Prove or reject a uniform real-energy gap.

The local shell became unavailable during this campaign. Use executable
dependency-free JavaScript in the available V8 runtime; report precisely
which checks ran. Do not claim unrun Python, large CP-12 trajectories,
interval arithmetic or independent validation. Publish the code and
saved evidence on an isolated branch.

Dimension is still supplied by the single-end grammar. Original charged
CP-12 matter, generic perturbations, equilibrium/thermodynamic phases,
universal gravity and distinctive empirical predictions remain unearned.
