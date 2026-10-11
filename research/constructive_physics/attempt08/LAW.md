# CP-08: reversible causal-order dynamics with matter transport

Frozen before executable tests, 2026-10-11. One new candidate, separately accounted. Authorization: owner “keep going work it out,” after CP-07's reported failures. Stage 3 inactive; Card 3 unopened. CP-07 is not modified.

## Supplied finite physical hypothesis

For a fixed N, let P_N contain ALL labelled N-element strict partial orders, not just orders compatible with a preferred numerical birth order. Geometry Hilbert space has orthonormal basis |C>, C in P_N. Matter is N qubits with ordinary complex quantum mechanics and Born probabilities. External evolution parameter t, Hamiltonian dynamics, fixed event inventory, tensor products, qubit axes and scales J, Omega, gamma are supplied; no spacetime interpretation of t is earned.

Join C to C' when adding one oriented relation i≺j changes exactly that relation and C' is still a strict partial order. Reversal removes it. Distinct orders are physical basis states, not irreversible records of every route taken to reach them. For addition define U_ij=exp(-i gamma Z_i X_j). Select the single Hamiltonian

    H = -J sum_additions [ |C'><C| tensor U_ij
                         + |C><C'| tensor U_ij^dagger ]
        + I_geometry tensor Omega sum_i X_i.

Parameters for the sole executed control: hbar=1, J=1, Omega=1/4, gamma=pi/4. Evolve exp(-i H t). Initial state is the empty order tensor |0...0>, unless a predeclared interference-path diagnostic starts at C={0≺2}. No Gaussian history weight, uniform physical probability measure or manifold preference is appended.

Time-independent Hermitian H conserves its total supplied energy and has a finite-dimensional positive-norm unitary interpretation. Addition/removal includes a global transitivity admissibility projector: microscopic locality is NOT presumed. Spectral lower bounds at each finite N do not establish a stable thermodynamic limit.

## Predeclared tests

- Enumerate all posets for N=3 and 4; verify exact counts, acyclicity/transitivity and that every Hamiltonian transition stays inside this space.
- Hermiticity, norm preservation, inverse evolution, energy conservation and affine mixture consistency.
- Full event-relabelling covariance, transforming matter factors as well as posets. This is not spacetime diffeomorphism symmetry.
- Derive and test interference of two additions with noncommuting matter transports. Compare gamma=0 at the same J, Omega, t without fitting. Scalar phases equal to differences of configuration potentials must be identified as a removable basis choice.
- Derive retarded response and symmetric fluctuations of a specified causal-relation projector from this same H and initial state. Test reversibility of records, and their possible effect on later probabilities when interfering routes are distinguished.
- Compute finite-size distributions of chain height. No inference of 4D geometry, Lorentz invariance, spin-2 gravity or universal coupling from two small N values.

A healthy finite joint dynamics is a legitimate conditional success; it is not a working theory of gravity. Do not add degree penalties, target-dimension terms or a tuned graph measure after seeing the results. No multiple candidates or parameter scan is hidden inside this packet.
