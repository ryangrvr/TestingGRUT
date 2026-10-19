# Candidate 1 two-object evolution

Parent: `96a7e7e6a49f49dc30e1c92f0855522af41d9641`. Isolated branch: `codex/grut-candidate1-two-object-20261019`.

The owner authorized trying the two-object calculation on a new branch. The ordinary shell workspace was unavailable at execution; the dependency-free JavaScript kernel can run both in the available V8 runtime and with Node.js. This package is a bounded classical axisymmetric (three spatial dimensional, head-on) calculation at Candidate 1's unchanged couplings. It does not derive spacetime, quantum theory, a formation ensemble, gravity or consciousness.

## Reproduction

From this directory, run `node run.cjs`, or `node run.cjs medium_pair_phase_0` for one frozen case. No package installation is required. The source profile is a documented decimation of the parent's archived CSV, not a new boundary-value solve. `protocol.json` defines the full run list and measurement/consistency criteria before those evolutions. A short single-object time-1 implementation pilot and a time-10 forward/reverse preflight preceded this freeze; it is retained in `pilot.json`.

Evolution uses the same real `phi` and complex `X` equations, with `D=0` and the flat-space benchmark approximation. Radial finite volumes preserve regularity at the symmetry axis; a weighted discrete Hamiltonian matches the spatial force. Kick-drift-kick Verlet preserves the global U(1) charge to roundoff and has a finite time-discretization error in energy. Boundaries have fixed vacuum ghost values. Domain and timestep controls are separately prescribed.

A superposition of Lorentz-boosted one-object profiles supplies the initial two-object state. It is not an exact two-object stationary solution; its interaction energy is included in diagnostics. All preparations, dimension, axisymmetry, couplings and initial resources are supplied. Opposite charge means opposite frequency rotation, not a change of coupling. Relative phase pi has an exact antisymmetry node; phases pi +/- 0.2 are included before evaluation to avoid treating a symmetry-enforced separation as a generic result.

Classification must use the full convergence evidence. Multiple density peaks alone do not prove persistent particle identity; sparse phase samples do not establish a populated phase or an open interval of successful initial data.

The full-run outcomes are initially unexecuted. A later report will state exactly what completed and what the numerical controls support.
