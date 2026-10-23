# Verification performed and boundaries

Same author as construction. No subagents or independent replication.

Five complete 3D NR onset runs were executed: seed 17 at 64^3, 96^3,
128^3 with maximum dt=.02, seed 17 at 96^3 with maximum dt=.01,
and seed 29 at 96^3 with maximum dt=.02. All ran until the declared
spectral limit. Each archive contains diagnostic CSV, summary/component
JSON and compressed snapshots. NPZ snapshots retain a spectrally truncated
32^3 complex field (complex64) and full-resolution middle-plane density
(float32) at the common saved states and the first stopped state. These
are diagnostic reductions, not full 3D checkpoints. Dynamics use complex128.

197/197 checks in `VERIFICATION_RESULTS.json` passed:

* Charge/Hamiltonian drift and last-safe spectral fractions in all five runs.
* Common-state grid and time-step comparisons at tau=0 and 6. Equal phases
  are compared directly, without an aligning fit. The common fields are
  truncated diagnostics; this is not a dense-core convergence assertion.
* Shared-grid agreement of the identical continuous initial Fourier series,
  periodic-face connected-component handling, and analytic growing-mode
  controls in an actual 3D cube for three wavevectors. Halving dt reduces
  the linear-growth error by about four, as expected for Strang splitting.
* 48 debris kernels: 16,000-to-32,000-point heat refinement, integrated decay
  probability and comoving particle-energy accounting (decay + cosmological
  redshift + remaining energy). Maximum relative heat refinement error is
  1.36e-6. This is kernel convergence, not validation of deposition assumptions.
* Pure-radiation cold-decay closed form, threshold amplitude versus full
  angular integration, angular quadrature refinement, and a separate
  solve_ivp route for the stated irreversible-channel analytic result.

The initial pilot inspected high momentum only at .1 output intervals;
accelerating collapse could require too many adaptive steps before the next
inspection. Two such in-progress pilots were interrupted. The finalized
code checks every five steps once peak density exceeds 8, as well as at
outputs; all five archived runs used this stop implementation. A filename
suffix bug was corrected before the final campaign completed: decimal dt
strings now append extensions rather than replacing `.02`. Completed
nominal files were renamed accordingly; the half-step run was restarted.
Neither correction changes the PDE, initial state or physical parameters.

No full-field dense-core simulation, final bound-object count, nonlinear
emission spectrum, complete hidden reaction network, BBN nuclear rerun,
updated FIRAS likelihood or multi-target detector forecast was executed.
The earlier wall and detector calculations were not rerun. The existing
197 checks validate the specified reduced calculations only.

The standalone figure was generated with Matplotlib and visually inspected.
`MANIFEST.json` records bytes and SHA256 for every published file beside it,
except the manifest itself; generated Python bytecode is excluded.
