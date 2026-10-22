# Verification actually performed

22 October 2026. One author; no subagents, new installations or outside
replication. This packet does not claim that earlier checks were all rerun
or that the proposed new fields have been observed.

## Scientific construction and numerical controls

The initial exploratory specification preceded execution. Its first run
passed 45 conditions and exposed an unexpected constant first-Born/mass
ratio. The stronger finite-charge identity was then derived and checked;
this after-result strengthening and subsequent numerical/scope additions
are explicit in AMENDMENT_01.md. Physical couplings, charges and matching
were unchanged. No empirical held-out test is claimed.

The final `calculate.py` run passes **51 internal numerical conditions**:
five fixed-integer-charge stationary profiles; larger-domain and tighter-
tolerance profiles at the three main charges; three scalar scattering
spectra at two resolutions; six finite-grid radial amplitude Hessians;
neighboring-charge energy differences; integrated field/virial and
mass/readout identities; zero-energy scattering and no-contrast controls;
finite-wavepacket discrimination and its momentum refinement.

| Quantity | Actual comparison |
|---|---:|
| Main profile mass refinement | largest relative difference 1.41e-11 |
| Main phase refinement | largest absolute difference 1.79e-11 rad |
| Wave equation versus variable-phase variables | largest difference below 5e-15 rad; common adaptive engine/steps |
| Refined mass–depletion identity, three main charges | relative residual below 2.2e-12 |
| Central d[E/(gv)]/dQ versus Omega | absolute difference 1.12e-10 |
| Fixed-Q radial amplitude Hessians | all six lowest sampled eigenvalues positive; smallest 0.217294 |
| Finite scalar scattering intervals | all three exact main lengths inside their analytic I<1 intervals |
| Packet momentum quadrature, 61 versus 121 points | squared-overlap difference 1.56e-15; both integer .95 counts unchanged |

These finite-grid and adaptive controls do not constitute interval-rigorous
continuum proofs. Radial positivity is not an all-channel lifetime. Probe
gravity is estimated by its small finite-momentum parameter, not evolved.
The full gravitational long-range scattering problem is excluded from
the short-range scattering-length claim.

## Read-only verifier

`code/verify.py` checks the full packet file set and hashes, the unchanged
source seed hash and the fixed physical input tuple. It reconstructs charge,
energy and the mass sum rule directly from all three archived profiles.

It then uses an independently organized **Numerov finite-difference**
scattering calculation from those raw profiles at nine charge/momentum
combinations. With a step refinement, the largest s-wave phase difference
from the adaptive results is 3.05e-9 rad, below the declared 2e-7 tolerance.
This is an additional algorithm by the same author.

An explicit 122-mode outgoing-state assembly and density-matrix eigenspectrum
give single-probe optimum .507233160200914, agreeing with the overlap
formula. A monochromatic one-mode density-matrix control has no distinguishable
global phase. Explicit summation of the binomial outcome probabilities
checks the fixed-measurement .95 count and the preceding unsuccessful count;
this does not rely on the calculator's likelihood-threshold implementation.

The final manifest covers **25 packet files** and excludes itself. Python
syntax, JSON, relative document links and git whitespace were checked.
The scientific figure was rendered and visually inspected. Earlier
scientific packets and the canonical GRUT-RAI implementation are unchanged.

## Fresh replay

`OPENBLAS_NUM_THREADS=1 python code/verify.py --replay` makes an isolated
temporary tree containing the new source, inputs and exactly hashed parent
seed. It recomputes the science and compares the complete schema:
**464 scientific JSON leaves and eight raw CSV arrays**. Maximum scientific
relative deviation is zero in the recorded environment; every raw array
comparison passes. Evaluation counts are computational diagnostics excluded
from physical comparisons. Repeated leaves are not independent discoveries.

This is same-code reproducibility plus the separate numerical identities
above. No blind independent author, external laboratory, empirical data
fit or outside peer review validates the result.

## Not performed

No new scalar-loop determinant, full quantum ground-state construction,
finite-Q quantum lifetime, all fragmentation partitions, dynamical writer,
physical D source/detector, cumulative repeated-read disturbance, coupled
metric scattering, realistic cosmological population, new vacuum selector,
unique outcome mechanism or explanation of consciousness. The prior
cross-repository synthesis is an input, not a fresh audit of every source
file in all repositories.
