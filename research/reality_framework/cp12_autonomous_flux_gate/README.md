# CP-12 autonomous flux and internal graph-record gate

Extends `4e2b7a876a3a45ed554408ed803b572dd9cab580` without changing the
frozen Hamiltonian. Read REPORT.md for results and DERIVATION.md for proofs,
assumptions and numerical limits. Main claim: persistent partial protection
and internal graph records occur in the tested unforced states. Spatial
selection, extended causal propagation and spacetime emergence remain open.

## Reproduction

Python 3.12, NumPy and SciPy suffice for the computations. From the repository
root, set a PRIVATE trusted scratch cache (pickle files must not be taken
from untrusted sources). The cache is optional and is not a physical input.

```bash
export CP12_FORMATION_CACHE=/tmp/cp12_private_flux_cache
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python research/reality_framework/cp12_flux_protection_tradeoff/check.py
python research/reality_framework/cp12_autonomous_flux_gate/dynamics.py --n 4 --q 4 --budget 256 --end 20
python research/reality_framework/cp12_autonomous_flux_gate/dynamics.py --n 5 --q 5 --budget 128 --end 10
python research/reality_framework/cp12_autonomous_flux_gate/dynamics.py --n 6 --q 6 --budget 128 --end 10
python research/reality_framework/cp12_autonomous_flux_gate/dynamics.py --n 5 --q 6 --budget 128 --end 10
python research/reality_framework/cp12_autonomous_flux_gate/dynamics.py --n 6 --q 12 --budget 6 --end .4
python research/reality_framework/cp12_autonomous_flux_gate/dynamics.py --n 6 --q 12 --budget 10 --end 1
python research/reality_framework/cp12_autonomous_flux_gate/refine_counts.py
python research/reality_framework/cp12_autonomous_flux_gate/records.py --n 4 --q 4 --budget 256
python research/reality_framework/cp12_autonomous_flux_gate/records.py --n 5 --q 5 --budget 128
python research/reality_framework/cp12_autonomous_flux_gate/mobility.py
python research/reality_framework/cp12_autonomous_flux_gate/energy_audit.py
python research/reality_framework/cp12_autonomous_flux_gate/deletion_audit.py
python research/reality_framework/cp12_autonomous_flux_gate/summarize.py
python research/reality_framework/cp12_autonomous_flux_gate/verify.py
```

Default main evolution settings are time spacing 0.1 and eight leakage
subintervals per spacing. The mobility probe uses sixteen. Matrix exponentials
supply actual finite-time states; only leakage quadrature uses a Taylor
remainder bound. Electric budgets are numerical projections, not a change
of the physical law. All outgoing transitions are included in the residual.
The fast rewiring decomposition uses total vertex fermion count: every
hopping matrix element changes it by one, whereas rewiring and diagonal
terms preserve it. Independent raw-transition and full-space assemblies
check this optimization; graph-mask equality is not used to identify hopping.
The budget-ten dense count refinement uses two leakage subintervals per
spacing, explicitly recorded with its looser forward bound. This changes
the error quadrature, not the Hamiltonian or the finite-time observable.

The Q=N six-event calculation has 102,992 retained orbit states and 1,306,862
nonzero Hamiltonian entries. Its initial build and run took several minutes
on the recorded environment; caching saves the basis build. The five-event
record audit unfolds 748,410 amplitudes. The verifier freshly reruns the
four-event main trajectory, small independent matrices, source probes and
the four-event record audit. It validates hashes and invariants of the
larger saved runs; it does NOT silently rerun them all.

`--record-manifest` is for intentional source/data updates after reviewing
all fresh results. Running it is not a formal certificate of floating-point
roundoff. `requirements.txt` records the numeric environment; Matplotlib is
only needed to regenerate the optional plot.

## Files and limits

- `N*_Q*_B*_T*.json`: unforced trajectories. High-budget selected runs carry
  geometry TV and continuous-time protection bounds. Earlier regulator
  records have fewer metadata fields; the physical evolution and error
  method are the same. Late low-budget errors can be vacuous.
- `RECORDS_N*.json`: graph/environment partial traces, classical conditional
  support bounds and r=0 diagnostic. There is no record insertion or
  postselection in the main trajectory.
- `MOBILITY.json`: finite +/-pi/2 source probes and matched no-hopping
  diagnostic comparators. These are external diagnostic operations with
  explicit work, not an autonomous state-selection mechanism.
- `ENERGY_AUDIT.json`: full-law electric-energy intervals inferred from
  exact total energy minus a bounded complementary operator.
- `DELETION_AUDIT.json`: residual coherent-removal capacity in the actual
  state, from the exact incidence/rewiring identity. This is not a causal
  response or a physically applied dephasing operation.
- `REFINED_COUNTS.json`: source-free forward/backward detector error bounds
  resolving flux formation in the denser Q=2N sector at an early time.
- `EXACT_CONTROLS.json`, `LEAKAGE_CONTROLS.json`: combinatorial and independent
  physical-matrix checks. The algebraic proofs concern unbounded fields;
  finite matrices cross-check their implementation.
- `INDEPENDENT_EXTENSION_CHECK.json`: newly executed upstream checker and
  its source hash.
- `MANIFEST.json`, `VERIFICATION.json`: source/input hashes and fresh-check
  results. Fresh timing fields and ARPACK roundoff can differ on rerun.

Bounds are analytically justified but evaluated in ordinary floating point,
without interval arithmetic. Reduced-state diagnostics other than the saved
graph entropy are projected numerical values unless a separate bound is
given. A finite gate does not establish the long-time or thermodynamic phase.

The Q=N size comparison is a mechanism control: connected graphs have only
one independent cycle and do not satisfy the 3D gauge-capacity target. The
initial pure uniform ray, its coherent preparation and the conserved resource
sector remain supplied. Internal records of an arbitrary graph do not
constitute selection of space.
