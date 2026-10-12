# CP12 reproducibility and execution accounting

Python 3 with NumPy and SciPy. Use `OPENBLAS_NUM_THREADS=1` to keep small sparse calculations reproducible and avoid BLAS oversubscription. Run from this directory, or give the full script path. Scripts write JSON next to themselves. Times and final decimal digits are machine-dependent.

## Quick physical and independent controls

```bash
OPENBLAS_NUM_THREADS=1 python independent_controls.py
OPENBLAS_NUM_THREADS=1 python validate_reduction.py
OPENBLAS_NUM_THREADS=1 python record_dynamics.py
OPENBLAS_NUM_THREADS=1 python distance_control.py
OPENBLAS_NUM_THREADS=1 python gauge_shortcut.py
OPENBLAS_NUM_THREADS=1 python local_reconnection_control.py
OPENBLAS_NUM_THREADS=1 python verify_packet.py
```

The independent electric enumeration and explicit orbit-isometry checks use constructions different from the optimized eigensolver path. `verify_packet.py` audits saved outputs and compares separately executed cutoffs; it does not rerun expensive diagonalizations or serve as independent empirical confirmation.

## Full joint ground states and hostile conditional controls

```bash
OPENBLAS_NUM_THREADS=1 python check.py
OPENBLAS_NUM_THREADS=1 python check.py --n 5 --q 5 --L 2
OPENBLAS_NUM_THREADS=1 python check.py --n 5 --q 5 --L 3
OPENBLAS_NUM_THREADS=1 python check.py --n 5 --q 6 --L 2
OPENBLAS_NUM_THREADS=1 python check.py --n 5 --q 6 --L 3
OPENBLAS_NUM_THREADS=1 python check.py --n 4 --q 4 --L 3 --mask 45
OPENBLAS_NUM_THREADS=1 python check.py --n 4 --q 4 --L 3 --mask 15
OPENBLAS_NUM_THREADS=1 python cycle_controls.py
python geometry_summary.py
```

N=5,Q=6,L=3 has 2,866,255 states and nearly 30 million nonzeros; it requires substantial memory. `--mask` disables reconnection and is explicitly a conditional fixed-graph comparator, not the full model. `cycle_controls.py` compares conditional rings through N=8, all at cutoff 2. Their tiny finite-size energy differences should not be used to claim a full fluctuating-geometry phase.

## Symmetry-sector extension

```bash
OPENBLAS_NUM_THREADS=1 python symmetric_sector.py --n 4 --q 4 --L 2
OPENBLAS_NUM_THREADS=1 python symmetric_sector.py --n 5 --q 5 --L 2
OPENBLAS_NUM_THREADS=1 python symmetric_sector.py --n 6 --q 6 --L 2
OPENBLAS_NUM_THREADS=1 python symmetric_sector.py --n 6 --q 6 --L 3
```

The reduced matrix equals the exact vacuum-relative invariant-sector restriction. Its N=6 ground need not be the ground of the entire system. No untested extrapolation to N=7 or a continuum is included.

## Thermal and actual preparation controls

```bash
OPENBLAS_NUM_THREADS=1 python thermal_control.py --n 4 --q 4 --L 2
OPENBLAS_NUM_THREADS=1 python thermal_control.py --n 4 --q 4 --L 3
OPENBLAS_NUM_THREADS=1 python thermal_control.py --n 4 --q 4 --L 4
OPENBLAS_NUM_THREADS=1 python thermal_control.py --n 4 --q 4 --L 5
OPENBLAS_NUM_THREADS=1 python quench.py --L 3
OPENBLAS_NUM_THREADS=1 python quench.py --L 4
OPENBLAS_NUM_THREADS=1 python quench.py --L 5
```

Thermal diagonalization is dense; cutoff 5 has dimension 10,275. The beta=0.5 result retains a visible cutoff effect and is marked unresolved at high precision. Quench cutoffs 4 and 5 agree on the T=20 cycle probability at 1.6e-8. Neither a ground-state solve nor closed unitary evolution supplies thermalization or cosmological preparation.

## Partner verification

The CP11 magnetic script was fetched from commit 2bca657df1a3539449f0468d06280a9b56e079a8 and executed before CP12. The owner later supplied SCOUT-RG01 and RG02. The exact RG02 `full_joint.py` and `shortcut_control.py` were downloaded from bbac3f19ddf7cca12393cc158aad40380b4e969d and run verbatim. Their hashes and results are in PARTNER_EXECUTION.json. The code remains on its original branch; no scout file was edited. Dense N=3 and random-start N=3,4 ground-energy checks were separately executed. The independent random-start command used the original `assemble` function, so it independently checks eigensolving, not construction of that Hamiltonian.

## Accounting

One CP12 law, one frozen parameter tuple. Multiple Q sectors and cutoffs were declared controls. No failed candidate was hidden and replaced. No coupling was tuned after the outputs. New partner compatibility tests are analytical/operator controls, not a CP13 law. The remote law commit records the earlier local freeze honestly; computations preceded remote upload. Earlier CP packets, Stage-3 Card 3, canonical GRUT-RAI and draft PR12 were not changed.
