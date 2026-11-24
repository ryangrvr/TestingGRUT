# GRUT-D4 reproducible packet

D4 tests historical repair without an external erasure flag. Its fixed
Hamiltonian acquires two records, responds to local parity discrepancies,
and supports four separated single faults with the same finite reservoirs.
The [report](REPORT.md) includes separate verdicts and the simultaneous-fault
ambiguity.

| File | Content |
|---|---|
| [MODEL.md](MODEL.md) | Frozen Hamiltonian, symbols, prepared state and test definitions |
| [PROOF.md](PROOF.md) | Finite-time writing, repair, storage and recurrence theorems |
| [INFORMATION_BOUND.md](INFORMATION_BOUND.md) | Two-copy obstruction, multi-label scaling and environment capacity |
| [CONTROLS.md](CONTROLS.md) | Exact failures, energy, entropy, resources and timing |
| [CRITIC.md](CRITIC.md) | Separate proof and independent small-matrix review |
| [PROTOCOL.json](PROTOCOL.json) and [FREEZE.json](FREEZE.json) | Frozen design and hashes |
| [RESULTS.json](RESULTS.json) | Uniform certificates, 243 comparator histories, and controls |
| [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json) | Separate small physical matrix results |
| [VERIFICATION.json](VERIFICATION.json) | Frozen inputs and numerical replay |
| [EXECUTION_RECORD.md](EXECUTION_RECORD.md) | How the results were produced and checked |

With Python 3.12 and the pinned dependencies installed, run from this
directory:

```bash
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 independent_check.py
```

`verify.py` recomputes the archived results and checks the frozen input
hashes. The separate checker builds the physical Hamiltonian directly for
small reservoirs. `python3 run.py --output PATH` writes a new numerical
record. The main theorem certifies the large writer-on system. The large
numerical calculation is the exact ready-sector comparison for this proof.
