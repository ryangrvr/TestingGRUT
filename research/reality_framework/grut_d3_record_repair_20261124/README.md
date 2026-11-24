# GRUT-D3 reproducible packet

GRUT-D3 gives an exact finite model in which a locally written record
selects which of two configurations is restored after one flagged erasure.
The proof covers repair and record retention together. The memory controls
match the full conditional energy distribution. The model, results,
assumptions and limits are summarized in [REPORT.md](REPORT.md).

| File | Content |
|---|---|
| [MODEL.md](MODEL.md) | Physical Hamiltonian, initial states and interventions |
| [PROOF.md](PROOF.md) | Invariant sectors, uniform time bounds and exact counterexamples |
| [CONTROLS.md](CONTROLS.md) | State and energy matching, resources and limitations |
| [CRITIC.md](CRITIC.md) | Separate internal critical review and full-matrix checks |
| [PROTOCOL.json](PROTOCOL.json) and [FREEZE.json](FREEZE.json) | Registered choices and frozen hashes |
| [RESULTS.json](RESULTS.json) | Finite dynamics and all scored probabilities |
| [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json) | Separate physical-matrix run |
| [VERIFICATION.json](VERIFICATION.json) | Frozen-input and numerical replay record |
| [EXECUTION_RECORD.md](EXECUTION_RECORD.md) | Provenance and verification scope |

From this folder with Python 3.12 and the pinned dependencies available:

```bash
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 independent_check.py
```

`verify.py` runs the full finite calculation, checks the uniform rational
bounds and compares every archived field within the declared tolerance.
The critic's script independently constructs small physical Hilbert-space
matrices. It compares the model against those matrices without needing to
build the exponential physical Hilbert space of the large reservoirs.
To generate a new production record, use `python3 run.py --output PATH`.
