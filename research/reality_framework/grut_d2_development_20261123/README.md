# GRUT-D2 reproducible packet

The [REPORT](REPORT.md) gives the result, evidence and limits. A fixed local
Hamiltonian yields finite-window ordering, an internal reservoir count record
and specified disturbance recovery. The preferred relation and reservoir
preparation remain supplied; active developmental feedback and new GRUT
physics are unestablished.

| File | Purpose |
|---|---|
| [MODEL.md](MODEL.md) | Complete physical model, preparation, time and resource inputs |
| [PROOF.md](PROOF.md) | Exact reduction and analytic finite-window bounds |
| [CONTROLS.md](CONTROLS.md) | Counterfactuals, disturbance work, entropy and D1 corrections |
| [CRITIC.md](CRITIC.md) | Separate internal AI critique and independent small-matrix results |
| [CONTROL_CORRECTION_01.md](CONTROL_CORRECTION_01.md) | Preserved failed hot-bath record forecast |
| [PROTOCOL.json](PROTOCOL.json) and [FREEZE.json](FREEZE.json) | Published design freeze and exact hashes |
| [RESULTS.json](RESULTS.json) | Full registered distributions and controls |
| [VERIFICATION.json](VERIFICATION.json) | Same-code replay and hostile validation controls |
| [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json) | Full physical-Pauli oracle record |

```bash
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python3 verify.py --output /tmp/grut_d2_verification.json
OPENBLAS_NUM_THREADS=1 python3 independent_check.py > /tmp/grut_d2_independent.json
```

Run from this folder or use absolute script paths. The scripts resolve inputs
relative to their own location, so neither a Git checkout nor the original
working directory is required. The optional [figure](ordering.png) is rebuilt
by `plot.py` using `requirements-plot.txt`.
