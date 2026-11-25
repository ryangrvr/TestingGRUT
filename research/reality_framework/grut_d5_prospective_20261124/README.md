# D5 reproducible packet

D5 uses a fixed local Hamiltonian to acquire a record, generate a
word from blank sites, recover partial development, and regenerate
a damaged completed word. [Report and verdicts](REPORT.md)

| File | Content |
|---|---|
| [MODEL.md](MODEL.md), [PROTOCOL.json](PROTOCOL.json), [FREEZE.json](FREEZE.json) | Microscopic design, frozen checks and hashes |
| [PROOF.md](PROOF.md) | Exact path and endpoint identities, time windows, recurrence |
| [CONTROLS.md](CONTROLS.md) | Matched record and rule interventions, failure cases |
| [LIMITS.md](LIMITS.md) | Indistinguishability theorem and pilot-null |
| [RESOURCES.md](RESOURCES.md) | Energy, entropy, hardware and archive scaling |
| [CRITIC.md](CRITIC.md) | Separate critique and direct physical validation |
| [RESULTS.json](RESULTS.json) | All frozen construction, control and lab-pulse rows |
| [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json), [VERIFICATION.json](VERIFICATION.json) | Direct physical reference, replay and cross-check |
| [EXECUTION_RECORD.md](EXECUTION_RECORD.md) | How the packet was produced |

With Python 3.12 and the pinned dependencies installed, run from
this directory:

```bash
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 oracle.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify.py
```

`python3 run.py --output PATH` writes a fresh production ledger.
The exact history reduction avoids enumerating the exponential
tape Hilbert space. The independent physical graph enumerates
states reachable by local forward and reverse hops and the
prescribed tape pulse, then uses direct sparse propagation.
The main lab-pulse scores include the clock endpoint and every
failure to arrive.
