# Reproduce the restored packet

```bash
cd research/reality_framework/grut_action_freezing_20261021
python3 code/verify.py
python3 code/verify.py --rerun
```

Python 3.12.14,numpy 2.3.5,scipy 1.17.0; exact environment in ENVIRONMENT.json.
No network, sympy, GPU, private dataset or new dependencies required. Same
author, not an independent scientific implementation.

The manifest requires the exact payload leaf set and SHA256 identities.
RAW_DATA.tar.gz contains all **16 unabridged raw CSV trajectories** (four full,
twelve late including refinements); RAW_MANIFEST.json checks exact members,
sizes and hashes. The verifier reads them without unpacking and independently
reconstructs metric velocity, scalar potential, energy balances and endpoints.

`--rerun` copies source into a fresh temporary directory, reruns every integration,
compares the complete result leaf schema and 237 scientific leaves at 2e-6
relative tolerance with zero absolute allowance (tiny densities cannot pass by
being replaced by zero). Wall time, solver counters and near-zero residual
magnitudes need not match bitwise; all fresh original gates must pass. No
committed file is overwritten. The first pre-reset artifacts are not claimed
to survive; the published archive is the restored execution.

To inspect the raw CSV, safely extract the known archive into an empty directory:

```bash
mkdir raw_inspection
tar -xzf RAW_DATA.tar.gz -C raw_inspection
```

Use a **separate working copy** to regenerate, because this overwrites outputs:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 code/build.py
python3 code/make_manifest.py
python3 code/verify.py
```

The optional `python3 code/plot.py` uses matplotlib 3.10.8 and regenerates the
figure from the archived CSV. Run it before regenerating a working manifest
if you want a new figure. The restored headline-reference comparison uses
values recorded in the conversation; it is not another empirical test.

The generator archives and removes its loose CSV only after writing all raw
trajectories and their manifest. Compression has deterministic metadata and
does not alter the numeric content. Seven analytic controls plus 28 numerical
gates comprise eight full conservation, two full refinement, twelve late work,
six late refinement checks. Refined runs tighten tolerances by ten and halve
steps. Late Friedmann is algebraic and is not an extra verification count.
Integrity and rerunning do not establish quantum preparation or empirical truth.
