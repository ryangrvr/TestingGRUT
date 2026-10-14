# Reproduce this CP-12 continuation

Use a checkout containing parent `ae91c9beb0d8dda0f7692ea8bab0c5f808355530`
and this additive directory. No external data downloads are required.

```sh
python -m pip install -r research/reality_framework/cp12_local_circulation_resource/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_local_circulation_resource/verify.py
```

Python 3.11+; the recorded run uses the exact versions in VERIFICATION.json.
NumPy and SciPy are both declared dependencies. Isolated Python mode is supported.
The full check typically takes well under a minute on the research host. Peak
memory is governed by a 1,875-dimensional dense matrix, not the 125-vertex
many-body Hilbert space. The latter check streams only one Hamiltonian action.

The verifier checks hashes, then executes all four checkers in fresh processes.
New numerical outputs go to a temporary directory and are compared against the
committed records (elapsed time ignored; floating values compared with explicit
tolerance). It leaves the checkout unchanged. `--record` is a maintainer option
for writing VERIFICATION.json before regenerating MANIFEST.json.

- REPORT.md: physical conclusions and limits.
- DERIVATION.md: complete counting, Gibbs, preparation, and Wilson-response proofs.
- check.py: exact graph packing counts, rational local-loop ranks, and tori.
- wilson_check.py: full unregulated Hamiltonian moments on Q=2N and Q=3N hosts.
- quantum_check.py: fresh full small matrices, Gibbs identities, unitary quenches,
  and a finite-impulse gauge-response control.
- graph_results.json, wilson_results.json, quantum_results.json: executed outputs.
- VERIFICATION.json, MANIFEST.json: execution record and provenance.

Dependencies from the parent tree are SHA-checked. `frozen_model.py` and
`exact_destination.py` are imported from the published ground-connectivity
packet; neither is edited. The parent's exact holonomy checker is run verbatim.

The theorem concerns the complete labelled physical Hilbert space, not only its
permutation-invariant subspace. The small numerical calculations are regulated;
the unbounded-rotor theorem is analytical. Prepared-host gauge capacity is not
geometry selection, and a moment is not a propagating photon or graviton pole.
