# Reproduce the cubic topology discrimination test

Use a full checkout of the published branch, so the byte-pinned sibling modules are present.

```bash
python -m pip install -r research/reality_framework/cp12_topology_sensitive_preparation/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_topology_sensitive_preparation/verify.py
```

The verifier executes the exact parent checker and all four new scripts afresh. It verifies rational certificates, the non-invariant full-law correction, independent factorial/overlap formulas, resource matching and the analytic derivative asymptotics. Typical runtime on the recorded host is several seconds; hardware varies.

- `full_moments.py`: auxiliary N=6,Q=12, independent full unregulated Hamiltonian moments through T³. Projected contractions are exact; no time evolution or full diagonalization is asserted.
- `topology.py`: Q=3N throughout, degree-matched mask pairs at N=24,48,96,192, exact motif-response derivatives through T³.
- `symbolic_check.py`: formal polynomial identity for the cubic probability and leading coefficient, standard library only.
- `asymptotic_check.py`: independent factorial-identity controls and closed-form scaling checks through N=6144. Those sizes are formula evaluations, not many-body simulations.
- `REPORT.md`, `DERIVATION.md`: results, proof, supplied resources and limits.

Set `CP12_OUTPUT_DIR` to an existing directory to redirect generated JSON. Default script outputs go into this packet directory. The verifier normally uses a temporary directory and does not overwrite the saved records. `--record` explicitly refreshes them and writes VERIFICATION.json. MANIFEST.json records source and dependency hashes.

These are fixed-N derivative results and their separately derived asymptotic coefficients. They are not a certified spatial phase, finite-time preparation mechanism, uniform-in-N Taylor approximation, or propagation calculation. The auxiliary Q=2N control is not combined with main Q=3N results as though they belonged to the same phase.
