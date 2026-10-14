# Reproduce the equal-energy CP-12 preparation study

From a full checkout containing this directory:

```bash
python -m pip install -r research/reality_framework/cp12_energy_matched_preparation/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_energy_matched_preparation/verify.py
```

The verifier runs the pinned parent min-entropy checker and all three new calculation scripts afresh. Outputs are regenerated in a temporary directory, then audited for exact certificates and cutoff differences. Saved tables alone do not constitute a pass. The run uses about tens of seconds on the recorded environment; hardware varies.

- `exact_check.py`: independent unregulated full-Hamiltonian rational certificates through T⁴, N=4,Q=4.
- `preparation.py`: Q=3N initial derivatives through T², N=8..256, with independent exact factorial moments and numerical quotient controls.
- `quantum_check.py`: complete finite-cutoff Hamiltonians, N=4,Q=4, L=1..5, full pure and mixed propagation at T=0,0.4,2,8,20.
- `REPORT.md`, `DERIVATION.md`: results, proof, resource accounting, and limits.

Imports from sibling `cp12_ground_connectivity` are byte-pinned in the scripts and manifest. The parent checker is likewise pinned. No working-directory-specific import assumption is needed; `python -I` works.

To redirect generated files, set `CP12_OUTPUT_DIR` to an existing directory. Default output locations are this packet's directory. JSON contains machine-readable exact coefficients and numerical diagnostics. Raw historical files and frozen constructions are not edited by verification.

The Q=3N calculations are not 256-event full many-body time evolutions. The four-event simulation is an auxiliary control in a different connection density, never evidence for a high-density spatial phase. Matching mean energy does not match the complete energy distribution. The two pure masks additionally match variance; their dephased comparison does not.
