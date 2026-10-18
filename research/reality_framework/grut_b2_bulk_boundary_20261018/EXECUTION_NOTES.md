# Execution disclosure

Charter freeze: ecff27de1be31a3a383b89ea63dfdee271f91d67.
Proof/code freeze: ca9ce9a00f834e686a7c68e74883b876517670f7.
Eight executed source/proof files were hashed before primary execution
and remained unchanged through the replay. No fallback or variants.

Python 3.12.14 executed controls.py once and once unchanged as replay,
standard library only. All word amplitudes and constants are Fractions.
The enumeration uses the actual constituent hops, potential denominators,
Q-only intermediate states and the actual quadrant restriction. It compares
bulk and quadrant static kernels, not probabilities of classical paths.

The declared cases are primitive lengths 2 and 4, source co-pairs n=0,1,2,
and energies 2 and 2-5J. Every intermediate charge separation is retained
for each finite word. The output archives all four bulk/quadrant tables.
The proof, not enumeration at a few sizes, controls every omitted word.

| Worker | Child wall seconds | Peak child RSS KiB | Result |
|---|---:|---:|---|
| Primary | 0.030584197000052882 | 10368 | 66 exact checks pass |
| Replay | 0.03127154799994969 | 10368 | Byte-identical output |

Combined child wall time is 0.061855745... seconds, including startup and
output writing but not assistant derivation or connector publication.
RSS is the operating system's child high-water mark. No ceiling was reached.

CONTROLS.json SHA-256:
`df0f895600d8a3318b7e25355a094d80467c72d8f5ab9d87f68f47a29a4cc91e`.
Logs record exit codes, stdout/stderr, timings, version and hashes.
MANIFEST.json covers final files except itself; verify.py checks that
archive, executed sources, replay and the frozen resource envelope.

No physical source evolution, time grid, real-axis quadrature, interval
spectral root certificate or operating probability was calculated. The
off-real Toeplitz factorization is an author derivation, not a numerically
validated limiting-absorption theorem. One assistant shared proof and code
ancestry; no independent external validation or reviewer contact occurred.
