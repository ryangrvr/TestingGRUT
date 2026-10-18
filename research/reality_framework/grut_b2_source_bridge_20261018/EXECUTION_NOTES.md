# Execution disclosure

The charter was published at 09fef5b1459f25e9380a052fac7bd74715199f4d.
The proof, code and source hashes were published at
dabee3a04a42b0b5236437e92f843d002a3b1b41 before primary execution.
No frozen source changed after either run. No fallback or second measure.

Python 3.12.14 executed controls.py once and then once as replay. It used
only the standard library. Fraction calculations sum all 225 declared
positive weights exactly, checking moments, central weight, continuous
bound constants, band margins and the source-to-joint arithmetic.
The continuous time statements and exact recurrence rest on their
elementary analytical proofs, not a floating point grid.

Three two-dimensional saturation controls check the capture formula and
nonzero cross terms. A single floating point recurrence evaluation checks
the exact integer-phase result. Floating point numbers are diagnostic;
no outward-rounded general matrix certificate is claimed. The diagnostic
return modulus 1.0000000000000018 is ordinary rounding of the exact one,
not a physical probability exceeding one.

| Worker | Child wall seconds | Peak child RSS KiB | Result |
|---|---:|---:|---|
| Primary | 0.03997951900009866 | 10368 | 40 checks pass |
| Replay | 0.042305081999984395 | 10368 | Byte-identical output |

The total child wall time is 0.082284600... seconds. Timing includes child
interpreter startup and output writing, not human/assistant analysis or
connector publication. RSS is the operating system's child-process high
water mark. Both workers stayed far inside the frozen limits.

CONTROLS.json SHA-256:
`a3a3d675db13a811509f292c9218dcda0bbc1c271c59391290ccdda55de5a7d5`.
PRIMARY_WORKER.json and REPLAY.json archive exit status, timing, RSS,
Python version, stdout/stderr and digest. EXECUTED_SOURCES.json pins all
eight frozen source/proof files. MANIFEST.json covers every final file
except itself; verify.py checks it and primary/replay agreement.

No source matrix, wavefunction trajectory, favorable endpoint search,
geometry formation experiment, actual source violation, paid computation,
installation or background process ran. One assistant shared all roles;
there was no independent implementation, reviewer or external validation.
