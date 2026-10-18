# Execution disclosure

Actual work: one assistant; no subagents, independent reviewer, external
messages, background job, persistent service, paid resources or experiment.

The candidate and gate were frozen at
`6b4f5478598a2f17eca97ed0d2dfaa0d9fc019a3`. The explicit proof, model,
control source and source hashes were frozen before primary execution at
`e11510191a373d9bfe5d3f74ea4ccebb02b98298`. No candidate, source, threshold,
time, decoder or fallback changed after seeing outputs. Pre-execution
analysis motivated the candidate; it is not a blind discovery claim.

Primary command: `checks.py` under timeout 50 seconds with OPENBLAS,
OMP and MKL threads all set to 1. It completed without error. Program
elapsed 0.312143865 seconds, peak RSS 71.58984375 MiB, whole process
0.584990163 seconds. A single replay under timeout 8 seconds completed:
program elapsed 0.274395186 seconds, RSS 71.67578125 MiB, whole process
0.478444920 seconds. Nonperformance parsed JSON compared exactly equal.
Total numerical worker wall time 1.063435083 seconds; resource limits were
60 seconds combined, 1 GiB peak and 100 MiB new files. No failed execution,
budget expansion, fallback or parameter search occurred. This replay is
same-code execution consistency, not independent mathematical review.

Environment: Python 3.12.14, numpy 2.3.5, scipy 1.17.0. The physical matrices
retain every allowed charge configuration in an exactly invariant resource
sector. Dimensions: 5,65,257 at one resource; 773 at two resources on eight
events. The full physical receiver and moving geometry were not simulated;
their comparison is analytic, pending independent proof review.

Scalar residual bounds and the Bessel arrival inequality were checked with
integer/Fraction arithmetic. Displayed decimal Bessel endpoints coincide
after float conversion; they are not machine intervals. The exact rational
inequality was evaluated before converting. Matrix diagonalization, norms,
phases and populations use ordinary double precision. Tiny negative
reported matter values at initial time (~4.4e-16) are subtraction roundoff,
left visible in raw output; no physical negative probability is inferred.
No full floating-point certificate or experimental uncertainty is supplied.

All 16 inherited source pins and six pre-execution source hashes remain
unchanged. `MANIFEST.json` covers the final packet except itself. `verify.py`
checks these hashes and accounting without repeating trajectories. No
canonical file or historical result was edited, and the closed numerical
recovery branch was not reopened. This bounded cycle has ended.
