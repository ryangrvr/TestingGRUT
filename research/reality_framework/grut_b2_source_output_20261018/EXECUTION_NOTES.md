# Execution disclosure

The source/output operating protocol was NOT numerically run. There are no
new trajectory populations, large-host simulations or tests at later times.
Only frozen exact-scalar and small embedding/derivative identities ran.

| Resource | Primary | Replay |
|---|---:|---:|
| Program seconds | 0.042380423 | 0.037424206 |
| Whole worker seconds | 0.405835826 | 0.234547917 |
| Peak RSS MiB | 44.531250 | 44.128906 |

Combined worker time 0.640383743 seconds. Runtime: Python 3.12.14,
NumPy 2.3.5. BLAS thread variables were all one. Each worker had a
20-second timeout; RSS and file ceilings were monitored, not OS-enforced.
No paid allocation, installation or persistent/background execution.

CONTROLS.json was byte-identical on replay. The computations validate
the explicit Fraction inequalities and small matrix identities at L=4,6,8.
Floating sparse matrices have ordinary double-precision tolerances; there
is no outward-rounded matrix certificate or hardware arithmetic contract.
The code reuses B1's exact catalog and local entries, so this is not an
independent implementation. The weighted graph itself and infinite source
scattering problem were analyzed, not numerically solved.

PROOF.md and controls.py were hashed and committed before execution.
Their hashes are unchanged. SOURCE_ESTIMATE.md was written subsequently
to express the still-unproved source conditions with unfiltered observables;
it changes no target and is not a successful numerical result.

One assistant performed all roles. There has been no external review,
empirical comparison or independent validation. No other agent or process
continues the work after the cycle. Historical records remain unchanged.
