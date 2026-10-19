# Closure-fork calculations and decision

Read [REPORT.md](REPORT.md), then [DERIVATION.md](DERIVATION.md) for the
mathematics and [SOURCES.md](SOURCES.md) for source/read scope. The findings
are classical channel gaps and explicit thermality/curvature qualifications.
This packet does not supply a new cosmological theory.

The follow-up [CLOSURE_SELECTION_RESPONSE.md](CLOSURE_SELECTION_RESPONSE.md)
answers the C2 question against the archived P1, supplies Obadia's full
reference, and corrects the proposed B control's source map. It introduces
no extra closure postulate and executes no new cosmological calculation.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/grut_closure_fork_20261019/probe.py
python research/reality_framework/grut_closure_fork_20261019/finalize.py
```

The executed environment is Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0.
The high-precision elementary calculations use Python's built-in Decimal at
60 digits. No package installation was performed. The complete probe takes
about five seconds on the execution host; it overwrites only its own
`results.json`. The initial dependency failure and coarse-epsilon failure
are retained with explicit amendments.

Expected outputs are a tree emission gap about `-.19999915 gv`, a tree
equal-split gap about `-3067.39734 gv`, quadratic midpoint residuals after the
first slow-roll temperature shift, a first-order endpoint correction, and
closed geometries with arbitrarily small negative Omega_K. The geometric
temperature check imports standard Euclidean-state de Sitter thermality.

`NEXT_QUESTION.md` states the unresolved cosmological selection question.
No prior report, coupling, state or gravitational matching parameter is
modified by this packet.
