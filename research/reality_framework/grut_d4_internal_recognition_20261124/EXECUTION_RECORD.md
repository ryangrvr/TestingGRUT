# D4 execution and preservation record

24 November 2026. The branch starts from D3 commit
`0f18c546ab0b0ad2f0f5968c52bc02349f2c01e4`. The research packet contains
the construction and separate internal AI critique.

1. The model, parameters, forecasts, control cases and required solver
   comparisons were frozen in `5c59a2dce5858752a785f563dd779fa1034a9371`.
   The design preceded finite numerical production. The two original
   hashes are recorded in [FREEZE.json](FREEZE.json).
2. The production solver evaluated 243 four-pulse timing rows across 81
   ordered strings and three waits. It evolved the exact ready comparator
   with the true physical fault maps and coherent output probabilities.
   The proof supplies an analytical bound for the writer-on system, using
   exact fraction comparisons. The large calculation did not instantiate
   the exponential full spin Hilbert space. [RESULTS.json](RESULTS.json)
3. The separate reviewer built physical tensor models of dimensions 144
   through 576 with the writer present. It compared the ready comparator
   to a separate 512-dimensional model in nine scenarios. The largest
   complete joint probability error there was 4.22e-15. In the small
   actual models the writer-on separation was checked against its
   Duhamel estimate. [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json)
4. The exact failures of simultaneous corruption, the identical-state
   collision, finite recurrence and occupied-bath backflow were retained.
   The two-copy and d-label information bounds were derived separately.
   [CONTROLS.md](CONTROLS.md) · [INFORMATION_BOUND.md](INFORMATION_BOUND.md)
5. The archived model replay agrees on 4,648 numerical fields at stored
   precision. Both checks were run from a copied folder outside the
   checkout using the installed dependencies. [VERIFICATION.json](VERIFICATION.json)

The calculation used Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 on Linux
with one OpenBLAS and OpenMP thread. The main gates are rational bounds,
while residuals between the separately implemented solvers are numerical.
External RAI and human peer review remain open. The change set contains
the D4 packet and its coordination entry.
