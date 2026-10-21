# Coarser timestep fails an accuracy control; reduce it before production

2026-10-21, before any production run. This corrects the incomplete wording in
Amendment 01 item 5: the first validation had TWO failures. Besides the mirror
initialization defect (fixed without changing its threshold), the homogeneous orbit
at dt=0.08 had max absolute field error 0.00200337 over duration 48 against a frozen
0.002 limit. Halving dt reduced the error by 4.00001. The error is numerical phase
accumulation in this second-order scheme, not a failure of the homogeneous solution.

The second validation passes exact mirroring but still fails that dt=0.08 control.
Both failed records are retained. Do not loosen the error limit. Before production,
reduce dt=0.08 rows to 0.06, the half-timestep control from 0.04 to 0.03, and validate
the homogeneous orbit at 0.06 and 0.03. The grid256 timestep stays 0.06 and grid384
stays 0.04. Physics parameters, ICs, seeds, endpoints and comparison criteria stay fixed.

Execute the two spatial refinements before the larger-box rows, following main seeds,
cold, kappa=4 and half-dt. Thus the cost cap cannot hide a grid problem behind a
completed larger-volume campaign. Unexecuted rows will still be reported.
