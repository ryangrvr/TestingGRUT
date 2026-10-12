# Six-vertex execution boundary

The full N=6,Q=5 basis contains 2,210,625 states at electric cutoff L=2
and 3,180,993 at L=3. L=2 completed through T=20. The preserved L=3
checkpoint completed through T=10; its requested T=20 is unfinished.
The first memory-reduced attempt was killed with exit137; the subsequent
segmented attempt lost its running session. No unfinished result is reported.

At T=10:

| Cutoff | Attached | Connected | Spanning path | Cutoff boundary weight |
|---|---:|---:|---:|---:|
| 2 | .167996064 | .105845240 | .024956170 | .199875396 |
| 3 | .167894855 | .105835432 | .024944194 | .008802639 |

Selected geometry probabilities agree closely, but electric-square per vertex
changes from .811473937 to .831984329. These are not an infinite-cutoff
convergence demonstration. At T=20 the L=2 boundary weight is .207097994,
so that late-time result must not be promoted to converged physical evidence.

`partial_N6_Q5_L3.json` preserves all completed rows and diagnostics. Earlier
T=10 direct and segmented executions agreed at about 1e-14 on the displayed
geometry probabilities; this is an execution control, not cutoff control.
