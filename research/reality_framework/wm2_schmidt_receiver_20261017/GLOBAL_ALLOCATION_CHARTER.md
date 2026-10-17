# S8-G1: separately frozen global charge-rank allocation

Freeze before G1 implementation checks or trajectories. Original S8 charter:
0e42a5d3023c7ec45029f38ef918c47f62b6c603. The first original S8 component
has stopped at t=30.0654: one charge sector hit 16 while total bond rank
was only 55 of 96. The second original component retains its original
protocol. Neither original result is reclassified by G1.

## The one change

Allocate at most 96 Schmidt vectors in total across all charge blocks by
their globally largest singular values. Permit up to 96 in one charge
instead of reserving at most 16 separately per charge. The *total* bond
ceiling remains 96. Original S8's per-charge gate is deliberately changed
only in this separately named numerical protocol, not silently relaxed.

This is a test of rank allocation, not an assertion that the old method
passed, a physical-law change, an additional favorable disorder family or
an empirical prediction. SVD output remains a proposal checked with the
same positive discarded-norm, reconstruction and isometry accounts.

Keep every other S8 choice: physical receiver L_N=64N^2, exact ordered
histories, original root weights, neutral preparation, g=.9, +/-1/4
offsets, common detuning zero, exact photon charge, all reverse amplitudes,
declared computational cut and full-H error bound, boundary gates, field
audit, timestep algorithm, .60 full-state envelope and its component
ceilings, 2048-MiB soft RSS and 900-second fresh G1-worker budget.

## Controls, recovery and locked follow-up

1. Check the global allocation rule against an exactly constructed charged
   packet for which each conditional history needs one Schmidt vector but
   the joint state needs 18 or 33. This control is not claimed to be the
   actual prepared state or an observed formation result.
2. Run the full S8 independent controls with the global-allocation class,
   including every amplitude of N=3, L=576, T=36. Verify that the new
   numerical class differs only at rank selection from preserved S8.
3. Run the same N=5 components, (2,2) and (3,1), T=100, L=1600, weights
   (3/5,2/5). Both must complete every gate and their weighted lower bound
   for m>=3 must be >=1/2; m=5 is separate. Preserve any stop.
4. Fresh immutable S7 workers on those same two components provide a
   common measured baseline, using their unchanged 600-second budgets.
   Different protocols are not an equal-error speed comparison.
5. Only if the *G1* recovery gate passes, attempt exactly N=6, counts
   (3,3), T=144, L=2304, with original equal root weights and the proved
   mirror reuse. No probability at a stopped endpoint replaces T=144.

Original S8's recovery gate stays failed if either original component
stops. G1 has its own frozen gate. No additional allocation/rank/budget
variant will be launched in this packet after G1 outcomes. Mathematical
analysis may identify a next construction but cannot retroactively change
any numerical account. All-size formation and independent disorder remain
open even if this finite G1 test passes.
