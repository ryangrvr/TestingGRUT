# Numerical-method ledger

Production monotone-sweep residual implemented; all prerecovery controls pass.
Full physical reference error .01359811 is inside .19647972.

The bounded recovery **fails** and this numerical branch is closed.
Exactly the two frozen workers ran; no fallback, new allocation, larger
budget or six-event trajectory follows. No target-time interval is
inferred from a stopped state. See RESULT_TABLES.md and summary.json.

* Counts [2, 2]: UNRESOLVED, model t=65.032680, wall 1202.430s, peak 68.305 MiB; whole-worker wall-time ceiling.
* Counts [3, 1]: UNRESOLVED, model t=76.746725, wall 1208.667s, peak 59.840 MiB; whole-worker wall-time ceiling.

The new account is mathematically checked and small-state validated,
but the direct residual implementation is too costly for this frozen
finite benchmark. Source pins, costs and failures are preserved.
The separate all-N analytic proof neither validates the numerical
method nor fills the missing five-event target. S5 remains the completed
finite-formation frontier; S9 adds a distinct asymptotic derivation.
