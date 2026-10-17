# Rank/resource interpretation; no gate changes

This clarification was recorded during G1, without changing either charter,
implementation, outcome or numerical gate.

At N=5, M=4 and the exact charges are 0,...,4: there are five charge
blocks. Original S8's cap of 16 per charge therefore implies an effective
maximum total rank of **80**, even though its additional total ceiling is
96. G1 removes the 16-per-charge restriction and can use all **96** of
its nominal total allowance. At N=6 there are six charges, so original
S8's per-charge caps could in principle use all 96.

G1 changes the admissible rank distribution and, at five events, the usable
maximum total rank. It is **not an equal-rank or equal-error performance
comparison**. The nominal total ceiling, wall/RSS ceilings and error gates
stay as predeclared, but actual memory and time can grow with concentrated
sector ranks. Report measured whole-process resources and completion status.

The observed original stops (total ranks 55 and 53) left capacity unused
even below the implied ceiling 80. The exact 18/33-history witness establishes
a sector-allocation capacity difference; it does not isolate how much of a
trajectory's improvement is due to allocation versus additional usable rank,
nor prove that any tested packet forms dynamically.

No further numerical allocation/rank/budget variant is introduced here.
SWEEP_BOUND.md is a separate analytical result and ordinary-double control;
it is not applied retroactively to any trajectory or stopped ledger.
