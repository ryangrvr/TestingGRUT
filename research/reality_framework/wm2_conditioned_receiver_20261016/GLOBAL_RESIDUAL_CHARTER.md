# S6-G1: separately frozen cumulative residual allocation

Freeze before G1 trajectories. Original S6 charter commit:
e627bb5f8287efb47d2ab4c79efbad6fd2984132. Preserve the original S6
charter, every original trajectory and every per-step quota stop.
An S6-G1 completion is not a pass under the original local quota policy.

The original S6 disordered runs stop under their local residual quota
while the global .20 residual allowance is mostly unspent. THEORY.md's
Duhamel/telescoping theorem requires only the cumulative sum, not an
equal share for each time step. Test that theorem without the additional
per-step quota. This is an openly revised numerical allocation policy,
not a new Hamiltonian or an improved full-state error envelope.

Keep the same exact history-conditioned representation, contact and
neighbor enrichment, projections, all reverse transitions, physical
L_N=64N^2 receiver, original neutral preparation and root weights, N=4,5,6
balanced +/-1/4 offsets, common detuning 0 or 1/2, g=.9, T_N=4N^2 and
macroscopic/full assembly events. No favorable initial packet or time
is chosen. The success probability target remains at least one half.

Keep the total state-error envelope .60, physical-H Strang <=.15,
word-conditioned projection <=.20, complete projected-emission residual
integral <=.20, and other contributions <=.05. Only replace the rule
"candidate step residual <= unspent residual / remaining steps" by
"candidate step residual <= unspent residual". The full candidate-step
cost and all preceding costs are still paid. Stop rather than cross .20.

Retain rank ceiling 10 per history, total coefficient ceiling 2,000,000,
cached dense overlap ceiling 8,000,000 entries, degree 8, three residual
coefficients, and 180-second step-boundary compute ceiling. All other
implementation choices are unchanged. Floating-point arithmetic is
not outward-rounded certification.

Reuse completed original S6 components. Rerun every unresolved primary
root component from the original seed under this policy, without choosing
by observed intermediate formation probability. Reuse only the same
proved sign identity. Unknown roots keep zero lower bound at their original
weight. Archive both ledgers, the original data hash, every target-time
stop and all error accounts.

Finite completion remains separate from improvement over S5, efficient
complexity, all-pattern independent disorder, or an all-size formation
theorem. None may be inferred from a cumulative-budget success.
