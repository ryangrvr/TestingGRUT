# S5-P1: separately frozen probability-certificate extension

Freeze before P1 trajectories. Preserve the original S5 charter and
every S5 result, including rank/error-budget stops at the target time.
Those stops remain UNRESOLVED under the original state-error envelope.
This extension changes a numerical certificate, not the microscopic law
or the physical formation test.

## Why this extension is opened

S5's full-state error target 0.20 is sufficient but unnecessarily strong
for a nearly certain event. For normalized physical/numerical states
with state-norm error <=e and computed event probability p_hat,

    p_physical >= max(0, 1-[sqrt(1-p_hat)+e]^2).

This follows by applying the triangle inequality to the complementary
event. Combine it with S5's original reverse-triangle lower bound.
The theorem is derived before the new trajectories. It does not assume
that an uncomputed final probability equals an intermediate one.

## Physical target is unchanged

Use the same original uniform-root seed, physical L_N=64N^2 receiver,
all type-word histories, no extra photon cutoff, full reverse dynamics,
macroscopic threshold max(2,ceil(N/2)), and T_N=4N^2. The primary finite
families remain balanced +/-1/4 offsets at N=4,5,6, bar_delta=0 or 1/2,
g=0.9. All original root weights are retained. A finite positive result
still requires a lower bound at least 1/2 at that fixed target time.

Reuse completed S5 instances without rerunning them. Attempt every
unresolved primary root component; no choice based on its observed
formation probability. Mirrored components may reuse the proved sign
identity. The secondary families remain S5 controls, not new P1 targets.

## New error envelope, frozen before P1 runs

Allow total analytically accounted state error <=0.60, with product
formula contribution <=0.20 and receiver projection contribution <=0.35;
remaining analytic/normalization/basis-change contributions <=0.05.
Use the maximum of zero, the original reverse-triangle probability
bound and the complementary-event bound above. Keep rank ceiling 24,
coefficient ceiling 2,000,000, degree-8 exponential and per-instance
120-second step-boundary compute budget. Ordinary double precision is
not outward-rounded certification.

An S5-P1 pass is not an S5 state-error-envelope pass. Report both ledgers.
No observation time, Hamiltonian parameter, disorder realization,
geometry event or success threshold is changed to obtain a pass.
Unknown root components contribute zero lower bound at original weight;
they are never dropped and renormalized away.

Large-N all-pattern formation, continuous-time guarantees, nonshrinking
all-size robustness and an emergent spacetime phase remain open unless
actually proved. Small successful instances are not an all-size proof.
