# GRUT-D7 · Learned phase compensation preserves quantum information

25 November 2026. D7 uses the internal sign record from D6 to
configure a physical compensator before loading a new payload.
The closed model transfers one half of a Bell pair from its
input site to a receiver through a sign-dependent phase
disturbance. The other half is isolated from the actuator.
The same Bell projector scores both hidden conditions.

With the weak assay and shallow phase error, calibrated
compensation yields **0.895284707521** entanglement fidelity.
The best physically allowed *blind* recovery is **0.853553390593**.
A full-size correction chosen from the same noisy record gives
0.875. The learned bit is useful even after allowing the blind
comparator any downstream quantum operation at the same interface.

| Weak assay, balanced source | Shallow phase π/8 | Deep phase π/4 |
|---|---:|---:|
| Calibrated correction with intact record | 0.895284707521 | 0.75 |
| Best blind recovery | 0.853553390593 | 0.50 |
| Full correction with intact record | 0.875 | 0.75 |
| Fixed replacement output | 0.25 | 0.25 |

At the deep phase, averaging the unknown sign without a record
completely dephases the transmitted qubit. The weak record raises
the Bell fidelity to 0.75, above the separable-state ceiling of
0.50. Thus it preserves entanglement across a physical link
under the model's access conditions. [Proof](PROOF.md)

## What the exact bound settles

For any accessible programme states $\rho_b$ conditionally
independent of the payload, define
$C=\sum_b\pi_b e^{-2i(-1)^b\phi}\rho_b$.
The optimum over all downstream recovery channels acting jointly
on Q and that programme is $(1+\|C\|_1)/2$. The physical phase
correction in the intact arm reaches this value for the tested
family. Keeping the original quantum probe coherent would not
raise the optimum here. The decoded classical bit already contains
the information needed for this recovery task.

The weak probe identifies the correct sign with probability 0.75.
Applying the full inverse for that sign overshoots in the shallow
case. The calibrated adjustment is $\tfrac12\arctan(1/2)$,
whereas the disturbance magnitude is $\pi/8$. The stage gains
4.1731 percentage points over optimal blind recovery and 2.0285
points over the full correction. The stronger assay reaches
unit fidelity. The magnitude, prior, actuator law and route are
supplied. Only the sign is acquired. [Model](MODEL.md) ·
[Proof](PROOF.md)

## Negative controls and timing

At the same shallow calibrated angle, scrambling, disabling or
using an independent passive selector gives 0.816227766017 after
averaging both source values. Inverting the record gives
0.737170824513. These same-angle ablations expose which input
the actuator uses. The best blind baseline remains the stronger
0.853553390593, which avoids applying a statistically unsupported
correction. An intentionally open direct-source bypass reaches
1 with its full correction and disabled memory.

The receiver measurement counts every trial. For the shortest
clock, the payload is ready with probability 0.922992943237 at
the window flanks, below the frozen 0.95 threshold. At L≥16
the complete-window bound on unfinished release is
$6.5536\times10^{-14}$ or less. The exact recurrence at 2T
returns the Bell pair to its preparation location and reduces
the receiver score to 1/4. The packet keeps both clock failures.
The interaction geometry, clean resources, Bell preparation and
Schrödinger time remain supplied. [Audit](CRITIC.md) ·
[All rows](RESULTS.json)

## Reproduction and physical status

The design and thresholds were committed before production at
`55217a742cb97978c713323e9cf90c6b23281aa0`. The frozen checks pass across the 648
reduced rows. A separate implementation
builds the local physical transition matrix without using the
spin reduction. It propagated 80 conditional rows over 320 coin
blocks. The 800 comparisons agree within $3.56\times10^{-15}$.
The energy-law residual in that implementation is below
$9.55\times10^{-15}$. [Freeze](FREEZE.json) ·
[Reproduce](README.md) · [Direct runs](ORACLE_RESULTS.json)

The comparison is experimentally measurable with the three
Pauli correlations XX, YY and ZZ on the receiver and reference.
The values are predictions of the specified quantum model.
There has been no hardware run. A gate-based test would check
the compensation law. Realizing the supplied exchange
Hamiltonian simultaneously would separately test its timing
and autonomous execution.

D7 establishes a supervised mechanism by which information
about an internal coupling can improve a physical function.
Entanglement fidelity and error correction with information
about the environment are established physics. GRUT-specific
physical novelty remains unestablished.
[Schumacher](https://arxiv.org/abs/quant-ph/9604023) ·
[Gregoratti and Werner](https://arxiv.org/abs/quant-ph/0403092)
