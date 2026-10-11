# CP09 finite autonomous clock-pointer apparatus

2026-10-11. Explicit measurement-apparatus extension, not a replacement fundamental law or CP10. Frozen CP09 remains unchanged. One declared apparatus with omega=g=1, N=2; no optimization or parameter scan. This supplements the ideal reference-record construction with time-independent joint dynamics. The coupling remains GLOBAL and supplied, so this does not solve locality.

## Hamiltonian and preparation

Let F be the projector onto geometry containing 0<1. At N=2 this identifies the nonisolated minimal/maximal roles 0 and 1. Supply a clock qubit c and record qubit r. In clock/system/record tensor order postulate

H_total = H09 + (I-X_c)/2 + |1><1|_c tensor F tensor (I-X_r).

All three terms are positive. This explicit apparatus coupling, clock Hamiltonian and pointer are extra assumptions. The clock gate is dynamical rather than externally pulsed, but external Schrödinger time remains in the calculation. Event-labelled F is a calibrated reference assumption, not an intrinsically localized or permutation-invariant observer. It is a coarse one-bit role record, not the entire f(C) register.

Prepare |0>_c tensor |ground09> tensor |0>_r. Initial total energy .5947552837. Full dense evolution of dimension 48 preserves norm and total energy to 1e-11. Minimum full energy .0947552837. Inverting the complete evolution at t=2 restores all subsystems with error 8.2e-16; reduced recording is not fundamental irreversibility.

## Energy and actual record information

At external t=2: system energy .16570827, clock energy .39202669, interaction energy .03702032. Their sum is the original total. System energy increases .07095299 while clock energy decreases .10797331. No unexplained external energy injection occurs during this time-independent evolution; preparing and measuring the apparatus itself has not been costed.

Condition on the declared clock measurement c=1, whose probability is .64636422. Joint role/record probabilities are

| | record 0 | record 1 |
|---|---:|---:|
| F=0 | .76531956 | .02415070 |
| F=1 | .15877763 | .05175211 |

The conditional mutual information is .04324332 nats. Role-record agreement is .81707167, but agreement alone would be misleading because F=0 is common; the nonzero mutual information explicitly establishes an imperfect informative record. These are joint endpoint measurements, not evidence that the record perfectly tracks a persistent earlier role while geometry changes.

At t=1 the conditional mutual information is .05111228 nats, so record quality is not monotonic. The pointer is reversible and has no dissipative memory stabilization.

## What this does and does not construct

This demonstrates a functioning, finite-energy autonomous record process with explicit apparatus-system energy exchange. It replaces an instantaneous ideal copy by a calculable interaction. It does not implement ideal relational source operations exactly or establish detector-response recovery; those claims are not credited.

The clock is a driven-backreacted two-state oscillator, not a monotonic time observable. Even an uncoupled qubit has recurrent readings, so conditioning on c=1 does not define a unique elapsed time or establish invariance under a change of clock. The supplied global F access, clock initial state and interaction remain physical costs. No local oracle, spacetime light cone, continuum scaling or gravity is derived.

The next substantive apparatus construction needs reference accessibility and a clock with a declared resolution and usable interval, with back-action included. Merely calling an ancilla a clock is insufficient. This calculation retains that failure rather than asserting that external time has disappeared.

Run python research/reality_framework/cp09_autonomous_apparatus/check.py (NumPy/SciPy). Expected CP09_AUTONOMOUS_APPARATUS_PASS. Results retain all energy components and joint tables at t=0,.4,1,2. Author finite quantum control; not a new GRUT law or experimental evidence.
