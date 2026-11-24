# GRUT-D3 shows that a written record can select historical repair

24 November 2026. A system writes its current configuration into a local
memory under a fixed Hamiltonian. After a common reversible disturbance
erases the system's configuration, the very same Hamiltonian restores the
original configuration and retains its record with probability at least
**99.1568%** at any chosen time in the specified windows. Scrambling the
record cuts that joint probability exactly in half. The arms have identical
present active states and identical energy distributions even conditional
on each historical label.

D2 showed that a reservoir can passively record how many defects it
absorbed. D3 uses the same finite-chain transport mechanism for two
exchanges. The first writes a local pointer and the second uses that pointer
to select a repair transition. Both allowed configurations are treated
symmetrically. Either can be restored, according to the actual earlier
configuration and whether its record is accessible. [D2 report](https://github.com/ryangrvr/TestingGRUT/blob/675d991927edd520daa4b000680a15e8a5200fa7/research/reality_framework/grut_d2_development_20261123/REPORT.md)

## Exact model and measured contrast

S is a three-level site with two ordinary configurations and one erased
state. M holds either configuration or a blank. Two initially empty spin
chains absorb excitations emitted during writing and repair. A timed SWAP
with a prepared, inert register D erases S and carries away the historical
label. D scores the terminal outcome and has no coupling to the Hamiltonian
during repair. A fair coin, present in every arm, controls the scrambling
operation on M. [MODEL.md](MODEL.md)

The post-fault amplitudes evolve on the product of the writer and repair
paths with the two edges at the common origin removed. At that origin
S is erased and M is blank. Either conditional exchange can proceed
elsewhere when the required label is present. This exact finite reduction
keeps both reservoirs and has 174,785 amplitudes for one selected label.
The half-line comparison bounds the finite-path returns. A Duhamel bound
on the missing edges certifies the whole timing rectangle, including times
between numerical samples. [PROOF.md](PROOF.md)

| Arm | Joint historical recovery and retention | Value at the weakest sampled intact timing |
|---|---:|---:|
| Intact record | Above 99.15% throughout both windows | 99.956368% |
| Coin-scrambled record | Exactly half the intact value | 49.978184% |
| Inverted record | Zero historical successes | 0% |

The arrangement ensures that the probability of reaching *either*
ordinary configuration is identical across the arms. The gain is recovery
of the earlier configuration. Flipping M after the fault acts the same as
the global symmetry on the prepared state, because S has the erased value.
That symmetry commutes with H. It preserves every energy probability,
even conditional on the historical bit and the coin. Fair averaging gives
the same density for the entire evolving system plus coin. Scrambling
removes information about the past from the evolving system. The isolated
coin recovers the association when read with a nonblank record.
[CONTROLS.md](CONTROLS.md)

## Resources and failure cases

In the certified experiment, writing lasts between `640/g` and `656/g`
and repair lasts between `8/g` and `16/g`. The large model uses 2688
writer spins and 64 repair spins and starts both chains in vacuum. The
one common fault supplies mean work `Delta=4g`. All arms share the complete
conditional energy distribution immediately after the fault and control.
The fault register, ready memory, external pulse timing and stable local
interaction rules are explicit resources. The probabilities concern one
terminal read at any allowed pair of times.

Two controls delimit the result. With one spin per chain, a specified
timing restores S with certainty while returning M to blank. Perfect
restoration can consume the record. An unflagged flip of S leaves the
repair charge at zero, so the repair interaction never acts and the
historical configuration fails to return. The positive result is for
the declared flagged erasure. It leaves repeated repair and general
disturbance correction open. [CONTROLS.md](CONTROLS.md)

## Reproduction and scope

The model and protocol were recorded before finite-system production at
[`53fd7c9`](https://github.com/ryangrvr/TestingGRUT/commit/53fd7c9e13185e79f0983932f32807f543591e2a).
The bounds were derived during design. Production implementation and the
separate physical-matrix check followed that design record. The minimum
joint probability across the 27 sampled timings is `0.9995636769842553`.
Sparse exponential propagation and a separate Chebyshev expansion agree
within `1.62e-14` in state-vector norm across nine comparisons. The
polynomial truncation bound concerns exact arithmetic. Numerical solver
agreement checks the floating-point implementation. The finite-window
certificate is an independent rational calculation.
[PROTOCOL.json](PROTOCOL.json) · [RESULTS.json](RESULTS.json)

A separately tasked AI constructed full physical tensor Hamiltonians for
six small reservoir pairs. The checker compared 54 timing pairs across
those cases. The largest discrepancy across its tests is `4.69e-14`.
The complete conditional energy probabilities agree within `7.63e-15`.
Independent external review remains open. Replay reproduced all 27 time
rows and 435 numerical fields at stored precision. [CRITIC.md](CRITIC.md) ·
[VERIFICATION.json](VERIFICATION.json)

The model supplies a classical record that guides one restoration and
persists through the scored window. The family of configurations, the rule
that reads M, the fault flag and the preparation are built into the
experiment. It singles out no absolute configuration. Standard quantum
mechanics already permits autonomous correction and correlation-dependent
feedback. D3 gives a finite example with exact energy matching. Deriving
the local rule, an observer model or a cosmological selection rule
remains open.
[Cohen and Mirrahimi](https://arxiv.org/abs/1409.6759) ·
[Chapman and Miyake](https://arxiv.org/abs/1506.09207) ·
[Najera-Santos et al.](https://arxiv.org/abs/2001.07445)
