# GRUT-I3 — Environmental predictability and the physical necessity of memory

**Verdict B.** In the screened single-port assembly model, a physical
record raises bond formation to .759906 at a deadline where every
screened no-record controller falls below .70. The comparison grants
the null arbitrarily fast port allocation and cost-free latching.
An additional bounded cue increases its ceiling to .669591.
The model includes acquisition, actuator transfer, thermal leakage,
finite fuel and ideal reset. Standard kinetics, decision theory and
stochastic thermodynamics give the entire effect.
[Verdict](VERDICT.md) · [Frozen protocol](REALITY_TEST_PROTOCOL.md)

## Reliability and causal controls

A fuel-consuming receptor observes a bistable reagent upstream. The
reagent can switch during transit, while the record can decay. A
screened reversible transfer sets one catalyst orientation before
reaction. High forward and reverse rates require a match between
the orientation and the reagent. A physical material bond is scored
at the fixed end time.

| Primary policy or bound | Deadline bond probability |
| --- | ---: |
| Predictive record, optimal locked response | .759906 |
| Disabled sensing, passive record or shuffled readout | .478354 |
| Inverted readout | .196801 |
| Optimal fixed pure orientation | .478354 |
| Optimized time-dependent open loop, relaxed port | .640458 |
| Optimized material feedback, includes free latch | at most .644591 |
| Same feedback with extra cue of TV .025 | at most .669591 |
| Free perfect current cue, locked orientation | .839075 |

The proof guarantees at least .716434 for the primary predictive
arm before exact propagation. Finite-state replay gives .949099 match
after acquisition, .907175 before transfer, and .890263 for the
delivered actuator. After charging both fuel preparations, the
material and orientation preparation, the error bit, and the common
overhead allowance, the conservative discard-and-reprepare budget
is $9.946752\,k_BT_b$ per trial. The stock contains one sensor and one
assembly molecule. Ideal predictive device reset costs $5.877154$
jointly or $6.443609$ with separate resets, in thermal units. Lower
yield arms can spend less. The comparison establishes a reliability
advantage under a common work cap.
[Model](MODEL.md) · [Proof](FEASIBILITY.md) · [Thermodynamics](THERMODYNAMICS.md)

## Robustness and scaling

The moderate noise row reaches .720035, above its no-record ceiling
of .662519 including the weak cue. Severe noise brings predictive
success to .557014 and defeats the reliability target. Fourteen of the
30 switching-rate and delay grid points meet .70. At the primary
switching rate, extending the delay from one to four time units gives
.686933 and misses. Faster switching erodes the decision contrast.
All conditions and arm outcomes are preserved.
[Negative results](NEGATIVE_RESULTS.md) · [Results](RESULTS.json)

For isolated compartments with independent environmental paths, the
probability of at least 180 bonds out of 256 is .984510 using records
and at most .141058 for the no-record policies with the declared
independent weak cues. The rate and work budgets are held fixed per
compartment. Each compartment needs its own fuel capacity and port.
This is independent Bernoulli concentration.
[Scaling derivation](DERIVATIONS.md)

## Interpretation and review

The acquired information is a necessary capability in the sense of
the frozen reliability test, with screened access, a single port, a
fixed energetic affinity and a fixed deadline. The same physical
prediction can reside in another register or an actuator orientation.
Supplying the present conformation bypasses the need for the past
record. The apparatus overhead allowance and the screening assumptions
remain to be validated for a real system.

Two implementations agree across 2,624 scalar checks to
$6.28\times10^{-13}$ or better. The second implementation and
critical audit were prepared by a separate internal reviewer.
External RAI review and laboratory validation are pending.
[Reproduction](REPRODUCTION.md) · [Internal critique](INTERNAL_CRITIC.md)
· [External review packet](RAI_REVIEW_PACKET.md)

Consolidate D1–D7 and I1–I3 as conventional effective physics, and close
the developmental series as a route to fundamental novelty. The
finite-fuel comparison of I2 and the positive screened-access result
of I3 address different physical tasks. The comparison illustrates
when prediction earns its cost under explicit constraints.
[Novelty audit](NOVELTY_AUDIT.md)
