# Reproduction

From the repository root at the F5 draft head, with Git and Python 3:

~~~
git merge-base --is-ancestor 1177b1690d6118f81c4f95706ef0ce3a9bda4ec2 HEAD
git diff --exit-code 1177b1690d6118f81c4f95706ef0ce3a9bda4ec2 HEAD -- \
  research/reality_framework/grut_f5_closure_consistency_20261127/FROZEN_PROTOCOL.md \
  research/reality_framework/grut_f5_closure_consistency_20261127/PARAMETERS.json
cd research/reality_framework/grut_f5_closure_consistency_20261127
python3 verify.py
~~~

The production program propagates signed kinematic basis states using the dressed shift, extracts their clock sections and evaluates the joint effects, phase transport and reduced densities. The independent program builds the histories directly from the cumulative gate table and relative coordinate, reconstructs the signed cycles, and separately evaluates the effects and controls. Both use exact integers and rational arithmetic. The verifier replays both archives byte-for-byte and checks every registered value without tolerance or external dependencies. The universal finite-prefix statement and finite-dimensional probability criterion require the analytical proofs. The tests audit a single finite witness and its declared controls.

| File | SHA-256 |
| --- | --- |
| FROZEN_PROTOCOL.md | 5acdc075ad38e7cc11456f48d83fa8618bf7c24aedd19aa4bc07a0f9c795b9a1 |
| PARAMETERS.json | 16426776b9239325a90bb3a31fffdd0e9e137c28e076d78d6217d015faad941a |
| RESULTS.json | 9a69db0c422be3227606cc8edf5369d883d9a3a5b849fe860cd822fa8878dbff |
| INDEPENDENT_RESULTS.json | ed3aec34c0f74530f9662c7f630e2f05b69716588446dfed0a9025c109d06dd6 |

The two outputs agree on all 31 registered quantities. Successful replay includes the expected negative outcomes, such as a failed B chart on the stopped branch, identical matter predictions for the distinct lost-reference seeds, and an empty invariant sector when the wrap is twisted. External replication is pending. [Proof and derivations](DERIVATIONS.md) · [Separate internal review](INTERNAL_CRITICAL_REVIEW.md)
