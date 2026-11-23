# Frozen control forecast corrected after separate review

The original [PROTOCOL.json](PROTOCOL.json) is preserved byte for byte.
Its maximally mixed control anticipated no new count record as well as no
new order. The order statement is right. The record statement is false if
the declared reference R is retained as an initial past-state tag.

Take one wall and one bath spin. Initially wall W and bath B are independent
fair bits and R=K=W. The active WB state is I/4, stationary under every
unitary, but RWB is **not** maximally mixed: R is correlated with W.
The exchange Hamiltonian at gt=pi/2 swaps W and B up to irrelevant phases
for this diagonal preparation. The bath then records K perfectly, while
the active WB state is still I/4 and the aligned probability is still 1/2.
Guessing K from B improves from 1/2 to 1, using the initially supplied
reference correlation. At intermediate time its success is

\[
 \Pr(B=K)=\tfrac12+\tfrac12\sin^2(gt).
\]

This is an unexpected failure of the frozen **control interpretation**, not
a change to the main reservoir-vacuum model, its equations, scored sizes,
finite-window thresholds or numerical outcomes. RESULTS.json and the separate
full-Hilbert check retain this counterexample. The full-state counterexample
rho_(RWB)=I/d would have no correlations to transfer, but would also have no
valid initial tag for a past-state decoding test. Those are different controls.

It is therefore legitimate to say "no new order for an initially maximally
mixed active whole," or "no new correlation for a fully maximally mixed
whole including every reference." It is not legitimate to infer the latter
from a maximally mixed active marginal with a correlated reference attached.
