# External RAI review packet

**Status.** Prepared for external review. The project includes a separate internal mathematical criticism and exact implementation. Outside replication and laboratory verification remain pending.

## Proposition to challenge

Physical closure is modeled by putting the clocks, matter, records and admitted readout apparatus inside one constrained quantum system. The probability rule and observable algebra are supplied. Test the candidate condition that different internal descriptions of the same operational event agree, even with interacting clocks and finite precision. The claimed obstruction is that this condition permits arbitrary finite matter histories and states within a specified cyclic class.

Take $a,b\in\mathbb Z_N$ for two clocks with $N$ an odd prime, a physical matter qubit $S$ with $s=0,1$, strides $\nu=(1,2)$, and $r=b-\nu_s a$. For each label set $V_{s,r,0}=I$ and define
\[
T|a,b,s,\psi\rangle
=|a+1,b+\nu_s,s,V_{s,r,a+1}V_{s,r,a}^\dagger\psi\rangle .
\]
Use the eigenvalue-one sector of the cyclic shift, equivalently the zero sector of its supplied spectral constraint. The physical state space is the image of
\[
J|s,r,\psi\rangle=
N^{-1/2}\sum_a |a,r+\nu_s a,s\rangle V_{s,r,a}|\psi\rangle .
\]
The reductions $U_A(a)=\sqrt N\langle a|J$ and $U_B(b)=\sqrt N\langle b|J$ are claimed unitary on their complete retained spaces. A change of ideal chart transforms both state and physical effects. Joint conditional probabilities use a separately transported numerator and denominator. The matter prefixes determine genuine differences at fixed local clock and record access, since removing the writer changes the observable distribution.

## Required challenges

1. Derive the spectrum, physical inner product, exact kernel and two charts. Check that each is a description of the complete physical state, that it retains the other clock, and that the anchors add no experimental conditioning to a matched event.
2. Challenge whether clock coupling and matter backreaction are physical within the supplied algebra. Check the signed phase minor at fixed matter. Identify which assertion would become mere conjugacy if the readout algebra were also allowed to change. Reproduce the no-write comparison with the inverse wrap updated by the same rule.
3. Verify the POVM, gauge average and conditional filter. Check for hidden postselection, lost coherences, a supplied external trigger and inappropriate multiplication of separately averaged effects. Establish exactly what a gauge-preserving pointer dilation supplies.
4. Derive the noise and access controls independently. At joint noisy reports 2 the filter and joint positive-record probabilities are $3/40$ and $1/20$. Separately conditioning on the same numeral gives different probabilities. Verify the signed phase effect, identical reduced states after losing B, and loss of injectivity for a stopped branch.
5. Test whether the obstruction survives for an arbitrary finite prefix list and whether any locality, energy, winding, capacity, constraint algebra or vacuum claim goes beyond that theorem. Check the sign-twisted wrap at fixed gauge character. Scrutinize why a different quantum cosmology or a new notion of closure could have additional physical input.
6. Compare against Page–Wootters, complete observables, group averaging, clock changes, interacting and periodic reference frames, history-state encoding and the multiple-choice problem. The key comparator is the same state and physical event in one ordinary Dirac theory with the same measurement access. Evaluate whether any genuinely independent observable restriction survives. [Novelty audit](NOVELTY_AUDIT.md)
7. Reproduce all frozen fixtures, document disputes and state an A–D verdict with the smallest valid scope. Keep any failed controls.

## Packet and checkpoints

- [Parent draft PR #27](https://github.com/ryangrvr/TestingGRUT/pull/27) · [Pre-production freeze](https://github.com/ryangrvr/TestingGRUT/commit/1177b1690d6118f81c4f95706ef0ce3a9bda4ec2)
- [Frozen criteria](FROZEN_PROTOCOL.md) · [Exact fixtures](PARAMETERS.json)
- [Foundations](MATHEMATICAL_FOUNDATIONS.md) · [Derivations](DERIVATIONS.md) · [Physical comparison](REALITY_TEST.md)
- [Reproduction](REPRODUCTION.md) · [Internal critical review](INTERNAL_CRITICAL_REVIEW.md) · [Negative results](NEGATIVE_RESULTS.md)
