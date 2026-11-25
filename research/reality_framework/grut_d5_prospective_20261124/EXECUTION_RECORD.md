# D5 execution and preservation record

25 November 2026. The D5 branch starts at the head of D4 draft PR #17,
`aa6ef7c1cfdc58ad941d8f37c413168cbd6acc02`. The construction and
independent physical check were separately tasked within Work.

1. Commit `9767724c9ee8302df9c29c8545744728989db136` froze the
   microscopic model, observables, sizes, controls and thresholds
   before numerical production. [FREEZE.json](FREEZE.json)
2. Production evaluated ordered histories for four sizes and every
   registered preparation. The external probe calculation grouped
   complete terminal work amplitudes, retaining nonarrival and
   incorrect outcomes. Failure probabilities use sums of disjoint
   case weights in logarithms. The sampled p values and model files
   retain their fixed design. The analytical certificate also applies
   the same tail inequality to p intervals around these samples.
   [RESULTS.json](RESULTS.json)
3. An independent implementation constructed the reachable physical
   graph with both forward and adjoint local hops, and included
   pulse-induced states at every head position. It directly propagated
   27 external pulse cases and 18 internal control cases at small
   sizes. Its failure case with the late external timing on the short
   clock is retained. [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json)
4. Both frozen files and every archived production field replay.
   All 4,656 numbers reproduce at stored precision in the recorded
   environment. Production and direct graph probabilities agree
   within 3.34e-15 on the small external cases. The complete packet
   also runs from a standalone copy. [VERIFICATION.json](VERIFICATION.json)

Runtime was Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Linux. Replays
set one OpenBLAS and OpenMP thread. Large-clock guarantees are
analytical. Direct graph comparisons use floating-point arithmetic
and deliberately tractable small systems. The independent referee
was a separate internal agent. Standalone external RAI and human
review remain open. Only the new D5 packet and its coordination
entry are in the change set.
