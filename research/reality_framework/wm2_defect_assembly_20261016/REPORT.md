# WM2-S1: generated adjacency and a protected end-defect mechanism

WM1 reference: `8741b6e74d571e643e9208187f5d93b5a8e2471c`.
This is a new WM2 construction, not another WM1 audit or a change to CP-12.

**S1 generates event adjacency without fixed pair affinities and
protects it against every virtual assembly excursion. Its closed
seeded dynamics fail appreciable formation.** The one-dimensional
class is also supplied by the single-end grammar. The full WM2
success theorem has therefore not been achieved.

## The explicit construction

The primitive registers are directed bond bits, one head bit per
event, free-event bits and one neutral probe species. A four-qubit
term appends any free event to the current head; its adjoint pops
that event. Every label pair has the same amplitude. Probe hopping
uses only currently occupied bonds. No coordinates, ordering of
event labels, supplied pair adjacency or Euclidean distance occurs
in the Hamiltonian. MODEL.md gives the exact operators.

The entire component reached from a bondless, permutation-symmetric
one-head seed is the set of all distinct-label partial words.
No intermediate geometry is dropped. The energy operator U K
counts free events. Its zero-defect sector contains every full
ordering, and an isolated low spectral band is derived separately.

Within a recorded full ordering, protected primitive transitions
define a path metric: degree at most two, diameter N-1, linear
ball growth and exact segment additivity. Which labels become
neighbors is generated; the path-like class follows from the
supplied grammar and is not a dynamically selected dimension.

## Protection is constructive, including changed final geometry

The prefix lemma proves, for every source ordering and every final
full ordering, that a virtual transfer over initial protected
distance d requires at least d+2 primitive moves. Activating a new
bond spanning d>1 also requires an intermediate defect count at
least d: the diagonal activation energy reaches at least U d.
The latter barrier follows from the dynamics, not a primitive
energy term containing the measured distance.

All append branches are counted by an exact operator identity.
With the supplied epsilon/sqrt(N) normalization their combined
kinetic norm is at most 2 epsilon, independent of the factorial
number of words. This normalization remains an input.

For U=16,epsilon=1,J=1/4, the complete Q-resolvent yields, for
every E<=1/2 and d>1,

\[\|R_j^P H_{\rm eff}(E)P_{i,d;j}\|_2\le
\frac1{13N}\left(\frac5{31}\right)^d.\]

Here the source projector includes arbitrary coherent full-geometry
packets with separation at least d, and the receiver includes
**all** final full geometries. This is stronger than a bound on
one selected configuration pair. The exact same tail applies to
the Schur correction at d=0,1. Direct neighboring hopping remains:
each same-word nearest effective element has magnitude at least
1/4, and arbitrary neighboring geometry packets have receiver
action at least 1/4-1/(13N).

The protected spectral gap is at least 13. Full-H propagation
also has a distance-weighted Dyson bound after summing every
partial and full receiver geometry. Thus the protection statement
does not rely on deleting virtual states, an arbitrary nonlocal
eigenbasis, or an unweighted global norm alone.

## Exact and numerical checks

* 2,364 finite-register configurations agree with the independent tuple transitions, reversibility and label covariance.
* 40,318 initial-source/final-protected-state pairs pass complete configuration-graph and Q-only shortest-path checks through N=7.
* Every virtual walk through length 10 is counted; no fixed-length distant shortcut appears.
* 7,830 long-bond configurations pass a complete minimax activation-barrier calculation.
* 115 receiver/source operator blocks pass the full Schur bound through N=6, retaining all final geometry records.
* Small full spectra check the gap and bare/exact-projector distance; full-H source pulses have variance 1/4 and agree with the +/-pi/2 identity.

These finite checks support the analytic all-N proofs; they are
not substituted for scaling arguments. Floating-point matrix
checks are not outward-rounded certificates or external review.

## The neutral preparation fails the formation gate

The preparation is the original CP-12 neutral ray, tensor the
explicit bondless one-head seed and probe vacuum. In the declared
product benchmark embedding, no energy flows from CP-12 into S1.
The old law is complete and unchanged. This embedding does not
make its original charged matter follow S1 geometry.

The seed has auxiliary energy U(N-1) and width
epsilon sqrt((N-1)/N). Conservation of the residual norm gives
the all-time formation bound

\[p_{\rm full}(t)\le
\frac{\epsilon^2(\sqrt{N-1}+1)^2}{N U^2(N-1)^2}.\]

Assemblies containing at least half the events have a separate
all-time upper bound that also vanishes as 1/N^2. This failure
therefore persists without requiring perfect full assembly.

| N | Full assembly: all-time upper bound | At least half assembled: all-time upper bound |
|---:|---:|---:|
| 8 | 0.0001324488064 | 0.001407854417 |
| 12 | 5.012818951e-05 | 0.0004869319633 |
| 16 | 2.576602289e-05 | 0.0002439974178 |
| 32 | 5.479258555e-06 | 5.18109603e-05 |
| 64 | 1.228307651e-06 | 1.1985431e-05 |
| 128 | 2.848336089e-07 | 2.884918227e-06 |

At fixed epsilon/U these bounds vanish; parameter retuning within
the proved gap regime U>2epsilon+4J still cannot rescue the full
assembly scaling. The exact radial reduction is checked against
small full Hamiltonians. Tiny sampled probabilities at large N
are not reported as resolved physical signals: the analytic
all-time bound is what establishes the failed gate.

The record includes all partial words, all successful words and
all roots. Its mean energy change and outcome entropy are supplied
in RESULTS.json. A particular full ordering has probability
p_full/N! under the uniform-root preparation. Conditional probe
tests are never promoted to an appreciable preparation success.

## What is supplied and what remains open

| Requirement | S1 result |
|---|---|
| Explicit finite-support, permutation-covariant primitive law | CONSTRUCTED |
| No fixed spatial adjacency or coordinates | SATISFIED |
| Generated event adjacency and protected operational distance | CONDITIONAL CONSTRUCTION |
| Dimension selected beyond the assembly grammar | NOT EARNED; path class supplied |
| Constant-length virtual shortcut exclusion | PROVED in the full reachable S1 class |
| Complete distance-weighted virtual sum | PROVED, including coherent source/all-record receiver blocks |
| Nonzero nearby transport | PROVED |
| Appreciable formation from declared neutral seed | FAIL, all-time scaling obstruction |
| Common geometry for original CP-12 matter or finite-density species | NOT ESTABLISHED |
| Full WM2 phase, continuum, gravity, distinctive empirical prediction | NOT ESTABLISHED |

As a sensitivity control, a different six-qubit block-reconnection
primitive permits rotation -> probe hop -> inverse rotation in
three steps across growing old S1 distance. Adding it changes the
protected metric itself; one must recompute that metric rather
than preserve the old distance by definition. Finite term support
alone does not imply a spatial configuration geometry.

## Next construction obligation

The prefix mechanism can be retained as a proved protection module.
The autonomous closed S1 variant is not a full WM2 success. A new
formation mechanism must exchange or carry the released binding
energy while preserving the primitive prefix barrier. Its reservoir,
fuel or preparation schedule must be explicit. Resonant channels
can remove the gap, so the existing Schur theorem cannot simply
be copied into that extension. No dissipative rate, cooling bath
or favorable low-energy state has been silently added here.

No general WM2 no-go is claimed. WM1, frozen CP-12 and SCOUT-1
are unchanged. No earned quotient value or empirical prediction
has been added.
