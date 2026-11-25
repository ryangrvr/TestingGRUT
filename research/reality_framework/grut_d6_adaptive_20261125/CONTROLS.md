# D6 controls and failed alternatives

All production arms use the same Hamiltonian, Hilbert space, initialized
head and tape, witness, dump, disable sink, and two independent mixed
coins. The arm's stationary intervention F and mobile port μ are
different known settings prepared independently of the condition. The
three strength settings S change the internal phase gate or close the
port. The full initial energy distribution is identical at fixed n
in all arms, including the physically different information paths.

No forward response gate makes a measurement or receives a timing pulse.
An evaluator scores both b values and both independent stage challenges.
A full word score has no per-cell averaging. A missed arrival or incomplete
bank counts as failure when it is required for a score.

## Weak probe at the terminal time

| Port and intervention | Pre-intervention record accuracy | Retained M accuracy | Both structures correct | Both structures and W forecast correct |
|---|---:|---:|---:|---:|
| Main intact | 0.75 | 0.75 | 0.75 | 0.75 |
| Main scrambled | 0.75 | 0.50 | 0.50 | 0.50 |
| Main inverted | 0.75 | 0.25 | 0.25 | 0.25 |
| Main disabled | 0.75 | 0.50 | 0.50 | 0.50 |
| Passive intact | 0.75 | 0.75 | 0.50 | 0.50 |
| Direct disabled | 0.75 | 0.50 | 1.00 | 1.00 |

The probabilities average b=0 and b=1 equally. Each result for the
intact arm holds for each b and each z. The disabled arm always writes
as if b were 0 and fails every b=1 structural trial. Its prepared
disable sink E receives the displaced memory. The scrambled arm XORs
the real memory with its own inaccessible stationary fair coin. The
stationary coin and other discarded records can hold correlations
that are absent from the response input. The intervention does not
erase global information.

Inversion leaves the accessible trace distance unchanged and reverses
the usual meaning of the pointer. The fixed response now scores 0.25.
A decoder informed of the inversion could recover the original 0.75.
Scrambling and disabling yield trace distance zero at the main port's
response cut. In the strong probe the intact, scrambled, inverted and
disabled structural scores are respectively 1, 1/2, 0 and 1/2.

The passive writer selects one independent fair bit and reuses it for
both tasks, even as M accurately registers the phase in the weak probe.
For equally likely source states its display ensemble is exactly the
same as the intact writer's. Correlations with the source and the
response to intervention separate the two. The direct-port writer
instead uses a literal copy of the source in A, and continues to score
1 through all memory interventions when that port is open. This bypass
uses an extra allowed interaction at the source. Equal energy laws
do not match information access. No such copy is loaded in the main arm.

## Indistinguishable and timing controls

With S sealed all interactions, including the direct port and witness,
are independent of b. The two complete accessible states coincide even
for another preparation of the head. Chance success (1/2) is optimal
for the two equally weighted structural targets. The closed-port
witness has a fixed Y outcome at T. A disabled writer that always
guesses 0 achieves forecast agreement of 1 against this constant
and scores only 1/2 on the unknown structural target. Forecast
agreement on this control supplies no evidence of acquisition.

The diagnostic A-only interval has a frozen 0.90 probability threshold.
It is missed for n=2, 4 and 8, with probabilities 0.660416206084,
0.755378503859 and 0.866130205885. The n=2 miss and the n=64 pass were
registered. Every intermediate size is retained. Ordered history
gates can be exact while the intended stage is weakly localized in
laboratory time. The full tape also returns to blank at 2T, which
excludes indefinite maintenance or retention in this finite model.

[Exact model](MODEL.md) · [Frozen gate definitions](PROTOCOL.json) ·
[Every conditional row](RESULTS.json)
