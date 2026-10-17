# Exact coefficient/workspace counts are not whole-memory measurements

For the frozen six-event root component, remaining type counts are (2,3).
The number of exact type words at photon count k=0,...,5 is

\[
a_k=(1,2,4,7,10,10),\qquad
a_k=\sum_{r=\max(0,k-3)}^{\min(2,k)}{k\choose r}.
\]

For a common rank R, *all* charge-allowed receiver occupations give

\[
D(R)=\sum_{k=0}^5 a_k{R+k-1\choose k}.
\]

The explicitly constructed one-body-density workspace has

\[
F(R)=R\sum_{k=1}^5a_k{R+k-2\choose k-1}
\]

complex entries. These counts are exact combinatorics, not measured RSS or
a claim about the cost of a compressed history-conditioned representation.

| R | Coefficients D(R) | Density workspace F(R) |
|---:|---:|---:|
| 10 | 28,951 | 97,770 |
| 20 | 525,251 | 2,110,040 |
| 26 | 1,686,959 | 7,093,658 |
| 28 | 2,358,511 | 10,029,768 |

At rank 28 the density workspace alone is about 153 MiB in complex128.
It excludes coefficients, polynomial work vectors, orbital/occupation maps,
Python objects, cached neighboring layouts, Fock rotations and library
overheads. Measured whole-worker peak RSS belongs in REPORT.md and result
JSON. The native kernel's compact maps can reduce resident memory even
when its retained coefficient count grows; it cannot remove the complete
Fock-space cost.

Each native B action visits occupied-mode incidences (at most N-1 per child
coefficient), including both append and exact reverse contributions. Natural
orbital proposals require the density contraction. A sequence of two-mode
Fock rotations also touches the coefficients repeatedly. Thus complete action
enclosure does not make propagation free: density and basis work can dominate
after the cancellation-prone residual has been eliminated structurally.

There is no all-size rank, history-count, runtime or memory theorem here.
No coefficient ceiling, approximate eigenvalue or favorable observed rank is
substituted for the paid physical-state error account.
