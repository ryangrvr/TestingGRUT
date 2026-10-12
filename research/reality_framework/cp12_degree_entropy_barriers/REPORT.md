# CP12 follow-up: degree defects and component entropy are separate barriers

2026-10-12. Based on research commit 7be83fc1f3d6a5799f11a587d3df5abc023dda4f.
No frozen Hamiltonian changes. The first switching argument below is supplied
by the owner's independent review; the enumeration independently checks it.
All probabilities in this note refer to labelled simple undirected graphs
with Q=N edges and classical weight exp[-beta lambda W], W=sum_i binom(d_i,2).
They are NOT the full CP12 quantum Gibbs state. Finite-degree defects do not
by themselves rule out approximate locality or a geometric phase.

## 1. The supplied finite-temperature bound is correct

Every 2-regular graph has W=N, the minimum at fixed sum d_i=2N.
Choose an oriented occupied edge (a,b) and a vertex c distinct from a and
not adjacent to a. Replacing ab with ac yields degrees d_b=1,d_c=3 and
leaves all other degrees two. The increase in W is exactly one.
There are 2N(N-3) operations per regular graph. In the result b and c are
identified by their exceptional degrees. A reverse operation chooses one of
at most three neighbors a of c. Thus there are at most three regular
predecessors, counted with operations, per resulting defective graph.
Writing R_N for the regular count, there are at least
(2/3)N(N-3)R_N distinct such defective graphs. Keeping just those two
energy shells in the partition function gives

    P(regular) <= [1 + (2/3)N(N-3) exp(-beta lambda)]^(-1).

This argument applies for N>=4; the N=3 graph is the single triangle.
For P(regular)>=p>0 it is necessary that

    beta lambda >= log[(2/3)N(N-3) p/(1-p)].

Hence beta lambda must grow at least as 2 log N + O(1) to retain a fixed
nonzero probability of exact regularity. For P(regular)->1 it is necessary
that beta lambda - 2 log N -> +infinity. Neither statement is sufficient;
other defect sectors have not been counted in that bound.

## 2. Even perfect regularity does not select one extended ring

Conditioning on regularity makes every graph's energy exactly lambda N.
The conditional distribution is therefore uniform, independently of temperature
and penalty strength, including size-dependent choices of both.
A labelled 2-regular graph is a disjoint union of undirected cycles of length
at least three. Its exponential generating function is the standard identity

    sum_N R_N z^N/N! = exp[sum_(l>=3) z^l/(2l)]
                     = exp(-z/2-z^2/4)/sqrt(1-z).

The single spanning-ring count is C_N=(N-1)!/2. Coefficient asymptotics give

    R_N/N! ~ exp(-3/4)/sqrt(pi N),
    P(one ring | regular) = C_N/R_N
                          ~ exp(3/4) sqrt(pi)/(2 sqrt(N)).

Thus P(one ring)<=C_N/R_N ->0 for ANY nonnegative sequence beta_N lambda_N,
including the zero-temperature limit of this degree-only model. Cooling can
remove defects; it cannot break the exact energy degeneracy among cycle
partitions. This is a counting result, not a new gravitational principle.
For fixed finite beta lambda the two bounds combine to

    P(one ring) <= (C_N/R_N)/[1+(2/3)N(N-3) exp(-beta lambda)]
                = O(N^(-5/2)).

This is an upper bound, not an asymptotic equality. It concerns a single
spanning ring, not the absence of macroscopic cycles or of a giant component.

## 3. Executed independent controls

check.py exhaustively enumerates every Q=N graph for N=4 through7, verifies
the regular counts, all forward switch counts, maximum reverse multiplicity,
and the Gibbs probability bound at beta lambda=0,1,5,10. It independently
counts regular graphs through N=500 using the component containing vertex1:

    R_N = sum_(l=3)^N binom(N-1,l-1) (l-1)!/2 R_(N-l), R_0=1.

The exact conditional ring probabilities at N=20,100,500 are approximately
.433017, .188794 and .084009. The asymptotic ratio at N=500 is 1.001252.
Run `python check.py`; results.json is the actual output. No CP12 quantum
thermal simulation was performed for this note.

## 4. Physical interpretation and the next falsifiable target

The positive fourth-order CP12 term supplies degree regularization in a
conditional weak-hopping comparator. It does not supply connectivity or
spatial dimension. The earlier exact loop-energy result adds a different
obstacle: the first tested topology-sensitive conditional coefficient favors
two triangles over one six-cycle. These statements have different scopes;
neither can be substituted for the full quantum free energy.

A meaningful next target is the effective free-energy difference between
one extended component and competing cycle partitions after integrating out
matter and electric fields, retaining rewiring-induced mixing. It must include
both their statistical multiplicity and their physical propagation. A finite
penalty for each disconnected component would itself be an additional physical
interaction, not a consequence of the degree term derived so far.

Exact degree two is too strong to be a universal requirement for physical
space. A successful spatial phase may carry finite-density defects. Therefore
these results reject exact regularity / spanning-ring selection by this
specific degree-only ensemble; they do not reject all emergent geometry.

## Literature and attribution

Cycle-set enumeration is standard labelled combinatorics. See Manuel Bodirsky,
Combinatorics lecture notes, generating-functions treatment of 2-regular graphs:
https://wwwpub.zih.tu-dresden.de/~bodirsky/Combinatorics.pdf .
The generating function above is rederived explicitly; no originality claim
is made for it. The CP12 interpretation and combination with the owner's
switching bound are analytical follow-up work, not empirical confirmation.
