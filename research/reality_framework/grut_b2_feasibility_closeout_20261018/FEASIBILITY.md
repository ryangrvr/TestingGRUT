# Quantitative screen, not another dynamical reduction

All constants below come from the frozen B2-RA/RB/IS laws and thresholds.
The arithmetic is reproduced by feasibility.py; FEASIBILITY.json records
exact fractions and string-encoded large integers. Decimal square roots
are display values, not interval certificates of a physical observable.

## Route A: what a positive access bound would have to deliver

With output already bounded below by p_*=1/10, B2-IS gives
D_cond<=e/(2sqrt(p_*)). The 1/100 quality threshold therefore follows
from e^2<=1/25000. The positive crossing inequality is

    e<=2jp integral_0^T sqrt(P_X(s)) ds, jp=1/64.

Its integral must be <=sqrt(128/3125), about 0.202386. If a uniform
population majorant P_X<=zeta is used, a sufficient requirement is
zeta<=128/(3125 T^2). This is not an endpoint-population substitute for
coherence: it is a conservative amplitude bound integrated over the full
history, including worst-case returns. It is not a necessary gate.

The ONLY currently available size/time/receiver-uniform occupation bound
used in this screen is 0<=P_X<=1. It gives e<=T/32=20L. At L=2016
this is 40,320; at the smallest inherited whole-word record it is
2,684,354,880. Neither is quantitatively usable. Initial P_X(0)=0 and
zero instantaneous accepted current do not improve this full-window bound.

Keeping the signed output kernel could avoid this majorant. Its quantitative
phase-sensitive estimate is precisely absent; writing that integral again
would not resolve it. Furthermore information access alone leaves common
production, arrival and retention unproved. No adequate Route A was selected.

## Route B: absolute current control cannot decide the probabilities

The mass, electric, fuel, density and degree diagonals commute with R,T;
the receiver-only Hamiltonian also commutes because its final state is
unobserved. Geometry-changing norm is <=g+4epsilon<=21/10. The two-resource
constituent hopping norm is <=4(jp+jm)=5/48: at most two of each species,
each on a degree-two path, with capacity gates only removing edges.
The two enabled conversion terms have combined norm <=2J=1/1920.
Thus for the self-adjoint off-diagonal K,

    ||K||<=21/10+5/48+1/1920=1411/640.

For an orthogonal projection O, [K,O] has off-diagonal block norm
||O K(1-O)||<=||K||; no extra factor two is needed. Since the exchange
term annihilates R and T, this also bounds the common current operators
i[A_c,R] and i[A_c,T]. Their signed prepared-state matrix elements remain
unknown. From zero initial output/retention, absolute integration gives
only 0<=p_O(t)<=min(1,(1411/640)t), with no positive lower bound.

At T=640L the upper majorant is 1411L. Already at L=2016 it is 2,844,576;
at the first inherited whole-word length it is 189,381,236,784. Hence it
cannot certify output >=0.1, retention >=0.25 or the coherent cross entry.
These are limitations of this route, not physical violations.

For scale, the fixed output margin is 114089/4000000=0.02852225 and the
retention margin is 184281/1000000=0.184281. A sufficient unnormalized
success-map norm budget for the 1% information condition is
514089/404000000, about 0.0012725. No such moving map estimate is known.
Direct use of the exact instrument criteria also needs actual values;
it cannot turn the available trivial ranges into a decision.

## Size and computation: two distinct scopes

L=2016 is the formal minimum of the prepared FIXED-PATH B2-RA theorem.
It is NOT the whole accepted moving record when the inherited preparation
requires N>=2^28 and L>=N/2. The first compatible multiple of 24 there is
L=134,217,744, giving T=85,899,356,160. This is an eligibility intersection,
not a selected favorable record or a new preparation probability.

64N^2 receiver modes are at least 2^62. Dense one-boson complex128 storage
alone is 16*2^62=2^66 bytes, compared with 2^30 authorized. The actual
receiver contains L-1 bosons, not one. This rules out DENSE enumeration
within this envelope; it is not a memory lower bound for every algorithm.
Sparse special states, analytical kernels or certified compression might
avoid it, but no current moving-window error theorem supports such a
calculation here. No truncation, surrogate receiver or small-host trajectory
was substituted. Finite sampling would also require separate time/size
coverage before it could settle this task.

## Feasibility conclusion

No currently justified estimate reaches an operating bound within the
existing computational envelope: at most 25 seconds per child, 60 combined
child seconds, 1 GiB memory and 25 MiB new files. No new A/B research
execution was launched. The exact signed currents, exchange influence and actual
receiver correlations remain the specialist problem, not new definitions
to commission. Stop is a method/resource decision; the physical verdict
remains unresolved.
