# What the constraints enforce: charged position transfer, to all orders

**Author-derived in the declared class; independent proof review pending.**
This is not a universality theorem for all matter observables. The internal
counterexample in INTERNAL_OBSTRUCTION.md lies inside the same class.

## The existing witnesses, computed as operators

For a species of nonzero charge z, an undressed term b_sj^* b_si changes
G_i and G_j by opposite nonzero amounts. The old bare all-pairs and star
terms therefore fail the Gauss constraint. Adding any constraint-preserving
local hopping term cannot cancel their commutators.

Single-link dressing b_sj^* U_(i->j)^z b_si fixes Gauss but fails C on
an unoccupied bond: starting with E_ij=0 it creates |E_ij|=z while
B_ij=0. Thus **Gauss alone is insufficient**. The additional supplied
flux-support constraint is essential. Again adding an invariant local
term leaves the nonzero C commutator unchanged.

Let P be the exact physical projector. For this elementary transfer,

    P T_s(i->j) P = P B_ij T_s(i->j) P.                          (T1)

This is a consequence of the constraint, not a derivation of the
constraint. The primitive implementation in CONSTRUCTION.md uses B_ij
explicitly and commutes with the constraints; a global projector is not
a permitted cheap Hamiltonian gate.

| Old witness | Gauss alone | Gauss plus flux support | Plus nonzero local hopping |
|---|---|---|---|
| Bare normalized all-pairs | Excluded | Excluded | Still excluded |
| Bare normalized root star | Excluded | Excluded | Still excluded |
| Single-link dressed all-pairs | Allowed by Gauss | Off-bond part excluded | Still excluded |
| Single-link dressed root star | Allowed by Gauss | Off-bond part excluded | Still excluded |
| Long occupied-path Wilson transfer | Gauge invariant | Allowed on that path | Requires as many link registers as its path length |
| Internal two-species interaction | Allowed | Allowed | Survives; full gate fails |

The finite-matrix controls check both charge values, with unphysical
electric configurations retained when computing the commutators. They
do not infer invariance merely from working in a zero-constraint basis.

## Complete words and virtual geometry

Consider an accepted record of length m, retain its first ell labels,
and source support at its first two positions. A receiver label has
recorded root distance D with 2<=D<ell. In the declared particle-number
sectors, each species is absent or present once, and has charge 1 or 2.

As long as the prefix has not been erased, its cut equations hold on
the original edges. If an A carrier crosses a cut, its charge change
is +/-1; B's possible change is 0 or +/-2. These cannot cancel. If B
crosses, its +/-2 similarly cannot be cancelled by 0 or +/-1 from A.
Consequently moving either carrier from an initial position i in {0,1}
to D changes at least D-i original electric registers, regardless of
the other carrier's endpoints, internal states, or coherent amplitudes.

A primitive geometry-diagonal term touching at most r electric registers
cannot change more than r of them. A gauge-compatible append/pop changes
none. Configuration-preserving D changes no electric values or positions;
receiver motion remains exact. Therefore any such primitive word must
contain at least ceil((D-1)/r) insertions. Alternatively, any history
that erases the first ell labels first needs m-ell+1 end pops. This
includes histories that restore the prefix later. Thus for

    k=min(ceil((D-1)/r), m-ell+1),                               (T2)

every interaction-picture word of order n<k has zero source-to-receiver
block, including all intermediate storage, electric, geometry and
receiver configurations. This is an operator-block statement for coherent
packets, not a count of one shortest route. In the concrete model r=1.
For any fixed r and D,ell,m proportional to N, k grows linearly in N.

A Wilson transfer along d occupied edges indeed evades first-order
single-link exclusion, but acts on d electric registers. It cannot be
a bounded-r primitive at growing d. Sequential Wilson segments need
growing word length. Next-nearest motion for r=2 is compatible with a
shared metric at a fixed range factor; strict nearest-bond hopping is
not claimed for every fixed-r class.

Suppose the complete changing interaction obeys ||V||<=W uniformly,
as required for the allowed class and explicitly proved for C4. The
entire Dyson expansion gives, without discarding any physical path,

    ||R_D exp(-itH) X_source||<=min(1,tail_k(W|t|)),
    tail_k(x)=sum_(n>=k) x^n/n!
             <=exp(exp(mu)x-mu k), mu>0.                         (T3)

The norm of the total interaction pays for path proliferation and all
interference. At |t|<=k/(4eW), T3 is <=exp(-3k/4). An extensive
receiver/diagonal norm is not substituted for W.

There are no virtual excursions into constraint-violating states because
P H (1-P)=0 exactly. This statement is a feature of a supplied physical
state space, not energetic suppression by a gap. Virtual configurations
**inside P**, including all storage errors and geometry returns, enter
T3. For a resolvent version, if eta>W and z=E+i eta, the complete
convergent expansion about the configuration-preserving D gives

    ||R_D (H-z)^(-1) X_source||
       <=(W/eta)^k/(eta-W).                                     (T4)

Each D resolvent preserves the same support count and has norm <=1/eta.
This is a positive-resolution bound; no eta->0 limit, real-energy gap
or uniformly bounded real-energy Schur operator is earned.

## Two distinct viable matter sectors

On every occupied nearest bond, the elementary constrained matrix element
for species s is exactly -j_s when the shift is within the physical kernel.
It always is for an allowed neighboring position of the declared carriers:
all physical fields have magnitude at most 3. Here j_A=1/64 and j_B=1/96,
with different charges and internal dimensions. The constraint does not
select these coefficients, equal speeds, or relativistic dispersion.

There is also a size-uniform nearby-response statement. In a single-species
experiment, inject the gauge-invariant adjacent pair of position/flux
configurations, with the same conditional storage and receiver field:

    (|r, E=0>+|w_1, E=the charge-z first-bond string>)/sqrt(2).

Both include the root compensator and the same internal level. A local
phase on n_sr has variance 1/4. With the other species absent, C3 vanishes
and ||V_s||<=g_max+2epsilon+2j_s<1.7. D gives identical receiver/diagonal
evolution to the two initial columns. The first cross term at receiver
w_1 is j_s t; all remaining terms are bounded by the complete double
Dyson tail:

    chi_(s,w1)(t)>=j_s t-[exp(3.4t)-1-3.4t].                    (T5)

At t=10^-4, T5 is >9*10^-7 for both species. This is uniform over the
original g/storage disorder and conditional receiver fields, for m>=2.
It can be prepared by root charge-neutral injection and a local first-link
rotation; no favorable remote spectral packet is required. The A/B
charges, internal spaces and coefficients remain supplied.

## Scope that cannot be dropped

The no-cancellation step uses at most one carrier of each positive charge.
Two charge-1 carriers could cancel a charge-2 transfer across a cut, and
neutral composites or gauge-neutral species need a separate analysis.
This is not a finite-density universality theorem. Already within the
declared two-carrier sector, unconstrained internal observables provide
the direct counterexample in INTERNAL_OBSTRUCTION.md. T3 constrains
position transfer, not the whole matter response algebra.
