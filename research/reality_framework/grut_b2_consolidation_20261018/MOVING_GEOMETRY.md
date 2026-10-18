# Composition with the actual moving-geometry law

**Verdict:** the existing persistence/comparison hypotheses do not match
the two-endpoint B2 setup. No moving-geometry operating pass or physical
failure is proved. This is a compatibility diagnosis, not a successor law.

## 1. The actual law being compared

Retain the archived B1 extension of the P1/S9 geometry law, not an
invented graph. A configuration is a rooted simple word w of length m
among N labels, with membership a_i, storage q_i, root R_i and a physical
64N^2-mode receiver. Its relative diagonal is

    16(Q-m)+64 sum_i(q_i-a_i)^2
      +sum_i delta_i(q_i-R_i)+dGamma(h_receiver), Q=sum_i q_i,
    n_receiver=Q-1, |delta_i|<=3/4, .9<=g<=1.1.

The receiver single-particle hopping is -1. Assisted append of a free,
empty child v raises q_v and creates a contact boson with coefficient
-g/N; reverse pop is its exact adjoint. Both directions have zero-field
and empty CHILD/removed-HEAD gates, including fuel in B1. The old head
need NOT be empty for a forward append. Bare append has coefficient
-(n_head+n_child)/(4 sqrt(N)), changes neither q nor receiver, and retains
its adjoint and the same support gates. Here n counts matter, not fuel.
All storage errors and receiver states remain in the physical law.

Matter, three-level fields, mass, density, fuel and conversion are exactly
the B1 terms. In particular the degree term is

    -A sum_i [2-deg_current_word(i)] n_+,i n_-,i.

The program enabling conversion is static on selected labels. Select the
two recorded endpoint labels, prepare their fuels, and do NOT make the
program follow a moving head. A following-head program would change the
law/protocol. Vacuum before fuel loading still has the old formation block,
but that fact does not establish operation after two fuels are loaded.

Initial host assumptions for an attempted integration are those of the
actual conditional record: m>=N/2, q=a, the actual receiver with m-1
bosons, no reset, vacuum matter/fields, and independent arbitrary input
colors in the two selected fuels. Recording and loading are supplied.

## 2. A well-defined common-space comparison

For a proposed B2 segment length L<=m, retain the first L recorded labels.
Let K retain that prefix and its storage bits, contain all matter/fuel
inside it, and have vacuum matter/fuel and zero flux beyond it. On K,
define the mathematical reference as the EXACT H_L of THEOREM.md plus
the genuine tail/storage/receiver dynamics reflected at length L. Omit
bare geometry and matter crossing the retained boundary in this reference;
use the endpoint degree -A at the retained right endpoint. On K-perp use
the corresponding compression of the actual Hamiltonian, with reference
cross blocks zero. This is a self-adjoint comparison on the SAME space.
The physical Hamiltonian is never cut or switched.

Use one common readout extension: success requires the original first L
recorded labels, bonds and prefix storage bits intact, both selected fuels
empty, and co-pairs at the same recorded X_L,Y_L. Decode their plus colors
in recorded left/right order and trace all actual tail/receiver registers.
Every other trial is erasure. The reference always retains this prefix,
so its success instrument is exactly the fixed-path instrument. For the
retention bound below, use the joint event of this prefix/storage survival
and total k=2, which includes separated charges. No intermediate measurement
or favorable receiver selection is introduced.

As in B1, a decoupled retained matrix would factor from the actual
conditional receiver. But the complete actual-minus-reference drive
must still be bounded. It includes prefix-breaking assisted pop, bare
end motion, charge hopping across the cut and the degree discrepancy.

The old proof controls these latter drives by resource occupation near
the cut. Its locally conserved resource is

    r_x=n_+,x+n_-,x+2n_f,x.

For support only at x<=b, retained length ell and boundary layer width 2,
it requires an EMPTY RESOURCE BUFFER

    d_buffer=ell-2-b>=0

and pays, for v=12(e-1), a term proportional to
exp(-d_buffer/2)[exp(v|t|/2)-1]. A separate assisted-tail penetration
term pays tau_(m-ell+1)(1.1|t|). No real-energy gap substitutes for either.

## 3. Endpoint placement and the buffer are incompatible here

For the exact B2 reference ell=L and fuels at x=0,L-1, b=L-1 already at
t=0. Thus d_buffer=-1. The support hypothesis fails BEFORE evolution.
Fuel cannot be dropped from r_x: conversion would then create resource
locally, invalidating the conservation step used in the old estimate.
Nor can the buffer be manufactured by choosing a smaller b while keeping
the actual initial weights. The right fuel alone gives
W(0)>=2 exp(L-1-b); with d_buffer=L-2-b this implies
sqrt(W(0)) exp(-d_buffer/2)>=sqrt(2e). The exponential suppression
therefore cancels in the inherited resource estimate.

There are two natural embeddings; neither earns inheritance:

| Embedding | Exact mismatch |
|---|---|
| B2 path is the entire accepted word, m=L | Right fuel sits at the movable head; there is no tail buffer m-L. It blocks pop while occupied but does not block forward append. |
| B2 path is an initial segment of a longer word, m>L | Right fuel is already an interior side register; its occupied outward bond and actual degree differ from the terminated B2 endpoint. The retained-cut resource buffer is still absent. |

Making N or the remote tail larger does not move the second fuel away
from the B2 cut. Keeping a much larger prefix instead puts that fuel in
the INTERIOR of the reference matter Hamiltonian. That Hamiltonian is no
longer H_L and requires its own source/channel analysis. These observations
do not refute an eventual new estimate for the unchanged full law.

## 4. Explicit local witnesses of the mismatches

### Interior embedding: a two-operation transfer across the B2 cut

Take any occupied record with m>=L+1 and b=L-1. Keep the actual receiver
state unchanged, prepare both selected fuels, and define the target with
the left fuel still excited, right fuel empty, plus at w_L, minus at w_b,
and E_(w_b,w_L)=-1. Gauss and the capacity-one constraint hold.

Exactly two matter operations reach it: right conversion (-J), followed
by outward plus hopping (-jp). Diagonals or geometry changes cannot
produce that matter target in two operations while returning to the same
word/storage/receiver state. Therefore

    <target|(H_full^2-H_ref^2)|initial>=J jp=1/245760.       (1)

The analogous outward-minus target has coefficient J jm=1/368640.
H_ref has zero such amplitudes because it preserves K. These coefficients
have no distant-head or N suppression. They are exact local identities,
not claimed lower bounds on a late-time probability. Separately, on a
co-pair at w_b with empty outside matter, the actual diagonal exceeds
the fixed-path endpoint diagonal by A: actual degree is 2, reference 1.

### Actual-head embedding: fuel does not prevent extension

If m=L<N and v is any free event, assisted append has matrix element

    -g/N sqrt(n_contact+1),                                (2)

even while the old head carries fuel. The new child is empty and its
new link has zero field, so every forward gate is satisfied. The old
fuel stays at its selected event, which becomes interior. This is
permitted dynamics of the old law, not a proposed modification.
If m=N there is no free child at that instant, but neither a tail nor
a resource buffer is thereby obtained; the old theorem still does not
apply. No appreciable full-assembly preparation is assumed to repair it.

### Information mechanism: exact order identity is no longer inherited

The archived B2 ten-operation color-exchange word is available when an
enabled fuel is an interior side register. On consecutive events 0,1,2
with empty fuel at 1 and co-pairs at 0 and 2: move both constituents
0->1, absorb; move both 2->1 then 1->0; emit at 1 and move both 1->2.
Every field remains capacity one; the plus colors are exchanged.

This invalidates the GLOBAL endpoint no-passing argument after extension,
even when the recorded prefix itself survives. It does not establish a
1% information failure from the specific B2 preparation in its window.
That would require controlling the actual identity/SWAP channel Gram
coefficients and their preparation weights, or a complete state comparison.

## 5. The operating time also needs a new comparison

The B1 integrated theorem used t<=15360 with N>=2^28 and a macroscopic
empty buffer. B2 requires T_L=L/(6J)=640L, already 1,290,240 at L=2016,
84 times that fixed horizon. A longer time is not covered by simply
reusing its conclusion. For illustration, if ell=L=N/4, then T_L=160N;
the older D1 size-scaled use window N/4096 is smaller by factor 655360.
Neither this timescale comparison nor saturation of a loose Dyson bound
proves physical destruction. They identify missing hypotheses/estimates.

## 6. Concrete margin for a FUTURE complete comparison

There is a precise, reviewable target without constructing another law.
Let V embed the fixed-path experiment and its actual conditional receiver
into K. A new proof would need a state comparison, uniform over all
allowed input/reference states and the stated host/receiver class,

    ||[exp(-it H_full)-exp(-it H_ref)] V psi||<=epsilon_geo,
                          0<=t<=T_L,                      (3)

where the reference retains genuine tail/receiver evolution as above.
Equation (3) must include color exchange and every geometry/storage/field
excursion and coherent return. It is NOT established in this packet.

An achieved epsilon_geo<=1/1000 would be sufficient. Applying the common
success projection and decoder to a purification gives

    p_out,full > (717/2000-1/1000)^2
                 =20449/160000=0.12780625,
    p_ret,full > (659/1000-1/1000)^2
                 =108241/250000=0.432964,
    conditional trace-distance error <=(1/1000)/(717/2000-1/1000)
                 =2/715<1/100.                            (4)

For the last line, if subnormalized successful vectors x,y differ by
<=epsilon, with reference norm >=sqrt(p0), their normalized pure-state
trace distance is <=epsilon/(sqrt(p0)-epsilon); partial trace cannot
increase it. This is uniform for entangled inputs by using the same
purification. It does not multiply independent readability/arrival events.

The margin is available; the comparison is missing. The actual local
witnesses show why the existing empty-buffer proof cannot furnish (3).
They are not a theorem that no future proof can do so.

## Decision and stop

B2-RA remains the archived author-derived fixed-path pass. Automatic
moving-geometry composition is NOT established; the exact blockers are
endpoint/buffer placement, exposed outward matter coupling and degree,
loss of the endpoint order identity, and the longer operating horizon.
No new law, protected boundary, parameter search or successor run is added.
The next specialist problem is (3) for the actual preparation, or the
corresponding actual two-source output instrument if that reference is
too different. It is not another scalar-source estimate or a gravity claim.
