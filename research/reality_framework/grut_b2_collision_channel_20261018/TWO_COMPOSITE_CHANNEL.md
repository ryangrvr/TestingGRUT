# Two-composite law, collision observable and information channel

The law and all coefficients are pinned B1. B2 adds a frozen experiment,
not a new microscopic Hamiltonian. The actual numerical experiment is the
eight-event FIXED-RECORD matter/fuel reference, not full moving geometry.

## Same complete matter law

Enable B1 conversion at both ends, x=0 and x=7. Initial fuel colors are the
two input qubits, initially arbitrary rho_ABZ; all charges and E are vacuum.
The exactly conserved sector F+N_+=F+N_-=2 retains every field-compatible
charge configuration and both fuel bits. Capacity-one links exclude only
states absent from the microscopic physical space, not virtual channels.
Positive field energy, density and hard-core constraints remain present.

For a scalar configuration (p,m,f), Gauss fixes link field to the cumulative
charge sum. At most two plus carriers, and two neutral fuels at the ends,
carry the qubits. Embed just these COLOR-CARRYING registers on the ordered
line fuel_left, event0,...,event7,fuel_right. This is bookkeeping for the
fixed record, not a proposed primitive spatial metric.

Every plus hop moves a color to an empty adjacent plus register. Conversion
moves a color between an end fuel and its adjacent endpoint plus register,
also requiring/creating the minus there. Fuel and plus are each hard-core.
Therefore two color carriers cannot pass on this extended line. All
Hamiltonian terms preserve their left-to-right color order exactly.

In this ordered basis the full colored matrix is

    H_colored = H_scalar tensor I_AB.

This is a new two-resource statement with an explicit no-passing proof,
not tensorization of B1's single-composite spatial dynamics. H_scalar still
contains ALL collisions, absorption and correlated returns. Opposite-charge
constituents and which minus forms a pair need not preserve an artificial
particle label; output composites are identified by co-location and order.

## One joint outgoing channel

At theta_f=6, R requires empty fuels and co-pairs at x in {1,2}, y in {5,6}.
Decode their colors in left/right order, independent of their unmeasured
spatial history. The joint successful instrument is exactly

    E_R(rho_ABZ)=p_R rho_ABZ.

Failure is an orthogonal erasure with weight 1-p_R and reference marginal
rho_Z. Conditional information error is ZERO, for arbitrary two-qubit
coherence or reference entanglement, if this model's spin symmetry holds.
There is no multiplication of matter success by an independently estimated
readability. Collisions change p_R and spatial states but cannot entangle
the ordered qubits in this restricted endpoint architecture. A decoder that
instead labels particles by continued direction would be a different map.

## Collision diagnostics, with no hidden measurement

The true contact term V_c=(1/2)sum_x d_x d_(x+1) costs 1/2 on adjacent
co-pairs, much greater than J. A repulsive encounter may reflect before
large physical contact occupancy. At theta=3 inspect C_2: co-pairs at
(2,4) or (3,5), the near-encounter shell with one empty event between them.
Its probability is a separately prepared-ensemble diagnostic, NOT a
projective measurement inserted in the outgoing protocol or a joint
probability of a classical trajectory.

Compare the same final R with H_ref=H-V_c. This control changes no
one-resource matrix entry and cancels only a specific two-composite
contact energy. It is not an independent-gas model: field constraints,
hard-core collisions and all other interactions remain. A resolved final
probability difference demonstrates sensitivity to this explicit contact
term. It does not by itself yield an asymptotic S matrix or identify the
history of each outgoing trial. This finite preparation-and-readout task
has no monochromatic incoming scattering state or asymptotic free limit.

## The geometry-extension obstruction

If the right fuel becomes an interior side register because the geometry
extends beyond its recorded event, the extended carrier graph branches.
No-passing on a line no longer applies. On three consecutive occupied
events 0,1,2, with fuel at 1 empty and co-pairs at 0 and 2, these ten allowed
primitive operations exchange the plus colors:

1. plus 0->1; 2. minus 0->1; 3. recombine into fuel at 1;
4. plus 2->1; 5. minus 2->1;
6. plus 1->0; 7. minus 1->0;
8. convert fuel into co-pair at 1;
9. plus 1->2; 10. minus 1->2.

All intermediate fields are capacity-one, all charges obey Gauss, and
every coupling is a nonzero B1 term. The local word is enumerated by the
frozen control. This proves that the endpoint ORDER-IDENTITY mechanism
cannot simply be inherited after branching. It is local color exchange,
not superluminal signaling or a failed universal locality bound.

More generally with exactly two colors, color-blind primitive motions act
on sorted input colors as I or S=SWAP. The full operator has form
H=A tensor I+B tensor S in a declared configuration encoding. Geometry
and receiver states belong to A,B, not to discarded noise. Set H_+=A+B,
H_-=A-B and U_+=exp(-itH_+), U_-=exp(-itH_-). Then

    U=U_a tensor I+U_b tensor S,
    U_a=(U_++U_-)/2, U_b=(U_+-U_-)/2.

For a pure prepared spatial/environment state psi and color-blind successful
projection R, put a=RU_a psi, b=RU_b psi. The actual successful channel is

    E_R(rho)=||a||^2 rho+||b||^2 S rho S
               +<b,a> rho S+<a,b> S rho.                 (1)

For mixed conditional receiver states sum these Gram coefficients over a
purification. They include coherent exchange histories. Success generally
depends on the input singlet/triplet content:

    p_R(rho)=||a||^2+||b||^2+2 Re(<b,a>) Tr(S rho).

Thus one cannot multiply an input-independent arrival probability by
perfect readability for a growing two-fuel geometry without computing or
bounding these terms. In the fixed endpoint path B=0 and b=0 exactly,
recovering the identity/erasure channel above. No full-geometry value or
uniform error bound for (1) is established by the eight-event experiment.
