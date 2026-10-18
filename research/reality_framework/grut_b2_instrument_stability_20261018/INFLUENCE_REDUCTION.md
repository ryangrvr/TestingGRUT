# Which boundary effects matter to the instrument?

**Tier:** author-derived structural reduction of the ACTUAL pinned law;
long-time estimates unresolved, independent mathematical review pending.
Preparation, readout and notation are fixed in INSTRUMENT_THEOREM.md.

## 1. Final environments, not invented history records

For each actual final spatial/geometry/receiver basis state alpha accepted
by R, the color Kraus operator is a_alpha I+b_alpha SWAP. Each coefficient
is the COMPLETE coherent sum of primitive histories ending in that state
with its corresponding permutation. Histories ending in the same alpha
are not made orthogonal by an artificial clock or measurement. Distinct
orthogonal final environment states disappear across a partial trace in
the ordinary way. In a general representation the same accounting is the
inner product of outgoing environment vectors.

Consequently a collection of histories is not harmless merely because
it changes unused labels, contains only color-preserving primitives, or
has a small/zero contribution to sector leakage. It can change later
production or interfere with other histories. Only the resulting full
Gram entries determine the specified instrument.

There is a precise harmlessness statement: a common final environment
isometry, independent of the input colors, leaves the instrument unchanged.
If the moving successful dilation is (I_color tensor W) times the fixed
successful dilation, then their operations are equal. W is an analytical
identification, NOT an extra decoder/resource. For this law's two-parity
description exact equality to p_f id is equivalent to u=v=c=p_f.
Large changes of the common environment are allowed by that condition.
No such equality has been proved for the actual moving experiment.

More generally, equality of the two outgoing Gram matrices is equivalent
to a common isometry between their two-vector spans. One sufficient
Hamiltonian situation is a common auxiliary generator commuting with
BOTH H_+,H_- and R. Then its unitary is the same W on both outgoing
vectors. A shared geometry term need not satisfy those commutators.

## 2. The only primitive that exchanges colors

Choose the carrier encoding specified in the theorem: root fuel before
root plus; static right fuel b after plus at b; all other fuel slots
vacuum. Then:

- Geometry appends empty registers and pops only empty registers. It
  never permutes occupied colors. Storage/receiver and diagonal terms
  are color identities; minus particles carry no color.
- On-event fuel conversion crosses no other occupied color register,
  because its plus slot must be empty. Its adjoint obeys the same fact.
- A plus hop across an ordinary consecutive word edge preserves order,
  EXCEPT across b and its CURRENT child while fuel_b is occupied. That
  hop passes the side-register color. Its adjoint also swaps the order.

This last point is the complete source of B_c in H=A_c tensor I+B_c
tensor SWAP. Let T_X be the scalar sum of the allowed crossing hop and
its adjoint, with unit matrix coefficients, including n_f,b=1, occupied
bond, Gauss shift and capacity/hard-core gates. Then

    B_c=-jp T_X, jp=1/64.                                     (1)

When fuel_b is occupied, the conserved two-resource sector has N_+<=1.
There is a unique current child of b. Each allowed crossing configuration
therefore has at most one partner, and the reverse hop recovers it.
T_X is a self-adjoint PARTIAL MATCHING, not a general extensive hopping
operator:

    T_X^2=Pi_X, ||B_c psi||=jp ||Pi_X psi||.                   (2)

Pi_X projects onto the configurations incident to an actually allowed
exchange hop. If b has no child, is absent, has empty fuel, or the field
gate forbids the hop, that configuration contributes zero. The invariant
inactive-fuel vacuum is essential; arbitrary extra internal fuels would
add exchange locations outside this experiment.

Initially both fuels are occupied and N_+=0, so Pi_X xi_w=0. Reaching
Pi_X requires occupied right fuel AND a plus at b or its current child,
with an allowed crossing. The other resource must already be matter.
This is a concrete prepared-state access problem. Initial vanishing does
not bound it over times proportional to L; no frozen-distance argument
is imported into a changing word.

## 3. Separate common production from relative information

H_+=A_c+B_c is exactly the usual scalar color-erased moving law: every
allowed primitive has its original coefficient, all receiver and geometry
states retained. H_- differs only by the signs of the exchange hops.
These are exact spin-symmetry blocks, not modified candidate laws.

Their complete difference satisfies

    (U_- - U_+)xi_w
       =2i integral_0^t U_-(t-s) B_c U_+(s)xi_w ds.             (3)

The actual SUCCESS-vector difference is the same signed integral with
R on the left. Define this output-restricted coherent kernel as

    K_X^R(t,s)=R U_-(t-s) T_X U_+(s)xi_w.                     (4)

Keeping the integral retains every return and cancellation. A sufficient,
potentially loose, restricted-access estimate is

    ||chi_--chi_+|| <=2jp integral_0^t ||Pi_X U_+(s)xi_w|| ds.  (5)

This is NOT a new absolute bound on global geometry motion. It is paid
only at the single exchange gate. Its right side is not bounded here;
replacing the access probability by one would make it useless on the
growing window. If (5) is too loose, the coherent integral (4), rather
than another global leakage estimate, is the exact information target.

With a proved output minimum p_*>0 and outgoing-vector difference e_R,
INSTRUMENT_THEOREM equation (6) gives D_R<=e_R/(2sqrt(p_*)).
Equations (3)-(5) do not establish production, arrival or retention.
Those require actual probabilities of R and T under this same moving
law, not the fixed-path band or absence of exchanges alone.

## 4. Zero instantaneous current is not harmlessness

For a color-erased state psi_s=U_s xi_w and an outcome projection O=R
or T, exact probability currents are

    d ||O psi_s||^2/dt = i <psi_s,[H_s,O]psi_s>.               (6)

A geometry primitive whose two endpoints both have the same final
acceptance status contributes no DIRECT current across O. Tail growth
preserving the record is such an edge. Growth of the initial two-fuel
state also leaves both endpoints outside success. But these transitions
can change later source dynamics; they are not discarded as harmless.
Retreat through a required recorded label, restoration of such a label,
matter routing, conversion and field gates can cross the final event.
The partition is by the actual measurement, not by the old P projector.

For clarity, a general same-space reference H_0 gives the exact identity

    U(t)* M U(t)-U_0(t)* M U_0(t)
      =i integral_0^t U(s)*[H-H_0,M_0(t-s)]U(s) ds,
    M_0(a)=U_0(a)* M U_0(a).                                 (7)

Here M can be the lifted decoded-output observable, including a reference.
This preserves signed expectations and future observable structure. An
auxiliary difference commuting with ALL these backward effects has zero
instrument influence. Commutation merely with M at time zero does not
suffice. The actual source defect from B2-RB shows why growing beyond
the right source cannot be treated as an independently evolving auxiliary
by assertion: the new degree and outward hop terms change future effects.

The actual exchange gate makes the warning especially concrete:

    R T_X=T_X R=0, and T T_X=T_X T=0.                          (8)

Exchange requires occupied fuel_b; R and T require both fuels empty.
It has ZERO direct current for both events, but need not have zero later
influence. Equation (4) contains R U_-(t-s) T_X, not R T_X; R must not
be moved through that propagation.

In fact the whole actual success Gram matrix G_ss'=<chi_s',chi_s> has
the common current operator J_R=i[A_c,R]:

    d G_ss'(t)/dt=<U_(s')(t)xi_w,J_R U_s(t)xi_w>, G(0)=0.      (9)

The same identity holds for the retention Gram with J_T=i[A_c,T].
The exchange terms vanish from the CURRENT OPERATOR but still alter
the two states and their cross expectation. Equations (4) and (9) are
explicit observable-level kernel targets, retaining their signs/phases.

## 5. Precisely unresolved quantities

An observable-level completion needs either:

1. The actual success Gram entries u(t),v(t),c(t) and retention values
   u_T(t),v_T(t), controlled on the complete frozen windows; or
2. A controlled common-mode moving success/retention calculation plus
   the output-restricted exchange integral (4), or another justified
   estimate of the same relative channel.

No raw sector-survival probability is a substitute. All quantities use
the ORIGINAL xi_w/sigma_w. No history postselection, extra tail, moved
fuel, geometry freeze, receiver reset or new physical parameter is used.
This cycle isolates the access/coherence and probability-current problems;
it does not solve their long-time values or trigger a successor campaign.
