# Is there a tractable analytical route?

**Decision: a concrete route is identified; the scaling operating statement
is mathematically unresolved. No B3 or further simulation is opened.**

The route uses an exact uniformly isolated *fixed-record matter band*, a
quantitatively mobile exact bulk band, and a boundary spectral-measure
problem. It requires no larger gap or change of the B2 law. Whether the
boundary and two-source estimates actually close remains a real risk.

## 1. A further analytic upper velocity estimate

This estimate was added after the derivation/source freeze and the identity
controls, as analytic reasoning only. It does not alter a prospective
threshold, source, time, law, or control run. No new trajectory is involved.

Use the exact fiber and positive even ground vector in DERIVATION.md.
With d=Delta-2b>0,

    ||u_Q|| <= sqrt(2)b u0/d,   u1<=||u_Q||/sqrt(2),
    <S+S*> = 4u0u1 + <u_Q,Q(S+S*)Q u_Q>
             <=4b/d+4b^2/d^2
             =4b(Delta-b)/(Delta-2b)^2.

The function (Delta-b)/(Delta-2b)^2 increases on the declared b interval.
Together with the earlier lower bound this gives, for 0<k<pi,

    (16928/8825) J sin k < e'(k)
        <=(1128/529) J sin k < (11/5) J.              (18)

The last rational comparison is 5640<5819. In particular the exact mobile
band's maximum bulk group speed is less than 2.2J, while momenta away from
0,pi have a size-independent lower speed. This is NOT a strict full-state
causal speed; the full model also has high spectral components and tails.

The prospective receiver/time rectangle has an intentional geometric margin:
velocities in [1.1J,2.2J] yield displacements between L/8 and 3L/8 throughout
t in [L/(8J),L/(6J)], with strict endpoint margins at least L/80 (lower) and
L/120 (upper), before integer rounding. These are ray displacements, not
probabilities for the prepared experiment. At the same times neither such
outgoing ray reaches the center. This suggests a **pre-encounter** route.
It would demonstrate dilute two-source operation, not scattering or a
repeatable interacting medium. No collision gate is retrofitted into B2-SA.

## 2. First load-bearing estimate: the actual fuel source

Define H_ch on the half-line as the exact one-pair charged Hamiltonian,
including its physical boundary degree counterterm and all separations.
Let b=|0,0>, fuel energy Omega, and conversion -J. There is no external sink.
For Im z>0, block inversion gives the exact boundary/fuel Green function

    m(z)=<b,(z-H_ch)^(-1)b>,
    G_f(z)=1/[z-Omega-J^2 m(z)].                      (19)

This is an identity, not a Markov or independent-emission approximation.
If the boundary spectral density rho_b(E) exists, then its absolutely
continuous contribution to the fuel spectral density is

    rho_f(E)=J^2 rho_b(E) /
       [(E-Omega-J^2 Re m(E+i0))^2
                           +(pi J^2 rho_b(E))^2].    (20)

Any source-bound eigenvalue off the continuous spectrum must also be
counted. Away from poles of m its condition is

    E-Omega-J^2 m(E)=0,
    residue=[1-J^2 m'(E)]^(-1).                      (21)

Here m'(E)=-<b,(E-H_ch)^(-2)b><=0. Thresholds, embedded modes, and singular
spectral parts cannot be ruled out simply by displaying the bulk band.

The outstanding source estimate must provide:

* Total prepared weight in outgoing channels with velocities fitting the
  frozen receiver, including any trapped/source residues in the complement.
* A uniform late-time control, including all source returns, to complement
  B1's short launch interval rather than extrapolate it.
* Quantitative finite-half-host and edge corrections through O(L/J) times.

Equations (12)-(18) make this a specific boundary scattering problem with a
known channel, rather than an unspecified new mechanism. They do not solve
it. Positive bulk velocity alone gives no lower bound on its fuel overlap.

## 3. Second load-bearing estimate: do two distant launches compose?

The needed comparison is to two one-source preparations before their
outgoing fronts meet, with the same physical observable dictionary. An
independent-particle ansatz is not a proof. In the physical B2 Hamiltonian
all hard-core constraints, field strings, pair contacts, both fuels, and
possible source recollection remain.

One possible proof strategy is an exact-band connected expansion: a term
that couples the two source regions must include a word spanning their
separation; use the all-length contour representation to control its weight.
Band-limited stationary-phase estimates can then exploit the velocity margin
in (18), rather than the much faster crude constituent propagation bound.

This strategy has an explicit open input: a connected comparison for the
exact band's input, generator AND output event. The contour distance bound
by itself is not that comparison. A naive product of canonical band maps
does not automatically identify the physical two-body band. Cutting a
contact can change which adjacent-pair states lie in the low band, so
equality of projector ranks must never be assumed across that reference.
Provide a common-Hilbert-space construction or a different justified map.

High-band norm can be separated once, uniformly in time; charging a fixed
high-band norm at EVERY time step would produce a useless O(L) bound.
Conversely, omitting high components inside a Duhamel integral without a
new justification is invalid. All omitted effects need a controlled account.

### A transparent sufficient interface, not a result

Suppose the actual original-trial probabilities can be compared uniformly
to the two one-source reference probabilities with additive event error
eta_L<=1/100, and suppose one proves for those source preparations

    p_one_matter(t)>=5/9 for 2<=Jt<=L/6,
    p_one_receiver(t)>=3/5 for L/8<=Jt<=L/6.

Then p_two_matter>=25/81-1/100>1/4 and
p_R>=9/25-1/100>1/10. Exact endpoint color factorization already supplies
zero conditional two-qubit error. These sufficient source/comparison
hypotheses are **not proved** here; their constants are not changed
acceptance gates or inferred favorable data. No success is recorded.

## 4. What the route does and does not decide

The analysis establishes an inspectable way to formulate long-time bounds
without exponentiating an energy-dependent Schur matrix or accumulating a
fixed-order phase error. It isolates a mobile physical channel and a source
resolvent whose spectral measure can settle whether fuel actually reaches it.
That justifies a focused mathematical challenge to (19)-(21), not a broad
solver or architecture campaign.

The physical task remains unresolved. No claim is made that the proposed
probability thresholds hold, that they fail for all L, or that their failure
time grows linearly. Even a positive pre-encounter result would leave
collisions, repeated emission, positive density, interior-fuel decoding,
actual moving geometry, and autonomous geometry/fuel preparation open.

The fixed-record spectral gap here does not restore the previously absent
uniform real-energy storage-consistency gap of the full WM2 receiver model.
There is a further timescale warning: B1's moving-geometry comparison used
a fixed maximum physical time 15360 and a much larger recorded buffer. If
one sets the whole record length proportional to L while taking t of order
L/J, its old tail estimate need not be small at all. That is a limitation
of the estimate, not a proof of geometry destruction. A different scaling
or proof is necessary, and growth can also make the right fuel an interior
source. Neither the old buffer guarantee nor endpoint color ordering is
automatically inherited by the full law in this proposed scaling family.
Supplied locality, one-dimensional grammar, quantum structure, and fuels
remain supplied. No empirical distinction from conventional quantum
dynamics has been obtained.

Stop with this assessment. B2 remains failed at theta=6 on eight events.
A1 and the numerical recovery branch remain closed. Independent review of
the present mathematics and earlier geometry results remains pending.
