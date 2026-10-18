# B2-SA: all-length identities and the remaining lower bound

Author-derived mathematics; independent proof review pending. No new law or
successful operating experiment. Throughout, the path is a **fixed recorded
path**, not the complete dynamical geometry/receiver. All its matter, fuel,
electric states and reverse transitions are retained.

## 1. The family and its exact sector

There are L events, endpoint fuels, hard-core plus/minus constituents, and
capacity-one electric links. The B1/B2 constants are

    j+=1/64, j-=1/96, kappa=1/4,
    A=(j+^2+j-^2)/(5/4)=13/46080,
    J=lambda=2j+j-/(5/4)=1/3840, Omega=2-2A.

The supplied occupied-bond interaction rule, color degeneracy, degree
counterterm and constraint encoding remain inputs. The exact resource sector
has q=N+=N-=2-F_left-F_right. Its scalar Hamiltonian is H=D+V, V=T+C,
with D the actual physical diagonal, T constituent hopping, C conversion.
Gauss determines the field; reject only states violating its capacity, not
large allowed separations. The read-only B1 `PathModel` defines every entry.

Let P contain: both fuels full; one full fuel and a co-located pair at ANY
event (including the occupied fuel's event); two nonadjacent co-located
pairs and both fuels empty. Adjacent pairs and every separated charge
configuration lie in Q. For L>=4,

    M_L=rank P=1+2L+(L-1)(L-2)/2=(L^2+L+4)/2.       (1)

P is a coordinate subspace used for analysis. It is NOT the exact invariant
low band; C can couple P to Q. No physical projection is performed.

Every P diagonal lies in [4-4A,4]. Every Q diagonal is >=4+1/2-A:
two adjacent co-pairs have contact energy 2kappa=1/2 and at most one endpoint
counterterm; other two-pair Q states have at least one nonzero field link,
giving energy >=5-2A; one-pair Q states have separated constituents and are
higher still. The endpoint counterterms sum to at most 2A.

At q<=2, each constituent has at most two hops. Both row and column absolute
sums give ||T||<=4(j++j-)=tau_2=5/48. Each fuel conversion term has norm
lambda, hence ||C||<=2J. Thus

    ||V||<=v=tau_2+2J=67/640,                         (2)

uniformly in L and the number of allowed configurations. This bound pays
for the number of primitive moves; it does not assume independent histories.

## 2. An exact, uniformly convergent low-band construction

Set r=1/4 and use the counterclockwise circle |z-4|=r. Its distance to D's
spectrum is at least delta=r-4A. Put a=v/delta. The exact fractions obey

    a<3/7,    delta>v,
    beta:=r v/[delta(delta-v)]<3/4,
    r/(delta-v)<7/4.                                  (3)

For every point of this contour the FULL resolvent has the norm-convergent
series

    (z-H)^(-1)=(z-D)^(-1) sum_{n>=0}[V(z-D)^(-1)]^n.   (4)

No virtual charged/fuel word is dropped. Diagonal propagation is exact;
all multiplicities and interference are included before the norm bound.
This is possible uniformly in L because only two resources are present.
It is not an extensive-density statement or a sum over moving geometries.

Define the exact spectral projector

    Pi=(1/(2pi i)) integral_contour (z-H)^(-1) dz.      (5)

The same contour works for D+sV, 0<=s<=1, so rank Pi=rank P=M_L. Integrating
the resolvent difference gives ||Pi-P||<=beta<3/4. In particular P Pi P has
eigenvalues >=1-beta^2>7/16 on P. The canonical band isometry

    W=Pi P (P Pi P)^(-1/2),    K=W^* H W               (6)

is well-defined and obeys W^*W=P, WW^*=Pi and H W=W K EXACTLY.
Equation (6) is an analysis of the same H, not a new pair-only law.

For a normalized initial double-fuel vector f, let

    w=||(1-Pi)f||^2, chi=W^*Pi f/sqrt(1-w),
    R_tilde=W^*R_L W.

The exact normalized low-band event is

    p_low(t)=<chi, exp(itK) R_tilde exp(-itK) chi>.     (7)

Input AND observable dressing matter: replacing chi by bare f, R_tilde by
P R_L P, or K by a second-order hopping matrix needs a new error bound.
There is no |t| times truncation residual in (7). It is not a claim that
computing K or proving its long-time transport is inexpensive.

### Complete virtual words, with their distance stated honestly

Let d_conf(X,Y) be the minimum number of nonzero V transitions between
configuration sets X,Y, INCLUDING intermediate Q states. For d>=1, the
matrix block of (4) vanishes through n=d-1. Consequently

    ||1_Y Pi 1_X|| <= r/(delta-v) a^d
                     <(7/4)(3/7)^d,
    ||1_Y (H-4)Pi 1_X|| <= r^2/(delta-v) a^d
                          <(7/16)(3/7)^d.            (8)

These are coherent source-set/receiver-set operator bounds, independent of
their cardinality. They concern CONFIGURATION distance, not a newly derived
spatial metric or a lower bound on successful matter arrival. They are not
a proof of exponential locality of K in an arbitrary band eigenbasis.

## 3. Uniformly small preparation weight above the band

T has zero P-P block, so lambda_max(PHP)<=4+2J. Conversely

    QHQ >= h_* Q,
    h_*=4+1/2-A-tau_2-2J > 4+2J.                    (9)

Min-max separates exactly M_L low eigenvalues from the rest; all high
eigenvalues are >=h_*. This is the same band as (5). The prepared fuel
state has E0=2Omega and variance 2J^2 exactly (the two possible first
conversion images are orthogonal). Its spectral measure therefore gives

    w <= 2J^2/(h_*-2Omega)^2
       =32/37027225 <1/10^6.                         (10)

For any projection/event, pure-state trace distance between full evolution
and its normalized low-band part is sqrt(w)<1/1000, for ALL times and L.
Thus |p_physical(t)-p_low(t)|<1/1000. High-band population is not the missing
order-one effect. This does not control recurrence WITHIN the low band.

The conserved energy also gives an all-time resource price:

    <sum E^2> <= 2Omega+2A+tau_2+2J < 411/100.         (11)

Mass, fuel and density energy are nonnegative; H_degree>=-2A and the other
two interaction norms are covered by (2). Energy and initial resource
number are O(1), path registers O(L); density tends to zero. Prepared
geometry and loaded fuels are supplied, not autonomously produced here.

## 4. The exact infinite-bulk neutral band is mobile

This paragraph concerns the translation-invariant ONE-resource bulk,
with no endpoint fuel, no boundary counterterm, and all relative charge
separations. It identifies an asymptotic channel, not its source weight.

Use center coordinate X=(i+j)/2 and relative coordinate s=j-i. Fourier
transform simultaneous translation and gauge the phase of relative hopping.
The k fiber on ell^2(Z) is a real Jacobi operator

    h(k)=diag[2+|s|+kappa 1_(|s|=1)]
                         -b(k)(S+S^*),
    b(k)=sqrt(j+^2+j-^2+2j+j- cos k),
    1/192 <= b(k) <= b_max=5/192.                     (12)

The potential tends to infinity as |s| grows, so each fiber has discrete
spectrum. Its unique strictly positive, reflection-even ground vector u
has energy e(k)<2. Compressing away s=0 gives spectrum >=2+Delta-2b(k),
Delta=5/4. Min-max therefore isolates the ground band by at least

    Delta-2b_max=115/96.                              (13)

No relative-coordinate truncation is involved. The band is analytic because
it is simple, uniformly isolated and b(k)>0. Hellmann-Feynman gives

    e'(k)=(j+j- sin k/b(k)) <u,(S+S^*)u> >0,
                                      0<k<pi.        (14)

Here positivity of u makes the displayed expectation strictly positive.
It is useful to bound its size rather than merely calling it nonzero.
Write u0 for its s=0 component. The Q resolvent and the s=1 equation imply

    u0^2 >= a0:=1/[1+2 b_max^2/(Delta-2b_max)^2],
    u1 > b(k) u0/(Delta+2b(k)).

Since <u,(S+S^*)u> >=4u0 u1, it follows that

    e'(k) > c sin k,
    c=4j+j- a0/(Delta+2b_max) > (19/10) J.            (15)

This is an all-order physical-time velocity statement for the exact bulk
band. Negative k has the opposite velocity. A bare co-located bulk pair
has projection weight at least a0 in this band. It is NOT the double-fuel
seed of the prospective task, and a useful lower bound on source emission
does not follow from (15).

For a smooth k wave packet confined to this band, multiplication by
exp[-it e(k)] gives the usual ballistic displacement directly: in Fourier
space X=i partial_k, so

    X exp[-it e(k)]a(k)u(k)
      = exp[-it e(k)] [t e'(k)a(k)u(k)
                         +i partial_k(a(k)u(k))].    (16)

After division by t the second term has L2 norm O(1/t). This proves the
velocity-operator limit on these smooth band packets. The phase gauge's
bounded, smooth derivative on the localized band changes only the O(1/t)
term. It is not an arrival lower bound for boundary fuel preparation or
two colliding pairs on a finite path.

## 5. Endpoint information: exact, but conditional on matter output

For every fixed L, place fuel color carriers just outside the two endpoints.
Hard-core plus carriers cannot pass on the path, and endpoint conversion
moves the extreme carrier between the endpoint and its adjacent fuel.
Emission/absorption is blocked whenever it could cross an occupied carrier.
Minus carriers carry no color. Thus ordered colors never exchange, and

    H_colored = H_scalar tensor I_4.                  (17)

For any color-blind outgoing event and left/right matter decoding, the
successful instrument is exactly p_R times the identity on two qubits,
including entanglement with an untouched reference. Failure is an erasure.
Color error on success is zero; success probability remains unproved.
No triplet recoding, fuel readout or global recovery is substituted.

This statement does NOT survive automatically when growth turns an endpoint
fuel into an interior fuel. B2's exchange word and B2-DX's channel analysis
remain valid warnings. A full moving-geometry theorem needs its own buffer,
geometry comparison and accessible output-channel argument.

## 6. Why the old approximation is not a scaling proof

B1 provides a valid short-time bound e0+rho t with a size-independent
rho>0. At t=L/(6J), its certified allowance is e0+[rho/(6J)]L. It eventually
exceeds any fixed useful error budget. This is failure of that bound to
settle scaling, NOT evidence the physical error or transport grows that way.
Even the one-resource estimate is insufficient; two-resource collisions
prevent tensorizing it in the first place. Using another fixed-order
effective Hamiltonian has the same problem unless its omitted phase/action
is controlled through O(L/J) times.

Equations (5)-(7) provide a controlled way around that particular analytical
obstacle, by retaining the band exactly. They leave the decisive tasks:

* Determine the boundary fuel's spectral measure, especially emitted vs
  localized weight, with a quantitative lower bound on relevant velocities.
* Control two-resource scattering, endpoint returns, and finite-host effects
  over the ENTIRE prospective window; do not assume two independent emitters.
* Prove the dressed event in (7) has its prescribed lower bound, and the
  q=2 retention probability stays above 1/4 from preparation to the window.
* Only then compose with the actual moving geometry and its receiver.

Neither a propagation upper bound, a mobile bulk band, nor a gap proves
these lower bounds. No finite-size recurrence estimate has been extrapolated.
