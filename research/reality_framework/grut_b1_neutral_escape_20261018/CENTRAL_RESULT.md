# B1: matched preparation, lifetime and information transport

Author-derived within the declared candidate; independent mathematical review
pending. The estimates below do not depend on S9. Its formation theorem is
needed only for the final original-seed acceptance corollary.

## 1. Exact one-resource reduction

On a recorded path of m events, the complete one-resource space has fuel
|f> and every charged position |i,j>, 0<=i,j<m. Color factors as I_2 and is
temporarily suppressed. No separation, orientation, or reverse process is
discarded. Gauss determines the unique capacity-one electric string. The
matrix is exactly

    H|i,j>=[2+|i-j|+(1/4)1_(|i-j|=1)
                    -A(2-deg(i))1_(i=j)]|i,j>
             -j_+ sum_(i'~i)|i',j>-j_- sum_(j'~j)|i,j'>;
    <f|H|f>=Omega;  <0,0|H|f>=-lambda.

All other fuel-matter elements vanish. These are the primitive matrix
elements of `GRUT_CANDIDATE.md`, not a postulated pair-only evolution.

Let P span fuel and co-located pairs, Q=1-P. Subtract 2I and write

    H-2I = [ A_P  B^* ; B  D+T ].

D is the complete diagonal |i-j|+(1/4)1_(|i-j|=1) on Q, so D>=Delta=5/4.
T is the Q compression of constituent hopping. Both ||B||,||T||<=tau=5/96,
uniformly in m. The diagonal degree correction acts only on P. At the
nominal parameters Omega0=2-2A and lambda=J,

    K:=A_P-B^*D^(-1)B=-2A I-J A_(m+1),                    (1)

where A_(m+1) is unit adjacency on the path fuel,0,1,...,m-1.
The two virtual charge-hop orders give J=2j_+j_-/Delta=1/3840.
The explicitly supplied degree correction cancels the endpoint self-energy
defect. It is necessary for this particularly simple reference channel.

Thus the *matched* effective description contains fuel energy, conversion
fuel<->pair, moving composites, their binding virtual states, and density
interaction. No conversion term was set to zero. In the one-resource sector
there is no pair-pair collision term; the original density term still acts
in every virtual charged state. At two resources it generates nontrivial
collisions, and the present effective reduction does not extend by addition.

## 2. A time-domain construction, not exponentiation of a Schur matrix

Define k=2A+2J=||K|| upper bound, r=tau/Delta, and

    x1=-D^(-1)B,
    x2= D^(-1)T D^(-1)B,
    x3=-D^(-1)T D^(-1)T D^(-1)B,
    xK=-D^(-2)B K,
    X=x1+x2+x3+xK,   Z=[I;X].                            (2)

Z is an approximate intertwiner, NOT an assumed unitary or exact spectral
isometry. Hopping flips parity of i+j, and co-located states have even
parity. Consequently B^*x2=0 exactly. Direct block multiplication, using
(1), gives residual R=(H-2I)Z-ZK with

    R_P=B^*(x3+xK),
    R_Q=T(x3+xK)-x2 K-(x3+xK)K.                          (3)

This identity is the load-bearing check. It keeps the whole Q space and
never assumes its dynamics is irreversible. Set

    u=r^3+tau*k/Delta^2,
    p=tau*u,
    q=p+r^2*k+u*k,
    rho=sqrt(p^2+q^2).

Then ||R||<=rho. Since B|f>=0 and K|f>=-2A|f>-J|0,0>,

    ||Z|f>-|f>|| <= e0=J*tau/Delta^2=1/115200.

Exact unitary Duhamel proves, for the physical one-resource evolution,

    || exp(-itH)|f>-Z exp[-it(2I+K)]|f> || <= e0+rho|t|.   (4)

This is a bound against a generally unnormalized dressed vector. For every
projection R contained in P, RZ=R, so (4) controls its physical amplitude
against the ordinary normalized effective walk without adding ||X||.
For matter *survival* use the complementary fuel amplitude; separated
matter in Q counts as matter, not failure. This is why no spectral-weight
or eigenvalue estimate has been substituted for an operating-time bound.

The rational constants obey rho<1/105000. At t<=4/J=15360,

    e0+rho*t <147/1000.                                  (5)

Both squared inequalities can be checked with exact fractions. No outward
floating-point matrix certificate is claimed.

## 3. End effects and nominal walk

Put theta=Jt. The half-line unit-adjacency walk from its endpoint has

    phi_n(theta)=i^n (n+1) J_(n+1)(2theta)/theta,
    n=0 fuel; n>=1 is pair position x=n-1.                 (6)

This follows by the sine transform or the image method applied to the
full-line walk. The scalar -2A in K supplies only a global phase.
In particular |phi_0|<=1/theta for theta>0, since |J_1|<=1.
For a finite m-event record, cutting the far edge changes the endpoint
state by at most

    e_end <= (1/2) tail_(m+1)(2theta),
    tail_d(z)=sum_(ell>=d) z^ell/ell!.                     (7)

Indeed the finite evolution reaches its far boundary only through at
least m adjacency factors; integrate that boundary drive in Duhamel.
For m>=64 and theta<=4, e_end<10^(-18), using
tail_d(z)<=exp(e*z-d). This is a comparison of full evolutions, not an
absorbing boundary in the physical candidate.

At theta=4, the unarrived half-line probabilities are

    1-p_(x>=4)=sum_(n=0)^4 [(n+1)J_(n+1)(8)/4]^2 <1/8.    (8)

The inequality is checked by an exact rational enclosure of each Bessel
power series: sum through k=80, then bound the decreasing alternating tail
by its first omitted term. The identity (6) supplies the normalization;
no normalization is fitted to the computed success event.

## 4. Declared parameter interval and bounds

Relative to nominal H, delta_Omega |f><f| and the fuel-pair coupling change
have combined operator norm at most 2J/1000 in the one-resource sector.
Their complete finite-time contribution is <=2theta/1000, hence <=.008
through theta=4. The same physical Hamiltonian is used before and after
production; this comparison is mathematical only.

Combining (4)-(8), e<31/200=.155, including far-end error. Thus for every
theta in [2,4],

    p_matter >=1-(1/2+31/200)^2
              =22839/40000=.570975 >1/2.                 (9)

At theta=4,

    p_arrival(x>=4) >=(sqrt(7/8)-31/200)^2 >3/5.           (10)

The strict margins leave room for the exponentially small full-geometry
comparison below. Statement (9) is continuous-time uniformity on an
interval, not a grid interpolation. Statement (10) is arrival at a receiver
set, not arrival at one chosen event or a relativistic propagation speed.

## 5. Complete geometry, storage and receiver

This extension requires the actual accepted initial record m>=N/2, not an
invented low-energy geometry packet. Retain ell=N/4, let b=N/16-1, and
L=ell-2-b=3N/16-1. Reference dynamics reflects the tail at ell and cuts
matter crossing the retained boundary; it is embedded on the SAME full
Hilbert space, with all actual receiver and storage states. Boundary degree
counterterm differences are included among the cut drives, not ignored.

Assign q_x=n_+,x+n_-,x+2n_f,x<=4. Conversion conserves q_x on its event;
each hop conserves q_x+q_(x+1). Diagonal terms commute with the weighted
resource W=sum exp(x-b)q_x. Initially W<=C=4/(1-exp(-1)). Keeping the
electric registers as local variables, bond and contact norms are bounded
by the same conservative local sum J0=6: kinetic <=2(j_++j_-), electric
<=1, density <=1/2, degree terms <=2A per incidence, and conversion <.001.
On-event mass/fuel terms can be treated in the local interaction picture.
The one plus's color increases dimension but not these norms.

The weighted quadratic-form/Gronwall argument gives

    <W(t)><=C exp(v|t|), v=12(e-1).

The only uncanceled reference drives are tail penetration by assisted
append/pop (norm <=1.1), resource-dependent end motion, and the occupied
boundary terms. The additional degree counterterm is <.001 per boundary;
the old conservative drive coefficient 14 has slack for it (it bounds
twice the local bound 6 plus the assisted/bare boundary terms). Exact
Duhamel, keeping every full-H excursion and return on the physical side,
therefore gives

    D_B(t)=min(2,tail_(m-ell+1)(1.1|t|)
             +(14 sqrt(C)/v)exp(-L/2)[exp(v|t|/2)-1]).      (11)

Empty-head and zero-field gates ensure reflection changes cannot orphan
fuel/charge or expose an uncounted neutral particle. In the reflected
reference the retained word, storage and link support are fixed; tail and
actual conditional receiver continue evolving and factor from the retained
one-resource Hamiltonian. Its far boundary degree convention differs from
the full law only by a resource-supported boundary drive already counted.

At N>=2^28 and t<=15360, (11) is <exp(-N/5)+exp(-N/12)<10^(-18).
For example v<24 and 1.1e t<46000. The two logarithmic exponents are at
most 46000-N/4 and 184323-3N/32, respectively, which imply the claimed
bounds at this N. This is a fresh evaluation on the B1 time window, not
inheritance of A1's shorter declared use window.

Add D_B to the amplitude error in (9),(10). It still totals <31/200.
Arrival now includes intact recorded prefix and all conditional receiver
fields. Prefix/storage error probability is <=D_B^2. A matter source does
not receive an unrecorded reset or access to a chosen receiver eigenmode.

## 6. Information, times and energy in this same experiment

H_B1 factorizes as H_spatial tensor I_2 on the one-resource space: color is
carried by excited fuel or the one plus. Every geometry/field gate is
color-blind. The successful arrival instrument therefore gives EXACTLY

    E(rho)=p_arr rho direct-sum (1-p_arr)|erasure><erasure|.

This proves arbitrary-qubit readability on success and pays its probability
without postselection normalization. The instrument supplies a color
decoder; equal color energies/couplings are assumptions, not a derivation
of matter universality or a protected qubit at finite density.

Matched times: by tau_prepare=2/J=7680 the probability of matter exceeds
1/2; the operation ends at tau_transport=4/J=15360 with distance>=4 and
success >3/5. No population threshold crossing occurs in that interval.
Thus the declared useful-lifetime window lasts at least 2/J=7680 after
preparation; a first subsequent loss time is >=4/J if it exists. We do not
identify an asymptotic decay constant or assert that loss never happens.
Finite-path recurrence ultimately precludes literal irreversible discharge.

In the retained reference, added energy Omega is conserved. Since matter
mass, density and fuel energies are nonnegative, ||H_hop||<=tau,
||H_conversion||<=lambda and H_degree>=-2A,

    <H_E(t)> <=Omega+tau+lambda+2A <2.06.                  (12)

Full-H field expectation adds at most 2(N-1)D_B. No extensive root charge,
silent bath or discarded recombination is used. The price is a supplied
extended path and near-band fuel; the released resource becomes a moving
neutral excitation rather than an outgoing unobserved sink.

## 7. What composes, what does not

Matter/fuel vacuum equals the old geometry-formation block exactly. If S9
is valid and recording/fuel loading are supplied, p_form>1/2 combined with
(10) gives original-trial accepted-and-arrived probability >3/10. This is
dependent on S9, not a new validation or reduced receiver resource demand.

Two resources access pair-pair density, hard-core exchange, additional
electric strings and fuel reabsorption. A second-order co-pair description
already contains +2kappa d_x d_(x+1), an interaction much larger than J.
The single-resource norm/gap argument cannot be tensorized over cells.
The prescribed two-cell full-state controls are diagnostics only. There is
no analytic many-cell preparation/lifetime or collision-information theorem
in this packet. That is a boundary of the result, not another solver task.
