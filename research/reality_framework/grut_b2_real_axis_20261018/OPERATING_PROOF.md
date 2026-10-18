# Original preparation to the frozen operating event

Author-derived, with internal challenge and exact constant checks;
independent mathematical review pending. Dependencies: B2-BB's exact
kernel, B2-SS's invariant graph and B2-SO's finite-host comparison.
No moving geometry, new Hamiltonian or replacement endpoint.

## 1. Precise claim

Use the unchanged nominal endpoint-source law and its original two
excited fuels, matter/field vacuum, any two-qubit input with any reference.
For every L in 24N, L>=2016, let I_L=[L/(8J),L/(6J)]. The physical
receiver event is two co-pairs at x in [L/8,3L/8-1] and y in
[5L/8,7L/8-1], with both fuels empty. The derivation below gives

    p_out(t)>514089/4000000=0.12852225, t in I_L,
    p_(q=2)(t)>434281/1000000=0.434281,
                                      2/J<=t<=L/(6J).              (37)

The endpoint-law successful two-qubit/reference channel is exactly
p_out times the identity channel. Thus (37) passes the existing 10%
output, 25% retention and 1% conditional information-error gate at
AUTHOR-DERIVED status. It is not an independently validated theorem.
B2's original eight-event failed endpoint remains failed.

## 2. Return: all spectral contributions included

REAL_AXIS.md proves (or, for the minimal route, derives at regular
boundary points) the actual band density w satisfying
|w-w0|<=w0/25, where w0=(2/pi)sin^2 k and integral w0=1. All remaining
source spectral mass, of ANY spectral type, is <=1/25. Hence

    |G(t)-I_e(t)|<=2/25,
    I_e(t)=(2/pi) integral_0^pi sin^2 k exp[-it e(k)] dk.             (38)

This comparison changes the density, not the phase. It never replaces
e(k) by a nearest-neighbor dispersion in a growing-time evolution.
Use e'=2J sin k/R and integrate by parts. The endpoint ratio
sin^2 k/e'(k) vanishes. Equations (24)-(25) imply

    ||d(sin k R)/dk||_1
       <=2Rmax+5J R_E+4p R_v < pi.

Thus |I_e(t)|<=1/(J|t|) and

    |G(t)|<=2/25+1/(J|t|)<=29/50<2/3, J|t|>=2.                    (39)

No early-time extrapolation, deletion of a bound state, or eventual-
decay-only assertion enters this bound. In the physical finite system
recurrence is controlled only up to the inherited comparison horizon.

## 3. Central preparation weight is earned from the actual density

Let C=[Omega-1.5J,Omega+1.5J] and

    c(t)=P_ac 1_C(h) exp(-ith)f, w=||c||^2.

This is an orthogonal component of the ORIGINAL state, not a new
preparation. Equation (3) implies that k in [3/4,pi-3/4] lies in C.
The elementary alternating-series bounds for cos(3/4) and sin(3/2),
with pi>31/10, give

    integral_(3/4)^(pi-3/4) (2/pi)sin^2 k dk >5/6.

Therefore the actual spectral density gives w>(24/25)(5/6)=4/5.
Using P_ac explicitly keeps any unclassified component in c's
orthogonal complement. Nothing removes it from the physical state.

## 4. Continuous-window spatial capture, with endpoint terms

Set theta=Jt and K=L/8. For t in I_L, K<=theta<=4K/3 and theta>=252.
On the central interval, the group velocity in units of J is between
1.3 and 2.01 and has a single maximum. Equation (34) retains both
phases, with x=n+2. If a phase derivative has fixed sign and two
monotonic pieces, direct integration by parts, including its endpoints,
gives

    |integral b(k) exp[i phi(k)]dk|
                   <=(4||b||_infinity+TV(b))/min|phi'|.             (40)

The variation of 1/phi' is bounded separately on its two monotonic
pieces. This proves (40), rather than assuming oscillatory cancellation.
Equation (36) bounds its numerator by 8/5 for both terms in (34).

For the outgoing term, missed sites n<=K-1 have
|phi'|>=1.3theta-n-2; missed sites n>=3K have
|phi'|>=n+2-2.01theta. Summing the inverse squares and bounding the
decreasing series by integrals gives

    sum_(n<=K-1) min|phi'|^-2 <=1/(.3theta-2),
    sum_(n>=3K) min|phi'|^-2 <=1/(.24theta+1).                     (41)

For the incoming term, over ALL n>=0,
sum min|phi'|^-2<=1/(1.3theta+1). The integer offsets in these
formulas keep the physical fuel-site phase and detector endpoints.
Minkowski's inequality, not incoherent addition, now gives

    ||missed central co-pairs||
      <=(8/5)sqrt[1/(.3theta-2)+1/(.24theta+1)]
        +(8/5)/sqrt(1.3theta+1)+1/200.                             (42)

The last term is the GLOBAL source/boundary remainder in (33). It
includes near-source sites; no claim of exponential smallness there.

The fuel part of c has amplitude <=10/theta. The existing exact low
invariant graph gives separated-charge squared norm <=w/362<=1/362.
All three complements are orthogonal. The right side below decreases
with theta>=252; exact outward rational square-root bounds give

    ||(1-R_left)c(t)||^2
      <= [right side of (42)]^2+(10/theta)^2+1/362 <3/20.           (43)

This covers every time in the frozen window and every length in the
family. It uses no time grid or favorable endpoint search.

## 5. Coherent complement is paid, not discarded

For a normalized state psi and an orthogonal component c of weight w,
if a projection R misses squared norm ell from c, the sharp projection
geometry inequality from B2-SB gives

    ||R psi|| >= [sqrt(w-ell)-sqrt((1-w)ell/w)]_+.

With w>4/5 and ell<3/20, the first square root is >4/5 and the second
is <1/5. Consequently, for the ACTUAL uncut original source state,

    A_L(t)>9/25=.36, every t in I_L.                              (44)

The rest of the original state may interfere destructively. The second
term above is precisely the price for that possibility. Neither energy
components nor two sources are treated as independent probabilities.

## 6. Joint finite-host reconstruction

Use B2-SO unchanged: Pi=1_(-infinity,2.5)(h), u=exp(-ith)Pi f,
omega<=36/27613^2, sigma<=6/27613, M=L/2-2, and its same-Hilbert-space
product reference cut at M. For 0<=t<=L/(6J),

    delta_L=12/27613+1500 exp(-L/105)<1/1000,
    tau=||(1-C_M)u||<1/1000,
    a=||R_left u||^2>=A_L-2sigma,
    m=||P_matter C_Mu||^2=1-omega-|<f,u>|^2-tau^2,
    |<f,u>-G|<=omega,
    p_out >=[a-delta_L]_+^2,
    p_q2  >=[m-delta_L]_+^2.                                      (45)

The product reference's JOINT output probability is a^2, not a.
The physical evolution retains collisions, both fuel returns, all
charge separations and finite-host reflections. Those effects are
paid by the norm defect delta_L, not assumed absent.

Equations (39),(44),(45), 2sigma<1/2000 and delta_L<1/1000 imply
a-delta_L>717/2000. Also
1-omega-(29/50+omega)^2-(1/1000)^2>33/50, so m-delta_L>659/1000.
Squaring gives (37). Exact Fraction checks validate this probability
composition separately from the spectral proof.

## 7. Scope and research consequence

This is an author-derived scaling regime for two resources on a fixed,
prepared path under the unchanged endpoint law. J is fixed; times and
receiver distances grow with L. It supplies neither extensive matter
nor autonomous geometry preparation. Arbitrary distributed fuel sites
have a different information-routing problem. Locality of permitted
generators, recorded geometry, fuel energy and quantum dynamics remain
supplied. Medium gravity and buoyancy play no role in the argument.

The constructive consequence, subject to independent review, is that
the escape mechanism need not remain a one-composite device. The next
legitimate integration question is compatibility with the actual moving
geometry law; no such composition or successor run is performed here.
Independent review should first challenge (4)-(16), the derivative
series, the global remainder, and the inherited finite-host comparison.
