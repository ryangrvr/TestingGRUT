# Exact source decomposition and quantitative closing conditions

Author-derived, independent mathematical review pending. Conditional
lemmas below are not established hypotheses for the physical source.
The B2-SO joint comparison is inherited unchanged, at its existing status.

## 1. One experiment, three different probability levels

Use the exact half-line source h and original fuel f of B2-SS. On
f and |i,j>, i,j>=0, its entries are

    h_ff=Omega=2-2A, h_f,00=-J,
    h_ij,ij=2+|i-j|+(1/4)1_(|i-j|=1)-A 1_(i=j=0),
    plus hops=-1/64, minus hops=-1/96,
    A=13/46080, J=1/3840.

All charged separations and both conversion directions remain. Set
v(t)=exp(-ith)f, theta=Jt, and let R=R_left,L project onto the physical
co-pairs |n,n> with L/8<=n<=3L/8-1. Define

    A_L=||R v||^2,       G=<f,v>.

A_L is ONE-source arrival. G is a scalar spatial amplitude, and the whole
one-resource fuel map is G I_2, including reference-entangled inputs.
Thus |G|^2 is its full fuel occupation, without output conditioning.

Let Pi=1_(-infinity,5/2)(h), u=Pi v, and omega=||(1-Pi)f||^2. Use

    sigma=sqrt(omega)<=6/27613, omega<10^-6,
    a_L=||R u||^2,       g=<f,u>,
    M=L/2-2, C=C_M,      tau_L=||(1-C)u||,
    m_L=||P_matter C u||^2.

u is NOT normalized. B2-SO gives, for L in 24N, L>=2016 and
0<=t<=L/(6J), a common-space reference

    Phi_L=E_L(Cu tensor Cu),
    ||Psi_L-Phi_L||<=delta_L,
    delta_L=12/27613+1500 exp(-L/105)<1/1000,
    tau_L <= (101/100)exp(3/10-L/105).                  (1)

E_L embeds reflected, disjoint source packets. It is an exact isometry,
not an assertion that the physical evolution factorizes. For the actual
joint outgoing projector and q=2 projector,

    ||R_joint Phi_L||=a_L,       ||P_q2 Phi_L||=m_L.    (2)

The reference probabilities are a_L^2 and m_L^2, respectively. If
r=Psi_L-Phi_L, the exact joint event decomposition is

    p_out=a_L^2+2 Re<R_joint Phi_L,R_joint r>
                    +||R_joint r||^2.                (3)

Hence, without independence or removal of collisions/returns,

    p_out >= [a_L-delta_L]_+^2,
    p_q2  >= [m_L-delta_L]_+^2.                      (4)

Here [x]_+=max(x,0). Since R is inside C,

    a_L >= [sqrt(A_L)-sigma]_+^2 >= A_L-2sigma,
    |g-G|<=omega,
    m_L=1-omega-|g|^2-tau_L^2.                       (5)

The second bound is omega, not sigma: the discarded fuel amplitude is
the Fourier transform of a positive measure of mass omega. Consequently

    p_out >= [[sqrt(A_L)-sigma]_+^2-delta_L]_+^2,
    p_q2 >= [1-omega-(|G|+omega)^2
                         -tau_L^2-delta_L]_+^2.      (6)

No errors in (5) replace the separate physical reconstruction error (1).
They are different legs: full source to filtered source, then product
reference to original finite two-source evolution. None is renormalized.

The unchanged sufficient conditions A_L>=17/50 on I_L and |G|<=2/3 on
[2/J,L/(6J)] imply

    p_out > (17/50-1/2000-1/1000)^2
           =458329/4000000=0.11458225,
    p_q2>1/4.                                      (7)

The exact endpoint color factor I_4 applies to this SAME outgoing event;
the successful channel is p_out rho, including an external reference.
Quality conditional on success does not supply the missing probability.

## 2. Exact coherent energy and spectral-type decompositions

Write C_en=[Omega-3J/2,Omega+3J/2], wholly within Pi, and

    c=exp(-ith)1_C_en(h)f,
    b=exp(-ith)(Pi-1_C_en(h))f,
    z=exp(-ith)(1-Pi)f.

These vectors are orthogonal BEFORE detection. Their norms are constant:

    ||c||^2=w_C>=5/9,
    ||b||^2=1-omega-w_C,     ||z||^2=omega.            (8)

But the exact detector probability is

    A_L=||Rc||^2+||Rb||^2+||Rz||^2
             +2 Re(<Rc,Rb>+<Rc,Rz>+<Rb,Rz>).        (9)

The mixed terms need not vanish. The high state z costs at most sigma
in detector AMPLITUDE; its fuel contribution costs omega. The remainder
b can carry O(1) weight. Neither its receiver contribution nor its phase
relative to c has yet been controlled. It is incorrect to sum the three
displayed diagonal probabilities as though detection dephased energy.

For the positive fuel measure nu(B)=||1_B(h)f||^2, define nu_C, nu_b,
nu_z by restriction to the three energy sets. Then exactly

    G=G_C+G_b+G_z, G_S=integral exp(-itE)dnu_S(E),
    |G_z|<=omega.                                  (10)

Within each energy set one may further split nu into pure point,
singular-continuous and absolutely continuous parts. Equivalently, with
P_pp,P_sc,P_ac the commuting spectral-type projections,

    A_L=sum_(S,T in {C,b,z}; a,b in {pp,sc,ac})
            <R exp(-ith)P_a 1_S(h)f,
             R exp(-ith)P_b 1_T(h)f>.               (11)

Equations (9)-(11) include all cross terms. Finite matrices do not classify
these infinite-source spectral types. The no-eigenstate theorem at/below
the bottom leaves possible embedded or higher localized contributions.

In a spectral partition limit, the position-sensitive kernel is

    K_L(B,B')=<1_B(h)f,R 1_B'(h)f>,
    A_L=double integral exp[it(E-E')]dK_L(E,E').      (12)

The scalar measure nu does not give K_L without the physical detector's
representation. In the cyclic spectral representation, nu determines h
and f up to unitary equivalence; it does not determine the extra operator
R. Here R IS fixed physically, so its intertwining/eigenfunction data must
be computed, rather than chosen after inspecting propagation.

## 3. A single central-packet estimate sufficient for arrival

The following lemma uses spectral orthogonality and pays for interference
without classifying every noncentral mode.

Lemma. Let psi be normalized, c an orthogonal component of psi with
||c||^2=w>0, and R an orthogonal projector. Set eta=||Rc||^2/w. Then

    ||R psi||^2 >=
      [sqrt(w eta)-sqrt((1-w)(1-eta))]_+^2.          (13)

Proof. Write psi=sqrt(w)chi+sqrt(1-w)zeta with chi,zeta orthonormal.
For eta>0 let e=R chi/sqrt(eta). It is a unit vector in range R and
|<e,zeta>|<=||(1-|chi><chi|)e||=sqrt(1-eta). Thus
||R psi||>=|<e,psi>|>=sqrt(w eta)-sqrt((1-w)(1-eta)).
For eta=0 the asserted zero lower bound is trivial. The two-dimensional
choice chi=(sqrt(eta),sqrt(1-eta)), zeta=(-sqrt(1-eta),sqrt(eta)) saturates
the nonzero branch. Global orthogonality has NOT been mistaken for
orthogonality after projection.

Apply (13) to the ACTUAL c in (8). Define the unresolved physical quantity

    eta_C(L,t)=||R_left,L exp(-ith)1_C_en(h)f||^2/w_C. (14)

A sufficient missing estimate is

    inf_(L in 24N,L>=2016; t in I_L) eta_C(L,t)
                      >=24/25.                    (15)

Because the right side of (13) is increasing in w and eta in this range,
(8) and (15) imply

    A_L >= [sqrt(8/15)-2/15]^2
          > (179/300)^2=32041/90000 >17/50.          (16)

The strict inequality uses sqrt(8/15)>73/100, verified by squaring
positive rationals. Thus (15) would close the ORIGINAL one-source arrival
condition, including worst-case coherent interference of every other
spectral component. It is a sufficient reduction, not a claim that 96%
capture is necessary or already proved.

Why bulk velocities make (15) plausible but do not establish it: B2-SS
places the central energies at outward bulk speeds 1.1J<v<2.2J. Their rays
over I_L lie strictly within [L/8,3L/8], with minimum continuous-coordinate
margin L/120; using integer detector sites costs at most one site.
But energy compatibility gives neither outward scattering weight nor
finite-time capture. Embedded/singular states, physical co-pair residue,
boundary scattering, slow/transient components and the finite-time
oscillatory remainder all enter (14). They have not been bounded by 4%.
No outgoing wave operator or asymptotic completeness is silently assumed.

## 4. A late-return estimate sufficient for retention

Push the LOW fuel measure forward by x=(E-Omega)/J. Write exactly

    dmu_low=sum_j alpha_j delta_xj+dmu_sc+r(x)dx,
    beta=sum_j alpha_j+mu_sc(R),
    g(t)=exp(-iOmega t)[sum_j alpha_j exp(-i theta xj)
             +integral exp(-i theta x)dmu_sc
             +integral exp(-i theta x)r(x)dx].        (17)

All symbols describe the actual measure, not independently adjustable
mixture weights. If D(theta) bounds the last Fourier integral, then

    |G(t)|<=omega+beta+D(theta).                    (18)

One explicit, separately UNPROVED sufficient regularity condition is
that r, extended by zero outside its support, is of bounded variation,
with V=|Dr|(R). Integration by parts as a distribution gives

    |integral exp(-i theta x)r(x)dx|<=V/|theta|.      (19)

All edge jumps are included in V. This is a dimensionless density; using
the variation of the dimensional E-density without the J rescaling
would give the wrong constants. The L1 and BV hypotheses justify the
boundary limit. Singular measures are retained in beta, not differentiated.

For example, the missing source conditions

    r in BV(R),     omega+beta+V/16<=2/3             (20)

would give |G|<=2/3 for ALL theta>=16. Combined with B2-SS's continuous
theta in [2,16] result, they would cover the complete retention interval
for every L, without a growing-time Duhamel error. This is only one
sufficient route; a suitable direct D(theta) can replace BV. A failure
of (20) would not refute retention. In particular threshold singularities
may require another oscillatory estimate instead of BV.

No current bound controls beta, V, or D(theta) strongly enough. First-band
mass >79.9% and central mass >=5/9 are ENERGY weights, not bounds on
the singular mass beta. They cannot be inserted into (18) as decay rates.

## 5. Exact error and obligation ledger

| Term | Current control | Required missing control | Dependence / contribution |
|---|---|---|---|
| Physical two-source reconstruction r | ||r||<=delta_L<.001 | None added here; inherited proof review pending | delta_L=12/27613+1500e^(-L/105), uniform through L/(6J); subtracts in joint event amplitude |
| High source state z | norm sigma<=6/27613; measure mass omega<=sigma^2 | None added here | Constant in L,t; costs sigma in source detector amplitude, omega in fuel amplitude |
| Source outside cut | tau_L<1.01e^(.3-L/105) | None added here | Uniform whole use horizon; subtracts tau_L^2 from m_L |
| Central receiver component | w_C>=5/9; compatible bulk energies | eta_C>=24/25 on every t in I_L is sufficient | Receivers and time both grow with L; (13)-(16) imply A_L>.34 |
| Noncentral detector amplitudes and mixed terms | Global orthogonality; no small individual detector bound | (13) handles worst case IF eta_C is bounded | O(1), coherent; cannot discard or add incoherently |
| Low singular fuel measure | No full spectral-type classification | beta plus continuous Fourier remainder small enough | May be O(1), no guaranteed decay; enters (18) |
| Continuous low fuel measure | Existence of its Lebesgue part only | D(theta), e.g. BV plus (20) | Must hold theta>=16 through unbounded L/6, not merely to 16 |
| Successful two-qubit channel | Exact I_4 for fixed endpoint law | Probability of SAME output event | Conditional quality is exact; does not increase acceptance |

## 6. Verdict of this decomposition

The gap has been converted to two quantitative sufficient estimates:
central-packet capture (15), and late fuel Fourier control such as (20).
Neither is proved here for the actual h. The separately proved logical
witness in SCALAR_LIMITATION.md shows why the banked scalar summaries
cannot replace them. This is an analytical-method limitation with an exact
scope, not an actual source violation or a uniform operating bound.
