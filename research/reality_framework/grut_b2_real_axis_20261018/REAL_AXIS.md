# Actual-source real-axis construction

Author-derived; independent proof review pending. All resolvents below
use **(z-h)^-1**. Inputs are pinned in SOURCE_PINS.json. No new law.

## 1. Operator, band and notation

The source Hilbert space is C f + ell2(N_0^2). The diagonal on |i,j> is
2+|i-j|+(1/4)1_(|i-j|=1)-A1_(i=j=0); constituent hops are -jp and -jm.
Fuel energy is Omega and f couples only to |0,0> with -J. Set

    jp=1/64, jm=1/96, s=jp^2+jm^2, p=jp jm,
    a0=5/4, A=s/a0=13/46080,
    J=2p/a0=1/3840, Omega=2-2A,
    v(k)=s+2p cos k, v*=v(0)=25/36864,
    u=4v*/a0^2=1/576, a_max=a0+5J.

Let e(k), 0<=k<=pi, be the already constructed first bulk band. It is
strictly increasing in the interior and lies in I_b=[2-5J,2]. The actual
closed-arm self-energy is sigma(E,k)=-f(E,v(k)), where f is defined by
the B2-BB continued fraction. Its positive closed-word expansion is

    f(E,v)=2v/a1 + sum_(r>=1)c_r(E)v^(r+1),
    a1=13/4-E, 0<=c_r(E)<=2*4^r/a1^(2r+1).       (1)

Each word has 2r hops and 2r+1 positive denominators. Dyck paths are
at most 4^r; differentiating denominators also bounds every positive
mixed E,v derivative by that of

    F(E,v)=2v/a1 / (1-4v/a1^2).                   (2)

This is a coefficient majorant, not an effective Hamiltonian. In
particular f_v>=2/a1 and f_E>=0. All constants below are evaluated at
the worst point a1=a0,v=v* unless explicitly stated otherwise.

The band equation is E-2+f(E,v(k))=0. A sharper useful bound is

    |e(k)-(Omega-2J cos k)| < J/250.              (3)

To prove it, write d=2-E>=0. The difference in (3) is exactly
2v d/[a0(a0+d)] - [f(E,v)-2v/(a0+d)]. The two bracketed contributions
are nonnegative. The first is <=(2v*)^2/[a0^3(1-u)]=J/276.
For the second use the symmetrized arm adjacency C: ||C||^2<=u and
(C^(2n))11=C12^2(C^(2n-2))22 <=v/(a1 a2) u^(n-1), a2=4-E>=2.
It is <=2v*^2/[a0^2*2*(1-u)]=(5/4416)J.
Take their maximum, not their sum. Both are <J/250.

## 2. Extract the actual band zero before inversion

Fix E=e(k0), v0=v(k0), x0=cos k0. Define the divided difference

    a_E(w)=[f(E,v(w))-f(E,v0)] / [(2/a1)(v(w)-v0)],
    v(w)=s+p(w+w^-1), J_E=2p/a1.                 (4)

Use its removable value at v(w)=v0. Then exactly

    kappa(E,w)=E-2-sigma(E,w)
              =J_E(w+w^-1-2x0) a_E(w).          (5)

This extracts a zero of the ACTUAL infinite-range kernel; no scalar
nearest-neighbor model replaces that kernel. The divided-power formula
and (1) give nonnegative, even Fourier coefficients for a_E-1 and

    ||a_E-1||_W <= eps=(1-u)^(-2)-1,
    1<=a_E(exp(ik))<=1+eps.                       (6)

Here ||.||_W is the sum of absolute Fourier coefficients. Since eps<1,
the convergent log series is in W, real and even. Let a_+,a_- be its
one-sided exponential factors with half of the constant in each, and
L=T_+(1/a_+). Write its real coefficients l_j, j>=0, and
c0=l_0=exp[-(log a_E)_0/2]. Then

    ||L||_W <=(1-eps)^(-1/2)<501/500=:l,
    (1+eps)^(-1/2)<=c0<=1,
    960J/961<=J_E<=J.                            (7)

Triangular coefficient convolution, with entries a_(n-m), proves

    T_+(kappa)=T_+(a_-) T_+(kappa0) T_+(a_+),
    g_T=L g0 L^T,
    (g0)_nm=-exp[ik0(max(n,m)+1)]
                  sin[(min(n,m)+1)k0]/(J_E sin k0). (8)

Transpose, not adjoint, is correct. At both thresholds take the
continuous removable limit of the finite geometric sum. Consequently

    |(g_T)_nm|<=l^2(min(n,m)+1)/J_E.              (9)

This estimate is uniform INCLUDING k0=0,pi. It is not an unweighted
ell2 inverse on the propagating band.

## 3. The actual boundary is inverted in localized spaces

Keep the exact U(E)=W(E)-A|0><0| of B2-BB. Use entrywise moments,
not its trace norm alone:

    M0 >= sum_nm |U_nm|,
    M1 >= sum_nm (m+1)|U_nm|,
    M11>= sum_nm (n+1)(m+1)|U_nm|,
    Mr >= sum_(m>=1,n>=0) m |U_mn|.               (10)

The actual real-energy kernel is symmetric, U_nm=U_mn. Thus the same
M1 controls either the row or column weight; the X estimate uses the
column weight on its input. SCATTERING_BOUNDS.md gives explicit
complete word-sum majorants. Put
C=l^2, j0=960J/961, eta=C M1/j0. Exact arithmetic gives eta<1/300.
On BOTH ell-infinity and X={x:sup |x_n|/(n+1)<infinity}, (9) gives

    ||g_T U|| <=eta<1,
    S=(I-g_T U)^-1 g_T.                          (11)

For X, use min(n,m)+1<=n+1; for ell-infinity use <=m+1.
All virtual boundary words are retained. The Neumann series converges
uniformly in these weighted/localized senses through the thresholds.
In particular

    (g_T)00=-c0^2 exp(ik0)/J_E,
    |S00-(g_T)00|<=C^2 M0/[j0^2(1-eta)].        (12)

The second bound uses the row-0 estimate and ||S .0||_infinity
<=C/[j0(1-eta)], rather than multiplying trace and ell2 norms.

## 4. These are genuine upper-half-plane boundary values

The implicit equation z-2+f(z,v0(z))=0 has an analytic solution near
each point of the closed band, since f_v>0 there. Near this compact
band, Im v0(z)<0 for Im z>0 sufficiently small. The solution r(z) of
r+r^-1=(v0(z)-s)/p with |r|<1 therefore has boundary r=e^(ik0),
including the square-root thresholds. The divided difference (4) is
analytic in the same neighborhood, stays close to 1 in W, and the
triangular factorization gives the actual off-real Schur inverse.

The geometric-sum bound (9) holds for |r|<=1, with slightly relaxed
constants nearby. Localized U(z) has uniformly summable moments there.
The product Q=g_T U maps X uniformly into ell-infinity (use M11);
truncate its output and the exponentially localized input sums to
deduce operator-norm continuity on X, even though it need not be
continuous on ell-infinity. Thus (11) and each source matrix element
converge to the displayed limits. No interchange of an uncontrolled
eta->0 limit with an unweighted ell2 inverse is used.

The fuel denominator and resolvents are

    D_f=E-Omega-J^2 S00,
    F_f=1/D_f, F_n=-J F_f S_n0.                  (13)

Combining (3), (12), c0^2 in [1/(1+eps),1], and (7) gives

    D_f=-J exp(-ik0)+J d_f(k0),
    |d_f|<1/100, hence |D_f|>99J/100.            (14)

See the exact constants in SCATTERING_BOUNDS.md. The strict inequality
persists in a thin upper neighborhood of the compact band. F_f is
therefore uniformly bounded there. Stieltjes inversion implies that
the fuel spectral measure has NO point or singular continuous part
on this band. The endpoints have zero source mass too. This is a
source-visible conclusion, not classification of every transverse
state or equality with the finite-host spectrum.

## 5. Actual source density, with unclassified weight paid

Set u_n=sum_(q=0)^n l_(n-q) sin[(q+1)k0],
psi=(I-g_T U)^-1 u. The free jump is rank one. Resolvent jump identities
give exactly

    -Im S00=|psi_0|^2/(J_E sin k0),
    rho_f(E)=J^2 |psi_0|^2/
                         [pi J_E sin k0 |D_f|^2]. (15)

The removable endpoint limits are intended. Put psi_0=c0 sin k0(1+zeta).
The X version of (11) bounds |zeta| by eta*l/[c_min(1-eta)],
c_min=1/(1+eps). Differentiation and the remaining constants are
specified, not assumed, in SCATTERING_BOUNDS.md.

The pushforward density w(k)=rho_f(e(k))e'(k) satisfies

    (24/25) w0(k) <= w(k) <= (26/25) w0(k),
    w0(k)=(2/pi)sin^2 k, 0<=k<=pi.               (16)

This is an actual spectral-density estimate, not a moment-matched
surrogate. Since integral w0=1, all source weight outside this band is
at most 1/25. Every unclassified atom, singular part and higher band is
included in that remainder; none is dropped from physical evolution.
The sharper inherited high-energy tail remains 36/27613^2.

The coherent central component used later may conservatively be defined
as P_ac 1_[Omega-1.5J,Omega+1.5J](h) exp(-ith)f. The band result above
makes P_ac redundant here, but avoids relying on any broader type claim.

## References and inheritance boundary

Wiener-Hopf factorization is standard methodology (M. C. Camara,
arXiv:1710.11572). The actual coefficient bounds, real-axis branch,
weighted inverse and denominator estimates are derived above; no
general scalar finite-defect Jacobi scattering theorem is imported.
The closed-channel word expansion and finite-host comparison remain
author-derived dependencies pending independent review.
