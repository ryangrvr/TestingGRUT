# Exact bulk symbol plus exponentially localized boundary kernel

Author-derived; independent mathematical review pending. This analyzes
the actual source, all charge separations and all reverse paths.

## 1. Translation-invariant bulk is derived from the whole-line operator

Extend the SAME constituent hopping and relative potential to (i,j) in Z^2,
without fuel or root counterterm. Let P_infty span |n,n>, n in Z, and
Q_infty its separated-charge complement. Simultaneous translation of both
coordinates commutes with the bulk operator and its Q resolvent. Hence

    (Sigma_infty(z))_(n,m)=sigma_(n-m)(z),
    sigma(z,k)=sum_(d in Z) sigma_d(z) exp(ikd).       (1)

This is a derived whole-line self-energy, not an assumed lattice law.
After Fourier transforming n=min(i,j), a phase gauge on each relative
arm leaves e_(r=0) unchanged and makes the arm hopping magnitude

    v(k)=jp^2+jm^2+2jp jm cos(k), b(k)=sqrt(v(k)).

For real k, v(k)>0. Each separated arm has diagonal q_1=13/4,
q_r=2+r (r>=2) and hopping -b(k). Its boundary resolvent gives exactly

    sigma(z,k)= 2 v(k) /
       [z-13/4 - v(k)/[z-4 - v(k)/[z-5 - ...]]].     (2)

The continued fraction is the actual arm-resolvent limit, with its branch
selected by the closed operator. The closed gap and convergent path
expansion justify it for Re z<=5/2, including real energies. It is not
an arbitrary scalar model sharing the source's moments.

The bulk low-energy equation is

    kappa(E,k)=E-2-sigma(E,k)=0.                    (3)

The relative fiber's open P space has rank one; its Q fiber is >=307/96.
Below 5/2 the monotonic scalar Schur function has at most one zero, and
the already constructed ground fiber band supplies that zero e(k).
Consequently this region has ONE bulk propagating branch, with infinitely
many eliminated closed transverse states. Its derivative is exactly

    e'(k)=partial_k sigma(e(k),k)/
                         [1-partial_E sigma(e(k),k)].

The existing exact e(k) and velocity bounds are retained. This says nothing
by itself about how the fuel populates that branch.

The word bound in 9ef64dc also gives |sigma_0|<=C_* and, for |d|>=1,

    |sigma_d(z)|<=C_* q^(2|d|-2),
    q=5/72, C_*=75/20636.                           (4)

Thus the Fourier series is absolutely convergent, with exponential
weights exp(mu|d|) whenever q^2 exp(mu)<1. It supplies an analytic bulk
symbol in a strip about real k, rather than erasing position phases.

## 2. Actual half-line differs by missing boundary-crossing words

Let Sigma_+(z) be the actual co-pair block of the quadrant self-energy,
and T_+(sigma) the Toeplitz compression of (1) to n,m>=0. Define

    W(z)=Sigma_+(z)-T_+(sigma(z)).                  (5)

Every quadrant virtual path occurs in the whole-line expansion with the
same amplitude and denominators. The difference consists exactly of paths
that leave the quadrant. For real E<=5/2 the entries of W(E) are
nonnegative; this is ENTRYWISE, not a general operator-positivity claim.

A missing path from (n,n) to (m,m) visits a negative coordinate. That
coordinate requires at least n+m+2 hops; the other requires at least
|n-m| hops. Thus it has at least 2max(n,m)+2 hops. The closed-word index
k counts all but the two endpoint B hops, so k>=2max(n,m).

Only even k survive. Absolute-value majorization by the positive real
path expansion at the maximal allowed real part yields, also at complex z,

    |W_nm(z)|<=C_* q^(2max(n,m))<=C_* q^(n+m),
                             Re z<=5/2.            (6)

This proves a size-uniform localized boundary kernel. Summing rank-one
matrix units gives the trace-class bound

    ||W||_1 <=sum_nm |W_nm|
       <= C_* (1+q^2)/(1-q^2)^2.                   (7)

No finite-rank assumption is made for the exact W. The two missing-path
conditions (closed-word length and reaching the boundary) control it.

## 3. Explicit first boundary coefficients

Put s=jp^2+jm^2, p=jp jm, a_1(z)=13/4-z, a_2(z)=4-z.
The k=0 boundary contribution is

    W^(0)(z)=s/a_1(z) |0><0|.                     (8)

The k=2 contribution, supported exactly within indices 0,1, is

    W^(2)(z)=1/[a_1(z)^2 a_2(z)]
          [[s^2+6p^2, 2ps],
           [2ps,       s^2-2p^2]].                 (9)

For example, for a single relative arm, writing its shift coupling as
L=aI+bS* gives the fourth-hop term (L*)^2 L^2. Compared with its
whole-line Toeplitz counterpart, the missing part is
(4a^2b^2+b^4)P_0+b^4P_1+2ab^3(|0><1|+|1><0|).
Summing arms (a,b)=(jp,jm),(jm,jp) gives (9). The primary control
independently enumerates the actual two-hop/four-hop coordinate paths.

For real E<=2, the matrix (9) is positive definite: s^2>=4p^2,
its diagonal is positive and its numerator determinant is
s^4-12p^4>=4p^4>0. Its trace norm therefore equals its trace.

The exact bulk first terms, useful for independent checking, are

    sigma^(0)(z,k)=-2v(k)/a_1(z),
    sigma^(2)(z,k)=-2v(k)^2/[a_1(z)^2 a_2(z)].       (10)

## 4. Small actual co-pair boundary defect on the first-band interval

The supplied microscopic root counterterm must be included. The co-pair
effective Hamiltonian is exactly

    2I+T_+(sigma(z))+U(z),
    U(z)=W(z)-A|0><0|.                             (11)

On I_b=[2-5J,2], the closed denominators permit sharper constants

    q_b=1/24, t=q_b^2=1/576, C_b=1/460.

The existing bound |e(k)-(Omega-2J cos k)|<J/12 and A/J=13/12 put the
entire first bulk band in [2-17J/4,2-J/12], contained in I_b. This does
NOT assert that all low-source spectral weight lies in I_b.

The remainder W-W^(0)-W^(2) requires closed index k>=4 and a boundary
crossing. Hence its entries are <=C_b t^max(2,n,m), and

    ||W-W^(0)-W^(2)||_1 <= C_b S_1(t),
    S_1(t)=9t^2+t^3[7/(1-t)+2t/(1-t)^2].          (12)

The nine pairs max(n,m)<=2 give the first term; the remaining shell r>=3
has 2r+1 entries. This bounds the COMPLETE higher-word sum.

At E=2, W^(0)=A P_0, so the physical root counterterm cancels it.
Throughout I_b,

    |s/a_1(E)-A| <=4AJ,
    ||W^(2)(E)||_1 <=(2s^2+4p^2)/[(5/4)^2*2].

Combining these exact bounds gives the proposed quantitative result

    ||U(E)||_1
      <=4AJ+(2s^2+4p^2)/[(5/4)^2*2]+(1/460)S_1(1/576)
       <J/500,               every E in I_b.       (13)

The final strict comparison is checked with exact Fractions. This is
the actual half-line co-pair boundary defect, with every remaining path
paid. It is NOT a source-coupling probability or an inverse bound. Fuel
is retained separately in SOURCE_KERNEL.md; it has not been absorbed
into an artificial bulk site or switched off.
