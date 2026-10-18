# Complete majorants and phase-sensitive packet bounds

Author-derived; independent review pending. REAL_AXIS.md supplies the
actual operator and factorization. All numerical assertions below are
strict rational comparisons in controls.py, not sampled spectral data.

## 1. Boundary word moments

Write t=u=1/576, Cb=1/460, d2=a0^2*2, and define

    S0=9t^2+t^3[7/(1-t)+2t/(1-t)^2],
    S1=18t^2+(t^3/2)[3t(1+t)/(1-t)^3
                          +23t/(1-t)^2+44/(1-t)],
    S11=36t^2+t^3[t(1+4t+t^2)/(1-t)^4
                  +12t(1+t)/(1-t)^3+48t/(1-t)^2+64/(1-t)].

They sum the COMPLETE max(2,n,m) higher-word majorant from B2-BB,
with weights 1, m+1 and (n+1)(m+1). For max(n,m)=r>=3 the shell
weights are respectively 2r+1, (r+1)(3r+2)/2 and (r+1)^3. Therefore

    M0 =4AJ+(2s^2+4p^2+4ps)/d2 +Cb S0,
    M1 =4AJ+(3s^2+2p^2+6ps)/d2 +Cb S1,
    M11=4AJ+(5s^2-2p^2+8ps)/d2 +Cb S11,
    Mr =(s^2-2p^2+2ps)/d2 +Cb(S1-S0).             (17)

These are bounds on the moments named in (10). In particular Mr is
derived by weighting each entry by m, NOT by subtracting unrelated
upper estimates of two physical norms.

Let V=(201/100)J. The majorant for sum(m+1)|dU_mn/dk| is

    Mdot1=V {s/a0^2
      +(3s^2+2p^2+6ps)/d2*(2/a0+1/2)
      +(Cb/a0)[(1+2t/(1-t))S1+2t S1']}.           (18)

Differentiate every actual positive denominator before majorizing.
The root derivative is s/a1^2; the W2 derivative factor is
2/a1+1/a2. The higher geometric sum Cb t^R has logarithmic derivative
[1+2R+2t/(1-t)]/a0. This proves (18), including all higher words.

## 2. Derivative majorants: no norm differentiation of plane waves

Let eps=(1-u)^(-2)-1, l=501/500, C=l^2, cmin=1/(1+eps),
j0=960J/961, eta=C M1/j0. Define

    A_E=(2/a0)[(1+u)/(1-u)^3-1],
    A_v=u/[v*(1-u)^3],
    A_m=(2p) A_v,
    A_dot=A_E V+A_v(2p),
    L_dot=l A_dot/[2(1-eps)],
    L_m=l A_m/[2(1-eps)],
    c_dot=A_dot/[2(1-eps)], jrel=V/a0.           (19)

Here v* (not v(k)) is the v in the denominators of (19). A_dot bounds
the W-norm parameter derivative of a_E, and A_m its first Fourier
moment. The factor 1/2 comes from even real Fourier coefficients of
log a. Thus L_dot bounds sum|dot l_j| and L_m bounds sum j|l_j|.

Use X=sup|x_n|/(n+1), and let Q=g_T U. The exact finite geometric sum
for g0 gives

    |d(g0)_nm/dk| <= (n+1)(m+1)/J_E
                       +jrel(min(n,m)+1)/J_E.

Lower triangular convolution only decreases these indices. Hence the
coordinatewise derivative matrix has the following uniform X bound:

    B=[(2l L_dot+C jrel) M1+C M11+C Mdot1]/j0.     (20)

It is NOT asserted that Q is norm differentiable on ell-infinity or X.
For each fixed Neumann order and fixed output coordinate, first truncate
all matrix indices and differentiate the finite sums. The common
ABSOLUTE kernel/moment majorants for Q,Q',tvec,tvec' bound the sums of
the differentiated terms, including their index tails, uniformly on
the central interval. Such coefficientwise majorants are obtained by
replacing every closed-word denominator by a0 and v0 by v*, then taking
absolute powers in the convergent log/exponential series. Their zeroth
and first moments are bounded by (19); (17)-(18) give the corresponding
boundary majorants. Thus this is a common summable bound, not only a
pointwise supremum of an operator norm. Dominated convergence removes those index cutoffs
and establishes coordinatewise C1 for each finite-order term. Only then
sum the Neumann orders. The r-th differentiated term is bounded by
r eta^(r-1) B ||input||_X +eta^r||input derivative||_X. Its sum is
uniform in k, so the scalar derivatives used below are justified.

On the CENTRAL interval sin k>=s0=13/20, put tvec=u/sin k. Then

    ||tvec||_X<=l,
    ||dot tvec||_X<=L_dot+l(1/s0+1/s0^2)=:t_dot,
    x=l/(1-eta), x_dot=(t_dot+B x)/(1-eta),
    z0=eta*x/cmin,
    z1=(B*x+eta*x_dot)/cmin+z0*c_dot/cmin.         (21)

These bound |zeta| and |zeta'| for psi0=c0 sin k(1+zeta). Only z0 is
used at thresholds; z1 and the derivative argument are central-band
statements. This distinction avoids a false threshold derivative bound.

## 3. Fuel denominator and density constants

Let

    eerr=1/276,
    alphaerr=max(eps/(1+eps),4J),
    R0=C^2 M0/[j0^2(1-eta)],
    derr=eerr+alphaerr+J R0.                      (22)

Here eerr is the sharper proved dispersion bound, normalized by J.
Thus |d_f|<=derr<1/100 in (14). Put fmax=1/(1-derr).

For the exact velocity write

    e'(k)=2J sin k/R(k),
    R(k)=2[1+f_E(E,v)]/[a0 f_v(E,v)].             (23)

All positive mixed derivatives of f are bounded by F in (2). If F_ab
denotes partial_E^a partial_v^b F at a0,v*, then

    Rmin=1/(1+eps),
    Rmax=(a_max/a0)(1+F10),
    P=1+F10, Bmin=2/a_max,
    R_a <= (2/a0)[P_a/Bmin+P B_a/Bmin^2],
    R_ab<= (2/a0)[P_ab/Bmin
               +(P_a B_b+P_b B_a+P B_ab)/Bmin^2
               +2P B_a B_b/Bmin^3].             (24)

In (24), P_E=F20,P_v=F11,P_EE=F30,P_Ev=F21,P_vv=F12;
B_E=F11,B_v=F02,B_EE=F21,B_Ev=F12,B_vv=F03. These are absolute
majorants, not claims that the derivatives of the ratio have one sign.

Using |v'|,|v''|<=2p, (24) gives

    R1=R_E V+R_v(2p)<1/100,
    E2=2J(1/Rmin+R1/Rmin^2),
    R2=R_EE V^2+2R_Ev V(2p)+R_vv(2p)^2
                              +R_E E2+R_v(2p)<1/100. (25)

It follows that 1.99J sin k<=e'<=2.01J sin k. Also e'' is positive
on [0,pi/3], negative on [2pi/3,pi], and e'''<0 between them:
differentiate 2J sin k/R twice, use sin k>=4/5 in the middle and
the displayed Rmin,Rmax,R1,R2. Thus e' has exactly ONE maximum.
This hypothesis of the later nonstationary-phase bound is proved,
not inferred from positivity or a speed bound alone.

The scalar fuel correction derivative is bounded as follows:

    Cder=2l L_dot+C jrel,
    Rdiff_dot=(B*C/(1-eta)^2
                   +eta*(C+Cder)/(1-eta))/j0,
    alphadot=(1+4J)(jrel+A_dot/(1-eps)),
    banddot=2 max(1/Rmin-1,1-1/Rmax),
    dder=banddot+alphaerr+alphadot+J Rdiff_dot,
    fprime=(dder+derr)/(1-derr)^2.                (26)

Specifically ftilde=-J F_f exp(-ik)=1/[1-d_f exp(ik)] has modulus
<=fmax and derivative <=fprime. The exact free cancellation is used
BEFORE differentiation; differentiating unrelated factors would miss it.

Put T=a_E(exp(ik))/(1+f_E). Then

    1/(1+F10)<=T<=1+eps,
    Tprime=A_dot+A_m+(1+eps)(F20 V+F11(2p)),
    Zmax=(1+eps) fmax l (1+z0),
    Zprime=Tprime fmax l(1+z0)
       +(1+eps) fprime l(1+z0)
       +(1+eps) fmax(L_dot+L_m+c_dot*l)(1+z0)
       +(1+eps) fmax l z1.                       (27)

These bound the actual outgoing amplitude factor

    Z=T ftilde L(exp(-ik)) c0(1+zeta).             (28)

No derivative follows merely from small ||Z-1||. Each derivative in
(27) has its own complete word or Wiener moment account.

The exact density ratio r=w/w0 equals

    r=T c0^2 |ftilde|^2 |1+zeta|^2.

Consequently

    rmin=(1-z0)^2/[(1+F10)(1+eps)(1+derr)^2]>24/25,
    rmax=(1+eps) fmax^2(1+z0)^2<26/25.           (29)

For later fuel-central variation, an explicit bound on |r'| is

    rprime=Tprime fmax^2(1+z0)^2
       +(1+eps)*2c_dot*fmax^2(1+z0)^2
       +(1+eps)*2fmax*fprime*(1+z0)^2
       +(1+eps)*fmax^2*2(1+z0)z1.                (30)

The controls verify rprime<1/8. These estimates include the Jacobian
e'(k); source density is not equated with energy-band occupancy.

## 4. Actual outgoing kernel and its GLOBAL local remainder

For 0<k<pi the exact identity is

    (g0)_nm=-exp[ik(n+1)] sin[(m+1)k]/(J_E sin k)
              +1_(m>n) sin[(m-n)k]/(J_E sin k).  (31)

Apply L on both sides and the boundary Neumann inverse. Since S is
complex symmetric, u^T(I+US .0)=psi0. Thus

    S_n0=-exp[ik(n+1)] L(exp(-ik)) psi0/(J_E sin k)+r_n.

The first term has its correct physical phase. The remainder has a
global ell2 bound over ALL n>=0, including sites near the source:

    ||r||_2 <= [L_m(1+z0)/j0]
                 +[C/(j0 s0)] [C Mr/(j0(1-eta))] =:Rloc. (32)

For the first term, the l2 norm of the omitted L-tail is at most
sum_(j>=1)sqrt(j)|l_j|<=L_m. For the second use (31),
||R_free y||2<=sum sqrt(m)|y_m|/(J_E sin k), and sqrt(m)<=m.
The m=0/root term contributes ZERO to that remainder; Mr pays the
remaining rows. Root-local corrections are not called exponentially
small at every n. The central energy width is exactly 3J, so after
F_n=-J F_f S_n0 and the spectral integral,

    ||central local packet remainder||2
                 <= (3J/pi) fmax Rloc <1/200.    (33)

Componentwise inversion plus the Fourier ell2 bound justifies this
vector statement without treating a generalized plane wave as an
ell2 vector at one energy.

The central co-pair amplitude is therefore, exactly,

    v_n^C(t)=integral_(k_a)^(k_b) exp[-it e(k)]
       {exp[ik(n+2)] b(k)+exp[-ik(n+2)] conjugate(b(k))} dk
       + local remainder,
    b(k)=sin k Z(k)/(pi i).                      (34)

The sign of the conjugate term is fixed by -Im F_n; conjugation of
b already includes 1/i. No coherent term is removed. Equation (34)
contains both directions and all hard-cut endpoint terms.

## 5. Central packet constants used by the operating proof

For |e(k)-Omega|<=1.5J, (3) gives |cos k|<=94/125, hence sin k>13/20.
Equations (23)-(25) give

    1.3 <= e'(k)/J <= 2.01.                     (35)

Using pi>31/10, integral_C|cos k|<=2-2s0 and integral_C sin k<=2,

    M_b=Zmax/(31/10),
    TV_b=[Zmax(2-2s0)+2Zprime]/(31/10),
    4M_b+TV_b <8/5.                             (36)

The integer n+2, not n, is retained in the phase. Fuel is the extra
endpoint site. For fuel-central w(k), (29)-(30) give a bounded
variation amplitude and (35), yielding |G_C(t)|<=10/(J|t|).
The conservative constant is checked with the same two-monotonic-piece
integration-by-parts estimate as (36).

Every step above concerns the actual source resolvent and its actual
spectral kernel. controls.py checks rational majorants; it does not
numerically solve a surrogate source or certify the all-size proof.
