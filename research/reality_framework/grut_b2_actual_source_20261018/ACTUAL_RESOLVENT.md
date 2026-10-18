# Actual source: exact block operator, spectral measure and receiver kernels

Author-derived; independent mathematical review pending. This is the
actual nominal h, not the B2-SB scalar countermodel or a supplied packet.
The off-real identities below are complete. Their real-axis estimates
needed for operation remain unresolved.

## 1. Exact reindexing with all charge separations retained

Write n=min(i,j), r=i-j. This is a bijection from the charge quadrant to
N_0 x Z: for r>=0, (i,j)=(n+r,n); for r<0, (i,j)=(n,n-r).
Let K=ell2(Z), e=e_0, and q_r=2+|r|+(1/4)1_(|r|=1).
Define the self-adjoint on-layer operator D on dom(q) by

    D_rr=q_r,
    D_(r+1,r)=D_(r,r+1)=-jp,  r>=0,
    D_(r-1,r)=D_(r,r-1)=-jm,  r<=0.                  (1)

Define the bounded, generally non-self-adjoint inter-layer operator T by

    T_(r+1,r)=-jm, r>=0,
    T_(r-1,r)=-jp, r<=0.                             (2)

There are no other entries. With psi_-1=0, the bare charged operator is

    (H_b psi)_n=D psi_n+T psi_(n+1)+T* psi_(n-1).     (3)

Every primitive charge hop appears exactly once with its adjoint. In
particular r=0 connects to both arms in (1) and (2); the n=0 boundary
removes exactly the unphysical negative-event hops. It inserts no spatial
interaction or new assembly grammar. If b=|n=0> tensor e, then

    H_ch=H_b-A|b><b|,
    h=Omega|f><f|+H_ch-J(|f><b|+|b><f|).             (4)

The diagonal q is self-adjoint, and everything added to it is bounded.
Finite-support vectors are a core. This establishes domains and
self-adjointness, rather than treating the infinite block as a finite
matrix. It includes all r, not a relative-distance truncation.

An important exact obstruction to ordinary transfer matrices is

    T*T e_r = jm^2 e_r              (r>0),
    T*T e_0 = (jp^2+jm^2)e_0,
    T*T e_r = jp^2 e_r              (r<0),
    ||T||=sqrt(jp^2+jm^2)=sqrt(13)/192,
    ker T* = span{e_0, jp e_1-jm e_-1}.              (5)

T is injective and bounded below, but not onto; T inverse cannot be used.
The infinite transverse space also falls outside a theorem whose
hypotheses require fixed finite nonsingular blocks. Those hypotheses
would need verification in a different representation, not assertion.

## 2. The physical matrix Weyl function

For Im z>0, let P_0 identify layer 0 with K, and define

    M(z)=P_0(z-H_b)^-1 P_0*,
    m_0(z)=<e,M(z)e>.                               (6)

Removing layer 0 leaves an exact translated copy of H_b. Schur inversion
therefore gives

    M(z)=[z-D-T M(z) T*]^-1.                        (7)

The definition (6), not merely the algebraic equation (7), selects the
physical branch. It has Im M(z)<=0 and ||M(z)||<=1/Im z in the convention
(z-H)^-1. No sign convention (H-z)^-1 is silently substituted.

An equivalent selection is the strong-resolvent limit of inherited
compressions Gamma(N,R)={0<=n<=N, |r|<=R}, as N,R grow cofinally.
Extend each finite compression by the uncoupled diagonal q outside it.
For any finite-support psi, the compressed operator agrees with H_b psi
once its support and its neighbors fit. On the dense set (z-H_b)core,
the compressed resolvents then converge to the true resolvent; their
uniform 1/Im z bound extends this to the whole Hilbert space.
This proves a strong off-real limit, not a rate uniform as Im z goes to
zero, not norm convergence of spectral measures, and not convergence
uniform in n proportional to L.

## 3. The actual scalar fuel spectral measure

The rank-one root counterterm and fuel Schur complement give exactly

    m_ch(z)=<b,(z-H_ch)^-1b>=m_0(z)/(1+A m_0(z)),
    F(z)=<f,(z-h)^-1f>
        =[1+A m_0(z)]/Delta(z),
    Delta(z)=(z-Omega)[1+A m_0(z)]-J^2 m_0(z).       (8)

The first fraction is used where its factors are finite; the full
resolvent definition gives the continuation. Do not count a zero of
1+A m_0 as a fuel pole without checking the numerator and denominator.

The actual positive measure nu_f is recovered by Stieltjes inversion:

    nu_f({E})=lim_(eta down 0) [-eta Im F(E+i eta)],
    rho_f(E)=-(1/pi)lim_(eta down 0) Im F(E+i eta)
              at its Lebesgue density points.       (9)

The remaining singular continuous measure is not zero by definition.
If the finite boundary value m_0(E+i0)=a(E)-i b(E) exists, then

    rho_f(E)=J^2 b(E)/[pi |Delta(E+i0)|^2].         (10)

This is a direct fractional-linear imaginary-part identity: writing
x=E-Omega, Delta=x+(Ax-J^2)m_0 makes the determinant J^2 explicit.
If b(E)>0 at a regular energy, Delta cannot vanish. Its imaginary part
would force Ax=J^2, where Delta=x=J^2/A is nonzero. Nevertheless,
positivity for almost every E does not exclude singular measure on an
exceptional null set. Uniform boundary control or an applicable theorem
is required there, including thresholds.

On a real interval where m_0 is analytic and real, a simple zero E_b
of Delta with 1+A m_0(E_b)!=0 is a fuel-visible atom with weight

    w_b=[1-J^2 m_0'(E_b)/(1+A m_0(E_b))^2]^-1.      (11)

Here d/dz (z-H)^-1=-(z-H)^-2. Formula (11) is a candidate residue
test where its analytic hypotheses hold, not a solved catalogue of atoms.
Inside a band, boundary regularity and normability must be established
before making a pole claim. Resonances require an actual continuation;
no resonance locations or absence theorem is supplied by (8).

B2-SS's no-eigenvalue-at/below-e0 result is retained. Fuel-visible
embedded components, possible above-band atoms, singular continuous
weight, edge behavior and the numerical density in (10) remain
unclassified. A band included in essential spectrum does not fix any
of those source weights. Nominal parameters are fixed; no robustness
over a new parameter interval is claimed.

## 4. The actual spatial resolvent, not only the scalar measure

Tail elimination at each n gives the off-diagonal bare resolvent block

    P_n(z-H_b)^-1P_0*=(M(z)T*)^n M(z).

Combining this with (8), the EXACT physical co-pair resolvent is

    F_n(z)=<n,n|(z-h)^-1|f>
          =-J <e,(M(z)T*)^n M(z)e>/Delta(z).        (12)

All source return, boundary scattering, relative excitations and both
propagation directions are retained. The physical co-pair event is
precisely r=0 at layer n; no velocity-based detector has been substituted.

Define the signed cross-spectral measure
kappa_n(B)=<n,n|1_B(h)|f>. The real matrix coefficients make it real;
where finite boundary limits exist, its ac density is

    kappa_n,ac'(E)=-(1/pi)Im F_n(E+i0).             (13)

Residues and singular measures must also be included. Exactly,

    v_n(t)=integral exp(-itE)d kappa_n(E),
    A_L(t)=sum_(n=L/8)^(3L/8-1) |v_n(t)|^2.        (14)

The integrals are coherent amplitudes; energy components are not
dephased. The general expression uses the upper/lower resolvent jump;
the imaginary-part form (13) uses this real source basis. Fuel return is

    G(t)=sum_b w_b exp(-itE_b)
          +integral exp(-itE)rho_f(E)dE
          +integral exp(-itE)dnu_f,sc(E).           (15)

Equations (8)-(15) put both missing physical functions under the SAME
actual resolvent. None has been evaluated or uniformly bounded here on
the late window. The unresolved operation is their quantitative inversion
with n proportional to L and t proportional to L/J, including exceptional
energies. Scalar m_0 controls (15); (14) additionally needs the matrix
sequence in (12). That distinction survives all moment estimates.

## 5. What finite operators do and do not represent

The finite Gamma(N,R) matrices are inherited compressions of h, with root
-A only. At fixed R their bare tail functions are DIFFERENT by layer:

    M_N=(z-D_R)^-1,
    M_n=[z-D_R-T_R M_(n+1)T_R*]^-1.                (16)

For their site resolvents replace (MT*)^n M in (12) by

    M_n T_R* M_(n-1) T_R* ... M_1 T_R* M_0.        (17)

Using a single fixed M in finite depth is an error. The controls test
(16)-(17) against direct coordinate assembly with exact arithmetic.

A finite PHYSICAL path has its original far-end energy and two-source
constraints/interactions; it is neither this finite compression nor the
infinite bare tail. It has a discrete spectrum. Its relation to (14)-(15)
on the frozen horizon is B2-SO's full same-space comparison, not equality
of spectra or omission of reflected paths. That comparison already pays
finite-host effects on its stated time range, if its proof survives review.
No finite diagonalization can itself exclude half-line singular spectrum.

## 6. Why the elementary off-real route stops quantitatively

For two trial Weyl functions in the resolvent norm ball, the Schur map
in (7) has the generic Lipschitz estimate

    ||Phi_z(X)-Phi_z(Y)||
        <= (||T||^2/eta^2)||X-Y||, eta=Im z.        (18)

This estimate proves contraction only when eta>||T||, that is
eta/J>20 sqrt(13)>72. It does not prove that the true branch fails to
converge at smaller eta, merely that THIS estimate gives no control.

The Cauchy-smoothed fuel measure
nu_eta(E)=-(1/pi)Im F(E+i eta) obeys exactly

    integral exp(-itE)nu_eta(E)dE=exp(-eta|t|)G(t).  (19)

Using only |G|<=1, its uncorrected Fourier error is at most eta|t|.
To guarantee that generic error <=1/1000 through T_L, one can take
eta/J<=6/(1000L). Already at L=2016 this is <=1/336000, incompatible
with the certified contraction range >72. Alternatively undamping (19)
amplifies absolute errors by exp(eta|t|); in the contraction range it
exceeds exp(18144) at the earliest target time theta=252.

Thus small finite off-real checks plus this norm bound cannot provide
the required real-axis measure or uniform late inversion. More specific
limiting-absorption/scattering or validated spectral approximation
estimates could succeed. No theorem of algorithmic impossibility, no
claim that external methods fail, and no physical recurrence are inferred.

## 7. A named part of the capture budget is already paid

For the central packet c=exp(-ith)1_C(h)f, w_C=||c||^2, the existing
invariant graph has ||X||<=1/19. Since c=(p,Xp), its separated-charge
fraction is at most

    ||Q_sep c||^2/w_C <=(1/19)^2/[1+(1/19)^2]=1/362.

The exact detector loss is the sum of three ORTHOGONAL physical sectors:

    1-eta_C = [|<f,c>|^2
       +sum_(n outside R_left) |<n,n|c>|^2
       +||Q_sep c||^2]/w_C.                       (20)

Therefore the previously sufficient eta_C>=24/25 follows if the ACTUAL
resolvent inversion proves

    [|integral_C exp(-itE)dnu_f(E)|^2
      +sum_(n outside R_left)
        |integral_C exp(-itE)d kappa_n(E)|^2]/w_C
            <=1/25-1/362=337/9050                 (21)

uniformly on I_L. This isolates fuel/source occupation and missed spatial
capture; transverse detector dressing is not an unconstrained O(1) loss.
The sharp central spectral cutoff can change initial spatial tails, so
the weighted bound for exp(-ith)Pi f cannot simply be copied to c.

Equation (21), or a direct bound (14)>=17/50, together with (15)<=2/3
in modulus on the full retention interval, is the focused unsolved
actual-source problem. No new scalar-only proof route is opened.
