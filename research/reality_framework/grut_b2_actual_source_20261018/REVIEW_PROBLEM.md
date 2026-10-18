# Specialist-review problem: this source, these detectors, these times

**Decision requested:** prove the actual operating gate, certify its
failure, or isolate a genuinely missing estimate. No new Hamiltonian.
This is a prepared request for mathematical assistance; nobody has been
contacted and no external review has occurred.

## Self-contained operator and target

Let h act on C f plus ell2(N_0 x N_0). Set

    jp=1/64, jm=1/96, kappa=1/4,
    A=13/46080, J=1/3840, Omega=2-2A,
    h_ff=Omega, h_f,00=h_00,f=-J,
    h_ij,ij=2+|i-j|+kappa 1_(|i-j|=1)-A 1_(i=j=0),
    h_(i+1,j),(i,j)=-jp, h_(i,j+1),(i,j)=-jm

and their adjoints. There are no other entries. Its domain is the
diagonal-potential domain; hopping and fuel terms are bounded. The seed
is f. No charged state, fuel return or spectral component is discarded.

For L in 24N, L>=2016, t in I_L=[L/(8J),L/(6J)], define

    A_L(t)=sum_(n=L/8)^(3L/8-1) |<n,n|exp(-ith)f>|^2,
    G(t)=<f|exp(-ith)f>.

One sufficient request is A_L>=17/50 throughout I_L and |G|<=2/3 for
2/J<=t<=L/(6J), uniformly for this family. These are sufficient SOURCE
conditions, not the physical gates and not necessary truths.

The unchanged finite two-source H_L operating gate is joint outgoing
co-pairs >=1/10 on I_L and q=2 >=1/4 on [2/J,L/(6J)], for original
double-endpoint-fuel/matter-vacuum preparation and any two-qubit/reference
state. The successful information channel under this fixed endpoint law
is exact. Moving geometry and new disorder/parameter classes are outside
this task.

## Known bridge to the physical operating experiment

Use Pi=1_(-infinity,2.5)(h), u=exp(-ith)Pi f, omega=||(1-Pi)f||^2,
sigma=sqrt(omega), M=L/2-2, C_M retaining max(i,j)<=M and fuel. Define
a=||R_left u||^2, m=||P_matter C_Mu||^2. B2-SO derives, pending external
review, a same-Hilbert-space product comparison with

    sigma<=6/27613,
    delta_L=12/27613+1500 exp(-L/105)<.001,
    tau=||(1-C_M)u||<1.01 exp(.3-L/105),
    a>=[sqrt(A_L)-sigma]_+^2,
    m=1-omega-|<f,u>|^2-tau^2,
    |<f,u>-G|<=omega,
    [a-delta_L]_+^2<=p_out<=min(1,(a+delta_L)^2),
    [m-delta_L]_+^2<=p_q2<=min(1,(m+delta_L)^2).

The norm estimate is uniform on the complete [0,L/(6J)] horizon and
includes finite-host scattering, collisions, return and interference.
Thus the sufficient source conditions give joint p_out>0.11458225 and
p_q2>0.3075. Do not transfer ONE-source 34% as joint 34%.

To certify a physical failure, use an UPPER event bound, e.g. a+delta_L
strictly below sqrt(.1), or m+delta_L below .5, at an allowed t. Failure
of the 34% or 96% sufficient source condition alone does not refute the
physical task. A source-only candidate violating point is not sufficient
unless its implication for the physical event is controlled.

## Exact actual-resolvent representation available for the proof

The coordinates n=min(i,j),r=i-j give a half-line block operator on
ell2(N_0;ell2(Z)). D has diagonal 2+|r|+kappa 1_|r|=1, positive-arm
off-diagonals -jp and negative-arm off-diagonals -jm. T has only
T_(r+1,r)=-jm for r>=0 and T_(r-1,r)=-jp for r<=0, with
(H_b psi)_n=Dpsi_n+Tpsi_(n+1)+T*psi_(n-1), psi_-1=0.
Let e=e_0,

    M(z)=P_0(z-H_b)^-1P_0*=[z-D-TM(z)T*]^-1,
    m_0=<e,M e>,
    Delta=(z-Omega)(1+A m_0)-J^2 m_0,
    F(z)=<f,(z-h)^-1f>=(1+A m_0)/Delta,
    F_n(z)=<n,n|(z-h)^-1f>
           =-J<e,(MT*)^nM e>/Delta.

These are off-real identities with the branch fixed by the actual
resolvent, not a freely chosen nonlinear solution. T is not onto and
has a two-dimensional adjoint kernel; D is unbounded with an infinite
transverse space. A finite nonsingular-block theorem needs new hypotheses.

At regular m_0(E+i0)=a-i b,

    rho_f(E)=J^2 b/[pi|Delta(E+i0)|^2].

Find or bound every fuel-visible atom and any singular continuous part,
including exceptional/threshold energies, then invert both F and F_n
coherently with constants valid at the required growing positions/times.
Eventual decay alone does not give a quantitative late-window bound.

## A narrower spatial estimate, if useful

For C_en=[Omega-1.5J,Omega+1.5J], w_C=||1_C_en(h)f||^2>=5/9. The
existing invariant graph pays separated-charge fraction <=1/362. Thus
the sufficient 96% central capture follows from

    [|integral_C_en exp(-itE)dnu_f(E)|^2
      +sum_(n outside receiver)
          |integral_C_en exp(-itE)d kappa_n(E)|^2]/w_C
                  <=337/9050,

where kappa_n(B)=<n,n|1_B(h)f>. This bound is NOT established. Its hard
energy cutoff can alter initial position tails. Do not infer it from
bulk velocities or a cone bound for the uncut preparation.

## What has been rejected and why assistance is needed

- Moment/energy/early-time scalar information alone is insufficient:
  B2-SB's exact recurrent countermodel refutes that inference, not h.
- No point eigenstate at/below the mobile-band bottom is known, but this
  does not classify embedded/above-band or singular preparation weight.
- A generic Weyl contraction bound certifies only Im z/J>72. To keep
  generic Cauchy smoothing error <=.001 at L=2016 one would use
  Im z/J<=1/336000. The norm-only estimate cannot bridge this gap.
- Exact small-compression resolvent equalities cannot replace a boundary
  limiting-absorption estimate, source scattering data or uniform inversion.

The assistance needed is a quantitative spectral/scattering theorem or
validated real-axis computational method for THIS singular block operator,
with a finite-time remainder adequate for the displayed probabilities.
No assertion is made that all such methods fail. No further automated
candidate construction or moving-geometry composition is justified here.

## Method reference, not an inherited theorem

Matthew J. Colbrook, *Computing Spectral Measures and Spectral Types*,
Commun. Math. Phys. 384, 433-501 (2021),
[primary publication](https://link.springer.com/article/10.1007/s00220-021-04072-4),
[author preprint](https://arxiv.org/abs/1908.06721).
It develops controlled resolvent/spectral methods for sparse operators and
distinguishes truncation from real-axis limiting questions. It is a
possible methodological starting point, not a proof of our source bounds
or an impossibility theorem for this specific problem. No implementation
of its algorithms or verified applicability claim is supplied here.
