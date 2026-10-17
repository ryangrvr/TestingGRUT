# S9 formation: an adjacent-charge rederivation and challenge ledger

Same Builder, fresh derivation and code. No failed load-bearing inequality
was found in the scope below. This is not independent external validation.
The source at `82bcab3` remains unchanged and author-derived.

## 1. Check the extension before manipulating a field

Fix a root, M=N-1. A geometry basis is the complete set of ordered words
of distinct non-root labels, including the empty word. Write K=length,
F=M-K, and S_v for append-v. A child has a unique parent and last label,
so S_v^* S_z=0 for v!=z and S_v^* S_v=P_(v free). Consequently

    S^* S=F, ||S||<=sqrt(M),
    A_u=sum_v exp(-i delta_v u) S_v,
    A_u^* A_u=F, ||A_u||<=sqrt(M).

This does not make append operators commute: S_v S_z and S_z S_v
produce different words. Both commutators [S,A_u] and [S^*,A_u]
must be bounded, each by 2M.

Use charge c=n_R-K. The physical preparation occupies c=0. For fixed
c>=-M, 0<=n_R<=M+c, so each finite-N, finite-L charge block is finite.
H_c is the restriction of

    D_word+dGamma(h_L)-lambda(S a_0^*+S^* a_0), lambda=g/N.

Define the ambient self-adjoint H by the orthogonal direct sum of H_c.
Annihilation and append both map c to c-1; creation maps c to c+1.
There is no artificial highest photon sector in the calculation.
The relative ambient Hamiltonian need not be bounded below. Adding the
physical conserved-charge offset restores the declared positive receiver
frequencies and leaves the prepared-sector evolution unchanged. The
relative charge-block calculation supplies no additional initial energy.

## 2. Derive the contact identity in Schrodinger space

Let Phi(t)=exp(-i H_0 t) psi_0, receiver vacuum initially, and
z_r(t)=a_r Phi(t), a vector in charge -1. Direct multiplication gives

    a_r H_0-H_-1 a_r=sum_s (h_L)_(rs) a_s-lambda S.

Thus i z'=H_-1 z+h_L z-lambda S Phi. The receiver-mode matrix acts on
the external r index and commutes with H_-1. Variation of constants gives

    a_0 Phi(t)=i lambda integral_0^t C_L(u)
                    exp(-i H_-1 u) S Phi(t-u) du.                 (F1)

Multiplying by exp(i H_-1 t) yields exactly S9's Heisenberg identity.
The free-field term vanishes because of the initial vacuum, not because
the field is reset or later factorized. In particular this argument
does not assume that different histories emit independently.

The fresh block code tests the operator intertwiner and (F1) separately,
including the entire charge -1, 0 and +1 blocks. It also verifies
||a_0^* Phi||^2=1+||a_0 Phi||^2. Dropping this +1 is a rejected mutation.

## 3. Recheck both receiver estimates

For the half-line path, the boundary spectral measure obtained from the
sine transform is rho(x)=sqrt(4-x^2)/(2 pi), -2<x<2. At u>=1 put
a=2-1/u. The two removed endpoint intervals have total mass at most
4/(3 pi) u^(-3/2). Two integrations by parts on [-a,a] give

    2 rho(a)/u+[2|rho'(a)|+integral_-a^a |rho''|]/u^2
       <=[2/pi+4/(pi sqrt(3))] u^(-3/2).

Here rho''<0, integral |rho''|=2|rho'(a)|, rho(a)<=1/(pi sqrt(u)),
and |rho'(a)|<=sqrt(u)/(pi sqrt(3)). The total coefficient is <2.
For all u, therefore |C_infty(u)|<=min(1,2u^(-3/2)). The four integral
consequences in S9 follow by elementary integration; the L1 bound 5 is
deliberately loose. No oscillatory cancellation is needed in these bounds.

A closed walk from mode 0 first distinguishes the finite L chain from
the half-line at order 2L. Both one-particle norms are <=2. If
L=64N^2 and u<=8N^2, the absolute Taylor remainder is bounded by

    2 sum_(n>=128N^2) (2u)^n/n!
      <=(16/7)(e/8)^(128N^2)<2^(-120N^2)=epsilon_N.                (F2)

The ratio of successive Taylor terms is <1/8. Using e<3 proves the
last inequality already at N=1; its margin increases with N.
The contact estimate from (F1) is consequently

    ||a_0(t)psi_0||<=lambda sqrt(M)(5+8N^2 epsilon_N)
                     <=6 lambda sqrt(M)=b_N.                    (F3)

This compares only free kernels. It never substitutes an infinite-bath
many-body evolution for the finite physical evolution.

## 4. Challenge the retarded replacement, retaining creation

At fixed t and u<=t differentiate the vector
exp(iH(t-v)) A_(u-v) exp(-iH(t-v)) psi_0. Its explicit phase derivative
cancels the D_word commutator, leaving -i[V,A_(u-v)](t-v) psi_0.
There are two terms, proportional to [S,A] a_0^* and [S^*,A] a_0.
The latter is controlled by b_N; the former by sqrt(1+b_N^2).
Both estimates are evaluated on the actual evolved initial state at
time t-v. Geometry and field commute at equal times, including between
the adjacent charge blocks. Therefore

    ||(S(t-u)-A_u(t))psi_0||<=min(B_N u,2sqrt(M)),
    B_N=2lambda M(sqrt(1+b_N^2)+b_N).                             (F4)

For N>=64, b_N<=33/40 and B_N<187/40<5. With tau=sqrt(N), integration
of the first bound up to tau and the second afterward gives

    integral_0^t |C_infty(u)| min(B_N u,2sqrt(M)) du
       <=4B_N sqrt(tau)+8sqrt(M)/sqrt(tau)<=28N^(1/4).             (F5)

This remains valid when t<tau. Replacing [S,A] by zero fails on the
initial seed with distinct detunings; the fresh control exhibits it.

## 5. Recompute the sign, rate and every remainder

The exact current is

    K'=2lambda Im <S(t)psi_0,a_0(t)psi_0>.

Inserting (F1), using (F2),(F5), and extending only the instantaneous
kernel integral to infinity yields three errors in the field vector:

| Vector error | Bound | Contribution to the current |
|---|---:|---:|
| Full history retardation | 28lambda N^(1/4) | 56g^2 N^(-5/4) |
| Initial missing kernel tail | 5lambda sqrt(M)/sqrt(t) | 10g^2/(N sqrt(t)) |
| Finite receiver kernel | 8lambda sqrt(M)N^2 epsilon_N | 16g^2 N epsilon_N |

Each current column multiplies the vector column by 2lambda sqrt(M).
The exact instantaneous product is S^* A_u=sum_v exp(-i delta_v u)
P_(v free). Abel regularization of the spectral measure and its L1
kernel imply

    2 Re integral_0^infty C_infty(u) exp(-i delta u) du
       =sqrt(4-delta^2).

The result is precisely S9's current and remainder (14). In particular
the rate has positive sign for every |delta|<=3/4. This is an equation
for one expectation with a controlled remainder, not a full-state
Markov approximation.

## 6. Expectation to probability: check the direction of every inequality

Since sqrt(4-delta^2)>=sqrt(55)/4>gamma_0=9/5, f=<F> obeys

    f'<=-gamma_0(g/N)^2 f+|r_N|, f(0)=M.

Keeping exponential damping on the constant history error gives
(56/gamma_0)N^(3/4)=(280/9)N^(3/4). Integrating the other errors
absolutely through 8N^2 gives bounds <73 and <155N^3 epsilon_N.
The t^(-1/2) singular upper bound is integrable. Start at positive time
and use continuity at zero. There is no neglected initial layer.

At t>=4N^2, gamma_0 g^2 t/N^2>=729/125 and exp(729/125)>256.
Macro failure, m<ceil(N/2), implies F>=floor(N/2)+1>=N/2.
The spectral distribution of the positive observable F therefore gives

    1-p<=2f/N
      <1/128+(560/9)N^(-1/4)+146/N+310N^2 epsilon_N.               (F6)

For N>=2^28 this is <1/2. The fresh exact-rational check verifies a
lower bound 0.5060758449965053 at that endpoint. The proof is pointwise
throughout [4N^2,8N^2], for every root and each allowed deterministic
disorder assignment. Original root averaging and the frozen uniform
time sample preserve the bound without conditioning on a favorable root
or instant. No all-size small-N guarantee follows.

## Outcome of this challenge

The adjacent-charge derivation reproduces the load-bearing identities
and estimates without importing the earlier diagnostic implementation.
The endpoint tail, vacuum assumption, noncommuting creation term and
negative-damping step survived this check. No formation claim is upgraded
to external verification. An independent mathematical reviewer should
still check (F1)-(F6) and their operator-domain interpretation. The fresh
finite diagnostics test execution and algebra; they do not prove all N.
