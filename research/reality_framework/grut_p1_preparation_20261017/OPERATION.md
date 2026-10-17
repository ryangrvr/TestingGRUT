# Positive-electric-energy operation of the prepared state

Author-derived; independent proof review pending. This proof is for the
specified controlled P1 law. Local generator support and finite electric
capacity are primitives. No real-energy consistency-sector gap is assumed.

## 1. Keep the links when proving locality

D1's matter-only tree isometry cannot be reused as a local tensor model:
E_x^2 becomes a square of a cumulative signed-charge sum after solving Gauss.
Instead, on a recorded prefix retain the three-level electric link as a
local tensor factor, together with the two hard-core matter slots. Put the
outgoing link at the cell's left endpoint. Every bond term then acts on
two consecutive cells, with the electric energy assigned to that bond.
The Gauss physical subspace is invariant under these explicit ambient
operators, so a bound on the ambient local algebra also bounds its physical
restriction. No nonlocal projection is used to implement the Hamiltonian.

The bond norm is bounded by

    A=1+j_++j_-+2kappa=293/192<2.

Each species' Hermitian hopping on one bond is a partial isometry plus its
adjoint between disjoint occupation sectors, with norm <=j_sigma. Partial
electric shifts have norm one. The field term has norm one and density
term norm <=2kappa. Thus the local interaction sum at any cell is bounded
by J_0=4, range R=1, uniformly over all particle sectors. Conserved mass
and passive fuel energies commute with the operating dynamics and can be
removed exactly from the propagation calculation.

For K_xy=sum_(I contains x,y)||h_I||, weighted row sum at mu>0 is at most
kappa_mu=2e^mu J_0. Summing connected overlapping-support chains with their
time simplex factors gives, for disjoint bounded local physical S,O,

    ||[S,U_in(t)^* O U_in(t)]||
      <=2||S||||O|| min{1,s e^(-mu d)[e^(2kappa_mu|t|)-1]},    (19)

s=min(|supp S|,|supp O|), distance between cell supports. This is the usual
finite-range Lieb-Robinson argument applied to a supplied local generator
class, not a new derivation of universal locality. An electric observable
on bond x,x+1 is counted at cell x; support-based distances differ from the
corresponding event-set distance by at most one. The statement includes
diagonal and information observables, not only particle motion.

## 2. Matter does not reach the empty end arbitrarily quickly

For prefix ell and initial matter support x<=b, n_x<=2, define
F_mu=sum_x e^(mu(x-b))n_x and C_mu=2/(1-e^(-mu)). Its initial quadratic
form is <=C_mu independently of K. Each bond conserves its local N_I;
the electric energy and density terms commute with F_mu. Within a fixed
N_I block, subtract w_min N_I from F_mu,I. Its positive residual has norm
<=N_I(w_max-w_min), hence

    |<i[h_I,F_mu]>|<=2||h_I||(e^mu-1)<F_mu,I>.

Summation and Gronwall give

    <F_mu(t)><=C_mu e^(v_mu|t|), v_mu=2J_0(e^mu-1).            (20)

This form argument also holds with local gauge factors and in the invariant
physical subspace. All coherences in the prepared state are allowed. If
L=ell-2-b, the amplitude for any carrier in the last two prefix cells is
at most sqrt(C_mu) exp[(v_mu|t|-mu L)/2]. It has no extensive K prefactor.

## 3. Full Hamiltonian and backreaction, not only a fixed graph

On the full physical space define K_prefix to retain the recorded ell
prefix and all matter within it. Every such state has zero net prefix
charge and zero fields in the tail and through the cut. On this subspace
use a reference with all prefix matter/link dynamics, all assisted tail
geometry/receiver/storage dynamics reflected at ell, no bare geometry,
and no cut-crossing matter. On its orthogonal complement use the actual
Hamiltonian compression, reference cross blocks zero. This is self-adjoint
on the **same** finite physical Hilbert space as H_op.

The prefix and tail factors commute. Prefix storage is a spectator, fuel
is passive, and the actual conditional receiver state is kept. Inconsistent
storage states are still in the physical space; they are not erased by
this comparison. For an initial word length m let a=m-ell. The entire
changing tail norm is <=g_max=1.1, so its end amplitude is <=tau_a(g_max t).

The exact uncanceled boundary drive comprises:

* Assisted pop breaking the retained prefix: <=g_max tau_a(g_max t).
* Bare geometry: norm <=4epsilon=1, and on the reference it acts only when
  head/parent contains matter in the last two prefix cells.
* Cut-crossing matter: norm <=R J_0=4 and vanishes without a boundary
  carrier. Positive electric energy on the cut is zero on the reference,
  since its cut flux is zero. Tail matter and field terms vanish there.

Everything else, including the extensive electric/mass/storage/receiver
diagonal and the extensive bulk matter Hamiltonian, cancels. Exact Duhamel
with the physical full H_op on its left therefore gives

    ||(U_op(t)-U_ref(t))X_b||<=D(t),
    D(t)=min{2,tau_(a+1)(1.1|t|)
       +(10 sqrt(C_mu)/v_mu)e^(-mu L/2)[e^(v_mu|t|/2)-1]}.    (21)

This controls all geometry excursions, matter backreaction, storage errors,
receiver returns and interference. The reference is not the physical law.
In particular this does not treat positive electric energy as an omitted
perturbation, or turn a fixed-energy Schur bound into dynamics.

At mu=1, b=N/16-1, ell=N/4, m>=N/2, J_0=4,

    L=3N/16-1, a>=N/4, v=8(e-1).

For 0<=|t|<=N/65536 and N>=256,

    D(t)<=exp(-N/5)+3exp(-N/12).                             (22)

One elementary check uses tau_k(x)<=exp(e x-k), e<3, a>=N/4.
For the second term use C<4, v>8, v<16. Its prefactor times e^(1/2)
is <(5/2)e^(1/2)<5, and its exponential rate is
-3N/32+N/8192. To obtain the stated constant 3 without numerical rounding,
use C<16/5 and v>40/3 (e>8/3), giving prefactor <27/20;
(27/20)e^(1/2)<9/4<3. The rate is <=-N/12. These elementary rational
bounds suffice uniformly over N; (21) is the sharper evaluated expression.

## 4. Compose with the actual generated matter state

The reference preparation has all matter in x<=b, but the actual preparation
is only within e_p of it. Thus the complete preparation-plus-operation
state error is

    E(t)<=e_p+D(t).                                         (23)

No matter measurement or buffer postselection is added. Prefix loss and
storage inconsistency probabilities are each <=min{1,E(t)^2}. Bounded
observable error is <=2||O||E(t). The joint conditional probability of
extensive matter and intact prefix is at least

    1-exp(-289N/65536)-2e_p-E(t)^2.                          (24)

For the original-trial bound multiply by p_form, leaving rejected trials
at zero; the use of S9 and its review status must be stated.

For pulse exp(-i theta S), chi=i<[S,U_op(t)^*OU_op(t)]>, sources preserving
the original support (or including their bounded images in b) obey

    |chi-chi_ref|<=4||S||||O||[e_p+D(t)].                     (25)

Together with (19) this is a complete information-response bound for the
supplied local operating law. It includes neutral charge-density products
and local physical field observables. Arbitrary constraint-preserving
nonlocal interactions remain excluded by the supplied generator rule.
For unit norms, supports <=2, distance >=N/32, the same use window yields

    |chi|<=4exp(-N/40)+4exp(-N/4)
                         +4exp(-N/5)+12exp(-N/12).           (26)

This is not a relativistic light cone, a zero-temperature phase, or an
all-time global-geometry stability theorem.

## 5. Nonzero response of the particles actually created

No external particle insertion is needed. Let e=(0,1) be the first matching
bond. For the plus species define its gauge-covariant forward hop K_+, and
S_+=(K_++K_+^*)/2, O_+=n_+,1. For minus use the reverse hop K_- from 1 to
0, S_-=(K_-+K_-^*)/2, O_-=n_-,0. These are local **bond sources**, not
single-event number phase sources. Their norms are <=1/2 and the prepared
reference state has Var(S_sigma)=p_e/4>=49/256.

At time zero chi=0. Since only the same-species hopping contributes to
[H,O], the local algebra gives

    chi_ref'(0)=-p_e j_sigma.                               (27)

Other same-species two-hop expectations vanish: minus occupations determine
the initial plus positions and vice versa, and those operators leave the
opposite species unchanged. Pair/vacuum coherence is not dropped. Electric
and density terms commute with O and so do not change the linear coefficient.

The local interaction sum J_0=4 gives ||ad_H^2(O)||<=8J_0^2=128.
The Taylor integral remainder for O(t) is <=64t^2; since ||S||<=1/2,

    -chi_ref(t)>= (49/64)j_sigma t-64t^2.                    (28)

At the fixed t_0=10^-5 the right side for the slower species is greater
than 7.3e-8. For N>=256, (8) and the sharp (21) at this t_0 make the
correction in (25) less than 2e-14. Thus both created species have
|chi(t_0)|>7.3e-8. This is an analytic, conservative local-response bound;
there is no claim of finite-density ballistic lower propagation over O(N)
distance, universal speed, free plasma, or new empirical law.
