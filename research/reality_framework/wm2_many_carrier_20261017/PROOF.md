# Finite-density recorded-interior theorem

**Author-derived, independent mathematical review pending.** Assumes D1 and
its supplied relational local interaction class from MODEL.md. No real-energy
gap, irreversible receiver approximation, or global extensive-H Dyson norm
is used. All statements are conditional on a recorded word and matter setup.

## 1. Uniform connected propagation inside a fixed prefix

On a retained prefix, the physical tree isometry gives a local hard-core
matter space, dimension 12 per site, with separate species-number conservation.
Write H_in=D_on+sum_I h_I. On-site D_on rotates only that site; interaction-
picture h_I(t) has the same support and norm as h_I. Define the nonnegative
matrix K_xy=sum_(I contains x,y)||h_I||. For range R and bound M5,

    sup_x sum_y exp(mu|x-y|) K_xy
       <= (R+1) exp(mu R) J_0 =: kappa_mu.                    (P1)

The usual Duhamel/Jacobi commutator recursion sums chains of overlapping
interaction supports. At order n, choose one intersection site at each step;
the nonnegative sum is bounded by sum_(x in X,y in Y)(K^n)_xy. This overcounts
chains rather than omitting them. The simplex integral supplies (2|t|)^n/n!.
For disjoint X,Y the zeroth term vanishes. P1 then gives, with s=min(|X|,|Y|),

    ||[S,U_in(t)^* O U_in(t)]||
      <=2||S||||O|| min(1, s exp(-mu d)
                                 [exp(2kappa_mu|t|)-1]).       (P2)

The symmetry follows by exchanging source/receiver and reversing time.
This is a connected bound on the full local site algebra, uniform over all
hard-core number sectors. It is not a bound using ||sum h_I||. The recursion
is the finite-range specialization of the established Lieb-Robinson argument;
the new moving-end obligation is addressed below, not absorbed into P2.

## 2. An empty-buffer estimate without a particle-number prefactor

All carriers initially lie at positions x<=b. Define n_x=n_Ax+n_Bx<=2 and

    F_mu=sum_(x in prefix) exp(mu(x-b)) n_x,
    C_mu=2/(1-exp(-mu)).

On the entire initial support subspace, F_mu<=C_mu in quadratic form,
independently of the number of occupied sites. The negative exponent towards
the root makes the geometric sum converge even for a fully occupied interior.

For an interval I, [h_I,N_I]=0, where N_I=sum_(x in I)n_x. In each N_I=n
block subtract w_min N_I from F_mu,I. The residual is positive and has norm
<=n(w_max-w_min). Cauchy-Schwarz gives, as an expectation bound for arbitrary
states (and hence quadratic forms),

    |<i[h_I,F_mu]>|
       <=2||h_I||(exp(mu R)-1)<F_mu,I>.

Summing with M5 and applying Gronwall yields

    <F_mu(t)> <= C_mu exp(v_mu |t|),
    v_mu=2J_0[exp(mu R)-1].                                   (P3)

On-site internal energies commute with every n_x and contribute zero.
If A is the boundary layer x>=ell-c, c=max(R,2), let
L=ell-c-b. The indicator that any carrier lies there is <=sum_(x in A)n_x.
Thus, uniformly over arbitrary initial matter/receiver entanglement,

    ||Pi_A U_in(t) X_b||
      <= min(1,sqrt(C_mu) exp([v_mu|t|-mu L]/2)).              (P4)

P3–P4 depend on local hard-core capacity and M5, not P_N. No classical
independent-particle or independent-emission assumption was introduced.

## 3. Common-space reference and the complete boundary drive

Record w of length m, retain its first ell labels, with m>=ell and
L=ell-c-b>=0. Let K be the physical subspace retaining that prefix and all
carriers inside it. The actual conditional storage/receiver state is kept.
On K define H_ref as H_in plus all assisted tail geometry/receiver/storage
dynamics reflected at length ell, with all bare geometry terms removed.
Matter terms across the cut are omitted in this **reference only**. On K-perp
use K-perp H K-perp, with no reference cross blocks. This is a self-adjoint
operator on the same full physical Hilbert space as H.

The tail fields are exactly zero by M1. Prefix storage bits are unchanged
spectator blocks; the prefix part of the storage diagonal is a matter-
independent constant in each block. Consequently H_ref on K factors into
prefix matter and tail/receiver dynamics, even for initial correlated fields.
The enlarged root compensators are fixed in each conserved number sector
and introduce no receiver coupling. No new generator reads them globally.

Let a=m-ell. Receiver motion and storage diagonal phases change no word
length. The reference tail's entire changing norm is <=g_max=11/10. Reaching
its lower reflecting boundary first requires a assisted end transitions, so

    ||Pi_(length=ell) U_ref(t) X_m||
          <=min(1,tau_a(g_max|t|)), tau_q(x)=sum_(n>=q)x^n/n!. (P5)

This is a complete Dyson sum for the **tail only**, whose norm is uniform.
The extensive matter evolution commutes with that factor and is not counted
in P5. This also covers a=0, though no growing head buffer then follows.

All omitted-action drives on K are now accounted for:

* An assisted pop breaking the retained prefix acts only at length ell;
  its norm is <=g_max.
* Bare append/pop requires a carrier at the head or parent. With head length
  at least ell it vanishes unless a carrier lies at x>=ell-2; its full norm
  is <=4epsilon, independently of particle number.
* Every cut-crossing matter term has an occupied interval meeting both sides.
  It annihilates its interval vacuum. On K the outside has no carriers, so
  this term vanishes unless a carrier lies at x>=ell-R. Summing its norms
  costs at most R J_0, charging each term to a prefix site near the cut.
  Terms wholly outside vanish on K; terms inside cancel exactly.

Putting B_0=4epsilon+R J_0, P4–P5 give the complete drive estimate

    ||(H-H_ref)U_ref(t)X||
       <=g_max tau_a(g_max|t|)
         +B_0 sqrt(C_mu) exp([v_mu|t|-mu L]/2).                (P6)

The common extensive storage/receiver diagonal cancels. The common extensive
bulk matter Hamiltonian cancels too. Bounds are uniform over all matter and
actual receiver states on the declared initial support.

Exact Duhamel, with **full H on the physical side**, integrates P6:

    ||(U(t)-U_ref(t))X|| <= Delta(t),
    Delta(t)=min(2,
       tau_(a+1)(g_max|t|)
       +(2B_0 sqrt(C_mu)/v_mu) exp(-mu L/2)
                                      [exp(v_mu|t|/2)-1]).     (P7)

This includes every physical storage error, matter excursion, geometry change,
complement visit and coherent return; the reference does not delete transitions
from H. No real-energy resolvent is exponentiated. If initially storage-consistent,
prefix loss and storage inconsistency each have probability <=min(1,Delta^2).
For any bounded receiver observable O_R, the difference from the stated
reflected-tail reference expectation is <=2||O_R||Delta. These are controls
of a specified finite-time backreaction experiment, not absence of backreaction.

## 4. Actual information response

Let S be a bounded local physical matter source preserving the initial support
subspace X (or include S X when defining b), and O a local physical receiver
observable. Their fixed-prefix representatives have disjoint supports at
recorded distance d. For a normalized state in X, pulse exp(-i theta S) and
chi=d<O(t)>_theta/dtheta at zero. The two bilinear terms and P7 imply

    |chi-chi_ref|<=4||S||||O||Delta.

Using P2,

    |chi(t)|<=min(2||S||||O||,
       2||S||||O|| s exp(-mu d)[exp(2kappa_mu|t|)-1]
                                      +4||S||||O||Delta(t)).  (P8)

This controls internal, number, diagonal-algebra and gauge-dressed local
motion responses in the supplied interaction class. It does not use UC1's
few-carrier cut-charge noncancellation argument, which does not extend to
multiple equal-charge carriers. Local generator support supplies the bound.

## 5. Explicit positive-density resource sequence and growing window

Take N=32n,n>=8, an accepted record m>=N/2, b=N/8-1 and ell=N/4. Use the
primary R=1,J_0=2,mu=1,epsilon=1/4. Then

    a>=N/4, L=N/8-1, B_0=3,
    v_mu=4(e-1), kappa_mu=4e.

Occupy the first N/8 positions with A at 4q, B at 4q+2, and vacancies
at 4q+1,4q+3, q=0,...,n-1. Then

    N_A=N_B=N/32, P_N=N/16.                                  (P9)

Both species have a nonzero limiting density. The theorem also covers any
other hard-core state supported in this block, including coherent or strongly
correlated states, without changing its constants. The remaining empty buffer
is a predeclared preparation restriction, not a density-averaged inference.

For every 0<=|t|<=N/4096, elementary 2<e<3 bounds in P7 give

    Delta(t) <= exp(-N/5)+3 exp(-N/17).                        (P10)

For example, the head tail is <=exp(e g_max t-a-1)<=exp(-N/5).
The second prefactor is <3, and its exponent is at most
1/2-63N/1024<=-N/17 for N>=256. For local supports of size <=2 and
d>=N/32, P8 consequently gives, for unit source/readout norms,

    |chi(t)|<=4 exp(-N/40)+4 exp(-N/5)+12 exp(-N/17).           (P11)

The constants and the useful O(N) window are independent of P_N. The global
matter norm is nevertheless Theta(N). Matter filling the entire path without
a buffer is outside this particular full-state comparison theorem. Failure
of that buffer hypothesis is not itself a no-go for local bulk observables.

## 6. Nonzero nearby response in an extensive background

Use P9, replacing A's particle at positions 0,1 by the equal coherent adjacent
pair (all other carriers unchanged), or B's particle at 2,3 by that pair.
The gauge fields in each branch are the exact M1 fields; compensators and
the actual conditional receiver are retained. Both coherent source branches
use the same internal level and background/receiver state. Source S=n_si
has variance 1/4; it and the readout O=n_s,i+1 have unit norm. The first response coefficient
is exactly j_s t; internal bond terms commute with both number observables.

Include the unit on-site clock in a local strength bound J_tot<=3. For
nearest-neighbor terms and single-site O,

    ||[H_in,[H_in,O]]||<=8 J_tot^2||O||.

This follows by summing first supports through O's site and then second
supports through at most two sites. The integral Taylor remainder therefore
gives chi_ref(t)>=j_s t-8J_tot^2 t^2, independent of background size.
At t_0=5*10^-5, N>=256, P7 gives Delta(t_0)<10^-9. Thus

    chi_s(t_0)>=j_s t_0-72t_0^2-4Delta(t_0)>3*10^-7           (P12)

for both species. To check the conservative Delta claim, use C_mu<4,v_mu<8,
L>=31 and exp(-15.5)<10^-6; the boundary integral is <6*10^-10 and the
head term is smaller than exp(-64). The positivity is for this explicit
phase-pulse preparation, not for every state or instant. Fully jammed slots
can block hopping, as expected.

These arguments establish a conditional finite-density regime. They do not
derive that preparation from the original neutral seed, select dimension,
restore a real-energy gap, establish relativistic propagation, or predict
a new empirical quantity. S9 formation review and IG1 proof review stay pending.
