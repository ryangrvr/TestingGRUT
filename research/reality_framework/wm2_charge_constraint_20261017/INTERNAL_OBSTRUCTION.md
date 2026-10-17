# The constraint survives; universal matter locality does not

**Primary gate: failed in the frozen full operator class U.** The charged
position theorem survives at its scoped, author-derived status. The
counterexample below is itself author-derived pending independent review.

## A physical bounded term requiring no flux transfer

Take C3 from CONSTRUCTION.md:

    H_int=kappa sum_(i!=j) Z_Ai X_Bj, kappa=1/4.                 (I1)

Each summand acts on two event/internal registers and zero electric
registers. It changes neither species' position, number, charge or
compensator. Thus it commutes with every G_i and C_ij, including on the
ambient electric space. It also conserves the original storage charge.
It is root-blind and covariant under all event relabelings.

In the one-A, one-B sector each distinct position pair is an orthogonal
block, carrying kappa Z_A tensor X_B. Its norm is exactly kappa,
independent of N. No large all-pairs norm or small-N artifact is used.
Giving the constraints a larger energy penalty would not exclude I1:
it already acts entirely within their exact kernel.

## An actual gauge-invariant source-to-receiver experiment

Fix a formed word w of m>=ceil(N/2), N>=8. Set ell=floor(N/4),
j=ell-1. Prepare A at the root and B at w_j, with their compensators
at the root, and the unique physical flux string. The original ell-label
prefix is recorded before the signal experiment. Keep an arbitrary
admissible conditional storage/receiver state; do not reset the field.

A begins in internal state |+y>, B in level |0>. Pulse the local,
gauge-invariant source

    exp(-i theta X_Ar/2), Var(X_Ar/2)=1/4.                       (I2)

Measure the success-weighted observable

    O=Pi_(original ell-prefix retained) Y_B,wj,                 (I3)

where Y acts on B levels 0,1 and is zero if B is absent there or in
level 2. It has norm <=1 and commutes with the constraints. The response
chi=d< O(t)>/dtheta at zero is a local internal-source response, not a
charged arrival probability or a phase derivative of number density.
The unconditioned **local** receiver Y_B,wj obeys the same lower bound
below. Prefix weighting is an additional check, not the source of a
nonlocal postselection signal.

For a mathematical reference only, set the two hopping coefficients to
zero while retaining every receiver, storage and geometry transition
in H_0=H_UC1-H_move,A-H_move,B. Actual H_UC1 retains both nonzero hoppings.
Under H_0 the particle positions are fixed. B's charge-2 flux string
prevents a pop through w_j, so the prefix is retained exactly. Geometry
and receiver evolution may continue beyond it. They do not act on the
internal levels and commute with I1 in this fixed-position component.
The third B level has a separate on-event energy and is initially empty.

The internal two-level calculation is exact:

    <Y_B>(t;theta)=-sin(theta) sin(2kappa t),
    chi_0(t)=-sin(2kappa t).                                    (I4)

This gives a constant response across recorded distance j proportional
to N, with no changed bond, transferred gauge charge or virtual activation.

## Both species still move in the physical law

Let h=||H_move,A+H_move,B||<=2j_A+2j_B=5/96. Exact Duhamel, with all
geometry and receiver states present, gives ||U-U_0||<=h|t|. For the
source S=X_Ar/2, ||S||<=1/2, and observable I3,

    |chi(t)-chi_0(t)|<=2h|t|.                                  (I5)

For example, write chi=i< [S,U^* O U] > and use
||U^* O U-U_0^* O U_0||<=2h|t|. This controls the whole response,
not just leakage at the final time. The actual prefix-loss probability
is at most (h|t|)^2 because U_0 never loses that prefix.

At t*=pi/(4kappa)=pi, independent of N, pi<4 implies

    |chi(t*)|>1-5/12=7/12,
    Pr(prefix retained at t*)>1-25/576>15/16.                    (I6)

I6 holds for both Y_B,wj and its prefix-weighted version: both have norm
at most one and coincide under H_0. There is no renormalization on rare
surviving records. Both hopping
coefficients are fixed, nonzero, independent of N and different between
species. This is not a frozen-geometry or zero-mobility physical model.
Arbitrary configuration-preserving diagonal electric energies independent
of the internal states also leave the argument intact: H_0 still factors
and the same hopping norm bounds the difference.

## Preparation cost and the formation dependency

I6 is conditional on a formed word **and a separated charged receiver**.
It does not inherit S9's root-pair injection protocol or its original-trial
probability. Preparing B at w_j with its root compensator requires the
charge-2 string along j links. This is an explicit O(j) support/preparation
cost if done by a string injection, or requires sequential gauge-respecting
transport/preparation. A controlled protocol under the always-on full law
has not been established here. It is not a free remote injection or a
supplied probability that the experiment reaches that state autonomously.

Once the receiver is prepared, the source pulse and readout are genuinely
local, preserve all constraints and have the response I6. A proposed
state-uniform locality theorem for the full class must cover such physical
receivers, as ordinary communication tests do. No S9 formation estimate is
needed for this conditional counterexample. No accepted-original-trial
lower bound is claimed without paying this extra preparation stage.
The preparation is the same for every later choice of source pulse, and
the initial internal state is a product. Its long setup does not itself
explain the fixed-time receiver dependence on a new local source choice.

## Why this is a structural obstruction to commutant-only selection

**Product-closure lemma.** If physical local observables A_i and B_j
commute with every constraint (or each preserve P and its orthogonal
complement), so does their product. Indeed [G,A_i B_j]=[G,A_i]B_j+
A_i[G,B_j]=0, for each constraint; the same holds for C and number.
If their bounded register supports have bounded union, their product
remains a finite-support operator regardless of recorded distance.

Therefore preserving constraints alone cannot exclude all distant
interactions while admitting independent gauge-invariant internal
observable algebras. I1 gives an explicit root-blind, norm-uniform
covariant witness, avoiding an unbounded extensive sum. I2-I6 exhibit
its operational consequence with noncommuting source/readout observables.

This is not a no-go for every microscopic theory or every possible
constraint. It is an obstruction to the **commutant-only Hamiltonian
class with independent local internal degrees of freedom** used here.
One must additionally restrict the Hamiltonian's permitted products,
change the physical observable encoding, or derive a dynamical suppression
mechanism. Imposing occupied-bond or connected-support interactions on
the internal algebra would be an additional locality input. No such
restriction has been derived by this construction, and it is not added
after the negative result.
