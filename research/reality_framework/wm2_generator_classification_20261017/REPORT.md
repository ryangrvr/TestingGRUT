# WM2-IG1: physical observables do not select the interaction law

The bounded task is complete at **author-derived, independent-review-pending**
status. The classification supports a firm negative and a conditional positive:

1. UC1's present constraints, independent local internal degrees of freedom,
   covariance and bounded few-carrier norm do **not** force information locality.
2. A separately supplied universal relational generator rule can bound both
   charged and internal-state responses under complete evolution in a recorded,
   buffered few-carrier interior. That rule is an input, not an emergence result.

No new conservation patch, architecture campaign or formation solver was opened.
UC1's negative full gate and S9's pending independent proof status are preserved.

## The general obstruction

The constraint commutant and the physical-kernel observable algebra are not
the same object; neither specifies which observable products are elementary
Hamiltonian interactions. Independent physical qubits at separated events admit
the product kappa Z_A X_B whenever all physical two-event products are allowed.
The exact pulse response is -sin(2kappa t), independent of separation. UC1's
occupation-tagged covariant sum has total norm kappa even as N grows.

A **diagonal** Z_A Z_B product also signals. Thus restricting where particles
move, or counting only registers whose basis value changes, is insufficient.
The necessity test derived here is

    Gamma_H(i,j)=sup ||[A_i,[H,B_j]]||
      <= C v exp(-mu d(i,j))

whenever a state-uniform exponential commutator front with differentiable
zero-time prefactor is claimed. The present conservation constraints impose
no such bound. The theorem says the broad class **contains** a witness;
it does not say that every permitted Hamiltonian is nonlocal.

## A sufficient restriction, with its price visible

RI(R,w) requires every matter generator, including internal and diagonal
interactions, to have bounded support on a consecutive occupied interval.
It also controls coefficients, receiver contacts and the total few-carrier
changing norm. This is stronger than finite event support and is explicitly
**supplied universal relational locality of the interaction law**.

For a recorded length-m word, carriers initially at positions <=b, and a
retained prefix of length ell, define

    a=ceil(d/R),
    k=min(m-ell,ceil((ell-max(R,2)-b)/R)).

For k>=1, source S preserving the setup subspace and a local receiver O,
the full physical response obeys

    |chi(t)| <= 2||S||||O|| tau_a(2w|t|)
                 +4||S||||O|| Delta_k(t),
    Delta_k(t)=min(2,(5/2)tau_(k+1)(2|t|)),
    tau_q(x)=sum_(n>=q) x^n/n!.

This controls information, not just charged arrival. With distance and buffers
proportional to N it gives an O(N) observation window with exponentially
small corrections. A common-Hilbert-space reference is used only in the proof;
the physical Hamiltonian retains every permitted geometry/receiver/storage
process and every return. No real-energy gap is invoked or restored.

The class is nonempty: the original two different charged hoppings, B's
third-level clock and a nontrivial occupied-bond Z_A X_B interaction satisfy
the norm cap. Both nearby single-species responses remain nonzero. In matter
vacuum the law is exactly the previous formation law, so S9 is inherited only
with its existing proof status and resource cost.

The guarantee is conditional on a recorded interior and declared carrier setup.
It does not cover arbitrary density, unrecorded coherent geometry, arbitrary
small nonlocal perturbations, or free separated charged-receiver preparation.

## Finite verification and its limits

The final suite passed in 1.59 s with peak resident memory 89.97 MiB. It uses
ordinary double precision and two-mode diagnostic receivers, **not** a new
simulation of S9's physical receiver or its formation regime.

* Independent two-qubit controls reproduce both signaling identities and
  check a deliberately commuting readout boundary.
* Full physical-kernel matrices through N=4 preserve Hermiticity and the
  exact old matter-vacuum block. The two nearby responses exceed their
  analytic lower bounds 1.50469e-6 and 9.83860e-7 at t=10^-4.
* A common-space N=4 reference has full-state error 3.06141e-5, below
  the evaluated bound 0.00202694. Its nonzero omitted-action amplitude
  exercises the shortest return path; the response bound in that particular
  co-located-carrier control is not a positive propagation measurement.
* History-retaining bounded-depth support enumeration checks first-action
  counts; independent seven/eight-site prefix matrices check low-order
  commutators. Adding a forbidden distant product produces an evolved
  commutator about 0.015, well above the local bound. Tiny positive-control
  commutator values below the floating tolerance are not numerical certificates.

These are diagnostics of implementations and finite identities. They do not
independently establish the all-size mathematics. The initial diagnostic
setup correction is disclosed in EXECUTION_NOTES.md; no budget was relaxed.

## Decision

Stop the constraint-only repair sequence. For controlled candidate building,
universal relational generator locality can be used only as an explicit
structural input until a deeper selector is derived. Changing the physical
observable encoding is another possible route, outside this bounded task.

The load-bearing proofs are in OBSTRUCTION.md and SUFFICIENT_RULE.md; the
assumption price is in ASSUMPTIONS.md. There is no new derived universality,
selected dimension, gravity result or empirical prediction.
