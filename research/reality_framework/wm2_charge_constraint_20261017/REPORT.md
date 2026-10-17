# WM2-UC1: charge motion can be constrained; universal locality still fails

**Primary architecture verdict: negative in the frozen class U.**
Explicit hard constraints enforce a growing word cost for charged
position transfer and preserve viable formation. They still allow a
bounded direct interaction between distant internal matter states.
Gauge invariance and flux support do not select a local Hamiltonian
for the whole physical observable algebra.

All mathematics is Builder-authored, independent review pending.
The charter was frozen at `3077093` before the new diagnostic runs.
WM2 and S9 remain unchanged, and the numerical recovery branch stays closed.

## What was constructed

CONSTRUCTION.md specifies seven-level pair electric registers, root
compensators, two charged species (charges 1,2; internal dimensions 2,3),
Gauss operators G_i and flux-support operators C_ij=(1-B_ij)E_ij^2.
The physical space is their joint kernel. The explicitly given H commutes
with every constraint. Its append/pop primitives carry zero-flux gates;
its charged hops shift the appropriate electric register. No global
P h P implementation or unproved finite energy gap is used.

In matter vacuum all added registers are zero, every new gate reduces
to one, and the formation component is exactly the old law. This inherits
S9's conditional status and resources rather than providing new formation
data or independent validation.

## Old witnesses are excluded, including their single-link dressings

Undressed charged all-pairs and star terms violate Gauss. Dressing one
link repairs Gauss but violates the additional flux-support constraint
off the occupied graph. Adding nonzero local hopping leaves those
commutators nonzero. Thus simply saying "gauge invariant" would have
failed; the extra supplied flux/graph relation does substantive work.

CHARGED_TRANSPORT.md proves that all physical primitive words transferring
a carrier across D recorded positions require at least

    min(ceil((D-1)/r), m-ell+1)

changing operations, including prefix erasure/restoration routes. The
complete Dyson tail and a convergent positive-resolution resolvent
expansion pay for every virtual path and its multiplicity. No real-energy
gap follows. Both species retain a uniform nearby response greater than
9*10^-7 at t=10^-4 under the declared single-species preparation.

## The surviving bypass lies inside the physical sector

The same class permits

    H_int=(1/4) sum_(i!=j) Z_Ai X_Bj.

It changes internal states, not charge positions or flux. It commutes
with **every** constraint, has two-event primitive support, is root-blind
and label-covariant, and has norm 1/4 on the declared two-carrier sector.
The physical law also retains nonzero hoppings j_A=1/64 and j_B=1/96.

Conditional on a formed word and an explicitly prepared separated
charged receiver, let j=floor(N/4)-1 and pulse A's internal state at
the root. INTERNAL_OBSTRUCTION.md derives at fixed time t=pi, for N>=8,

    |d< Pi_prefix Y_B,wj >/dtheta|_(theta=0)>7/12,
    Pr(original prefix retained)>15/16.

The same response bound holds for the unconditioned local Y_B,wj readout.
The source has variance 1/4. Its distance to the receiver grows as N.
The response already counts unsuccessful prefix retention with zero
readout. Every geometry, storage and receiver process is included by
the exact Duhamel comparison; the physical hopping is not removed.

This conditional counterexample does not depend on S9's formation proof.
The far receiver preparation requires an O(N)-link string or sequential
preparation and is explicitly paid as an extra setup; no new original-trial
probability is claimed. The later source choice is local, with initially
product internal states, so the response is not due to a prearranged message.

## Checks and their limits

| Diagnostic | Outcome |
|---|---|
| Bare/dressed old witnesses, charges 1 and 2, full electric space | Expected nonzero Gauss or flux-support commutators; local addition cannot cancel |
| Complete physical-kernel matrices at N=2,3,4 | Dimensions 141, 1,776, 24,780; Hermitian; zero-flux geometry and single-link matter changes checked |
| Seeded formation blocks versus frozen S4 implementation | Maximum evaluated entry difference 1.12*10^-16; zero leakage out of that component |
| Nonuniform detunings, field superposition and full geometry dynamics at N=3, t=pi | Local and weighted internal response about -0.998997; prefix loss about 0.00014236 |
| N=4 tail-growth/return control at t=.7 | Local and weighted response about -0.342879; prefix loss about 3.75*10^-7 |
| Independent +/- source-pulse propagation | Agrees with the propagated tangent response |
| Both charged nearby protocols | Nonzero response consistent with the uniform analytical lower bound |
| Integer cut changes, relabeling/orientation and rational constants | All recorded controls pass |

These use complete **two-mode diagnostic receivers**, not the physical
L=64N^2 receiver for a new formation trajectory. The all-size claims rest
on the displayed analytic arguments. Ordinary double arithmetic is not
an outward-rounded hardware certificate. CHECKS.json and executed-source
hashes preserve the values and scope.

## What this decides

The supplied constraints enforce a finite-range charged transport class;
they do not derive why nature selects those constraints, or shared locality
for the whole matter algebra. The product of two allowed independent local
physical observables again preserves all constraints. The counterexample
shows that this algebraic closure can implement a distant interaction with
uniform norm and an actual local-source response.

Therefore a commutant-only Hamiltonian class with these independent
internal degrees of freedom cannot establish universal geometry. One must
add or derive a restriction on allowed interactions, change the observable
encoding, or provide a dynamical suppression mechanism. No particular such
principle is inferred for nature. This one bounded construction closes at
the negative full gate, preserving its narrower charged-motion result.
