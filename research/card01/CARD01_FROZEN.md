# Card 1 — single-operator phase/decay law

**DRAFT — author science, pending independent review.** This is a genuine candidate
attempt charged to the three-card budget, not a standard-control exemption.
The law and scope below are frozen before computational hostile evaluation.
The law will not be repaired after a failure.

## C1. Identity, provenance and search disclosure

Slot 1. Authorization is recorded verbatim in `AUTHORIZATION.md`. Base
`9b07210445b17fa3cef99324e92cd2df9fd09ee3`; governing charter blob
`8dd7b05ab1ccec191e1ee8c54fcb41fbfd85c052`, verified live against
`ryangrvr/vsctestinggrut:grut2-stage3`. Evaluator envelope
`c32f672288fcad2a9d08772b8ef2bf46ddee70a2`, subject source
`3a4bd12d9b269a3a74adaf994db50121b8e51838`, not externally cleared.
Freeze hash/commit will be in the separate freeze receipt and attempt ledger.

This is the first law chosen in this branch. No program searched candidate laws.
The author has read the kinetic-freedom, fixed-lift, admissibility and selection
audits and the three standard positive controls. They supplied the demand for a
joint restriction; none is evidence for this law. The sole design choice is to
try linking the coherent generator and dissipative generator to one operator.
A preliminary literature lookup already returned Milburn (1991) and Adler (2001);
they are disclosed now, not treated as an independent discovery. The exact
assembled equivalence will be audited. No lock data are accessed. No prior
GRUT law card is claimed in the available record. Global draft-count completeness
is not independently certified. Conservatively log two local IP-12 entries:
the symbolic family and the one numerical instance; local selection-tax lower
bound is log2(3). Further time points, input states and basis controls do not
change this same parameter assignment. Comparator models are tests, not repairs.

## C2. Complete mathematical law

POSTULATED commitments C: a finite complex Hilbert space H of dimension d >= 2;
an external real time t >= 0 with its unit; a constant tau > 0 of dimension time;
the usual positive trace-one density operators. No Gibbs state, thermostat,
spatial observer, privileged tensor factorization or initial memory is supplied.
Quantum structure and the clock are supplied, not derived.

Xi = (A, rho), where A is a finite Hermitian operator of dimension inverse time
and rho:[0,infinity)->End(H) is continuously differentiable. Initial conditions
are arbitrary A(0)=A0=A0^dagger and rho(0)=rho0>=0, Tr rho0=1. No boundary
condition at infinity. All A0, including scalar and degenerate operators, and
all initial density operators are admitted. Tau is shared within an instance;
the proposal makes no derived numerical prediction for its absolute value.

For every t >= 0:

    K[C,Xi] = (dA/dt, d rho/dt + i[A,rho] + tau/2 [A,[A,rho]]) = (0,0).

[M,N] := MN-NM. Multiplication is operator composition; dagger is the Hilbert
adjoint; Tr is the ordinary finite trace. This is a bounded linear ODE in rho
with constant A. Reference implementation: `law.py`.
No operator named "self-consistency" is left undefined.

The physical energy interpretation H_phys=hbar A supplies hbar in addition;
it is not part of the generating derivation. Eigenvalues of A are angular
frequencies. A shift A->A+cI has no effect on the law. Basis changes conjugate
all objects. The continuum scope is finite d, not QFT or gravity. No size-,
sector- or protocol-dependent sub-law is allowed.

The mathematical equation is complete; a full audited L0 unfolding of its
continuum imports and proposed chain is NOT supplied by merely writing it this
way. That is an explicit certification obligation, not a macro-price claim.
L0 matrix arithmetic will be itemized in the price report. If the chain cannot
be completed without protected readout or partition inputs, the card dies.

## C3–C7. Ontology, decidability, inventory, price and extraction

Finite matrices, their adjoint/product, continuum time and a given Hilbert
inner product are postulates. Exact rational fixtures and spectral formulas
give controlled computations. Degeneracy is algebraically decidable on exact
algebraic input; no uniform decision procedure for arbitrary unencoded real
input is claimed. A generic real-domain DC certificate is therefore pending.

Supplied vocabulary to price includes VB-2, VB-3, VB-7 when records are computed,
VB-8, VB-9, VB-10, VB-11 when energy is used. Tensor factors in composition tests
are conventional comparator resources, not a generated subsystem claim.
No quantum selection credit is claimed from supplied Hilbert/Born structure.
A0, rho0, d, tau, units, preparations, POVMs, initial apparatus states and any
proposed record/interface class are supplied data. Their entire price, including
tables and protected items, must be disclosed; menu price is not their full price.

Proposed X_Pi: spectral blocks of A (equivalently candidate fixed-algebra blocks,
if that equivalence is proved). Proposed Z: their weights Tr(P_j rho).
This construction is tested for input encoding and law dependence. It does
not assume these blocks qualify as generated physical subsystem identities.
Persistent populations are a proposed memory, not memory-writing dynamics.

There is no claimed unique detector map. Candidate symmetry stabilizers may
exist but must NOT be identified with T_Pi without a physical accessible-record
derivation. If X_h or X_T is unavailable, return bottom. The operational test
uses explicitly SUPPLIED Ramsey preparation/POVMs. Gamma is then the finite
protocol-dependent Born record family. Epsilon_R is NOT specified separately,
and will not be computed without a certified carrier and admissible T_Pi.
No numeric epsilon is fabricated to fill a broken chain.

## C8–C13. Restriction, baseline, locks and measurement

Let a_j be distinct eigenvalues identified by the spectral blocks. Define
omega_ij=a_j-a_i independently from the phase of a Ramsey coherence and Gamma_ij
from its amplitude decay. The primary proposed restriction, for distinct
nonzero-frequency pairs (i,j),(k,l), is

    R = Gamma_ij omega_kl^2 - Gamma_kl omega_ij^2 = 0.

At least three distinct levels are needed to give a nontrivial multi-pair test;
the rest of the law's domain is not silently removed. Claimed consequence to
prove: Gamma_ij=tau omega_ij^2/2 for all pairs. It is a dynamical restriction,
not an epsilon definition. Whether it genuinely couples GENERATED identities
and accessible interfaces is a separate mandatory gate.

Hostile baseline at fixed A and the identical Ramsey carrier: ordinary CPTP,
trace-preserving, unital, time-homogeneous pure-dephasing quantum dynamics

    L_std rho = -i[A,rho] - 1/2 sum_mu [D_mu,[D_mu,rho]],

with arbitrary commuting Hermitian diagonal D_mu. Preserve the same populations,
energy, preparations, POVMs, time unit and record statistics definitions.
No Gibbs preparation is needed. Check weighted detailed-balance compatibility
where applicable; do not claim a full microscopic KMS/FDT result from a finite
semigroup alone. Compare with exact finite-ancilla implementations on a finite
intervention grid. No relativistic/QFT/EFT claim is made for this finite model.

One declared numeric instance only: d=3, A=diag(0,1,3) in a supplied inverse-time
unit, tau=1/5 in its time unit. The comparator has one D=sqrt(1/5)diag(0,2,1).
Matrix-unit, conservation, semigroup, positivity and basis-change identities
are controls. No coefficient is adjusted after observing the result.

Operational test: prepare each (|i>+|j>)/sqrt(2); record the complementary
Ramsey observables X_ij=|i><j|+|j><i| and Y_ij=-i|i><j|+i|j><i| at known times.
Extract z=<X>+i<Y>, its phase frequency and amplitude decay. Calibrate tau on
one pair and predict the others without further adjustment. Resolved unequal
ratios 2Gamma_ij/omega_ij^2 refute the stated intrinsic-channel restriction.
No named platform, residual-noise separation, independently fixed tau,
power/sample-size certificate or sealed lock is established. A programmable
simulator may test arithmetic but cannot establish a universal intrinsic law.

## C14. Frozen kill conditions

KILL if any of the following holds; no post-test repair is authorized:

1. The full chain requires a supplied partition, readout/interface class or
   observer/apparatus structure; stable populations alone do not discharge it.
2. An input-only decoder recovers the proposed identity/interface, or the same
   proposed X_Pi/X_T is returned when the law is ablated.
3. Nontrivial generated identity fails anywhere in the stated full domain.
4. The relation is a pre-existing consequence of the assembled single-operator
   law, rather than a distinctive GRUT result; acknowledge any exact equivalence.
5. Standard physics with the same operation/resource assumptions reproduces the
   entire proposed evidence and no intrinsic/ordinary-noise discriminator is found.
   This does not logically negate an extra universal postulate; it kills the
   claim that the displayed reduced channel establishes it.
6. The restriction is only a fitted numerical coincidence or fails an exact
   identity; or a cheap dissipation clause forces it while the identity chain
   is independent of that clause (joint-necessity failure).

All universal charter kills remain acknowledged. This is an author investigation,
not a substitute for the independently administered charter batteries.

## C15–C20. Remaining mandatory fields and non-claims

Selectivity: supplied quantum supports; zero earned selectivity. AI instances:
the exact qutrit above; scalar A; a repeated-eigenvalue composite; deletion of
the double-commutator clause; an arbitrary spectral-basis conjugation.
Representation: simultaneous unitary conjugation, spectral relabeling,
energy-zero shifts; record calibration moves are not generated by these.
Tradeoff: no presumed monotone identity/response rule. FR/GY tables: endpoint
requirements unfulfilled, no pending theorem borrowed to claim a pass.
No cross-sector transfer, quantum emergence, consciousness mechanism,
thermodynamic emergence, gravitational limit, universal constant determination,
confirmed prediction, full price/compression certificate or external sign-off
is claimed. These are gaps to expose, not promises that survive a kill.

The law is frozen as written. The external-prerequisite waiver permits this
attempt; it neither approves the evaluator nor supplies certificates.
