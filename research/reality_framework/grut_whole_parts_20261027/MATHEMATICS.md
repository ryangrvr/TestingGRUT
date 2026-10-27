# Exact whole/parts construction

All spaces and controls below are supplied finite quantum mechanics. They are not a
quantization of the continuum GRUT action. Charge is in arbitrary logical units.

## 1. Kinematics, constraints and one chosen state

Take K = C⁴_T ⊗ C²_X ⊗ C²_C ⊗ C²_R ⊗ C²_P. On the nonclock factors define

\[
Q_X=6Z_X,\qquad Q_C=6Z_C,\qquad G=Q_X+Q_C.
\]

G is **global neutrality only**. There are no local Gauss operators or gauge links.
Let U=diag(1,i) act on P. Define W_t = U_P^t for t=0,1 and
W_t = CNOT_(X→R) U_P^t for t=2,3. Tensor identities are implicit. Put

\[
S=\sum_{t=0}^{3}|t+1\rangle\langle t|\otimes W_{t+1}W_t^\dagger,
\qquad C=2I-S-S^\dagger,
\]

with cyclic clock indices and W_4 = W_0. S is unitary, S⁴=I, C is positive, and
[C,G]=0. Its spectrum is 0,2,2,4 with each clock eigenvalue repeated 16 times.
The constraint C=0 selects history states. Global neutrality then leaves dimension 8.

The isometry Vψ = (1/2)Σ_t |t⟩W_tψ maps the nonclock space into ker C. On ker G,

\[
\mathcal H_{\rm phys}\simeq
\mathbb C^2_{XC,\,\rm neutral}\otimes\mathbb C^2_R\otimes\mathbb C^2_P.
\]

The neutral pair's basis is |0_X1_C⟩, |1_X0_C⟩. The 64×8 logical isometry J is
stored in HISTORY_MODEL.npz; J†J=I and JJ† is the physical projector.

Choose the one actual history

\[
|\Psi\rangle=\frac12\sum_t|t\rangle W_t|\psi_0\rangle,
\quad |\psi_0\rangle=\frac{|010\rangle_{XCR}+|100\rangle_{XCR}}{\sqrt2}
\otimes|+\rangle_P.
\]

It obeys C|Ψ⟩=G|Ψ⟩=0 and has purity one. The physical space still has dimension
eight. The constraint's zero eigenvalue is not a measured zero sum of matter and
gravitational energy, nor a thermodynamic ground-state assertion about every subsystem.

## 2. Physical internal observations

For a nonclock effect O commuting with G and clock reading t define

\[
D_{O|t}=\sum_{n=0}^{3}S^n(|t\rangle\langle t|\otimes O)S^{-n}.
\]

It commutes with C and G. Positivity and Σ_a O_a=I imply positivity and
Σ_a D_(O_a|t)=I. The physical Born rule is

\[
\langle\Psi|D_{O|t}|\Psi\rangle
=\langle\psi_t|O|\psi_t\rangle,
\qquad |\psi_t\rangle=W_t|\psi_0\rangle.
\]

The four equal clock weights cancel in conditional probabilities. Raw clock-conditioned
kinematic projectors are not substituted for physical effects. Probe X,Y,Z dress to
a full Pauli algebra; their physical compression gives the third logical factor. X's
own charge is physical; its charge-flipping Pauli operator alone fails neutrality.
The compensated transfer |01⟩⟨10|+h.c. is neutral and acts on the pair.

At t=0,1 the record is blank. At t=2,3 it records X's charge sign perfectly, with
probabilities 1/2 each and conditional charges ±6. The probe has Bloch vectors
(1,0,0), (0,1,0), (−1,0,0), (0,−1,0). These are internal correlations of one stationary
history. The cycle's wrap reverses record copying. No irreversible arrow is implied.

At t=0 the relational dressing of the compensated transfer has expectation one in Ψ
and zero in the dephased mixture of the two charged histories. This distinguishes
global quantum coherence from a classical mixture with the same local charge statistics.
The joint operator is permitted by global neutrality. Local operational access to a
corresponding dressed operator in a continuum gauge theory would require its own proof.

An alternative base |010⟩|+⟩ has X charge +6 and C charge −6, still satisfying both
constraints. A blockwise-dressed X_X X_C exchanges its charge-conjugate histories.
This is a symmetry of this finite charge model; it is not asserted to be relativistic CP
or the time-reversal-type conjugation of the earlier no-antiparticle toy.

## 3. Logical clock resolution

Replace perfect clock readout by the supplied POVM

\[
\widetilde E_t=(1-\eta)|t\rangle\langle t|+\eta I_T/4,
\qquad 0\le\eta\le1.
\]

For the probe, ρ_(P|t)=(1−η)|p_t⟩⟨p_t|+ηI/2, so

\[
\operatorname{Tr}\rho_{P|t}^{\,2}=1-\eta+\eta^2/2.
\]

The globally pure state and constraint have not changed. The ideal logical readout
η=0 gives exact conditional unitarity. This refutes an implication from purity and
no external clock alone to a universally positive decoherence floor. It does not
establish arbitrary clock precision in a gravitating laboratory. Physical clock
resource/noise bounds require a clock Hamiltonian and operational assumptions.

## 4. Dimension and information

A unital *-representation of the qubit algebra on a one-dimensional ordinary Hilbert
space would send Pauli X,Z to scalars x,z. X²=Z²=I gives x,z=±1; XZ+ZX=0 would give
2xz=0, impossible. Thus a faithful ordinary internal qubit requires dimension at least
two. k fully independent ordinary logical bits require capacity at least 2^k. This
model supplies three, while choosing one actual global state.

This does not refute all one-dimensional fundamental-Hilbert-space proposals: an
effective observer dictionary with different assumptions must be specified and checked.
Nor does spatial compactness make a continuum field theory finite-dimensional.

In finite dimension, Tr[A,B]=0 forbids [q,p]=iI exactly. The explicit four-state oscillator
has [a,a†]=diag(1,1,1,−3), trace zero. A low-state effective oscillator can approximate
the canonical relation; a finite exact fundamental algebra must differ somewhere.

## 5. What makes an individual physical in the field theory

The example establishes subsystem observables and records. A material object's
individuality additionally needs identifiable local properties, a dynamically sustained
configuration, and a lifetime longer than its relevant internal times. Approximate
autonomy can coexist with energy exchange and correlations. Those properties have to
come from the action and state; a convenient tensor factorization is insufficient.

For Candidate 1 the relevant tests are E(Q), perturbation stability, exchanges with the
surrounding fields, and actual formation/survival. The finite factors here do not derive
the X/phi soliton, replace the core/wave bookkeeping, or certify a relic abundance.
