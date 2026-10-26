# Exact internal clocks, records and the missing preparation law

This is a finite supplied history-state construction using standard quantum
mechanics. It makes the architecture's observation rule executable. It is
not a microscopic theory of spacetime or a derivation of quantum outcomes.

## Hilbert space and the one fixed constraint

Let K=C^8_clock tensor C^2_S tensor C^2_R1 tensor C^2_R2. Use the eight gates
U_0 through U_7:

1. H on S.
2. CNOT from S to R1.
3. CNOT from S to R2.
4. exp(−i .37 Z_S).
5. The inverse phase.
6. The inverse second copy.
7. The inverse first copy.
8. The inverse H.

Their ordered product is I. Write W_0=I and W_(t+1)=U_t W_t. The cyclic
propagation constraint is the single time-independent positive operator

\[
C=\frac12\sum_{t=0}^{7}\big[
 (|t\rangle\langle t|+|t+1\rangle\langle t+1|)\otimes I
 -|t+1\rangle\langle t|\otimes U_t
 -|t\rangle\langle t+1|\otimes U_t^\dagger\big],
\]

with t+1 understood modulo eight. The clock reading t is a label of an
internal degree of freedom. No external sequence of applied controls is
used once C and the state are supplied. This circuit is deliberately
engineered; its existence does not show that the universe selects it.
The finite clock cycle is not a claim that cosmological time is periodic.

Define the isometry

\[
V=\frac1{\sqrt8}\sum_t|t\rangle\otimes W_t,
\quad P=VV^\dagger.
\]

Then C V=0 and V^dagger V=I. Dressing by
D=sum_t |t><t| tensor W_t maps C to the eight-cycle Laplacian tensor I_8.
The spectrum is 1−cos(2 pi k/8), repeated eight times. The physical kernel
has dimension eight and the first gap is .292893218813.

For any density operator rho_0 on S,R1,R2, rho_phys=V rho_0 V^dagger is a
stationary physical state. Conditioned on a clock reading, the state is
rho_t=W_t rho_0 W_t^dagger. A stationary whole can therefore contain
changing relational observations. The model's C is not the GR Hamiltonian
constraint and its eigenvalue is not an observed total cosmic energy.

## Physical observation effects

The raw clock projector need not commute with C. For an outcome projector
Pi_a on the matter/record bits define

\[
E_{a|t}=8P(|t\rangle\langle t|\otimes\Pi_a)P.
\]

Restricted to the physical subspace these are positive effects with
sum_a E_(a|t)=P and [C,E_(a|t)]=0. Their probabilities are

\[
\operatorname{Tr}(\rho_{\rm phys}E_{a|t})
 =\operatorname{Tr}(\rho_0W_t^\dagger\Pi_a W_t).
\]

This gives a physical-subspace dictionary rather than merely multiplying
unphysical clock projectors. The Born rule and the meaning of the pointer
basis are supplied. No new probability law is inferred from stationarity.
An autonomous reader implementing every compressed effect has not been
constructed; the declared observable algebra is not itself that apparatus.

## Records and retained coherence

For rho_0=|000><000|, clock reading 3 has
(|000>+|111>)/sqrt(2) on S,R1,R2. Each record agrees with S, and the two
records share one bit of classical pointer information. The reduced S state
has entropy one bit. At readings 3–5 the correlations are retained; the
phase changes the global XXX witness but not the pointer probabilities.

The coherent state gives <X_S X_R1 X_R2>=1 at reading 3, whereas its
diagonally dephased mixture gives zero. Both have identical individual
pointer probabilities. Thus copying and local decoherence do not select a
unique global outcome. At reading 7 the copies have been undone, the
records are blank and S is pure again. A global arrow of time is not derived
from these temporary records. The low-correlation input and clock orientation
were chosen. This model proves no macroscopic record lifetime or human
consciousness mechanism.

## Closure leaves the state free

For any p in [0,1], choose
psi_0(p)=H[sqrt(1−p)|0>+sqrt(p)|1>] tensor |00>. Every resulting history
is pure, obeys the same constraint and starts with blank records. At reading
3 it predicts record probability p. The numerical examples p=0,.1,.5,.9,1
illustrate an analytic whole interval, not a universal claim inferred from
five grid points.

The maximum-entropy state on the physical subspace is P/8. Its conditional
state is W_t(I/8)W_t^dagger=I/8 for every t. It generates no changing pointer
correlations. The proof holds for any unitary sequence in this isometry
class, not all possible quantum-gravity descriptions. Maximum ignorance
does not supply the low-entropy preparation required here. Global purity
also leaves the whole p-family available.

For comparison, the additional state condition Tr(rho_phys F_0)=0 with

\[
F_0=V(I-|000\rangle\langle000|)V^\dagger
\]

selects the blank-input example. F_0 is positive and commutes with C; on
the original physical space its kernel has dimension one. The p=.1 history
has expectation .2 and is excluded. This is a genuine strict restriction
of allowed preparations at the same observation interface. It is supplied
by a state law, not forced by closure. The physical Hilbert space remains
dimension eight and its observation operators are unchanged. A selected
actual state is different from a one-dimensional fundamental state space.

Page–Wootters and history-state constructions are prior art. This packet
combines an elementary finite cyclic example with the specified physical
POVM and state-selection controls. No general observer dictionary for
quantum gravity is claimed.
