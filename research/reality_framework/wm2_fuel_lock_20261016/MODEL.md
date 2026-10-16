# WM2-S2: explicit storage-assisted assembly

S2 is a new law, referenced to S1 at
`7b1d9234deaeaddea3f34c419bdb35b99efd944b`. It does not modify that
packet or frozen CP-12. The main construction retains the single-end
grammar. Branching is a separately declared sensitivity proposal.

## Hamiltonian and primitive support

Retain S1's bond bits b_uv, head bits h_i, free bits f_i and neutral
probe c_i. Let X_uv=b_uv^+ h_v^+ h_u^- f_v^- append a free v at the
current head u. Add one two-level **energy storage** bit q_i per event.
Define C_i=f_i+q_i-1, M=sum_i C_i^2 and F_uv=n_cu+n_cv. The complete
time-independent candidate is

\[
\begin{split}
H_\delta={}&U\sum_i f_i+(U+\delta)\sum_i q_i
+\Lambda\sum_i C_i^2\\
&-\frac g{\sqrt N}\sum_{u\ne v}
 (X_{uv}q_v^++q_v^-X_{uv}^\dagger)\\
&-\frac\epsilon{\sqrt N}\sum_{u\ne v}
 F_{uv}(X_{uv}+X_{uv}^\dagger)\\
&-J\sum_{u\ne v}b_{uv}(c_u^\dagger c_v+c_v^\dagger c_u).
\end{split}
\]

Here b in the last line denotes its occupation. In the one-probe sector
F_uv is 0 or 1 and commutes with the geometry operators. A storage-assisted
term acts on five qubits. Each summand of the probe-assisted mismatch
interaction acts on five registers. Matter hopping uses three. C_i^2
uses two on-event bits. All ordered pairs have identical coefficients;
no coordinates, fixed pair affinity or spatial adjacency is supplied.

S2 **replaces** the independent bare S1 append with the storage-assisted
append. Unsynchronized append/pop terms remain as an explicit probe
interaction. Their vanishing in probe vacuum is a property of this
declared law, not a numerical deletion of transitions. Once a probe is
injected every reachable storage-mismatch state is included.

The construction conserves head number, probe number and K plus bond
number, where K=sum f_i. Storage-assisted terms conserve each C_i;
the probe-assisted terms can change them. Reverse terms are always
present. The bits called spent or charged are **higher-energy storage
bits**, not an irreversible waste sink.

## Initial registers and energy ledger

Use the original CP-12 neutral ray tensor the following uniform-root
seed and probe vacuum. For each root r: no bonds; h_r=1; f_r=0,q_r=1;
all other h=0, f=1,q=0. Superpose these N root choices with equal phase.
The root and its initially charged storage bit are supplied seed resources.

The benchmark embedding is the full product Hamiltonian
H_CP12 tensor I + I tensor H_delta. Original CP-12 modes are not removed,
and their gauge law is unchanged. No energy flows from that factor into
S2; the new spins provide the explicit binding-energy storage channel.
This product embedding does not establish a common metric for CP-12's
original charged a,b excitations or a global spectral band for the
unbounded combined system.

The auxiliary initial energy is UN+delta: free-event energy U(N-1),
root-storage energy U+delta, zero mismatch penalty and zero mean kinetic
energy. Code subtracts the constant UN+delta only as an energy origin.

| Primitive process | Free-event energy change | Storage energy change | Penalty change at a locally consistent event |
|---|---:|---:|---:|
| Assisted append | -U | U+delta | 0 |
| Assisted pop | U | -(U+delta) | 0 |
| Unsynchronized append | -U | 0 | Lambda |
| Unsynchronized pop | U | 0 | Lambda |

Thus at delta=0 binding energy is transferred into q exactly. There is
no added source of energy, bath, Markov approximation, external driving
or ignored reverse absorption. At detuning the kinetic sector pays the
remaining delta per net attachment. All mismatch signs are included in
the code; a pre-existing mismatch can release its penalty when removed.

## Full configuration component and protected sector

Geometry remains a distinct-label rooted word w. Root q_r=1 is invariant.
With one probe, the full reachable component consists of every partial
word, every other storage bit pattern, and every probe position, including
free events. Complete finite graphs verify this characterization.
Independent bit-register rules check the tuple representation.

At delta=0, the bare protected sector is P0=1_(M=0), so q_i equals
membership of i in w. It includes **all assembly stages**, not only
full geometries. This is the crucial change from S1's K=0 sector.
Free-event count K and storage mismatch M must never be conflated.
The isolated exact low spectral band Pi is derived for the auxiliary
one-probe sector in MATHEMATICS.md, not assumed.

In probe vacuum the mismatch interaction vanishes. The consistent seed
therefore stays exactly in P0 and admits an exact N-state radial
description. This reduction does not discard a physically reachable
mismatch. In the later one-probe experiment, all such states return.

For any recorded word s and its attached events i,j, define operational
distance using **all protected primitive configuration moves**, allowing
assembly and return, before computing a response. The prefix lemma
shows that its minimum is still separation along s. This is stronger
than defining distance by freezing the geometry. The choice of a
one-dimensional class remains a supplied grammar input.

## Parameters, resonance and locking limits

Primary values are U=16,Lambda=64,g=1,epsilon=1/4,J=1/4,delta=0.
Units set hbar=1; these are supplied model energy units, not measured
fundamental constants.
Every value, the 1/sqrt(N) rule, the new species, storage bits and seed
are inputs. The local, permutation-symmetric perturbation delta sum q_i
is allowed by the declared symmetries. Exact resonance is not derived.

Spent storage can be converted back into free-event energy by an assisted
pop. Thus the model has no one-way fuel lock. Conditional bulk stability
comes instead from end-only access and the explicit normalization. It
must not be described as an irreversible transition into one fixed graph.

## Separate branching proposal

The branching control replaces the single head by a local outgoing-child
counter d_i in {0,1,2}. A finite-support term may attach a free v to any
bound u with d_u<2, raise d_u, add b_uv, and charge q_v. Its adjoint can
remove v only when v's local child counter is zero. The zero-child test
uses a local qudit, not an implicit whole-graph oracle. This generates
rooted two-child trees without a fixed event adjacency.

It is a **different law**, not included in single-head S2 matrices.
Explicitly, let d_i^+ be the truncated shift on the three-state child
counter and P_(d_v=0) its zero-child projector. One attachment term is

`Y_uv = b_uv^+ (1-f_u) f_v^- d_u^+ P_(d_v=0) q_v^+`.

Use H_branch=U sum f+U sum q+Lambda M
-(g/sqrt(N)) sum(Y+Y^dagger)-J sum b(c^dagger c+h.c.). This control has
no unsynchronized probe interaction. Every accessible state in its
consistent initial component is retained; the diagonal after the UN
shift is zero there. Its raw tree is supplied only as a counterexample
record, never as a successful formation result.

Leaf deletion, attachment elsewhere and restoration permit a five-step
matter transfer across arbitrarily large old tree distance while staying
storage-consistent. The proper protected operational metric must include
that route. This exposes a dimension/locality obstruction for the naive
branching class; it does not rule out every multi-defect grammar.
