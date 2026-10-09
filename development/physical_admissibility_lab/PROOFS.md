# Physical embeddings, exclusions and resource bounds

**DRAFT — same-author proofs and controls; external scientific review pending.**
All results below are consequences of established comparator premises, not
GRUT laws. No candidate K, Stage-3 price, score or Card was produced. Sources
P01–P12 are identified in [SOURCES.json](SOURCES.json). Finite arithmetic
fixtures are in [RESULTS.json](RESULTS.json); they check instances, not the
universal quantifiers in these proofs.

## Operational contract

A setting a is an external intervention, y an accessible record, z retained
state. W(z',y|z,a) is a normalized nonnegative causal kernel. A policy may
choose a from all earlier settings/records and independent randomization.
Equivalence means equality of every record distribution under the **same**
preparations, policy and readout. A comparator may add hidden state, but may
not let it know a future randomized setting or change what the intervention
does. Geometry, physical elapsed time, energy units, accessible outputs,
resource boundaries and error tolerances must be independently frozen when
those premises are invoked. Descriptive state enlargement does not establish
their physical implementability.

POSTULATED throughout: the specified state spaces, preparation, operators,
intervention and apparatus. DERIVED below: consequences under each listed
premise. There is no derivation of a physical subsystem identity or detector.

## T1 — finite classical kernels have a fixed CPTP instrument

Let Z,A,Y be finite. On H_Z=C^|Z| define, for each a,

$$K^a_{z',y,z}=\sqrt{W(z',y|z,a)}\,|z'\rangle\langle z|,
\quad \mathcal I^a_y(\rho)=\sum_{z,z'}K^a_{z',y,z}\rho K^{a\dagger}_{z',y,z}.$$

Each outcome map is completely positive. Normalization gives
Σ_y,z,z' K†K=I, so Σ_y I_y is trace preserving. Starting with a diagonal
state, the unnormalized conditional output is diagonal with weights
Σ_z p_z W(z',y|z,a), exactly the classical update. Induction on record length
proves equality under **every** adaptive policy; randomized policies follow
by convex mixing. No new intervention or relabelled readout is used.

If the record is retained as an orthogonal output register, the joint channel
has Kraus operators sqrt(W)|z',y><z|. Its Choi operator is diagonal on
|z>⊗|z',y>, with nonzero eigenvalues precisely the positive W entries. Let
r be their count. Its Choi rank, minimum Kraus count and minimum pure
Stinespring environment dimension equal r [P01]. Explicitly

$$V_a|z\rangle=\sum_{z',y}\sqrt W\,|z',y\rangle\otimes|z,z',y\rangle_E$$

is an isometry; its images are orthogonal for distinct z and normalized.
For a chosen common finite output/input ambient space it can be extended
to a unitary with initialized output/environment registers. A Hermitian
logarithm supplies a finite-dimensional Hamiltonian for a clocked stroke.
This statement supplies neither a local Hamiltonian, an energy-conserving
thermal implementation, nor a controller-free autonomous apparatus.

For PR #6's eight-state equations the two successor weights are 3/4 and 1/4.
There are 16 nonzero Choi eigenvalues for either setting, hence r=16 for
**this explicitly dephasing extension on arbitrary quantum inputs**.
This is NOT a minimum over every coherent channel agreeing only on
classical preparations and records. Such extensions have extra freedom.
The 1,024 depth-three feedback/preparation comparisons (8,192 cells)
match over exact rationals. The induction, not this enumeration, proves
all finite record lengths.

Price of this comparator construction: a supplied eight-dimensional
state, two fixed instruments, their amplitudes, preparation, readout,
external setting/clock and initially independent ancillas. Fresh environments
of dimension r per stroke give the crude finite-horizon upper bound r^N.
Resetting or recycling them requires further resources. No perpetual
bounded closed apparatus, computable implementation of arbitrary real W,
relativistic completion or origin of these primitives is established.

## T2 — the actual four-word register needs d≥4 quantum memory

The inherited control writes two input bits into (m1,m2), with
(m1,m2)→(m2,a), and reports old m1. Four preparations 00,01,10,11 and the
common probe 00 give the response matrix I4.

In any fixed quantum comparator let ρ_i be the retained state after
preparation i and {E_j} the POVM induced by the whole common probe. Fresh
independent ancillas, adaptive probe operations and discarded outputs can
all be absorbed into that POVM. Exact agreement requires
Tr(E_jρ_i)=δ_ij. Since 0≤E_j≤I, support(ρ_i) lies in the eigenvalue-one
space of E_i and in the kernels of E_j for j≠i. These supports are mutually
orthogonal; thus d≥Σ_i rank(ρ_i)≥4. Four orthogonal states and a common
measurement attain the bound. A fixed four-state classical register also
attains it. This is standard state discrimination [P01,P13].

The finite-error version for M equiprobable preparations is particularly
simple:

$$p_{\rm success}=\frac1M\sum_i\mathrm{Tr}(E_i\rho_i)
\le\frac1M\sum_i\mathrm{Tr}E_i=\frac dM.$$

It uses ρ_i≤I and ΣE_i=I. For M=4, a certified success probability greater
than 3/4 excludes all retained d≤3; greater than 1/2 excludes a qubit.
A one-sided confidence bound must exceed the threshold, including declared
leakage/readout error. Generic response rank≤d² is weaker: four linearly
independent qubit density operators exist, but cannot be perfectly decoded
by one common probe. There is no contradiction with PR #6's rank warning.

All preparation-dependent retained environment, classical labels,
controller/clock information and side channels accessible to the probe
must count toward d. An uncounted copy of the preparation trivially defeats
the bound. The device's capacity cannot be inferred just by calling it a
qubit. This is a measurable resource exclusion, not a universal GRUT law.

## T3 — reversible finite control and irreversible memory cost

For any finite deterministic update f(z,a,b), the map

$$(z,a,b,s)\mapsto(z,a,b,s\mathbin{\mathrm{xor}}f(z,a,b))$$

on bit-string scratch is a permutation, indeed an involution. Initialize
s=0 and it computes the update while retaining input and noise. The
eight-state fixture has 256 distinct domain points and 256 images. Its
finite Boolean functions can be decomposed into reversible gates; nearest
neighbor swaps route them on a finite lattice, with supplied gate schedule,
latency and scratch. Probability 1/4 can be generated by two fresh unbiased
bits (b=1 only for 11). Reversible embeddings do not produce free garbage
disposal. This is an exact finite control construction; it is not a proof
of arbitrarily precise smooth classical Hamiltonian hardware with bounded
time/energy, or a local, autonomous, perpetual realization of arbitrary W.

For exact reset on arbitrary d-dimensional inputs with pure environment
|e0>, a closed unitary must send |i>|e0> to |0>|ei>. Unitarity makes the
|ei> orthonormal, so dim E≥d. For N independent input registers, all
reset exactly, all cyclic controller registers returned to the same pure
state, and no other retained information, the same argument on d^N
orthogonal inputs gives

$$\dim E\ge d^N.$$

Output records or garbage carrying the inputs are part of E. If only a
restricted sequence rather than d^N independently preparable inputs is
admissible, this bound does not follow. Pure ancillary swaps attain the
one-reset bound. The conclusion concerns closed reversible realization,
not a universal heat cost for every memory update.

**Finite faithful thermal bath obstruction.** For d>1, initial I_d/d⊗τ_E,
finite dim E=D, and a finite-temperature finite-energy Gibbs state τ_E>0,
the joint state has rank dD. Unitary evolution preserves rank. A pure
reset memory requires joint output |0><0|⊗σ_E (a pure marginal necessarily
factorizes), of rank at most D. Contradiction. Exact pure reset is therefore
impossible under these resources. Approximation, infinite baths, pure
ancillas or additional work/entropy reservoirs change the premise and must
be counted. A mixed faithful extra finite reservoir does not repair rank;
an initially pure extra system can. This is an instance of the standard
finite-bath erasure boundary [P02], not a new third law.

## T4 — Landauer equality with a priced finite example

Let ρ_SE=ρ_S⊗τ_E, τ_E=e^(-βH_E)/Z, β>0, finite dimensions, and a global
unitary U. Define S(ρ)=−Trρ logρ, ΔS=S(ρ_S)−S(ρ'_S),
Q=TrH_E(ρ'_E−τ_E), I'=S(ρ'_S)+S(ρ'_E)−S(ρ'_{SE}), and
D(ρ||τ)=Trρ(logρ−logτ). Entropy conservation and the initial product give
S(ρ'_E)−S(τ_E)=ΔS+I'. Substituting logτ_E=−βH_E−logZ yields

$$\boxed{\beta Q=\Delta S+I'+D(\rho'_E\Vert\tau_E)\ge\Delta S.}$$

This is Reeb–Wolf's equality [P02], rederived here. Its independence,
thermal-state, closed-unitary and reservoir accounting premises matter;
initial correlations or nonthermal ancillas can alter the bound.

Concrete stroke: β=1, H_S=0, H_E=diag(0,ln3), initial memory I2/2 and
thermal bath diag(3/4,1/4), U=SWAP. Final memory is that mixed bath state,
final bath I2/2, and I'=0. Exactly,

$$\beta Q=\tfrac14\ln3,\quad
\Delta S=\ln2-h(1/4),\quad
D=\tfrac12\ln(4/3),\quad \beta Q=\Delta S+D.$$

The stored floating residual checks that exact log identity numerically.
It is evidence, not a proof of a universal entropy inequality. With H_S=0,
SWAP does not conserve bare H_S+H_E; external work Q is supplied and
priced. The example is not a free thermal operation or exact pure erasure.
Choosing a pure reset ancilla makes exact SWAP reset possible but abandons
the faithful thermal premise; it does not refute T4.

## T5 — stationary detailed-balance lifts cannot hide oriented pair current

Suppose a classical hidden chain is stationary μ and obeys
μ(z)P(z,z')=μ(z')P(z',z). Use a fixed time-reversal-even passive readout
h_i(z), with independent emission noise absorbed into the hidden state or
conditional factors. Then

$$p(i,j)=\sum_{z,z'}\mu(z)P(z,z')h_i(z)h_j(z')=p(j,i).$$

The equality follows by swapping z,z' and detailed balance. Thus no such
lift, whatever its hidden dimension, realizes the stationary visible
qutrit kernel P=(I+C)/2, C a forward three-cycle: π=uniform and
π_0P_01−π_1P_10=1/6. Extra memory cannot remove this obstruction while
retaining the passive stationary time-even readout premise.

**Hostile physical countermodel to an overbroad conclusion.** Take degenerate
H_S=0 on a qutrit and H_E=0 on a bath qubit, τ_E=I2/2. The unitary
U=I⊗|0><0|+C⊗|1><1| commutes with H_S+H_E. Discarding E gives
Φ(ρ)=(ρ+CρC†)/2 and realizes exactly the same P on labelled diagonal
states (transpose C if using row rather than column vector convention).
This is a conventional energy-conserving thermal operation [P03]. It
preserves Gibbs but fails classical detailed balance. Repeated strokes
use refreshed bath bits and a supplied timing/controller, not an autonomous
equilibrium DB semigroup. Hence Gibbs preservation, thermal operations,
passivity of a state and dynamical detailed balance are different claims.

A driven three-cycle CTMC with forward rates 2, backward 1 and uniform π
has entropy production Σ_edges(2−1)ln(2)/3=ln2 per unit time. Standard
nonequilibrium driving provides circulation once affinities/resources are
included. It is not an exact finite-time realization of the zero-entry
discrete P. Generalized equilibrium time reversal can reverse momenta or
magnetic fields; T5's plain DB/time-even hypothesis must not be applied to
those observables without that transformation. Quantum invasive records
also need their own instrument-level reversibility premise.

## T6 — bounded-rate CTMC excludes exact finite-time NOT

For any finite generator Q, Q_ij≥0 off diagonal and Q1=0,
det(e^(tQ))=e^(t TrQ)>0. The visible two-state NOT matrix has determinant
−1, so it is not e^(tQ) at finite t. A bounded time-dependent generator has
detP=exp∫TrQ(s)ds>0 by Liouville's formula [P04].

Enlarging a finite CTMC while keeping a fixed initially certain bit
readout also cannot implement an exact NOT: a no-jump trajectory retains
the original hidden state and readout with positive probability. If every
exit rate is at most Λ, its probability is at least e^(−ΛT). Therefore
the finite-time flip error is at least e^(−ΛT), even with hidden mixtures.
This assumes a fixed certain readout, bounded rates and no external
deterministic reset or continuous drift that changes the readout on the
no-jump path. Such changes leave the comparator class.

The exclusion is not a ban on physical NOT. H=πσ_x/(2T) gives
U(T)=−iσ_x and flips the computational basis exactly. Clocked reversible
gates also implement NOT. A CTMC premise needs justification; a hidden
clock with deterministic continuous drift is not a finite bounded-rate CTMC.

## T7 — positivity/causality alone do not establish quantum admissibility

Qubit transposition is positive and trace preserving, but its action on
half of a Bell pair produces eigenvalues (1/2,1/2,1/2,−1/2); the singlet
expectation is −1/2. It is not CP [P01]. This exclusion invokes the same
operation on independently supplied entangled ancillas. Reduced maps
defined only on an initially correlated restricted domain need not meet
that extension premise. T1's classical kernels always have a CP extension.

A PR box has p(a,b|x,y)=1/2 for a xor b=xy, otherwise zero. Marginals are
uniform, but CHSH=4. For separated tensor-factor quantum ±1 observables,
B=A0⊗(B0+B1)+A1⊗(B0−B1) satisfies
B²=4I−[A0,A1]⊗[B0,B1]. Operator norms of the commutators are at most 2,
so ||B||≤2sqrt2 [P05]. A PR box is not a quantum realization, in any
dimension, under these setting-independence and no-communication premises.
Sixteen deterministic local assignments give |CHSH|≤2. These are standard
exclusions; neither no-signalling nor an abstract fixed stochastic lift
implies quantum realizability.

## T8 — locality: strict circuit cones versus Hamiltonian tails

A radius-r classical lattice update has site-j state after N steps
depending only on initial inputs within rN, by induction. A fresh local
setting outside that cone cannot affect the remote record when the common
initial state is independent of that setting. Shared initial correlations
do not convey a new randomized setting. Finite-depth quantum circuits have
the analogous exact gate cone. A global copy-setting kernel can violate
these constraints; renaming its remote dependency as hidden memory does
not make it local.

Do not transfer that exact cone to continuous finite-range quantum lattice
Hamiltonians. In units hbar=1 let
H=J(X1X2+Z2Z3), A=X3, B=Z1. The two terms anticommute, H²=2J²I, and
with θ=sqrt2 Jt direct multiplication gives

$$A(t)=\cos^2\theta X3-\sqrt2\cos\theta\sin\theta Z2Y3
-\sin^2\theta X1Y2Y3,$$
$$[A(t),B]=2i\sin^2\theta Y1Y2Y3,\quad
\|[A(t),B]\|=2\sin^2\theta=4J^2t^2+O(t^4).$$

At θ=π/2 the identities H A H/(2J²)=−X1Y2Y3 and
[A(t),B]†[A(t),B]=4I are checked exactly. The tails are nonzero at
arbitrarily early nonzero times despite only nearest-neighbor terms.
Lieb–Robinson bounds control their exponential envelope, e.g.
C||A||||B||e^(−μ(distance−v|t|)), with constants depending on the supplied
interaction norms/geometry [P06]. They are not exact relativistic cones.

In a relativistic local algebra, spacelike commutativity does give no
unconditional signalling by localized trace-preserving operations:
Σ_r K_r† B K_r=B Σ_r K_r†K_r=B. Selective conditioning can change
conditional remote distributions but requires conveying the outcome to
use it. A lattice embedding is not a Lorentz/QFT completion.

## T9 — standard physical constraints leave measurable kinetic parameters

For a positive frozen stationary vector π on N labelled states and a
complete graph, reversible positive CTMC generators are in one-to-one
correspondence with E=N(N−1)/2 positive symmetric conductances:
Q_ij=c_ij/π_i (i≠j), Q_ii=−Σ_j≠i Q_ij. Conversely DB gives c_ij=π_iQ_ij.
The open positive cone has dimension E after normalization, conservation
and DB. Knowing rates, Hamiltonian barriers or bath couplings removes
freedom; those are additional physical information. This is a baseline
dimension count, not another kinetic no-go family.

A quantum control makes the same boundary concrete. Take H=Δ|1><1|,
βΔ=ln3, σ_-=|0><1|, D[L](ρ)=LρL†−{L†L,ρ}/2, and

$$\dot\rho=-i[H,\rho]+\gamma D[\sigma_-](\rho)
+\tfrac\gamma3 D[\sigma_+](\rho)
+\tfrac{\gamma_\phi}{2}D[\sigma_z](\rho),\qquad
\gamma>0,\ \gamma_\phi\ge0.$$

This fixed GKLS semigroup is CP, has Gibbs populations (3/4,1/4), thermal
up/down ratio 1/3 and positive relaxation. Its population decay rate is
4γ/3; its coherence decay rate is 2γ/3+γ_phi. The dissipative part obeys
the usual thermal detailed-balance relation; Hamiltonian evolution is
included separately. Rates depend on ordinary reservoir spectra/couplings
[P07,P29]. Quantum positivity and the thermal ratio do not determine them.
Microscopic weak-coupling derivations have their standard validity regime;
no exact finite-bath all-time semigroup or relativistic completion is claimed.
The prior three-state kinetic manuscript similarly leaves projector
orientation free within its reversible coarse class. This audit does not
enlarge that theorem family or claim all microscopic baths realize it.

## T10 — conditional EFT positivity removes real freedom

For forward identical massive scalar scattering, define ν=s−2m² and the
pole-subtracted crossing-even amplitude F(ν). Assume a mass gap, real
analyticity, cuts |ν|≥ν_th=2m², optical-theorem positivity ImF≥0 on the
right cut, and sufficient high-energy growth control for a twice-subtracted
dispersion relation. Cauchy's formula, deforming the contour to both cuts
and using crossing, gives

$$a_2=\tfrac12F''(0)=\frac2\pi\int_{\nu_{th}}^\infty
\frac{\operatorname{Im}F(\nu+i0)}{\nu^3}\,d\nu\ge0.$$

It is strict if the absorptive part is nonzero on positive measure. Thus a
negative forward ν² coefficient is excluded by this UV comparator [P08].
The sign is a restriction on independently measurable scattering data,
not a definition. An arbitrary negative polynomial is not a countermodel
obeying the full hypotheses. This conditional result is known physics,
not new GRUT reduction. Massless/gravitational forward poles and failed
growth bounds require separate treatment; no gravity positivity verdict
is inferred here. This theorem has no finite numerical proof fixture.

## Hostile boundaries retained

1. T1 survives every classical finite feedback setting, but fails to prove
   bounded autonomous thermal/local/relativistic resources.
2. T2 disappears if the probe receives an uncounted preparation label.
3. T3's d^N bound requires independently preparable inputs and exact reset;
   it cannot be attached to every correlated history or memory write.
4. T5's DB obstruction is broken by an ordinary thermal-operation stroke;
   using Gibbs preservation as a synonym for DB is false.
5. T6's CTMC obstruction is broken by an ordinary qubit Hamiltonian.
6. T8's exact lattice light-cone claim is killed by a three-spin Hamiltonian;
   the correct continuous-time comparator uses an LR envelope.
7. T9's remaining rates are ordinary constitutive/bath inputs, not evidence
   of nature lacking a law. Complete specified microphysics can determine them.
8. None of T1–T10 derives subsystem identity, readout, clock or a new joint
   identity/response generator. All obtained restrictions have a standard
   framework source or direct conventional derivation.
