# CARD 2 — FUNDAMENTAL GENERATING-LAW RESULT

**AUTHOR VERDICT: KILL — PENDING INDEPENDENT SCIENTIFIC REVIEW.**

One generating law was frozen at
[`b871073fc2c7e02097dd82657480c6a667a012c4`](https://github.com/ryangrvr/TestingGRUT/commit/b871073fc2c7e02097dd82657480c6a667a012c4)
before numerical evaluation. The law has not been amended. It produces a genuine
conditional nonlinear phase prediction, but an exact symmetry obstruction
prevents its intrinsic interaction from changing any elementary local state.
Initially independent blocks also remain independent. It therefore fails the
record/interaction-generation task and does not generate the required
partition–carrier–readout–interface chain.

The formal first terminal is **S0 CARD-INCOMPLETE**: the frozen submission
discloses absent mandatory certificates, experimental lock, complete price and
fresh sealed kill coverage. The public calculations below document the failed
unchanged proposal. They are not a fresh passing S1–S9 evaluation, a replacement
candidate or independent certification. No sealed holdouts were drawn.

The owner's limited §18.1 development waiver is recorded in
[AUTHORIZATION.md](AUTHORIZATION.md). External evaluation-kit review remains
unestablished. **Two slots consumed; one remains; Card 3 is reserved and
unopened.** Previous scientific packets, Card 1's ledger and canonical records
are unchanged. PR #12 remains a draft.

## 1. THE LAW

### Complete finite model

POSTULATED kinematics: for every integer (n\ge2),
\(\mathcal H_n=(\mathbb C^2)^{\otimes n}\), its elementary tensor factors,
Pauli operators \(Q_i^a\), density matrices, the Born trace rule and a continuous
oriented time coordinate. These primitives are supplied, not derived quantum
structure, observers or physical identities.

For \(\rho\ge0\), \(\operatorname{Tr}\rho=1\), define

\[
r_i^a=\operatorname{Tr}(\rho Q_i^a),\quad
C_{ij}^{ab}=\operatorname{Tr}(\rho Q_i^aQ_j^b)-r_i^ar_j^b,
\quad
F(\rho)=\frac\kappa2\sum_{i<j}\sum_{a,b=x,y,z}(C_{ij}^{ab})^2.
\]

The **new physical postulate** is that this scalar generates the entire
intrinsic Hamiltonian:

\[
\boxed{
H_C(\rho)=\frac{\delta F}{\delta\rho}
=\kappa\sum_{i<j,a,b}C_{ij}^{ab}
\left(Q_i^aQ_j^b-r_j^bQ_i^a-r_i^aQ_j^b\right),
\qquad
\dot\rho=-i[H_C(\rho)+V(t),\rho].}
\]

Here \(V(t)=\sum_{i,a}v_i^a(t)Q_i^a\) is a supplied bounded,
piecewise-continuous local intervention. The parameter is \(\kappa>0\), with
units of inverse time in the supplied \(\hbar=1\) convention. Initial conditions
are **all** density matrices \(\rho(0)=\rho_0\). There are no spatial boundaries:
the model describes one finite causally connected cell, not a relativistic field
theory. The time horizon and control waveform are experimental inputs. The
Hamiltonian has no added bare interaction, thermostat, Gibbs state, bath,
memory register, observer or response kernel.

The choice of a quadratic connected-correlation scalar is the speculative
commitment. Local-basis and permutation invariance motivate it, but do not
uniquely imply it. It is neither derived from data nor an established origin of
physical forces. A correlation is a state property; postulating that it sources
an interaction is additional physics.

### Information price and accounting

[CARD02_FROZEN.md](CARD02_FROZEN.md) contains the unchanged C1–C20 submission.
[L0_CORE.json](L0_CORE.json) explicitly transcribes the finite real/imaginary
matrix derivative arithmetic. Its supplied tensor/quantum/time/domain imports
are outside that arithmetic transcription.

| Disclosed price component | Author calculation |
|---|---:|
| Pointwise arithmetic under the pinned, uncleared oracle | 9,005 bits |
| Explicit Pauli real/imaginary literal tables | 306 bits |
| Local selection tax: 14 Card-2 entries plus 7 prior Card-1 entries | \(\log_2(22)=4.4594316\) bits |
| Sum of those components | 9,315.4594316 bits |
| Complete frozen-language information price | **NOT CERTIFIED** |
| Awarded joint compression \(b_J\) | **0** |

The sum is a **reported partial-price lower bound**, not the total Price or a
compression pass. The oracle is not externally cleared; normal-form/import
unfolding, full axiom encoding, state/protocol specification, range/window costs,
selection maxima and IP-13 exclusions remain unpriced obligations. The global
post-G2-11 search inventory is not certified. No missing component is charged
as zero. A mechanical bit count does not replace a premise-sensitivity audit.

The [attempt ledger](ATTEMPT_LEDGER.json) logs the symbolic draft and 13 public
candidate fixtures before their execution. Identical replay does not append a
new fixture; changed inputs under an existing identifier are rejected. This
documents local fixture accounting, not an independently certified complete
global BU-5/IP-12 inventory. No alternate functional, sign, coefficient sweep,
candidate rescue or Card-3 law was tested.

## 2. THE DERIVATION

### A. The generator is the stated gradient — DERIVED

For a Hermitian tangent \(\delta\rho\),

\[
\delta C_{ij}^{ab}=
\operatorname{Tr}(\delta\rho Q_i^aQ_j^b)
-r_j^b\operatorname{Tr}(\delta\rho Q_i^a)
-r_i^a\operatorname{Tr}(\delta\rho Q_j^b).
\]

Consequently \(\delta F=\operatorname{Tr}(H_C\delta\rho)\). Adding a scalar
multiple of the identity to a constrained gradient would not change the
commutator. The displayed representative fully defines the law.

### B. A finite positive trajectory exists — DERIVED

At fixed \(n\), the autonomous part of the vector field is polynomial in the
matrix entries and locally Lipschitz. With bounded piecewise-continuous controls,
the corresponding piecewise/Carathéodory initial-value problem has a unique
local solution. Hermiticity of \(H_C+V\) gives a unitary propagator along that
solution:

\[
\dot U=-i(H_C(\rho(t))+V(t))U,\qquad
\rho(t)=U(t)\rho_0U(t)^\dagger.
\]

Trace, positivity and the entire spectrum of \(\rho\) are preserved. The density
matrix domain is compact; finite-dimensional boundedness on that domain prevents
finite-time blow-up, so the solution extends over every finite control horizon.
All \(\operatorname{Tr}\rho^m\) are constant. Undriven,
\(\dot F=-i\operatorname{Tr}(H_C[H_C,\rho])=0\).

This proves positivity **of each trajectory**. It does not prove a linear CPTP
map, operational mixture consistency or a causal composite extension.

### C. Exact local-marginal obstruction — DERIVED

Under an independent active local unitary \(W=\bigotimes_iW_i\), Pauli
expectations rotate by real orthogonal matrices \(O_i\):

\[
C_{ij}(W\rho W^\dagger)=O_iC_{ij}(\rho)O_j^{\mathsf T}.
\]

Each Frobenius norm is unchanged. Hence \(F(W\rho W^\dagger)=F(\rho)\).
This is an **active symmetry of the frozen state dynamics**, not just a harmless
coordinate relabeling.

Choose any Hermitian local \(A_i\) and differentiate the symmetry at zero:

\[
0=\left.\frac{d}{ds}F(e^{-isA_i}\rho e^{isA_i})\right|_{s=0}
=-i\operatorname{Tr}(H_C[A_i,\rho]).
\]

Cyclicity of trace therefore gives
\(\operatorname{Tr}(A_i\dot\rho)=0\) for every local \(A_i\) in undriven
evolution. Since these operators span the local Hermitian algebra,

\[
\boxed{\operatorname{Tr}_{\bar i}[H_C(\rho),\rho]=0,
\qquad \dot\rho_i=0\quad\text{for every state and every site}.}
\]

With local interventions the only local-state change is
\(\dot\rho_i=-i[V_i(t),\rho_i]\); all other controls vanish under the partial
trace. Local spectra still cannot change. In particular, the intrinsic law
cannot transfer an excitation, change a local energy expectation, change local
purity, or decohere an elementary local system. The scalar \(F\) also remains
constant under these local interventions by the same symmetry.

This proves a broader obstruction than one bad product fixture. The symmetry
that was meant to avoid a preferred coordinate instead creates too many
conserved physical quantities.

### D. Initial organization cannot bootstrap across a product cut — DERIVED

Let \(\rho=\bigotimes_B\rho_B\) for any partition of the supplied elementary
factors. Every cross-block connected correlator vanishes. The frozen generator
then decomposes as

\[
H_C(\rho)=\sum_B H_C(\rho_B)\otimes I_{\bar B}.
\]

The tensor product of the block solutions solves the full equation with the
same initial state and local controls. Uniqueness proves that this product-block
manifold remains invariant. No correlation or coupling is created across such
a cut. For a fully product state, every \(C_{ij}^{ab}=0\), so \(H_C=0\) and
the undriven state is exactly stationary. Local controls merely rotate its
independent factors.

There are also correlated states invisible to this source. The three-qubit
mixture assigning weight \(1/4\) to \(000,011,101,110\) has zero one- and
two-site correlators but \(\langle Z_1Z_2Z_3\rangle=1\). Its intrinsic
Hamiltonian is zero. Thus “no pair correlation” is not “no physical structure.”

The required failed theorem is explicit:

> The frozen pair-correlation-gradient law can generate interaction-connected
> organization and a record from an initially independent admissible system and
> environment using only its intrinsic dynamics and local controls.

**Counterexample:** every fully product initial state; more generally every
product-block initial state across the proposed new connection. This failure
cannot be repaired by adding initial correlations, seed noise, a bare Hamiltonian
or a new measurement rule to the frozen card. Each supplies new physical content.

### E. A conditional phase formula survives — DERIVED, insufficient

For the undriven pure family

\[
|\psi\rangle=\sqrt p|00\rangle+e^{i\phi}\sqrt{1-p}|11\rangle,
\quad z=2p-1,\quad c^2=1-z^2,
\]

the squared pair-correlation norm is \(2c^2+c^4\), independent of \(\phi\),
and \(F=\kappa(2c^2+c^4)/2\). The local-marginal theorem fixes \(p\).
On this subspace, direct application of \(H_C\) gives the relative phase rate

\[
\boxed{\dot\phi=-4\kappa z(2-z^2).}
\]

For example, at \(\phi=0\), the off-diagonal Hamiltonian element between
\(|00\rangle\) and \(|11\rangle\) is \(2\kappa c\), while the diagonal
difference is \(-4\kappa c^2z\). The relative phase equation is therefore
\(-4\kappa c^2z-4\kappa z=-4\kappa z(2-z^2)\); covariance under local
phase rotations extends it to arbitrary \(\phi\).

This relation genuinely follows from the law. It is not an inserted response
definition. It nevertheless uses supplied entanglement, preparation, qubits,
clock and joint tomography. It establishes no law-generated persistent identity,
accessible memory, carrier or interface.

### F. The required chain fails

| Required structure | What this attempt actually establishes |
|---|---|
| Primitive quantum relata | POSTULATED tensor factors and Pauli algebra |
| Persistent physical differentiation \(\Pi\) | UNESTABLISHED; a correlation graph decodes supplied initial organization, cannot designate a generated S/E split |
| Independently certified common carrier | UNESTABLISHED; reference extractor returns failure |
| Accessible readout class \([h]_{T_\Pi}\) | UNESTABLISHED; diagnostic Pauli tomography is supplied |
| Interface freedom \(T_\Pi\) | UNESTABLISHED; basis covariance is not a generated observation interface |
| Observable response \(\Gamma_\Pi\) | Conditional supplied-observable trajectories exist; the required chain does not |
| R1 reciprocity \(\epsilon_R\) | NOT COMPUTED because the upstream chain fails |

The numerical graph output is illustrative, not an exact-zero/persistence
certificate. No favorable threshold, observer tie-break or replacement metric
is inserted. Bottom outputs are failures, not successful structures.

## 3. THE RESTRICTION

The strongest exact restriction of this formal law is conservation of every
local marginal under intrinsic evolution, together with invariant independence
across every initially product cut. These are nondefinitional, measurable
restrictions on dynamics. They exclude ordinary record-writing and exchange
processes, but for the wrong reason: the proposed interaction is unable to
perform the task it was designed to generate.

The frozen primary diagnostic relation is also quantitative. At
\(p_1=3/4\) and \(p_2=5/8\),

\[
\omega_1=-\frac72\kappa,\quad
\omega_2=-\frac{31}{16}\kappa,\quad
\boxed{56\omega_2-31\omega_1=0,\quad \omega_2/\omega_1=31/56.}
\]

This removes actual observable freedom in a **formal replacement-dynamics
model**. To see the baseline freedom explicitly, restrict an ordinary fixed
Hamiltonian to the \(|00\rangle,|11\rangle\) subspace:

\[
H_{b,c}=c(|00\rangle\langle00|-|11\rangle\langle11|)
+b(|00\rangle\langle11|+|11\rangle\langle00|).
\]

At initial \(\phi=0\), its instantaneous phase rate is

\[
\omega(p)=2c-\frac{2bz}{\sqrt{1-z^2}}.
\]

The two chosen \(z\)'s give distinct coefficients of \(b\), so
\((b,c)\mapsto(\omega_1,\omega_2)\) has nonzero determinant. Even the small
resource region \(|b|+|c|\le\kappa\) has an observable image with nonzero area;
the candidate's rates lie on one fixed ray when \(\kappa\) varies, or one point
when it is independently fixed. The unrelated microscopic coordinates need not
be uniquely selected. This argument concerns initial rates; the conventional
off-diagonal Hamiltonian need not preserve \(p\) at later times.

**It is not an earned GRUT exclusion theorem over the complete baseline.**
The full conventional operational baseline also includes consistent randomized
preparations and linear CPTP operations. The candidate fails that extension
below. It has not supplied a physically admissible candidate family satisfying
all the baseline requirements, nor generated the charter's joint coordinates.
Discarding those constraints would change the scientific comparison.

## 4. THE COUNTERMODEL

### Same operational preparation and finite resources

The conventional comparison uses the same two qubits, product or Schmidt
preparations, time units, local diagnostics, finite horizon and norm ceiling.
For Pauli pairs, covariance Cauchy–Schwarz gives \(|C_{ij}^{ab}|\le1\), and
each Hamiltonian bracket has norm at most 3. Thus

\[
\|H_C\|\le27\kappa\binom n2.
\]

For two qubits the shared ceiling is \(27\kappa\). All displayed conventional
Hamiltonians lie below it; no ancilla, bath, state-feedback controller, refitting
or unbounded resource limit is needed. The baseline permits conventional
interactions; the candidate's claim is to generate their physical role from
state structure instead. Giving it an extra entangling control would evade that
claim, not test it.

The comparators are ordinary finite unitary evolutions, positive on all states
and extensions, with conserved time-independent Hamiltonian energy. They reside
in one causally connected cell. Relativistic remote-control claims are absent.
KMS/FDT/Onsager constraints are not imposed on these coherent nonequilibrium
preparations without their equilibrium hypotheses. Conservation, causality and
resource assumptions are matched where applicable.

| Fixed conventional model | Common preparation and duration | Conventional result | Frozen candidate |
|---|---|---|---|
| \(H=gZ_AZ_B\), \(g=\kappa\) | \(|++\rangle\), \(t=\pi/(4\kappa)\) | Local visibility \(V=0\), conditional E record distinguishability \(D=1\) | \(V=1,D=0\) |
| \(H=g(X_AX_B+Y_AY_B)/2\), \(g=\kappa\) | \(|10\rangle\), \(t=\pi/(2\kappa)\) | \(|01\rangle\), B excitation probability 1 | \(|10\rangle\), B excitation probability 0 |
| \(H=-(7\kappa/4)Z_A\) | The two frozen Schmidt preparations | Both relative phase rates \(-7\kappa/2\); ratio 1 | Rates \(-7\kappa/2,-31\kappa/16\); ratio \(31/56\) |

For record writing, the conditional E vectors are
\(e^{\mp i\pi Z/4}|+\rangle\), which are orthogonal; their overlap gives zero
system visibility. For exchange,
\(|10\rangle\mapsto\cos(gt)|10\rangle-i\sin(gt)|01\rangle\).
These are physical changes at a fixed readout, not coordinate transformations.
Here S/E designation and readouts are explicit comparator diagnostics, not
claimed generated interfaces or a calculation of \(\epsilon_R\).

### An additional operational obstruction: non-affinity

Let \(\rho_0=|00\rangle\langle00|\),
\(\rho_1=|+1\rangle\langle+1|\), and
\(\rho_m=(\rho_0+\rho_1)/2\). Both pure product inputs are stationary.
For their mixture the frozen formula gives

\[
H_C(\rho_m)=\frac\kappa2(Z_A-X_A)Z_B,
\qquad
\left.\frac d{dt}\langle Y_AZ_B\rangle\right|_{\rho_m}=\kappa,
\]

whereas the average of the two evolved product-branch derivatives is zero.
Therefore no linear CPTP evolution agrees with the frozen law on all three
preparations. This is an exact witness, not a negative-Choi calculation.

Usual branchwise randomized-preparation predictions and evolution of the density
matrix give incompatible answers. A different ontology could try to deny that
identification, but this card does not supply such a measurement/composite
theory, and it may not add one after failure. The established no-signaling
arguments have explicit entanglement, trace-rule and mixture/local-extension
hypotheses [1]; alternative nonlinear extensions must be analyzed on their own
assumptions [2,3]. **This witness alone is not a proved superluminal signal.**

### Novelty audit of the assembled proposal

| Framework | What was supplied or reused; what is not earned |
|---|---|
| Dynamical systems | Finite polynomial, isospectral Hamiltonian flow; existence and Noether conservation follow from that structure, not new physics |
| Nonlinear quantum/mean-field dynamics | State-dependent Hamiltonian evolution belongs to an established broad class [4,5]; no literature-wide originality claim for this functional |
| Predictive-state theory | State, intervention coordinates and tomography are inputs; no minimal physical memory or observer is generated |
| Causal inference | Connected correlation is a supplied-state statistic, not an independently certified causal interaction selector |
| Information theory | Covariance norms and trace distinguishability are conventional diagnostics; defining them creates no physical record |
| Open-system physics | The fundamental model is closed and isospectral; partial traces do not supply an operational CPTP extension, bath or irreversible record dynamics |
| Statistical mechanics | A conserved scalar is not an attracting free energy; no thermal state, thermostat, equilibration or selection of initial organization is derived |

The speculative content is the universal correlation-gradient interaction
postulate. Combining established equations does not establish its truth or
originality. Its exact new restrictions are reasons to reject this attempt,
not evidence that nature obeys it.

## 5. THE TEST

### Reproducible computation

From the repository root, using the pinned requirements:

```bash
bash research/card02/reproduce.sh
```

This verifies the frozen SHA-256 hashes, logs/replays the unchanged fixtures,
uses a separate exact Gaussian-rational matrix implementation, runs the frozen
NumPy reference and integrates the two prescribed phase trajectories. It emits
[results.json](results.json), [CONTROL_TABLES.md](CONTROL_TABLES.md) and
[price_report.json](price_report.json).

| Check | Observed author result |
|---|---:|
| Exact rational control assertions | 111 |
| Local-marginal derivative cells within those controls | 66 exactly zero |
| Exact/reference RHS discrepancy on rational fixtures | 0 |
| Separate numerical local-drive comparison | 1 passes |
| Phase trajectories | Two, nine output times each, \(0\le t\le1\), \(\kappa=1\) |
| Largest phase density-matrix discrepancy | \(1.82656\times10^{-11}\) |
| Local logged Card-2 entries including symbolic draft | 14 |

The phase integrator uses DOP853 with rtol \(10^{-11}\), atol \(10^{-12}\).
Its tolerance checks and small eigenvalue roundoff are numerical evidence, not
certified interval enclosures. The continuum proofs do not depend on a finite
test grid. The exact assertion count includes local conservation, rational
gradient/phase/covariance identities, block closure and explicit conventional
unitary/writing checks; it is not a count of independent discoveries.

The failure-control harness was corrected to classify its floating local-drive
comparison separately, use a nonzero gradient tangent, and explicitly verify the
exchange comparator's unitary.
These audit corrections did not modify a frozen law file, change a candidate
fixture or add a replacement candidate. Reproduction checks all six frozen
source hashes before evaluating anything. Output replay is deterministic in
this environment; timestamps in the ledger are retained, not regenerated.

### Operational falsification, not an experimental lock

If the supplied elementary factors could independently be identified with two
physical qubits, local preparations and Pauli tomography could measure the
two initial phase rates. Use the first rate to calibrate \(\kappa\), then
predict the second without refitting. With certified absolute error bounds
\(e_1,e_2\), reject the phase relation when

\[
|56\widehat\omega_2-31\widehat\omega_1|>56e_2+31e_1.
\]

The bounds must include state preparation, finite-time rate extraction, timing,
leakage, tomography, residual conventional interaction and model discrepancy.
Certified hypothesis-test coverage and power would still be required; they are
not provided by this deterministic inequality. Record-writing and excitation
transfer are separate, simpler falsifications of the product-state prediction.

Published superconducting-qubit work provides a conventional controlled-phase
platform precedent [6]. It is not our run-specific calibration, a mapping from
device qubits to fundamental relata, a held-out result or an LT-1 certificate.
No laboratory data, apparatus access, sample-size/power certificate, independent
calibration or generated interface exists here. Experimental feasibility and
the full rival discrimination are **not certified**. Programming a simulator
to implement \(H_C(\rho)\) would test a controller, not fundamental dynamics.

The conditional phase-rate ratio is not known to distinguish the candidate from
every possible conventional preparation-dependent controller or hidden
environment. The explicit fixed Hamiltonian is a hostile admissible witness;
it is not a complete lock-side conventional panel.

## 6. THE VERDICT

**KILL at author level.** The primary finite phase consequence is true, but the
generating mission fails. The exact local symmetry preserves all local states;
product blocks cannot acquire new interactions; record writing requires
organization already encoded in the input. The required identity/carrier/readout
chain, operational composite extension and qualifying joint-observable
restriction are absent. The full price and experimental lock are also absent.

Evidence grade: author analytic proofs and exact rational controls; supplementary
numerical trajectories; **no independent scientific review and no experiment**.
The evaluation-kit waiver permitted development, not a certification. The
formal first terminal remains S0 CARD-INCOMPLETE; missing prerequisites are
hostile failures, never favorable unknowns. No BANK or validated-new-physics
claim is made, and no canonical KILL/bank/state record is edited.

The failed theorem and counterexamples apply to this one frozen law, not to
every possible theory of interactions. **Card 2 is consumed and author-killed;
two total attempts consumed; one remains. Card 3 is unopened.** No repair or
third attempt is included in this packet.

### Primary sources and read scope

1. C. Simon, V. Bužek and N. Gisin,
   [The no-signaling condition and quantum dynamics](https://arxiv.org/pdf/quant-ph/0102125).
   Full four-page argument read: its mixture/trace/entanglement assumptions and
   local identity-extension condition are relevant; it is not asserted as a
   hypothesis-free exclusion of every nonlinear ontology.
2. J. Polchinski,
   [Weinberg's nonlinear quantum mechanics and the Einstein–Podolsky–Rosen paradox](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.66.397).
   Publisher abstract only; contextual warning about nonlinear extensions, not
   a substituted proof of a signaling protocol for this law.
3. M. Czachor,
   [Nonlocal-looking equations can make nonlinear quantum dynamics local](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.57.4122).
   Publisher abstract only; confirms that composite-extension assumptions matter.
4. S. Weinberg,
   [Testing quantum mechanics](https://www.sciencedirect.com/science/article/pii/0003491689902765),
   Ann. Phys. 194, 336–386 (1989). Publisher abstract retrieved through search;
   direct page access returned 403. Used only for the existence of a prior
   nonlinear Hamiltonian framework, not a claim of exact equivalence.
5. M. Czachor and J. Naudts,
   [Microscopic foundation of nonextensive statistics](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.59.R2497).
   Publisher abstract only: nonlinear Lie–Poisson density-matrix evolution is
   already an established mathematical construction.
6. R. Barends et al.,
   [Logic gates at the surface code threshold](https://arxiv.org/pdf/1402.4848).
   Main controlled-phase architecture and gate sections read; finite conventional
   platform precedent, not fresh GRUT apparatus evidence.

The package's mathematical claims are proved above and controlled in its code.
No paper is quoted at length, redistributed, or credited as validating this law.
