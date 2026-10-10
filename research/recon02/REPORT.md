# Recon02: a measurable response/record target, without a new selector

10 October 2026. Authorized development branch; author audit. **Card 2 remains
unopened.** One frozen slot consumed, two remain; external evaluator/CR-5 review
is still unmet. This is a bounded §15.6 NULL-PiT conventional comparison, not a
fourth positive-control mechanism, a new foundations lane or a GRUT law attempt.
See [SCOPE.md](SCOPE.md).

The operational pair is system coherence retained after an interaction and
information about the system's conserved binary label in one fixed accessible
memory fragment. Standard quantum theory permits the entire region

\[
S_P=\{(V,D_F):0\le V\le1,\quad0\le D_F\le\sqrt{1-V^2}\}.
\]

The upper boundary is known complementarity. The interior is real conventional
freedom, even with the same pure preparation, two memory qubits, conserved
system populations, fixed readout and bounded interaction resources. The current
GRUT record provides no additional equation selecting a proper subset of this
region. **Target found; GRUT bridge not established; HOLD Card 2.**

## 1. Operational baseline P, supplied in full

There are three labelled qubits S, F, R, with Hilbert space C²⊗C²⊗C². F is the
fixed accessible record; R is a retained additional memory outside that readout.
Their names and access status never change between comparators. Initial memory
state is |00⟩FR. For visibility runs S starts in |+⟩; for record calibration the
intervention sets S to |0⟩ or |1⟩ with equal prior probability. These are separate
repeated runs with identical memory preparation and interaction settings.

Born probabilities, tensor factors, computational bases, clock, controls and
tomography apparatus are POSTULATED conventional inputs. The system's bare
Hamiltonian is H0=ℏω ZS/2; memory bare Hamiltonians are zero. Controlled unitaries
preserve ZS and the bare energy. In the interaction frame, admissible evolutions
have the form U=P0⊗U0+P1⊗U1, where Pi=|i⟩⟨i| and U0,U1 are unitaries on FR.
Their microscopic values are baseline freedoms. No thermal equilibrium state,
KMS bath, FDT spectrum or Onsager transport regime is invoked. Thus these are
not omitted constraints on an equilibrium calculation.

The common resource ceiling is two initially pure memory qubits, no fresh
ancillas, no feedback or postselection, and local gates on the S-F-R chain.
Allow six entangling two-qubit gates, eight single-qubit gates and total time T.
Allocate each elementary pulse at least T/16 and use Hermitian pulse logarithms
with eigenphases in [-π,π]; interaction norm is at most 16πℏ/T. Identity padding
is allowed. Both the product family and environmental encoding below fit this
same ceiling. This fixes capacity and a finite control budget, not equality of
each microscopic pulse or an autonomous origin for the controller. Zero-energy
memories and their purity are supplied physical resources, not free memory
generation. Controller preparation/reset costs are not derived or claimed zero.

Nearest-neighbour synthesis is explicit: a controlled rotation on R uses
SWAP(F,R), controlled rotation on F, SWAP(F,R). A controlled rotation on F uses
one further two-qubit gate. The optional F-R encoding uses one CNOT and local
gates. Sequential finite-time gates do not claim relativistic spacetime/QFT
completion or arbitrarily fast distant influence. All record-forming and
environment-encoding interactions commute with H0. System preparation and
analysis pulses are supplied interventions whose work is part of the control
budget. The bare-energy conservation check does not certify the controller's total
thermodynamic work budget. Stronger microscopic restrictions would define a
different P and require a fresh admissibility audit.

Let |ei⟩=Ui|00⟩. The reduced channel obeys
ρ01↦qρ01, q=⟨e1|e0⟩, and ρii↦ρii. Remove a calibrated coherent phase and define

\[
V=|q|,\quad \rho_{F|i}=\operatorname{tr}_R|e_i\rangle\langle e_i|,
\quad D_F=\tfrac12\|\rho_{F|0}-\rho_{F|1}\|_1.
\]

D_F is an independently measured distinguishability, not a definition in terms
of V. Equal-prior optimal guessing probability is (1+D_F)/2. Full F tomography
is admitted in both comparators; D_F=0 means **no** measurement on F can identify
the label better than guessing. It does not mean a deliberately bad detector
was selected. D_F concerns a prepared classical label, not arbitrary quantum
state copying, phenomenal experience or a derived Born outcome.

## 2. Exact observable image: both inclusion and realization

For pure conditional states on FR, their difference has nonzero eigenvalues
±√(1−|⟨e1|e0⟩|²). Hence D_FR=√(1−V²). Partial trace cannot increase trace
distance: the variational formula D(ρ,σ)=max(0≤M≤I) tr[M(ρ−σ)] embeds every F
effect as M⊗IR. Therefore D_F≤D_FR. This proves the upper inclusion for the
whole stated controlled-unitary family, not just numerical fixtures.

To attain every point, define

\[
R(\theta)=e^{-i\theta Y}=
\begin{pmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{pmatrix},
\quad U_0=I,\quad U_1=R(\theta_F)\otimes R(\theta_R),
\quad0\le\theta_F,\theta_R\le\pi/2.
\]

Direct reduction gives

\[
V=\cos\theta_F\cos\theta_R,\qquad D_F=\sin\theta_F.
\]

For any admissible (v,d) with d<1 choose
cosθF=√(1−d²), cosθR=v/√(1−d²). The latter lies in [0,1]. For d=1 the region
requires v=0; choose θF=π/2. Thus every point of S_P is realized using the same
resource ceiling and unchanged interface. The associated constant generator
in units ℏ=T=1 is H=P1⊗(θF YF+θR YR); exp(−iH) gives this U, is unitary and
commutes with H0. The direct generator is an algebraic control; the nearest-
neighbour pulse synthesis establishes the declared interaction graph.

An exact fixed-response slice illustrates the freedom:

| V | cosθF | cosθR | D_F | optimal guessing | 1−V²−D_F² |
|---|---|---|---|---|---|
| 3/5 | 1 | 3/5 | 0 | 1/2 | 16/25 |
| 3/5 | 4/5 | 3/4 | 3/5 | 4/5 | 7/25 |
| 3/5 | 3/5 | 1 | 4/5 | 9/10 | 0 |

All three implement the same reduced S channel on **every input state** at the
declared observation time. All conserve the same population labels; all begin
with the same pure memories. Pure initial environmental preparation alone
therefore does not force the local-readout boundary to saturate. This is an
exact conventional witness, not a new impossibility theorem about nature.

## 3. Hostile access control: information in correlations

Start from the last row, |e0⟩=|00⟩, |e1⟩=v|00⟩+√(1−v²)|10⟩. Apply the same
fixed, system-independent environmental unitary W to both branches, whose
columns in the ordered FR computational basis are

\[
W=(|\Phi^+\rangle,|\Psi^+\rangle,i|\Phi^-\rangle,i|\Psi^-\rangle),
\quad|\Phi^\pm\rangle=(|00\rangle\pm|11\rangle)/\sqrt2,
\quad|\Psi^\pm\rangle=(|01\rangle\pm|10\rangle)/\sqrt2.
\]

It is CNOT(F→R)·(HadamardF⊗IR)·(diag(1,i)F⊗IR). Conditional states become
Φ+ and vΦ++i√(1−v²)Φ−. Their two single-qubit marginals are all I/2 because
|v±i√(1−v²)|²=1. Thus D_F=D_R=0 while V=v and D_FR=√(1−v²). The labelled
readout is unchanged. The different microscopic evolution has moved the
classical information into correlations, not destroyed it or changed q.

Any environmental unitary leaves the reduced S density matrix unchanged.
If applied after the final S-environment interaction, with subsequent
decoupling, it also leaves every later S-only intervention/record probability
unchanged. This is not a claim of process equivalence if S subsequently
recouples to the encoded environment. The example supplies no lifetime law
for the stored record; storage isolation is an explicit supplied condition.

## 4. Existing bridges and the precise missing step

Englert's 1996 visibility/which-way inequality already supplies the upper
boundary; renaming D_F “accessible memory” cannot create new physics [1].
For pure conditional states, **complete** access to FR already gives saturation
by the eigenvalue calculation above. Complete access is a different, stronger
premise than the fixed fragment F. Supplying it would buy saturation from
ordinary quantum theory; GRUT has not generated that access premise.

In continuous superconducting measurement, Eddins et al. independently extract
dephasing and information-acquisition rates and discuss unmonitored channels;
their rate convention is ηmeas=Γmeas/(2Γφ) [3]. That ratio is not this packet's
finite-state D_F. Unden et al. prepare and read conditional nuclear-spin
records with system Ramsey measurements; their reported Holevo information
upper-bounds communicated classical information, rather than equalling D_F [2].
These sources establish operational methods, not certification of this exact
three-qubit protocol or a full Stage-3 experimental lock.

The existing GRUT books and working synthesis explicitly retain a supplied
split, access seed, preparation and noise/response structure. They do not
establish a further response/record restriction for this P. Persistent Z labels
here are invariant because the investigator supplied a diagonal interaction.
The partition, physical access, equivalence class and interface freedom are
not generated. Computing an R1 quotient would not repair that failure.

The surviving freedom is **allocation and encoding of response information
among accessible fragments**, rather than an unknown absolute noise scale.
This is a possible operational comparison domain, not an independently
established gravity/matter cross-sector bridge. No GRUT K or narrower allowed
region is proposed or tested. An arbitrary restriction on S_P would merely
postulate the answer; the record presently supplies no derivation of one.

## 5. Reproducible conventional test and falsification scope

From repository root: `bash research/recon02/reproduce.sh`. All fixtures,
runtime versions, exact rational slice controls, matrix errors and input hashes
are in [results.json](results.json) and [MANIFEST.json](MANIFEST.json). The full
8×8 unitary, Hamiltonian exponential, matrix-unit channel, normalized Choi
matrix, conditional partial traces and Helstrom measurements are checked
independently of the displayed scalar formulas. This is numerical evidence
supporting exact algebra; floating-point tolerances are not certified enclosures.

At fixed settings, estimate V via normalized S Ramsey contrast, and determine
ρF|0 and ρF|1 by F tomography after separate label interventions. Calibrate
SPAM, phase drift, population leakage and readout loss. For the three slice
settings above, conventional predictions are V=0.6 in all runs and
D_F=0,0.6,0.8, respectively. The encoding control predicts D_F=D_R=0,
D_FR=0.8; checking that last number requires an explicitly additional joint
tomography diagnostic, not changing F in the main comparison.

If independently certified intervals force V²+D_F²>1, the stated quantum
baseline or its operational premises fail. That is not a GRUT prediction or
evidence for a GRUT selector. For example, V∈[0.79,0.81] and D_F∈[0.69,0.71]
give the exact lower value 1.1002>1, whereas intervals about (0.6,0.8) straddle
the boundary. The script rejects invalid probability intervals and reports
only conservative decisions. No real measured intervals are supplied.

## 6. Verdict

**BASELINE TARGET ESTABLISHED; NO NEW GRUT RESTRICTION; HOLD CARD 2.** Exact
finite-model proofs plus reproducible synthetic matrix controls; author grade,
pending independent review. No candidate information price or scientific
credit is claimed. No new law was formulated/evaluated; no slot or IP-12
candidate assignment was added. Card 1 remains killed on the authorized branch,
with two cards remaining and external scientific review still pending.

## Primary sources

[1] B.-G. Englert, PRL 77, 2154 (1996), DOI
[10.1103/PhysRevLett.77.2154](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.2154).
Publisher abstract inspected; full PDF unavailable. Proof for this packet's
controlled finite model is given above rather than attributed to unseen text.

[2] T. K. Unden et al., PRL 123, 140402 (2019),
[primary manuscript](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=926782).

[3] A. Eddins et al., PRX 9, 011004 (2019),
[primary paper](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.9.011004).

Read scopes and local record references are pinned in [SOURCES.json](SOURCES.json).
