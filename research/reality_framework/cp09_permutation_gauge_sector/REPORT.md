# CP-09 permutation-invariant quantum sector: a constructive mode-identity correction

**Date:** 2026-10-11. Independent, additive research construction on \`codex/grut-cp09-permutation-gauge-sector\`, forked from [CP-09 mode identity e135cef](https://github.com/ryangrvr/TestingGRUT/blob/e135cef7c60b27cbc15f0f1e49f83d8f5c9c838e/research/reality_framework/cp09_mode_identity/REPORT.md), using the [independently assembled N=5 full-poset Hamiltonian](https://github.com/ryangrvr/TestingGRUT/blob/476019ef2f0664a7e6d061f29fafc33106770e18/research/reality_framework/cp09_five_event_geometry/check.py). Frozen CP-09 at 1dadbd6, earlier constructions, canonical GRUT-RAI and Stage-3 Card 3 are UNCHANGED. This does not open CP-10.

## Central physical choice: permutation symmetry versus gauge identification

The CP-09 Hilbert space is spanned by **labelled** causal configurations \(|C\rangle\) tensor N labelled matter qubits. Its frozen Hamiltonian is covariant under *simultaneous* permutation of geometry and matter labels:
\[
[H_{09},U(\pi)]=0,\qquad \pi\in S_N.
\]
Covariance is an established exact property of the supplied finite law. It does not, by itself, say whether event permutations are (a) a global symmetry that can carry distinct physical charge sectors, or (b) a redundancy of event labels that should be quotiented as a gauge.

**We explicitly test option (b) as a declared additional physical interpretation**, rather than quietly treating it as derived. Define
\[
P_{\mathrm{inv}}=\frac{1}{N!}\sum_{\pi\in S_N}U(\pi),\quad
\mathcal H_{\mathrm{inv}}=P_{\mathrm{inv}}\mathcal H,\quad
H_{\mathrm{inv}}=H_{09}\big|_{\mathcal H_{\mathrm{inv}}}.
\]
Since the group commutes with H, \(\mathcal H_{\mathrm{inv}}\) is an exactly closed unitary sector without any Hamiltonian modification, penalty, new coupling, finite-neighbor cap or loss of causal histories. This is a **construction of a possible physical-state space**, not a derivation that Nature gauges event relabelling. Other gauge treatments that include reference frames and nontrivial sectors would require their own physical postulates.

## Constructive computation: exact orbit basis through N=5

Construct simultaneous permutation orbits of every computational basis pair \((C,b)\), where b is an N-qubit matter bit string. For each orbit \(\mathcal O\), use the normalized orbit sum
\[
|\mathcal O\rangle=\frac{1}{\sqrt{|\mathcal O|}}\sum_{(C,b)\in\mathcal O}|C,b\rangle .
\]
The columns of isometry T are precisely these orthonormal orbit vectors. Then
\[
T^\dagger T=I,\quad
H_{\rm inv}=T^\dagger H_{09}T,\quad
H_{09}T=T H_{\rm inv}.
\]
All three identities are explicitly checked in finite arithmetic; maximum subspace-intertwining residual at N=5 is about \(4.4\times10^{-13}\). The original 135,392-dimensional five-event system becomes a **1,490-dimensional exactly invariant subsystem**, with NO truncation of allowed causal orders. Note this does not reduce the full physical Hilbert space unless the gauge interpretation is supplied.

| Events N | Total labelled posets | Full matter×geometry dimension | Simultaneously invariant dimension | Full first gap | Invariant-sector first gap |
|---:|---:|---:|---:|---:|---:|
| 2 | 3 | 12 | 7 | 0.318298542 | 0.388191625 |
| 3 | 19 | 152 | 32 | 0.185712205 | 0.288781726 |
| 4 | 219 | 3504 | 192 | 0.096503864 | 0.209244875 |
| 5 | 4231 | 135392 | 1490 | 0.041881530 | **0.149418632** |

The ground energies remain the same as the original full model: 0.094755284, 0.310979227, 0.603757645, 0.953282910, respectively. That is because the finite-N ground state is already permutation invariant, numerically and analytically under the nondegeneracy assumption. Projecting onto this proposed physical sector **does not change the ground-state geometry distribution**: it cannot, by itself, establish spacetime selection. The invariant gap also decreases over these four sizes. No continuum limit, critical exponent or finite excitation mass has been derived.

N=5 full eigensolver's first gap is the published 0.04188153 and is independently reproduced; as earlier reported, a small \`eigsh\` call can miss the complete degenerate full shell. We do **not** certify its full N=5 representation multiplicity from that truncated call. The fully invariant gap at N=5 is separately obtained by diagonalizing the exactly closed invariant sector.

## Why the first low-energy shell was invisible to collective scalar probes

For N=3, the first full excited shell is two-dimensional; for N=4 it is three-dimensional. After carefully orthonormalizing the numerically returned degenerate shell, the action of a transposition has trace 0 at N=3 and 1 at N=4, as required for the **standard representation** of S_N (dimension N-1, transposition character N-3). The shell has zero projection into the invariant space. This is an **actual representation-theoretic identification for N=3,4**, stronger than the original "consistent with a nontrivial sector" statement. For N=5 the full shell is not exhaustively resolved here.

If the ground state |g> is invariant and an observable O is invariant, an excited state |e> in a nontrivial representation obeys
\[
\langle e|O|g\rangle
=\langle e|P_{\mathrm{inv}}O|g\rangle=0.
\]
Therefore total causal-relation number and collective mean matter X must have vanishing first-shell transition amplitudes by a symmetry selection rule. Their vanishing is **not evidence that this entire shell carries no causal-geometry information**.

To test it directly, independently predeclare the *event-specific* geometry-vector operator
\[
D_i(C)=\frac{n_i^{\mathrm{out}}(C)-n_i^{\mathrm{in}}(C)}{N-1}.
\]
The normalized summed spectral transition weights into the original lowest shell are:
\[
|\langle \mathrm{first\ shell}|D_0|g\rangle|^2=
0.122740428\ (N=3),\qquad
0.110641539\ (N=4).
\]
The nonisolated-minimal-role indicator of event 0 likewise couples with weights 0.053913724 and 0.059568092. These operators transform with an event label, so they are **not gauge invariant by themselves**; a physical event-specific measurement would require a relational identifying reference. The result shows that the low shell is sensitive to local causal role while still not demonstrating spatial propagation or physical localization.

### Numerical correction to the original mode-identity note

The earlier [e135cef shell-weight script](https://github.com/ryangrvr/TestingGRUT/blob/e135cef7c60b27cbc15f0f1e49f83d8f5c9c838e/research/reality_framework/cp09_mode_identity/check.py) summed \( |V^\dagger O g|^2 \) across raw nearly degenerate \`scipy.sparse.linalg.eigsh\` eigenvectors. For complex sparse Hermitian matrices SciPy's algorithm can return vectors **not mutually orthogonal within a degenerate numerical shell**. Summing raw squares is then not the squared norm of an orthogonal eigenspace projector; the values can change with how many modes were requested.

Use QR to form an orthonormal shell basis Q and compute \( \|Q^\dagger O g\|^2\) instead. Independent runs with 10, 12, 16, 20 requested eigenvectors give stable corrected \`local_X0\` first-shell weights:
\[
N=3:\quad\boxed{0.002600147453},\qquad
N=4:\quad\boxed{0.004987763689}.
\]
The original published numbers 0.00261678656 and 0.00386211812 are **not reliable shell-projector weights**; do not cite them as independently stable values. The original zero scalar couplings remain true by the exact symmetry argument. Frozen previous files are deliberately not overwritten. The new test includes this correction as an assertion.

## An actually invariant collective geometry–matter observable responds

As a physical observable candidate predeclared independently of the eigenmodes, take
\[
\boxed{\mathcal G_N=\frac{1}{N(N-1)}
\sum_i[n_i^{\mathrm{out}}(C)-n_i^{\mathrm{in}}(C)]Z_i.}
\]
It is diagonal in geometry/computational matter basis and **invariant under simultaneous permutations**. It measures the correlation between a site's causal past–future orientation and its matter Z-polarization. It is a global collective operator (requiring geometry-dependent access), not a localized detector, metric, wavevector, gravitational field or derived stress tensor.

This observable has a **nonzero** transition to the first invariant-sector excitation. Its first-mode weight and finite-system stationary retarded self-response
\[
R_{\mathcal G\mathcal G}(t)=i\theta(t)\langle g|[\mathcal G_N(t),\mathcal G_N(0)]|g\rangle
\]
are obtained from the same reduced H with no fitted kernels:

| N | First invariant-mode squared matrix element | R_GG(0.4) |
|---:|---:|---:|
| 2 | 0.014207570 | 0.377611372 |
| 3 | 0.016808506 | 0.133426761 |
| 4 | 0.016461917 | 0.070987432 |
| 5 | 0.014221339 | 0.043643907 |

Independent source kicks \(\exp(\pm i\,10^{-5}\mathcal G_N)\) followed by full invariant-sector unitary evolution reproduce each Kubo response within \(10^{-8}\). The decreasing values partly reflect the explicit 1/[N(N-1)] normalization; interpreting that trend as a propagation law or universal scaling exponent would be unjustified.

The calculation **constructs a symmetry-respecting joint geometry–matter quantum observable with actual finite-N dynamics**, and distinguishes it from low modes associated with event-label representation. It still does **not** produce physical separation, light cones, velocity, a continuum dispersion law, spacetime-class preference or Einstein dynamics.

## The next physical hypothesis to test—not to assume as earned

Treat labels as gauge only if the theory supplies a consistent relational-reference interpretation. A useful next construction would **couple two individually identifiable matter reference systems to geometry**, define their separation entirely through gauge-invariant conditional correlations, and then calculate their response and the *jointly generated geometry distribution*. Those references must have an explicit apparatus Hamiltonian and cannot use a global role oracle as an unpriced free resource. The same physical criterion must be applied to both manifoldlike and layered competitors.

Do not open a new candidate solely because a quotient gap is smaller, larger or decreasing. The genuinely missing physics is a law connecting *reference interactions, physical separation and causal propagation* while selecting the background structure from the very same dynamics. Full experimental gravity and massless spin-2 degrees are still open.

## Reproduce and verification scope

Run from any working directory:
\`\`\`bash
python research/reality_framework/cp09_permutation_gauge_sector/check.py
\`\`\`
Requires NumPy and SciPy. The executable reuses the independently assembled full Hamiltonian from \`cp09_five_event_geometry/check.py\`, with a **separate** orbit/isometry/gauge-reduction algorithm, verifies unitary-closure of invariant sectors, independently computes full and invariant spectral gaps, corrects the nonorthogonal degenerate shell, evaluates the causal alignment observable and compares its retarded response with a source-impulse derivative. Components and equivalent calculations were run in the analysis environment. **The exact committed script itself is not yet separately certified verbatim and remains for partner execution.** Sparse eigensolvers and numerical residuals are not proofs of asymptotic completeness.

Scientific status: an exact finite quantum symmetry-sector construction with a tested relationally meaningful *collective* response, and a correction to spectral-mode identification. No CP-10, no externally confirmed physical law, and no claim of gauge principle uniqueness.
