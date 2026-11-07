# GRUT completion review

6 November 2026 · Research synthesis and finite mathematical results · Base commit 6f5fb8f43fe7680c04158585bc961646ab15732d

## The working theory

The strongest current formulation of GRUT is a closed-system research framework. One physical state and one set of dynamics must account for the probabilities of records, interventions, and effective histories inside the system. Wherever different descriptions predict the same experiment, their probabilities must agree to the stated accuracy. The construction uses quantum probabilities and adopts general relativity and the Standard Model as low-energy inputs. Its candidate matter sectors and cosmological preparations are separate hypotheses.

The [existing working theory](../grut_working_theory_v2_20261101/WORKING_THEORY.md) expresses the inputs as \((\mathcal K,\{\hat C_A\},S_{\rm EFT},\vartheta,\mathcal B,\{E_{a|T=t}\})\). The [emergence ledger](../grut_emergence_synthesis_20261105/CLAIMS.md) leaves physical access, state selection, and a distinctive relation open. This review takes the foundational question literally. How much of a subsystem, its interactions, and its effective state can be inferred from the closed whole?

**Result.** A finite commuting stationary datum cannot select a concrete proper quantum subsystem by a deterministic, basis-equivariant rule. With specified local controls, the support of interactions can be recovered exactly from commutators. With specified downstream observables, the smallest sufficient linear state description can also be found exactly in a finite model. These are separate steps. They identify which ingredients a new principle would have to supply and which consequences can already be computed using quantum mechanics.

## What completion would require

| Layer | What the theory has to specify or recover | GRUT status |
|---|---|---|
| Physical states and probability | A physical state space or observable algebra, positive probabilities, and consistent joint predictions for allowed operations. Define the gauge constraints, physical inner product where applicable, and meaning of an internal record. | Quantum structure and the Born rule are adopted. Finite constructions exhibit records relative to supplied access. The physical algebra and state-selection rule remain open. |
| Dynamics and symmetry | A well-posed dynamics or constraint system, with controlled limits and regulator dependence. Recover causal propagation, tested GR, and quantum matter in the applicable regime. A fundamental claim requires a defined domain beyond the low-energy EFT. | Einstein gravity and quantum matter provide the working low-energy theory. Candidate actions and matching rules are inputs. The finite graph program has yet to derive continuum locality or Einstein dynamics. |
| Clock, geometry, and boundary state | Operational clocks and spatial access, justified initial or boundary conditions, and a mechanism for the low-entropy conditions used to generate an arrow. Predictions must specify which choices are calibrated. | Conditional clock and record examples exist. Their clock, access, state, and blank records are supplied. Stationarity and closure leave these choices open. |
| Fields, scales, and vacuum | Explicit field content, symmetries and couplings, controlled quantum corrections, and physical renormalization conditions. A claimed explanation of the vacuum energy must predict the gravitating renormalized value or constrain it independently. | Candidate content and parameters are supplied. The vacuum energy remains a renormalized input. A flat-space energy reference has no power to fix the gravitational term. |
| Effective physical histories | Derive production, feedback, and the observables used at each transition from compatible states and dynamics. Verify that reductions preserve the specific quantities used downstream. | Finite controls and conditional transfer results exist. H55 and the vector reheating realization have different unresolved physical inputs. Their numerical results have distinct scopes. |
| Empirical closure | Fix a prediction for an independent observable or a relation between observables. Declare calibration, uncertainty, comparator, and an outcome that would reject the claim. | GRUT has benchmark and failure results. There is no independently tested prediction that selects a GRUT-specific law over its adopted physics. |

Parameters and boundary conditions may be calibrated in a complete effective theory. Predictivity requires the calibrated construction to constrain additional results. GRUT's stronger claim that the whole forces its parts needs a law selecting or constraining inputs that would otherwise be independent. The program's assumptions of finite total capacity, compact space, or a specific cosmic preparation are conjectures. Internal records are operational correlations. A further claim about subjective experience would require its own definitions and evidence.

Low-energy quantum GR has an established EFT formulation. Internal-clock evolution in a stationary description is also established as a conditional construction. Neither fact determines the physical clock or the GRUT state. [Donoghue](https://arxiv.org/abs/gr-qc/9512024) · [Page and Wootters](https://doi.org/10.1103/PhysRevD.27.2885)

## 1. Stationary input cannot select a concrete quantum factor

Let \(\mathcal H=\mathbb C^D\). The input \(\mathbf M\) is a tuple of pairwise commuting Hermitian matrices. It can include a Hamiltonian \(H\) and a stationary density matrix \(\rho\), or a constraint and the projector onto a vector satisfying it. A rule \(F\) returns a concrete unital embedded algebra and is equivariant under unitary change of basis,

\[
F(U\mathbf M U^\dagger)=UF(\mathbf M)U^\dagger.
\]

A proper quantum factor here means an algebra unitarily equivalent to \(M_a(\mathbb C)\otimes I_b\) with \(a,b>1\) and \(ab=D\). The multiplicity \(b\) expresses the presence of a complementary subsystem.

**Proposition.** No such rule can select a proper quantum factor from this input. The same holds if the output is a finite unordered collection of concrete proper quantum factors.

**Proof.** Simultaneously diagonalize \(\mathbf M\). Every diagonal phase rotation \(D(\theta)\) fixes the input, so the selected algebra \(\mathcal A\) must be invariant under conjugation by this torus. If every element of \(\mathcal A\) is diagonal then \(\mathcal A\) is commutative. Otherwise some \(A_{ij}\ne 0\) with \(i\ne j\). The phase-weighted average of \(D(\theta) A D(\theta)^\dagger\) with the factor \(\exp[-i(\theta_i-\theta_j)]\) isolates \(A_{ij}|i\rangle\langle j|\) and belongs to \(\mathcal A\). Adjoint and multiplication give the rank-one projector \(|i\rangle\langle i|\) in \(\mathcal A\). Every nonzero projector of \(M_a\otimes I_b\) has rank divisible by \(b>1\). Both cases contradict the proposed output. For a finite collection, the continuous torus orbit of each algebra lies in a finite set and hence is a single point. The same argument applies to each algebra. \(\square\)

This closes a specific route to the access problem in the [closure-forcing audit](../grut_closure_forcing_20261027/SYNTHESIS.md). Covariance alone on stationary commuting inputs supplies neither a concrete clock factor nor a concrete spatial tensor factor. Conditioning on a clock with its algebra already specified is consistent with the proposition. The proof allows outputs such as a family with a continuous orbit, an equivalence class, or an algebra with a center. It also leaves noncommuting physical inputs, gauge-invariant algebra nets, and the continuum outside its scope. Generic spectral locality arguments infer a tensor structure up to dualities under a locality prior. Those inferences concern a different output and use a restriction on the admissible descriptions. Gauge theories can instead require regional algebras with centers. [Cotler, Penington and Ranard](https://arxiv.org/abs/1702.06142) · [Casini, Huerta and Rosabal](https://arxiv.org/abs/1312.1183)

## 2. Operational input recovers an interaction graph

Supply an unordered set of Hermitian involutions. Assume each has exactly one anticommuting partner and commutes with all other controls. Each pair \(X_i,Z_i\) generates \(M_2\) with \(Y_i=iX_iZ_i\). The different pairs commute. Their joint multiplication represents \(M_{2^L}\) for \(L\) pairs. The representation is unital and faithful by simplicity. If \(D=2^L\), it has multiplicity one and the recovered algebras form a full tensor decomposition up to ordering. Otherwise, when a representation exists, its multiplicity must be accounted for separately. The noncommutation graph of the controls determines the grouping. The input controls encode this operational structure, as in observable-induced tensor decompositions. [Zanardi, Lidar and Lloyd](https://arxiv.org/abs/quant-ph/0308043)

There is an exact test for interaction supports once the factors are known. Let \(\|A\|_{2,n}^2=\operatorname{Tr}(A^\dagger A)/D\), and expand \(H=\sum_P h_P P\) in orthonormal Pauli strings. For a set of distinct sites \(S=\{i_1,\ldots,i_s\}\), define

\[
W_S^2=\frac{1}{8^s}\sum_{\mu_1,\ldots,\mu_s\in\{x,y,z\}}
\left\|\left[\cdots[H,\sigma_{i_1}^{\mu_1}],\ldots,\sigma_{i_s}^{\mu_s}\right]\right\|_{2,n}^2
=\sum_{P:\ S\subseteq\operatorname{supp}P}|h_P|^2.
\]

For fixed \(\vec\mu\), a Pauli string survives precisely when its factor at each selected site is a nonidentity Pauli different from the chosen \(\sigma^\mu\). The map on surviving strings is injective up to signs and the factor \(2i\) per site. Orthogonality removes cross terms, and each site gives two surviving choices and a squared factor of four. This yields the stated normalization and exact formula. The sum is invariant under orthonormal changes of local Pauli basis. Descending subtraction of contributions from strict supersets gives the weight at each exact support. A pair weight alone includes every many-body interaction containing that pair.

There is an operational readout. If the sender applies \(\exp(-i\epsilon B_j)\) at time zero and the receiver reads \(A_i\) at time \(t\), the first-order kick response is

\[
\chi_{ij}(t)=\left.\partial_\epsilon\langle A_i(t)\rangle_\epsilon\right|_0
=i\operatorname{Tr}\big(\rho[B_j,A_i(t)]\big),
\qquad A_i(t)=e^{iHt}A_i e^{-iHt}.
\]

Its expectation depends on the chosen state and readout. Across the three Pauli axes at both sites, define the operator response
\(R_{ij}(t)=(\sum_{\alpha,\beta}\|[\sigma_i^\alpha(t),\sigma_j^\beta]\|_{2,n}^2)^{1/2}\).
For distinct sites its first-order coefficient is \(8W_{\{i,j\}}\). This follows from the expansion of \(A_i(t)\) and the two-site identity above. Draw an edge wherever a support contains both vertices. An \(n\)-fold commutator with \(H\) can extend the operator support at most \(n\) steps in this graph. At distance \(d\), all derivatives of the commutator before order \(d\) vanish. Particular interactions and states can delay or eliminate the response. The test concerns early-time operator support in a specified finite Hamiltonian. Size-uniform propagation bounds require locality and norm hypotheses. [Nachtergaele and Sims](https://arxiv.org/abs/1004.2086)

### Finite replay

In a five-qubit chain with XY bonds \(J=(0.31,0.43,0.37,0.29)\) and \(Z\) fields \((0.13,-0.19,0.07,0.11,-0.17)\), jointly rotate the Hamiltonian and controls by a seeded complex unitary and shuffle the controls. The evaluator groups controls using their matrices alone. It consults the planted labels only to score the answer. It recovers all five factors and the four bond weights \(\sqrt{2}|J_i|\) with maximum residual \(8.2\times10^{-16}\). The endpoint time-series slopes use eigendecomposition, separate from computing the nested commutators.

| Setup | Inclusive endpoint weight | First nonzero endpoint derivative order | Finite-time log slope |
|---|---:|---:|---:|
| Nearest-neighbor chain | Roundoff | 4 | 3.99979 |
| Add \(0.11Z_0Z_4\) | 0.110000 | 1 | 0.99980 |
| Add \(0.09Z_0X_2Y_4\) | 0.090000 | 1 | 0.99968 |
| Keep \(H\) and its Gibbs state fixed and use CNOT-conjugated endpoint access | 0.552449 | 1 | 0.99931 |

For the three-site control, support subtraction gives squared weight \(0.0081\) on \(\{0,2,4\}\) and zero on the exact endpoint pair. A deliberately inserted joint control in place of a primitive fails the grouping test. A separate diagonal-phase symmetry of \(H,\rho\) preserves them to roundoff while moving a norm-one site operator a distance \(0.891\) outside its original local algebra in normalized Hilbert–Schmidt norm. The two fixed-\(H\) checks address different questions. The phase symmetry changes the concrete algebra within the input's stabilizer. The CNOT changes operational access and with it the graph measured at fixed \(H,\rho\).

One measurable comparison uses \(B=X_4\), \(A=X_0\), and prepares the endpoints in \(+Y\) with mixed intermediate qubits. The linear-in-time coefficient of \(\chi(t)\) is zero for the chain. With the direct link it is \(-4(0.11)=-0.44\) in the simulation's units with \(\hbar=1\), since the direct term gives \(-4(0.11)Y_0Y_4\) in the derivative of the kick response. The matrix replay gives the same coefficients. This could discriminate the two Hamiltonians on a calibrated simulator. It is a prediction of their specified ordinary quantum dynamics. The example uses supplied time, preparation, and controls in a finite graph. A physical clock and a relativistic metric are separate requirements.

## 3. The exact finite criterion for an effective state

Fix a finite Hamiltonian and the Hermitian readouts \(\mathcal O\) to be predicted at arbitrary later times. In a model with unitary evolution and no further interventions, put \(\mathcal L(O)=i[H,O]\). Define the real operator space

\[
\mathcal V=\operatorname{span}_{\mathbb R}\{I,\mathcal L^n(O):O\in\mathcal O,\ n\ge0\}.
\]

It can be computed by closing the readout span under \(\mathcal L\) until it stops growing. The expectations of a basis of \(\mathcal V\) suffice to predict all the stipulated future expectations. Their dynamics close linearly because \(\mathcal L\mathcal V\subseteq\mathcal V\). This is the smallest such span among summaries by fixed linear expectations valid for every initial density matrix. To see necessity, the span of the evolved readouts and \(I\) is \(\mathcal V\) by the finite-dimensional time expansion and its derivatives. If the summary misses a direction in that span, one can choose a traceless Hermitian \(D_0\) orthogonal to all summary operators but with nonzero inner product with an evolved readout. Both \(I/D+\epsilon D_0\) and \(I/D-\epsilon D_0\) are states for small enough \(\epsilon\). They agree on the summary and disagree on that readout. The trace is fixed, leaving \(\dim\mathcal V-1\) independent real expectations.

For an allowed sequence of quantum instruments, the corresponding branch probabilities are expectations of effects obtained by applying the adjoints of the branch maps in reverse time order. Their span gives the same criterion for that declared family of experiments. Restricting the class of allowed preparations can also reduce what is required. The criterion gives neither a guaranteed small state space nor a particle interpretation.

**Example.** For \(H=X\otimes X\) and both qubits' single-qubit readouts, six Bloch components describe the two initial marginals. Exact prediction of both future marginals for arbitrary preparations requires ten expectations. Four additional ones correspond to \(Y\otimes X,\ Z\otimes X,\ X\otimes Y,\ X\otimes Z\). The local operators and these four span the nonconstant directions under commutation with \(H\). Both states
\(\rho_\pm=(I\pm Y\otimes X)/4\) have identical initial marginals. The first-qubit trace distance after evolution is \(|\sin(2t)|\), giving \(0.707107\) at \(t=\pi/8\) and \(1\) at \(t=\pi/4\). A common predicted first-qubit state must incur error at least half their separation in one of the preparations, by the triangle inequality. The replay checks the rank of the predictive span and the trace distances.

This gives a general handoff test. With the candidate retained state fixed, use two physically allowed matched preparations and evolve them with the parent dynamics. If an intended downstream observable separates by \(\Delta\), a single prediction from the retained state alone makes an error of at least \(\Delta/2\) on one of the pair. This turns a disagreement into a lower bound on the information lost by the proposed reduction. A usable tolerance has to be set for the actual downstream observable. The finite example does not identify the needed moments of the continuum H55 state.

## What to carry forward

| Research lane | Defensible result and decision |
|---|---|
| H55 and sterile transport | Preserve homogeneous and propagating-source controls and the evidence that projected occupations can lose work, pressure, or channel memory. Feedback and abundance remain suspended pending a controlled physical source and handoff. The historical 0/135 ledger applies to the separate higher-scale M2-X acceptance set. [Research reset](../grut_research_reset_20261106/AUDIT.md) |
| ALP cogenesis and minimal thermal scalar | The bounded ALP conversion calculation closed as unresolved, so there is no derived physical late conversion factor from this campaign. The reported 40 GeV minimal scalar exclusion is conditional on that frozen portal, freeze-out history, and all-dark-matter assumption. [Comparator records](../grut_research_reset_20261106/) |
| Inflationary 1 meV vector | Preserve the replicated 4.95-fold transfer enhancement at matched primordial spectrum for the specified backgrounds. The ordinary causal gas realization of the logistic conversion failed the seed-population test. The \(r_{0.05}=0.002765\) numerical fit retains assumed scalar calibration. [Vector claim ledger](../grut_research_reset_20261106/vector_audit/CLAIM_LEDGER.md) |
| Foundations | Preserve operational records and the tests of spatial transport. Close deterministic concrete factor selection from commuting stationary data within the theorem's scope. Use the response and sufficiency calculations as finite criteria for claims of geometry or effective state emergence. |

The first foundational decision is now explicit. GRUT can **supply** physical access and aim for a predictive effective theory, or **derive** it from additional physical structure and test that extra premise. A claim of derivation must declare the new input, explain which alternatives it rules out, and show that the target access or interaction graph was absent from the inputs. The simple paired-control example has its access structure encoded in the control commutators and therefore serves as a benchmark. The same-\(H,\rho\) example is a rejection control for rules that infer a concrete response graph from stationary data alone. A cosmological model can separately qualify as a physical theory of a mechanism when it earns its own preparation and evolution and passes a discriminating test. It need not settle the foundations.

The new results here are exact finite-model statements with a numerical replay of the examples. The observable-induced decomposition and the support argument use established quantum mechanics. Literature priority for the specific propositions has not been assessed exhaustively. The calculation has not received independent external review or been tested on hardware. There is presently no derived GRUT-specific law selecting the observed universe.

## Replay

Run the following from this directory with Python and NumPy installed.

```sh
OPENBLAS_NUM_THREADS=1 python -I replay.py
```

The reference [output](results.json) records the seed, thresholds, residuals and time samples. Near-zero residuals may vary with numerical libraries. The assertions check recovery, the rejected control, three distinct interaction supports, the fixed-input access comparisons, the scalar kick slopes, and the predictive dimension and trace-distance tests.
