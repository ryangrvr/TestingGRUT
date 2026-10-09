# Padgett / Quantum Information Holography audit

**DRAFT — external review pending. Cutoff: 2026-10-09 UTC.** This is an
equation-level audit of five retrieved PDFs, including the 149-page April
2026 revision and the May 2026 manuscript. It is not a judgment based only on
the early two-page proposal. Source identifiers and file hashes are in
`SOURCES.json`; exact controls and scope qualifications are in `PROOFS.md`.
The October 2023 asymmetric-primes record Q06 was found but its PDF was not
retrieved; no body-level verdict is assigned to that document.

## What the strongest available version supplies

The later presentation supplies a complex qubit description, a two-channel
radial/angular operator family, a proposed screen encoding, a metric ansatz,
and spectral/counting prescriptions. This is more concrete than a verbal
claim that information produces gravity. A defensible reduced interpretation
would be a **standard quantum/open-system model with a stipulated code and
calibrated geometric observables**. That could be studied once the domain,
coefficients, state preparation, boundary data, and measurement channel are
specified. It would initially be an effective model, not a derivation of all
physics or a demonstration of a unique microscopic selector.

The April text itself acknowledges the standard status of the qubit Born
probabilities. Its valid trigonometry and operator identities should be
retained. The question is what extra, internally consistent restriction the
assembled proposal supplies. The audit found no reproducible quantitative
prediction distinguishing that assembled proposal from its supplied QM/GR
components. This is a finding about the inspected corpus, not proof that no
future QIH model can be completed.

## Equations and obligations

Rows describe mathematical content; **IDENTITY** does not mean new physics.
**POSTULATE** may define a possible model, but does not establish its physical
selection. **GAP** means the stated implication has not been demonstrated.

| Item / source location | Classification | Decisive check or missing obligation |
|---|---|---|
| `E=hf=ħomega`, mass-energy identification (Q01) | Standard identity with scope | Photon energy is not a nonzero photon rest mass; no gravity dynamics follows. |
| Reciprocal angular projection (Q01) | Failed stated identity | C02 gives 5/3 for real unit vectors. |
| Time/velocity/potential from phase derivatives (Q02) | Definitions / unsupported physical identifications | Need dispersion, scales, coordinate convention and field equation. SI velocity and potential dimensions fail as written (C02). |
| Qubit amplitudes and Z probabilities (Q03 opening) | STANDARD QM | Hilbert space, basis and Born measurement rule are already supplied; phi matters for other readouts (C01). |
| Planar velocity and its derivative (Q03 opening) | POSTULATE + correct kinematic identity | Bloch-to-space identification and theta dynamics are additional inputs. Later derivative has correct acceleration units. |
| Phase ticks and derivative conversion (Q03) | Reparameterization | A clock variable and frequency conversion do not select dynamics or solve absence of a preferred clock. |
| Bulk/screen tensor product and `U*U=I` (Q03 §15) | POSTULATE + capacity obligation | Isometry requires source dimension no larger than target; code restriction is needed for a finite screen (C05). |
| Horizon Hilbert space and block H (Q03 §117) | Model family | Complex QM, radial domain, sphere, two channels and operator coefficients are supplied. Closed-domain block construction can be self-adjoint; coefficients still do not follow from it. |
| Vanishing boundary current / self-adjointness (Q03 §144) | GAP | A zero-current statement does not replace a specified maximal self-adjoint operator domain. Deficiency/domain conditions must actually be checked. |
| Functional calculus and unitary evolution (Q03 §129) | STANDARD RESULT | Conditional on a valid self-adjoint H; does not select H. |
| Gibbs trace (Q03 §129) | Missing existence hypothesis | Self-adjointness alone does not imply trace-class `exp(-beta H)`; chiral unbounded examples fail (C04). |
| Prime trace/potential prescriptions (Q03 §§155–157) | Target supplied / inverse spectral gap | Requiring prime periods and coefficients to match zeta inserts the desired target. It does not prove the unique operator or orbit dynamics. |
| Plus-sign spectral determinant target (Q03 §171, PDF pp.140–141) | **Exact obstruction in stated class** | Q03 uses exponent -z/2. The normalized finite determinant is ≥1 for real omega, whereas the xi target has real zeros (C03). Infinite claim needs a valid regularized construction with the stated scope addressed. |
| Metric from phase/state derivatives (Q03) | POSTULATE | Signature, covariance of supplied fields, causal structure and coefficient units need proof; naming the tensor g does not prove these properties. |
| Geodesic variation / Einstein-form equation (Q03) | Standard conditional result / imposed matching | Geodesic motion follows from a supplied metric action. Defining curvature does not derive Einstein dynamics, matter content or a stress tensor. |
| EM mode filter and alpha counting recipe (Q03 §§159–160) | Proposed prescription, not calculated prediction | Projector, coherence measures, thresholds, spectrum and maximization need fixed operational definitions and a unique solution with uncertainty. “No tuning” does not follow from writing a maximum. |
| Spectral counting and gauge beta functions (Q03) | Parameterization / standard-theory matching | Counts and weights required to reproduce known field-theory coefficients are supplied targets; SM representations and independent couplings are not derived. |
| Paired frequency sum / vacuum residual (Q03) | GAP | Exact block grading pairs ± frequencies; signed finite trace vanishes. A physical nonzero vacuum stress needs an explicit, justified breaking/observable/regularization (C04). |
| Modulo-30 prime/dark classification (Q05) | Arithmetic + unsupported physical inference | Four primes multiply to 210, not 30. Residue counts are not cosmic mass fractions; composite 49 belongs to a permitted residue class (C12). |
| Symmetry/asymmetry protects coherence (Q05) | Noise-dependent claim | X-dephasing preserves symmetric plus state and mixes biased zero state; a universal ordering fails (C12). A specified coupling can change that ordering. |
| Unit-circle velocity then `v=omega` (Q04 §2) | Failed simultaneous equations | Norm is always 1, while omega=1/2 demands norm 1/2. Distinct variables need a replacement map (C13). |
| Microtubule Hilbert factors, coupling and switching (Q03 §137, Q04) | Supplied biological model / untested bridge | Finite qubits, local/environment Hamiltonians, thresholds and gravitational collapse times are premises, not neural mechanisms derived from H. |

## Dynamical, dimensional, and covariance assessment

There is an operator **family**, not a completely fixed action/Hamiltonian
and boundary-value problem that yields all claimed sectors. Constraints on
high-frequency speed, symmetry and a large-radius potential do not determine
all coefficient functions and higher-order terms. The latest half-line
operator (§169) also requires its weighted domain and asymptotic data.
Compressing an operator with P does not automatically preserve
self-adjointness unless the domain and invariance/reduction conditions are
controlled. No calculation here assumes those unspecified conditions.

The valid geometric qubit expressions are representation-dependent until
state, measurement and spatial frame are related by a physical channel.
Quantum reference-frame covariance is not recovered merely by writing a
scalar phase or a metric-shaped expression. Conversely, an SI dimensional
objection must not be applied unchanged to a consistently defined natural-
unit formulation. C13 remains a contradiction even at c=1; C03 is an
algebraic sign/zero obstruction independent of unit convention.

No assembled derivation of the SM spectrum, interacting field content,
quantitative Newton/Einstein limit, measured dimensionless couplings, or an
independently selective new observable was established in these PDFs. Many
of their target equations would be useful **checks of a future completed
model**, but imposing them as matching conditions is not generating them.

## Peer review, citation, and experimental status

Zenodo deposition and a deposit's “journal article” resource type do not
establish journal peer review. The searches did find an independent 2026
Zenodo critical/propositional comparison (Q07); its landing page later
returned 410, so only retrieved metadata is used. This is not an independent
replication. A conference abstract and commercial claims are also different
evidence categories. No independent experimental QIH validation or
reproducible gravity/constants benchmark was found in the bounded search
record. The company's retrieved live page said “Coming Soon”; cached
engineering/AI claims do not substitute for published test protocols. [Q07,Q08]

The 2013 Padgett-associated neuroscience case report concerns mathematical
imagery and neural activity, not a measurement of QIH horizon physics.
Subjective experience, savant abilities, and proposed biological coherence
are not evidence that the determinant, field dynamics or constants are
correct. This assessment does not depend on any judgment about the author's
abilities or experiences. [N01]

## Hostile evaluator conclusion

**Do not promote QIH to a GRUT law or use it to fill unspecified operators.**
Retain valid QM/kinematic identities as controls, and the five precise
obstruction families as audit fixtures: angular/unit/velocity errors,
encoding capacity, Gibbs existence, determinant zeros, and noise-basis
dependence. The determinant target fails the stated finite class. The
broader unification claims remain **NOT ESTABLISHED**, not a theorem that
every corrected information-based theory is impossible. Any later revision
must first supply a fully defined model and independent predictions; matching
well-known equations by relabelling them will not meet that obligation.
