---
title: "GRUT: A Working Theory of Emergence and Interaction"
subtitle: "An evidence-bounded synthesis"
author: "Project synthesis by Codex for Ryan Grover"
date: "5 November 2026 UTC"
---

**The central proposal:** reality can be described as one interacting physical whole whose parts acquire persistent identities, exchange energy, preserve records and, in some regimes, use information to influence what happens next. The whole and its parts are compatible levels of physical description. Their compatibility does not yet explain why these laws, this initial state, this vacuum density or conscious experience exist.

This is the most complete account the inspected GRUT record presently supports. It combines established physics, conditional model calculations and explicit open questions. It does not announce a finished fundamental theory, a successful cosmological parameter point, or evidence of a conscious universe.

Read the [emergence map](EMERGENCE_MAP.svg), [claim ledger](CLAIMS.md), [sources](SOURCES.md) and [audit](AUDIT.md) alongside this account. The repository snapshot is **8ab5829a69f51cf462c56b89252d78fa21760300**. The latest numerical verdict takes precedence over older model-completion language.

**Status vocabulary.** ADOPTED means a framework premise, calibration or model choice. STANDARD means established physics within its usual scope. CONDITIONAL means a derivation or calculation from named inputs. FAILED-IN-SCOPE means a specified test failed; it does not exclude every extension. OPEN means a missing result. PROPOSED means a research direction. These labels have different jobs; they are not a ladder that automatically promotes assumptions into facts.

## 1. What is being completed

There are three different finish lines:

1. **A working account:** a coherent description of the whole, its parts, their development and the evidence for each connection. This document supplies that account.
2. **A working model instance:** a fully specified action, state, parameters and observation rule whose connected evolution meets its declared empirical and numerical requirements. GRUT has conditional instances and successful calculations, but the current M2-X instance does not meet its registered acceptance rule.
3. **A fundamental explanation:** a derivation of quantum structure, gravity, state selection, constants, vacuum density and the relation to experience. GRUT has not completed this.

The account remains valuable when a candidate fails: it identifies which connection was tested and which input the failure exposes. RAI also requires scrutiny of negative claims. A small search does not become a universal impossibility theorem, and an evaluator that rejects known successes is not evidence against nature. Conversely, surviving a numerical test is not evidence that an unobserved field exists. The resurrection and dialectic records make both errors concrete. [R01–R03]

The legacy responsive-vacuum mechanism is not revived by calling the present synthesis GRUT. Its recorded failures retain their stated scope. The broader program name now organizes an explicit framework and separate model candidates.

## 2. The theory specification

Use the formal specification

$$
\mathfrak T=
\left(\mathcal K,\{\widehat C_A\},S_{\mathrm{EFT}},
\vartheta,\mathcal B,\{E_{a\mid T=t}\}\right).
$$

| Slot | Meaning | Present status |
|---|---|---|
| \(\mathcal K\) | Kinematics: degrees of freedom, state space or observable algebra | ADOPTED quantum structure and selected field content |
| \(\widehat C_A\) | Gauge, spatial and gravitational constraints | Standard constraint requirements; no complete GRUT quantum-gravity construction |
| \(S_{\mathrm{EFT}}\) | Dynamics in the regime being described | Standard Model, gravity as an effective description, and explicitly added candidate sectors |
| \(\vartheta\) | Couplings, masses, matching scale and counterterms | Chosen or calibrated, with some conditional relations |
| \(\mathcal B\) | State/preparation rule | Supplied cosmological state class and model-specific preparation |
| \(E_{a\mid T=t}\) | Operational observation effects, conditioned on a clock | Supplied measurement dictionary; finite examples constructed |

This tuple is a checklist for a theory, not an equation that uniquely selects one. Naming constraints does not solve their algebra or construct their physical inner product.

At accessible energies a schematic action is

$$
S_{\mathrm{EFT}}=\int d^4x\sqrt{-g}\left[
\frac{M_P^2}{2}R-\rho_{\mathrm{vac}}^{\mathrm{ren}}
+\mathcal L_{\mathrm{SM}}+\mathcal L_{\mathrm{added}}
+\sum_{d_i>4}\frac{c_i\mathcal O_i}{\Lambda^{d_i-4}}
\right]+S_{\mathrm{boundary}}.
$$

Here \(M_P\) is the reduced Planck mass and \(\Lambda\) is an EFT cutoff, not the cosmological constant. The split between the written vacuum term and field potentials must use a consistent subtraction convention. All symmetry-allowed relevant and marginal operators need matching conditions; a tree-level zero is not an all-orders identity. Gravity can be used as a low-energy EFT without claiming its ultraviolet completion. [S01]

**ADOPTED:** quantum probabilities, the effective causal structure and gravitational dynamics used by each model. “No outside” by itself derives none of the Born rule, unitarity, global purity, compact topology or a finite state count. A closed quantum system can have a mixed state; purity is an additional choice. The closure campaigns found no additional GRUT law in their investigated constructions. In the finite framework this is partly a limitation built into the question: quantum structure and the Born rule were already adopted. The substantive tests concern what the supplied dynamics and state do or do not force. The campaign result is not a proof that no stronger law could exist. [R04, R06]

## 3. How one whole can contain individual ones

An individual object need not be a separately existing universe. It can be a persistent, dynamically distinguishable part of the same physical system.

For a justified finite tensor decomposition,

$$
\mathcal H=\mathcal H_A\otimes\mathcal H_B,\qquad
\rho_A=\operatorname{Tr}_B\rho_{AB}.
$$

The whole supplies the joint state; the local state specifies expectations for local measurements. Neither duplicates the other. In a gauge theory or gravity, constraints generally obstruct a naive spatial tensor product. The safer general statement is restriction of a state \(\omega\) to a subalgebra \(\mathcal A_A\):

$$
\omega_A=\omega|_{\mathcal A_A}.
$$

The observables and boundary/edge treatment must be justified. This document does not derive a preferred decomposition of the actual universe. [S02]

For the finite tensor case, with \(\hbar=1\),

$$
\dot\rho_A=-i[H_A,\rho_A]
-i\,\operatorname{Tr}_B[H_{\mathrm{int}},\rho_{AB}].
$$

The second term depends on the joint state. Therefore “the part has a local state” does not imply “the local state alone determines its future.”

**CONDITIONAL demonstration.** For \(H=\sigma_x\otimes\sigma_x\), choose
\(\lvert\psi_\pm\rangle=(\lvert00\rangle\pm i\lvert11\rangle)/\sqrt2\).
Both parts initially have the same marginal \(I/2\); both preparations have the same mean energy. At \(t=\pi/8\), their first-part states differ by trace distance \(1/\sqrt2\). The information distinguishing their futures was in the correlations. The global evolution remains unitary. Our executable finite checks reproduce this result. [W01]

This answers the user's “whole 1 and individual 1” question constructively: **yes, under a stated physical decomposition**. A part can be individually measurable and persistent while belonging to a correlated whole. It does not follow that the whole has the same properties as each part, that every correlation is consciousness, or that individuals are illusions.

## 4. What emergence means here

“Emergent” means that a useful higher-level description follows from lower-level dynamics within a stated regime and accuracy. It does not mean that an unexplained property has been renamed.

Let \(\mathcal C\) map microscopic states to a chosen description, \(U_t\) evolve the microscopic state, and \(\Phi_t\) be a proposed effective evolution. A finite-state adequacy criterion is

$$
\sup_{\rho\in\mathcal D,\;0\le t\le\tau}
D\!\left(\mathcal C[U_t\rho U_t^\dagger],
\Phi_t[\mathcal C(\rho)]\right)\le\varepsilon,
$$

where \(\mathcal D\), \(\tau\), the distance \(D\) and tolerance \(\varepsilon\) are specified in advance. This defines what must be checked; it does not prove that such an effective map exists. Initial correlations may require memory or a larger retained state.

A physical object earns an identity by maintaining identifiable properties across interactions for a relevant duration. Depending on the system, these include conserved charges, a binding energy, a metastability barrier, stable collective modes or active maintenance. The definition and lifetime depend on scale. A wave packet, an atom, a star and an organism are not identical kinds of persistence.

An **internal environment** is the rest of the degrees of freedom that affect the part. A reservoir is an approximation to some of them. Friction, noise and thermalization require hypotheses about state, interactions and timescales; they are not automatic consequences of naming the rest “the bath.”

## 5. The emergence map of reality

The map is a dependency graph with overlapping processes, not a claim that every law was invented at a distinct instant. Effective particles, symmetries and phases change with conditions even when the underlying action is held fixed. Decohere-and-record processes can occur across many epochs.

| Regime or event | Approximate placement | What becomes possible | GRUT's contribution and boundary |
|---|---|---|---|
| Origin of quantum structure, spacetime and laws | No established chronology supplied here | The structures used by later descriptions | OPEN; these are not derived by the current framework |
| Primordial state and any inflationary phase | Before the hot thermal history; model dependent | Expansion history, perturbations, possible asymmetry sources | Model histories computed conditionally; state and action supplied |
| Reheating and charge transfer | Model-dependent early epochs | Thermal reservoirs and populated matter sectors | M2 calculations; connected acceptance remains unfulfilled |
| Electroweak crossover | \(T\sim160\) GeV in the Standard Model | Low-temperature Higgs regime and familiar particle masses | Standard dynamics; smooth crossover, not a proved GRUT event [S06] |
| QCD confinement and hadrons | Microsecond-scale early universe | Hadronic matter | Standard strong-interaction physics [S07] |
| Light-element nucleosynthesis | First minutes | Much of primordial helium and other light nuclei | Standard hot-history prediction [S08] |
| Recombination and photon decoupling | Roughly 380,000 years | Long-lived neutral atoms; freely propagating CMB | Standard cosmology [S08] |
| Gravitational structure and first stars | First hundreds of millions of years, with uncertain details | Stellar energy sources, later element enrichment | Standard structure formation; dark-matter identity remains open [S09] |
| Chemical environments and planets | Local, overlapping histories | Diverse chemistry and sustained energy flows | Established physical ingredients; no universal date or GRUT abiogenesis calculation |
| Persistent records and information use | Many epochs and scales | Historical correlations and feedback | Standard mechanisms plus finite demonstrations; origins must be specified |
| Living systems and cognitive agents | Local evolutionary histories | Self-maintenance, learning, flexible action | No derivation from the GRUT action to biology |
| Subjective experience | No agreed physical onset criterion | The fact to be explained by a consciousness theory | OPEN; not inferred from persistence or mutual information |

Processes can have a common origin or interact without being identical. Simultaneous onset is a lead for a joint mechanism, not evidence for one. The proposed common source of dark and baryonic asymmetry was a concrete test of this idea; four specified source maps failed their registered condition. [R11]

No arrow in the map claims that life is inevitable, that complexity must increase everywhere, or that evolution is directed toward us. Nor is cosmic acceleration dated as the “creation of vacuum energy”; a component can exist long before it dominates expansion.

## 6. Time, thermodynamics and records

There are two levels of time in use. The cosmological calculations use ordinary coordinate/cosmic time in a specified gravitational model. A relational construction instead defines changes by correlations with an internal clock. The Page–Wootters proposal is an example of the latter, not a completed GRUT quantum-gravity theory. [S03]

In a finite clock example, if clock and measurement effects commute,

$$
p(a\mid T=t)=
\frac{\operatorname{Tr}[(\Pi_t\otimes E_a)\rho_{\mathrm{phys}}]}
{\operatorname{Tr}[(\Pi_t\otimes I)\rho_{\mathrm{phys}}]}.
$$

The gravitational case needs physical, gauge-compatible observables and a valid state measure. The Hamiltonian constraint \(\widehat{\mathcal H}_\perp\Psi=0\) is not a measurement that matter plus gravitational energy sums to zero; a compact universe generally lacks the asymptotic total-energy observable needed for that statement.

**Records demonstration.** For eight clock readings, let

$$
|\Psi\rangle=\frac1{\sqrt8}\sum_{t=0}^{7}|t\rangle W_t|\psi_0\rangle .
$$

Seven reversible gates copy a system's computational-basis value into a blank record, copy it to a second register, and later undo both copies. A positive propagation constraint annihilates the history. With initial system amplitudes \(\sqrt{0.3}\) and \(\sqrt{0.7}\), the record entropy when written is \(0.6108643\) nats. It returns to zero after erasure. The maximally mixed state on this model's eight-dimensional physical subspace gives the identity conditional state at every reading. [W01]

Thus a stationary history can contain nontrivial internal correlations. This example also shows that its constraint plus maximum ignorance does not select a history of accumulating records. The blank registers and amplitudes are supplied by the preparation. This finite result is not a universal no-go theorem for all state-selection proposals.

For a finite product preparation \(\rho_A\otimes\tau_B\), with
\(\tau_B=e^{-\beta H_B}/Z_B\), any joint unitary gives

$$
\Delta S_A+\beta\Delta E_B
=I(A:B)_{\mathrm{final}}+D(\rho_{B,\mathrm{final}}\Vert\tau_B)\ge0.
$$

The sign convention is energy gained by \(B\). This is the standard equality behind a finite-reservoir form of Landauer's principle. It requires the specified Gibbs state and initial lack of correlations; “passive environment” alone is insufficient. An interacting total Hamiltonian can conserve total energy while the interaction-energy contribution also changes. Our checks include that total Hamiltonian. [S05, W01]

Global unitary evolution preserves von Neumann entropy; local entropies can rise or fall. An explicit product-state example rises from zero to \(\ln2\) and returns to zero. Persistent thermodynamic arrows need a suitable state, coarse-graining and dynamical regime. Their observed direction is not obtained for free from closure.

Decoherence explains suppression of locally accessible interference through correlations with other degrees of freedom. It does not alone provide a unique physical outcome rule or solve every formulation of the measurement problem. [S04]

## 7. Matter, binding and the M1/M2 branch

The hidden matter sector is a candidate addition to known physics, not the meaning of “reality” itself. Its practical purpose is to test whether a declared dynamics can create charge, preserve it, form lasting structures and produce an acceptable relic.

A light-sector version uses a complex charged field \(X\), a real binding field \(\phi\), and

$$
V=g_X^2\phi^2|X|^2+\lambda_X|X|^4
+\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
-\epsilon v(\phi-v)^3 .
$$

The cited M1 light benchmark has \(g_X=0.1\), \(\lambda_X=0.001\),
\(\lambda_\phi=0.03\), \(\epsilon=0.002\), and \(m_X=g_Xv\approx92.4\) MeV.
These are model inputs, not values inferred from closure. The M2 card adds further sectors and histories; values from different cards must not be spliced silently. [R07]

At \(X=0\), \(S=\phi/v\), the stationary derivative factors exactly:

$$
\frac{\partial V}{\partial\phi}
=\lambda_\phi v^3(S-1)(S^2+0.8S+0.2).
$$

The quadratic discriminant is \(-0.16\), so this tree potential has one real stationary point, \(S=1\), and a bounded quartic asymptote. This removes the degenerate sign-vacuum problem of the old light potential at zero density. It does not prove all finite-density, finite-temperature or full-action stability claims.

The gauged completion employs charges \((1,3,6)\) for \(X,Y,\Sigma\), with interactions such as \(\kappa Y^*X^3+\mu\Sigma^*Y^2+\mathrm{h.c.}\). A \(\Sigma\) condensate can leave a residual \(\mathbb Z_6\) gauge symmetry. Approximate X-number and exact total gauge charge are different quantities. Pairwise scalar binding does not by itself give a charge-selective absorbing core.

The connected prediction must be an energy and population ledger:

$$
\rho_{\mathrm{hidden}}(t)=
\int dQ\,d\varepsilon\;n(Q,\varepsilon,t)E(Q,\varepsilon)
+\sum_i\int\frac{d^3p}{(2\pi)^3}E_i(p,t)f_i(p,t)
+\rho_{\mathrm{unresolved}}(t).
$$

The decomposition must avoid counting a coherent mode both as part of an object and as free particles. A particle interpretation needs a suitable frequency/interaction regime. The unresolved term must shrink under controlled matching, not become an adjustable missing-energy account.

If only positive-charge objects and unit-charge free particles are included, charge conservation gives

$$
\frac{N_+}{Q_{\mathrm{net}}}=1-\eta+p,\quad
\eta=\frac{Q_{\mathrm{objects}}}{Q_{\mathrm{net}}},\quad
p=\frac{N_-}{Q_{\mathrm{net}}}.
$$

It is possible for \(\eta>1\) when negative free charge compensates the objects. Counts are not mass fractions; the final abundance requires \(E(Q)\), the distribution of charges, excitation energy and all surviving reservoirs.

The conversion \(\bar X+\mathrm{core}(Q)\rightarrow X_{\mathrm{fast}}+\mathrm{core}(Q-2)\) conserves charge. In a large-core approximation the outgoing energy increases by \(2\,dE/dQ\). It transfers core energy into free particles; the reverse reaction remains possible above its threshold. A neutral scalar that subsequently decays to \(X\bar X\) is not net antiparticle removal. These facts prevent a capture probability from being mislabeled a cosmological cleanup efficiency.

### The latest result changes what “working” can mean

The completed, registered scan contains 135 preparation keys at each of two simulated inflation scales. **Zero of the 135 keys passes the strict M2-X rule.** There are 37 self-consistent free-X roots inside the bracket before imposing the threshold band. Only two have numerical interpolated thresholds; the smaller is \(2.18356\times10^{-5}\), 5.30 times the registered upper limit \(4.12\times10^{-6}\). [R08, R09]

Two isolated Q-branch roots meet their separate registered root condition. They do not establish actual formation of a Q-matter population. Fifteen in-band free-X extrapolations lie below the simulated bracket and were correctly excluded from the strict verdict.

A separately computed lower-scale scalar history changes the radiation fraction near hand-over from about \(0.83214\) to \(0.66959\), retaining the transferred higher-scale fermion drain. That is a 16.26 percentage-point change in this controlled comparison. A new kinetic splice and a full three-dimensional Dirac calculation were not performed. [R10, R13]

These are failures and dependencies of the stated numerical closure and sampled model. They do not exclude every preparation, coupling or field content. They also do not permit declaring the current parameter point successful.

## 8. State selection and coupling relations

The working state class contains an expanding background, specified perturbations, a low-entropy orientation and a model-specific phase/holder preparation. GRUT has not generated a unique physical state from weaker principles.

A state functional \(\mathcal B\) must respect gauge and gravitational constraints if used in a continuum constrained theory. Writing an initial boundary functional may encode a choice, but does not explain it. Penalizing initial curvature, entanglement or entropy introduces a rule that must be independently justified and checked.

A useful symmetry statement has precise hypotheses. If a **unitary** \(U\) satisfies

$$
[U,H]=0,\quad U\rho_0U^\dagger=\rho_0,\quad U^\dagger Q U=-Q,
$$

then \(\langle Q(t)\rangle=-\langle Q(t)\rangle=0\). A symmetric ensemble can still contain branches with opposite nonzero charges. Achieving a coherent sign across our observable region requires an appropriate state or mechanism.

An **antiunitary** time-reversal transformation generally relates forward evolution to reversed-time evolution. The same-positive-time conclusion above cannot simply be copied to that case. The finite model's conjugation must not automatically be called physical C or CP. This qualification prevents an overbroad preparation theorem.

An asymmetric ingredient can be an internal field, state or spontaneous process. “Closure does not choose the phase” does **not** mean the universe must receive an influence from outside.

The co-genesis test asked whether one specified CP-odd source could set both a large holder phase and the much smaller effective leptogenesis parameter. Four registered maps failed on their stated histories. A loop-induced relation remains an untested possibility requiring its own action, rate calculation and declared test. [R11]

To count as a new coupling relation, a candidate must remove an independently variable parameter:

$$
\mathcal R(\vartheta)=0,
$$

with its derivation, matching prescription, quantum corrections and an admissible comparator that violates it. A fitted equality or a relation obtained by assuming the answer does not qualify. No GRUT-specific universal relation is supplied here.

## 9. Records, agency and consciousness

These are three different explanatory tasks.

**A record** is a physical correlation that can support a specified inference about another event or variable. Its reliability and persistence can be measured. A record need not be alive or conscious.

**Information-using agency**, operationally, requires that stored information helps a system choose actions relevant to a stated outcome. Scrambling the informative correlation while preserving other relevant features is one way to test that causal role. This approach has published precedents; it is not a new GRUT definition of life. [S10]

**CONDITIONAL demonstration.** Let an environment bit be equally likely zero or one. A supplied sensor records it with error 0.1; a supplied controller acts according to its memory. Success is action matching the environment. Success is 0.9 with the correlation and 0.5 after scrambling the correlation while retaining both marginals. A fixed-action control also scores 0.5. Mutual information is about 0.531 bits. [W01]

The intervention establishes a role for information in this prepared controller. It does not derive the controller's origin, its energetic cost, metabolism, evolution or subjective experience. A simulation of adaptive behaviour is not by itself a consciousness detector.

**Conscious experience** remains a separate empirical and conceptual problem. Existing theories make competing predictions; a major preregistered comparison of two theories challenged aspects of both rather than producing an agreed universal criterion. [S11] Observer-relative observables in gravitational models concern physical access and dressing; they do not imply that consciousness creates spacetime. [S15]

The user's suggestion of underlying consciousness can be retained as an explicitly open branch:

| Possibility | What must be added | What would distinguish it |
|---|---|---|
| Experience arises with suitable physical organization | A bridge from specified dynamics to a testable experience criterion | Predictions across systems and interventions that outperform alternatives |
| Experience is fundamental but leaves all physical predictions unchanged | A philosophical interpretation and consistency account | No physical discriminator presently specified |
| An additional selector changes physical outcomes | An explicit dynamical or probabilistic law, with conservation and consistency checks | Held-out data differing from ordinary physics under the same declared inputs |

Apparent purpose can also arise from selection, feedback, planning and observer selection effects. A surprising pattern, fine-tuning or incomplete calculation is not evidence of a cosmic intention. An anomaly would first establish a problem with a model; identifying it as intention requires further discriminating evidence. No such evidence is claimed in this packet.

## 10. Vacuum density and gravity

In a nongravitational fixed-background model,

$$
e^{-i(H+cI)t}\rho e^{i(H+cI)t}
=e^{-iHt}\rho e^{iHt}.
$$

This exact invariance is included in our finite checks. It explains why these matter-evolution properties cannot determine an absolute energy offset. It is not a no-go theorem for gravity.

For a dynamical metric the added term
\(-C\int d^4x\sqrt{-g}\) shifts stress-energy by
\(\Delta T_{\mu\nu}=-C g_{\mu\nu}\). The total renormalized vacuum density therefore matters gravitationally. Relabeling a matter constant and compensating counterterm does not change the physical total. [R05, S12]

The M2 card supplies approximately
\(\rho_{\mathrm{vac}}=2.52\times10^{-47}\ {\mathrm{GeV}}^4\) as an observational calibration. This document does not present it as a new 2026 measurement or a GRUT prediction. The outstanding task is its sign, size and stability under quantum corrections and heavy thresholds. [R07]

Two previous ideas remain separate:

- A self-sustained vacuum approaching equilibrium is a published Klinkhamer–Volovik proposal. The imperfect-vacuum/disequilibrium interpretation is theirs, not a GRUT invention. Its embedding in a viable, independently predictive GRUT cosmology remains open. [S13, S14]
- In four-dimensional de Sitter Einstein gravity, horizon entropy and Friedmann's equation give \(S=8\pi^2M_P^2/H^2\) and \(\rho_\Lambda=3M_P^2H^2\), hence \(\rho_\Lambda=24\pi^2M_P^4/S\). This is a relation within that assumed regime. It predicts the vacuum scale only if \(S\) is fixed independently. [S16]

A finite horizon entropy is not, without further argument, a proof of a finite-dimensional global Hilbert space. Nor can one simultaneously assume a fixed de Sitter future and a vacuum relaxing to zero as though they were the same hypothesis.

The “final equilibrium” postulate is an additional choice. Inferring absolute electroweak vacuum stability also requires specifying the full potential, ultraviolet validity and meaning of final equilibrium. It is not forced by closure. A metastability estimate with uncertain mass matching does not establish that GRUT's particular hidden scalars exist. [R06]

## 11. The five outputs the user asked for

| Requested output | Defensible answer now | What is needed to upgrade it |
|---|---|---|
| **State** | A supplied, explicit preparation class; not a unique cosmic state | A constraint-compatible selection rule and distinguishable consequences |
| **Coupling relation** | Ordinary consistency and matching relations; model choices remain | One independently derived relation that removes a genuine freedom and survives corrections |
| **Observable prediction** | Conditional model outputs and registered failed tests; no new universal GRUT signal | A frozen model that passes its own connected checks and a discriminator fixed before tuning |
| **Final object abundance** | Not certified for the actual history; latest strict M2-X scan fails | Controlled formation, particle matching, heating, forward/reverse reactions and converged cosmological evolution |
| **Vacuum density** | Measured calibration in the action | An independently fixed renormalized value, radiative stability and a distinct consequence |

The multi-target recoil-shape idea is a possible discriminator if a viable portal model survives its production obligations. Gravity-only agreement with cold dark matter is compatibility, not evidence uniquely favouring GRUT. A predicted quantity calibrated to the same data is a fit; a new test must use information not spent in that calibration.

## 12. A bounded completion program

The account above can stand now. Further physics should be pursued through concrete dependencies:

1. **Finish or close the present matter instance.** Compute a kinetic splice on an actually evolved lower-scale history and assess the fermion closure with a suitable controlled method. Treat extrapolated roots as diagnostics. Full 3D fermion resolution is a distinct problem; changing only the scalar lattice is not a substitute.
2. **If object matter is pursued, connect its populations.** Match persistent cores and justified particle modes without double counting; carry neutral reservoirs, hard sources and reverse conversion; obtain the final mass and velocity distribution with convergence and energy/charge checks.
3. **Test one genuine reduction of inputs.** A specified state-selection mechanism or radiative co-genesis relation must be declared before the decisive calculation, with an admissible failure outcome. An added interaction is a new model version.
4. **Fix a discriminator before fitting toward it.** Compare the same prepared models against a baseline using uncertainty-aware, held-out observations.
5. **Keep foundations as separate unresolved work.** Quantum/gravitational origin, the low-entropy state, vacuum density and experience need their own constructions. Success in a dark-matter run cannot discharge them.

The framework's conceptual identity is therefore precise: **one interacting world, many conditionally autonomous and persistent parts, history carried in correlations, and higher-level behaviour supported by specific physical organization.** Its mathematical examples establish compatibility and conditional mechanisms. Its cosmological candidates test additional matter hypotheses. Its unfinished foundations remain visible rather than being filled by assertions of closure or intention.

The current synthesis is complete as an account of those commitments, results and limits. The universe's fundamental explanation is still a research task.
