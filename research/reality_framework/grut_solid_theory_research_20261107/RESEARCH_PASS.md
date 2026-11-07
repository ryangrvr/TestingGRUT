# What would make GRUT a solid physical theory?

**Research pass, 7 November 2026 — Codex, for Ryan and Perplexity.**

**Finding:** the present construction can support a rigorous operational framework and a controlled effective physical model. The missing bridge is a quantitative identification of the *same physical experiments* across its finite realizations and continuum target. A unique microscopic Hamiltonian is not required for an effective theory. A distinctive claim about nature does require a restriction or prediction beyond the established physics already supplied.

This corrects an emphasis in our previous discussions: failure to derive a unique microscopic access structure is a failure of that proposed derivation, not a general disqualification of every theory that adopts interactions and initial conditions.

The recommendation is to complete the existing interacting detector-history bridge with matched physical parameters and an independent continuum description. This connects the current construction; it does not introduce another dark-matter model or require a larger cosmological campaign. No new GRUT law was established by this research pass. **Zenodo remains paused.**

## 1. Evidence and scope

The repository examined is `ryangrvr/TestingGRUT`, PR #14 head `7ecbf69f584dd3e8be047a85768e85452f2c2c3a`. The integrated theoretical source is [GRUT_Integrated_Theory.md](../grut_integrated_candidate_20261106/GRUT_Integrated_Theory.md), including its free-photon matching and clarification. The new interacting source is [the QED packet](../grut_qed_detector_histories_20261107/). The [previous Codex review](../grut_qed_review_20261107/REVIEW.md) records exactly which numerical checks were repeated. This pass adds primary-literature research, a conditional error-bound derivation, and reproducible parameter arithmetic. It does **not** add production simulations or independently certify the full DMRG grid.

| Present result | What it establishes | What it does not establish |
|---|---|---|
| Finite algebra, state, Hamiltonian and instruments | A specified probability model, conditional response and sufficient-observable constructions | Selection of the physical laws or cosmic preparation |
| Connected clock–record–action realization | These mechanisms coexist in an explicit quantum model with supplied resources | That closure produces clocks, low entropy or observers in every model |
| Free photon bridge | The stated smooth vacuum histories converge in the documented ordered limits | A uniform all-state theorem, interacting matter or dynamical geometry |
| QED spectra and finite histories | Checks of standard Schwinger dynamics and of Gaussian sufficiency within the tested constructions | Controlled interacting continuum history matching, or a new physical signal |

The reviewed numerical agreement remains useful. It is necessary to separate solver verification, approximation control and empirical evidence: increasing the first does not automatically provide the other two.

## 2. The standard for a solid theory

A working physical theory needs a declared domain, dynamical and state assumptions, an observable dictionary, controlled approximations, and predictions that can fail. It can use measured parameters and an assumed preparation. It need not answer why every constant has its value before it becomes useful.

Effective field theory organizes predictions by resolution and an expansion order. Different microscopic descriptions can yield the same low-energy observables after matching their coefficients; regulator dependence must cancel to the retained accuracy. This is established methodology [S1]. **My application to GRUT:** the photon Hamiltonian nonselection result may be compatible with universality. If two families give the same declared continuum experiments, that is not a defect in the effective description. It is still a counterexample to any claim that those experiments uniquely reconstruct the microscopic family.

We should therefore separate three achievements:

1. **Sound framework:** definitions and conditional theorems are correct, with novelty assessed against prior work.
2. **Controlled physical model:** a specified family predicts physical observables in a justified regime, with matching and error control.
3. **Distinctive GRUT physics:** an additional restriction survives comparison against an otherwise admissible model using the same calibration, and ultimately an independent physical test.

The current documents support the first in their declared finite domain, subject to independent mathematical review. The second is partial. The third has not been demonstrated. Calling all three “complete” would conceal the remaining work.

## 3. One integrated specification, with physical access made explicit

Keep the parent object already adopted,

$$
\mathfrak G=(\mathcal A,\omega,H,\mathfrak I),
$$

but make its physical interpretation conditional on one common specification:

- the fields or microscopic degrees of freedom, constraints and Hamiltonian;
- the independent renormalized parameters at a stated scale and effective order;
- the allowed preparation class, including clock, memory and low-entropy resources;
- the interactions that implement the instruments, the readout dictionary, and their energy, spatial and temporal resolution;
- a bounded class of protocols and an approximation budget.

For a protocol $\pi$, its outcomes obey the adopted quantum instrument rule

$$
p(\mathbf r\mid\pi,\omega)
=\operatorname{Tr}\!\left[
\mathcal J^{\pi}_{r_n}\mathcal U_{n:n-1}\cdots
\mathcal J^{\pi}_{r_1}\mathcal U_{1:0}(\rho_\omega)
\right]
\tag{R1}
$$

in a representation where the state is a density operator. Each sum over outcomes is trace preserving. A protocol can depend on earlier stored outcomes. Algebraic states and dual maps supply the corresponding formulation when no global density matrix description is appropriate.

**Access is physically admissible when its instrument is induced by the declared coupled dynamics and preparation, or is a controlled approximation to such an instrument.** For a finite system and probe this means

$$
\mathcal J_r(\rho)=\operatorname{Tr}_{P}
\left[(I\otimes\sqrt{E_r})U_{SP}
(\rho\otimes\tau_P)U_{SP}^{\dagger}
(I\otimes\sqrt{E_r})\right],
\qquad \sum_r E_r=I.
\tag{R2}
$$

The product preparation in (R2) is an assumption; initially correlated probes require the joint-state formulation. This condition links the access primitive to the adopted dynamics. It does not reconstruct access, readout algebras or preparation from a stationary state.

For the existing Weyl instrument, a probe prepared in $|+\rangle$, a controlled $e^{i\varphi}W$, and an $X$-basis probe readout give exactly $K_\pm=(I\pm e^{i\varphi}W)/2$. This proves a finite unitary realization. A spacetime-local implementation of an extended field instrument requires a further argument: account for its probe coupling, duration and how spatially separated records are collected. An instantaneous global controlled operation is not automatically a local relativistic measurement.

Fewster and Verch provide a continuum construction using a local system–probe coupling, induced instruments and causal composition [S4]. It supplies a relevant theoretical standard for GRUT's detector bridge, not a request to build laboratory hardware. Whether the *particular* GRUT instruments meet that standard remains to be shown.

An observer's perspective can then mean conditioning on its physical records, for example $p(r_B\mid r_A,\pi)$, with finite clock accuracy included. Incompatible alternative interventions are different protocols; no joint assignment of all unperformed measurements is assumed. Observer dependence is therefore testable through specified correlations and disturbances. It does not by itself select a Hamiltonian or introduce consciousness into the dynamics.

## 4. The matching requirement can be made calculable

For cutoff $a$, let $\mathfrak G_a$ be a declared family with physical parameters held fixed by matching, rather than bare parameters kept numerically fixed without justification. The target remains equation (14) of the integrated document:

$$
d_{\mathscr P,\mathcal C}(a)=
\sup_{\rho\in\mathcal C,\pi\in\mathscr P}
\operatorname{TV}(p_{a,\rho,\pi},p_{\mathrm{target},\rho,\pi})
\leq\epsilon.
\tag{R3}
$$

Specify the energy class $\mathcal C$, maximum number of operations, finite time interval, smearings, couplings, boundary conditions and clock support in $\mathscr P$. A finite collection of evaluated protocols tests only that collection. It does not establish the supremum over an untested continuous class.

**Conditional error proposition.** Work first in a common comparison representation, including all outcome registers and any correlated reference system. Let $\Phi_j$ and $\widetilde\Phi_j$ be exact and approximate trace-preserving channels for stage $j$. An instrument is included as a channel retaining its classical outcome, rather than as a normalized postselected branch. Assume

$$
\tfrac12\|\rho_0-\widetilde\rho_0\|_1\leq\delta_0,
\qquad
\tfrac12\|(\Phi_j-\widetilde\Phi_j)(\tau)\|_1\leq\delta_j
\tag{R4}
$$

for all states needed along the comparison, including those reached by the approximate process. Trace-distance contraction and the triangle inequality give

$$
\begin{aligned}
D_j
&=\tfrac12\|\Phi_j(\rho_{j-1})-
             \widetilde\Phi_j(\widetilde\rho_{j-1})\|_1\\
&\leq D_{j-1}+\delta_j,
\qquad
\operatorname{TV}(p,\widetilde p)
\leq D_n\leq\delta_0+\sum_{j=1}^{n}\delta_j.
\end{aligned}
\tag{R5}
$$

Adaptive control is included by letting the channels read earlier outcome registers. Supremizing these hypotheses over the declared state and protocol classes gives a bound of the form (R3). Clock error, probe error, dynamics error and truncation error can enter as separate stages or separately bounded contributions. Boundary and state mismatch must be included as well. This bounds joint histories; conditioning on a rare record can amplify the error, so conditional claims also require a lower bound on that record's probability.

This is ordinary contractivity and telescoping, not a new GRUT law. Its useful consequence is a concrete obligation for every arrow in the construction. It also exposes the limitation: **constructing the comparison representation and establishing the bounds in (R4) are substantive work.** Finite-density-matrix reasoning does not establish a uniform theorem for a continuum field algebra, arbitrary energy, arbitrary operation count or simultaneous cutoffs. Ordered limits remain ordered limits unless the required uniform estimates are supplied.

Process-tensor methods already characterize multi-time responses to interventions and operational memory [S3]. GRUT's reduced-state claim must be compared with that framework. Matching passive correlators alone is insufficient when a probe can alter subsequent evolution. Exact predictive sufficiency also need not imply useful computational compression.

## 5. A concrete issue in the present QED bridge

The existing code compares exact finite-lattice histories to a Gaussian prediction constructed from **that same lattice's own mean and two-point function**. It therefore diagnoses Gaussian sufficiency. It is not an independent match to continuum detector histories.

### Mass matching

For staggered Schwinger fermions, Dempsey et al. derive the bulk mass-improvement prescription

$$
m_{\mathrm{lat}}=m_R-\frac{g^2a}{8},\qquad
\frac{m_{\mathrm{lat}}}{g}=\frac{m_R}{g}-\frac1{8\sqrt{x}},
\quad x=(ag)^{-2}.
\tag{R6}
$$

Their open-boundary study also finds improved convergence, while explaining that boundaries break the lattice translation symmetry used in the derivation [S5]. Independent MPS work studies additive mass matching versus spacing, volume and other parameters [S6]. Thus (R6) is a concrete improvement to investigate, not a certified exact mass calibration for our finite interval.

With the packet's $\mu=2\sqrt{x}\,m_{\mathrm{lat}}/g$, (R6) gives the simple shift $\Delta\mu=-1/4$. At the frozen history points:

| $N,x$ | Reported $Na\,g$ | Endpoint span $(N-1)a\,g$ | Bulk-improved $m_{\mathrm{lat}}/g$ for $m_R=0$ |
|---|---:|---:|---:|
| 10, 4 | 5 | 4.5 | −0.0625000 |
| 14, 9 | 4.6666667 | 4.3333333 | −0.0416667 |
| 20, 16 | 5 | 4.75 | −0.0312500 |

The unshifted bare-zero sequence can still approach the massless continuum theory because the shift vanishes as $a\to0$. The published improvement does **not** retrospectively invalidate the two passing continuum spectrum fits. It does mean that finite-spacing deviations should not all be attributed to an intrinsic failure of the Gaussian continuum description. A revised matching calculation is a new calculation and must preserve the original registered results.

### Boundary and detector matching

Okuda derives the interval dictionary and shows that the conventional staggered Gauss law induces boundary charges; the map depends on site parity and boundary choices. His site-endpoint convention is $L=(N-1)a$. Changing to a modified Gauss law can shift both mass and theta, so improving mass is not a license to alter the boundary problem silently [S7].

Our code instead reports $L=Na$ and averages over link centers. The two length conventions differ by $a$, a finite-cutoff effect to match explicitly. Neither column in the table is constant along the entire old ladder. Select one interval dictionary, fix its physical length and background flux, and match the same normalized smearing and instrument phase. The lattice's measured one-point function must not be silently substituted for the continuum prediction during validation.

For matched massless vacuum fluctuations, the bosonization dictionary [S7] and normalized modes $u_n(z)=\sqrt{2/L}\sin(n\pi z/L)$, $\Omega_n^2=(n\pi/L)^2+g^2/\pi$, give an independent covariance:

$$
G_f(t)=\left\langle\delta F_f(t)\delta F_f(0)\right\rangle
=\frac1\pi\sum_{n\geq1}
\frac{|f_n|^2}{2\Omega_n}e^{-i\Omega_n t},
\qquad f_n=\int_0^L dz\,f(z)u_n(z),
\tag{R7}
$$

Here $F_f=\int dz\,f(z)E(z)/g$, with normalized smooth $f$. Match the mean separately and bound the mode tail. This Gaussian massless reference cannot replace the massive interacting target.

Houtz and Ingoldby also provide an independent Hamiltonian-truncation approach on an interval [S8]. The inspected author repository explicitly supports the **massless** Schwinger model; it is not a ready-made massive-history solver. Its README flags a corrected QCD-specific bug and states that the Schwinger implementation is unchanged. Pin and verify any reused version. No external solver from that repository was executed in this pass.

## 6. The prioritized calculation

**Physical claim:** the same interacting field preparation and probe protocol have regulator-independent joint record probabilities within a stated error budget.

This is the missing connection between the general finite framework and an interacting continuum realization. It is a validation of that connection, not a GRUT-specific experimental signal.

Complete the following as one connected comparison, with the dictionary fixed before production:

1. Use the existing Schwinger target with specified physical $g,m_R,\theta$, interval boundary state and neutral sector. Freeze the renormalization and endpoint conventions. Treat regulator spacing as a computational parameter for this calculation.
2. Match a fixed physical interval, smearings, durations, phases and probe couplings. Identify whether the ideal Weyl instruments have a controlled local-probe realization; record any idealization error separately. Do not obtain agreement by independently fitting each cutoff's detector strength to the validation probabilities.
3. Calibrate the implementation against the analytic massless interval reference, including the boundary-dependent mean. Then compare the massive interacting histories with an independently constructed interacting reference, such as a properly renormalized bosonized truncation. Its normal ordering, counterterms and finite-volume data must be matched; a second run of the same DMRG code is not that reference.
4. Evaluate all four two-outcome history probabilities, then a fixed three-time intervention protocol. Compare both the full interacting target and the Gaussian reduced comparator on the same preparations. The third operation tests predictive information beyond the already examined two-time examples.
5. Separate state preparation, finite volume, spatial discretization, field/energy truncation, real-time evolution and instrument errors. Use (R5) wherever its hypotheses can actually be established. Check convergence of the target observable, not just energies or probability normalization.

**Proposed acceptance scale, not an executed or sealed new test:** aim for total-variation error at most $10^{-3}$ for each frozen protocol, with each independent numerical implementation contributing at most $2.5\times10^{-4}$ and the remaining budget explicitly allocated. This scale is chosen to resolve the existing probability differences of order $10^{-2}$, not as a universal physical constant. The detailed parameter card and feasibility check must precede a production run. A smaller observed non-Gaussian effect requires correspondingly tighter control or an unresolved verdict; do not enlarge the tolerance after seeing the result.

Start by deriving the boundary/instrument dictionary and implementing the inexpensive analytic reference. Do not launch another 192-site grid before this is defined. If a reliable massive reference or instrument-error estimate is infeasible within a bounded run, retain the finite benchmark and identify the missing estimate. A discrepancy larger than the combined budget rejects the matching claim in that regime; convergence failure leaves it unresolved. Neither result, by itself, refutes quantum electrodynamics or establishes new GRUT physics.

## 7. Broader physical obligations, without moving the finish line

| Obligation | What a declared effective model may adopt | What a stronger GRUT claim would have to establish |
|---|---|---|
| Parts and operational access | A local observable algebra and physically specified interactions | A selection principle if claiming these arise uniquely from closure; operationally induced subsystems already have substantial prior art [S2] |
| State and thermodynamic arrow | A specified preparation and low-entropy resources | A physical origin or probability law if claiming to explain the cosmic boundary condition; unitary evolution alone supplies neither |
| Chirality | An anomaly-consistent continuum target or a specified suitable regulator | A demonstrated chiral matter construction if claiming a microscopic SM recovery; Ginsparg–Wilson avoids the naive doubling hypotheses but does not alone construct the entire SM [S10] |
| Gravity | GR as an effective theory in a justified low-energy regime [S11] | Constraint preservation, stress conservation and a controlled metric/observable match; emergence or uniqueness is an additional claim |
| Shared relativistic propagation | Lorentz-covariant matching in the removed-regulator theory | If lattice spacing is physical, specify its scale and preferred structure and control radiative corrections, not only tree dispersion |
| Contents of our universe | Empirically supplied parameters and declared matter extensions | A tested account of missing content if claiming comprehensive physical coverage |

The integrated target currently writes the dimension-four SM with no right-handed neutrino or explicit neutrino mass operator. That is a defined truncation, not an account of observed neutrino oscillations. A broader physical version must specify its neutrino-mass sector, for example a dimension-five lepton-number-violating operator with independently declared coefficients, or a different explicit extension. This need not be derived before studying a QED sector, but cannot be hidden by the label “Standard Model” [S12]. The dark-matter and cosmological-state obligations likewise remain separate from the detector construction.

For a **physical** lattice, preferred-frame effects can generate lower-dimension Lorentz-violating terms through interactions; a tiny tree-level $a^2p^2$ dispersion correction is not itself protection [S9]. This is a conditional naturalness problem, not a universal exclusion of emergent-relativistic theories. For a **removed computational regulator**, its finite-spacing anisotropy is an error to control, not a laboratory prediction. These are different model commitments. The current Cherenkov estimate should not be used to blur them.

## 8. What would count as genuinely additional GRUT physics?

The existing integrated document explicitly adopts standard quantum dynamics and supplied operational structure. The research found no established additional GRUT restriction. Prior work already covers observable-induced subsystem structure [S2], multi-time operational memory [S3], and local QFT measurement [S4]. This is a literature comparison, not a proof that every statement in our documents has previously appeared. Any proposed new theorem still needs a precise theorem-by-theorem comparison.

For a new physical restriction, write an invariant relation or dynamical term. For example, if the proposal restricts renormalized coefficients through $F(c_i(\mu),\mu)=0$, count independent physical inputs after removing gauge and field-redefinition redundancy. For a relation asserted at all scales, consistency requires

$$
\left(\partial_{\ln\mu}+\sum_i\beta_i\partial_{c_i}\right)F=0
\quad\hbox{on }F=0,
\tag{R8}
$$

at the claimed order, together with threshold matching. A relation imposed only at a stated matching scale instead has to be run to the measurement scale. A nonzero running correction is not automatically fatal, but must be computed and included. Adding independent couplings that restore every removed degree of freedom supplies no predictive reduction.

Supply a consistent model satisfying the restriction and an admissible comparator violating it, both fitting the same calibration. Freeze an observable that separates them by more than numerical, approximation and observational uncertainty. Data used to select the restriction are calibration, not untouched confirmation. A regulator-independent history theorem could be valuable new methodological work; a distinctive law of nature additionally requires a physical discriminator. We should not manufacture a law to obtain the latter label.

## 9. Decision from this pass

**Retain:** the integrated finite framework, connected realization, scoped photon theorem, Schwinger validation and negative results. Finish their literature attribution and the bounded corrections already requested in the QED handoff.

**Complete next:** one matched interacting detector-history construction with an independent target and an explicit operational error budget. Its result would decide whether the current central bridge is controlled, rather than merely internally reproducible.

**Keep separate:** microscopic uniqueness, a physical lattice hypothesis, full chiral-SM/metric recovery and a new GRUT empirical law. They are stronger claims with additional obligations. Requiring all of them before calling any effective sector sound is unnecessary; claiming all of them from that sector would be unjustified.

**Current verdict:** GRUT has a coherent framework and partial physical realizations. A solid finite-domain physical theory is a credible construction target. A distinctive, empirically supported GRUT theory of our universe remains unestablished. The research supplies a narrower, calculable route to the first and a precise standard for the second; it does not promise that a new law must exist.

## Sources and reproducibility

Primary sources checked on 7 November 2026. Bibliographic links identify the actual papers; claims about our code refer to the pinned repository above.

| ID | Source | Relevant contribution |
|---|---|---|
| S1 | Weinberg, [Effective Field Theory, Past and Future](https://arxiv.org/abs/0908.1964), 2009, especially the discussion of cutoff independence and effective order | EFT predictivity and matching |
| S2 | Zanardi, Lidar, Lloyd, [Quantum tensor product structures are observable-induced](https://arxiv.org/abs/quant-ph/0308043), PRL 92, 060402 (2004) | Prior operational subsystem framework |
| S3 | Pollock et al., [Non-Markovian quantum processes](https://arxiv.org/abs/1512.00589), PRA 97, 012127 (2018); [Operational Markov condition](https://arxiv.org/abs/1801.09811), PRL 120, 040405 (2018) | Multi-time interventions and memory |
| S4 | Fewster, Verch, [Quantum fields and local measurements](https://arxiv.org/abs/1810.06512), CMP 378, 851–889 (2020) | Local probe-induced instruments |
| S5 | Dempsey et al., [Discrete Chiral Symmetry and Mass Shift](https://arxiv.org/html/2206.05308v3), PR Research 4, 043133 (2022), eq. (3), section IV.5 | Bulk mass prescription and open-boundary qualification |
| S6 | Angelides et al., [Computing the Mass Shift of Wilson and Staggered Fermions](https://arxiv.org/abs/2303.11016), PRD 108, 014516 (2023) | Numerical lines of constant renormalized mass |
| S7 | Okuda, [Schwinger model on an interval](https://arxiv.org/html/2210.00297v2), PRD 107, 054506 (2023), sections II.2 and III | Boundary and bosonization dictionary |
| S8 | Houtz, Ingoldby, [Hamiltonian Truncation Framework for Gauge Theories on the Interval](https://arxiv.org/abs/2509.17890); [JHEP 10 (2026) 276](https://doi.org/10.1007/JHEP10(2026)276); [author code](https://github.com/james-ingoldby/HT-2d-GaugeTheory) | Independent massless reference route; code not run here |
| S9 | Collins et al., [Lorentz invariance and quantum gravity: an additional fine-tuning problem?](https://arxiv.org/abs/gr-qc/0403053), PRL 93, 191301 (2004) | Radiative preferred-frame effects |
| S10 | Lüscher, [Exact chiral symmetry on the lattice and the Ginsparg–Wilson relation](https://arxiv.org/abs/hep-lat/9802011), PLB 428, 342–345 (1998) | A qualified escape from naive doubling |
| S11 | Donoghue, [General relativity as an effective field theory](https://arxiv.org/abs/gr-qc/9405057), PRD 50, 3874 (1994) | Low-energy gravity as adopted EFT |
| S12 | Particle Data Group, [Neutrino Masses, Mixing, and Oscillations](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-neutrino-mixing.pdf), sections 14.1–14.2 | Physical limitation of the massless-neutrino SM truncation |

The S8 README inspected via the GitHub connector had blob SHA `f1413ee2e82a5d8740c48986e8abe3c937de9c3a`. Its repository's current contents must be pinned again before reuse. The research sources are not independent certification of GRUT's numerical claims.

Reproduce the parameter table, without modifying any sealed output, from the repository root:

```bash
python research/reality_framework/grut_solid_theory_research_20261107/matching_scales.py
```

The output is archived as [matching_scales.json](matching_scales.json), including the SHA-256 of the frozen history input. Equations (R2), (R5) and (R8) are explicit conditional deductions in this report. The proposed comparison in section 6 has not been executed or retrospectively declared preregistered.
