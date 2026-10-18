# GRUT Candidate 1

**18 October 2026 · explicit effective physical proposal · version 1**

GRUT Candidate 1 is Einstein gravity with ordinary quantum matter and a specified additional scalar sector that can support charge-stabilized organized matter. One real field changes the mass and binding environment of a populated complex field. A separate reporter field reads that organization through an explicit local interaction. All fields couple to the same dynamical metric.

The proposed operating domain is a weakly coupled, low-curvature effective field theory (EFT). A complete nonperturbative quantum-gravity state space and ultraviolet completion are **missing**. The worked calculation concerns a narrower semiclassical scalar configuration, a weak reporter, and linearized gravity. The proposal does not derive spacetime dimension, quantum probabilities, universal locality, or Einstein dynamics.

## Architecture decision

Exactly two architectures were compared before full specification.

| Requirement | Continue the moving-word architecture | Selected covariant field architecture |
|---|---|---|
| Populated interacting matter | Requires a new many-composite and positive-energy lifetime argument | A local interacting scalar action with an exactly conserved charge is specified |
| Shared propagation | Supplied relational generator rule; current medium is nonrelativistic and one-dimensional | Common metric and covariant principal kinetic terms are postulated |
| Gravity | No demonstrated compatible gravitational sector | Explicit Einstein action and universal stress coupling are postulated |
| Demonstrated regime | Moving two-source delivery/retention remains unresolved | Solve a charge-stabilized configuration and its reporter/gravitational observables |
| Main cost | Several unresolved composition problems before a full model exists | Much of known physics is adopted, rather than derived |

The second architecture is selected. This is a separately named successor. It does not modify CP-12, B1, B2, MC-E1, or their frozen evidence. The moving-channel bottleneck takes **option 3: replace the mechanism and calculate the replacement's obligations**. No historical transport, density, or geometry-formation theorem is inherited.

## 1. Postulates and their roles

| Postulate | Mathematical content | Physical purpose and obstruction addressed | Dependent conclusions | Motivation/status |
|---|---|---|---|---|
| P1: spacetime | An oriented, time-oriented, globally hyperbolic Lorentzian four-manifold equipped with a spin structure, with signature `(-,+,+,+)` and dynamical metric `g` | Supplies dimensional/causal structure that the single-end grammar and graph selection work did not derive | Covariant dynamics, clocks, shared characteristic cones | Standard framework adopted; dimension and orientation are supplied |
| P2: quantum framework | Complex amplitudes, unitary perturbative EFT evolution, Born probabilities and completely positive instruments | Supplies probabilities and a physical quantum lift without relabeling decoherence as outcome selection | All quantum measurement claims | Standard quantum postulate; no consciousness/outcome mechanism derived |
| P3: local field content | Standard Model with three right-handed neutrinos, plus real `phi,D` and complex `X`; local covariant kinetic terms | Gives actual interacting matter and detector variables; excludes arbitrary distant internal-state products at leading order | Causal field propagation and local feedback | Ordinary matter imported; new scalar content is a construction choice |
| P4: symmetries | Standard gauge symmetry; `U(1)_X × Z2_phi × Z2_D`; global `B−L` at EFT scope for Dirac neutrinos | Protects X charge and states which interactions are admitted | Charge-stabilized matter; particle content and counterterms | Explicit assumptions; no quantum-gravity protection of global symmetries claimed |
| P5: gravity | Einstein action with reduced Planck mass `M_P`; all stress is obtained by metric variation | Supplies the gravitational dynamics and universal leading coupling absent from the graph route | Metric response, free fall, clock shifts and lensing | Established gravitational law adopted, not an emergent GRUT derivation |
| P6: scalar interactions | The potential below, with positive couplings and `g_X^4 > lambda_phi lambda_X` in the binding class | Allows mutual matter/organization feedback and positive-energy binding | Scalar mass shifts, finite-charge organization, reporter scattering | Known soliton mechanism instantiated; values are supplied parameters |
| P7: matching and preparation | Renormalized matching prescription, an initial state, physical apparatus, and measured parameter vector are specified below | Prevents importing incompatible modules, an unmodeled reservoir, or a favorable postselection | Domain of every prediction | Necessary specification; preparation is supplied in the worked example |

The word “relational” refers to covariantly defined relations between physical fields and apparatus. It does not make supplied locality or the spacetime manifold emergent.

## 2. Complete leading action

Use natural units `hbar=c=1`. Write `g_X` for the new scalar coupling, reserving `g_{mu nu}` for the metric. The leading action is

\[
S_0=\int d^4x\sqrt{-g}\left[\frac{M_P^2}{2}R-\rho_\Lambda
+\mathcal L_{\rm SM+\nu_R}
-\frac12\nabla_\mu\phi\nabla^\mu\phi
-\nabla_\mu X^*\nabla^\mu X
-\frac12\nabla_\mu D\nabla^\mu D-U\right]+S_{\rm boundary},
\]

\[
U=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
+g_X^2\phi^2|X|^2+\lambda_X|X|^4
+\frac12(m_{D0}^2+h_D^2\phi^2)D^2+\frac{\lambda_D}{4}D^4.
\]

`S_boundary` is the Gibbons–Hawking–York term with the chosen boundary orientation, together with the appropriate asymptotic reference subtraction for an ADM energy. The asymptotically flat worked benchmark uses renormalized `rho_Lambda=0`; a nonzero cosmological density is a separately calibrated parameter of the same EFT family. The benchmark is not a prediction that the observed cosmological constant vanishes.

The displayed scalar potential is nonnegative for positive couplings. In the vacuum `phi=+v, X=D=0`, the new masses at leading order are

\[
m_X=g_Xv,\qquad m_\phi=\sqrt{2\lambda_\phi}\,v,\qquad
m_D=\sqrt{m_{D0}^2+h_D^2v^2}.
\]

The choice of the `+v` asymptotic vacuum is boundary/preparation data. Domain walls between the two Z2 vacua are allowed by the action but are outside the demonstrated regime.

### Ordinary matter is defined, not generated

The gauge group is `SU(3)_c × SU(2)_L × U(1)_Y`. For each of three generations use `Q_L:(3,2,1/6)`, `u_R:(3,1,2/3)`, `d_R:(3,1,−1/3)`, `L_L:(1,2,−1/2)`, `e_R:(1,1,−1)`, and `nu_R:(1,1,0)`, with Higgs `H:(1,2,1/2)`. The new scalars are gauge singlets. In tetrad notation the adopted matter Lagrangian is

\[
\mathcal L_{\rm SM+\nu_R}
=-\frac14G^a_{\mu\nu}G_a^{\mu\nu}
-\frac14W^I_{\mu\nu}W_I^{\mu\nu}
-\frac14B_{\mu\nu}B^{\mu\nu}
-(D_\mu H)^\dagger D^\mu H
-\lambda_H(H^\dagger H-v_H^2/2)^2
+\sum_\psi i\bar\psi e^\mu_a\gamma^a\mathcal D_\mu\psi
-[\bar Q_LY_dHd_R+\bar Q_LY_u\widetilde H u_R
+\bar L_LY_eHe_R+\bar L_LY_\nu\widetilde H\nu_R+\mathrm{h.c.}].
\]

Use the conventional particle-physics gamma matrices with `{gamma^a,gamma^b}=−2 eta^{ab}` for the stated `(-+++)` metric: `gamma^0` is Hermitian, the spatial matrices are anti-Hermitian, and `bar psi=psi† gamma^0`, with the corresponding tetrad spin connection. Gauge field strengths and covariant derivatives have their usual definitions: `F=dA+g A wedge A`, and `mathcal D` includes the spin connection and the displayed representations. Gauge couplings, Yukawa matrices and Higgs parameters are independently calibrated inputs. The strong theta term is set to zero as an explicit matching choice in this candidate; it is not derived. The hypercharge assignment and ordinary gauge anomaly cancellations are imported. Dirac neutrinos are a model choice consistent with the imposed EFT `B−L`; their nature has not been established by this calculation. No claim to derive Standard Model content or masses is made.

## 3. Quantum and EFT completion rule

The finite action above is the leading truncation, not a quantum-closed fundamental action. Define the EFT through the following explicit prescription:

1. Use a fixed covariant gauge, dimensional regularization, and modified minimal subtraction for perturbative calculations; physical outputs must be gauge independent at the claimed order.
2. Include every local counterterm consistent with diffeomorphism, Standard Model gauge invariance, and the stated global/discrete symmetries. This includes curvature-squared terms, nonminimal terms `xi_phi R phi²`, `xi_X R|X|²`, `xi_D R D²`, `xi_H R H†H`, a bare `m_X0²|X|²`, the `|X|²D²` portal generated by phi loops, and mixed scalar/Higgs portals.
3. At the benchmark matching scale `mu_0=g_X v`, set the renormalized `xi` coefficients, `m_X0²`, mixed portals absent from `U`, curvature-squared coefficients, and higher-dimension Wilson coefficients to zero. Their counterterms and running away from `mu_0` are retained. These zeros are matching choices, not symmetry protection or naturalness results.
4. At any finite EFT order, retain all operators and loop terms contributing at that order; higher-derivative terms are expanded perturbatively, not resummed into extra exact ghost degrees of freedom.
5. Restrict use to momenta, field amplitudes, and curvature below the cutoff `Lambda_EFT`; the illustrative prescription is `Lambda_EFT=0.1 M_P`. Sensitivity to unknown alternative matching data is a physical limitation of the proposal. No all-orders convergence or ultraviolet completion is asserted.

This is one formal EFT matching prescription, rather than an unspecified choice of arbitrary coefficients after data. The session executes only the leading classical scalar and linearized-gravity terms, with a weak-reporter calculation. Loop and higher-order quantum numerical predictions remain open.

## 4. The theory datum

\[
\mathcal T_{\rm GRUT1}=(\mathcal H_{\rm EFT},\mathfrak A_{\rm phys},
\mathcal C,H_{\rm ADM/EFT},\mathcal S_0,\mathcal O,\Theta).
\]

### State space and observable algebra

Within perturbation theory about the specified asymptotically flat vacuum, use an infrared-regulated wave-packet Fock space of the two physical graviton helicities and elementary matter fields as a perturbative starting construction. Massless photons/gravitons require infrared dressing or inclusive scattering probabilities when the regulator is removed. The full nonperturbative QCD state space and bound-state physics are imported as a Standard Model obligation, not constructed by this elementary-field Fock prescription. In covariant gauge include the gauge/ghost complex and take physical states in BRST cohomology order by order. This is a perturbative definition; it is not an exact nonperturbative Hilbert space over every geometry.

`mathfrak A_phys` consists of gauge-invariant scattering observables, asymptotic charges, and relational or gravitationally dressed apparatus observables, constructed to the same EFT order. Products of matter fields at coordinate labels are computational representatives, not exact local Dirac observables in nonlinear gravity. In the weak-field apparatus limit, locally smeared fields and detector probabilities approximate these physical quantities.

The complex scalar charge is

\[
j_X^\mu=i(X\nabla^\mu X^*-X^*\nabla^\mu X),\qquad
Q_X=-\int_\Sigma d^3x\sqrt h\,n_\mu j_X^\mu,
\]

where `n` is the future-directed unit normal. This convention gives `j^0=omega f²` and positive charge for the displayed rotating soliton in flat coordinates. Charge conservation follows from the postulated `U(1)_X`, not from a readout convention.

### Constraints and evolution

Classical states satisfy the Einstein Hamiltonian/momentum constraints, all Standard Model Gauss constraints, the boundary conditions, and the matter equations. In ADM variables `h_ij,pi^ij`, lapse `N`, shift `N^i`,

\[
H=H_{\rm ADM}+\int d^3x\,[N\mathcal H_\perp+N^i\mathcal H_i
+\lambda^a\mathcal G_a],
\]

\[
\mathcal H_\perp=
\frac{2}{M_P^2\sqrt h}(\pi^{ij}\pi_{ij}-\pi^2/2)
-\frac{M_P^2\sqrt h}{2}\,{}^{(3)}R+\mathcal H_{\perp,\rm matter}=0,
\]

\[
\mathcal H_i=-2h_{ik}D_j\pi^{jk}
+p_\phi\partial_i\phi+p_D\partial_iD
+p_X\partial_iX+p_{X^*}\partial_iX^*
+\mathcal H_{i,\rm SM}=0.
\]

The new scalar part of the matter Hamiltonian constraint is

\[
\mathcal H_{\perp,\rm scalars}=
\frac{p_\phi^2+p_D^2}{2\sqrt h}+\frac{p_Xp_{X^*}}{\sqrt h}
+\sqrt h[\tfrac12h^{ij}(\partial_i\phi\partial_j\phi+
\partial_iD\partial_jD)+h^{ij}\partial_iX^*\partial_jX+U].
\]

The quantum evolution used for calculations is the gauge-fixed perturbative Hamiltonian/S-matrix with the physical-state constraint prescription above. An exact quantum implementation of the full ADM constraint algebra is not delivered here.

**Time is supplied at the framework level:** the manifold is time-oriented, and ADM evolution uses a chosen foliation. Coordinate time is gauge dependent; the worked asymptotically flat benchmark uses the asymptotic proper-time parameter. This is not a derivation of time or its arrow.

## 5. Coupled field equations and gravity

At leading order,

\[
M_P^2G_{\mu\nu}+\rho_\Lambda g_{\mu\nu}
=T_{\mu\nu}^{\rm SM}+T_{\mu\nu}^{\phi,X,D},\qquad
T_{\mu\nu}=-\frac{2}{\sqrt{-g}}\frac{\delta S_{\rm matter}}{\delta g^{\mu\nu}}.
\]

\[
\Box\phi-\lambda_\phi(\phi^2-v^2)\phi-2g_X^2\phi|X|^2-h_D^2\phi D^2=0,
\]

\[
\Box X-g_X^2\phi^2X-2\lambda_X|X|^2X=0,
\]

\[
\Box D-(m_{D0}^2+h_D^2\phi^2)D-\lambda_DD^3=0.
\]

No fundamental semiclassical rule `G = <T>` is added. The classical/semiclassical metric used in the worked example is a leading saddle of the same metric-plus-matter EFT. Branch-dependent metric fluctuations, graviton loops, and quantum back-reaction beyond this saddle remain uncalculated.

All species share the metric characteristic cone through their leading kinetic operators. Massive group velocities need not be equal. Nonminimal and higher-order EFT terms require their own propagation analysis; no exact all-frequency universal cone beyond the stated leading regime is proved.

A freely falling sufficiently weak localized test body obeys the geodesic equation at the monopole/test limit, a clock measures `d tau²=−g_mu nu dx^mu dx^nu`, and electromagnetic rays follow the metric null characteristics at geometric-optics order. Scalar reporter D also feels its explicit organization-dependent mass; its trajectory cannot be called pure gravitational free fall inside a varying phi background.

## 6. Allowed preparation and physical organization

`mathcal S_0` contains finite-energy, constraint-satisfying EFT states in the chosen asymptotic vacuum. The worked subclass is a large-charge, approximately spherical scalar preparation near a nodeless rotating solution, plus a weak D probe and an ordinary prepared readout apparatus.

At classical level,

\[
X(t,r)=e^{-i\omega t}f(r)/\sqrt2,\quad \phi(t,r)=\varphi(r),\quad D=0.
\]

At quantum level one may use a large-mean-charge semiclassical state with charge fluctuations, or project onto a fixed integer charge sector. A fixed-charge state has no nonzero one-point `X` expectation; the rotating classical profile represents its leading large-charge description. Neither alternative is an exact all-orders coherent quantum solution.

The energy and charge of the prepared matter are supplied initial resources. The calculation does not create charge from the vacuum, convert absolute nothing into matter, or prove autonomous formation from a diffuse generic seed. No reservoir or controller acts to sustain the stationary leading configuration after preparation. Apparatus operation and probe injection remain supplied interventions whose energy is accounted separately.

Organization occurs through feedback in the equations: a high X density lowers the preferred magnitude of phi; smaller phi lowers the X mass; this can bind the charge and sustain the altered field profile. The observable `v^−2` times a physically smeared `phi²` is an occupation-correlated structural record. In a quantum calculation use the renormalized composite `[phi²]_R` with a stated vacuum/reference and apparatus calibration, not an unregulated pointwise product. Persistence after removal of X, classical error correction, learning, consciousness, and universal self-organization are not established.

## 7. Measurements and calibration

For any finite apparatus record `r`, use

\[
\mathcal I_r(\rho)=\operatorname{Tr}_{A}
[(I\otimes P_r)U(\rho\otimes\sigma_A)U^\dagger(I\otimes P_r)],
\quad p_r=\operatorname{Tr}\mathcal I_r(\rho).
\]

The apparatus state `sigma_A`, interaction region and final record projectors are physical preparation/readout data. A fully implemented apparatus must obtain its interaction Hamiltonian and unitary `U` from this same action with prepared apparatus/control degrees of freedom, respecting the constraints and total resource accounting; an arbitrary extra detector coupling cannot be chosen for free. The ideal instrument limit is adopted for the worked readout. The session does not derive a macroscopic detector from the full Standard Model. The reporter propagation phase is calculated in `WORKED_REGIME.md`. Its conversion to a probability uses a separately specified ideal interferometric instrument. An action-compatible microscopic source, beam splitter, mode recombination and detector have not been constructed; their laboratory realization and costs remain an explicit operational gap.

`Theta` includes `M_P,rho_Lambda`, ordinary matter parameters, `v,g_X,lambda_phi,lambda_X,m_D0,h_D,lambda_D`, the matching prescription and scale, apparatus calibration, and the prepared charge/state. New-sector masses and local scattering amplitudes calibrate its couplings; state charge is measured independently of the observables used for the held-out test. Predictions are not independently calibrated observations merely because numerical values are printed.

The benchmark uses `g_X=.1`, `lambda_phi=.01`, `lambda_X=.001`, `m_D0=g_X v`, `h_D=.01`, `lambda_D=.001`, and `v/M_P=10^−5`. These are a declared illustrative point, not fitted or measured parameters of nature. The value of the prepared charge is specified by the numerical state in the execution records.

## 8. Status and companion files

The final integrated requirement/status table is in `COMPATIBILITY.md`. The worked regime, execution results, and limits are in `WORKED_REGIME.md`; calibration and rejection are in `PREDICTIONS_AND_TESTS.md`; outstanding claims are in `OPEN_CLAIMS.md`. The code implements the stated leading equations and does not certify their quantum or empirical validity.

Known mechanisms are adopted deliberately: the scalar organization sector is a quartic extension of the Friedberg–Lee–Sirlin/non-topological-soliton construction, and gravity/ordinary matter are imported. Candidate 1 is a proposed synthesis with new-sector parameters, not an originality claim for GR, Q-balls, quantum mechanics, or scalar detectors.

