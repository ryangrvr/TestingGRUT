# Predictions, calibration, and rejection

**GRUT Candidate 1 · 18 October 2026**

This document separates calculated synthetic outputs from future empirical tests. The benchmark couplings are declared construction inputs, not independently measured parameters. No dataset was supplied for empirical calibration or validation.

## 1. What must be calibrated first

The main candidate fixes the action, symmetries, matching scale/rules and an illustrative cutoff. A physical realization would require independent measurements of the following quantities before a held-out comparison.

| Input | Independent calibration route | What that alone does not determine |
|---|---|---|
| `M_P`, ordinary matter constants | Ordinary gravitational and particle measurements in the adopted framework | New scalar couplings |
| `rho_Lambda` and boundary state | Cosmological/background calibration appropriate to the experiment | A claim that the flat benchmark describes cosmology |
| `v,lambda_phi,g,lambda_X` | Scalar pole masses plus independent scattering amplitudes and mediator response, in a declared renormalization scheme | Vacuum masses alone leave parameter degeneracies |
| `m_D0,h_D,lambda_D` | Reporter mass in independently characterized backgrounds and reporter scattering | Target profile or ideal instrument efficiency |
| Matching coefficients | Independent operator-sensitive data at a specified EFT order, or retain the declared benchmark prescription | Symmetry does not protect the chosen zero portals/masses/xi values |
| Prepared state/charge | Charge and spatial-state characterization independent of reporter phases to be tested | Generic formation probability, exact lifetime |
| Apparatus | Source spectrum, arm phase, selection, visibility and record efficiency | A computed microscopic realization from the present package |

Calibration must not adjust the target profile, couplings or arm phase separately for each reported momentum. Once a physical Q and parameters are specified, solve the same stationary branch and propagate all measured input uncertainties. Numerical central amplitude was used to find this session's preparation; Q is an output of that search. In a charge-controlled experiment Q would be the input and the boundary-value solve would have to impose it directly.

## 2. Calculated benchmark targets

| Output | Candidate prediction at the declared preparation | Role |
|---|---:|---|
| `E/(m_X Q)` | 0.91108652 | Binding-energy benchmark |
| `R90 g v` | 5.1034322 | Matter profile benchmark |
| `phi(0)/v` | 0.049079979 | Organization response |
| `delta_0(k=3 gv)` | 0.0067765110 radians | Reporter momentum dependence |
| `delta_0(k=5 gv)` | 0.0040665116 radians | Reporter momentum dependence |
| `delta_0(k=7 gv)` | 0.0029047695 radians | Reporter momentum dependence |
| Central clock redshift | 6.68813×10^−10 | Full-stress weak-gravity benchmark |
| Light bending, `b g v=20` | 2.89601×10^−10 radians | Exterior metric consistency |
| Selected packet output at `tau=40`, finest grid | `p_select≈.992`, joint `P+≈.492`, conditional `P+≈.49591` | Approximate finite-bandwidth instrument example |

The momentum set and couplings were written before reporter evaluation. They are reproducibility targets, not empirical held-out successes. Every value must be read with the approximation/error inventory below. The first six rows could provide correlated new-sector tests if such particles, states and instruments were physically realized. Exterior Einstein lensing and Born probabilities reproduce adopted laws; they do not distinguish GRUT from those same laws.

## 3. Comparison models

| Baseline | Predicted difference | Purpose and limit |
|---|---|---|
| Uniform phi, no lump | `V_D=0`, reporter phase zero | Tests the structure-dependent interaction; does not model all backgrounds |
| Rigid `phi=v`, same positive X quartic | No mediator-induced attraction; leading scalar energy is at least `m_X Q` | Isolates the feedback needed for this nongravitating bound configuration |
| Original FLS potential `lambda_X=0` | Different energy/charge and large-charge behavior | Tests the additional self-interaction after independent calibration; no old numerical bound transfers |
| Density-only Newtonian source | Sets interior `Phi=Psi` and misses pressure-sensitive clock response | Internal diagnostic; not the candidate's Einstein prescription |
| Einstein + ordinary matter without new scalar sector | No corresponding X/phi/D object or reporter channel | A physical comparison requires the new sector to exist and be accessible |

Do not claim a novel modification of gravity from a scalar source obeying standard Einstein equations. Do not claim a universal new prediction merely by giving a new name to a familiar scalar soliton.

## 4. Numerical rejection rules used now

These acceptance checks assess the implemented leading regime, not nature. `code/verify.py` records the actual pass/fail results.

| Claim | Reject the executed calculation if… |
|---|---|
| Nontrivial stationary solution | A refinement fails, gives a vacuum/nodal inadmissible profile, has bad boundary falloff, or fails the radial residual tolerance |
| Resolved global outputs | The last two domains disagree in Q or E by more than 10^−6 relative, or in `E/(m_X Q)` by more than 10^−7 |
| Classical binding comparison | `E/(m_X Q)>=1` at the same couplings/state |
| Conserved stress | Dimensionless stress residual exceeds 10^−7 or relative virial error exceeds 10^−7 |
| Weak-gravity expansion | Potentials/compactness are not small; the present implementation threshold is 10^−4 |
| Metric source/correspondence | Active mass and energy disagree by more than 10^−7 relative, or the exterior bending benchmark fails at the declared numerical tolerance |
| Tested local amplitude minimum | A nonsymmetry negative mode persists under mesh/domain refinement; treating the observed translation artifact as harmless requires its identified symmetry overlap and h² scaling |
| Phase calculation | Independent wave/variable-phase methods disagree above 10^−7 radians or the exact-profile Born bound is violated beyond numerical errors |
| Finite packet bookkeeping | Norm fails 10^−10 conservation; output probabilities leave [0,1], do not sum to selection, or visibility exceeds one beyond roundoff |

Thresholds were consolidated for final reporting after execution and are not statistical preregistration. Direct numerical outputs, failed starts and full prescriptions are retained so these checks can be challenged. Eigenvalue convergence is assessed on actual spectra, not inferred from code assertions alone.

## 5. Future tests that can reject central claims

1. **Global bound state:** compute the fixed-charge minimum across competing branches and all relevant fragment partitions. Reject global stability if a lower-energy allowed partition is found. A local Hessian or `E<m_XQ` cannot override that result.
2. **Nonlinear persistence:** evolve perturbed three-dimensional configurations with conserved charge/energy and resolved radiation boundaries. Reject any proposed robust operating class if perturbations in that class reproducibly disperse or fragment it. Specify amplitude, duration and mesh convergence before claiming a lifetime.
3. **Quantum persistence:** calculate renormalized fluctuations and allowed decay channels of this quartic model. Reject a claimed semiclassical long-lived regime if corrections remove binding, add an instability, or produce a lifetime shorter than its stated operational time. Results for a different FLS potential do not settle this test.
4. **Self-consistent gravity:** solve the Einstein–scalar constraints/profile and perturbations beyond linear response. Reject the claimed approximation where corrections are not controlled. Strong-field extrapolation of the flat solution is prohibited.
5. **Physical observation:** exhibit an action-compatible source, arm splitter/recombiner and detector, including losses, recoil and energy accounting. Reject an unconditional instrument claim if its required coupling or record dynamics is absent. Until then, the reported readout remains conditional on an ideal instrument.
6. **Empirical discrimination:** after independent calibration, compare reporter phases and profile/energy predictions across held-out momenta/charges without refitting. Use a registered covariance/error model and competing explanations. A statistically significant systematic disagreement outside the complete error budget rejects that parameterized EFT/preparation claim; it does not justify adding a fitted mechanism without a new version.

No physical discovery or experimental result is asserted by this list. It defines how the proposal can lose.

## 6. Error inventory

| Error/source | What is controlled now | What is not quantified |
|---|---|---|
| Stationary numerical solution | Four tolerance/domain runs, residuals, identities | Rigorous outward-rounded solution enclosure |
| Local fluctuation spectrum | Three grids, known translation/phase zero modes, tested angular sectors | Rigorous global spectral theorem; nonlinear evolution |
| Reporter phase | Two independent ODE forms, tolerance refinement, analytic Born remainder for exact profile | Complete interval error including profile interpolation |
| Finite packet | Three grids; norm/selection/visibility bookkeeping | Rigorous continuum error or packet domain refinement; fine digits are not claimed |
| Weak gravity | Compactness/potentials about 10^−9; all leading stresses retained | Nonlinear correction coefficient and gravitating branch response |
| Quantum/EFT truncation | Explicit weak-coupling/subcutoff regime and counterterm prescription | Loop shifts, running amplitudes and lifetime at this state |
| Preparation/probe | State and small test-probe hierarchy explicitly supplied | Laboratory acceptance, target recoil/deformation and inelastic loss |
| Empirical/calibration | None; all current numbers are synthetic | Measured inputs, physical systematic errors, empirical confidence |

Tree loop-counting factors such as `g²/(16 pi²)≈6.3×10^−5` suggest a perturbative regime but are not a computed uncertainty on a soliton observable. Small gravity and a well-converged classical solver do not make all omitted physical effects smaller than the printed numerical digits.
