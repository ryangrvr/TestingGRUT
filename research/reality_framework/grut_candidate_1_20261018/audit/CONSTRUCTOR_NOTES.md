# Constructor notes: GRUT Candidate 1, covariant architecture B

Date: 18 October 2026. Internal construction/review notes for the integrating worker. These notes do not replace the five requested deliverables, do not modify TestingGRUT, and make no empirical validation claim.

## Assessment

The selected architecture is a coherent narrower **effective model proposal**: Einstein gravity, the Standard Model with three Dirac neutrinos, and a specified scalar sector supporting self-bound charge configurations and a local reporter. It adopts the structures that the previous microscopic route did not produce. A real scalar profile correlates with charge occupation because both fields satisfy coupled local equations. This is a physical organization/occupation record in a restricted sense. It is not a general theory of memory, measurement outcomes, consciousness, dimension, the arrow of time, or emergent gravity.

Its main scalar mechanism is a quartic extension of the Friedberg–Lee–Sirlin/non-topological-soliton mechanism. GR, ordinary quantum matter, quantum probabilities, and detector instruments are imported. Specifying these together can produce a calculable candidate without constituting an original fundamental-physics derivation. The conceptual viability of this construction must be separated from existence and stability of the **chosen numerical solution**, accuracy of its approximations, and empirical viability of its new particles. The numerical worker must establish the first two; neither proves the third.

## 1. Covariant action and explicit symmetries

Use natural units, metric signature \((-+++)\), an oriented, time-oriented globally hyperbolic four-manifold, and the reduced Planck mass \(M_P^{-2}=8\pi G\). At leading order,

\[
S_0=\int d^4x\sqrt{-g}\left[\frac{M_P^2}{2}R-\rho_\Lambda
+\mathcal L_{\rm SM+3\nu_R}
-\frac12(\nabla\phi)^2-\nabla_\mu X^*\nabla^\mu X
-\frac12(\nabla D)^2-U\right]+S_{\rm GHY,ref},
\]

\[
U=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
+g_X^2\phi^2|X|^2+\lambda_X|X|^4
+\frac12(m_{D0}^2+h_D^2\phi^2)D^2+\frac{\lambda_D}{4}D^4.
\]

All displayed scalar couplings are positive. \(S_{\rm GHY,ref}\) is the usual boundary term and asymptotic reference prescription needed for the chosen variational problem/ADM charge. The new fields are Standard Model gauge singlets. Impose exact \(U(1)_X\), \(\phi\mapsto-\phi\), \(D\mapsto-D\) **as EFT symmetries**. Impose unbroken global \(B-L\) on the ordinary matter sector, including three right-handed neutrinos, to exclude Majorana and Weinberg mass operators. No ultraviolet or nonperturbative quantum-gravity protection of these global symmetries is supplied.

The asymptotically flat worked boundary condition is \(\phi\to+v\), \(X,D\to0\), the ordinary electroweak vacuum, and \(\rho_{\Lambda,R}=0\). Choice of the positive scalar vacuum, time orientation, and asymptotic flatness is supplied boundary/preparation data. \(\rho_{\Lambda,R}=0\) is not a prediction for our cosmology. The two disconnected \(\phi\) vacua admit domain walls; none are needed in the worked regime.

### Ordinary matter

The full renormalizable ordinary sector is \(SU(3)_c\times SU(2)_L\times U(1)_Y\) with three copies of

| Field | Representation \((SU(3),SU(2),Y)\) |
|---|---|
| \(Q_L\) | \((3,2,1/6)\) |
| \(u_R\) | \((3,1,2/3)\) |
| \(d_R\) | \((3,1,-1/3)\) |
| \(L_L\) | \((1,2,-1/2)\) |
| \(e_R\) | \((1,1,-1)\) |
| \(\nu_R\) | \((1,1,0)\) |
| \(H\) | \((1,2,1/2)\), one Higgs doublet |

Its terms are gauge kinetic terms, positive-energy fermion kinetic terms with the tetrad/spin connection, \(-(D_\mu H)^\dagger D^\mu H\), the Higgs potential, and

\[
-[\bar Q_LY_dHd_R+\bar Q_LY_u\widetilde H u_R
+\bar L_LY_eHe_R+\bar L_LY_\nu\widetilde H\nu_R+{\rm h.c.}].
\]

Gauge couplings, Yukawa matrices, Higgs parameters, and CP parameters are imported calibrated constants. The new model derives none of them. Dirac neutrinos are a choice, not a measured conclusion of this model. If using the written fermion convention \(i\bar\psi\gamma^ae_a^\mu\mathcal D_\mu\psi\) with \((-+++)\), specify a consistent Clifford convention, e.g. \(\{\gamma^a,\gamma^b\}=-2\eta^{ab}\). A strong theta term must either be included with its coefficient or explicitly set to zero matching; the integrator chose the latter. A future empirical claim about strong CP must then acknowledge this assumption.

Classical positivity of the new potential does not establish the all-scale stability of the full running Standard Model Higgs potential. Use the prepared electroweak vacuum within its applicable lifetime and field range. This demonstration does not calculate Standard Model vacuum decay.

## 2. Equations, signs, charge, and positive new-sector energy

At the displayed order,

\[
M_P^2G_{\mu\nu}+\rho_\Lambda g_{\mu\nu}
=T^{\rm SM}_{\mu\nu}+T^{\phi,X,D}_{\mu\nu},
\]

\[
\Box\phi-\lambda_\phi(\phi^2-v^2)\phi
-2g_X^2\phi|X|^2-h_D^2\phi D^2=0,
\]

\[
\Box X-g_X^2\phi^2X-2\lambda_X|X|^2X=0,
\qquad
\Box D-(m_{D0}^2+h_D^2\phi^2)D-\lambda_DD^3=0.
\]

The new scalar stress is

\[
T^{\phi,X,D}_{\mu\nu}=
\nabla_\mu\phi\nabla_\nu\phi+\nabla_\mu D\nabla_\nu D
+\nabla_\mu X^*\nabla_\nu X+\nabla_\nu X^*\nabla_\mu X
-g_{\mu\nu}\left[\tfrac12(\nabla\phi)^2+\tfrac12(\nabla D)^2
+\nabla X^*\cdot\nabla X+U\right].
\]

One transparent charge convention is

\[
j_X^\mu=i(X\nabla^\mu X^*-X^*\nabla^\mu X),\quad
\nabla_\mu j_X^\mu=0,\quad
Q=-\int_\Sigma d^3x\sqrt h\,n_\mu j_X^\mu,
\]

where \(n^\mu\) is the future unit normal. For \(X=e^{-i\omega t}f(r)/\sqrt2\), \(\omega>0\), flat space gives \(j_X^0=\omega f^2\), hence \(Q=4\pi\omega\int_0^\infty dr\,r^2f^2>0\). Changing the signs of both current and surface orientation is also consistent, but a single convention must be used in text and code.

For the stationary \(D=0\) configuration,

\[
E=4\pi\int_0^\infty dr\,r^2\left[
\tfrac12\varphi'^2+\tfrac12 f'^2+\tfrac12\omega^2 f^2
+\tfrac{\lambda_\phi}{4}(\varphi^2-v^2)^2
+\tfrac12g_X^2\varphi^2f^2+\tfrac14\lambda_Xf^4\right].
\]

Every term is nonnegative. A nonzero finite solution therefore has strictly positive scalar energy. Binding means \(E<m_XQ\), not negative total energy. A positive-energy classical saddle is not automatically a quantum-stable state or a global fixed-charge minimum.

Vacuum masses are \(m_X=g_Xv\), \(m_\phi=\sqrt{2\lambda_\phi}v\), and \(m_D=\sqrt{m_{D0}^2+h_D^2v^2}\), at tree level with the specified matching conditions.

## 3. Nontrivial binding condition and normalized equations

Positive couplings alone do **not** ensure the desired Q-ball. Let \(y=|X|^2\), set \(D=0\), and minimize the potential over \(s=\phi^2\ge0\):

\[
s_*(y)=\max\left(v^2-\frac{2g_X^2}{\lambda_\phi}y,0\right).
\]

For \(y\le\lambda_\phi v^2/(2g_X^2)\),

\[
\frac{U_{\min}(y)}y=g_X^2v^2+
\left(\lambda_X-\frac{g_X^4}{\lambda_\phi}\right)y.
\]

An energy-per-charge advantage over dilute particles requires

\[
g_X^4>\lambda_\phi\lambda_X.
\]

In that case the homogeneous interior minimum occurs at

\[
\phi=0,\quad y_0=\frac{v^2}{2}\sqrt{\frac{\lambda_\phi}{\lambda_X}},\quad
\omega_0^2=\min_{y>0}\frac{U_{\min}(y)}y
=v^2\sqrt{\lambda_\phi\lambda_X}<m_X^2.
\]

The familiar flat-space frequency range is \(\omega_0<\omega<m_X\). It describes a candidate branch, not proof that every charge or every shooting guess has a finite bound configuration. At equality \(g_X^4=\lambda_\phi\lambda_X\), this strict bulk binding advantage disappears.

Set

\[
x=g_Xvr,\quad h=\varphi/v,\quad F=f/v,\quad
\Omega=\omega/(g_Xv),\quad a=\lambda_\phi/g_X^2,\quad b=\lambda_X/g_X^2.
\]

The exact leading flat radial equations are

\[
h''+\frac2x h'=a h(h^2-1)+hF^2,
\qquad
F''+\frac2x F'=(h^2-\Omega^2)F+bF^3,
\]

with \(h'(0)=F'(0)=0\), \(h(\infty)=1\), \(F(\infty)=0\), and regular finite central values. Primes here mean \(x\) derivatives. Fixing \(F(0)=2.5\) selects a requested shooting condition; it is not a dynamical law or generic-state claim.

Define

\[
I=\int_0^\infty dx\,x^2F^2,
\]

\[
\varepsilon=\int_0^\infty dx\,x^2\left[
\tfrac12h'^2+\tfrac12F'^2+\tfrac12\Omega^2F^2
+\tfrac a4(h^2-1)^2+\tfrac12h^2F^2+\tfrac b4F^4\right].
\]

Then \(E=(4\pi v/g_X)\varepsilon\), \(Q=(4\pi/g_X^2)\Omega I\), and \(E/(m_XQ)=\varepsilon/(\Omega I)\). These factors provide direct normalization checks between analytic notes and code.

At the frozen illustrative point, \(a=1,b=0.1\), \(\Omega_0=(ab)^{1/4}=0.5623413252\), and the infinite-bulk \(\Omega_0\) core has \(F=(a/b)^{1/4}=1.7782794100\). A finite configuration need not have either of those central/bulk values.

## 4. Constraints and scope of the quantum construction

In ADM variables with Einstein action \(M_P^2R/2\), the gravitational Hamiltonian constraint is

\[
\mathcal H_\perp=\frac{2}{M_P^2\sqrt h}
(\pi^{ij}\pi_{ij}-\tfrac12\pi^2)
-\frac{M_P^2\sqrt h}{2}\,{}^{(3)}R
+\mathcal H_{\perp,\rm matter}=0.
\]

The momentum constraint is

\[
\mathcal H_i=-2h_{ik}D_j\pi^{jk}
+p_\phi\partial_i\phi+p_D\partial_iD
+p_X\partial_iX+p_{X^*}\partial_iX^*
+\mathcal H_{i,\rm SM}=0.
\]

The new matter contribution is

\[
\mathcal H_{\perp,\rm new}=
\frac{p_\phi^2+p_D^2}{2\sqrt h}+\frac{p_Xp_{X^*}}{\sqrt h}
+\sqrt h\left[\tfrac12h^{ij}(\partial_i\phi\partial_j\phi+
\partial_iD\partial_jD)+h^{ij}\partial_iX^*\partial_jX+U\right].
\]

The full Hamiltonian includes lapse/shift multipliers, all Standard Model Gauss constraints and their multipliers, and the asymptotic ADM surface Hamiltonian. Gauge/spinor formulations may additionally impose the usual nondynamical and second-class constraints before the reduced description. Initial data must satisfy these constraints. Flat radial data followed by a computed weak metric is only a perturbative construction of such initial data; it is not exact nonlinear constraint satisfaction unless those constraints are solved to the claimed order.

Use a fixed covariant gauge and perturbative BRST cohomology for physical state representatives at finite EFT order. The formal gauge-fixed matter/graviton oscillator/Fock construction is a computational starting point with a regulator. QED/gravity require infrared dressing or inclusive observables, while confined quarks/gluons are not exact asymptotic physical particles. Hadronic observables require ordinary low-energy matching. No exact nonperturbative Hilbert space across arbitrary geometries, complete nonlinear Dirac-observable algebra, or quantum ADM constraint representation is produced.

At EFT scope the physical algebra comprises gauge-invariant scattering/relational apparatus observables, asymptotic charges, and gravitationally dressed matter operators to the chosen order. Strict compact-region tensor factors and exactly commuting local gravitational subalgebras are not assumed. Coordinate-labelled \(\phi(x)\), \(X(x)\) are useful gauge-fixed representatives, not exact diffeomorphism-invariant observables.

The quantum/Born framework is **postulated**: positive normalized states on the applicable algebra, standard unitary perturbative dynamics, and completely positive instruments. The classical Q-ball can approximate a large-charge state. A fixed charge eigenstate has \(\langle X\rangle=0\); its stationary charge density and correlators, rather than a literal rotating one-point function, correspond to the classical background. A coherent large-charge representative has charge fluctuations. Charge projection, loop corrections, and finite-charge corrections are distinct calculation obligations.

## 5. Local reporter, records, instruments, and clocks

Because \(D=0\) on the background, its added action leaves the leading \((\phi,X)\) radial solution unchanged. A weak incoming \(D\) wave obeys

\[
[\Box-m_{D0}^2-h_D^2\varphi^2(r)]D=0.
\]

In flat space take an asymptotic frequency \(\nu\) and momentum \(k^2=\nu^2-m_D^2>0\). Each partial wave satisfies

\[
u_\ell''+\left[k^2-\frac{\ell(\ell+1)}{r^2}
-h_D^2(\varphi^2(r)-v^2)\right]u_\ell=0.
\]

Thus the readout interaction is specified by the same local action, with no extra record-reading force. The observable phase shift is measured against a physically prepared reference, not an absolute unobservable phase. For \(h_D/g_X=0.1\), the dimensionless potential is \(0.01(h^2-1)\). Its negative sign in the suppressed core is an attractive mass-squared well. Since \(m_{D0}^2>0\), the leading \(D\) quadratic operator has no negative-frequency-squared instability on this static scalar background, even if it admits positive-frequency bound levels. Probe amplitude, pulse shape, and apparatus bandwidth are preparation/calibration inputs. Backreaction scales with probe energy/occupation and must be checked if the weak limit is relaxed.

Classically, a smeared normalized \(\varphi^2/v^2\) diagnoses organization. At quantum order \([\phi^2]_R\) is a renormalized composite operator and can mix with other permitted dimension-two scalar operators and the identity/curvature. A calibration must define its scheme and vacuum reference. The physically read quantity is the reporter response to the renormalized action, not an uncalibrated coordinate field value. This issue does not alter the displayed tree-level reporter calculation.

For a finite apparatus with state \(\sigma_A\), a local interaction generating \(U\), and complete record projectors \(P_r\),

\[
\mathcal I_r(\rho)=\operatorname{Tr}_A\{(1\otimes P_r)
U(\rho\otimes\sigma_A)U^\dagger(1\otimes P_r)\},\qquad
p_r=\operatorname{Tr}\mathcal I_r(\rho).
\]

The outcome probabilities and conditional-state rule are the adopted quantum framework. A fully microscopic Standard Model detector is not calculated by writing this instrument formula. Strict system-apparatus tensor products are ideal weak-gravity approximations, with the gravitational dressing caveat above. Preparation of the apparatus, reference arm, pointer stability and decoherence must be supplied or separately modeled. No theorem selecting one outcome from a unitary state is claimed.

The scalar profile is an **occupation-correlated record**, not established permanent memory. Removal of charge removes the stabilizing condition for the suppressed core; a finite local domain of the opposite Z2 vacuum is also not automatically stable. Neither the double-well drawing nor charge stability proves autonomous formation from a generic seed, self-repair, replication, or learning.

Coordinate time is a foliation/chart choice; asymptotically flat energy and stationary frequencies are referred to asymptotic proper time. An actual clock follows \(d\tau^2=-g_{\mu\nu}dx^\mu dx^\nu\) and uses ordinary matter transition/oscillation phase calibrated to \(\tau\). No scalar is assumed to provide a globally monotonic internal clock. The rotating \(X\) phase can serve as a local phase reference only with an operational comparison and its quantum limitations specified. This construction postulates temporal geometry and does not derive time or an arrow from records.

## 6. Gravity consistency checks

For a stationary spherical solution write

\[
ds^2=-\alpha(r)^2dt^2+a(r)^2dr^2+r^2d\Omega^2,
\qquad a^2=(1-2Gm(r)/r)^{-1}.
\]

At \(D=0\),

\[
\rho=\frac{\varphi'^2+f'^2}{2a^2}
+\frac{\omega^2f^2}{2\alpha^2}+U,
\quad
p_r=\frac{\varphi'^2+f'^2}{2a^2}
+\frac{\omega^2f^2}{2\alpha^2}-U,
\]

\[
m'=4\pi r^2\rho,\qquad
\frac{\alpha'}\alpha=
\frac{G(m+4\pi r^3p_r)}{r^2(1-2Gm/r)}.
\]

These formulas show why simply inserting energy density into a Newtonian potential inside relativistic scalar material can miss stress terms. The exterior asymptotic \(1/r\) mass and clock/lensing limits must agree with the ADM energy at the computed order. Interior metric potentials require the relevant stress components or the consistent linearized Einstein equations. The calculator owns the actual weak-field treatment and error estimate.

For the benchmark \(v/M_P=10^{-5}\), gravity is parametrically weak at moderate dimensionless charge/radius, but the compactness \(2Gm(r)/r\) must actually be evaluated; the small value of \(v/M_P\) alone is not a uniform bound at arbitrarily large charge. Free scalar gradients and the common metric produce the same leading principal light cone. Organization-dependent scalar forces and higher-derivative EFT corrections must not be called a proof of identical trajectories or an exact all-frequency cone for every species.

No separate fundamental law \(G_{\mu\nu}=\langle T_{\mu\nu}\rangle/M_P^2\) is introduced. The displayed classical geometry is a leading large-occupation saddle of the same EFT. Stress fluctuations, branch-dependent metrics, quantum backreaction and graviton loops beyond that saddle are uncalculated.

## 7. Parameter, operator, and matching ledger

Matching values refer to \(\mu_0=g_Xv\), for example in dimensional regularization and \(\overline{\rm MS}\) with a stated tadpole convention. Tree-level \(v\) and mass relations receive radiative corrections; they are not automatically on-shell definitions.

| Datum/operator | Frozen matching/preparation value | Status |
|---|---|---|
| \(M_P\) | Reduced Planck mass, independently calibrated if applying to nature | Imported gravitational parameter |
| \(\rho_{\Lambda,R}\) | 0 for the asymptotically flat worked regime | Boundary-sector choice; not a cosmological prediction |
| \(v/M_P\) | \(10^{-5}\) | Illustrative dimensionful scalar scale |
| \(g_X\) | 0.1 | Illustrative scalar coupling |
| \(\lambda_\phi,\lambda_X\) | 0.01, 0.001 | Illustrative values; satisfy strict binding inequality |
| \(m_{D0},h_D,\lambda_D\) | \(g_Xv,0.01,0.001\) | Explicit reporter parameters |
| \(m_{X0,R}^2|X|^2\) | \(m_{X0,R}^2=0\) at \(\mu_0\) | Allowed relevant operator; tuned matching input, not forbidden |
| \(\phi^2H^\dagger H, |X|^2H^\dagger H,D^2H^\dagger H\) | Zero renormalized coefficients at \(\mu_0\) | All allowed portals; no gravity-loop protection asserted |
| \(|X|^2D^2\) | Zero coefficient at \(\mu_0\) | Allowed portal; scalar loops involving \(\phi\) can generate it |
| \(R\phi^2,R|X|^2,RD^2,RH^\dagger H\) | All \(\xi_R=0\) at \(\mu_0\) | Common minimal leading coupling at matching, not an RG fixed requirement |
| \(R^2,R_{\mu\nu}R^{\mu\nu}\), curvature basis equivalents/boundaries | Zero finite renormalized coefficients at \(\mu_0\) | Counterterms still retained; basis and topological terms handled consistently |
| Remaining symmetry-allowed higher-dimension operators | Zero matching coefficients in declared operator basis | A proposed EFT matching point, not a derived UV result |
| SM gauge, Higgs, Yukawa/CP parameters | Conventional independently calibrated values; strong theta set to zero by chosen matching | Imported; not needed numerically for the scalar vacuum calculation |
| Right-handed Majorana and Weinberg operators | Absent by imposed unbroken EFT \(B-L\) | Additional symmetry assumption; Dirac nature not established |
| Initial charge/profile | Fixed by selected finite numerical solution, initially aiming at \(F(0)=2.5\) | Preparation/state datum, not a coupling fit |
| Reporter frequency, amplitude, reference phase and readout | Explicit values in reporter execution | Apparatus/preparation data |
| Cutoff | Integrator's illustrative \(\Lambda_{\rm EFT}=0.1M_P\) | Scope assumption; no demonstrated new threshold there |

The polynomial scalar list is complete at dimension \(\le4\) given the stated discrete/global symmetries, up to vacuum/mass/quartic/kinetic renormalization already in the displayed action and the absent operators listed above. Terms odd in \(\phi\) or \(D\), and terms violating \(X\) number, are excluded by the declared symmetries. Universal metric coupling does not forbid the allowed portals. In particular, the \(g_X^2\phi^2|X|^2\) and \(h_D^2\phi^2D^2\) vertices generate an \(|X|^2D^2\) counterterm even with its matching coefficient zero.

All necessary local counterterms, running and nonlocal loop effects are retained at the perturbative order actually calculated. Zero finite Wilson coefficients in a specified scheme and basis are legitimate matching inputs, but they are not evidence of naturalness, a unique microscopic model, or decoupled physics at every scale. Higher derivatives are treated as order-by-order EFT insertions; one must not resum the truncation to claim physical high-energy ghost modes. A finite-order leading calculation needs only the coefficients contributing at that order, not a nonperturbative definition of the entire infinite tower.

Fields and gradients in the worked solution must remain within the declared EFT domain. \(F(0)=2.5\) corresponds to \(f(0)=2.5\times10^{-5}M_P\), well below the illustrative cutoff; characteristic inverse lengths are of order \(g_Xv\) multiplied by numerical profile factors. This is an initial scale check; actual profile gradients, curvature, the classical fluctuation spectrum and loop estimates still need verification. Unknown alternative matching coefficients are physical model sensitivity, not numerical discretization error.

## 8. Main-file audit and evidence boundaries

Reviewed `/workspace/scratch/9a57d63bee22/GRUT_Candidate_1/GRUT_CANDIDATE_1.md` during integration. The displayed scalar field equations, scalar masses, Hamiltonian/momentum constraints and scalar Hamiltonian factors were consistent. Sent the integrator the binding inequality, normalized equations, allowed mass/portal operators, Dirac-neutrino symmetry requirement, spinor-convention issue and state-space/infrared caveats. The integrator reports patching the remaining convention and scope points. This is not a claim that the final exported file or numerical results have passed an independent full audit.

No fatal internal obstruction to this **leading EFT proposal** was found. The following remain conditional until supplied by the execution workers:

1. A converged finite-energy profile with the frozen couplings and declared preparation, rather than a diverging or finite-box artifact.
2. Charge, energy and boundary normalization verified from that profile.
3. A fixed-charge fluctuation calculation with clearly stated angular sectors, symmetry zero modes, continuum threshold, discretization/tail control, and distinction between linear/local stability, fission resistance, tunneling and full nonlinear longevity. \(E/Q<m_X\) only forbids complete decay to free unit-charge X quanta; it does not by itself exclude energetically allowed fission into smaller bound objects.
4. A consistent weak-gravity calculation including the relevant stress and a measured compactness bound.
5. Reporter scattering on the same numerical background, with phase/probability conventions and weak-probe assumptions checked.
6. A specified comparison/testing protocol. Numerically displaying an EFT prediction is not a new experimental observation, independent calibration, or evidence favoring this sector over standard GR+SM and conventional scalar alternatives.

If the workers deliver only radial linear stability, the final demonstrated claim must remain radial linear stability. Translation and global-phase modes are symmetry zero modes, not instabilities. A nonlinear formation experiment, finite-temperature lifetime, complete quantum spectrum and phenomenological exclusion study require separate work. Do not transfer any previous graph/density/record recovery theorem to this unrelated scalar architecture.

## Primary literature actually consulted

The formulas above are the constructor's direct derivations for the displayed action; the sources below identify antecedents and framework limitations rather than validate this particular benchmark.

- R. Friedberg, T. D. Lee and A. Sirlin, *Class of scalar-field soliton solutions in three space dimensions* (1976), original scalar nontopological-soliton construction: https://doi.org/10.1103/PhysRevD.13.2739
- S. Coleman, *Q-balls* (1985), conserved-charge soliton construction: https://doi.org/10.1016/0550-3213(85)90286-X . Original article copy hosted by the research institute: https://www.ift.uam-csic.es/sites/default/files/QBALLS.pdf
- J. F. Donoghue, *General relativity as an effective field theory: The leading quantum corrections* (1994), low-energy quantum-gravity EFT: https://arxiv.org/abs/gr-qc/9405057
- R. Arnowitt, S. Deser and C. W. Misner, *The Dynamics of General Relativity*, original chapter republication: https://arxiv.org/abs/gr-qc/0405109
- W. Donnelly and S. B. Giddings, *Observables, gravitational dressing, and obstructions to locality and subsystems* (2016): https://arxiv.org/abs/1607.01025
- C. J. Fewster and R. Verch, *Quantum fields and local measurements* (2018/2020), local probe/instrument framework on a supplied curved spacetime: https://arxiv.org/abs/1810.06512
- J. Heeck, *Unbroken B−L Symmetry* (2014), primary model paper discussing unbroken B−L and Dirac neutrinos: https://arxiv.org/abs/1408.6845
- B. S. Merzlikin, I. L. Shapiro, A. Wipf and O. Zanusso, *Renormalization group flows and fixed points for a scalar field in curved space with nonminimal F(phi)R coupling* (2017), for the narrower point that nonminimal scalar curvature couplings run: https://arxiv.org/abs/1711.02224 . No numerical beta function from that paper is used here.
