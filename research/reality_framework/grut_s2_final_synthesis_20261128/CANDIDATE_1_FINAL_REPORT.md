# Candidate 1 final scientific report

**C1/C2/C3/C4 verdicts A/A/C/C. Unchanged Candidate 1 investigation closed.** Candidate 1 is a mathematically supported example of a localized charged state in a specified conventional scalar matter EFT. The results apply to prepared configurations at the original matching, as qualified below. Formation and an accessible distinctive signal are unsupported. The physical existence of this matter in nature remains unresolved. [Program verdict](FINAL_PROGRAM_VERDICT.md) · [Pins and replay](REPRODUCTION_INDEX.md)

## 1. Original covariant EFT and classical binding

For signature $(-+++)$ and $\hbar=c=1$, the leading action is
$$
S_0=\int d^4x\sqrt{-g}\left[\frac{M_P^2}{2}R-\rho_\Lambda+\mathcal L_{\mathrm{SM}+\nu_R}
-\tfrac12(\nabla\phi)^2-\nabla_\mu X^*\nabla^\mu X-\tfrac12(\nabla D)^2-U\right]+S_{\rm boundary},
$$
$$
U=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
+g_X^2\phi^2|X|^2+\lambda_X|X|^4
+\frac12(m_{D0}^2+h_D^2\phi^2)D^2+\frac{\lambda_D}{4}D^4.
$$
Einstein gravity and ordinary quantum matter with three right-handed neutrinos enter as inputs. The complex $X$ carries $U(1)_X$, and the binding field $\phi$ and reporter $D$ are real, with discrete parities in the declared EFT. The vacuum choice is $X=D=0,\phi=+v$ and renormalized $\rho_\Lambda=0$ for the asymptotically flat benchmark. The reduced Planck scale is $M_P=(8\pi G)^{-1/2}$. The leading masses are $m_X=g_Xv$, $m_\phi=\sqrt{2\lambda_\phi}v$ and $m_D=(m_{D0}^2+h_D^2v^2)^{1/2}$. [Original action and matter definitions](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md)

The matching fixes $g_X=0.1$, $\lambda_\phi=0.01$, $\lambda_X=0.001$, $\lambda_D=0.001$, $h_D=0.01$, $m_{D0}=g_Xv$, $v/M_P=10^{-5}$ and $\mu_{\overline{\rm MS}}=g_Xv$. Absent renormalized portals, curvature couplings, bare carrier mass parameter and higher-order coefficients are zero at that matching under the original prescription. Symmetry-allowed counterterms and running are retained at the working order. C1 includes only the leading flat saddle and new-scalar loops. C2 uses the leading flat scalar equations. C3 and C4 use their own probe and gravity orders. Formation-sector extensions have separate actions and evidence. [Matching and running](../grut_candidate_1_quantum_20261018/RG_DERIVATION.md) · [C1 protocol](../grut_c1_quantum_channels_20261127/PROTOCOL.json)

At $D=\phi=0$ and $X=f/\sqrt 2$, minimizing $2U/f^2$ over positive $f^2$ gives $v^2\sqrt{\lambda_\phi\lambda_X}<m_X^2$ for $g_X^4>\lambda_\phi\lambda_X$. This is the conventional bulk binding criterion from lowering the carrier mass inside. The finite charged object must also pay gradient and quartic energy. Set $\rho=g_Xvr$, $\tau=g_Xvt$, $X=vF(\rho)e^{-i\Omega\tau}/\sqrt2$, $\phi=vS(\rho)$, $D=0$. The Noether current $j_X^\mu=i(X\nabla^\mu X^*-X^*\nabla^\mu X)$ gives
$$
I=\int_0^\infty \rho^2F^2\,d\rho,\qquad \bar Q=\frac{g_X^2 Q}{4\pi}=\Omega I.
$$
Regular radial data have $F'(0)=S'(0)=0$. At infinity $F\to0$ and $S\to1$. At fixed $\Omega$, $K_0=E_0-\omega Q$ with $\omega=g_Xv\Omega$. Stationarity gives
$$
-\frac{(\rho^2F')'}{\rho^2}+(S^2+0.1F^2-\Omega^2)F=0,\qquad
-\frac{(\rho^2S')'}{\rho^2}+(S^2-1+F^2)S=0 .
$$
The leading energy is $(4\pi v/g_X)\int_0^\infty \rho^2(K_t+K_r+W)d\rho$ with $K_t=\Omega^2F^2/2$, $K_r=(F'^2+S'^2)/2$, $W=(S^2-1)^2/4+S^2F^2/2+F^4/40$. Prepared binding is assessed by actual energy differences. The earlier fixed-charge full breakup comparison gave $E(39944)-39944m_{X,\mathrm{pole}}=-3553.70573404\,g_Xv$ with a roughly $0.11\,g_Xv$ registered allowance. C4's C1-based same-order central mass defect is about $-3553.7132\,g_Xv$. The two central values are from different registered numerical dossiers at the same matching. [Original worked branch](../grut_candidate_1_20261018/WORKED_REGIME.md) · [Earlier full breakup report](../grut_candidate_1_energy_20261019/REPORT.md) · [C4 mass](../grut_c4_gravitational_accessibility_20261127/REPORT.md)

## 2. C1 fixed-charge quantum comparisons

In angular momentum $\ell$, use measure $d\rho$ by multiplying each field perturbation by $\rho$. The rotating $X$ amplitude and $\phi$ perturbations have stiffness
$$
L_{+,\ell}=-\partial_\rho^2+\frac{\ell(\ell+1)}{\rho^2}+
\begin{pmatrix}
 S^2+0.3F^2-\Omega^2 & 2FS\\
 2FS & 3S^2-1+F^2
\end{pmatrix}.
$$
The phase form is $L_{-,\ell}=-\partial_\rho^2+\ell(\ell+1)/\rho^2+S^2+0.1F^2-\Omega^2$. The amplitude–phase velocity coupling has skew entries $\pm2\Omega$. The reporter is decoupled at quadratic order and has $H_{D,\ell}=-\partial_\rho^2+\ell(\ell+1)/\rho^2+1+0.01S^2\ge I_{\rm op}$. Its vacuum block has dimensionless mass squared $1.01$. [Rotating and reporter operators](../grut_candidate_1_energy_20261019/ENERGY_DERIVATION.md)

With $u=\rho F$ and normalized $e=(u,0)/\sqrt I$, the quadratic cyclic phase coordinate has canonical momentum $\dot\theta-2\Omega\langle e,\eta\rangle$ for amplitude perturbation $\eta$ in the archived convention. Its removal at zero first-order charge gives $L_Q=L_++4\Omega^2|e\rangle\langle e|$ in $\ell=0$. Its Euclidean Schur complement gives $\det P_{\rm full}(\zeta)=\zeta^2\det P_{Q_1=0}(\zeta)$ in an orthonormal kinetic basis. The full background determinant divided by the complete original vacuum block retains that factor. The phase and three translations are collective coordinates with fixed charge and rest momentum. Their four low-frequency logarithms survive the relative determinant, including when the singular part is integrated analytically. [Routhian and collective normalization](../grut_candidate_1_energy_20261019/ENERGY_DERIVATION.md)

With the reporter included, the dimensionless ratio $K_1/(g_Xv)$ is the integral over $\zeta>0$ and sum over $\ell$ of
$$
\frac{2\ell+1}{2\pi}
 \left(\log\frac{\det P_\ell}{\det P_{\ell,0}}-B_{1\ell}-B_{2\ell}\right)
$$
plus the two restored $\overline{\rm MS}$ graphs in $g_Xv$ units. The first and second Born terms are $\operatorname{tr}(G_0V)$ and $-\operatorname{tr}(G_0VG_0V)/2$ for the same vacuum-relative operator. The restored graphs keep the laboratory $\Omega$ and $2\Omega$ harmonics, and the reporter enters every applicable part at the fixed matching. For $K=K_0+\hbar K_1+\cdots$, $E(Q)=K(\omega)+\omega Q$ at the stationary $\omega=\omega_0+\hbar\omega_1+\cdots$. At first order $\omega_1$ multiplies $K_0'(\omega_0)+Q=0$. Thus the energy correction is $K_1(\omega_0)$. The pole uses the strict first-order expansion with reporter and vacuum shift. [Renormalization and pole](../grut_candidate_1_energy_20261019/ENERGY_DERIVATION.md) · [Pole thresholds](../grut_candidate_1_fluctuations_20261019/POLE_THRESHOLDS.md)

| $Q$ | $E_0/(g_Xv)$ | $K_1/(g_Xv)$ at retained cutoffs | Allowance on summed energy |
| ---: | ---: | ---: | ---: |
| 39,944 | 36392.42286136 | 2.78783262 | 0.02668966 |
| 39,943 | 36391.62286052 | 2.78777941 | 0.02668966 |
| 19,972 | 19729.91009896 | 1.65541564 | 0.02668972 |

The strict expanded $m_{X,\mathrm{pole}}/(g_Xv)$ is $1.0001232694040532$. The archived sign convention is $\Delta_1=E(39944)-E(39943)-m_{X,\mathrm{pole}}$ and $\Delta_{1/2}=E(39944)-2E(19972)$. A negative interval gives energetic protection in the channel at the stated order and allowances. A positive interval gives energetically accessible final states. An interval containing zero is inconclusive.

| Outcome | Energy difference in $g_Xv$ units |
| --- | ---: |
| Paired one-carrier threshold | $-0.2000692023\pm0.0032516234$ |
| Independent numerical integrated chemical-potential interval | $[-0.201140826,-0.198997579]$ |
| One-carrier threshold adding separate absolute allowances | $[-0.253448528,-0.146689876]$ |
| Equal fission into the tested daughter branch | $-3067.92033523\pm0.08006910$ |

The threshold comparison uses correlated neighboring-charge subtraction on matched grids. The direct and chemical-potential centers differ by $3.71\times10^{-11}$ while sharing the determinant and calibration. The paired allowance combines an analytical padded-spline tail estimate $0.0000516234$, empirical allocations totaling $0.0032$ and $10^{-8}$ for the pole. The separate check adds both absolute-energy allowances and the pole allowance. The equal split spends the half-charge error twice. Discrete budgets use numerical refinement and analytical inequalities evaluated on numerical profiles. The C2 certificate supplies the leading exact scalar background premise, leaving a rigorous all-cutoffs continuum determinant bound and higher-loop/gravitational corrections outside C1. The reported curvature diagnostic has no separate verifier gate. [C1 report and full budget](../grut_c1_quantum_channels_20261127/REPORT.md) · [C1 critique](../grut_c1_quantum_channels_20261127/INTERNAL_CRITICAL_REVIEW.md)

## 3. C2 certified continuum linear stability

Use $(u,w)=(\rho F,\rho(S-1))$ on the half-line, regular at the origin and decaying to vacuum. The archive bounds a bordered inverse for the stationary equations and normalized charge. Bernstein cell enclosures bound the residual and curve derivatives. An exact inverse of the predictor derivative defines a Newton fixed-point map on the graph norm with components $\|\eta\|_2$, $\|\eta''\|_2/16$, $|\delta\Omega|$. For inverse norm bound $\gamma$ and residual bound $\epsilon$, the accepted radius $r=2\gamma\epsilon$ and $360\gamma r<1/2$ establish a contraction in the recorded tube. The exterior ODE gives a decaying bound $\|(u,w)(\rho)\|_{\mathbb R^2}\le3r\,e^{-\kappa(\rho-R)}$ beyond the taper, with controlled derivative. This trailing predictor taper is a certification tool. The exact solution retains its tails. [C2 method](../grut_c2_continuum_stability_20261127/CERTIFICATE_METHOD.md)

The continuum spectral lower form uses Neumann separation and a positive vacuum exterior for the predictor, a cellwise constant lower form and a Poincaré estimate for the endpoint-zero remainder. A verified negative shifted trial below $-g$ and lower-form index at most one below $+g$, transferred through the exact profile tube, certify a unique negative radial direction and positive gap. The original propagated interval LDL failed at block 123. The frozen amendment instead encloses point-factor residuals and inverse norms at the same shifts, mesh and precision, verifying that the homotopy to the actual interval form stays invertible. This transfers the point inertia. [Failed result](../grut_c2_continuum_stability_20261127/INITIAL_ATTEMPT.json) · [Amendment](../grut_c2_continuum_stability_20261127/AMENDMENT_01.md) · [Critical review](../grut_c2_continuum_stability_20261127/INTERNAL_CRITICAL_REVIEW.md)

Differentiating the stationary system at fixed frequency gives $L_+\partial_\Omega(u,w)=(2\Omega u,0)$. Since $\bar Q=\Omega I$, one obtains $1+4\Omega^2\langle e,L_+^{-1}e\rangle=\bar Q'/I$. Consider the block form with diagonal $(L_+,-1)$ and off-diagonal $2\Omega e$. Eliminating $-1$ leaves $L_Q$ and one negative scalar. Eliminating invertible $L_+$ leaves one negative direction in $L_+$ and the scalar $-\bar Q'/I>0$. The two reductions imply index zero and nullity zero for $L_Q$. Its positive essential floor persists under the rank-one perturbation. This reasoning uses the verified gap and nonzero negative slope. [Rank-one derivation and archived enclosures](../grut_c2_continuum_stability_20261127/REPORT.md)

| $Q$ | Enclosed exact $\Omega$ | Profile radius upper | $g$ lower | $d\bar Q/d\Omega$ |
| ---: | --- | ---: | ---: | --- |
| 39,944 | [0.79999952449, 0.79999972818] | $1.019\times10^{-7}$ | 0.03121578 | [−326.52232, −324.92872] |
| 39,943 | [0.80000196660, 0.80000217233] | $1.029\times10^{-7}$ | 0.03121543 | [−326.51768, −324.90815] |
| 19,972 | [0.88507412869, 0.88507426173] | $6.652\times10^{-8}$ | 0.03122765 | [−99.23759, −98.77926] |

For each configuration the radial phase kernel is $\rho F$, and the amplitude $\ell=1$ kernels arise from $\rho(F',S')$ times the three translation harmonics. The positive complement gaps are at least $0.01559078$, $0.01559043$ and $0.01560265$ respectively in both symmetry channels. Higher sectors use centrifugal ordering plus a positive exterior floor. The reporter positivity bound is one at every $\ell$. On fixed charge and zero rest momentum after quotienting the symmetry coordinates, the skew gyro preserves positive quadratic energy. This is continuum linear energy stability of the scalar system, with control at each of the three charges. [C2 sector verdicts](../grut_c2_continuum_stability_20261127/REPORT.md)

## 4. C3 reporter scattering and empirical feasibility

On the flat static $D=0$ background the regular reporter wave obeys
$$
u_\ell''+[k^2-\ell(\ell+1)/\rho^2-0.01(S^2-1)]u_\ell=0,\qquad k=p/(g_Xv).
$$
Equal incident and elastic outgoing Klein–Gordon flux factors $p/E_D$ give $\sigma_\ell=4\pi(2\ell+1)\sin^2\delta_\ell/p^2$. The reference S-wave at $k=5$ was $0.004066511575$ rad. At exact parent charge, the predictor gives $0.004066517304$ rad with a separate $1.52\times10^{-8}$ rad exact classical profile allowance. At $k=3,7$ it gives $0.006776520524,0.002904773628$ rad. At $k=5$, $\sigma_0=5.46\times10^{-62}\,\mathrm m^2$ and the retained angular sum is $1.15\times10^{-59}\,\mathrm m^2$. Empirical angular cutoff refinements support that retained sum. Recoil and vacuum pole diagnostics have distinct scopes and the full inhomogeneous quantum elastic and inelastic amplitudes remain open. [C3 phase, angular data and error sources](../grut_c3_empirical_feasibility_20261127/REPORT.md)

Here $g_Xv\simeq2.4353\times10^{12}\,\mathrm{GeV}$. The target rest energy is about $1.42\times10^7\,\mathrm J$ with 39,944 units of hidden net charge. Assembly from the corresponding cold free carrier population would have to dispose of about $1.39\times10^6\,\mathrm J$. Starting from ordinary uncharged matter requires compensating hidden charge elsewhere. A $k=5$ reporter quantum carries about $1.99\times10^3\,\mathrm J$. The hidden reporter interaction provides a test field equation. A calibrated source, prepared population and ordinary readout are absent. $D$ parity and the stated portal matching constrain proposed connections while gravity remains available. In the ideal dilute calculation $N_{\rm sig}=\mathcal F N_T\sigma_{\rm acc}\epsilon$, so with one target and perfect collection the retained total would require fluence about $8.66\times10^{58}\,\mathrm m^{-2}$ for mean one event. Background and exposure for an actual experiment are unassigned. [C3 source, target and detector analysis](../grut_c3_empirical_feasibility_20261127/REPORT.md)

A matched scalar well with the same vacuum-relative reporter mass profile yields all the same leading elastic phases. The angular dependence can distinguish a restricted isotropic comparator under calibrated preparation, and is equally a prediction of the matched conventional extended well. A correctly matched rejection experiment would test the theory and uncertainty in a calibrated regime. Its physical access requirements were unmet, giving C3 Verdict C. [C3 comparator and novelty](../grut_c3_empirical_feasibility_20261127/NOVELTY_AUDIT.md)

## 5. C4 gravity and the scope of closure

Einstein's equation with shared matter stress already gives a gravitational connection. In orthonormal source components in $g_X^2v^4$ units, $e=K_t+K_r+W$, $p_r=K_t+K_r-W$, $p_t=K_t-K_r-W$. Conservation gives $p_r'+2(p_r-p_t)/\rho=0$. The boundary term $[\rho^3p_r]_0^\infty$ vanishes, giving $\int\rho^2(p_r+2p_t)d\rho=0$ and the integrated active monopole from $e$. With $\alpha=v^2/(2M_P^2)$ and $\mathcal M(\rho)=\int_0^\rho s^2e(s)ds$, the areal lapse at leading order obeys $\nu'=\alpha[\mathcal M(\rho)/\rho^2+\rho p_r]$, $\nu(\infty)=0$. [C4 analytical derivation](../grut_c4_gravitational_accessibility_20261127/ANALYTICAL_GATE.md)

| Leading benchmark and bounds | Value |
| --- | ---: |
| Classical mass | $1.57993\times10^{-10}\,\mathrm{kg}$, about 158 ng |
| 90% energy radius | $4.13516\times10^{-28}\,\mathrm m$ |
| 90% charge radius | $3.82104\times10^{-28}\,\mathrm m$ |
| Compactness $GM/(R_{E90}c^2)$ | $2.83732\times10^{-10}$ |
| Central lapse perturbation | $-6.68814243\times10^{-10}$ |
| Exterior exact-profile mass fraction at $\rho=60.125$ | $\le3.636\times10^{-13}$ |
| Leading exterior scalar-tail lapse there | $\le8.738\times10^{-25}$ |
| Static material tidal lapse for the registered tide/readout | $\le2.546\times10^{-133}$ |

The classical mass and small metric perturbations support the leading weak-field treatment. The scalar one-loop rest mass estimate is $1.58005\times10^{-10}\,\mathrm{kg}$. The full gravitating correction and renormalized local stress at that order remain open. A spherical vacuum comparator at the same gravitating mass and external time normalization has the Schwarzschild exterior. The actual scalar tails are nonvacuum and are enclosed at leading gravitational order. A static, asymptotic monopole measurement alone detects mass and cannot identify composition. The stationary rotating phase gives stationary spherical stress. The other dynamical options were screened, with no physical mode-excitation, scattering or radiation pathway established. [C4 report and controls](../grut_c4_gravitational_accessibility_20261127/REPORT.md)

The sole bounded dynamical calculation concerns material rearrangement under a static prescribed weak quadrupole. The active source $J=e+p_r+2p_t=4K_t-2W$ gives the forcing $b_2=\rho^3\chi\nabla_{F,S}J$. Its response $\langle b_2,L_{+,2}^{-1}b_2\rangle$ is bounded by $\|b_2\|^2/\lambda_{*,2}$ with $\lambda_{*,2}=1/(144R^2)$. For $10^8\,\mathrm s^{-2}$ and a $1\,\mu\mathrm m$ readout the listed lapse is at most $2.170\times10^{-102}$ of the monopole. Coupled metric, contact and matching terms of a full relativistic tidal response were outside that bound. A conventional scalar soliton with matched action, state and Einstein coupling gives identical predictions. Gravity can see a hypothetical target's mass, but this investigation found neither a prepared target nor a distinctive measurable structural effect. [C4 tidal scope](../grut_c4_gravitational_accessibility_20261127/ANALYTICAL_GATE.md) · [C4 negative results](../grut_c4_gravitational_accessibility_20261127/NEGATIVE_RESULTS.md)

## 6. Final scope by scientific claim

| Claim | Final status |
| --- | --- |
| Specified asymptotic decay-channel energetics | Negative C1 intervals at the calculated scalar one-loop order and allowances for one-carrier emission and equal fission, plus the earlier complete breakup comparison. Other partitions and higher orders remain open. |
| Local fixed-charge linear stability | Certified at all three charges in the leading flat continuum scalar theory with identified symmetry modes removed. |
| Nonlinear orbital stability, absolute stability and lifetime | Unresolved. A stable linearized energy and favorable specified channels do not establish these properties. |
| Physical formation and abundance | Unresolved. The original action has no demonstrated prepared target or physical population. Conditional histories in other matter variants do not supply one. |
| Experimental accessibility | Unsupported for C3 and C4. Matched conventional comparators reproduce the tested responses. The gravitational coupling is physical within the adopted theory. |
| Scientific disposition | The authorized route has ended. The benchmark EFT remains in the archive with no empirical identification and no physical exclusion claimed. |

The C-series preserves the mathematically valid results and retains the failed tests. Foundational Verdict B/C governs the program. External replication is pending. [Evidence map](FINAL_EVIDENCE_MAP.md) · [Open questions](OPEN_PHYSICAL_QUESTIONS.md) · [RAI review packet](EXTERNAL_RAI_REVIEW_PACKET.md)
