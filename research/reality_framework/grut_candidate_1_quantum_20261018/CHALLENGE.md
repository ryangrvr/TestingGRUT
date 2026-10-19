# Adversarial audit: Candidate 1 bounded scalar quantum check

Audit date: 19 October 2026. The experiment folder retains its 20261018 identifier. This audit fixes the same action and benchmark matching; no new mechanism, fitted coupling, or finite counterterm is authorized by a favorable desired result. The owned file is this challenge only.

## Acceptance boundary

The calculation can establish a one-loop pure-new-scalar flat-space renormalization check, constant-background Coleman–Weinberg (CW) diagnostics, and the stated order-by-order renormalization-scale cancellation. It cannot establish a quantum Q-ball energy, fixed-charge binding threshold, particle pole mass, quantum stability/lifetime, or quantum gravitational consistency.

In particular, **small CW changes at the vacuum and a selected homogeneous bulk point do not validate the rotating, inhomogeneous soliton quantum mechanically.** The absence of the relevant finite-frequency fluctuation determinant is a substantive missing calculation, not a numerical precision issue. A negative local static Hessian on the wall is also not by itself a rejection of the soliton: the omitted gradients and charge constraints are essential to that stability problem.

## Normalization and four-mode counting

Write `X=(x+i y)/sqrt(2)` with four canonically normalized real fields `(phi,x,y,D)` and `r^2=x^2+y^2`. The explicit convention read from `rg_protocol.json` and `effective_potential.py` is

\[
V_0=C+\tfrac12m_\phi^2\phi^2+\tfrac12m_X^2r^2+\tfrac12m_D^2D^2
+\tfrac14\lambda_\phi\phi^4+\tfrac14\lambda_Xr^4+\tfrac14\lambda_DD^4
+\tfrac12c\phi^2r^2+\tfrac12d\phi^2D^2+\tfrac12\alpha r^2D^2.
\]

Thus alpha multiplies `|X|^2 D^2` in complex notation. Its normalization must not be confused with a convention containing an extra half. The signed mass coefficients at matching are `m_phi^2=-0.01`, `m_X^2=0`, `m_D^2=0.01`, in units of v squared. They are not the positive vacuum curvature masses.

The scalar quartic tensor has `lambda_iiii=6 lambda_i` and mixed `lambda_iijj=2 portal_ij`, including the O(2) relation `lambda_xxyy=2 lambda_X`. The transverse y fluctuation is compulsory; counting X as only one real scalar gives incorrect counterterms and CW multiplicities.

## Independently checked beta coefficients

Define `B_z=16 pi^2 d z/d ln(mu)`. With the convention above, contraction of the canonical real-scalar one-loop tensor beta gives

\[
B_{\lambda_\phi}=18\lambda_\phi^2+4c^2+2d^2,
\quad B_{\lambda_X}=20\lambda_X^2+2c^2+2\alpha^2,
\quad B_{\lambda_D}=18\lambda_D^2+2d^2+4\alpha^2,
\]

\[
B_c=c(6\lambda_\phi+8\lambda_X+8c)+2d\alpha,
\]

\[
B_d=d(6\lambda_\phi+6\lambda_D+8d)+4c\alpha,
\quad
B_\alpha=\alpha(8\lambda_X+6\lambda_D+8\alpha)+2cd.
\]

The signed mass and constant betas are

\[
B_{m_\phi^2}=6\lambda_\phi m_\phi^2+4c m_X^2+2d m_D^2,
\]

\[
B_{m_X^2}=2c m_\phi^2+8\lambda_Xm_X^2+2\alpha m_D^2,
\]

\[
B_{m_D^2}=2d m_\phi^2+4\alpha m_X^2+6\lambda_Dm_D^2,
\qquad
B_C=\tfrac12[(m_\phi^2)^2+2(m_X^2)^2+(m_D^2)^2].
\]

Pure quartic scalar anomalous dimensions vanish at this loop order. This does not remove finite derivative terms in an inhomogeneous effective action.

At the declared matching point, `B_alpha=2e-6`, `B_mX2=-2e-4`, and `B_C=1e-4`. Consequently the restricted tree potential with alpha and the independent X mass held identically zero is not closed under scalar loops. Their renormalized zeros are allowed once at the specified scale, with the admitted counterterms and running. A negative running signed X mass coefficient does not alone imply a vacuum tachyon: the actual vacuum Hessian also contains `c phi_vac^2` and other terms. Likewise the tiny negative alpha obtained below the matching scale is not alone a boundedness violation.

The decisive algebraic cross-check is the full polynomial identity

\[
\sum_z B_z\,\frac{\partial V_0}{\partial z}
=\tfrac12\operatorname{Tr}H^2,
\qquad H_{ab}=\partial_a\partial_bV_0.
\]

It includes the constant, quadratic, and quartic coefficients. Omitting C, the generated X mass, or the generated X–D portal must leave an explicit residual. A finite numerical agreement at a few selected field points is weaker than the proposed exact coefficient comparison. This is closure of the stated flat scalar potential family only. Curved-background scalar loops also require nonminimal and curvature counterterms even without graviton loops; Standard Model, gravitational, and derivative-sector checks are not included here.

## CW source audit

`effective_potential.py`, `potential_protocol.json`, and `potential_results.json` were read directly. The implementation uses all four real Hessian eigenvalues and

\[
V_1=\frac1{64\pi^2}\sum_a h_a^2[\log(h_a/\mu^2)-3/2],
\quad h_a=\text{eigenvalues of }H\text{, i.e. mass squared}.
\]

The notation `h_a^2` here is a squared mass-squared eigenvalue, not a second independent square of a physical mass. The factor, scalar constant, and multiplicities are consistent with the MSbar scalar formula. For positive eigenvalues the implemented Frechet derivatives of the trace use the divided differences of `g(h)=h(log(h/mu^2)-1)`; the degeneracy limit is `g'(h)=log(h/mu^2)`. Those matrix derivatives and the strict first-order vacuum-shift expansion are consistent. In particular, off-diagonal fluctuation mixing cannot generally be omitted by differentiating only diagonal entries of H.

One bookkeeping error was identified and sent to the root: the initial version's `bulk.tree_background_masses_squared` used `mb` after it had been overwritten by the Hessian at the truncated optimum. The correct eigenvalues at the tree bulk field `F=10^(1/4)` are approximately `[0.00316227766017,0.00948683298051,0.01,0.02162277660168]`. The initial JSON instead labeled the nearby optimum eigenvalues as tree values. The independently stored tree CW value and ratio-shift calculations were unaffected. The final source/results were reread and this label is now corrected; the initial output and a correction record are retained. No coupling or physical result was adjusted by the correction.

## Vacuum and constant matching

The fixed-C scalar-loop vacuum energy at the tree point is

\[
V_1(1,0,0,0)=-1.22652632365\times10^{-6}\,v^4.
\]

The strict field shift is `delta phi/v=1.21937398077e-4`. The first-order curvature-squared values are approximately `0.02001277843` for phi, `0.01000277531` for each X component, and `0.01009983063` for D, in units of v squared. These are local potential curvatures in the specified renormalized field convention. They are not momentum-dependent particle pole masses or complete physical screening masses. The latter require the relevant two-point/kinetic information.

The numerical minimum of `V0+V1` is separately near `phi/v=1.00012188180`, with vacuum value `-1.22667495222e-6 v^4`. Its difference from the strict first-order energy is approximately `-1.49e-10 v^4`. Keeping this shifted minimum is an admissible truncated-potential diagnostic, but its induced higher-order pieces are not a computed two-loop result.

At fixed original matching the nonzero raw vacuum value means that exact one-loop asymptotic flatness does not follow from the tree calibration. A constant change of approximately `+1.22652632365e-6 v^4` would restore zero at strict first order, but it is a new explicitly stated total-constant matching condition. It cannot silently be fitted or called a prediction. In the full gravity-plus-matter action only the appropriately matched total field-independent vacuum density gravitates; splitting it into scalar C and a separately written cosmological term does not create independent observable freedoms.

Subtracting the vacuum value in a relative bulk ratio cancels a constant for a matter energy comparison. It does not impose the gravitational cosmological-constant condition or solve the cosmological-constant problem. Algebraic Einstein curvature scales inferred from the raw density are conditional diagnostics, not a computation of graviton loops or the self-consistent quantum metric.

## Static bulk ratio is not quantum finite-charge binding

The reported relative static ratio is `2[Veff(phi,F)-Veff(vacuum)]/F^2`. At tree order its minimum is `0.00316227766017 v^2`; the strict one-loop correction is `+1.17872091832e-7 v^2`. The local truncated minimum remains near `phi=0,F=1.77846602v`, with ratio `0.00316239555657 v^2`. The sampled rectangle lies in a positive-tree-Hessian region, and the declared tree convexity bound has a positive determinant lower bound. This makes the selected real local potential calculation well defined there. The sample does not certify a global quantum minimum.

The square root of that ratio divided by a tree or curvature X mass is a named static ratio, not a computed quantum emission threshold. A pole mass has not been calculated, and the charged state is rotating. The relevant fluctuation problem in a rotating frame contains finite omega terms and mixing of the two X components. At the actual classical reference omega is about `0.08v`, while vacuum masses are about `0.1v`; no parametrically small time-derivative expansion has been established. Finite charge, wall gradients, changes in the Noether-charge functional, zero modes, derivative operators, and the nonlocal determinant are absent from ordinary constant-field CW.

Therefore neither a small positive ratio correction nor a ratio below one proves quantum binding, survival, or a lifetime of the finite-Q state. No original classical energy-versus-free-particle comparison is upgraded to a quantum theorem by these numbers.

## Imaginary parts and the profile wall

The implementation retains `log(h<0)=log|h|+i pi` and records the sign convention. This is preferable to clipping negative eigenvalues to zero or retaining only `log|h|` while declaring the answer real. A branch-conjugate/Feynman convention can reverse the displayed sign; the nonzero magnitude and identified unstable constant background are the relevant diagnostic here.

The scan finds negative static Hessian eigenvalues at 84 of 1001 sampled classical-profile points, over dimensionless radii approximately 2.64 through 4.30; the minimum is about `-0.00357857 v^2`. These are constant-background Hessians evaluated on the wall values, not eigenvalues of the full charged soliton fluctuation operator. They are neither a proven decay channel nor evidence that the earlier finite-time classical tests were inconsistent. The exact ground-state Legendre effective potential and a complex perturbative homogeneous-background potential must not be conflated.

The protocol explicitly declines to integrate a local CW density over this rotating profile. That restraint is necessary: both the negative-mode region and the omitted finite-frequency/gradient effects prevent treating such an integral as the renormalized one-loop Q-ball energy. Similarly, the perturbative homogeneous-background imaginary part is not a computed Q-ball lifetime or tunneling rate.

## RG scale test: required interpretation

At fixed background and one scalar loop, all parameters are run from the one specified matching point. The running must include C, all signed mass coefficients, and generated alpha; a zero can be imposed at mu0, not reset at each evaluation scale. The identity above cancels the explicit `mu dV1/dmu=-Tr(H^2)/(32 pi^2)` through first loop order. The scalar background has no one-loop anomalous rescaling in this restricted theory.

The remaining variation of `V0(running)+V1(running)` is generally of the next loop order, not exactly zero. Integrating one-loop beta functions resums some logarithmic terms; it does not deliver the omitted genuine two-loop potential. A consistently applied loop-counting multiplier on both beta and V1 can test quadratic scaling of the residual. Modifying physical couplings or refitting C at each scale would test a different problem.

Small factor-two scale dependence is internal consistency evidence, not a universal uncertainty bound, all-orders renormalization, UV validity, or empirical evidence. No scale cancellation at a constant positive-Hessian point resolves the wall or finite-omega determinant problem.

## Primary sources inspected

- S. P. Martin, *Two-loop effective potential for a general renormalizable theory and softly broken supersymmetry*, Phys. Rev. D 65, 116003 (2002), [hep-ph/0111209v2](https://arxiv.org/html/hep-ph/0111209v2). The real-field counting, scalar one-loop convention, and RG organization provide the relevant normalization reference; no supersymmetric result is imported as this model's spectrum.
- M. Luo, H. Wang, Y. Xiao, *Two-loop Renormalization Group Equations in General Gauge Field Theories*, Phys. Rev. D 67, 065019 (2003), [hep-ph/0211440](https://arxiv.org/pdf/hep-ph/0211440). General mass/quartic renormalization conventions were checked against the explicit one-loop contractions used here.
- E. J. Weinberg, A. Wu, *Understanding complex perturbative effective potentials*, Phys. Rev. D 36, 2474 (1987), [DOI 10.1103/PhysRevD.36.2474](https://doi.org/10.1103/PhysRevD.36.2474). The primary abstract concerns the instability of a specified homogeneous state; it does not identify a Q-ball decay rate from a local wall Hessian.

## Final audit status

The final `rg_symbolic.py`, `rg_results.json`, `flow_check.py`, `flow_results.json`, corrected potential output, and `POTENTIAL_NOTES.md` were inspected. The independent tensor contraction and directly differentiated polynomial routes agree exactly: all 50 monomials on each side of the RG identity match, with zero residual. The reported 290 checks comprise deterministic arithmetic, tensor-component, symmetry, and identity checks within this derivation; they are not 290 independent physical experiments. The displayed beta coefficients match the independent contractions above. Omitting alpha, the independent X mass, or C leaves the recorded nonzero polynomial residuals.

In particular, the running tree minimum value is not identically zero. At matching,

\[
\frac{d}{d\ln\mu}\left[C-\frac{(m_\phi^2)^2}{4\lambda_\phi}\right]
=2.22276512902\times10^{-6}v^4
\]

at the stated one-loop order. This equals the appropriate vacuum Hessian trace contribution. Keeping `C=(m_phi^2)^2/(4 lambda_phi)` at every scale would contradict this flow rather than preserve the original matching.

The numerical flow uses all ten coefficients from the same fixed matching point at `mu/mu0=0.25,0.5,1,2,4`, at seven declared four-field positive-Hessian backgrounds. It does not reset the vacuum constant or change the field units with mu. The independently entered direct mass matrix agrees with the tensor conventions. The maximum pointwise RG identity relative residual is `1.04e-15`; the strict first-order scale spread is at most `8.67e-19 v^4`. The genuinely integrated one-loop running plus CW has residual scale spreads approximately `5.14e-9` through `7.96e-9 v^4`. Changing ODE tolerances changes parameters by at most `1.73e-18`, far below those residuals. These values support the intended algebraic and numerical scale-cancellation checks.

The eta bookkeeping test specifically uses the leading-log trajectory `theta0+eta log(mu/mu0) beta(theta0)` together with `eta V1`; it is not a separately integrated resummed trajectory for each eta. Its residual divided by eta squared changes by at most `8.99e-4` between adjacent tested eta values. This is the expected next-order scaling diagnostic and does not supply a computed two-loop coefficient or rigorous theory-error bar.

**Final verdict: the expanded flat scalar potential survives the one-loop closure and RG-consistency challenge. The finite-Q quantum-soliton claim remains untested.** The small local vacuum/bulk corrections, retained imaginary parts, fixed matching, and explicit scope are internally consistent. The corrected output does not contain a remaining normalization blocker.

| Proposed conclusion | Audit verdict |
|---|---|
| The displayed restricted tree potential can keep alpha and mX2 identically zero under scalar running | Rejected; both are generated at one loop |
| The admitted full ten-parameter flat scalar potential family closes at one loop | Supported by the tensor/polynomial calculation |
| The stated constant-field one-loop expression has the required leading RG cancellation | Supported at the specified order and backgrounds |
| The local vacuum curvatures remain positive at this matching point | Supported as a local potential-curvature diagnostic |
| Tree zero vacuum density remains automatically zero quantum mechanically with fixed C | Rejected; the explicit scalar loop gives a nonzero raw vacuum shift |
| The static bulk ratio supplies quantum binding, a pole threshold, or a Q-ball lifetime | Rejected; the necessary finite-frequency, momentum, and finite-charge calculations are absent |
| Local negative wall Hessians prove the finite soliton decays | Rejected; these are not the full constrained fluctuation operator |
| Full curved-space/gravitational/Standard Model EFT closure or ultraviolet completion has been demonstrated | Not established by this flat pure-new-scalar check |

The numerical results are useful execution of a previously admitted EFT obligation. They are not a quantum upgrade of the earlier classical persistence evidence, external validation, or a new empirical success.
