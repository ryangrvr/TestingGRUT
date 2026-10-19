# Candidate 1: one-loop pure-scalar RG closure

Protocol frozen 18 October 2026; derivation/exact calculation completed during the continuing 18–19 October session. This is a bounded one-loop \(\overline{\rm MS}\) check of the existing scalar action and matching prescription. Physical couplings are unchanged. No Standard Model, graviton, curvature, thermal or finite-density loop is included. No quantum-corrected finite-charge soliton energy, pole mass, decay rate or lifetime is calculated.

## 1. Potential and conventions

Use four canonically normalized real fields \(z_a=(\phi,x,y,D)\), with \(X=(x+iy)/\sqrt2\) and \(R^2=x^2+y^2\). The renormalizable scalar potential is

\[
V_0=C+\frac12m_\phi^2\phi^2+\frac12m_X^2R^2+\frac12m_D^2D^2
+\frac{\lambda_\phi}{4}\phi^4+\frac{\lambda_X}{4}R^4+\frac{\lambda_D}{4}D^4
+\frac c2\phi^2R^2+\frac d2\phi^2D^2+\frac\alpha2R^2D^2.
\]

Here \(m_\phi^2,m_X^2,m_D^2\) are signed polynomial **mass-squared coefficients**, not physical pole masses. In particular, \(m_X^2\) denotes the bare/matching quadratic coefficient previously set to zero, while the tree vacuum excitation has squared mass \(m_X^2+c v_{\rm cl}^2\). The scalar symmetries are \(O(2)\) in \((x,y)\), \(\phi\mapsto-\phi\), and \(D\mapsto-D\). These admit \(\alpha\) even though its matching value is zero.

Define \(\beta_p=dp/d\ln\mu\) and write \(B_p=16\pi^2\beta_p\). The code's `mp,mx,md,lp,lx,ld` mean \(m_\phi^2,m_X^2,m_D^2,\lambda_\phi,\lambda_X,\lambda_D\). At \(\mu_0=0.1v\),

\[
(\lambda_\phi,\lambda_X,\lambda_D,c,d,\alpha)
=(0.01,0.001,0.001,0.01,0.0001,0),
\]

\[
(m_\phi^2,m_X^2,m_D^2)=(-0.01,0,0.01)v^2,\qquad C=0.0025v^4.
\]

Thus \(c=g_X^2\), \(d=h_D^2\), and the original tree double-well is reproduced. All zeros are conditions at \(\mu_0\), rather than instructions to set a permitted operator to zero at every scale.

## 2. Tensor derivation and scalar anomalous dimension

Write \(V_4=\Lambda_{abcd}z_az_bz_cz_d/4!\). The nonzero independent symmetric tensor entries are

\[
\Lambda_{0000}=6\lambda_\phi,\quad
\Lambda_{1111}=\Lambda_{2222}=6\lambda_X,\quad
\Lambda_{3333}=6\lambda_D,
\]

\[
\Lambda_{0011}=\Lambda_{0022}=2c,\quad
\Lambda_{0033}=2d,\quad \Lambda_{1122}=2\lambda_X,\quad
\Lambda_{1133}=\Lambda_{2233}=2\alpha,
\]

with every permutation of each entry included. At one loop, the scalar quartic diagrams give

\[
B_{\Lambda_{abcd}}=\sum_{e,f}\left(
\Lambda_{abef}\Lambda_{efcd}+\Lambda_{acef}\Lambda_{efbd}
+\Lambda_{adef}\Lambda_{efbc}\right).
\]

For a general symmetric quadratic matrix \((m^2)_{ab}\),

\[
B_{(m^2)_{ab}}=\sum_{c,d}\Lambda_{abcd}(m^2)_{cd},
\qquad B_C=\frac12\operatorname{Tr}[(m^2)^2].
\]

In this pure quartic scalar theory the one-loop two-point ultraviolet divergence is momentum independent. Consequently the one-loop scalar anomalous-dimension matrix vanishes, \(\gamma_{ab}^{(1)}=0\). This does not state that every finite momentum-dependent self-energy or effective-action derivative correction vanishes. Expanding about a broken vacuum produces cubic vertices but does not change this one-loop ultraviolet wavefunction conclusion.

For example, \(B_{\Lambda_{0000}}=3[(6\lambda_\phi)^2+(2c)^2+(2c)^2+(2d)^2]\), which must be divided by six to obtain \(B_{\lambda_\phi}\). This explicitly accounts for both real components of \(X\).

## 3. Complete one-loop beta functions at this scope

The six quartic functions are

\[
\boxed{B_{\lambda_\phi}=18\lambda_\phi^2+4c^2+2d^2,}
\]

\[
\boxed{B_{\lambda_X}=20\lambda_X^2+2c^2+2\alpha^2,\qquad
B_{\lambda_D}=18\lambda_D^2+2d^2+4\alpha^2,}
\]

\[
\boxed{B_c=c(6\lambda_\phi+8\lambda_X+8c)+2d\alpha,}
\]

\[
\boxed{B_d=d(6\lambda_\phi+6\lambda_D+8d)+4c\alpha,}
\]

\[
\boxed{B_\alpha=2cd+\alpha(8\lambda_X+6\lambda_D+8\alpha).}
\]

The three mass-squared functions and constant are

\[
\boxed{B_{m_\phi^2}=6\lambda_\phi m_\phi^2+4c m_X^2+2d m_D^2,}
\]

\[
\boxed{B_{m_X^2}=2c m_\phi^2+8\lambda_Xm_X^2+2\alpha m_D^2,}
\]

\[
\boxed{B_{m_D^2}=2d m_\phi^2+4\alpha m_X^2+6\lambda_Dm_D^2,}
\]

\[
\boxed{B_C=\tfrac12[(m_\phi^2)^2+2(m_X^2)^2+(m_D^2)^2].}
\]

The mass tensor remains diagonal with equal \(x,y\) entries. No additional symmetry-breaking quartic or off-diagonal mass is generated. Including \(\alpha,m_X^2,C\) and their counterterms closes the renormalizable **pure-scalar** polynomial sector at this order. It does not close the entire gravity/Standard Model EFT, and it does not establish higher-loop or ultraviolet completion.

## 4. Independent polynomial RG identity

Let \(H_{ab}(z)=\partial_a\partial_bV_0\) be the field-dependent squared-mass matrix. Its diagonal entries are

\[
H_{\phi\phi}=m_\phi^2+3\lambda_\phi\phi^2+cR^2+dD^2,
\]

\[
H_{xx}=m_X^2+c\phi^2+\lambda_X(3x^2+y^2)+\alpha D^2,
\quad
H_{yy}=m_X^2+c\phi^2+\lambda_X(x^2+3y^2)+\alpha D^2,
\]

\[
H_{DD}=m_D^2+d\phi^2+\alpha R^2+3\lambda_DD^2.
\]

The off-diagonal entries are \(H_{\phi x}=2c\phi x\), \(H_{\phi y}=2c\phi y\), \(H_{\phi D}=2d\phi D\), \(H_{xy}=2\lambda_Xxy\), \(H_{xD}=2\alpha xD\), and \(H_{yD}=2\alpha yD\).

For positive eigenvalues \(h_i\) of \(H\), the scalar one-loop potential is

\[
V_1=\frac1{64\pi^2}\sum_{i=1}^4h_i^2\left[\ln(h_i/\mu^2)-\frac32\right].
\]

Its explicit scale derivative at fixed parameters/fields is \(\partial_{\ln\mu}V_1=-\operatorname{Tr}(H^2)/(32\pi^2)\). With \(\gamma^{(1)}=0\), one-loop RG consistency therefore requires

\[
\boxed{\sum_p\beta_p\frac{\partial V_0}{\partial p}
=\frac{\operatorname{Tr}(H^2)}{32\pi^2},}
\qquad
\boxed{\sum_p B_p\frac{\partial V_0}{\partial p}=\frac12\operatorname{Tr}(H^2).}
\]

Here \(p\) includes all six quartics, all three quadratic coefficients **and \(C\)**. The polynomial equality is valid independently of whether a particular background has positive Hessian eigenvalues. Negative eigenvalues obstruct interpreting the logarithmic expression as a real stable-vacuum potential there; analytic continuation and its imaginary part do not supply the lifetime of the rotating inhomogeneous soliton.

`rg_symbolic.py` derives the left side from explicit symmetric-tensor contractions, then independently enters \(V_0\), differentiates its Hessian and expands the right side. Sparse polynomials use exact rational coefficients; no floating-point tolerance is involved in the equality. There are 50 expanded monomials per side, counting parameter and field powers. The residual has **zero monomials**. All 256 tensor components, all 16 mass components, the projected beta formulas, \(O(2)\) identities and arithmetic preflights passed, for **290 exact checks** in total.

The code also checks that omission is detectable. After setting the corresponding matching-zero parameters to zero, suppressing their beta functions leaves

\[
\left[\sum_p B_p\partial_pV_0-\tfrac12\operatorname{Tr}(H^2)\right]_{B_\alpha\,\mathrm{omitted},\alpha=0}
=-cdR^2D^2,
\]

\[
\left[\cdots\right]_{B_{m_X^2}\,\mathrm{omitted},\alpha=m_X^2=0}
=-c m_\phi^2R^2,
\quad
\left[\cdots\right]_{B_C\,\mathrm{omitted}}=-\tfrac12[(m_\phi^2)^2+2(m_X^2)^2+(m_D^2)^2].
\]

Thus keeping these coefficients zero at every scale would fail the one-loop identity for the already specified interactions.

## 5. Numerical matching values and generated operators

Mass and constant beta values below are normalized by \(v^2\) and \(v^4\), respectively. The remaining entries are dimensionless.

| Parameter | One-loop \(\beta\) at matching |
|---|---:|
| \(\lambda_\phi\) | \(1.39317894023\times10^{-5}\) |
| \(\lambda_X\) | \(1.39316627508\times10^{-6}\) |
| \(\lambda_D\) | \(1.14112983077\times10^{-7}\) |
| \(c\) | \(9.37220948692\times10^{-6}\) |
| \(d\) | \(4.23015941707\times10^{-8}\) |
| \(\alpha\) | \(1.26651479553\times10^{-8}\) |
| \(m_\phi^2/v^2\) | \(-3.78687923863\times10^{-6}\) |
| \(m_X^2/v^2\) | \(-1.26651479553\times10^{-6}\) |
| \(m_D^2/v^2\) | \(3.67289290703\times10^{-7}\) |
| \(C/v^4\) | \(6.33257397765\times10^{-7}\) |

For an illustrative first-order leading log \(p(\mu)=p(\mu_0)+\beta_p(\mu_0)\ln(\mu/\mu_0)\),

| \(\mu/\mu_0\) | \(\alpha\) | \(m_X^2/v^2\) |
|---:|---:|---:|
| 1/2 | \(-8.77881159659\times10^{-9}\) | \(8.77881159659\times10^{-7}\) |
| 1 | 0 | 0 |
| 2 | \(8.77881159659\times10^{-9}\) | \(-8.77881159659\times10^{-7}\) |

These are local leading-log illustrations, not a new flow integration or a fit. The negative ultraviolet-running bare \(m_X^2\) is not a demonstrated physical tachyon: a physical mass also involves the running vacuum, \(c\), self-energy and the relevant renormalization conditions. None of those effects can consistently be replaced by reading a bare coefficient in isolation.

## 6. Polynomial constant versus physical vacuum energy

For \(m_\phi^2<0\), \(\lambda_\phi>0\), on the \(x=y=D=0\) branch,

\[
v_{\rm cl}^2=-m_\phi^2/\lambda_\phi,\qquad
V_{0,\min}=C-\frac{(m_\phi^2)^2}{4\lambda_\phi}.
\]

The original tree matching has \(V_{0,\min}=0\), but the relation \(C=(m_\phi^2)^2/(4\lambda_\phi)\) is not preserved as a running identity. Direct differentiation gives

\[
B_{V_{0,\min}}=B_C-\frac{m_\phi^2}{2\lambda_\phi}B_{m_\phi^2}
+\frac{(m_\phi^2)^2}{4\lambda_\phi^2}B_{\lambda_\phi}.
\]

An additional exact check verifies this equals \(\operatorname{Tr}(H_{\rm vac}^2)/2\); the code multiplies by \(\lambda_\phi^2\) to keep the comparison polynomial. At matching the squared-mass eigenvalues are

\[
(0.02,0.01,0.01,0.0101)v^2,
\quad B_{V_{0,\min}}=0.000351005v^4,
\quad \beta_{V_{0,\min}}=2.22276512902\times10^{-6}v^4.
\]

This running tree value is canceled by the explicit scale dependence of the one-loop potential at the same order. It must not be identified with a physical change of the cosmological constant. To one-loop order the scalar vacuum value uses \(V_{0,\min}+V_1(v_{\rm cl})\); the shift of the stationary field contributes to the value only at higher order because \(V_0\) is stationary at \(v_{\rm cl}\).

In the covariant action, a separate gravitational coefficient \(\rho_\Lambda\) and the scalar polynomial \(C\) contribute through their **sum** to the constant density. Their separation is bookkeeping. The scalar beta \(B_C\) calculated here must be assigned to that total constant once; it must not be counted again as a second independent vacuum contribution. Full Standard Model and gravitational contributions are outside this calculation.

If \(C_{\overline{\rm MS}}(\mu_0)=0.0025v^4\) and \(\rho_{\Lambda,\overline{\rm MS}}(\mu_0)=0\) are both literal coefficient inputs, the finite scalar \(V_1(v_{\rm cl})\) is generically nonzero. To impose a zero **physical total** vacuum energy at one-loop instead requires a stated finite vacuum-energy matching condition at a chosen scale, propagated consistently by the RG. It is not legitimate to reset the running minimum to zero separately at every scale while claiming unchanged \(\overline{\rm MS}\) matching. This is a normalization/calibration issue for the existing vacuum parameter, not a new prediction or a mechanism solving the cosmological-constant problem.

## 7. Reproduction, outputs and remaining work

`rg_protocol.json` was saved before the tests. `rg_symbolic.py` uses only Python's standard library, so the check requires no symbolic-package installation. From the repository root:

```bash
python research/reality_framework/grut_candidate_1_quantum_20261018/rg_symbolic.py
```

The script refuses to overwrite `rg_results.json`. Reproduce in a fresh copy/output directory while preserving the frozen protocol and implementation. `rg_results.json` contains protocol/code hashes, exact checks, beta values and leading logs; `rg_polynomials.json` contains the exact coefficient lists for \(V_0\), its Hessian, both RG-identity sides and the zero residual. `rg_beta_functions.txt` and `rg_execution.log` provide compact equation/check records.

The result establishes one-loop pure-scalar polynomial closure and a nontrivial exact consistency identity for the fixed action. It does not make the constant-background potential a controlled approximation to the rotating finite-charge lump. That calculation requires the inhomogeneous, charge-constrained fluctuation determinant, collective/zero-mode treatment, derivative and frequency dependence, matching/subtraction and a consistent quantum state. Pole masses and lifetimes require additional on-shell/real-time calculations. The independent RG flow and numerical scale checks are owned by the integrating worker's separate `flow_*` files and are not modified here.

## Primary source and derivation status

Stephen P. Martin, *Two-loop effective potential for a general renormalizable theory and softly broken supersymmetry*, hep-ph/0111209v2: https://arxiv.org/html/hep-ph/0111209v2 . The scalar one-loop expression and scale conventions are in section 3; the general RG operator is equation (7.1), and the explicit derivative of the one-loop scalar function is equation (B.1). The tensor contractions and all model-specific formulas above were independently expanded by the supplied exact-arithmetic code; no two-loop result from that paper is claimed to have been evaluated here. Abstract/DOI record: https://arxiv.org/abs/hep-ph/0111209 .
