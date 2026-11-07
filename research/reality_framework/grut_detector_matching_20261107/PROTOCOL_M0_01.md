# M0-01 — frozen massless interval history comparison

**Status at writing:** equations, parameters, gates and analysis fixed before running the new reference or lattice histories. This file and the implementation must be committed before production. Existing PR #14 results are known; this is a new comparison, not a claim of ignorance of prior physics. Author: Codex, 7 November 2026. Zenodo remains paused.

## Claim and limit

Compare the same ideal Weyl-instrument histories in a staggered-fermion lattice sequence and an independently computed massless bosonized interval vacuum. The question is whether a bounded finite-resolution ladder reaches the declared numerical matching window. A pass applies only to this finite protocol set. Neither a pass nor massless solver verification certifies massive histories, the supremum in integrated equation (14), or a causal physical realization of the instantaneous extended instrument.

The lattice implementation reuses the already reviewed `ed/schwinger.py` Hamiltonian. The independent reference must not read lattice means, covariances, energies or probabilities. It uses continuum mode sums and the boundary value problem below. Keep production and controls in new files; never overwrite the original QED packet outputs.

## Physical dictionary

Set $\hbar=c=g=1$; all lengths and times below are in $g^{-1}$. The target is the massless Schwinger model with continuum $m_R/g=0$, theta zero and the neutral sector. Choose even site counts, conventional staggered Gauss law, zero incoming link flux and $Q=\sum_n Z_n=0$.

Use the site-endpoint interval convention $L=(N-1)a=5$, sites $z_j=ja$, link midpoints $z_{j+1/2}=(j+1/2)a$. Thus

$$
a=\frac5{N-1},\quad x=a^{-2},\quad
m_{\rm lat}/g=-a/8,\quad \mu=2\sqrt{x}\,m_{\rm lat}/g=-1/4.
$$

This is the published bulk mass-improvement prescription, not an assumed exact finite-interval mass calibration. The primary sequence is $N=10,14,18,20$. Additional *predeclared controls* at $N=10,20$ retain bare $m_{\rm lat}=0$. They do not replace a failing primary sequence.

For this conventional even neutral chain, the continuum scalar boundary values are $\phi(0)=\phi(L)=\sqrt\pi/4$. The ground state is Gaussian around

$$
\frac{\bar\phi(z)}{\sqrt\pi}
=\frac{\cosh[\mu_S(z-L/2)]}{4\cosh(\mu_S L/2)},
\qquad \mu_S=1/\sqrt\pi.
$$

This follows by solving $(-\partial_z^2+\mu_S^2)\bar\phi=0$ with those boundaries. Its nonzero mean must be included. The field dictionary is $E/g=\phi/\sqrt\pi$ at theta zero. Sources: [Okuda, eqs. (10), (13), (52), (66)](https://arxiv.org/html/2210.00297v2); [Dempsey et al., eq. (3) and section IV.5](https://arxiv.org/html/2206.05308v3). No modified Gauss law or theta shift is introduced.

The observable is $F=\int_0^L f(z)E(z)/g\,dz$, with

$$
f(z)=\frac{\exp[-(z-L/2)^2/(2\sigma^2)]}
{\int_0^L\exp[-(y-L/2)^2/(2\sigma^2)]\,dy},\qquad \sigma=0.8.
$$

On the lattice use the same Gaussian evaluated at link midpoints and normalize its discrete weights to sum to one. This midpoint discretization is an explicit approximation to the continuum observable, not an additional fitted coupling. Prepare each model's ground state with the declared parameters. Do not replace the continuum mean by the measured lattice mean or retune an instrument to make the histories agree.

## Instruments and protocols

$W=\exp(icF)$, $c=2$, and $K_s=(I+s e^{i\varphi}W)/2$, $s=\pm1$, $\varphi=0$. Outcomes are retained, and the measured state is propagated, so measurement disturbance is included. Protocol times are exactly `(0)`, `(0,1)`, `(0,2)` and `(0,1,2)`. Compute every outcome string at full numerical precision. No conditioning on rare outcomes and no fitted time rescaling are used.

Also compute the second measurement without the first, and a wrong-reference control in which only the continuum mean is set to zero. The latter intentionally changes the preparation and is an implementation diagnostic, not a competing matched physical model. Its outcome is reported without inventing a separate numerical threshold after inspection.

## Independent reference and ultraviolet control

Let $u_n=\sqrt{2/L}\sin(n\pi z/L)$, $\Omega_n=\sqrt{(n\pi/L)^2+1/\pi}$, $f_n=\int f u_n$, and $h_n=f_n/\sqrt{2\pi\Omega_n}$. Then

$$
\bar F=\int f\bar\phi/\sqrt\pi,\qquad
G(t)=\sum_{n\ge1}h_n^2e^{-i\Omega_n t}.
$$

Derive ordered Weyl expectations from displacement-operator multiplication and the Gaussian vacuum. The implementation must retain commutator phases. Verify that algorithm against direct finite-oscillator state evolution at occupation cutoffs 32 and 48, with $\Omega=0.8$, $h=0.3$, mean $0.17$, times `(0,0.7,1.9)` and phases `(0.2,-0.3,0.4)`. This is a code verification case, not the lattice validation experiment.

Use $M=8192$ continuum modes; repeat at $M=16384$. Compute sine coefficients by weighted adaptive quadrature, absolute tolerance $10^{-12}$, and repeat the reference with $10^{-13}$. The quadrature estimates are numerical estimates, distinct from the analytic mode-tail bound.

For $B=|f(0)|+|f(L)|+\int|f'|=2f(L/2)$, integration by parts gives $|f_n|\le\sqrt{2/L}B/k_n$. Therefore

$$
V_{>M}=\sum_{n>M}h_n^2
\leq \frac{B^2L^2}{2\pi^4 M^2}.
$$

An ordered word of at most $2n$ Weyl factors has omitted exponent bounded by $2n^2c^2V_{>M}$. Interpolating its Gaussian expectation between retained and full modes preserves nonpositive real exponent. Summing the $4^n$ expansion terms and $2^n$ outcomes yields

$$
\operatorname{TV}(p_M,p_\infty)\leq 2^n n^2c^2 V_{>M}.
$$

This bounds mode truncation for the declared finite histories; it is not a uniform arbitrary-protocol result.

## Solver, refinement and gates

Lattice: sparse eigensolver for the lowest two states, deterministic starting vector seed `20261107+N`, tolerance $10^{-11}$; use the physical $H/g=W/(2\sqrt{x})$. Propagate with SciPy sparse `expm_multiply`. Independently rerun the $N=14$ primary point with eigensolver tolerance $10^{-13}$ and two exponential substeps per interval. Report eigenvector residuals, gap, state normalization and maximum probability normalization error. Numerical agreement and residuals are diagnostics, not a general certified eigenstate-error theorem.

Fixed gates:

| Gate | Registered criterion |
|---|---|
| Reference implementation | Direct-oscillator versus Gaussian histories differ by TV $<10^{-10}$ at cutoff 48; the 32-to-48 change is also $<10^{-10}$ |
| Continuum reference control | Three-time analytic mode-tail bound $<10^{-6}$; doubling modes and tightening quadrature each change probabilities by TV $<10^{-6}$ |
| Lattice numerical diagnostics | Eigenvector residuals $<10^{-9}$; probability normalization error $<10^{-9}$; tighter $N=14$ replay changes all protocol distributions by TV $<10^{-8}$ |
| Finite-resolution matching | Each primary $N=20$ protocol has TV $\le10^{-3}$ relative to the independent continuum reference |
| Resolution diagnostic | Each primary protocol's distance is nonincreasing from $N=14$ to 18 to 20; the 18-to-20 distributions differ by TV $\le2.5\times10^{-4}$ |

No continuum extrapolation, optimized endpoint offset, fitted detector normalization or target-dependent mass adjustment is allowed in M0-01. Report mean, covariance, gap and history differences separately to locate a disagreement. Refinement differences alone are not rigorous estimates of every omitted systematic. Even passing these gates would be numerical evidence for this protocol set, not a complete physical error certification: finite-interval parameter matching and ideal-detector realization still need their own bounds.

## Resources, outcomes and continuation

Run points sequentially with one BLAS thread. Per-point wall cap 300 seconds, total lattice production cap 1200 seconds, no sites above 20 in this card. Archive any timeout, failed gate or exception with inputs and code hashes. No result-dependent extension counts as this same test.

Classify separately: implementation failure; identified dictionary mismatch; inadequate resolution; or unresolved disagreement after controlled comparison. A persistent numerical discrepancy is not automatically new physics. An analytically controlled massless result does not bound massive-history errors. Build the massive reference under a separate card only after reporting this outcome and its remaining limitations.

Reproducibility here means replaying this implementation. Evaluator reliability on unfamiliar theories requires its own frozen external test set and calibrated alternatives; this packet makes no such certification.
