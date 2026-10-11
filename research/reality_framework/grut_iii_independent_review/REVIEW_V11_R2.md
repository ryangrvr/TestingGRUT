# GRUT III revision 2: independent first-order quasi-de Sitter check

2026-10-11. **Cross-branch review only**: author's baseline [GRUT III revision 2, 24ebc104](https://github.com/ryangrvr/TestingGRUT/tree/24ebc104ecd5aab1f591c5d4042a583ded5e2929/GRUT_III). This note is placed on \`codex/grut-causal-state-review\`; it does not change Claude's branch, CP-06, frozen constructions, or PR #12. Author independent verification of v7–v11 remains outstanding.

## Substantive correction to v11 / D11 / O5

The exact conformal power-law FRW worldline kernel in \`v11_quasi_ds_kms_defect.py\` yields a **pure shift in de Sitter H at first order in epsilon** when expanded at fixed midpoint H and finite (including noncoincident complex) relative time. The numerical frequency-dependence observed at finite epsilon cannot be called a demonstrated O(epsilon) nonthermal shape without an explicit epsilon->0 extrapolation *after* subtracting that exact first-order shift.

Let epsilon=1/p, x=H_mid (t-t'), a(t)=(t/tbar)^p, tbar=p/H_mid, and u±=1±epsilon*x/2. Direct substitution into the exact conformal vacuum Wightman function gives

\[
W_\epsilon(x)=-\frac{H^2(1-\epsilon)^2}
{4\pi^2 (u_+u_-)^{1/\epsilon}
\left(u_+^{1-1/\epsilon}-u_-^{1-1/\epsilon}\right)^2}.
\]

As epsilon->0, \`W_0(x)=-H^2/[16*pi^2*sinh^2(x/2)]\`, with an understood time ordering/i0 prescription.

At fixed finite x the expansion is

\[
\frac{W_\epsilon(x)}{W_0(x)}
=1+\epsilon\left[x\coth(x/2)-2\right]+O(\epsilon^2).
\]

For a de Sitter kernel at effective Hubble rate H_eff=H(1-epsilon), define
\[
W_{\rm dS}(x;H_{\rm eff})=
-\frac{H_{\rm eff}^2}{16\pi^2
\sinh^2[(H_{\rm eff}/H)x/2]}.
\]

Because \(-H\partial_H W_{\rm dS}/W_0=x\coth(x/2)-2\) with physical Delta fixed,

\[
\boxed{W_\epsilon(x)=W_{\rm dS}(x;H(1-\epsilon))+O(\epsilon^2).}
\]

Moreover direct expansion to the *next* order gives

\[
\boxed{
\frac{W_\epsilon(x)}{W_{\rm dS}(x;H(1-\epsilon))}
=
1+\epsilon^2\left[\frac{x^2}{4}
-\frac{x^3}{12}\coth\frac{x}{2}\right]
+O(\epsilon^3).
}
\]

The first additional time-separation-dependent nonthermal shape therefore starts at O(epsilon²) for this local conformal correlator. It is not a free parameter or new GRUT prediction; it is conventional QFT in power-law FRW. If a controlled Fourier/Wigner limit may be interchanged with the expansion, the corresponding **first-order** detailed-balance temperature is T_eff=H(1-epsilon)/(2*pi), with no O(epsilon) frequency-dependent defect after this shift. The full global Wigner statement still needs boundary/switching control: power-law FRW is defined for t>0; a midpoint Fourier transform over all Delta does not automatically define one physical stationary detector spectrum, and v11's contour-shifted finite interval can introduce edge contributions. Avoid claiming a universal exact thermal state.

## Independent numerical check

Accompanying [check_v11_powerlaw_expansion.py](check_v11_powerlaw_expansion.py) evaluates the **exact unexpanded** W at complex, noncoincident x = 1-i, 0.75-i*pi, 2-i*pi, 3-2i and epsilon=0.02, 0.01, 0.005. It compares the limiting O(epsilon) and O(epsilon²) coefficients with their independent analytic expressions. The O(epsilon²) coefficient residual decreases linearly with epsilon for every x tested.

Executed in a separate Python/mathematical environment before uploading. Expected final lines:

\`\`\`
PASS: first-order correction is an effective H shift;
      first non-temperature-shift shape appears at O(epsilon^2)
\`\`\`

This crosscheck was written independently of the \`v11\` script and does not reuse its quadrature implementation. It does not certify the global Fourier/Wigner extension or detector response, so those remain open.

## Other specific items requiring adversarial verification

1. **T6 realizability scope** (\`07_ACCESS_RESOLVED_RESPONSE.md\`): A stationary, causal, passive linear response with suitable regularity, dispersion subtractions and a positive oscillator spectral density can be represented by a harmonic bath at the linear-response level. A **local Lorentz-covariant field** does not follow from an arbitrary temporal kernel; Källén–Lehmann requires additional Poincaré covariance, positivity and spectral conditions. Mori–Zwanzig projection does not itself imply a unique physical bath. Contact terms and zero-mode effects must be tracked. Do not broaden the theorem to "every causal stable kernel has a local integrated-out field" without those premises.

2. **v10 diffusion/thermal statement**: The code demonstrates the \`D*min(t,t')\` behavior as a **large-separation asymptote** of the massless-minimally-coupled two-point function, and stationarity of differentiated correlators. It does not establish equality with an exact Brownian covariance at arbitrary separations. Imaginary-time periodicity of a \`log(sinh²)\` expression alone also requires a branch/i0 and analyticity check to certify KMS detailed balance.

3. **v9 Newton-constant accounting**: The reported 1/G induction coefficients are cutoff/scheme and field-content dependent. A cutoff inferred larger than Planck for minimal Standard Model fields shows failure of a one-loop dominated *particular* Sakharov scenario, not a universal no-go for emergence of G. This is largely acknowledged in its text.

4. **v7 smearing**: The stress response retains a UV 1/t edge with purely spatial smearing; it should be treated as a distribution with physically defined temporal switching, not an integrable phenomenological memory kernel by default.

## Research division and freeze

GRUT III is a valuable **conventional correspondence and falsifier baseline** for a microscopic causal theory. The causal-history program's unsolved problem is more fundamental: derive a common microscopic action, quantum history measure and physical records; no number of baseline correlator identities by itself selects those. Reproduce/repair the v11 O(epsilon) statement first, then test a physically dressed/finite-switching gravitational response. Keep all distinctive predictions as explicit residuals against a fully matched conventional baseline, with no refitting.

This review is a negative/precision result, not CP-07, a new framework, or a proposal to merge branches.
