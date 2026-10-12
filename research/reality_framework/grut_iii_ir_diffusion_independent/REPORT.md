# GRUT III v13/v14: independent asymptotic derivation of the sigma-free IR coefficient

2026-10-11. Branch \`codex/grut-iii-independent-ir-diffusion\`, forked from Claude's **GRUT III revision 2.4** [e123dc5](https://github.com/ryangrvr/TestingGRUT/commit/e123dc5d0fd10c8a986109dd5c38f4d44095d18c). This is an additive independent crosscheck, not a rewrite of GRUT III, not a new postulate, and not a Stage-3 Card 1. The frozen CP programs and other partner branches are untouched.

## Answer: the gamma-function coefficient has an independent derivation

The fixed-background mathematical claim of \`v13_zero_mode_temperature.py\` and \`v14_massless_sector_first_order.py\` is reproducible. Rather than expanding Claude's claimed gamma identity in isolation, derive it **from the large-separation asymptotic of the exact massless power-law hypergeometric two-point function**, with a different operation order.

Set H(0)=a(0)=1 and constant epsilon with 0<epsilon<<1. Define
\[
a(t)=(1+\epsilon t)^{1/\epsilon},\quad H(t)=(1+\epsilon t)^{-1},
\quad \tau(t)=-[(1-\epsilon)a(t)H(t)]^{-1},
\quad \delta=\frac{\epsilon}{1-\epsilon}.
\]
The same exact analytic-continuation scheme supplies
\[
W(t,t')=\frac{(1-\epsilon)^2 H(t)H(t')}{16\pi^2}
\Gamma(3+\delta)\Gamma(-\delta)
\,{}_2F_1(3+\delta,-\delta;2;w)+C_{\rm IR},
\]
\[
w=1+\frac{[\tau(t)-\tau(t')]^2}{4\tau(t)\tau(t')}.
\]
For t->infinity with earlier t' fixed, |w|->infinity and the leading Gauss-hypergeometric connection formula is
\[
{}_2F_1(3+\delta,-\delta;2;w)
\sim\frac{\Gamma(3+2\delta)}{\Gamma(3+\delta)\Gamma(2+\delta)}
(-w)^\delta.
\]
The physical **symmetric real part** across the timelike cut contributes \(\cos(\pi\delta)\). Because
\(H(t)|\tau(t)|^{-\delta}=H(t')|\tau(t')|^{-\delta}\),
\[
\lim_{t\to\infty}\Re W(t,t')
=\frac{H(t')^2(1-\epsilon)^2}{16\pi^2}
\,\frac{\Gamma(-\delta)\Gamma(3+2\delta)}
{\Gamma(2+\delta)4^\delta}\cos(\pi\delta)
+ C_{\rm IR}.
\]
Assume C_IR is a fixed, time-independent subtraction in the specified state so that differentiating with respect to t' removes it. Then \(\partial_{t'}H(t')^2=-2\epsilon H(t')^3\), giving
\[
\boxed{\frac{D}{H^3}
=-\frac{\epsilon(1-\epsilon)^2}{8\pi^2}
\frac{\Gamma(-\delta)\Gamma(3+2\delta)}{4^\delta\Gamma(2+\delta)}
\cos(\pi\delta).}
\]
This is the **independently derived form**, with no stochastic coarse-graining window or assumed diffusion parameter.

Use \(\Gamma(1-\delta)=-\delta\Gamma(-\delta)\), \(\delta=\epsilon/(1-\epsilon)\), and the **duplication identity**
\[
\frac{\Gamma(3+2\delta)}{\Gamma(2+\delta)}
=\frac{2^{2+2\delta}}{\sqrt\pi}\Gamma(3/2+\delta).
\]
The result becomes exactly the claimed closed form
\[
\boxed{
\frac{D}{H^3}=
\frac{\sqrt\pi}{2\pi^3}(1-\epsilon)^3
\Gamma(3/2+\delta)\Gamma(1-\delta)\cos(\pi\delta).
}
\]
Logarithmic differentiation at epsilon=0 uses
\(\psi(3/2)-\psi(1)=2-2\ln2\), so
\[
\left.\frac{d}{d\epsilon}\ln(D/H^3)\right|_0
=-3+\psi(3/2)-\psi(1)
=\boxed{-(1+2\ln2)}=-2.3862943611198906\ldots
\]
This closes the *mathematical re-derivation requested in Claude's revision-2.4 closeout*, under the explicit analytic-continuation/IR assumptions. It is not an independent proof that every physically prepared state has the same coefficient.

## Direct finite-time test of the unexpanded two-point function

An independent mpmath implementation evaluated
\( \partial_{t'}\Re W(t,t')|_{t'=0} \)
**without first taking** t->infinity or using the closed-form D for the differentiation. Taking t=10,000 (and H(0)=1) yields:

| epsilon | Direct derivative | Independently derived asymptote | Relative difference |
|---:|---:|---:|---:|
| 0.05 | 0.022204676618995 | 0.022204676618995 | below 1e-40 |
| 0.10 | 0.018882672488711 | 0.018882672488711 | below 1e-40 |
| 0.20 | 0.011654110728280 | 0.011654110728280 | below 1e-24 |
| 0.30 | 0.003304455543328 | 0.003304455543328 | about 5.5e-17 |

The independently derived hypergeometric-asymptotic and gamma-reduced expressions agree to ~45 digits at each of these values. The companion [check.py](check.py) implements both paths and numerical controls, **without importing any GRUT III calculation script**. Equivalent equations and the finite-time derivative tests were independently executed in the analysis environment. The exact committed script should additionally be run verbatim by a partner; do not overstate local equivalent execution as a repository CI run.

### Exact scope warning: this is not a positive diffusion constant for arbitrary epsilon

The exact expression contains \(\cos(\pi\delta)\). It vanishes at delta=1/2 (epsilon=1/3), and for epsilon=0.4 is **negative** (~-0.0089504573 at H=1). For 0<epsilon<1/3 the displayed branch is positive; the intended perturbative slow-roll regime epsilon<<1 is safely inside it.

A mathematical derivative of a symmetrized two-time function may be well-defined even when its stochastic *positive white-noise diffusion* interpretation is not. Therefore **do not extend the "sigma-free diffusion coefficient" headline to arbitrary constant epsilon or to an arbitrary state/IR prescription**. An altered initial zero mode or a time-dependent infrared subtraction could alter the derivative. This does not threaten its first-order expansion in the declared BD-like analytic-continuation control.

## What the independent check does NOT establish about thermality

- A state with \(\partial_{\bar t}W(\bar t,\Delta)\ne0\) is not globally stationary, so it cannot have a conventional global single-parameter KMS state simply because the **fixed-mean-time Wigner spectrum** has a detailed-balance relation to O(epsilon). Claude explicitly notes this limitation after audit. Keep "Wigner detailed balance" as the headline, not unqualified thermality.
- Analytic midline symmetry of a function of \(\cosh^2(H\Delta/2)\) is an implementation clue, **not a theorem that an arbitrary function has the requisite analytic strip, positivity and distributional boundary conditions of KMS**. Claude's exact-function strip check helps, but it remains a check in the chosen continuation.
- Polynomial contributions in \(\Delta\) correspond to distributions supported at omega=0 in the infinite-time spectrum. Removing them to compare nonzero frequencies is an explicit comparison convention, not an innocuous physical operation for an instrument with finite switching: its spectral resolution convolves zero-frequency contributions into the response. A real finite-time Unruh–DeWitt detector must be assessed with its switching and coupling, not asserted thermal by the Wigner result.
- The resulting \(-[1+2\ln 2]\) is standard free-field QFT on supplied power-law FRW. It supplies no new independent GRUT postulate, prediction, geometry selection or Born mechanism. It **should not** be cited as independent experimental confirmation.

## Literature crosscheck and priority

A targeted search of equations and publication abstracts on 2026-10-11 did not identify an already-published exact occurrence of this *specific* gamma combination or the v14 first-order \(r(\omega)\) function. This is **not a priority/novelty finding**: abstract/index searches do not exclude a result contained in a paper's appendices, another notation or earlier derivation. In particular the Wigner shape \(r(\omega)\) has *not* been independently reproduced in this review.

Established related prior art:
- Cho, Kim & Kitamoto, [*Stochastic Dynamics of Infrared Fluctuations in Accelerating Universe*](https://arxiv.org/abs/1508.07877) (2015): constant/variable slow-roll stochastic noise, state and drift; this is the precise underlying stochastic source discussed by v13.
- Kandrup, [*Stochastic inflation as a time-dependent random walk*](https://doi.org/10.1103/PhysRevD.39.2245) (1989): stochastic/IR behavior and power-law inflation, including limits to fluctuation–dissipation language.
- Garbrecht & Prokopec, [*Unruh response functions for scalar fields in de Sitter space*](https://doi.org/10.1088/0264-9381/21/21/016) (2004): detector response, IR influence and separate measurement conventions.
- Garbrecht & Prokopec, [*Energy density in expanding universes as seen by an Unruh detector*](https://doi.org/10.1103/PhysRevD.70.083529) (2004): a separately defined detector observable with different conclusions about thermal interpretation.
- Cruces & Germani, [*Stochastic inflation at all order in slow-roll parameters: foundations*](https://arxiv.org/abs/2107.12735) (2021): why mode-based stochastic approximations require their own dynamical scope.
- [*Renormalisation and matching of massless scalar correlation functions in Soft de Sitter Effective Theory*](https://doi.org/10.1007/JHEP07(2026)177) (2026): extended IR EFT matching framework and source-state qualifications.

## Decision for parallel research

Claude's GRUT III revision 2.4 can be retained as a **conventional QFT baseline**, with the stated fixed-mean-time convention and corrected infrared qualifications. This independent derivation closes the outstanding narrow *gamma-formula audit*, not the entire priority check for \(r(\omega)\).

A new KMS-breaking postulate is **not forced** by these conventional results and should not be invented simply to generate novelty. The constructive microscopic program (CP-11 and successors) still owes: dynamically selected geometry, physical local interactions, an extensive thermodynamic phase, and an actual novel tested consequence. Neither an effective temperature nor a first-order diffusion coefficient earns those.

No changes to any frozen research report or to Claude's branch.
