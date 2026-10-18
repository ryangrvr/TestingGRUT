# Worked regime: populated structure, propagation, and gravity

**GRUT Candidate 1 · 18 October 2026**

The displayed action in `GRUT_CANDIDATE_1.md` supports a numerically resolved, populated, charge-bound classical configuration. Its organization field changes a reporter's mass and scattering phase. Its conserved stress sources a weak gravitational metric with a calculated clock shift and light deflection. These are calculations within the proposed action, not observations of nature.

The demonstrated hierarchy is: stationary interacting scalar solution at zeroth order in gravity; linearized metric from that solution's complete stress; linear reporter propagation on the fixed scalar background; an explicitly supplied ideal readout instrument. No coupled quantum lifetime or microscopic laboratory apparatus was calculated.

## 1. Frozen physical parameters and supplied preparation

Let `g=g_X`, `rho=g v r`, `F=f/v`, `S=varphi/v`, `Omega=omega/(g v)`, `a=lambda_phi/g²=1`, and `b=lambda_X/g²=0.1`. Use

\[
g=0.1,\quad\lambda_\phi=0.01,\quad\lambda_X=0.001,\quad
h_D=0.01,\quad m_{D0}=gv,\quad\lambda_D=0.001,
\quad v/M_P=10^{-5}.
\]

The asymptotically flat benchmark has renormalized `rho_Lambda=0`, as a boundary/matching choice. Thus `Gv²=10^−10/(8 pi)`.

Prepare

\[
X(t,r)=f(r)e^{-i\omega t}/\sqrt2,\qquad \phi=\varphi(r),\qquad D=0,
\]

with the regular nodeless profile described below. This is supplied charged matter. A large-mean-charge semiclassical state is its interpretation; the noninteger classical charge is not an exact integer particle count or an exact coherent state of fixed charge. A charge-projected state requires a different quantum description with `⟨X⟩=0`.

The physical couplings were chosen before solving. Six preliminary high-amplitude boundary-value starts failed and are archived. A fixed-frequency `Omega=0.8` bag-like seed then converged. Its central amplitude `F(0)=1.979050307129842` was frozen for all verification refinements. This is a successful preparation search, not evidence of generic spontaneous formation, laboratory preparation acceptance, or statistical pre-registration. No physical coupling was changed to rescue the solve.

## 2. Why this potential admits a binding class

Write `y=|X|²`. Minimization over a spatially uniform phi gives

\[
\phi_{\min}^2=\max(0,v^2-2g^2y/\lambda_\phi).
\]

On the nondepleted branch,

\[
U_{\min}=g^2v^2y+(\lambda_X-g^4/\lambda_\phi)y^2.
\]

Allowing phi to respond creates an attractive effective quartic when `g^4>lambda_phi lambda_X`, while the full two-field potential stays nonnegative and bounded. On the depleted branch, minimizing `U/y` gives

\[
f_0^2=v^2\sqrt{\lambda_\phi/\lambda_X},\qquad
\omega_0^2=v^2\sqrt{\lambda_\phi\lambda_X}.
\]

Hence `Omega_min=(a b)^(1/4)=0.5623413<1`. The candidate localization/binding interval is `Omega_min<Omega<1`. This analytic argument supplies a binding class; the actual finite-charge solution and its stability require separate checks. Eliminating phi algebraically in the wall would discard its gradient energy and is not the calculation performed here.

The underlying two-field charge-supported soliton mechanism is established prior work. The positive `lambda_X |X|^4` term changes the original Friedberg–Lee–Sirlin potential, so original FLS numerical bounds and asymptotics are not imported. See [Heeck and Sokhashvili, 2023](https://arxiv.org/html/2303.09566v2).

## 3. Solve the same field equations

The flat, stationary equations are

\[
F''+\frac2\rho F'=(S^2+bF^2-\Omega^2)F,
\qquad
S''+\frac2\rho S'=[a(S^2-1)+F^2]S.
\]

Primes denote `d/d rho`. Conditions are `F'(0)=S'(0)=0`, `F(infinity)=0`, `S(infinity)=1`, and the frozen `F(0)` above. `Omega` is solved as a boundary-value parameter. `scipy.integrate.solve_bvp` handles the regular radial singular term explicitly; domain and tolerance refinements are separate.

Define

\[
\bar Q=4\pi\Omega\int_0^\infty\rho^2F^2d\rho,\qquad Q=\bar Q/g^2,
\]

\[
\bar E=4\pi\int_0^\infty\rho^2d\rho\left[
\frac{F'^2+S'^2+\Omega^2F^2+S^2F^2}{2}
+\frac b4F^4+\frac a4(S^2-1)^2\right],\quad E=\frac vg\bar E.
\]

`R90` encloses 90% of energy, not charge. Results from `code/results.json` are:

| Domain `rho_max` | BVP tolerance | `Omega` | `Qbar` | `Ebar/Qbar` |
|---:|---:|---:|---:|---:|
| 30 | 10^−5 | 0.7999999545 | 399.4385720 | 0.9110864932 |
| 30 | 10^−7 | 0.8000000014 | 399.4384639 | 0.9110865192 |
| 40 | 10^−8 | 0.8000000019 | 399.4384625 | 0.9110865196 |
| 50 | 10^−9 | 0.8000000020 | 399.4384624 | 0.9110865197 |

The final solve has maximum collocation RMS residual `1.0×10^−9` and negligible boundary residual. The asymptotic F tail decays with exponent `sqrt(1−Omega²)≈0.6`; the homogeneous S scale is `sqrt(2a)≈1.414`, while the F² source can drive a slower S deviation with exponent about 1.2. Domain comparison checks truncation empirically; these are not interval-certified errors.

| Calculated quantity | Value |
|---|---:|
| `S(0)` | 0.0490799789 |
| `Qbar` | 399.4384624 |
| `Ebar` | 363.9229986 |
| `R90bar` | 5.1034322 |
| Mean/classical charge `Q` at `g=.1` | 39,943.8462 |
| Energy `E/v` | 3,639.2300 |
| Radius `R90 v` | 51.03432 |
| `E/(m_X Q)` | 0.91108652 |

The last ratio excludes complete conversion of this leading classical energy into free X charge carriers at the same charge. It does not prove stability against every charged fragmentation or quantum decay.

## 4. Perturbation checks and their precise scope

At fixed charge minimize

\[
\bar E_{\bar Q}[F,S]=\frac{\bar Q^2}{2I}
+4\pi\int\rho^2d\rho\,[\tfrac12(F'^2+S'^2)+\bar U],
\quad I=4\pi\int\rho^2F^2d\rho,
\]

where `Ubar=a(S²−1)²/4+S²F²/2+bF⁴/4`. The amplitude Hessian, after radial unitary transformation, has local potential matrix

\[
V_{\rm amp}=\begin{pmatrix}S^2+3bF^2-\Omega^2&2FS\\2FS&a(3S^2-1)+F^2\end{pmatrix}.
\]

Add `−d²/d rho²+ell(ell+1)/rho²` on each diagonal. In the spherical sector include the essential positive rank-one term `4 Omega² |q><q|/I`, where `q=(sqrt(4 pi) rho F,0)` and the inner product is `d rho`. Omitting it gives a negative eigenvalue about `−0.188855`; that is the wrong fixed-charge problem.

| Amplitude sector | Lowest eigenvalue, `L=40,N=3200` | Interpretation |
|---|---:|---|
| `ell=0`, fixed charge | 0.36696654 | Positive finite-box value; continuum threshold is 0.36 |
| `ell=1` | −2.66211×10^−6 | Translation zero mode with overlap 0.999999999987 |
| `ell=2` | 0.28548987 | Positive tested nonradial mode |
| `ell=3` | 0.39051344 | Positive tested sector |

The translation artifact decreases proportional to mesh spacing squared across the three grids. Higher amplitude `ell` add a positive centrifugal term relative to `ell=1`; the numerical low sectors support the expected absence of an extra negative amplitude mode. The displayed `ell=0` number is not an infinite-space spectral gap.

The phase operator is

\[
H_{\rm phase}=-d^2/d\rho^2+S^2+bF^2-\Omega^2.
\]

It has the positive, nodeless zero solution `rho F`; the ground-state factorization gives continuum nonnegativity, with higher `ell` adding a positive centrifugal term. The computed `ell=0` eigenvalue approaches zero as `h²`: `−4.05×10^−6`, `−2.59×10^−6`, `−1.01×10^−6`. These are symmetry zero-mode discretization errors, not measured growth rates.

This is numerical evidence for local fixed-charge classical energetic stability, supported by symmetry identities. No full nonlinear perturbation evolution, all-partition fission minimization, quantum fluctuation determinant, or lifetime was executed. Small finite-grid errors are identified by refinement and analytic modes, not erased by declaring every negative number harmless.

## 5. Gravity from the complete stress

Use dimensionless stress in units `g²v⁴`. Set

\[
T_g=\tfrac12(F'^2+S'^2),\quad T_t=\tfrac12\Omega^2F^2,
\quad e=T_g+T_t+\bar U,
\]

\[
p_r=T_g+T_t-\bar U,\qquad p_t=-T_g+T_t-\bar U.
\]

The stationary identity is `p_r'+2(p_r−p_t)/rho=0`. The final maximum dimensionless residual is about `4.7×10^−10`. Integrated stresses satisfy

\[
T_{\rm grad}+3(U_{\rm total}-T_{\rm time})=0
\]

to `1.83×10^−12` relative to energy. Consequently the integrated active mass agrees with energy at that accuracy.

In the weak static metric

\[
ds^2=-(1+2\Phi)dt^2+(1-2\Psi)d\mathbf x^2,
\]

Einstein's linearized equations give

\[
\nabla_\rho^2\Phi=4\pi Gv^2(e+p_r+2p_t),\qquad
\nabla_\rho^2\Psi=4\pi Gv^2 e.
\]

The local active source is `2 Omega²F²−2 Ubar`; replacing it by energy density alone changes the interior answer. For either dimensionless source `s`, the code evaluates

\[
P_s(\rho)=-Gv^2\left[\frac{4\pi}{\rho}\int_0^\rho x^2s(x)dx
+4\pi\int_\rho^\infty xs(x)dx\right].
\]

| Observable | Calculated leading value |
|---|---:|
| `2GE/R90` | 5.67463×10^−10 |
| `Phi(0)` | −6.68813×10^−10 |
| `Psi(0)` | −5.23912×10^−10 |
| Central clock redshift to infinity `z≈−Phi(0)` | 6.68813×10^−10 |
| Light deflection at `b g v=20` | 2.89601×10^−10 radians |

Exterior stress conservation gives `Phi=Psi=−GE/r` at leading order. Integrating the transverse gradient of `Phi+Psi` gives the reported bending, agreeing with `4GE/b` at `b g v=20` to about `3×10^−10` relative within this computation. That agreement is a consistency check on adopted Einstein gravity, not a new gravitational prediction beyond it. Slow ordinary test bodies respond to `−grad Phi`; clocks and light depend on the same metric. D also feels its scalar mass gradient and is not a clean pure-gravity test inside the lump.

The scalar profile has not been recomputed in its gravitational field. This is a consistent first-order expansion: its metric source uses the zeroth-order solution; profile feedback changes the metric at next order. The very small potentials control the classical weak-field expansion, subject to the absence of a singular branch response. No exact nonlinear Einstein–scalar solution or rigorous coefficient bound was calculated.

## 6. Reporter: from this preparation to an ideal measurement

The vacuum reporter mass is `m_D,infinity²=m_D0²+h_D²v²`. Write its dimensionless momentum as `k`, with energy `epsilon=E_D/(gv)=sqrt(k²+1.01)`. The linear s-wave obeys

\[
u''+[k^2-V(\rho)]u=0,\qquad
V=0.01(S^2-1),\qquad u(0)=0.
\]

The well is attractive and the local reporter mass squared stays positive. With `u~sin(k rho+delta)`, the elastic channel is `S_0=exp(2i delta)`. The variable-phase equation

\[
\delta'(\rho)=-\frac{V(\rho)}k\sin^2[k\rho+\delta(\rho)],\qquad\delta(0)=0
\]

was integrated directly and checked against an independent radial wave fit. Their final differences are below `3×10^−11` radians at all three momenta. The first Born result and a controlled exact-profile remainder are

\[
\delta_B=-\frac1k\int V\sin^2(k\rho)d\rho,\qquad
|\delta-\delta_B|\leq\frac1{2k^2}\left(\int|V|d\rho\right)^2.
\]

The bound follows from `|sin²(x+y)−sin²x|≤|y|` and the integrated variable-phase bound. It is analytic for an exact potential; numerical profile/ODE/quadrature errors must still be added. Here `integral |V|=0.0406684994`.

For a balanced ideal reference arm with phase `beta=pi/2`, unit selected-mode probability, and otherwise identical modes,

\[
P_+=\frac{1-\sin(2\delta)}2.
\]

| Dimensionless `k` | `delta` | Absolute Born error | Born remainder bound | Ideal `P+` |
|---:|---:|---:|---:|---:|
| 3 | 0.0067765110 | 1.568×10^−6 | 9.188×10^−5 | 0.49322370 |
| 5 | 0.0040665116 | 3.384×10^−7 | 3.308×10^−5 | 0.49593353 |
| 7 | 0.0029047695 | 1.233×10^−7 | 1.688×10^−5 | 0.49709525 |

These are conditional classical linear/ideal elastic-channel predictions. Full finite-mass target scattering, recoil, entanglement and inelastic losses have not been solved.

### An actual finite packet was evolved

An incoming radial positive-frequency Gaussian packet with center `rho=22`, envelope width 3 and mean inward momentum 5 was evolved by

\[
i\partial_\tau u=\sqrt K\,u,\qquad K=-\partial_\rho^2+1+0.01S^2,
\quad \tau=gv t,
\]

on `[0,80]` with Dirichlet endpoints. This is the positive-frequency representation of the linear local Klein–Gordon problem; it does not supply an exact local relativistic position observable or compact wavefront. The supplied ideal instrument measures the radial L2 projector on `R=[12,60]` at `tau=40`. Let `p_s,p_ref` be its separate signal/reference selection probabilities and `C_R=<u_ref|P_R|u_s>`. Then

\[
P_\pm(R)=\tfrac14[p_s+p_{ref}\ \pm2\operatorname{Re}(e^{i\beta}C_R)],
\qquad p_{select}=\tfrac12(p_s+p_{ref}),
\quad P_{+|R}=P_+(R)/p_{select}.
\]

Thus phase, visibility, delivery and conditional output remain distinct.

| Radial grid points | Mean radius at `tau=40` | Selected probability | Joint `P+(R)` | Conditional `P+ given R` |
|---:|---:|---:|---:|---:|
| 700 | 15.6034 | 0.953610 | 0.472698 | 0.495693 |
| 1400 | 16.8074 | 0.988089 | 0.489962 | 0.495868 |
| 2800 | 17.1109 | 0.991960 | 0.491923 | 0.495910 |

Global norm error is below `10^−13`. That checks the chosen unitary grid evolution, not its continuum accuracy. The last grid change is `0.00387` in selection, `0.00196` in joint P+, and `0.0000421` in conditional P+. Quote this packet result only at modest precision: about 99% delivery and 49% joint plus output for the chosen time/region. No rigorous continuum error enclosure or domain refinement of this packet was performed. The converged monochromatic calculation above supplies the more accurate propagation benchmark.

The finest packet overlap phase is about `0.00818`; its finite bandwidth and selected readout make it different from the monochromatic `2 delta=0.00813302`. No retention probability for an internal message is inferred from packet delivery.

## 7. Information and resource accounting

For equiprobable hypotheses “uniform vacuum” and “this prepared lump,” the ideal normalized monochromatic arm states at `k=5` have overlap magnitude `|cos delta|`. Their trace distance is `|sin delta|≈0.00406650`. The optimal one-copy binary error is about `0.49796675`; the corresponding accessible mutual information is approximately `1.19×10^−5` bits per ideal probe. The Holevo upper bound is about `7.99×10^−5` bits. These standard quantum-information calculations use the supplied Born framework and ideal coherent arms; no experimental information rate was measured. Finite selection and imperfect visibility require the actual instrument probabilities instead.

| Quantity | What was supplied or calculated | What remains unknown |
|---|---|---|
| Preparation acceptance | A specific admissible stationary profile was supplied and solved | Probability and cost of physically preparing it |
| Matter energy/charge | `E=3639.23 v`, mean/classical `Q≈39943.85` | A loading protocol with all apparatus losses |
| Probe resource at `k=5` | One reporter quantum has energy `5.1 gv=0.51 v` | Source, beam splitter, recombination, detection and reset costs |
| Information recoverability | Ideal phase-dependent probabilities and binary distinguishability | Full target–probe quantum channel and detector implementation |
| Delivery | Finite packet reaches the selected outgoing region with about 99% probability in the stated reduced instrument | General source acceptance, finite-mass losses, other preparations |
| Retention | Stationary leading profile plus local fluctuation checks | Nonlinear/quantum lifetime and stored-message retention |

One `k=5` quantum carries about `1.4×10^−4` of the lump's energy. This suggests a test-probe hierarchy but is not a rigorous bound on back-reaction or inelasticity. No reservoir sustains the stationary leading solution. All preparation and readout interventions remain physical resources, with missing costs explicitly uncomputed.

## 8. Reproduction and claim boundary

Run `python code/qball.py`, then `python code/reporter.py`, then `python code/verify.py` from the package root. See `code/README.md` for environment, execution logs, failed starts, source hashes and output meanings. Numerical tolerances establish internal convergence and consistency; they are not experimental error bars. The preparation selection and all failed starts are retained in the archive.

**Achieved:** a populated leading classical regime, numerical support for local fixed-charge energetic stability, complete-stress weak gravity, reporter phase, and a finite-packet ideal readout calculation within the same action.

**Not established:** a universal microscopic theory, a quantum soliton lifetime, global fission stability, autonomous production, implemented apparatus, or agreement with measurements of nature.
