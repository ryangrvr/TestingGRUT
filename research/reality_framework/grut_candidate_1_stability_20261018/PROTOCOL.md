# Candidate 1: bounded nonlinear and fragmentation stability protocol

**18 October 2026 · physical action and couplings unchanged**

Continue the Candidate 1 preparation at `a=1,b=.1`, `Qbar=399.4384624223012`, `Ebar=363.9229985516867`, with its archived radial profile. X charge is supplied, global and not electric charge. This follow-up tests the leading nongravitating scalar mechanism. It does not simulate the reporter, apparatus, Standard Model or metric. Small benchmark gravity is not a substitute for calculating those sectors.

## Questions and finite limits

1. Do found lower-charge stationary branches permit a lower-energy fragmentation of the benchmark?
2. Does a smooth nonradial perturbation at the same charge retain a localized interacting core over the stated finite interval?
3. Are any changes resolved against independent time-step, spatial and domain effects?

The Calculator freezes its branch continuation and alternative starts separately in `branch_protocol.json` before the branch scan. Failed starts, turning points and incomplete charge ranges must remain visible. A sampled positive branch margin is not a proof of global minimization.

Root's nonlinear matrix is frozen below before the first physical trajectory. Local FFT/stencil timing benchmarks were performed to select finite resources; they used random arrays, not physical states. Calculations are finite local CPU work, no installs, paid execution or indefinite workers. The intended nonlinear matrix is approximately ten minutes of local execution, with hardware-dependent timing recorded rather than guaranteed.

## Exact continuum equations and Hamiltonian

Let `A=sqrt(2)X/v`, `S=phi/v`, `tau=g_X v t`, and `y=g_X v x`. Then

\[
\ddot A=\nabla^2A-(S^2+.1|A|^2)A,\qquad
\ddot S=\nabla^2S-[(S^2-1)+|A|^2]S.
\]

The charge is `Qbar=−integral Im(A* A_dot)`. The energy is

\[
\bar E=\int d^3y\left[\tfrac12(|\dot A|^2+|\nabla A|^2+\dot S^2+|\nabla S|^2)
+\tfrac14(S^2-1)^2+\tfrac12S^2|A|^2+.025|A|^4\right].
\]

Use Fourier collocation in a periodic cube and Strang splitting into the exactly evolved linear Klein–Gordon Hamiltonian and a local nonlinear potential kick. For `q=(Re A,Im A,s=S−1)`, the linear masses squared are `(1,1,2)`. The remainder potential is `Ubar−|A|²/2−s²`. Both subflows preserve the actual discrete global U(1) charge. Energy is approximately conserved by splitting and evaluated using the matching spectral `k²` gradient energy, not a different derivative stencil.

This is a finite-dimensional Hamiltonian approximation to the continuum action. No exactly local finite-grid light cone is assumed, and no dealiasing filter or artificial damping is applied. Spatial refinement tests the collocation/truncation effects.

## Preparation and perturbations

Evaluate the archived F,S profile on the actual Cartesian grid, with its vacuum continuation beyond the radial table. Set `R=4` and

\[
q_2=\frac{2z^2-x^2-y^2}{R^2}e^{-r^2/(2R^2)},\quad
\chi=\frac{xy}{R^2}e^{-r^2/(2R^2)}.
\]

Both are smooth at the origin. For amplitude epsilon,

\[
A_0=F(r)(1+\epsilon q_2)e^{i\epsilon\chi},\qquad
S_0=S(r)+\tfrac\epsilon2(1-S(r))q_2,
\]

\[
\dot A_0=-i\Omega_{grid}A_0,\qquad \dot S_0=0,
\quad \Omega_{grid}=\bar Q_{target}/\sum|A_0|^2h^3.
\]

This normalization fixes the same charge on every grid. The perturbation changes energy; the actual initial energy and excess must be printed. Initial angular gradients and scalar deformation are paid for. It is not a cost-free information perturbation.

## Nonlinear run matrix

All runs end at `tau=40` and record diagnostics every unit time.

| Name | Cube side | N per direction | dt | epsilon | Purpose |
|---|---:|---:|---:|---:|---|
| control | 64 | 96 | .1 | 0 | Unperturbed continuum-profile control |
| small | 64 | 96 | .1 | .05 | Smooth weak nonradial perturbation |
| stress | 64 | 96 | .1 | .20 | Larger fixed-charge perturbation |
| time | 64 | 96 | .05 | .20 | Independent time refinement |
| coarse | 64 | 64 | .1 | .20 | Spatial comparison |
| space | 64 | 128 | .1 | .20 | Independent spatial refinement |
| domain | 256/3 | 128 | .1 | .20 | Same h as stress, larger periodic box |

Compare only matched physical times. The finest spatial/time run is not silently substituted for every control. This duration is about five carrier rotation periods and a few internal dynamical periods, not a demonstrated long lifetime. Periodic radiation can return; initial tail/boundary mismatch can also matter. Print initial boundary amplitudes and compare stress versus domain. A continuum propagation estimate alone is insufficient to certify numerical isolation.

## Diagnostics, rejection and reporting

Record global charge and the matching Hamiltonian; charge/energy inside `r<=8` and `r<=12`; charge radius, central S, core peak amplitude, quadrupole norm, phase-aligned amplitude/profile error and central slices. Track a core directly: conserving global totals alone cannot show localization.

- Numerical validity target: global relative charge drift below `1e−9`, relative energy drift below `2e−3`. If these fail, report an unresolved numerical regime rather than a physical instability or success.
- Operating diagnostic: core charge within `r<=8` remains above 90% of the initial total charge for every sampled time. Failure is a failure of that particular finite-time operating diagnostic, not a no-go theorem.
- Report resolved growth, oscillation or dispersal relative to the unperturbed control. Report maximum differences of matched observables under dt, h and box refinements; no rigorous continuum error bound is presumed.
- Do not infer all-perturbation stability, every fission partition, quantum lifetime, full gravitational stability, autonomous formation or empirical agreement.

If a solver implementation or protocol must be corrected, preserve the first result and explicitly version the correction. Do not tune a physical coupling, erase a failure, or change the perturbation after seeing an unfavorable physical result.
