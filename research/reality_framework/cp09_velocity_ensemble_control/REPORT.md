# CP-09 bounded clock: an independently calibrated classical velocity-ensemble control

2026-10-11. Independent analysis of the **corrected unit-area** positive quadratic clock at [986394a](https://github.com/ryangrvr/TestingGRUT/blob/986394af97f64111af62dcd10ed2254ad6779f9e/research/reality_framework/cp09_clock_area_control/REPORT.md). Separate branch forked from 986394a. Frozen CP-09, CP-09 moving/bounded clock simulator source, Stage-3 Card 3 and prior constructions untouched. **This is not CP-10, a new clock law or a derivation of spacetime.**

## Why another comparator is needed

The corrected bounded-clock response at m=320, spatial-grid 2048 and time step .00125 is 0.0536830950652. At the same finite window width 0.15 and integrated couplings epsilon=zeta=0.001, the exactly translating clock gives 0.0523491926013; the instantaneous CP-09 Kubo target is 0.0525885135387.

The remaining discrepancy against the finite-width ideal mixes a quantum clock's momentum spread, phase coherence, recoil, wave-packet dispersion, spatial discretization, hard-edged interaction windows, periodic-circle effects and residual finite-pulse corrections. **A classical velocity mixture is a controlled way to account for one ingredient without claiming it explains all the others.**

## Independently derived classical transport ensemble

The originally prepared quadratic clock has an unchirped Gaussian wavepacket centered in position at q=-0.7 with sigma_q=0.08, mean momentum p_0=m, and sigma_p=1/(2 sigma_q)=6.25, at hbar=1. If one temporarily *replaces the quantum clock motion and clock-matter back-action by independent deterministic classical trajectories at sampled velocity v=p/m>0* with the same Gaussian initial momentum distribution, a trajectory traverses the source and detector spatial windows (width delta=0.15, separation tau=0.4) in durations
\[
t_{\rm pulse}=\delta/v,\qquad t_{\rm gap}=(\tau-\delta)/v.
\]
The physical spatially normalized pulse potential is epsilon f_delta(q), with f_delta=1/delta on each region. **At velocity v it therefore has integrated temporal strength epsilon/v, not epsilon**. Both source and detector get the same factor 1/v. The correct trajectory simulation does *not* simply rescale tau and leave the pulse strengths unchanged.

For each classical trajectory use the full supplied finite CP-09 system+pointer evolution
\[
U_s^\pm(v)=\exp[-i(H_s\otimes I_r\mp(\epsilon/\delta)B\otimes I_r)\delta/v],\\
U_0(v)=\exp[-i(H_s\otimes I_r)(\tau-\delta)/v],\\
U_d(v)=\exp[-i(H_s\otimes I_r-(\zeta/\delta)A\otimes Y_r)\delta/v].
\]
Compute the same pointer-X source-sign-difference estimator as the author, then average this trajectory response over the Gaussian velocity distribution \(v\sim {\cal N}(1,\sigma_v^2)\) with \(\sigma_v=1/(2m\sigma_q)=0.01953125\) for m=320. This is a **classical surrogate**, not the actual quantum calculation or a new physical clock proposal. The exponentially tiny negative-velocity Gaussian tail is numerically irrelevant for these parameters.

An independent weak-pulse calculation reconstructs the stationary CP-09 Kubo response directly from its exact twelve-eigenstate spectrum:
\[
R_{AB}(s)=-2\operatorname{Im}\sum_n
\langle 0|A|n\rangle\langle n|B|0\rangle e^{-i(E_n-E_0)s}.
\]
Integrating the response over two finite square windows yields the **linear-response** expression
\[
\boxed{R_\delta(v)=\frac{-2}{v^2}\operatorname{Im}
\sum_n\langle0|A|n\rangle\langle n|B|0\rangle
e^{-i\omega_n\tau/v}
\operatorname{sinc}^2\!\left(\frac{\omega_n\delta}{2v}\right).}
\]
The factor \(v^{-2}\) is essential. It is a clock-speed-dependent *pulse-area* effect. The two independently implemented methods—full 24x24 source/pointer exponentials, and spectral linear response—agree at each sample velocity to better than \(7\times10^{-8}\) for epsilon=zeta=.001, the expected weak-strength difference. Averaging with Gauss–Hermite quadrature at 8, 12 and 16 nodes agrees to better than 1e-10.

## Numerical benchmark at the SAME controls

| Quantity | Response or shift |
|---|---:|
| Instantaneous ideal Kubo \(R(\tau)\) | 0.052588513539 |
| Translation clock, spatial square width delta=.15, v=1 | 0.052349192601 |
| Same finite-window model with classical Gaussian velocities, m=320 | **0.052489125240** |
| Independent first-order spectral double-window Gaussian average | 0.052489174957 |
| Corrected positive quadratic clock, imported from 986394a | **0.053683095065** |

With the finite-strength 24x24 trajectory benchmark, the observed difference has the **exact arithmetical decomposition**

\[
\begin{aligned}
R_{\rm bounded}-R_{\rm instantaneous}
&=(R_{v=1,\delta}-R_{\rm instantaneous})\\
&\quad+(\overline R_{\rm class.\,velocity}-R_{v=1,\delta})\\
&\quad+(R_{\rm bounded}-\overline R_{\rm class.\,velocity})\\
&=\boxed{-0.000239320937+0.000139932638+0.001193969825}\\
&=\boxed{+0.001094581526}.
\end{aligned}
\]

Thus the classical velocity distribution accounts for about **10.5%** of the m=320 bounded-minus-translation finite-width gap, with 89.5% left **unattributed** by this classical surrogate. This 89.5% is **not evidence for a uniquely quantum mechanism**. In particular it also contains possible finite-grid/circle/sampling errors, coherent wavepacket effects, momentum exchange and the possibility that the packet has not fully cleared the interaction windows in all tails.

A separate small-spread calculation expands the trajectory response around v=1:
\[
\overline R_{\rm class.\,velocity}
=R(1)+\frac{\sigma_v^2}{2}R''(1)+O(\sigma_v^4).
\]
Numerically \(R''(1)\simeq0.73221\), giving a leading shift of approximately 0.000139659, in good agreement with the full Gaussian average +0.000139933. This makes the control sensitive to the velocity distribution's *variance* and its predicted \(m^{-2}\) scaling, rather than a fitted offset.

## What this does and does not certify

This calculation uses an **independent, unmodified N=2 CP-09 Hamiltonian** and retains exactly the source/detector operators of the frozen role-defined benchmark. It does **not** rerun the full corrected quantum clock lattice dynamics, prove its continuum-grid convergence, separate recoil from coherent propagation, or derive a physically local implementation of the global relation-dependent probes. The positive quadratic clock's actual finite-width response must not be represented as this positive classical mixture: they are distinct experiments.

The strongest next control is now well defined: retain the corrected discrete pulse-area normalization, use multiple masses and window widths with separate spatial-grid and timestep convergence, and compare the full quantum simulation against (a) the exact fixed-v translation baseline, (b) the **predeclared classical Gaussian velocity-mixture baseline**, and (c) source/pointer clock-arrival records. Smooth normalized window profiles would reduce grid aliasing from sharp indicators. Report the source/detector arrival distributions, clock momentum distributions before and after the interactions, and a complete Hamiltonian energy ledger. That would begin to separate dispersion, clock recoil, coherent quantum effects and numerical artifacts without changing CP-09.

The foundational gap is larger: the detector/source coupling still reads a supplied global causal-role operator and the low-energy causal ensemble is not demonstrated to favor physical spacetime. Clock refinement cannot earn locality, Einstein gravity or a new empirical prediction by itself.

## Reproduce

\`python research/reality_framework/cp09_velocity_ensemble_control/check.py\` (NumPy/SciPy). The script independently builds the CP-09 N=2 Hamiltonian, cross-checks the full finite pulse with spectral convolution, integrates the declared classical free-packet momentum distribution by Gaussian quadrature and compares the reported 986394a corrected number. It prints \`CP09_VELOCITY_ENSEMBLE_CONTROL_PASS\` when its internal assertions pass. The published corrected bounded-clock output is an **explicit source number**, not independently recomputed by this script.
