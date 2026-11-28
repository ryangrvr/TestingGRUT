# C3 original reporter and empirical feasibility

**Verdict C.** The unchanged Candidate 1 action gives a quantitative scattering observable for a prepared target. A calibrated source, target population, and ordinary-matter readout for that observable remain unsupported at this benchmark. The strongest conventional comparator has exactly the same scalar elastic predictions. The calculation advances the prediction and identifies its physical access requirements. [Original action](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md) · [Continuum certificate](../grut_c2_continuum_stability_20261127/REPORT.md) · [Raw results](RESULTS.json)

## Static elastic observable

On the leading flat scalar saddle with $D=0$, the linear reporter equation decouples. With $k=p/(g_Xv)$, physical on-shell frequency $\sqrt{k^2+1.01}\,g_Xv$ and radial measure $d\rho$, its vacuum-relative equation is
\[
u_\ell''+\left[k^2-\frac{\ell(\ell+1)}{\rho^2}-0.01(S^2-1)\right]u_\ell=0.
\]
Regularity fixes the inner solution. The asymptotic Riccati–Bessel form gives the real elastic phase $\delta_\ell$ and $S_\ell=e^{2i\delta_\ell}$. Incoming and elastically outgoing Klein–Gordon wave fluxes use the same $p/E_D$, so the partial cross section has the ordinary flux normalization
\[
\sigma_\ell=\frac{4\pi(2\ell+1)}{p^2}\sin^2\delta_\ell,\qquad
f(\vartheta)=p^{-1}\sum_\ell (2\ell+1)e^{i\delta_\ell}\sin\delta_\ell\,P_\ell(\cos\vartheta).
\]
The calculation requires a weak probe, an approximately stationary target, and the leading EFT valid at the specified momenta. The partial and total quantities below are for the short-range scalar potential in this flat, coherent elastic limit. Long-range gravitational phase and scattering, resolved recoil, inelastic channels and full quantum corrections require separate treatment.

| $k$ | Original reference $\delta_0$ (rad) | $Q=39,944$ predictor $\delta_0$ (rad) | $\sigma_0(g_Xv)^2$ | Retained $\sigma_{\rm elastic}(g_Xv)^2$ |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 0.006776510978 | 0.006776520524 | $6.41172\times10^{-5}$ | $4.88330\times10^{-3}$ |
| 5 | 0.004066511575 | 0.004066517304 | $8.31214\times10^{-6}$ | $1.75787\times10^{-3}$ |
| 7 | 0.002904769535 | 0.002904773628 | $2.16390\times10^{-6}$ | $8.96859\times10^{-4}$ |

The original $k=5$ phase is reproduced to $4.96\times10^{-12}$ rad against its displayed value. The target phase is checked by an independent variable-phase integration and a regular Numerov wave matched to free Bessel functions. The maximum S-wave method difference across the three momenta is $3.78\times10^{-10}$ rad. At $k=5$ the change in the retained total between the two radius, mesh and angular choices is $5.03\times10^{-13}$ in $(g_Xv)^{-2}$ units. The $\ell=65,\ldots,80$ contribution is $1.03\times10^{-13}$ in those units. These are empirical refinements. No full angular-tail enclosure is inferred from them. Some far-angle differential values are close to cancellation zeros and are sensitive to the angular cutoff. The full retained angular data and changes are in the archive. [Results](RESULTS.json) · [Reproduction](REPRODUCTION.md)

The exact C2 profile contributes a controlled S-wave sensitivity. If $r_{\rm C2}$ is its graph-norm error radius, then $\|\delta(S)\|_\infty\le12r_{\rm C2}$ on $[0,1]$ and $\|\delta(\rho S)\|_2\le r_{\rm C2}$ on the half-line. Cauchy–Schwarz with $\int_1^\infty\rho^{-2}d\rho=1$ gives $\|\delta(S)\|_1\le13r_{\rm C2}$. Since both profiles have $|S|\le2$, the change of scalar potential has $L^1$ norm at most $0.52r_{\rm C2}$. The variable-phase equation $\delta'=-(U/k)\sin^2(k\rho+\delta)$ and Grönwall's inequality give an exact-profile phase allowance $1.52\times10^{-8}$ rad at $k=5$. That allowance concerns the continuum classical profiles only. The numerical integration comparison is separate.

## Quantum and recoil effects

At fixed measured momentum the vacuum reporter pole shift changes the frequency used to set $p$. At fixed probe energy, the archived strict pole changes the inferred mass by $-8.39\times10^{-6}$ fractionally and the leading S-wave phase by approximately $-1.38\times10^{-9}$ rad. The two-body center-of-mass momentum correction alone shifts the phase by approximately $+5.70\times10^{-7}$ rad. These are kinematic estimates. A scale change to $5g_Xv$ gives a scalar leading-log $|X|^2 D^2$ coefficient of approximately $2.04\times10^{-8}$ from the prescribed zero-matched coefficient. Cancellation of scale dependence requires the remaining graphs. At this probe energy, excited composites or free scalar radiation can be kinematically relevant if the corresponding interactions and thresholds allow them.

A full soliton one-loop reporter response needs the renormalized inhomogeneous self-energy, correction to the prepared profile and consistent elastic/inelastic matching. It is unresolved in this packet. The vacuum pole effect is small for the current phase estimate. The identified absence of calibrated physical access persists independently of the numerical size of the inhomogeneous amplitude. [Archived pole](../grut_candidate_1_fluctuations_20261019/POLE_THRESHOLDS.md) · [Running conventions](../grut_candidate_1_quantum_20261018/RG_DERIVATION.md)

## Preparation and physical scales

Use the unchanged illustrative $v/M_P=10^{-5}$ and $g_X=0.1$, with the reduced Planck mass. CODATA conversion gives $g_Xv=2.4353\times10^{12}$ GeV. The required target carries 39,944 units of net $X$ charge. [CODATA constants](https://physics.nist.gov/cuu/Constants/)

| Quantity | Required or predicted scale |
| --- | ---: |
| Target rest energy and mass | $1.4201\times10^7$ J, $1.5800\times10^{-10}$ kg |
| Radius enclosing 90% of classical charge | $\rho_{90}\simeq4.716$, $3.8210\times10^{-28}$ m |
| Rest energy of the corresponding cold free carrier population | $1.5587\times10^7$ J |
| Energy to remove if assembling those cold carriers into the specified rest state | $1.3866\times10^6$ J |
| One $k=5$ reporter energy | $1.990\times10^3$ J per quantum |
| Reporter pair rest-energy threshold | $7.842\times10^2$ J |

The exact $U(1)_X$ imposes charge accounting. Ordinary initial matter carries zero $X$ charge. A prepared positive-charge target from such an initial state would require a compensating negative charge elsewhere, or a pre-existing charged sector. Combining cold carriers also requires a dissipation route. The rest-state capture of one carrier by the neighboring-charge target releases $0.200069\,g_Xv$, below the free binding-scalar mass $\sqrt 2\,g_Xv$. A single on-shell binding scalar cannot carry that release in the cold capture channel. Hot collisions, collective energy redistribution, other processes and gravitational radiation require rates and an actual environment. The action and certified saddle provide no characterized preparation yield in a lab or astrophysical setting. The applicable survival results are classical linear stability at the specified charges and the two C1 energetic comparisons. A survival probability or lifetime for a prepared population remains open. [C1 comparisons](../grut_c1_quantum_channels_20261127/REPORT.md) · [C2 certificate](../grut_c2_continuum_stability_20261127/REPORT.md)

## Source, detector and attainable event rate

The $\phi^2D^2$ interaction couples the reporter to the hidden scalar well. Mixed portals to the Standard Model have zero renormalized coefficients at the specified matching scale. Allowed counterterms and running remain part of the EFT. Both sectors couple to the Einstein metric. Gravitation can therefore mediate transitions, including pair production with sufficient energy. It supplies no calibrated beam source or characterized readout at this benchmark. Exact $D$ parity keeps vacuum production from ordinary matter even in $D$. The original ideal interferometer supplies arm, flux and measurement data externally. A physical realization needs reporter generation, energy and momentum calibration, preparation and monitoring of targets, angular acceptance or mode-resolved recombination, and a registered energy or timing signal in ordinary matter with characterized efficiency and background. No same-action implementation of this chain was demonstrated. [Original action and apparatus assumptions](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md)

At $k=5$ the S-wave cross section is $5.46\times10^{-62}\,\mathrm{m}^2$. The retained coherent scalar elastic total is $1.15\times10^{-59}\,\mathrm{m}^2$. Given reporter fluence $\mathcal F$, $N_T$ independent known targets, a calibrated acceptance cross section $\sigma_{\rm acc}$, and efficiency $\epsilon$, the dilute expected signal is $N_{\rm sig}=\mathcal F N_T\sigma_{\rm acc}\epsilon$. An ideal one-target experiment accepting the whole retained elastic total with unit efficiency needs about $8.66\times10^{58}\,\mathrm{m}^{-2}$ for mean one event. Using only the S-wave term gives $1.83\times10^{61}\,\mathrm{m}^{-2}$. These are conditional figures using supplied fluence and target. Even saturating the elastic S-wave unitarity cross section at this momentum gives only $3.30\times10^{-57}\,\mathrm{m}^2$. A measured background, physical flux, target population, exposure and phase-error model are unavailable. An empirical significance or rate cannot be assigned from the static cross section.

## Comparator and falsification

A conventional FLS-type scalar bound state with the same quartic term, matched reporter and the same static target profile has the identical mode equation and asymptotic flux. More generally any static scalar well with the same vacuum-relative mass profile reproduces all leading elastic phase shifts in the declared preparation. The predicted event distributions agree whenever the probe, target profile and detector calibration agree. The likelihood ratio is one for this matched comparison. No finite exposure separates identical distributions.

A more restrictive zero-range isotropic model matched only to the $k=5$ S-wave component has a total smaller by a factor of about 211.5 and differential cross section $\sigma_0/(4\pi)$. The finite-size comparator is sharply forward peaked. This is a conditional finite-size discriminator for a declared isotropic model. The same extended-well mechanism supplies it in conventional physics. At distances far beyond the soliton radius the leading gravitational monopole is that of a point mass with the same total energy. Other given matter models require their own profile, channel and scale matching. [Conventional comparisons](NOVELTY_AUDIT.md)

A future test of this benchmark would independently calibrate the mass scale, probe, charge, target profile class and acceptance. It would reject the predicted scattering law if a replicated angular or phase distribution disagreed beyond combined experimental and theoretical uncertainty in the matched elastic regime. The static results alone cannot specify an experiment-level rejection threshold. This packet establishes a conditional scalar scattering prediction and its missing physical links. It supplies no relic abundance or independent fundamental-law deviation. [External RAI review](RAI_REVIEW_PACKET.md)
