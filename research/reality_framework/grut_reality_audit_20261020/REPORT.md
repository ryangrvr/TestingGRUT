# GRUT reality audit: which hypotheses work, and what is still missing?

20 October 2026, America/Chicago · Separate full-field numerical tests and empirical-access audit · Author research; independent scientific review pending

**Verdict:** the published classical bulk-wave and mixed-branch growth predictions survive a separate full-field implementation. The record does not establish that the added fields exist in nature, that finite objects are robust in unrestricted 3D, or that GRUT explains the observed particles, cosmological abundance or consciousness. The displayed leading benchmark is a hidden scalar sector coupled to ordinary matter through the adopted metric; its ideal reporter is not yet an ordinary laboratory detector.

Claude's current tests are unpublished according to the owner. Their sources, preparations, output and error controls are unavailable here. This report neither imports their conclusions nor assigns them a pass/fail result. It uses the already published record and the new calculations preserved in this directory.

## 1. The new test is an evolution, not another matrix comparison

The [theory packet](../grut_theory_synthesis_20261020/THEORY.md) at [0cbb364](https://github.com/ryangrvr/TestingGRUT/commit/0cbb3644f3466c6dade3a00d919c62a97b24f087) predicts the coupled scalar bulk spectrum and a mixed-branch unstable band. Its numerical spectrum checks diagonalized the linear generator. This audit instead evolves the **full nonlinear equations for u,w,S and their momenta**, using newly written Python code, periodic Fourier spatial differentiation and explicit fourth-order Runge–Kutta time steps. The earlier collision solver used cylindrical finite volumes and Verlet.

No rotating-frame force, maintained background, imposed rotation frequency, tangent equation or damping acts in these evolutions. Only the initial coherent preparation is supplied. A frame derived from the evolving spatial mean is used afterward to measure small fluctuations.

The fields depend on one spatial coordinate: these are plane-dependent solutions of the 3D field equations around uniform bulk. Energy and charge are measured per unit transverse area. They are not finite droplets, collision trajectories, cosmological histories or unrestricted-3D stability experiments.

All ten prescribed evolutions completed. Before launching them, [LOCAL_FREEZE.json](LOCAL_FREEZE.json) recorded hashes of the solver, run list and supervisor. The frequency/growth targets, tolerances, estimator description and fitting windows were fixed in [protocol.json](protocol.json). That was a local pre-run freeze, not a public preregistration or a blinded empirical experiment. The frequency-estimation script was implemented after execution according to that specification. No physical coefficient was fitted to the evolution data.

The estimator starts with the largest separated nonzero Fourier peaks, then fits sinusoidal frequencies with free amplitudes and phases. Theoretical frequencies do not initialize the fit. Its number of branches and search interval are theory-informed. Growth uses the previously specified log-norm windows and published matrix predictions.

## 2. Outcomes

The dense background has a=1,b=.1, F=10^(1/4), S=0 and Omega=.1^(1/4). A charged-amplitude perturbation excites both charged branches; a mediator perturbation tests the third branch.

| Full-field test | Previously predicted frequency/rate | Measured frequency/rate |
|---|---:|---:|
| Dense charged, k=.2, lower branch | .11703821171 | .11703821153 |
| Dense charged, k=.2, upper branch | 1.40130962071 | 1.40130961896 |
| Dense charged, k=.8, lower branch | .53027194345 | .53027194139 |
| Dense charged, k=.8, upper branch | 1.70181616578 | 1.70181616196 |
| Dense mediator, k=.2 | 1.48400729788 | 1.48400728742 |
| Dense mediator, k=.8 | 1.67400049587 | 1.67400047760 |
| Mixed coherent background, k=.2, growing mode | .07581066444 | .07580819919 |
| Mixed coherent background, k=.4, growing mode | .11212539255 | .11211726772 |

The largest measured dense-frequency relative discrepancy is 1.10e-8. Mixed-branch growth discrepancies are approximately 0.0033% and 0.0073%, within the declared 2% limit. At k=.4, increasing the seed tenfold changes the fitted growth rate by 8.59e-8 relative. The larger-seed run ends with a mode norm about 1.27% of mean F, so these runs remain mainly in the small-perturbation regime. They do not establish the nonlinear fragmentation endpoint.

The k=.8 mixed-background control is above the predicted unstable band. Its mode norm gain stays below 1.381 during the tested interval. This is a reported finite-time observation; no post-run universal stability threshold is invented for it.

Halving the time step changes the two k=.2 charged frequencies by at most 1.66e-9. Doubling the spatial grid changes them by at most 4.87e-13. The largest sampled relative energy drift over all ten cases is 3.09e-9; the largest charge drift is 3.32e-9. RK4 does not exactly preserve charge, so its measured drift is retained rather than described as roundoff conservation.

All **22 source, numerical and prediction checks** pass: three source identities, ten numerical-consistency checks and nine frequency/growth comparisons. These are not 22 independent physical discoveries or experimental confirmations. Residuals are numerical comparisons in highly controlled synthetic preparations, not laboratory uncertainties. The same author constructed both implementations; this is separate-code support, not independent human review.

![Full-field wave and growth comparisons](AUDIT_FIGURE.svg)

The sound-speed statement is a long-wavelength limit. For example, the finite-k=.2 lower branch has nu/k about .58519; its difference from 1/sqrt(3) is the predicted dispersion correction. The new runs support the full dispersion formula, not an assertion that every finite-wavelength wave travels exactly at c/sqrt(3).

## 3. What the hypotheses actually say

| Hypothesis or larger interpretation | Audited standing |
|---|---|
| Scalar feedback permits binding | Supported by the restricted classical energy bound, existing stationary configuration and conditional assembly calculations. This instantiates a known FLS/Q-ball mechanism. |
| Dense bulk and the specified coherent mixed branch have the predicted waves/growth | Supported by the new full-field tests above, with small perturbations and stated geometry. |
| The condition ab<1 forces all preparations to organize | False as a universal claim. The exact classical vacuum X=D=0, phi=v remains a stationary solution. The binding inequality makes states energetically possible; it does not populate them. |
| Every diffuse charged cloud forms a durable object | Not established. The previously tested width-8 cloud is marginal; preparation matters. |
| The objects survive general 3D disturbances | Unresolved. The collision m=1 response was resolution-sensitive. Uniform bulk modes do not settle finite interfaces or trajectories. |
| A collision leaves an acquired charge distinction | Supported over the earlier registered finite windows with an ideal observable. No implemented readout/controller or human-scale retention is established. |
| The model derives quantum mechanics, the Standard Model or gravity | Those structures are adopted. Their successful reproduction is not a new empirical distinction. |
| These fields describe observed matter or dark matter | No empirical calibration, abundance prediction or new-sector detection is present in the inspected record. |
| The mechanism explains consciousness | Not established by the field dynamics, records or this audit. |

There is also a useful limit on the apparent universality of the sound speed. The EFT admits a bare X mass operator, whose finite matched value was set to zero at the benchmark scale. If its dimensionless coefficient is c0, the uniform S=0 branch has

\[
\Omega^2=c_0+bF^2,\quad n=\Omega(\Omega^2-c_0)/b,\quad
c_s^2=\frac{\Omega^2-c_0}{3\Omega^2-c_0}.
\]

The value 1/3 therefore follows from the selected leading quartic core with c0=0. It is not a fundamental constant derived without matching assumptions. This sensitivity formula describes an allowed changed matching parameter; c0 was not changed in the ten physical evolutions.

## 4. Put the benchmark into physical units

The original theory specified g_X=.1 and the illustrative v/M_P=1e-5, with **reduced** Planck mass M_P. Using NIST's 2022 CODATA G and SI defining constants gives M_P approximately 2.4353e18 GeV. [physical_scales.py](physical_scales.py) performs the conversion from the archived stationary profile; [PHYSICAL_SCALES.json](PHYSICAL_SCALES.json) retains inputs, formulas, outputs and uncertainty caveats.

| Quantity | Original illustrative benchmark |
|---|---:|
| v | 2.44e13 GeV |
| Free carrier mass m_X | 2.44e12 GeV |
| Existing soliton mass | 8.86e16 GeV, or 1.58e-10 kg |
| Radius containing 90% of its charge | 4.14e-28 m |
| Time represented by tau=240 | 6.49e-35 s |
| Duration represented by tau=160..240 | 2.16e-35 s |

These are consequences of an assumed scale, not measured properties. A finite simulated interval is not a lifetime estimate. The numbers do not establish human-scale memory, and do not disprove longer persistence that has not been calculated. The inherited compactness, about 5.67e-10, is recovered by the SI conversion.

A scale-independent description is useful: tau240 spans about 47.0 light-crossing times of the archived initial R90, and the late readout interval spans about 15.7. Rescaling tau240 to one second would simultaneously give that reference object a radius about 6.38e6 meters. This relates the present one-scale preparation and observation window; a composite protected memory could acquire additional timescales, but that construction and its lifetime have not been calculated.

Changing the free scale changes the interpretation dramatically. For illustration, choosing m_X=1 GeV at the same dimensionless scalar parameters would give a soliton mass about 36,392 GeV, radius about 1.01e-15 m, and tau=240 about 1.58e-22 s. No observation selects that alternative here. The flat dimensionless simulations cannot determine v.

At the specified scalar ratios and archived preparation, one can eliminate v:

\[
\frac{m_\phi}{m_X}=\sqrt2,\quad
\frac{m_D}{m_X}=\sqrt{1.01},\quad
\frac{M_{\rm sol}}{m_X}=36392.3,\quad
\frac{R_{90}m_X}{\hbar c}=5.10343.
\]

The last expression uses m_X as an energy in natural-unit notation. These are conditional correlated targets after fixing the couplings and state, not measured relations or parameter-free universal constants. A pole-mass/scattering calibration and quantum corrections are needed for a precision physical prediction.

## 5. The missing detector connection is structural

In the displayed leading action, X,phi,D are Standard Model singlets. The benchmark's direct scalar-Higgs/mixed portals are set to zero at matching. The leading ordinary sector and hidden sector communicate through their common metric. In particular, a nonzero scattering phase of D does **not** by itself predict a phase for an ordinary laser or an electron detector. D is another proposed new particle, not a synonym for light.

The EFT admits additional portal, higher-order and gravitational-loop effects. Their coefficients and observational response would have to be calculated or calibrated at the stated order. A direct laboratory portal cannot be silently added and called a tested prediction of the unchanged benchmark. A known ordinary apparatus can measure a metric response in principle, but the archived weak-gravity calculation has not supplied an experimentally achievable source or resolution at the scales above.

An exterior field agreeing with Einstein gravity at a given source mass does not uniquely identify the internal scalar composition. A distinguishing observation must constrain the additional scalar mechanism or its correlated population behavior.

The prior head-on collision cases also do not provide a halo-velocity transfer cross section. The scale file reports a geometric contact-area proxy only, explicitly labeled as such. It is not compared as a computed cross section against observational exclusion curves.

## 6. A concrete route toward nature

The present action most directly motivates a **candidate hidden matter sector coupled to GR**. That is a research interpretation, not evidence that the modeled objects are dark matter. Its next empirical specification must contain all of the following in one compatible model:

1. **Population production:** an initial hidden-sector state, temperature/reheating history, charge asymmetry and formation mechanism. The current Gaussian and homogeneous preparations do not predict a cosmic number density. A thermal transition requires the actual finite-temperature potential; a first-order transition cannot be borrowed from a different potential.
2. **Abundance and lifetime:** predict the resulting charge/mass distribution, survival and mergers. For an entirely one-sign population carrying charge density n_Q, rho_sol=(M/Q)n_Q. Adjusting an unconstrained n_Q to match an observed density is a calibration, not an abundance prediction. Mixed-sign populations require a separate treatment.
3. **Cosmological stress and structure:** a population of nonrelativistically moving stable droplets can have small center-of-mass pressure even when internal dense matter has sound speed c/sqrt(3). Its distribution, velocity history and interactions must be evolved; internal sound is a different observable from the population's effective pressure.
4. **Observable channel:** select the gravitational/population response allowed by the current matching, or explicitly version an extended action with a measurable portal. Then calculate source, propagation and detector response with backgrounds and noise.
5. **A held-out comparison:** fix parameters and preparation from separate information, predict a named data product with uncertainty and compare competing models. A single fitted number or a successful scalar simulation is insufficient.

The 2026 PDG dark-matter review provides a current map of abundance, structure, lifetime and interaction requirements; it is not evidence for this scalar model [R2]. A recent primary FLS formation study illustrates the amount of production modeling required, using a different thermal potential and fluid/dissipative dynamics [R3]. Its abundance or stability results cannot be imported into the present quartic Hamiltonian benchmark. Selected source sections were inspected; neither the PDG evidence base nor that paper's solver was independently reproduced here.

## 7. The next deciding work

The existing symmetry-generated translation/tilt controls and unrestricted-3D perturbation test remain the decisive finite-object question. The new uniform-bulk tests do not close it. After finite-object behavior is resolved, choose a stated cosmological/experimental regime and calculate a population plus observable without changing parameters separately for each target.

The project also needs a physical readout/control construction before describing the acquired charge distinction as an observer. Functional causal use of a record and subjective experience remain separate claims.

For coordination, [START_HERE.md](START_HERE.md) is a short durable entry point with exact source boundaries, what ran, what is unknown about Claude's current work, and the next required handoff contents. A new chat should read it before reinterpreting a large context summary as evidence. No messages were sent to Claude or any other person by this audit.

## Reproduce and inspect

- [Protocol](protocol.json), [frozen local source hashes](LOCAL_FREEZE.json), [full-field solver](full_field.py), [analysis](analyze.py), and [summary](SUMMARY.json).
- [Exact raw diagnostic archive](raw_outputs.tar.gz) and [entry hashes](RAW_OUTPUTS_MANIFEST.json). These are saved Fourier-mode, conservation and mean-field histories, not full spatial-state checkpoints.
- `python3 unpack_raw.py`, then `python3 analyze.py` reproduces the comparisons. `python3 physical_scales.py` reproduces unit conversions. `python3 plot.py` recreates the figure.
- `python3 run_suite.py` reruns all ten prescribed cases with two local numerical workers; individual cases use `python3 full_field.py CASE_ID`. NumPy/SciPy are needed for simulation/analysis, and Matplotlib for plotting. Reruns may replace local extracted histories; the archive preserves this execution.

R1. [NIST 2022 CODATA full constants table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt), retrieved 20 October 2026. SI defining constants and measured G were used for conversions. Reduced Planck mass was calculated with the 8*pi convention; the measured G uncertainty is retained conceptually and printed digits are not experimental precision.

R2. L. Baudis and S. Profumo, [Dark Matter, PDG 2026](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-dark-matter.pdf), revised August 2025. Selected opening, candidate-property and production sections inspected; no exclusion analysis of this benchmark was performed.

R3. Y.-J. Li, J. Liu and Z.-K. Guo, [Q-balls from thermal balls during a first-order phase transition: a numerical study](https://arxiv.org/abs/2601.19150) (2026). Abstract and selected dynamics/production sections inspected for prior context. Different thermal and dissipative model; no numerical result transferred.
