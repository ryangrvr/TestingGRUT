# Gravitational quantum information: apparatus feasibility decision

**Verdict: NOT YET FEASIBLE on the documented apparatus.**

Assessment date: 11 October 2026 UTC. This is one bounded literature-based
experimental assessment, not an experiment, a new law, or independent scientific
clearance. Stage 3 stays inactive; Card 3 is unopened and unauthorized. PR #12
stays draft. The accepted gravitational-information note and all earlier frozen
packets are preserved.

## Apparatus selection and evidence boundary

Select the **Leiden superconducting magnetic-levitation/SQUID apparatus** as the
strongest gravity-ready lead in this bounded comparison: it has actual
gravitational-force detection and a documented subsequent cooling implementation.
This is a task-specific choice, not a claim that it leads every platform metric.

Other current experimental evidence does not supply its missing capabilities:

| Independently documented experiment | Relevance to this decision |
|---|---|
| Bound atomic clusters up to 608 AMU [3] | Coherent preparation; no demonstrated mutual gravitational quantum interaction in that apparatus |
| 2026 one-milligram torsional pendulum [4] | Cooled to 240 microkelvin; reported occupation about 280,000. No completed gravitational entanglement measurement |
| Levitated oscillator–optical-field entanglement [5] | Actual mechanical quantum control; interaction/readout is optical, not evidence of gravitational mediation |

No parameter in those experiments is imported into the selected apparatus.

**Source ledger.** The 2024 gravity configuration [1] used a 0.43-mg detector,
26.7-Hz mode and kilogram rotating sources at roughly 0.5 m. It reported a
30-aN gravitational signal and 0.5-fN/sqrt(Hz) force noise. Its response/model
amplitude factor was 0.35 +/- 0.02; platform motion matters. Supplement A
documents electrical and magnetic shielding.

The 2026 cooling preprint [2], v2 dated 1 July, describes a changed configuration,
not a concurrent repeat of that gravity measurement. Modes 3/4 have frequencies
50.59/67.98 Hz, bare Q values 3.8/5.5 million and simultaneous cooled temperatures
7.1/6.6 mK. Detection noise is approximately 1 pm/sqrt(Hz). Its mass is internally
inconsistent: 0.35 mg in the setup versus 0.38 mg in the conclusion. Supplement A
also qualifies mode-axis identification. Treat these as author-reported results;
the retrieved record does not establish journal review of [2].

## What the measured numbers imply

The accompanying [arithmetic record](GRAVITATIONAL_FEASIBILITY_NUMBERS.json)
and [short script](gravitational_feasibility_numbers.py) reproduce these
calculations using exact SI constants. They contain no apparatus simulation or
new fixture campaign.

Define the energy-equivalent occupation and quantum temperature scale by

\[
n_E=\frac{k_B T_{\rm mode}}{h f},\qquad T_q=\frac{h f}{k_B}.
\]

These are energy diagnostics. An inferred mode temperature does not certify
a thermal quantum state, purity, or spatial superposition.

| Calculation using the 2026 configuration alone | Mode 3 | Mode 4 |
|---|---:|---:|
| Quantum temperature scale | 2.428 nK | 3.263 nK |
| Energy-equivalent occupation | 2,924,293 | 2,022,973 |
| Zero-point position scale, using the textual mass envelope | 0.661–0.688 fm | 0.570–0.594 fm |
| Bare mechanical energy ringdown, Q/(2 pi f) | 11,955 s | 12,877 s |

The mass envelope is a textual ambiguity, not a confidence interval. Long
mechanical ringdown is **not a measured quantum coherence lifetime**. Feedback
also changes damping; substituting the bare Q into a cooled-state lifetime
would conflate operating conditions.

For the near-ground-state oscillator-control route, the measured energy must
fall by roughly two to three million times to approach one quantum. This is not
a universal requirement for every possible gravitational experiment; it is the
size of the shortfall along the cooling route examined here.

The published linear-feedback model gives, for one-sided force and detection
PSDs,

\[
n_{\min}=\frac{\sqrt{S_F S_x^{\rm det}}}{2\hbar}.
\]

Its reported current floors are approximately 270,000 and 570,000 quanta. Within
that model, reaching a floor of one requires reducing the **product of force
and detection amplitude spectral densities** by those factors, or their PSD
product by about 7.29e10 and 3.25e11. More feedback gain alone cannot remove that
floor. This is a conditional engineering diagnostic, not a certified quantum
noise limit or a newly derived physical restriction. After such improvement,
back-action and state certification would still need independent measurement.

For comparison, use only the **2024** force configuration. A known-phase sinusoid
of amplitude F in locally white, stationary one-sided noise S_F has fitted
amplitude standard error sqrt(S_F/T). Its reported 30-aN signal and
0.5-fN/sqrt(Hz) noise therefore give

\[
T_{5\sigma}=\left(\frac{5(0.5\times10^{-15})}{30\times10^{-18}}\right)^2
=6,944\ {\mathrm{s}}\simeq1.93\ {\rm hours}.
\]

That is a conditional **ordinary force-detection** averaging estimate, not a
quantum-signature forecast. It assumes stable drive phase and noise over the
integration, independent calibration, and negligible bias. It does not absorb
the response-scale discrepancy or certify these assumptions. Combining this
force floor with the 2026 temperatures would invent an unmeasured configuration.

## Complete measurement requirements

| Requirement | Apparatus-specific assessment |
|---|---|
| Massive coherent source | No certified spatial coherence, branch separation or preparation fidelity in the selected experiment; rotating classical masses do not supply it |
| Gravitational interaction | Classical force established in the earlier configuration; no independently calibrated quantum source–probe interaction in the cooled configuration |
| Coherence during interaction | No measured protocol-specific coherence lifetime; classical Q cannot substitute |
| Joint quantum readout | SQUID displacement detection is not a demonstrated pair of calibrated quantum measurement channels or joint state tomography |
| Preparation/readout error | No measured quantum source preparation fidelity, measurement confusion matrix or quantum detection efficiency for this task |
| Statistical power | Quantum signal contrast, joint outcome distribution and full uncertainty budget are absent; no defensible finite shot count can be assigned |

The kilogram-source experiment can establish a gravitational response but
does not test the coherent-source versus matched-mixture comparison. Under the
accepted branch-preserving controlled interaction, their probe-only outputs
are equal. Repeating that measurement cannot recover a missing joint distinction.
The previous note already proves this limitation; no additional witness is
derived here.

## Background accounting

The following are **uncertified contributions**, not quantities assumed to be
zero. Their required limits must refer to the quantum signal in the actual
geometry; that signal has not yet been established.

| Contribution | Consequence for feasibility and required apparatus evidence |
|---|---|
| Electromagnetic | The sensor is a permanent magnet with SQUID, coils and active feedback. Existing shielding is useful, but a second nearby mass introduces new paths. Measure residual magnetic gradients, charge/patch forces, coil cross-talk and shared electrical feedback. Bound their force gradients and coherent coupling, not only their mean force; common control must not become an entangler mistaken for gravity |
| Casimir and surface forces | The distant-source measurement supplies no calibration for nearby quantum masses. Determine actual gaps, surfaces and intervening shields; measure geometry-dependent forces/gradients with gravity-insensitive controls. A conductive shield can itself introduce motion, electromagnetic noise or shared mechanical paths |
| Vibration | Measure source-synchronous support motion, cross-spectra, transfer functions and nonlinear upconversion in the final geometry. Cryostat isolation does not by itself bound communication through shared supports or controller channels |
| Thermal/environmental | Measure reheating with feedback disabled/enabled, gas-collision and photon noise, magnetic fluctuations, and the spectrum relevant to the intended preparation. Refrigerator temperature and inferred motional temperature are distinct; neither alone certifies coherence decay |
| Preparation | Verify initially separable source/probe states and an independently prepared coherent source and matched mixture. Specify measured branch geometry and phase stability; otherwise pre-existing entanglement or imperfect matching can create the contrast |
| Readout | Calibrate local gains, efficiencies, imprecision, back-action, drift and detector cross-talk independently of held-out outcomes. Include controller records and systematic bias in the joint uncertainty; classical displacement SNR does not certify quantum observables |

No quantitative gravitational-attribution certificate is available for these
combined effects. A plausible shielding story cannot close that gap.

## Conventional comparisons and interpretation

| Framework | Potential distinguishing observation | Present apparatus |
|---|---|---|
| Conventional low-energy quantum gravity [6] | Supplied weak-field interaction can generate correlations/entanglement from suitable quantum preparations; predictions require geometry, states and environmental couplings | Current mean-force response is compatible; relevant quantum preparation/readout absent |
| Mean-field or explicitly local classical communication comparison | Under a stated separability-preserving interaction model, independently certified entanglement from separable inputs would exclude that model | No such observation or preparation certificate |
| Time-local completely positive classical–quantum dynamics [7] | Independently measured back-reaction, decoherence and classical diffusion can test its joint trade-off with null-space compatibility | No complete joint coefficient measurement, noise subtraction or time-locality certificate |
| Broader classical–quantum models with memory | Must be specified and compared operationally; the preceding time-local trade-off is not a universal bound on this class | Existing force/cooling data do not select among these possibilities |

A positive future entanglement result would require nongravitational-coupling
exclusion before gravitational attribution, and explicit mediator assumptions
before a claim about quantization. A null result could reflect noise, preparation
or readout. Recoherence after an implemented inverse would test recovery of
accessible records; it would not alone certify reversibility of omitted field
and apparatus degrees of freedom. Neither outcome would automatically establish
a GRUT law.

## Evidence that would change this decision

Changing the verdict requires **one configuration** with independently certified
mass/geometry and gravitational coupling; suitable source coherence and initially
separable probe preparation; measured lifetime through the interaction; calibrated
joint quantum readout; and quantitative background bounds. It must produce a
predeclared signal separation exceeding total calibration/systematic/statistical
uncertainty with adequate power against specified conventional alternatives.

For the assessed cooling route, measured occupation near one and demonstrated
quantum control would address the quantified energy/noise deficit, but would not
alone complete that certificate. No single cooling factor can substitute for
the missing coherent source, joint readout or attribution controls.

**Stop at this decision.** No new law, candidate search, experimental lock,
canonical integration, bank claim, or issue queue is created. Independent kit
and failed-card reviews remain outstanding. Two cards consumed; Card 3 untouched.

## Primary sources and access

1. Fuchs et al., *Measuring gravity with milligram levitated masses*, Science
   Advances 10, eadk2949 (2024): [published full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10889343/).
   [Author manuscript including supplements A–F](https://arxiv.org/pdf/2303.03545v2),
   especially pp. 3–8 and 14–16. The publisher supplement link encountered an
   access challenge; the manuscript supplement was inspected instead. Raw
   data are linked by the paper; they were not independently reanalysed here.
2. Uitenbroek, Langendorff and Oosterkamp, *Picometer control of a levitating
   milligram gravity sensor*, [arXiv:2605.28479v2](https://arxiv.org/pdf/2605.28479v2),
   pp. 2–5 and supplement A. Full author manuscript inspected; no independent
   laboratory replication or review certificate supplied.
3. Zhang et al., *Scalable generation of massive Schrödinger cat states via
   quantum tunnelling*, Nature Physics 22, 1037–1043 (2026):
   [publisher abstract](https://www.nature.com/articles/s41567-026-03281-9).
   Abstract-level scope only; no inaccessible apparatus parameters used.
4. Agafonova et al., *One-milligram torsional pendulum toward experiments at
   the quantum-gravity interface*, Communications Physics 9, 80 (2026):
   [full paper](https://www.nature.com/articles/s42005-026-02514-w.pdf), pp. 1–4.
5. Deplano et al., *Stationary entanglement of a levitated oscillator with an
   optical field*, Science 394, 113–117 (2026):
   [author manuscript](https://arxiv.org/pdf/2602.03456). Used solely to distinguish
   optical entanglement from a demonstrated gravitational interaction.
6. Donoghue, *General relativity as an effective field theory: The leading
   quantum corrections*, Physical Review D 50, 3874 (1994):
   [author record](https://arxiv.org/abs/gr-qc/9405057). This assessment does not
   recompute quantum-loop corrections or claim that they are measurable here.
7. Oppenheim et al., *Gravitationally induced decoherence vs space-time
   diffusion: testing the quantum nature of gravity*, Nature Communications
   14, 7910 (2023): [full paper](https://www.nature.com/articles/s41467-023-43348-2.pdf).
   Comparison restricted to the stated completely positive time-local class.
