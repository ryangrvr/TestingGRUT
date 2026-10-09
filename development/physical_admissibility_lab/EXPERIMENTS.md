# Operational tests and experimental evidence

**DRAFT — no physical experiment performed, no Stage-3 lock certified.**
These are standard discriminators with explicit comparators. Historical
demonstrations are not advertised as the strongest October 2026 bounds.
Primary-source inspection depths are in [SOURCES.json](SOURCES.json).

## Finite, operational memory test from the inherited control

Choose each of the four input words equiprobably, write it via the specified
shift-register intervention, then apply the identical two-step 00 probe.
Success means decoding both bits correctly from the two records. All
preparation-dependent degrees of freedom available during the probe,
including controller labels, environment, leakage levels and correlated
timing, belong to the declared retained capacity.

Exact four-state comparator: p_success=1. Any d≤3 quantum retained-state
comparator with independent fresh probe ancillas: p_success≤3/4. A qubit:
≤1/2. This bound does not rely on response-matrix rank or a fitted rate.
It follows from T2. Software or classical register implementation is
straightforward, and known prepare/measure dimension demonstrations [P13]
support platform feasibility in principle. This package certifies no
particular apparatus resource boundary or detector interface.

If N trials are independently drawn under a frozen stable device model,
a conservative one-sided bound is
p_lower=p_hat−sqrt[ln(1/α)/(2N)] (Hoeffding). At α=0.01 and N=1,000 the
statistical subtraction is about 0.048; p_hat=0.90 would exclude d≤3
before separately bounded systematic errors. That is an illustrative
statistical design under independence, not achieved data or a universal
sample-size claim. With drift/device memory use a justified martingale or
block/randomization certificate instead; optional stopping is not covered.
Refutation of a **specific four-state ideal register implementation** is
an incompatible noisy record distribution after its error model is frozen.
Failure to cross the witness threshold does not establish low dimension.

The d≤3 exclusion is already quantum information theory. A larger ordinary
fixed register survives. It does not test a novel generative law and cannot
be promoted to a Card discriminator solely by renaming memory as Π.

## Other constrained classes that observations can distinguish

| Test / evidence | Same-operation comparator excluded or restricted | Operational assumptions and practical boundary |
|---|---|---|
| Bell correlations: Hensen et al., 2015 [P15], actual experiment | Bell-local setting-independent models obey CHSH≤2; ordinary quantum theory permits up to 2sqrt2 | Spacelike randomized settings, trial definition and memory-robust statistics matter. This does not reject fixed quantum dynamics or demonstrate evolving laws |
| Prepare/measure dimension: Hendrych et al., 2012 [P13], actual experiment | Specified bounded classical/quantum dimensions | Leakage and retained side information are part of dimension; the T2 four-word bound is an exact standard resource test, not a claim this prior experiment implemented this register |
| Contextuality: Kirchmair et al., 2009 [P23], actual experiment | The stated noncontextual model with its compatibility/operational assumptions | Disturbance/compatibility must be controlled. A classical controller retaining measurement context can evade the restricted comparator; contextuality does not exclude all fixed hidden-memory dynamics |
| Quantum erasure: trapped-ion single-atom experiment, 2018 [P16], actual experiment | Thermal product-bath erasure with correctly measured entropy, heat and correlations obeys T4 | Heat is not work; all side information and ancilla renewal count. It constrains a priced process, not a temperature-independent cost of every memory event |
| Non-Markovian characterization: quantum-processor process-tensor experiment, 2021 [P17], actual experiment | System-only memoryless instrument dynamics under the frozen interventions | Environmental memory, correlated control and leakage supply fixed enlarged comparators. Non-Markovianity alone is not a breakdown of fixed laws |
| Correlation propagation: trapped-ion long-range experiment, 2014 [P18], actual experiment | The specified interaction graph/range/norm and appropriate propagation bound | Its long-range dynamics cannot be compared to a short-range LR bound and called a locality violation. T8 supplies a three-spin control against an exact-cone claim; relativistic causality is different |
| Gravity-mediated entanglement proposals, 2017 [P19,P24], proposals | Classical local mediator/LOCC class that cannot entangle initially separable masses under the mediation assumptions | Must bound nongravitational channels, direct interactions and witness errors. Does not exclude every classical-gravity hybrid description or certify a full quantum-gravity theory. No performed experiment is claimed from these proposal papers |
| Collapse/noise tests [P25], established model framework | Specific CSL/collapse parameter ranges vs calibrated ordinary environmental decoherence | A mass-, separation- and time-dependent prediction plus nuisance bounds is required. An unrestricted extra noise process can imitate loss of visibility. No new collapse simulation or October 2026 strongest bound is claimed |
| MICROSCOPE final test, 2022 [P27], actual experiment | Differential acceleration beyond the reported uncertainty for the measured material pair, under calibrated gravity/apparatus model | Equivalence-principle tests constrain specified composition-dependent couplings. They do not select all gravitational actions or a universal Π; no latest-limit ranking is supplied |
| Optical Lorentz clock comparison, 2018 [P28], actual experiment | The specified orientation/time-dependent Lorentz-violating coefficients | Clock sensitivities, rotation/background fields and nuisances are independent inputs. A drift in a ratio needs a specified physical model, not an architecture interpretation |
| EFT scattering/positivity [P08], theorem and possible inference | UV comparator obeying T10 when a negative coefficient is reliably inferred | Cutoff, other operators, pole subtraction and dispersion hypotheses matter. No collider fit or empirical positivity violation is produced here |

## Why this is not yet an experimental lock

These sources establish that physically restricted comparator classes are
operationally useful. They do not certify the program's specific common
carrier, whole-chart error bound, generated readout/interface, non-equilibrium
restriction, guard-band or frozen platform. None gives a new relation
beyond the listed standard premises. Card 1 remains CLOSED. Software
fixtures certify exact model arithmetic; they cannot independently certify
apparatus capacity, locality, reservoir state or preparation independence.
