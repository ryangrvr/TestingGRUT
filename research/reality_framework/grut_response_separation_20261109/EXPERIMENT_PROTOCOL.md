# EX1-01: prospective physical test of the Construction 1 coherence prediction

Prepared 9 November 2026. **Protocol and deterministic ideal verification only. No hardware has been selected, no experimental observations exist in this packet, and no experimental pass is claimed.** This is a hardware-agnostic protocol, not a device-independent certification in the Bell-test sense. Its statistical conclusions depend on stated sampling and apparatus assumptions.

## Purpose and claim

Implement one small GRUT operational model on a real device and test a prediction on data withheld from calibration: two coherent disturbances give final signal-action agreement 13/16; deliberately discarding parity coherence gives 5/8. The difference is 3/16. Both descriptions have the same immediate registered readout and agree on zero- and one-disturbance calibration experiments.

Success would support a finite operational realization under standard quantum mechanics, with a calibrated preparation, operations, readout, and error budget. It would reject the registered population-only reduction at the achieved resolution. It would not establish a new GRUT law, autonomous physical clocks, a detector map for the Schwinger model, a continuum limit, gravitational dynamics, or empirical completeness of GRUT. EX1-01 is a separate experiment; N22, B2, D1, massive-transfer and interacting-parent results retain their existing verdicts.

## Immutable source specification

All repository references below use commit `7ba04fd35e47870cbdb41268b2d7ffd33f567f51`:

- [Construction 1 specification, microscopic model and preparation](https://github.com/ryangrvr/TestingGRUT/blob/7ba04fd35e47870cbdb41268b2d7ffd33f567f51/research/reality_framework/grut_construction_1_20261106/SPECIFICATION.md#L14-L45).
- [Exact gate implementation, model.py lines50-67](https://github.com/ryangrvr/TestingGRUT/blob/7ba04fd35e47870cbdb41268b2d7ffd33f567f51/research/reality_framework/grut_construction_1_20261106/model.py#L50-L67); [preparation, lines87-96](https://github.com/ryangrvr/TestingGRUT/blob/7ba04fd35e47870cbdb41268b2d7ffd33f567f51/research/reality_framework/grut_construction_1_20261106/model.py#L87-L96); [nonselective readout and joint probabilities, lines124-133](https://github.com/ryangrvr/TestingGRUT/blob/7ba04fd35e47870cbdb41268b2d7ffd33f567f51/research/reality_framework/grut_construction_1_20261106/model.py#L124-L133).
- [Frozen population comparator and two-pulse prediction, specification lines73-90](https://github.com/ryangrvr/TestingGRUT/blob/7ba04fd35e47870cbdb41268b2d7ffd33f567f51/research/reality_framework/grut_construction_1_20261106/SPECIFICATION.md#L73-L90).
- [Integrated candidate Eq.(8), intervention and contrast](https://github.com/ryangrvr/TestingGRUT/blob/7ba04fd35e47870cbdb41268b2d7ffd33f567f51/research/reality_framework/grut_integrated_candidate_20261106/GRUT_Integrated_Theory.md#L269-L280); [explicit scope, lines455-465](https://github.com/ryangrvr/TestingGRUT/blob/7ba04fd35e47870cbdb41268b2d7ffd33f567f51/research/reality_framework/grut_integrated_candidate_20261106/GRUT_Integrated_Theory.md#L455-L465).

The original numerical 1e-10 tolerance is a solver criterion. The physical tolerances below are new, prospective EX1-01 criteria and do not replace that historical result.

## Four-qubit work protocol

Use qubits in order `(S, M1, M2, A)`: signal, two memory carriers, actuator. The readout convention is `Z|b> = (-1)^b|b>`. Put `P = Z_M1 Z_M2`, `Rx(theta) = exp(-i theta X/2)`, `Rz(phi) = exp(-i phi Z/2)` and `RZZ(beta) = exp(-i beta P/2)`.

Start in `|0000>`. Apply H on S, H on M1, and CNOT(M1,M2). The prepared state is `|+>_S |Phi+>_M1M2 |0>_A`. State preparation errors belong in the apparatus model; this notation does not certify a laboratory state.

| Step | Gate or operation | Purpose |
|---|---|---|
| 1 | CNOT(S,M2) | Write the signal sign into memory parity |
| 2 | Rz(M1, pi/2) | Historical hold phase; retain it even though the target parity readouts are insensitive to it |
| 3 | Rx(M2, pi/3) if d1=1, matched idle if d1=0 | First disturbance |
| 3a | Identity slot or parity-dephasing channel specified below | Matched intervention comparison |
| 4 | RZZ(M1,M2, 2pi/3) | Relative parity phase |
| 5 | Rx(M2, pi/3) if d2=1, matched idle if d2=0 | Second disturbance |
| 6 | CNOT(M1,A), then CNOT(M2,A) | Set initially blank A to the parity bit without consulting S |
| 7 | Read S,M1,M2,A in Z | Record all four raw bits |

The exact RZZ decomposition, in temporal order, is CNOT(M1,M2), Rz(M2,2pi/3), CNOT(M1,M2). Device compilation may use a native equivalent only after verifying its angle, sign, global/relative-phase convention, timing and routing. Do not cancel the apparatus controls or alter schedules selectively after seeing outcomes.

The selector bits d1,d2 are classical settings fixed before each shot. In the original model they are conserved qubits in fixed computational-basis sectors. Replacing them by fixed settings leaves the ideal conditional work protocol unchanged. Removing the eight-level propagation clock changes the physical claim: this four-qubit version uses external scheduling and does not test the autonomous 512-dimensional Hamiltonian or its energy/clock statements. The full clock experiment would additionally require the eight-level clock, two selector qubits, the specified work-dependent clock links and their physical calibration.

## Parity dephasing without an unintended readout

At slot 3a, apply I or P with equal probability and ignore that selection in subsequent work and the primary analysis. The unconditioned channel is

    D_P(rho) = (rho + P rho P)/2
             = Pi_P,+ rho Pi_P,+ + Pi_P,- rho Pi_P,-.

Thus a classical random operation implements exactly the nonselective parity measurement channel on the work system. No fifth qubit is required. Within each 1,000-shot dephased block, use exactly 500 I and 500 P settings in a precommitted random order. Retain the selection log for checking implementation; do not condition, postselect or correct primary outcomes with it. This balanced design estimates the equal mixture of its two physical branches. The independence assumption applies to shot noise conditional on the schedule.

The nonselective projector channel is the rule in [Luders' original measurement paper, Eq.(8), English translation](https://arxiv.org/html/quant-ph/0403007v2). Here the equality to the randomized I/P channel follows directly by expanding `(I+P)rho(I+P)/4 + (I-P)rho(I-P)/4`; it is also verified in the supplied finite-matrix calculation. This mathematical identity does not certify that a particular device implements the channel.

Use Z on both memory carriers for P. Applying only Z_M1 is the wrong channel and leaves the final contrast intact. Measuring both memory Z values mid-circuit is also a stronger instrument than measuring their product: it can erase coherence inside a parity eigenspace. It is not the registered implementation. A parity ancilla implementation is an alternative only in a separately frozen device adapter, with additional gates and measurement backaction calibrated.

Match the duration and idle exposure of slot 3a in the coherent and dephased protocols. Include residual phases, crosstalk, leakage and any different pulse exposure in the physical model. A virtual Z implementation still needs calibrated phase-frame bookkeeping. Verification must include preservation of within-parity coherence, for example the response of a Bell-state input, as well as loss of cross-parity coherence. Immediate Z-readout equality alone cannot identify the full channel.

## Derivation and control predictions

The logical memory Pauli operators are `tau_x = X_M2`, `tau_y = Z_M1 Y_M2`, `tau_z = P`. Conditional on signal sign, the sign-corrected logical Bloch vector after writing is `(0,0,1)`. The two Rx disturbances and the parity Rz phase give

    z_final = cos(theta1)cos(theta2) - sin(theta1)sin(theta2)cos(beta),
    Pr(S = A) = (1 + z_final)/2.

At theta1=theta2=pi/3 and beta=2pi/3, `z_final=5/8`, so action agreement is `13/16`. Dephasing after the first disturbance deletes the transverse components: `z_final=cos(theta1)cos(theta2)=1/4`, giving `5/8`. The population-only comparator uses independent parity flips with probabilities `q_i=sin(theta_i/2)^2`; its two-pulse agreement is `1-q1-q2+2q1q2=5/8`. This comparator is a restricted memoryless reduction, not every possible classical explanation.

| Condition | Setting/intervention | Ideal action agreement |
|---|---|---:|
| Calibration 00 | Neither disturbance | 1 |
| Calibration 10 | First disturbance only | 3/4 |
| Calibration 01 | Second disturbance only | 3/4 |
| Writer-removal control | 00, writer replaced by matched idle | 1/2 |
| Actuator-removal control | 00, actuator gates replaced by matched idles | 1/2 |
| C11, held out | Both disturbances, identity slot 3a | 13/16 |
| P11, held out | Both disturbances, D_P at slot 3a | 5/8 |

Two additional held-out conditions I_C and I_P stop immediately after slot 3a, before RZZ, the second disturbance and the actuator. They independently read the joint `(S, parity, A)` distribution. They are separate preparations; no mid-circuit readout is inserted into C11. Ideally both have A=0 and probabilities `(s,p)=(0,0),(0,1),(1,0),(1,1)` equal to `3/8,1/8,1/8,3/8`. Therefore their registered TV distance is zero.

The prospective comparison includes levels, a contrast, and the immediate-control distribution. A contrast alone could be created by a bad dephaser; matching immediate readouts alone could hide errors in later operations.

## Two-stage freeze and bounded resources

**Current status: apparatus-independent design prepared.** Hardware execution requires the adapter and calibration freeze below. The absence of a device, pulse schedule and data is an explicit incomplete field, not evidence of a pass.

1. Before collecting calibration data, freeze this design, the device identity and configuration, physical-to-logical mapping, all compiled circuits and schedules, outcome/erasure rules, characterization circuits and counts, parameterized noise models, inference method, drift tests, randomization seeds, stopping rules and analysis code. Store SHA256 digests and a timestamp. The exact physical gates, durations and calibration observation model cannot be filled honestly before selecting a device.
2. Calibration stage: at most **1,000,000 attempted shots total**, inclusive of preflight, characterization, baseline controls and post-run stability sentinels. Reserve 320,000 shots for the 16 computational-basis readout preparations (20,000 each), 100,000 for the five control rows above (20,000 each), and at most 580,000 for independently specified preparation/gate/phase/dephaser characterization and stability sentinels. The device adapter allocates the latter budget before any data. Every attempt consumes the cap. Uncertainty in basis-state preparation must be included when inferring readout errors.
3. Fit only the frozen calibration model family. Allocate calibration failure probability 0.0025 to the five baseline control confidence intervals (0.0005 each), and 0.0025 to all other parameter/characterization and stability confidence objects together. Thus their joint failure probability is at most 0.005. Derive coverage of the latter set from the actual observation model rather than labeling a best fit or bootstrap error bar a certificate. Publish the calibration data, parameter set, prediction bands, numerical enclosure errors and whether qualification passed. Freeze all of them before opening held-out labels/outcomes.
4. Validation stage: **20,000 attempted shots each for C11 and P11**, plus **50,000 each for I_C and I_P**: exactly 140,000 validation attempts. Arrange 140 blocks of 1,000 shots (20,20,50,50 blocks respectively) in a precommitted randomized order, interleaved with calibration-budget sentinels. Block order, dephasing selections and stopping rules must be immutable before acquisition. Use masked condition labels for the primary analyst until qualification and analysis code are sealed.
5. Total cap is **1,140,000 attempts**. No stopping when significance is reached and no extra shots to rescue a failed gate. A device failure can terminate acquisition; an incomplete scheduled sample has the verdict `incomplete`, not an analyzed smaller successful sample. Report every aborted/failed job and its attempted-shot accounting. Provider jobs rejected before any physical attempt are not observations; record them separately.

This bounds data collection, not wall-clock queue time or money. Those require an identified device and a separately accepted booking. No device booking or paid job is authorized or executed by this document.

## Apparatus qualification and uncertainty ledger

The independent calibration must establish a usable physical dictionary for the specific gates and readout, including the sign and magnitude of beta. Zero- and one-disturbance action statistics cannot determine beta. Independent phase-sensitive probes are mandatory: e.g. prepare M1=0 or 1 and M2=+, apply the parity phase, and measure M2 X/Y; the phase changes sign with M1. Characterize Rx using phase-sensitive probes as well as Z populations. Characterize writer/actuator truth tables, Bell preparation and dephaser action on both within- and cross-parity coherence. Routing, idles, spectator interactions and gate composition must be covered by held-out calibration sequences distinct from C11/P11 and I_C/I_P.

For each of the five 20,000-shot calibration controls, require `abs(p_hat-target) + sqrt(log(2/0.0005)/(2*20000)) <= 0.04`, using its table target and the same unconditional valid-and-agree event as validation. These five confidence intervals consume the separate 0.0025 control allocation above. Failure is a qualification miss. Every remaining characterization or stability gate and its joint confidence coverage must be specified in the hardware adapter before calibration; no unnamed discretionary quality test may be introduced after outcomes are known. Post-run sentinels test the already frozen drift envelope; they do not refit or narrow the primary predictions. An out-of-envelope sentinel makes model qualification unresolved even when the raw statistical gates pass.

Use a declared preparation density matrix, physical gate/idle channels and a correlated 16-input by 17-output readout response (16 bitstrings plus erasure), or a more economical model with independently justified assumptions. The forward prediction for raw outcomes is the preparation propagated through the physical operations and then the physical readout response. Include loss/leakage explicitly. Do not apply inverse readout correction, clip negative reconstructed probabilities, discard errors, or renormalize onto detected shots for the primary analysis. Auxiliary corrected estimates may be labeled exploratory.

Confidence sets must include parameter uncertainty, preparation/readout ambiguity, model truncation, compiled-pulse error, dephaser imbalance, crosstalk, leakage, finite duration and drift. Independent validation of a model ansatz is required; a convenient local depolarizing model is not an error bound. Specify a conservative residual model-discrepancy allowance at most 0.01 in primary event probability and at most 0.01 in immediate registered TV. If it cannot be supported within the calibration cap, qualification is unresolved and validation is not opened. This allowance cannot be inferred from the held-out result.

Over the calibration confidence set, freeze forward prediction intervals `[L_C,U_C]`, `[L_P,U_P]`, and joint-distribution sets for the immediate arms. Each primary interval width must be at most 0.04; its allowance-expanded interval must lie inside the corresponding ideal target +/-0.04. The worst permitted immediate-arm TV difference, including the residual allowance, must be at most 0.02. These qualification limits leave statistical room for the validation criteria below. A hardware configuration that cannot qualify is a recorded hardware/model qualification failure; do not weaken the gates after observing it.

For a separately validated gate-error ledger, a useful sufficient bound is `|p_actual - p_target| <= delta_prep + sum(delta_gate) + delta_readout`, where each channel error is half its diamond-norm distance over the relevant input/ancilla class. Full diamond-norm certification is not presumed available. State-restricted errors require proof that the reached states lie in the certified class, including dephased preparations. For unitary finite-pulse implementations, Duhamel's identity gives `||U_actual-U_target|| <= integral ||H_actual(t)-H_target(t)|| dt`; convert this to an event-probability bound before budgeting it. A Rabi fit alone does not establish this operator bound. Known finite pulse profiles may instead enter the forward channel model, with their uncertainty carried into its confidence set. No impulse approximation or zero pulse duration is assumed here.

The population-only comparator must use the same independently calibrated preparation/readout nuisance information, not an intentionally ideal comparator against a noise-corrected quantum model. Its allowed transition/noise family, parameters and prediction band must be frozen using the calibration data alone. For the ideal model it is exactly the independent-flip formula above. With physical noise, a scalar q inferred directly from a raw single-pulse agreement is unjustified unless preparation/readout effects are controlled. Report both models' raw-outcome predictions. Reject only the registered population comparator if its entire prediction band, including its residual allowance, is disjoint from the observed coherent confidence interval. If calibration freedom can fit either held-out outcome, report that the comparator is not discriminated.

## Prospective statistical decision

Count one outcome per attempted shot. Valid outcomes are four raw bits in the fixed order S,M1,M2,A. Use a separate erasure category for an attempted shot without a valid outcome; never silently discard it. The primary Bernoulli success event is `valid AND S=A`. The immediate readout has nine categories: eight `(s,p,a)` combinations plus erasure. Save all 17 raw categories even though the preregistered immediate comparison uses nine.

This profile assumes independent bounded outcomes conditional on the fixed randomized schedule; the bound also accommodates different shot probabilities, in which case it addresses their scheduled average. Time dependence and dependence between shots are different issues. Interleaving and stability checks do not prove independence. If correlated noise invalidates the frozen sampling model, report the statistical verdict as unresolved; a block-dependence alternative must be specified before data, with new sample/error calculations, rather than substituted after seeing a result.

Reserve failure probability 0.005 for validation, split equally across four confidence objects: C11 and P11 event probabilities, and I_C and I_P registered distributions. For a scalar sample mean of n independent [0,1] variables, use `r(n)=sqrt(log(2/alpha_object)/(2n))`. For a k-category distribution, TV equals the largest probability error of an event/subset. A union bound over its `2^k-2` nontrivial subsets gives `r_TV(n,k)=sqrt(log(2*(2^k-2)/alpha_object)/(2n))`. Here `alpha_object=0.00125`, k=9. These are conservative bounds derived from [Hoeffding's bounded-sum inequality](https://www.tandfonline.com/doi/abs/10.1080/01621459.1963.10500830), not asymptotic normal error bars. The executable computes the exact numerical radii.

The four validation gates are:

1. `abs(p_hat_C - 13/16) + r(20000) <= 0.06`.
2. `abs(p_hat_P - 5/8) + r(20000) <= 0.06`.
3. `p_hat_C - p_hat_P - 2*r(20000) >= 0.10`.
4. `TV(p_hat_I_C,p_hat_I_P) + 2*r_TV(50000,9) <= 0.05`.

The bounds are shared confidence objects; reusing them for multiple gates does not require another union penalty. On their simultaneous coverage event, the gates certify the stated ideal tolerances, contrast and immediate-TV tolerance for the sampled device conditions. They are equivalence/separation decisions, not merely failure to reject an ideal null. Combine calibration and validation failure allocations by the union bound for at most 0.01 overall, provided both confidence constructions and their assumptions are valid.

Also compare the observation to the frozen physical prediction set. A primary model discrepancy is excluded at the declared accuracy only if `max(abs(p_hat-L),abs(p_hat-U)) + r + residual_allowance <= 0.06` for both arms; record failure otherwise. Include both the calibrated model and population-comparator comparisons in the report. Statistical-gate success without apparatus qualification or model support is not an experimental GRUT realization.

Calibration controls, dephaser implementation, uncertainty coverage, completed sample counts and stability qualifications must all pass. Report a missed gate as `failed at registered accuracy`; report unavailable/invalid controls, unsupported model bounds, incomplete acquisition or invalid sampling assumptions as `unresolved/incomplete` with the reason. Do not turn a simulator replay, a best-fit model or a single probability contrast into an empirical pass.

## Reproduction and data contract

`experiment_verify.py` uses only Python's standard library. It verifies exact four-qubit state-vector probabilities against an independent logical-rotation formula, the RZZ decomposition, randomized-dephasing/Luders identity, immediate-readout invariance, writer/actuator removal controls, and negative controls with missing parity phase or the wrong local dephaser. Its output explicitly says `empirical_status: not_run`. It neither generates simulated counts nor contacts hardware.

    python experiment_verify.py --output experiment_ideal.json

An optional future raw-count scorer accepts a JSON object with the sole top-level key `conditions`. That object must contain exactly `C11`, `P11`, `I_C`, `I_P`. Each condition maps all 16 four-bit strings `0000` through `1111` and the key `erasure` to nonnegative integer counts. Totals must equal the scheduled count for that condition. No blank count file is presented as data.

`experiment_precheck.py` separately exercises the scorer using deterministic synthetic fixtures in a temporary directory. `experiment_precheck.json` archives its checks for valid synthetic counts, invalid negative/boolean counts, missing conditions/categories, wrong totals and erasure handling. These fixtures are software tests, not simulated experimental observations. `EXPERIMENT_MANIFEST.json` records the exact payload hashes and the pending device adapter, with actual experimental data explicitly absent.

    python experiment_precheck.py --output experiment_precheck.json

    python experiment_verify.py --counts actual_counts.json --output scored_counts.json

The scorer computes the four statistical gates and hashes the count file. It deliberately does not certify apparatus qualification, model confidence-set coverage, the comparator, drift, independence or an overall empirical pass. Preserve raw shot order, timestamps, settings, erasures, pulse/dephaser logs, pre/post sentinels, calibration confidence-set files, inference code, model predictions and the two freeze manifests in addition to aggregate counts. Store immutable results even when a gate fails. Outputs refuse overwrite.

The first justified physical milestone is a successfully calibrated, prospectively tested realization of this finite work protocol. A further distinctive GRUT theory requires an additional fixed physical restriction or selection principle and independent tests; none is introduced by EX1-01.
