# Approach B: experimental boundary and joint measurement contract

2026-10-10. Starting branch revision:
`079e876eed80cf6eb7e0db2e07cd9723ebf1a6ab`.

**Status: BLOCKED — independent apparatus calibration and a certified energy
readout are absent.** A demonstrated platform is identified; no apparatus is
controlled by this workspace and no experiment was performed. This completes
the public-evidence assessment and measurement-contract work, not the requested
experimental calibration or validation milestone.

The design-stage oscillator packet at `19ffa04c` remains byte-for-byte fixed.
This is an additive experimental handoff, not another theoretical investigation,
reconnaissance lane or search for novelty in those equations. Card 2 is unopened;
two slots remain. There is no new physical principle, experimental lock,
canonical integration or independent review certificate.

## One demonstrated platform

Use the documented ETH single-40Ca+ axial-mode implementation as the engineering
reference. Fluhmann and Home's arXiv:1907.06478v1 reports approximately 1.9 MHz
motion, bichromatic state-dependent displacement and fluorescence readout.
Its supplement describes a separate Fock-state force calibration, a typical
shift-rate parameter `c = 36 +/- 0.4 ms^-1`, and motional-frequency recalibration
with about +/-20 Hz accuracy, taking 20 seconds every five minutes. It also
describes ensemble sideband population reconstruction. These are historical
published summaries, not fresh whole-domain error certificates. [P1]

The force convention must be checked against actual displacement data before
`c` becomes our branch-force amplitude. The published example states include
squeezing and superpositions; their occupations cannot be substituted for a
verified centered isotropic thermal preparation. A fitted holdout parameter
cannot be promoted into our calibration. [P1]

The publisher confirms publication in PRL 125, 043602 (2020). Its retrieved
landing page restricts access to the final supplement; the openly accessible
arXiv v1 supplement was inspected instead. No complete independent-calibration
plus held-out shot bundle was retrieved in the bounded search. This is not a
claim that no such data exist. [P2]

Haljan et al.'s separate 111Cd+ experiment provides a caution: its displayed
revival curves fit occupations, detuning/heating and contrast, and distinguish
inferred heating from directly measured heating. Those fits cannot become
independent calibration for our held-out test. Nor may its numbers be combined
with the calcium apparatus to manufacture one measured system. [P3]

## Evidence actually needed from this same apparatus

| Measurement | Available public evidence | Required before a prospective forecast |
|---|---|---|
| Carrier, mode and spin axis | Documented physical implementation | Stable apparatus identity, mode inventory and aligned spin/readout frame |
| Complex force amplitudes | Calibration method and a typical scale | Raw independent calibration runs, displacement normalization, amplitude/phase drift bounds and coverage |
| Physical frequency and detuning | Approximate mode frequency and reported recalibration accuracy | Run-specific values, independent drive detuning, clock error and whole-sequence bounds |
| Initial occupation and state | State-preparation methods | Independently measured occupation, Gaussian/isotropy/centering/correlation certificate and heating bounds for the actual preparation |
| Applied waveform | Demonstrated force control | Measured timing, switching and leakage bounds for both frozen histories |
| Spin detection | Demonstrated fluorescence | Independent `+X,-X` detection probabilities, phase error, trial definition and unbiased trial retention |
| Energy detection | Ensemble reconstruction procedure | Frozen bounded estimator, response calibration, range, mean bias, unresolved-tail moment bounds and sampling design |
| Model discrepancy | Known effective approximations | Independent force nonlinearity, extra-mode, Lamb-Dicke/RWA and heating error enclosures |
| Held-out outcomes | Published results from other protocols | Fresh, separately collected counts/estimator outcomes after this forecast is frozen |

`platform_readiness.json` records the available summaries separately from these
missing inputs. It emits no real calibration, prediction or power claim. Raw
calibration and holdout records have not been found in the current repository.
No investigator was contacted or asked to run an experiment.

This tool inherits the original forecast's disclosed floating-point allowance;
it is not a formally certified interval-arithmetic implementation.

The energy scale is naturally measured in phonons. For this nominal physical
frequency, one quantum is approximately `h * 1.9e6 = 1.259e-27 J`. A desired
uptake of `0.16` quanta would be approximately `2.014e-28 J`; **0.16 is a design
target, not a measured displacement or predicted outcome for that apparatus**.
The calculator does not treat these planning numbers as calibration.

## Energy sampling: an operational requirement, not a new oscillator result

The original packet already predicts `Delta n_j = |alpha_j|^2`. What was missing
was the measurement and sampling contract. Energy has an unbounded spectrum;
Hoeffding cannot simply be applied to the untruncated occupation with a claimed
range `0..C`. A finite-Fock software cutoff is not a detector certificate.

For a chosen mode, freeze an estimator `Y` with an actual per-trial range
`[L,U]` and a cap `C` on the target observable `Q=min(N,C)`. Establish independent
certificates, for both the undriven baseline B and driven sequence A, that

\[
|\mathbb E Y_s-\mathbb E Q_s|\le b_s,\qquad
0\le\mathbb E[N_s-Q_s]\le\tau_s,
\quad s\in\{A,B\}.                              \tag{1}
\]

An independently bounded energy-model discrepancy `b_model` is also required.
A bound on tail *probability* alone does not bound its contribution to mean
energy. If (1) is unavailable, the full mean-energy test stays blocked. A capped
observable can be tested separately, but must not be labelled the full energy.

The estimator need not report a phonon number on each shot. For a frozen
sideband response inversion, for example, randomized probe label `l` with
probability `q_l>0` and binary outcome `X_l` can give
`Y=b+ w_l X_l/q_l`. Its range follows from its realized values, including any
negative inversion weights, not from `C`. The calibration must establish (1)
for that decoder. Ill-conditioned inversion enlarges `U-L` and the required
sample count. This packet supplies no response matrix or such decoder
certificate for the chosen apparatus.

Let the frozen force forecast enclose the physical uptake as `[d_-,d_+]` in
phonons. Then the mean of the reported estimator difference is enclosed by

\[
\mathbb E(Y_A-Y_B)\in
[d_- -\tau_A-b_A-b_B-b_{\rm model},\quad
 d_+ +\tau_B+b_A+b_B+b_{\rm model}].              \tag{2}
\]

This follows by writing `N=Q+(N-Q)` in both preparations. Notice the asymmetric
tail terms: the driven tail lowers the capped uptake; the baseline tail raises
it. No Gaussian tail is inferred from the outcome being tested.

For `G=2K` baseline/driven groups covering `K` frozen energy tests, allocate total
energy sampling failure probability `a_E`. For group `s` with `N_s` independent
trials and frozen width `W_s=U_s-L_s`, use

\[
e_s=W_s\sqrt{\log(2G/a_E)/(2N_s)}.               \tag{3}
\]

Reject a test if its observed difference falls outside (2) widened by
`e_A+e_B`. Hoeffding and the union bound control simultaneous energy sampling
error by `a_E`. Cross-group independence is not required by that union bound;
the implementation nevertheless uses separately defined baseline/driven groups
and does not silently reuse a baseline as two independent samples.

Keep the original spin allocation `a_S` and the force-calibration simultaneous
coverage failure `a_C`. The energy certificate has its own simultaneous failure
bound `a_D`. The full spin/energy false-rejection bound is at most
`a_S+a_E+a_C+a_D`, conditional on the stated physical model, independent-shot
and discrepancy assumptions. Two tests cannot each spend the entire common
sampling budget and then claim its original confidence.

## Power must refer to a specified alternative

The companion tool reads the existing forecast implementation without changing
it. It reports whether every possible bright fraction is accepted. It also
evaluates frozen operational contrasts with a requested total miss probability
`beta_total`, split equally over the frozen contrasts. Below `beta` is the
allocated per-contrast miss probability; a union bound gives simultaneous power
at least `1-beta_total` if all contrast conditions pass.

For two predicted mean intervals, let `d_min` be their minimum possible
separation, zero if they overlap. For an operational equality null, let
`e_null` be the sum of the applicable group sampling allowances. Let

\[
e_{\rm alt}=\sum_{s=1}^{r} W_s
\sqrt{\log(2r/\beta)/(2N_s)}.                   \tag{4}
\]

If `d_min > e_null+e_alt`, the specified model, uniformly over its admitted
mean intervals, rejects equality with probability at least `1-beta`, conditional
on the certificates. Otherwise the calculator reports
`INSUFFICIENT_CERTIFIED_POWER`; it does not say the real experiment has no power.
For spin `r=2, W_s=1`. For two reported energy-uptake differences, `r=4` and the
estimator widths apply. The energy contrast concerns **reported estimator
means**, not an equality theorem for physical energies in the presence of
unknown decoder biases. The individual energy rejection rule (2) remains the
physical model test.

A narrow acceptance interval alone is not a general power certificate: an
alternative arbitrarily close to its edge can be hard to detect. Equations
(3)–(4) expose the tradeoff with sample counts and estimator conditioning; they
do not certify actual calibration, hardware stability or novelty.

## Executable handoff and stop condition

`joint_contract.py` accepts the same calibration and protocol schemas as the
fixed design packet plus an independent energy-measurement contract. Its
software-control fixtures are explicitly invented and never substituted into
`platform_readiness.json`. There are no fitting or outcome arguments. To run
the contract checks:

**28 author contract checks passed.** They exercise probability allocation,
asymmetric tail bounds, estimator conditioning, vacuous intervals and rejected
inputs. They add no experimental evidence or new oscillator-physics claim.

```bash
bash research/approach_b_experiment/reproduce.sh
```

For independently supplied files:

```bash
python research/approach_b_experiment/joint_contract.py \
  --calibration /absolute/path/calibration.json \
  --protocol /absolute/path/frozen_protocol.json \
  --measurement /absolute/path/energy_measurement_contract.json \
  --output /absolute/path/frozen_joint_contract.json
```

The actual platform remains **BLOCKED** until the table's run-specific evidence
and energy decoder/tail contract exist. Predicted numbers and minimum sample
counts from software controls are not an experimental forecast. Once an actual
bundle is supplied, freeze the input/code/evidence hashes and prospectively
collect all predeclared trials, including failed preparation/detection events
under the frozen retention rule. A change of calibration, decoder, mode
inventory or discrepancy budget after examining holdout data requires a new
holdout.

No further oscillator algebra or software agreement will unblock a missing
measurement. This handoff identifies the platform, implements the previously
absent measurement/power criteria, and records exactly why calibration and
validation cannot currently be claimed. The fundamental-theory question remains
separate; these statistical and experimental requirements add no GRUT law.

## Sources and retrieval scope

- [P1] Fluhmann and Home, arXiv:1907.06478v1, main text and supplement,
  https://arxiv.org/pdf/1907.06478v1 . Specific evidence: main p.2 and supplement
  calibration routines, Fock reconstruction discussion, and calibration/fit table.
- [P2] Publisher record for PRL 125, 043602 (2020),
  https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.043602 .
- [P3] Haljan et al., arXiv:quant-ph/0411068, figures 2–3 and their captions,
  https://arxiv.org/pdf/quant-ph/0411068 .
- The design-stage packet is pinned in `MANIFEST.json` by its existing hashes.
  No published fitted parameter, paper figure or inaccessible supplement is
  claimed to be a fresh measurement or an independently reanalysed raw dataset.
