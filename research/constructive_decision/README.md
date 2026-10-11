# Calibrated conventional effective-model forecast

Read [DECISION.md](DECISION.md) for the scientific decision, exact derivation,
calibration contract, conventional comparison and boundaries. This is not a
new GRUT law, new reconnaissance lane, Card 2 or a gravity prediction.

## Reproduce the software controls

From the repository root, with Python 3.11+, NumPy and SciPy (tested versions
are pinned in `requirements.txt`):

```bash
bash research/constructive_decision/reproduce.sh
```

The forecast CLI itself uses only the Python standard library. Its example
inputs are **software controls with invented values**, not apparatus data,
literature parameter extraction, an experimental lock or independently
calibrated evidence. `software_control_forecast.json` retains that label.
The verification script independently exponentiates finite-Fock Hamiltonians;
its truncation checks and analytic controls test the implementation.
They do not constitute independent scientific review.

## Make a prospective forecast from actual calibration

Copy the calibration and protocol schemas to separate files, replace the
software values with independently calibrated values, and set calibration
`purpose` to `INDEPENDENT_CALIBRATION_FORECAST`. Supply `evidence` entries with
`id`, `role`, local `path` (relative to the calibration file) and the SHA-256 of
each actual evidence file. Required references are:

| Input object | Field | Evidence role |
|---|---|---|
| Calibration | `clock_evidence_id` | `clock_and_waveform` |
| Every mode | `frequency_evidence_id` | `frequency` |
| Every mode | `occupation_evidence_id` | `occupation_and_state` |
| Every calibrated drive | `force_evidence_id` | `force` |
| Readout | `readout_evidence_id` | `readout` |
| Readout | `model_error_evidence_id` | `model_error` |

The role checks and byte hashes do not establish that the experiments or
certificates are scientifically adequate. Independent assessment is still
required. Claimed calibration intervals must cover the whole frozen chart
with the declared simultaneous failure probability. Model discrepancy must be
bounded independently; setting it to zero is a physical assertion.

The protocol selects a measured drive ID and its sign for each mode/segment,
with positive durations and absolute timing errors. An absent mode command
means **exactly zero force** in the model. Use a calibrated off-drive entry if
leakage needs a nonzero force-error bound. Frequency/detuning are radians per
second, force amplitudes are inverse seconds, occupation is dimensionless,
and segment durations are seconds. Physical mode frequency and drive detuning
are separate quantities. State and frame assumptions are in the decision.

```bash
python research/constructive_decision/forecast.py \
  --calibration /absolute/path/calibration.json \
  --protocol /absolute/path/frozen_protocol.json \
  --output /absolute/path/prospective_prediction.json
```

This command has no fitting or observation input. Freeze its input files,
evidence and output before collecting the holdout. The result predicts mode
energy uptake, qubit visibility and raw spin detection probability. Its
`accepted_observed_bright_fraction` is a predeclared rejection threshold
under the stated independent-shot, calibration and model-error assumptions.
No energy-readout sampling certificate is implied by that spin threshold.

A held-out bright count divided by the frozen shot count lying outside that
interval rejects the combined model/calibration/error contract. Record the
outcome even if the model fails. Refitting parameters or adding modes after
seeing the result requires fresh holdout data. Do not interpret a disagreement
as a new fundamental law or a refutation of quantum mechanics.

No actual calibration bundle is present in this repository packet.
**Experimental status: AWAITING INDEPENDENT CALIBRATION.**

`MANIFEST.json` pins this packet and its repository source bytes. Prior packet
manifests and the Card-1 ledger are left intact.
