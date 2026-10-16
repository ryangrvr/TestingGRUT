# WM2-S5: history-preserving formation

Read [REPORT.md](REPORT.md), then [MATHEMATICS.md](MATHEMATICS.md) and
[PROBABILITY_BOUND.md](PROBABILITY_BOUND.md). The two charters were
published before their respective trajectories. [MODEL.md](MODEL.md)
records the unchanged physical law, preparation, resources and record
cost; [HANDOFF.md](HANDOFF.md) preserves the next boundary.

The tested five-event balanced-disorder families now have controlled
finite formation lower bounds under P1. The original stricter S5 runs
and all six-event resource stops remain in the archive. All-size,
all-pattern independent-disorder formation remains OPEN.

## Quick verification

From the repository root with Python, NumPy and SciPy:

```sh
python research/reality_framework/wm2_history_formation_20261016/reproduce.py
```

This verifies frozen charter hashes, every packet hash, every stop and
error gate, original root weights, the unchanged original S5 data hash
recorded by P1, and the complete audited probability ledger.

## Full reproduction without changing the archive

```sh
OPENBLAS_NUM_THREADS=1 python research/reality_framework/wm2_history_formation_20261016/reproduce.py --recompute --output /tmp/grut_s5_rerun.json
```

This copies S5 and its S4 source dependency to a temporary workspace,
runs the frozen campaign and independent complete-matrix controls,
then runs P1 and exports the new evidence. Archived files stay unchanged.
It can take tens of minutes. Wall-time-limited instances may stop at
different steps on different hardware; that is a computational limit,
not an alteration of the physical formation criterion. No positive
status is inferred from an expected archived result.

## Data schema and numerical limits

* `RESULTS.json`: original S5 record, retained verbatim, including the
  pilot reuse note and every incomplete target-time instance.
* `PROBABILITY_RESULTS.json`: P1 record and its hash of original S5 data.
  The inherited `finite_instance_gate` is the original <=.20 envelope;
  the separate `P1_gate_pass` concerns the P1 certificate.
* `AUDITED_RESULTS.json`: physical interpretation with the conservative
  contact-threshold term included in every early partial ledger.
* `ASSESSMENT.json`: proved, evaluated and open claims.
* `EXECUTION.json`, `VERIFICATION.json`, `manifest.json`: execution
  disclosure, verification summary and content hashes.

Incomplete populations are intermediate only. Unknown roots keep zero
lower bound at their original weight and are not renormalized away.
The physical receiver is L_N, while the active numerical rank is at
most 24. Geometry histories and charge are exact in the representation;
field-mode and time discretization errors are explicitly bounded.

Ordinary double precision is used. Analytic reduction bounds are not
outward-rounded floating-point certificates or independent external
validation. The reported response protocol is an analytic combination
with S4; S5 did not propagate a giant joined one-probe matrix. The
model uses hbar=1 and supplied receiver/coupling units. No laboratory
scale, selected dimension or distinctive empirical prediction is derived.
