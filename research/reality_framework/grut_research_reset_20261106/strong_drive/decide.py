"""Apply the published comparison rule to saved strong-drive runs."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from propagation import BINS, CASES, FINE, Grid, HERE, Pulse, tangent, write_json

ARMS = ("propagating", "local")
AXES = ("time", "spacing", "longitudinal", "transverse")


def comparison(runs, metric):
    shifts = {}
    crossed, regulator, fine = {}, {}, {}
    for arm in ARMS:
        base = runs["base"][arm][metric]
        shifts[arm] = {case: runs[case][arm][metric]-base for case in AXES}
        fine[arm] = runs["fine"][arm][metric]
        crossed[arm] = fine[arm] - base - sum(shifts[arm].values())
        regulator[arm] = runs["regulator"][arm][metric]-fine[arm]
    error = sum(abs(delta) for arm in ARMS for delta in shifts[arm].values())
    error += sum(abs(crossed[arm])+abs(regulator[arm]) for arm in ARMS)
    work_closure = {
        arm: abs(runs["fine"][arm]["work"]-runs["fine"][arm]["energy"])
        for arm in ARMS}
    if metric == "work":
        error += sum(work_closure.values())
    scale = max(abs(v) for v in fine.values())
    gap = abs(fine["propagating"]-fine["local"])
    return dict(fine=fine, signed_one_axis_shifts=shifts,
                signed_crossed_remainder=crossed, signed_regulator_shift=regulator,
                fine_work_closure=work_closure, allowance=error, gap=gap,
                relative_gap=gap/scale, relative_allowance=error/scale,
                gap_over_allowance=gap/error,
                large=gap > max(.01*scale, 3*error),
                small=gap+3*error < .01*scale)


def verify(runs):
    worst = dict(sector_number_defect=0., spectrum_relative_error=0.,
                 work_partition_relative_error=0., work_trace_relative_error=0.)
    for case, row in runs.items():
        if row["case"] != case or row["config"] != vars(CASES[case]):
            raise ValueError("Run configuration mismatch for " + case)
        for arm in ARMS:
            r = row[arm]
            if not 0 < r["pressure_free"] < r["energy"]/3:
                raise ValueError("Free outgoing moments out of bounds")
            for observable in ("number", "energy"):
                e = abs(sum(r["spectrum_"+observable])/r[observable]-1)
                worst["spectrum_relative_error"] = max(e, worst["spectrum_relative_error"])
            e = abs((r["positive_work"]+r["negative_work"])/r["work"]-1)
            worst["work_partition_relative_error"] = max(e, worst["work_partition_relative_error"])
            e = abs(r["work_trace"][-1][-1]/r["work"]-1)
            worst["work_trace_relative_error"] = max(e, worst["work_trace_relative_error"])
            worst["sector_number_defect"] = max(r["sector_number_defect"],
                                               worst["sector_number_defect"])
    if max(worst.values()) > 1e-8:
        raise ValueError("Outgoing number or ledger inconsistency " + str(worst))
    return worst


def diagnostics(runs):
    """Post-run diagnostics, outside the registered decision rule."""
    pulse = Pulse()
    duration = pulse.duration

    def transform(frequency):
        out = np.zeros_like(frequency)
        for j, coefficient in enumerate((1/16, 1/4, 3/8, 1/4, 1/16), start=-2):
            carrier = pulse.frequency + j*np.pi/duration
            out += coefficient*duration*(np.sinc((frequency-carrier)*duration/np.pi)
                                         + np.sinc((frequency+carrier)*duration/np.pi))
        return out

    grid = Grid.quadrature(replace(FINE, bloch_nodes=32, radial_nodes=12))
    compact_tangent = {arm: tangent(grid, transform, local=arm == "local")
                       for arm in ARMS}
    opened, cyclic = (runs[case]["propagating"] for case in ("fine", "regulator"))
    first, second = np.asarray(opened["work_trace"]), np.asarray(cyclic["work_trace"])
    if opened["work_trace_dt"] != cyclic["work_trace_dt"] or not np.array_equal(
            first[:, 0], second[:, 0]):
        raise ValueError("Regulator traces must sample the same times")
    difference = np.abs(first[:, 1:]-second[:, 1:]).max(axis=0)
    return dict(compact_tangent_per_amplitude_squared=compact_tangent,
                max_sampled_regulator_power_shift=float(difference[0]),
                max_sampled_regulator_cumulative_work_shift=float(difference[1]),
                positive_work_regulator_shift=cyclic["positive_work"]-opened["positive_work"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    path = parser.parse_args().output
    if not json.loads((path / "qualification.json").read_text())["passed"]:
        raise SystemExit("Reserved run was unqualified.")
    runs = {case: json.loads((path / (case + ".json")).read_text()) for case in CASES}
    checks = verify(runs)
    metrics = {metric: comparison(runs, metric) for metric in ("number", "work")}
    relative_closure = {
        arm: abs(runs["fine"][arm]["work"]/runs["fine"][arm]["energy"]-1)
        for arm in ARMS}
    close = max(relative_closure.values()) < .002
    status = "unresolved"
    if close and all(m["large"] for m in metrics.values()):
        status = "large_difference"
    if close and all(m["small"] for m in metrics.values()):
        status = "small_difference"
    result = dict(status=status, metrics=metrics,
                  relative_work_energy_residual=relative_closure,
                  bookkeeping=checks, diagnostics=diagnostics(runs),
                  spectrum_edges=[value for value in BINS[:-1]]+["infinity"])
    write_json(path / "decision.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
