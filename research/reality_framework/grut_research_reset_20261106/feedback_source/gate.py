"""Apply the registered prescribed-force and live-feedback acceptance rules."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from source import CASES, HERE, write_json

LIMIT_FORCE = .03
LIMIT_ENERGY = .002
LIVE_CASES = ("base", "time", "fine", "tapered")


def load(directory, name):
    row = json.loads((directory / (name+".json")).read_text())
    with np.load(directory / (name+".npz")) as archive:
        trace = dict(archive)
    return row, trace


def curve_difference(first, second, key):
    t, s = first["time"], second["time"]
    if np.array_equal(t, s):
        return float(np.max(np.abs(first[key]-second[key])))
    # Use the denser samples inside the interval spanned by both grids.
    if len(t) > len(s):
        return curve_difference(second, first, key)
    keep = (s >= t[0]) & (s <= t[-1])
    interpolated = CubicSpline(t, first[key])(s[keep])
    return float(np.max(np.abs(interpolated-second[key][keep])))


def prescribed(directory, cases=CASES):
    rows, traces = {}, {}
    for case, config in cases.items():
        rows[case], traces[case] = load(directory, case)
        if rows[case]["config"] != vars(config):
            raise ValueError("Unexpected prescribed run configuration for "+case)
    fine = traces["fine"]
    peak = rows["fine"]["peak_force"]["local_match"]
    relative = {case: curve_difference(traces["base"], traces[case], "local_match")/peak
                for case in ("time", "spacing", "longitudinal", "transverse", "fine")}
    relative["tapered"] = curve_difference(fine, traces["tapered"], "local_match")/peak
    relative["linear_replacement"] = float(
        np.max(np.abs(fine["local_match"]-fine["linear_match"]))/peak)
    raw_peaks = {name: rows["fine"]["peak_force"][name]
                 for name in ("old", "raw", "local_match", "linear_match")}
    diagnostic_shifts = {
        name: curve_difference(fine, traces["tapered"], name)
        for name in ("old", "raw", "local_match", "linear_match", "nonlinear")}
    closure = rows["fine"]["work_relative_error"]["local_match"]
    result = dict(passed=max(relative.values()) < LIMIT_FORCE and closure < LIMIT_ENERGY,
                  force_peak=peak, relative_shifts=relative,
                  fine_peaks=raw_peaks, boundary_absolute_shifts=diagnostic_shifts,
                  transverse_linear_replacement_shift=curve_difference(
                      traces["base"], traces["transverse"], "linear_match")/peak,
                  transverse_relative_work_shift=abs(rows["transverse"]["work"]["local_match"] /
                                                     rows["base"]["work"]["local_match"]-1),
                  fine_work=rows["fine"]["work"],
                  fine_energy=rows["fine"]["outgoing"]["energy"],
                  fine_number=rows["fine"]["outgoing"]["number"],
                  fine_work_error=closure)
    write_json(directory / "prescribed_gate.json", result)
    return result


def live(directory, cases=CASES):
    rows, traces = {}, {}
    for case in LIVE_CASES:
        rows[case], traces[case] = load(directory, "live_"+case)
        if rows[case]["config"] != dict(vars(cases[case]), wave=1.3):
            raise ValueError("Unexpected live run configuration for "+case)
    fine = traces["fine"]
    signal = rows["fine"]["amplitude_departure"]
    shifts = {case: curve_difference(traces[case], fine, "amplitude")
              for case in ("base", "time", "tapered")}
    allowance = shifts["base"]+shifts["tapered"]
    closure = {case: rows[case]["relative_energy_defect"] for case in LIVE_CASES}
    result = dict(passed=max(closure["fine"], closure["tapered"]) < LIMIT_ENERGY
                  and signal > 3*allowance,
                  signal=signal, trajectory_shifts=shifts, allowance=allowance,
                  signal_over_allowance=signal/allowance, energy_defects=closure,
                  fine_final=rows["fine"]["final"],
                  fine_peak_total_energy=rows["fine"]["peak_total_energy"])
    write_json(directory / "live_gate.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gate", choices=("prescribed", "live"))
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    result = (prescribed if args.gate == "prescribed" else live)(args.output)
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("Registered gate missed.")


if __name__ == "__main__":
    main()
