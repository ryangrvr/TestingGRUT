"""Isolate the static and slope errors of the first time-split tangent."""
import json

import numpy as np

from gate import curve_difference, load
from source import CASES, HERE, System, write_json


def main():
    directory = HERE / "results"
    shifts, traces = {}, {}
    for name, config in CASES.items():
        _, trace = load(directory, name)
        system = System(config)
        dt = config.dt
        angle = system.pair_frequency*dt/2
        strength = system.pair_strength
        dt_value = system.contact-np.sum(strength*(dt/2)/np.tan(angle))
        dt_slope = -np.sum(strength*dt**3*np.cos(angle)/(8*np.sin(angle)**3))
        dv = dt_value-system.match["cutoff_static"]
        ds = dt_slope-system.match["cutoff_omega2_slope"]
        trace["corrected"] = trace["local_match"]-dv*trace["amplitude"]+ds*trace["acceleration"]
        traces[name] = trace
        peak = np.max(np.abs(trace["corrected"]))
        shifts[name] = dict(static_error=dv, slope_error=ds,
                            corrected_force_peak=peak, remaining_method_difference=
                            np.max(np.abs(trace["corrected"]-trace["linear_match"]))/peak)

    fine = traces["fine"]
    peak = shifts["fine"]["corrected_force_peak"]
    result = dict(tangent_shifts=shifts, residual_refinements={
        name: curve_difference(traces["base"], traces[name], "corrected")/peak
        for name in ("time", "spacing", "longitudinal", "transverse", "fine")})
    result["residual_refinements"]["tapered"] = curve_difference(
        fine, traces["tapered"], "corrected")/peak
    write_json(directory / "split_diagnosis.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
