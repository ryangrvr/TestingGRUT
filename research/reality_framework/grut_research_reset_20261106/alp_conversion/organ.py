"""Retrospective check using the pinned public tabulation of ORGAN Fig. 5."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model import HERE, MASS_EV, matching, write_json

DATA_HASH = "d9f34132793d1fcd49bffa2642a43d8aa64f607adda978defd98b22662e1f4ea"
DATA_URL = ("https://github.com/cajohare/AxionLimits/blob/"
            "7d375f4879b32406a239fe48d2615a4bfd9bc0bb/"
            "limit_data/AxionPhoton/ORGAN-1b.txt")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    args = parser.parse_args()
    directory = HERE / "results"
    conversion = json.loads((directory / "conversion_decision.json").read_text())
    if conversion["prospective_status"] != "unopened":
        raise ValueError("Expected a frozen conversion decision before this retrospective audit.")
    digest = hashlib.sha256(args.data.read_bytes()).hexdigest()
    if digest != DATA_HASH:
        raise ValueError("Unexpected ORGAN curve digest.")

    data = np.loadtxt(args.data)
    curve = data[data[:, 1] < 1e-6]  # Drop polygon closure used for plotting.
    if not np.all(np.diff(curve[:, 0]) > 0) or not curve[0, 0] < MASS_EV < curve[-1, 0]:
        raise ValueError("The archived limit curve does not bracket the registered mass.")
    i = np.searchsorted(curve[:, 0], MASS_EV)
    limit = float(np.interp(MASS_EV, curve[:, 0], curve[:, 1]))
    limit_log = float(np.exp(np.interp(np.log(MASS_EV), np.log(curve[:, 0]), np.log(curve[:, 1]))))
    expected = matching(2.)["g_agamma_per_gev"]
    crossing = 2*(limit/expected)**2
    result = dict(
        data_url=DATA_URL, data_sha256=digest, source="AxionLimits tabulation of ORGAN Phase 1b Figure 5",
        experiment_doi="10.1103/PhysRevLett.132.031601", nominal_confidence=.95,
        mass_eV=MASS_EV, frequency_GHz=MASS_EV*241798.9242,
        neighboring_rows_eV_and_per_GeV=curve[i-1:i+1].tolist(),
        nominal_limit_per_GeV=limit, log_interpolation_limit_per_GeV=limit_log,
        nominal_conversion_crossing=crossing,
        C2_coupling_over_curve=expected/limit,
        C2_power_over_curve=(expected/limit)**2,
        experiment_reported_relative_coupling_systematic=.061,
        local_density_gev_cm3=.45,
        physical_verdict="unresolved", prospective_interval="unopened")
    write_json(directory / "organ_retrospective.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
