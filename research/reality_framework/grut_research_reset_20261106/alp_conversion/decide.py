"""Evaluate the frozen conversion gate without importing detector data."""
from dataclasses import asdict
import json

import numpy as np

from model import CASES, HERE, matching, write_json


def main():
    directory = HERE / "results"
    if not json.loads((directory / "qualification.json").read_text())["passed"]:
        raise ValueError("The conversion ladder did not qualify.")
    rows = {}
    for name, config in CASES.items():
        row = json.loads((directory / (name+".json")).read_text())
        if row["config"] != asdict(config):
            raise ValueError("Configuration mismatch for "+name)
        with np.load(directory / (name+".npz")) as raw:
            late = raw["late"]
            spectra = raw["late_quadratic_spectrum"]
            if not np.allclose(spectra.sum(axis=1), late[:, 1], rtol=1e-12):
                raise ValueError("Spectrum and quadratic proxy disagree for "+name)
            if not np.all(np.diff(raw["ledger"][:, 0]) > 0):
                raise ValueError("Energy ledger has a nonmonotone scale factor.")
            row["quadratic_zero_mode_fraction"] = float(np.mean(late[:, 3]))
            row["quadratic_rms_p_over_m_final"] = float(late[-1, 4])
            row["quadratic_spectrum_over_charge"] = np.mean(spectra, axis=0).tolist()
            row["comoving_bin_edges"] = raw["comoving_edges"].tolist()
        rows[name] = row

    base = rows["base"]["late_quadratic_proxy"]
    differences = {name: row["late_quadratic_proxy"]/base-1 for name, row in rows.items()
                   if name != "base"}
    controlled_cases = ("base", "time", "spacing", "volume")
    refinement_pass = (all(rows[name]["particle_gate"] for name in controlled_cases)
                       and max(abs(differences[name]) for name in controlled_cases[1:]) < .02)
    physical = json.loads((directory / "physical_linear.json").read_text())
    onset = physical["fine"]
    inferred_hierarchy = physical["scales"]["m_over_H_star"]
    result = dict(
        reduced_hierarchy_particle_gate=bool(refinement_pass),
        physical_C=None,
        physical_coupling_band=None,
        prospective_status="unopened",
        hierarchy_gap_vs_slowest_run=inferred_hierarchy*CASES["slower"].h,
        onset_before_nonlinear_start={label: row["a_onset"] < CASES["base"].a_start
                                      for label, row in onset.items()},
        measured_proxy_shifts=differences,
        cases=rows,
        conditional_inference_examples=[matching(c) for c in (.5, 1., 2., 3.)])
    write_json(directory / "conversion_decision.json", result)
    print(json.dumps({key: value for key, value in result.items() if key != "cases"}, indent=2))


if __name__ == "__main__":
    main()
