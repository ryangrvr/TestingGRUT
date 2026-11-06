"""Apply the bounded route gate without reading detector data."""
from dataclasses import asdict
import json

import numpy as np

from match import CASES, DIRECTORY
from model import write_json


def decide():
    matched = json.loads((DIRECTORY / "match.json").read_text())
    rows, misses = {}, {}
    for name, config in CASES.items():
        path = DIRECTORY / (name+".json")
        if not path.exists():
            misses[name] = ["missing_endpoint"]
            continue
        row = json.loads(path.read_text())
        if row["config"] != asdict(config):
            raise ValueError("Unexpected lattice configuration for "+name)
        rows[name] = row
        with np.load(DIRECTORY / (name+".npz")) as raw:
            if not all(np.isfinite(raw[k]).all() for k in raw.files):
                raise ValueError("Non-finite trace for "+name)
            if not np.all(np.diff(raw["ledger"][:, 0]) > 0):
                raise ValueError("Non-monotone scale factor for "+name)
            if not np.allclose(raw["late_quadratic_spectrum"].sum(axis=1), raw["late"][:, 1],
                               rtol=1e-12):
                raise ValueError("Inconsistent spectrum for "+name)
        if not row["completed"] or "terminal_window" not in row:
            misses[name] = ["missing_endpoint"]
            continue
        limits = dict(drift=.02, dictionary=.02, anharmonic=.02,
                      outside_branch=1e-4, cutoff_tail=.01)
        late = row["terminal_window"]
        misses[name] = [key for key, limit in limits.items() if late[key] >= limit]
        if not late["fixed_branch"]:
            misses[name].append("rotating_mean")
        if row["energy_residual"] >= .01:
            misses[name].append("energy")
    reference = rows.get("base", {})
    base = reference.get("quadratic_proxy") if reference.get("completed") else None
    variation = {name: row["quadratic_proxy"]/base-1 for name, row in rows.items()
                 if base and row["completed"] and "quadratic_proxy" in row and name != "base"}
    pass_route = matched["passed"] and not any(misses.values()) and (
        len(variation) == len(CASES)-1 and max(map(abs, variation.values())) < .05)
    result = dict(route_passed=pass_route,
                  decision="continue_ALP_redshift_test" if pass_route else "close_ALP_lane",
                  physical_C=None, prospective_search="unopened",
                  linear_match_passed=matched["passed"], cases=rows, misses=misses,
                  quadratic_proxy_shifts=variation)
    write_json(DIRECTORY / "route_decision.json", result)
    return result


if __name__ == "__main__":
    print(json.dumps(decide(), indent=2))
