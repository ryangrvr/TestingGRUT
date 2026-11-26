"""Compare the two frozen realizations of the protocol without dropping any arm."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path


FIELDS = (
    "acquisition_match",
    "stored_match",
    "actuator_match",
    "bond_probability",
    "sensor_fuel_consumed",
    "repair_fuel_consumed",
    "heat",
    "entropy_full",
    "entropy_device",
    "entropy_device_separate",
    "entropy_production_full",
    "source_device_mi",
    "total_device_correlation",
    "preparation_work",
    "joint_reset_work",
    "separate_reset_work",
)


def read_json(path):
    return json.loads(path.read_text())


def index_rows(data):
    rows = {row["name"]: row for row in data["rows"]}
    if len(rows) != len(data["rows"]):
        raise ValueError("Duplicate rows")
    return rows


def compare(builder, critic, parameters):
    left = index_rows(builder)
    right = index_rows(critic)
    declared_count = parameters["production"]["parameter_rows"]
    if left.keys() != right.keys() or len(left) != declared_count:
        raise AssertionError("Incomplete registered coverage")
    arms = set(parameters["arms"])
    errors_by_field = defaultdict(float)
    worst = {"error": -1.0}
    count = 0

    def check(label, a, b, row_name, arm=None):
        nonlocal count, worst
        if not (math.isfinite(a) and math.isfinite(b)):
            raise AssertionError(f"Non-finite input for {row_name} {arm} {label}")
        delta = abs(a - b)
        errors_by_field[label] = max(errors_by_field[label], delta)
        if delta > worst["error"]:
            worst = {"field": label, "row": row_name, "arm": arm, "error": delta}
        count += 1

    for name, row in left.items():
        independent = right[name]
        if row["parameters"] != independent["parameters"]:
            raise AssertionError(f"Parameter mismatch for {name}")
        if row["arms"].keys() != arms or independent["arms"].keys() != arms:
            raise AssertionError(f"Incomplete arm coverage for {name}")
        for key in ("match", "mismatch"):
            check("yield_" + key, row["yields"][key], independent["yields"][key], name)
        for arm in parameters["arms"]:
            for key in FIELDS:
                check(key, row["arms"][arm][key], independent["arms"][arm][key], name, arm)

    tolerance = parameters["gates"]["solver_agreement_tol"]
    return {
        "registered_rows": len(left),
        "arms_per_row": len(arms),
        "scalar_comparisons": count,
        "tolerance": tolerance,
        "max_differences_by_field": dict(errors_by_field),
        "worst_disagreement": worst,
        "independent_method": critic["method"],
        "independent_audit": critic["audit"],
        "passed": worst["error"] <= tolerance,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("builder", type=Path)
    parser.add_argument("critic", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    parameters = read_json(Path(__file__).with_name("PARAMETERS.json"))
    result = compare(read_json(args.builder), read_json(args.critic), parameters)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"{result['scalar_comparisons']} comparisons, "
          f"maximum error {result['worst_disagreement']['error']:.3g}")
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
