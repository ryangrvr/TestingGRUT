"""Replay the frozen D4 calculation and verify all archived fields."""

from __future__ import annotations

import hashlib
import json
import math

from run import FREEZE_COMMIT, ROOT, run


def compare(saved, current, tolerance, path="results"):
    if type(saved) is not type(current):
        raise AssertionError(f"Type mismatch at {path}: {type(saved)} and {type(current)}")
    if isinstance(saved, dict):
        if saved.keys() != current.keys():
            raise AssertionError(f"Different fields at {path}")
        for key in saved:
            yield from compare(saved[key], current[key], tolerance, f"{path}/{key}")
    elif isinstance(saved, list):
        if len(saved) != len(current):
            raise AssertionError(f"Different lengths at {path}")
        for index, (before, after) in enumerate(zip(saved, current)):
            yield from compare(before, after, tolerance, f"{path}/{index}")
    elif type(saved) in (int, float):
        if not math.isclose(saved, current, rel_tol=0, abs_tol=tolerance):
            raise AssertionError(f"Numerical mismatch at {path}: {saved} versus {current}")
        yield abs(saved - current)
    elif saved != current:
        raise AssertionError(f"Different value at {path}")


def main() -> None:
    freeze = json.loads((ROOT / "FREEZE.json").read_text())
    assert freeze["protocol_commit"] == FREEZE_COMMIT
    for name, expected in freeze["sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        assert actual == expected, f"Frozen file changed: {name}"
    protocol = json.loads((ROOT / "PROTOCOL.json").read_text())
    tolerance = protocol["verification"]["tolerance"]
    saved = json.loads((ROOT / "RESULTS.json").read_text())
    differences = list(compare(saved, run(), tolerance))
    samples = saved["samples_in_ready_comparison"]["four_fault_rows"]
    assert len(samples) == 243
    assert saved["no_fault"]["false_repair_bath_activation"] == 0
    assert saved["ambiguous_pair"]["equal_prior_historical_success_ceiling"] == 0.5
    assert all(case["exact_writer_on_historical_codeword_at_nominal_write"] == 0
               for case in saved["simultaneous_faults"])
    print(json.dumps({
        "frozen_files_match": True,
        "repeated_timing_rows_replayed": len(samples),
        "numerical_fields_replayed": len(differences),
        "maximum_numeric_difference": max(differences, default=0),
        "simultaneous_failure_and_state_collision_retained": True,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
