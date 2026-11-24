"""Check frozen inputs, archived production results and required negative controls."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from run import FREEZE_COMMIT, run

ROOT = Path(__file__).resolve().parent


def compare(saved, actual, tolerance, path="results"):
    """Compare every archived field, allowing roundoff in numerical values."""
    if type(saved) is not type(actual):
        raise AssertionError(f"Type mismatch at {path}")
    if isinstance(saved, dict):
        if saved.keys() != actual.keys():
            raise AssertionError(f"Field mismatch at {path}")
        for key, item in saved.items():
            yield from compare(item, actual[key], tolerance, f"{path}/{key}")
    elif isinstance(saved, list):
        if len(saved) != len(actual):
            raise AssertionError(f"Length mismatch at {path}")
        for index, (item, value) in enumerate(zip(saved, actual)):
            yield from compare(item, value, tolerance, f"{path}/{index}")
    elif type(saved) in (int, float):
        if not math.isclose(saved, actual, abs_tol=tolerance, rel_tol=0):
            raise AssertionError(f"Numerical mismatch at {path}: {saved} versus {actual}")
        yield abs(saved - actual)
    elif saved != actual:
        raise AssertionError(f"Value mismatch at {path}")


def main() -> None:
    freeze = json.loads((ROOT / "FREEZE.json").read_text())
    assert freeze["protocol_commit"] == FREEZE_COMMIT
    for name, promised in freeze["sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        assert actual == promised, f"Frozen input changed: {name}"

    protocol = json.loads((ROOT / "PROTOCOL.json").read_text())
    tolerance = protocol["verification"]["replay_tolerance"]
    saved = json.loads((ROOT / "RESULTS.json").read_text())
    computed = run()
    differences = list(compare(saved, computed, tolerance))
    assert saved["short_chain_control"]["joint_gate"] == "fails"
    assert saved["verdicts"]["new_physical_law"] == "unestablished"
    assert saved["verdicts"]["cosmological_selection"] == "unestablished"
    for row in saved["time_rows"]:
        intact = row["intact_restore_and_retain_probability"]
        assert row["scrambled_restore_and_retain_probability"] == intact / 2
        assert row["wrong_restore_and_retain_probability"] == 0
    print(json.dumps({
        "frozen_inputs_verified": True,
        "time_rows_replayed": len(saved["time_rows"]),
        "numerical_fields_replayed": len(differences),
        "maximum_absolute_replay_difference": max(differences, default=0),
        "negative_control_retained": True,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
