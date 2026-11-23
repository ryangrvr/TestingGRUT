"""Replay saved physics and reject changes to frozen inputs/negative verdicts.

This verifies a numerical implementation, not physical novelty or empirical
truth. independent_check.py supplies the separate small full-matrix checks.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import sys
from run import run

ROOT = Path(__file__).resolve().parent


def compare(saved, fresh, path="", discrepancies=None):
    if discrepancies is None:
        discrepancies = []
    if isinstance(saved, dict) and isinstance(fresh, dict):
        if saved.keys() != fresh.keys():
            raise AssertionError(f"key mismatch at {path}")
        for key in saved:
            if path == "" and key == "runtime":
                continue
            compare(saved[key], fresh[key], f"{path}/{key}", discrepancies)
    elif isinstance(saved, list) and isinstance(fresh, list):
        if len(saved) != len(fresh):
            raise AssertionError(f"length mismatch at {path}")
        for i, (a, b) in enumerate(zip(saved, fresh)):
            compare(a, b, f"{path}/{i}", discrepancies)
    elif isinstance(saved, bool) or isinstance(fresh, bool):
        if saved is not fresh:
            raise AssertionError(f"boolean mismatch at {path}")
    elif isinstance(saved, (int, float)) and isinstance(fresh, (int, float)):
        if not math.isfinite(saved) or not math.isfinite(fresh):
            raise AssertionError(f"nonfinite result at {path}")
        difference = abs(saved - fresh)
        discrepancies.append(difference)
        if difference > 1e-10 * max(1, abs(fresh)):
            raise AssertionError(f"numerical mismatch at {path}: {difference}")
    elif saved != fresh:
        raise AssertionError(f"value mismatch at {path}")
    return discrepancies


def validate_science(result):
    assert all(v == "not_established" for v in result["negative_results"].values())
    assert result["maximally_mixed_active_control"]["frozen_no_new_record_forecast"] == "failed"
    assert len(result["cases"]) == 5
    for case in result["cases"]:
        cert = case["certificate"]
        assert cert["all_time_order_gate_pass"]
        assert cert["all_time_count_record_gate_pass"]
        assert cert["kick_recovery_gate_pass"]
        assert case["sampled"]["minimum_order_probability"] >= cert["order_lower_display"]
        for row in case["time_rows"]:
            assert 0 <= row["order_probability"] <= 1
            assert row["order_count_identity_error"] < 1e-10
            assert row["joint_normalization_error"] < 1e-10
            assert min(v for r in row["joint_K_Nbath"] for v in r) >= -1e-13
        assert case["disturbance"]["order_probability_union_lower"] >= cert["kick_recovery_lower_display"]


def verify(saved_path):
    freeze = json.loads((ROOT / "FREEZE.json").read_text())
    for name, expected in freeze["files"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        assert actual == expected, f"frozen file changed: {name}"
    saved = json.loads(saved_path.read_text())
    fresh = run()
    differences = compare(saved, fresh)
    validate_science(saved)
    hostile = []
    for kind in ["forged_physical_success", "corrupted_probability", "erased_control_failure"]:
        altered = deepcopy(saved)
        if kind == "forged_physical_success":
            altered["negative_results"]["new_GRUT_physics"] = "established"
        elif kind == "corrupted_probability":
            altered["cases"][0]["time_rows"][0]["order_probability"] = 1.1
        else:
            altered["maximally_mixed_active_control"]["frozen_no_new_record_forecast"] = "passed"
        try:
            compare(altered, fresh)
            validate_science(altered)
        except AssertionError:
            hostile.append({"case": kind, "rejected": True})
        else:
            raise AssertionError(f"failed to reject {kind}")
    return {"freeze_hashes_match": True, "physics_replay_pass": True,
            "replayed_time_rows": sum(len(c["time_rows"]) for c in saved["cases"]),
            "numerical_values_compared": len(differences),
            "max_absolute_replay_difference": max(differences, default=0),
            "hostile_controls": hostile,
            "scope": "Same-code replay and analytic gate arithmetic, not an independent physical replication or experiment."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "RESULTS.json")
    parser.add_argument("--output", type=Path, default=ROOT / "VERIFICATION.json")
    args = parser.parse_args()
    result = verify(args.results)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
