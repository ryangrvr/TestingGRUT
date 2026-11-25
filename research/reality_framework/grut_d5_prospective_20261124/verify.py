"""Verify frozen input, replay the production ledger and cross-check the independent graph."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
import math

from model import pulse_result, seed
from run import FREEZE_COMMIT, ROOT, run


def numeric_differences(first, second, tolerance, path="results"):
    if type(first) is not type(second):
        raise AssertionError(f"Type mismatch at {path}")
    if isinstance(first, dict):
        assert first.keys() == second.keys(), path
        for key in first:
            yield from numeric_differences(first[key], second[key], tolerance, f"{path}/{key}")
    elif isinstance(first, list):
        assert len(first) == len(second), path
        for i, (a, b) in enumerate(zip(first, second)):
            yield from numeric_differences(a, b, tolerance, f"{path}/{i}")
    elif type(first) in (int, float):
        assert math.isclose(first, second, rel_tol=0, abs_tol=tolerance), path
        yield abs(first - second)
    else:
        assert first == second, path


def main() -> None:
    manifest = json.loads((ROOT / "FREEZE.json").read_text())
    assert manifest["protocol_commit"] == FREEZE_COMMIT
    for name, digest in manifest["sha256"].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    protocol = json.loads((ROOT / "PROTOCOL.json").read_text())
    tolerance = protocol["thresholds"]["independent_oracle_absolute_tolerance"]
    saved = json.loads((ROOT / "RESULTS.json").read_text())
    differences = list(numeric_differences(saved, run(), tolerance))
    physical = json.loads((ROOT / "INDEPENDENT_CHECK.json").read_text())
    assert physical["all_oracle_checks_pass"]
    comparisons = []
    for case in physical["external_physical_graph_cases"]:
        stage = case["trial"] + 1
        p = Fraction(protocol["external_trials"]["probabilities_p"][case["trial"]])
        reference = pulse_result(
            seed(case["n"], case["h"], case["g"]), case["length"], stage, p,
            null=case["pilot_null"],
        )
        comparisons.append(abs(reference["joint_completion_and_archive_probability"]
                               - case["direct_target_and_archive_success"]))
        comparisons.append(abs(reference["correct_stage_probability"]
                               - case["direct_stage_probability"]))
        assert reference["orthogonal_endpoint_work_states"] == case["complete_endpoint_groups"]
    assert max(comparisons) < tolerance
    assert any(case["direct_stage_probability"] < 0.12 and
               case["direct_target_and_archive_success"] < 0.56
               for case in physical["external_physical_graph_cases"])
    print(json.dumps({
        "freeze_hashes_match": True,
        "numerical_fields_replayed": len(differences),
        "max_numeric_replay_difference": max(differences, default=0),
        "physical_graph_cases_cross_checked": len(comparisons) // 2,
        "max_physical_graph_probability_difference": max(comparisons),
        "small_clock_negative_diagnostic_retained": True,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
