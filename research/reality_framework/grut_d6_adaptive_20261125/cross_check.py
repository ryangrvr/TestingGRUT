#!/usr/bin/env python3
"""Compare physical propagation with the frozen D6 reduced calculation."""

import argparse
import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROW_KEYS = ("n", "strength", "port", "intervention", "source", "challenge")
METRICS = (
    ("a_only", "a_only_probability"),
    ("a_only_target", "a_only_target"),
    ("endpoint", "endpoint_probability"),
    ("endpoint_memory_zero", "endpoint_memory_zero"),
    ("endpoint_memory_correct", "endpoint_memory_correct"),
    ("dump_correct", "preintervention_dump_correct"),
    ("dump_correct", "preintervention_record_correct"),
    ("endpoint_a_target", "endpoint_A_success"),
    ("endpoint_b_target", "endpoint_B_success"),
    ("endpoint_joint_target", "endpoint_both_success"),
    ("endpoint_witness_zero", "endpoint_W_zero"),
    ("endpoint_forecast", "endpoint_forecast"),
    ("endpoint_joint_with_forecast", "endpoint_both_with_forecast"),
)


def identity(row):
    return tuple(row[k] for k in ROW_KEYS)


def load(path):
    return json.loads(path.read_text())


def compare(reduced_path, oracle_path):
    protocol = load(HERE / "PROTOCOL.json")
    reduced = load(reduced_path)
    oracle = load(oracle_path)
    reduced_rows = [dict(zip(reduced["columns"], r, strict=True))
                    for r in reduced["conditional_rows"]]
    if not oracle["all_oracle_checks_pass"]:
        raise ValueError("The oracle did not pass its own physical checks")
    reduced_by_key = {identity(r): r for r in reduced_rows}
    required = protocol["independent_oracle_required_cases"]
    expected = {
        (required["n"], *arm, b, z)
        for arm, b, z in itertools.product(required["arms_order"], required["source"],
                                           required["stage_challenge"])
    }
    oracle_rows = oracle["rows"]
    if len(oracle_rows) != len(expected) or {identity(r) for r in oracle_rows} != expected:
        raise ValueError("The oracle rows differ from the frozen grid")
    errors = []
    for row in oracle_rows:
        key = identity(row)
        n, strength, port, intervention, _, z = key
        main_row = reduced_by_key[key]
        for left, right in METRICS:
            errors.append(abs(main_row[left] - row[right]))
        for side, flank in (("left", "late_minus"), ("right", "late_plus")):
            errors.append(abs(main_row[f"late_{side}_joint_target"] - row[flank]["both_success"]))
            errors.append(abs(main_row[f"late_{side}_joint_with_forecast"]
                              - row[flank]["both_with_forecast"]))
        source_zero = reduced_by_key[n, strength, port, intervention, 0, z]
        source_one = reduced_by_key[n, strength, port, intervention, 1, z]
        memory_distance = abs(source_zero["memory_zero_at_cut"] - source_one["memory_zero_at_cut"])
        # The main port carries only M. The passive bit is independent, and
        # the open direct port carries a perfect additional source label.
        response_distance = 1.0 if port == "direct" and strength != "sealed" else memory_distance
        errors.append(abs(memory_distance - row["memory_cut_trace_distance"]))
        errors.append(abs(response_distance - row["complete_response_cut_trace_distance"]))

    stages = {c["n"]: c["a_only"] for c in reduced["clock"]}
    for stage in oracle["stage_diagnostics"]:
        errors.append(abs(stages[stage["n"]] - stage["A_only_probability"]))
    tolerance = protocol["numerical_gates"]["oracle_abs_error"]
    residual = max(errors)
    if residual > tolerance:
        raise ValueError(f"Reduced and physical propagation disagree by {residual}")
    return {"rows": len(oracle_rows), "scalar_checks": len(errors), "max_abs_error": residual}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reduced", type=Path, default=HERE / "RESULTS.json")
    parser.add_argument("--oracle", type=Path, default=HERE / "ORACLE_RESULTS.json")
    args = parser.parse_args()
    print(json.dumps(compare(args.reduced, args.oracle), indent=2))


if __name__ == "__main__":
    main()
