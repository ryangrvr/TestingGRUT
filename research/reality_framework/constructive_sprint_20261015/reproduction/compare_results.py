"""Compare new immutable-wrapper execution with published scientific outputs.

Residuals at floating-point noise scale are recorded, never required to match
bitwise. Physical probabilities and enclosure endpoints receive stated mixed
tolerances; dimensional integers and target counts must agree exactly.
"""
from pathlib import Path
import json
import math
import platform
import sys
import numpy as np
import scipy

P = Path(__file__).resolve().parent
PACKET = P.parents[1] / "cp12_charged_accessibility_gate"

def main():
    expected = json.loads((PACKET / "NUMERICAL_AUDIT.json").read_text())
    actual = json.loads((P / "python.stdout.log").read_text())
    node = json.loads((P / "node.stdout.log").read_text())
    checks = []
    def record(name, a, b, atol=0., rtol=0.):
        passed = (a == b) if not (atol or rtol) else math.isclose(float(a), float(b), abs_tol=atol, rel_tol=rtol)
        checks.append(dict(name=name, actual=a, published=b, atol=atol, rtol=rtol, passed=passed))
    for kind in ["node", "python"]:
        metadata = json.loads((P / f"{kind}.execution.json").read_text())
        record(f"execution.{kind}.exit_code", metadata["exit_code"], 0)
        record(f"execution.{kind}.cache", metadata["CP12_FORMATION_CACHE"], "unset")
        record(f"execution.{kind}.threads", metadata["threads"], 1)
    record("node.exact_core_status", node["status"], "CP12_CHARGED_ACCESS_EXACT_CORE_PASS")
    bounds = json.loads((P / "ANALYTIC_BOUNDS_RERUN.json").read_text())
    bounds_expected = json.loads((PACKET / "ANALYTIC_BOUNDS.json").read_text())
    bounds_differences = []
    def compare_bounds(a, b, path=""):
        if isinstance(a, dict) and isinstance(b, dict):
            if a.keys() != b.keys():
                bounds_differences.append(dict(path=path, passed=False, reason="keys"))
            for key in a.keys() & b.keys():
                compare_bounds(a[key], b[key], path+"."+key)
        elif isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                bounds_differences.append(dict(path=path, passed=False, reason="length"))
            for i, (x, y) in enumerate(zip(a, b)):
                compare_bounds(x, y, path+f"[{i}]")
        elif a != b:
            numeric = isinstance(a, (float, int)) and isinstance(b, (float, int))
            passed = numeric and math.isclose(a, b, rel_tol=2e-13, abs_tol=0.)
            bounds_differences.append(dict(path=path, actual=a, published=b, passed=passed,
                relative_difference=abs(a-b)/max(abs(a), abs(b)) if numeric and max(abs(a), abs(b)) else None))
    compare_bounds(bounds, bounds_expected)
    record("analytic_bounds.complete_saved_table_with_float_tolerance", all(x["passed"] for x in bounds_differences), True)
    onset = json.loads((P / "PYTHON_ONSET_CHECK.json").read_text())
    record("python.onset.status", onset["status"], "PASS")
    record("python.onset.signed_amplitude", onset["exact_signed_amplitude"], "-49583226880")
    record("python.onset.minimum_hops", onset["minimum_hops"], 22)
    for section, keys in {
        "spectral": ["N", "Q", "budget", "retained_dimension", "outside_orbits", "target_orbits", "trial_modes", "coefficient_pruning"],
        "capacity": ["N", "Q", "budget", "retained_dimension", "outside_orbits", "nnz", "boundary_nnz", "target_orbits_by_vertex_fermions", "time_steps"]
    }.items():
        for key in keys:
            record(f"{section}.{key}", actual[section][key], expected[section][key])
    for key in ["absolute_weighted_target_Gram_sum", "projected_infinite_time_average_non_degenerate"]:
        record("spectral."+key, actual["spectral"][key], expected["spectral"][key], atol=2e-11, rtol=2e-10)
    for section in ["spectral", "capacity"]:
        for i, (a, b) in enumerate(zip(actual[section]["rows"], expected[section]["rows"], strict=True)):
            record(f"{section}.rows[{i}].T", a["T"], b["T"])
            key = "p" if section == "spectral" else "projected_probability"
            record(f"{section}.rows[{i}].{key}", a[key], b[key], atol=2e-12 if section == "spectral" else 2e-24, rtol=2e-8)
            for j in range(2):
                record(f"{section}.rows[{i}].full_probability_interval[{j}]", a["full_probability_interval"][j], b["full_probability_interval"][j], atol=2e-11 if section == "spectral" else 2e-18, rtol=2e-7)
            if section == "capacity":
                record(f"{section}.rows[{i}].full_state_error_bound", a["full_state_error_bound"], b["full_state_error_bound"], atol=2e-15, rtol=2e-7)
    residual_keys = ["initial_reconstruction_norm", "full_H_weighted_trial_eigenvector_residual", "boundary_weighted_residual", "internal_weighted_residual"]
    residuals = {key: dict(actual=actual["spectral"][key], published=expected["spectral"][key]) for key in residual_keys}
    for key in residual_keys:
        # Eigensolver-dependent roundoff is assessed by absolute useful scale.
        checks.append(dict(name="spectral.roundoff_scale."+key, actual=actual["spectral"][key], maximum_allowed=1e-9, passed=actual["spectral"][key] < 1e-9))
    outcome = dict(status="PASS" if all(c["passed"] for c in checks) else "FAIL",
        executed_commit="1801afee27096c03c2a3b64cc811ac0be0d2a4a7", cache_used=False,
        node_exact_status=node["status"], checks=checks,
        analytic_table_float_rtol=2e-13,
        analytic_table_nonbitwise_entries=bounds_differences,
        residuals_not_required_to_match_bitwise=residuals,
        numeric_runtime_seconds=actual.get("seconds"),
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()),
        limitation="Floating-point reproduction of analytical residual bounds; no outward-rounded interval arithmetic or external mathematical certification.")
    (P / "COMPARISON.json").write_text(json.dumps(outcome, indent=2)+"\n")
    print(json.dumps({"status": outcome["status"], "checks": len(checks), "failed": [c for c in checks if not c["passed"]]}, indent=2))
    raise SystemExit(0 if outcome["status"] == "PASS" else 1)

if __name__ == "__main__":
    main()
