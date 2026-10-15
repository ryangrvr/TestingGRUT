"""Validate saved checks; --rerun reproduces both physical calculations first."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import os
import subprocess
import sys

P = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true")
    args = parser.parse_args()
    if args.rerun:
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        for script in ["run.py", "exact_onset.py", "neutral_record.py", "report.py"]:
            subprocess.run([sys.executable, str(P / script)], check=True, env=env)
    data = json.loads((P / "RESULTS.json").read_text())
    exact = json.loads((P / "EXACT_ONSET.json").read_text())
    assert data["status"] == exact["status"] == "PASS"
    assert [b["cutoff"] for b in data["boxes"]] == [6,10]
    assert len(data["boxes"]) == 2
    maximum_phase_error = maximum_H_error = 0.
    for box in data["boxes"]:
        maximum_H_error = max(maximum_H_error, box["independent_H_column_maximum_error"])
        assert box["outgoing_dimension"] > 0
        assert box["independent_H_column_maximum_error"] < 1e-12
        for row in box["results"]:
            for name in ["host", "protection", "original_circulation"]:
                lo, hi = row[name]["interval"]
                assert 0 <= lo <= hi <= 1
            for d in row["detector_responses"]:
                assert d["zero_excluded"]
                assert d["interval"][0] * d["interval"][1] > 0
                maximum_phase_error = max(maximum_phase_error, d["exact_phase_shift_formula_difference"])
                assert d["exact_phase_shift_formula_difference"] < 1e-12
    for check in data["cutoff_comparisons"]:
        assert check["maximum_detector_difference"] <= check["analytic_allowance"] + check["floating_point_tolerance"]
    assert exact["zero_edge_target_configurations"] == 93
    assert exact["maximum_absolute_target_field"] == 2
    for d in exact["minimum_field_distances"]:
        assert d["minimum_electric_L1_distance_to_zero_edge"] == 8*d["flux"]-5
    for d in exact["exact_unregulated_onsets"]:
        assert d["first_zero_edge_H_power"] == 8*d["flux"]-5
        assert d["first_changed_graph_H_power"] == 8*d["flux"]-4
        assert Fraction(d["graph_change_probability_coefficient"]) > 0
    neutral_checks=None
    if (P / "NEUTRAL_RECORD.json").exists():
        neutral=json.loads((P / "NEUTRAL_RECORD.json").read_text())
        assert neutral["status"] == "PASS"
        assert neutral["model_checks"]["quotient_dimension"] == 2242
        assert neutral["model_checks"]["physical_dimension"] == 13413
        assert neutral["model_checks"]["isometry_error"] < 1e-11
        assert neutral["model_checks"]["H_intertwiner_error"] < 1e-10
        assert abs(sum(x["probability"] for x in neutral["complete_record_outcome_probabilities"])-1)<1e-10
        for name in ["host_survival","circulation_survival","joint_host_and_circulation_survival"]:
            lo,hi=neutral[name]["conditional_interval"]
            assert 0 <= lo <= hi <= 1
        assert neutral["density_response"]["zero_excluded"]
        assert neutral["same_host_joint_density_response"]["zero_excluded"]
        neutral_checks=dict(event_probability=neutral["event_probability_at_wait"],
            density_response_zero_excluded=True,complete_channel_normalized=True,
            bounded_record_energy_error_observable=True)
    sha = {}
    for path in [P.parents[1]/"cp12_ground_connectivity"/"frozen_model.py",
                 P.parents[1]/"cp12_ground_connectivity"/"exact_destination.py",
                 P.parents[1]/"cp12_distance_history"/"evolution.py"]:
        sha[str(path.relative_to(P.parents[1]))] = hashlib.sha256(path.read_bytes()).hexdigest()
    out = dict(status="PASS", fresh_full_numerical_reproduction=args.rerun,
        scope="Saved-result consistency checks. --rerun first reruns every finite-time and exact-onset calculation.",
        maximum_independent_H_column_error=maximum_H_error,
        maximum_phase_identity_discrepancy=maximum_phase_error,
        maximum_electric_cutoff_response_discrepancy=max(x["maximum_detector_difference"] for x in data["cutoff_comparisons"]),
        neutral_record_checks=neutral_checks,
        source_sha256=sha)
    (P / "VERIFICATION.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
