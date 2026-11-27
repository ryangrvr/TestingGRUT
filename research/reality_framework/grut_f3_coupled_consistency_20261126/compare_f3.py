"""Compare the independent rational checks for the frozen F3 packet."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path


def compare(production, critique):
    counts = Counter()

    def verify(group, a, b):
        if isinstance(a, list):
            if len(a) != len(b):
                raise AssertionError(f"Dimension mismatch in {group}")
            for left, right in zip(a, b):
                verify(group, left, right)
        else:
            if F(a) != F(b):
                raise AssertionError(f"{group}: expected {a}, received {b}")
            counts[group] += 1

    rows = production["background_rows"]
    independent_rows = critique["fixtures"]
    if len(rows) != len(independent_rows):
        raise AssertionError("Incomplete fixture coverage")
    for row, reference in zip(rows, independent_rows):
        for left, right in (
            ("x", "x"), ("rho", "rho_over_lambda4"), ("pressure", "p_over_lambda4"),
            ("px", "p_x"), ("rho_x", "rho_x"), ("epsilon_h", "epsilon_H"),
            ("sound_speed_sq", "cs_squared"), ("flow", "F"),
            ("matched_kinetic", "canonical_kinetic_over_lambda4"),
            ("matched_potential", "canonical_potential_over_lambda4"),
        ):
            verify("background", row[left], reference[right])
        for mode in ("active", "canonical", "metric"):
            verify("mode_coefficients", row[f"pump_{mode}"], reference["pump"][mode]["A"])
            if mode != "metric":
                verify("mode_coefficients", row[f"w_{mode}"], reference["pump"][mode]["W"])
        if row["clock_running"] == reference["formal_limit_only"]:
            raise AssertionError("Stopped-clock fixture interpreted as a physical curvature mode")

    verify("underconstrained_symbols", production["hessian_rank"]["linear_map_rank"],
           critique["hessian_rank"]["rank"])
    permutation = (0, 2, 1)  # Critic uses X11, X22, X12 coordinates.
    for case, label in (("rank_one_control", "rank_one"), ("rank_two_control", "independent")):
        row = production[case]
        reference = critique[label]
        verify("underconstrained_symbols", row["x_entries"],
               [reference["X_coordinates_11_22_12"][i] for i in permutation])
        verify("underconstrained_symbols", row["determinant"], reference["determinant"])
        hessian = reference["det_squared_hessian"]
        verify("underconstrained_symbols", row["coordinate_hessian"],
               [[hessian[i][j] for j in permutation] for i in permutation])
        verify("underconstrained_symbols", row["projected_null_covector"],
               reference["null_gradient_contractions"])
        verify("underconstrained_symbols", row["principal_delta"], reference["null_principal_correction"])

    reality = production["reality_bounds"]
    reference = critique["band_benchmark"]
    bench = next(row for row in rows if row["x"] == production["bench_x"])
    verify("separation_bounds", reality["cone_gap"], reference["principal_gap"])
    verify("separation_bounds", reality["flat_matter_dispersion_error"],
           reference["flat_completion_error_upper"])
    verify("separation_bounds", reality["flrw_gap_lower_vs_matched"],
           reference["canonical"]["conservative_lower_from_caps"])
    verify("separation_bounds", reality["max_finite_k_flrw_error_vs_matched"],
           reference["canonical"]["curvature_error_upper_from_caps"])
    verify("separation_bounds", reality["max_finite_k_flrw_error_vs_metric"],
           reference["metric"]["curvature_error_upper_from_caps"])
    for control in ("canonical", "metric"):
        verify("separation_bounds", F(bench["pump_active"]) - F(bench[f"pump_{control}"]),
               reference[control]["pump_difference"])
    verify("separation_bounds", production["uniform_flrw_pump_caps"]["metric_cap"],
           critique["pump_bounds"]["metric"])

    c_squared = F(bench["sound_speed_sq"])
    flat_interval = list(map(F, reference["flat_dispersion"]["light_speed_squared_interval"]))
    # The independently bracketed exact completion interval is a subset of the proven cap.
    flat_consistent = c_squared <= flat_interval[0] <= flat_interval[1] <= (
        c_squared + F(reality["flat_matter_dispersion_error"])
    )
    if not (production["passed"] and flat_consistent):
        raise AssertionError("Registered or independent interval test failed")
    return {
        "exact_common_scalars": sum(counts.values()),
        "scalar_checks_by_group": dict(counts),
        "exact_agreement": True,
        "independent_dispersion_interval_inside_cap": flat_consistent,
        "stopped_clock_scope_agrees": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("main", type=Path)
    parser.add_argument("critical", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = compare(json.loads(args.main.read_text()), json.loads(args.critical.read_text()))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{result['exact_common_scalars']} exact cross-checks passed")


if __name__ == "__main__":
    main()
