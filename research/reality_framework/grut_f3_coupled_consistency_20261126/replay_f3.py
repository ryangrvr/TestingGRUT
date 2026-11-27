"""Exact arithmetic certificates for the frozen GRUT-F3 obstruction."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
FREEZE = "0452465ea9a3ba9b29b2543395f66549331d295c"
FIELD_PAIRS = ((0, 0), (0, 1), (1, 1))


def pump(w, w_x, flow, epsilon):
    return w**2 + (1 - epsilon) * w + flow * w_x


def background(x):
    px = 1 + 2 * x
    density_derivative = 1 + 6 * x
    density_factor = 1 + 3 * x
    rho, pressure = x * density_factor, x * (1 + x)
    speed = px / density_derivative
    epsilon = 3 * px / density_factor
    flow = -6 * x * speed

    active_denominator = density_derivative**2 * density_factor
    active_numerator = x * px
    active_w = 1 - 9 * active_numerator / active_denominator
    active_w_x = -9 * (
        1 + 4 * x - active_numerator * (12 / density_derivative + 3 / density_factor)
    ) / active_denominator
    canonical_denominator = density_derivative * density_factor
    canonical_w = 1 + 3 * x / canonical_denominator
    canonical_w_x = 3 * (
        1 - x * (6 / density_derivative + 3 / density_factor)
    ) / canonical_denominator

    kinetic = x * px
    potential = kinetic - pressure
    return {
        "x": x,
        "rho": rho,
        "pressure": pressure,
        "px": px,
        "rho_x": density_derivative,
        "epsilon_h": epsilon,
        "sound_speed_sq": speed,
        "flow": flow,
        "rho_loga": density_derivative * flow,
        "matched_kinetic": kinetic,
        "matched_potential": potential,
        "matched_rho": kinetic + potential,
        "matched_pressure": kinetic - potential,
        "w_active": active_w,
        "w_canonical": canonical_w,
        "pump_active": pump(active_w, active_w_x, flow, epsilon),
        "pump_canonical": pump(canonical_w, canonical_w_x, flow, epsilon),
        "pump_metric": 2 - epsilon,
        "clock_running": x > 0,
    }


def uniform_pump_caps(x_max):
    flow_cap = 6 * x_max
    active_w_lower = 1 - 9 * x_max * (1 + 2 * x_max)
    active_w_x_cap = 9 * (1 + 4 * x_max + 15 * x_max * (1 + 2 * x_max))
    canonical_w_upper = 1 + 3 * x_max
    canonical_w_x_cap = 3 * (1 + 9 * x_max)
    return {
        "active_w_lower": active_w_lower,
        "active_cap": 3 + flow_cap * active_w_x_cap,
        "canonical_cap": canonical_w_upper**2 + 2 * canonical_w_upper
        + flow_cap * canonical_w_x_cap,
        "metric_cap": F(1),
    }


def contracted_hessian(hessian, direction, i, k):
    result = F(0)
    for j, l in product(range(2), repeat=2):
        a = FIELD_PAIRS.index(tuple(sorted((i, j))))
        b = FIELD_PAIRS.index(tuple(sorted((k, l))))
        coordinate_factor = (1 + int(i != j)) * (1 + int(k != l))
        result += hessian[a][b] * direction[j] * direction[l] / coordinate_factor
    return result


def row_rank(matrix):
    echelon = [[F(entry) for entry in row] for row in matrix]
    row = 0
    for col in range(len(echelon[0])):
        pivot = next((i for i in range(row, len(echelon)) if echelon[i][col]), None)
        if pivot is None:
            continue
        echelon[row], echelon[pivot] = echelon[pivot], echelon[row]
        pivot_row = echelon[row]
        for other in echelon[row + 1:]:
            ratio = other[col] / pivot_row[col]
            if ratio:
                other[col:] = [a - ratio * b for a, b in zip(other[col:], pivot_row[col:])]
        row += 1
        if row == len(echelon):
            break
    return row


def null_hessian_rank(vectors):
    basis = []
    for i, j in combinations_with_replacement(range(3), 2):
        unit = [[F(0) for _ in range(3)] for _ in range(3)]
        unit[i][j] = unit[j][i] = F(1)
        basis.append(unit)
    rows = [
        [contracted_hessian(unit, direction, i, k) for unit in basis]
        for direction in vectors for i, k in FIELD_PAIRS
    ]
    return {"linear_map_rank": row_rank(rows), "variables": len(basis), "tested_rows": len(rows)}


def minkowski_inner(x, y):
    return -x[0] * y[0] + sum(a * b for a, b in zip(x[1:], y[1:]))


def determinant_fixture(gradients, null_k, coefficient):
    gradients = [[F(x) for x in covector] for covector in gradients]
    a, b, c = [-minkowski_inner(gradients[i], gradients[j]) / 2 for i, j in FIELD_PAIRS]
    determinant = a * c - b**2
    differential = (c, -2 * b, a)
    differential2 = ((0, 0, 1), (0, -2, 0), (1, 0, 0))
    hessian = [
        [2 * coefficient * (differential[i] * differential[j]
                             + determinant * differential2[i][j])
         for j in range(3)] for i in range(3)
    ]
    direction = [minkowski_inner(gradient, null_k) for gradient in gradients]
    delta = [[-contracted_hessian(hessian, direction, i, j) for j in range(2)] for i in range(2)]
    return {
        "x_entries": [a, b, c],
        "determinant": determinant,
        "coordinate_hessian": hessian,
        "projected_null_covector": direction,
        "principal_delta": delta,
    }


def reality_bounds(benchmark, cfg):
    h_over_k = F(cfg["max_H_over_k_physical"])
    k_over_gap = F(cfg["max_k_over_matter_gap"])
    speed_sq = benchmark["sound_speed_sq"]
    gap = 1 - speed_sq
    active = cfg["expansion_pump_caps_active"]
    metric = cfg["expansion_pump_caps_metric"]
    canonical = cfg["expansion_pump_caps_canonical"]
    flrw_bound = (active + metric) * h_over_k**2
    canonical_bound = (active + canonical) * h_over_k**2
    flat_bound = gap**2 * k_over_gap**2

    # Characteristic polynomial of the two rotating-state matter modes.
    discriminant_increment = 4 * gap * k_over_gap**2
    sqrt_lower = 1 + discriminant_increment / 2 - discriminant_increment**2 / 8
    sqrt_upper = 1 + discriminant_increment / 2
    return {
        "cone_gap": gap,
        "max_finite_k_flrw_error_vs_metric": flrw_bound,
        "max_finite_k_flrw_error_vs_matched": canonical_bound,
        "flat_matter_dispersion_error": flat_bound,
        "flrw_gap_lower_vs_matched": gap - canonical_bound,
        "flat_completion_gap_lower": gap - flat_bound,
        "positive_flrw_frequency_bound": speed_sq - active * h_over_k**2,
        "h_gap_squared_over_v_planck_squared": benchmark["rho"] / (6 * benchmark["rho_x"]),
        "illustrative_observational_error": F(cfg["illustrative_error_allowance"]),
        "sqrt_interval_verified": 0 < sqrt_lower and sqrt_lower**2 <= 1 + discriminant_increment
        <= sqrt_upper**2,
    }


def run(definition):
    model = definition["representative"]
    theorem = definition["rank_check"]
    cfg = definition["reality_test"]
    rows = [background(F(value)) for value in model["x_fixtures"]]
    benchmark = next(row for row in rows if row["x"] == F(model["benchmark_x"]))
    caps = uniform_pump_caps(F(model["upper_x"]))
    rank = null_hessian_rank(theorem["vectors"])
    metric_null = theorem["metric_null_covector"]
    coefficient = theorem["determinant_control_coefficient"]
    parallel = determinant_fixture(theorem["rank_one_gradients"], metric_null, coefficient)
    independent = determinant_fixture(theorem["independent_gradients"], metric_null, coefficient)
    reality = reality_bounds(benchmark, cfg)

    checks = {
        "matching_entire_stress_at_fixtures": all(
            row["rho"] == row["matched_rho"] and row["pressure"] == row["matched_pressure"]
            for row in rows
        ),
        "stress_conservation": all(row["rho_loga"] == -3 * (row["rho"] + row["pressure"])
                                   for row in rows),
        "positive_kinetic_and_nested_cones": all(
            row["px"] > 0 and row["rho_x"] > 0 and 0 < row["sound_speed_sq"] <= 1 for row in rows
        ),
        "whole_band_pump_bounds": caps["active_w_lower"] > 0
        and caps["active_cap"] <= cfg["expansion_pump_caps_active"]
        and caps["canonical_cap"] <= cfg["expansion_pump_caps_canonical"]
        and caps["metric_cap"] <= cfg["expansion_pump_caps_metric"],
        "null_hessian_rank": rank["linear_map_rank"] == theorem["expected_rank"],
        "homogeneous_control_blind": parallel["determinant"] == 0
        and any(value for row in parallel["coordinate_hessian"] for value in row)
        and all(value == 0 for row in parallel["principal_delta"] for value in row),
        "independent_gradient_control_distinguishes": independent["determinant"] != 0
        and any(value for row in independent["principal_delta"] for value in row),
        "metric_null_covector": minkowski_inner(metric_null, metric_null) == 0,
        "effect_floor_in_both_comparisons": min(reality["flrw_gap_lower_vs_matched"],
                                                 reality["flat_completion_gap_lower"])
        > F(cfg["minimum_squared_speed_gap"]),
        "positive_physical_frequencies": reality["positive_flrw_frequency_bound"] > 0,
        "flat_matching_bound": reality["sqrt_interval_verified"],
    }
    return {
        "protocol": definition["protocol"],
        "freeze": FREEZE,
        "background_rows": rows,
        "bench_x": benchmark["x"],
        "uniform_flrw_pump_caps": caps,
        "hessian_rank": rank,
        "rank_one_control": parallel,
        "rank_two_control": independent,
        "reality_bounds": reality,
        "checks": checks,
        "passed": all(checks.values()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    definition = json.loads((HERE / "PARAMETERS.json").read_text())
    result = run(definition)
    args.output.write_text(json.dumps(result, indent=2, default=str, allow_nan=False) + "\n")
    print("All frozen algebra, coupled-background and cone checks:", result["passed"])
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
