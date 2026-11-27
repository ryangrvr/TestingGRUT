#!/usr/bin/env python3
"""Independent exact log-derivative and principal-symbol replay for F3."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import signal


HERE = Path(__file__).resolve().parent


def log_terms(x, terms):
    """Log derivative and its derivative of prod_j (1+c_j x)^power_j."""
    first = sum(power * c / (1 + c * x) for c, power in terms)
    second = -sum(power * c**2 / (1 + c * x)**2 for c, power in terms)
    return first, second


def pump(x, flow, flow_prime, epsilon, log_factors):
    log_factors = [(c, F(sign, 2)) for c, sign in log_factors]
    ell, ell_prime = log_terms(x, log_factors)
    w = 1 + flow * ell
    dw = flow_prime * ell + flow * ell_prime
    return {
        "W": w,
        "W_x": dw,
        "A": w * w + (1 - epsilon) * w + flow * dw,
    }


def physical_row(x):
    rho = x + 3 * x**2
    pressure = x + x**2
    px = 1 + 2 * x
    rho_x = 1 + 6 * x
    c2 = px / rho_x
    epsilon = 3 * (1 + 2 * x) / (1 + 3 * x)
    flow = -6 * x * c2
    flow_prime = -6 * (1 + 4 * x + 12 * x**2) / rho_x**2
    kin = x * px
    pot = kin - pressure
    assert kin + pot == rho and kin - pot == pressure
    assert rho_x >= px > 0
    if x > 0:
        assert epsilon == F(3, 2) * (rho + pressure) / rho
        assert rho_x * flow == -3 * (rho + pressure)
    return {
        "x": x,
        "rho_over_lambda4": rho,
        "p_over_lambda4": pressure,
        "p_x": px,
        "rho_x": rho_x,
        "cs_squared": c2,
        "epsilon_H": epsilon,
        "F": flow,
        "F_prime": flow_prime,
        "canonical_kinetic_over_lambda4": kin,
        "canonical_potential_over_lambda4": pot,
        "formal_limit_only": x == 0,
        "pump": {
            "active": pump(x, flow, flow_prime, epsilon, ((6, 1), (3, -1))),
            "canonical": pump(x, flow, flow_prime, epsilon, ((2, 1), (3, -1))),
            "metric": pump(x, flow, flow_prime, epsilon, ()),
        },
    }


def pump_bounds(upper):
    # 0 <= x <= upper, |F| <= 6x, |F'| <= 6 and |1-epsilon_H| <= 2.
    flow_max = 6 * upper
    active_log_max, active_dlog_max = F(3, 2), F(27, 2)
    canonical_log_max, canonical_dlog_max = F(1, 2), F(5, 2)
    w_active_min = 1 - flow_max * active_log_max
    assert w_active_min > 0
    w_canonical_max = 1 + flow_max * canonical_log_max
    active_pump_bound = 3 + flow_max * (6 * active_log_max + flow_max * active_dlog_max)
    canonical_pump_bound = (w_canonical_max**2 + 2 * w_canonical_max
                            + flow_max * (6 * canonical_log_max + flow_max * canonical_dlog_max))
    return {
        "flow_abs_upper": flow_max,
        "flow_prime_abs_upper": 6,
        "active_W_floor": w_active_min,
        "canonical_W_upper": w_canonical_max,
        "active": active_pump_bound,
        "canonical": canonical_pump_bound,
        "metric": 1,
    }


def band_test(row, limits):
    beta = 1 - row["cs_squared"]
    alpha2 = F(limits["max_H_over_k_physical"])**2
    u = F(limits["max_k_over_matter_gap"])**2
    active = row["pump"]["active"]["A"]
    outputs = {"principal_gap": beta, "flat_completion_error_upper": beta**2 * u}
    for comp in ("canonical", "metric"):
        difference = active - row["pump"][comp]["A"]
        published = limits[f"expansion_pump_caps_{comp}"] + limits["expansion_pump_caps_active"]
        coarse = beta - published * alpha2
        fine = beta + min(F(0), difference) * alpha2
        assert 0 < coarse <= fine and fine > F(limits["minimum_squared_speed_gap"])
        outputs[comp] = {
            "pump_difference": difference,
            "lower_at_any_H_over_k": fine,
            "conservative_lower_from_caps": coarse,
            "curvature_error_upper_from_caps": published * alpha2,
        }
    # Coefficients, in units where m_gap^2=1, of the exact frequency
    # polynomial y^2 - (2u+1)y + u(u+1-beta), y=omega^2/m_gap^2.
    discr = 1 + 4 * beta * u
    root_lo, root_hi = sqrt_interval(discr)
    low_speed_lo = 1 - 2 * beta / (1 + root_lo)
    low_speed_hi = 1 - 2 * beta / (1 + root_hi)
    assert row["cs_squared"] < low_speed_lo <= low_speed_hi
    assert low_speed_hi - row["cs_squared"] <= outputs["flat_completion_error_upper"]
    outputs["flat_dispersion"] = {
        "scaled_wavenumber_squared": u,
        "frequency_polynomial_coefficients_descending": [1, -(2 * u + 1), u * (u + 1 - beta)],
        "light_speed_squared_interval": [low_speed_lo, low_speed_hi],
        "metric_front": "limit of each canonical branch at large k",
    }
    return outputs


def sqrt_interval(value, denominator=10**12):
    numerator2 = value.numerator * denominator**2
    lower = math.isqrt(numerator2 // value.denominator)
    return F(lower, denominator), F(lower + 1, denominator)


H_POSITIONS = list(itertools.combinations_with_replacement(range(3), 2))
POLARIZATIONS = ((0, 0), (0, 1), (1, 1))


def lifted(u):
    a, b = map(F, u)
    return [[a, 0], [0, b], [b / 2, a / 2]]


def rank_rows(directions):
    rows = []
    for u in directions:
        s = lifted(u)
        for a, b in POLARIZATIONS:
            row = []
            for i, j in H_POSITIONS:
                contribution = s[i][a] * s[j][b]
                if i != j:
                    contribution += s[j][a] * s[i][b]
                row.append(F(contribution))
            rows.append(row)
    return rows


def eliminate(rows):
    a = [row.copy() for row in rows]
    pivot_columns = []
    pivot = 0
    for j in range(len(a[0])):
        row = next((i for i in range(pivot, len(a)) if a[i][j]), None)
        if row is None:
            continue
        a[pivot], a[row] = a[row], a[pivot]
        scale = a[pivot][j]
        a[pivot] = [v / scale for v in a[pivot]]
        for k in range(len(a)):
            if k == pivot:
                continue
            g = a[k][j]
            a[k] = [left - g * right for left, right in zip(a[k], a[pivot], strict=True)]
        pivot_columns.append(j)
        pivot += 1
    return {"rank": pivot, "pivot_columns": pivot_columns, "rref": a}


def lorentz_pair(v, w):
    return -v[0] * w[0] + sum(x * y for x, y in zip(v[1:], w[1:], strict=True))


def fixture(gradients, k, coefficient):
    ga, gb = (tuple(map(F, vector)) for vector in gradients)
    a, b, c = (-lorentz_pair(i, j) / 2 for i, j in ((ga, ga), (gb, gb), (ga, gb)))
    det = a * b - c * c
    d = (b, a, -2 * c)
    d2 = [[0, 1, 0], [1, 0, 0], [0, 0, -2]]
    h = [[2 * coefficient * (d[i] * d[j] + det * d2[i][j])
          for j in range(3)] for i in range(3)]
    u = [lorentz_pair(grad, k) for grad in (ga, gb)]
    s = lifted(u)
    correction = [[-sum(s[i][a] * h[i][j] * s[j][b] for i in range(3) for j in range(3))
                   for b in range(2)] for a in range(2)]
    return {"X_coordinates_11_22_12": [a, b, c], "determinant": det,
            "det_squared_hessian": h, "null_gradient_contractions": u,
            "null_principal_correction": correction}


def run(spec):
    rep = spec["representative"]
    r = spec["reality_test"]
    rows = [physical_row(F(x)) for x in rep["x_fixtures"]]
    upper = F(rep["upper_x"])
    assert rep["kinetic_correction_coefficient"] == 1
    bounds = pump_bounds(upper)
    for mode in ("active", "canonical", "metric"):
        assert bounds[mode] <= r[f"expansion_pump_caps_{mode}"]
        assert all(abs(row["pump"][mode]["A"]) <= bounds[mode] for row in rows)
    benchmark = next(row for row in rows if row["x"] == F(rep["benchmark_x"]))
    band = band_test(benchmark, r)
    rankcheck = spec["rank_check"]
    rankmatrix = rank_rows(rankcheck["vectors"])
    rank = eliminate(rankmatrix)
    assert rank["rank"] == rankcheck["expected_rank"]
    k = rankcheck["metric_null_covector"]
    assert lorentz_pair(k, k) == 0
    multiplier = rankcheck["determinant_control_coefficient"]
    parallel = fixture(rankcheck["rank_one_gradients"], k, multiplier)
    independent = fixture(rankcheck["independent_gradients"], k, multiplier)
    assert parallel["null_principal_correction"] == [[0, 0], [0, 0]]
    assert any(map(any, parallel["det_squared_hessian"]))
    assert independent["null_principal_correction"] == [[F(-1, 8192), 0], [0, F(-1, 16384)]]
    return {
        "method": "independent Fraction log derivatives, null-symbol coefficient elimination, and interval check",
        "fixtures": rows,
        "pump_bounds": bounds,
        "band_benchmark": band,
        "hessian_rank": rank | {"matrix": rankmatrix},
        "rank_one": parallel,
        "independent": independent,
    }


def timeout(*_):
    raise TimeoutError("Frozen replay wall limit exceeded")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads((HERE / "PARAMETERS.json").read_text())
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(cfg["compute"]["wall_seconds_per_implementation"])
    result = run(cfg)
    args.output.write_text(json.dumps(result, indent=2, default=str, allow_nan=False) + "\n")
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ("method", "pump_bounds", "band_benchmark")},
                     indent=2, default=str))


if __name__ == "__main__":
    main()
