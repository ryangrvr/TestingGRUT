#!/usr/bin/env python3
"""Independent rational LDL replay of the frozen F2 intervention data."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import signal


HERE = Path(__file__).resolve().parent


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def spring_stiffness(n, edges, pin):
    result = [[pin * x for x in row] for row in eye(n)]
    for i, j, weight in edges:
        assert i != j and weight > 0
        result[i][i] += weight
        result[j][j] += weight
        result[i][j] -= weight
        result[j][i] -= weight
    return result


def ldl(matrix):
    """Factor an exact positive definite matrix as L D L transpose."""
    n = len(matrix)
    lower = [[F(x) for x in row] for row in eye(n)]
    diagonal = []
    for j in range(n):
        pivot = F(matrix[j][j]) - sum(lower[j][k] ** 2 * diagonal[k] for k in range(j))
        if pivot <= 0:
            raise AssertionError("Matrix lacks positive pivot")
        diagonal.append(F(pivot))
        for i in range(j + 1, n):
            residual = F(matrix[i][j]) - sum(lower[i][k] * lower[j][k] * diagonal[k]
                                             for k in range(j))
            lower[i][j] = residual / pivot
    return lower, diagonal


def solve_ldl(factor, unit_index):
    lower, diagonal = factor
    n = len(diagonal)
    forward = []
    for i in range(n):
        forward.append(int(i == unit_index) - sum(lower[i][j] * forward[j] for j in range(i)))
    scaled = [x / d for x, d in zip(forward, diagonal, strict=True)]
    backward = [F(0)] * n
    for i in range(n - 1, -1, -1):
        backward[i] = scaled[i] - sum(lower[j][i] * backward[j] for j in range(i + 1, n))
    return backward


def inverse(matrix):
    factor = ldl(matrix)
    columns = [solve_ldl(factor, i) for i in range(len(matrix))]
    return list(map(list, zip(*columns, strict=True)))


def response(v, z, ports=(0,), loads=(), mass=1):
    matrix = [row[:] for row in v]
    for i in range(len(v)):
        matrix[i][i] += mass * z
    for i, change in loads:
        matrix[i][i] += change
    factor = ldl(matrix)
    columns = [solve_ldl(factor, j) for j in ports]
    return [[column[i] for column in columns] for i in ports]


def root(v, z, loads=(), mass=1):
    return response(v, z, loads=loads, mass=mass)[0][0]


def sqrt_fraction(x):
    assert x >= 0
    a, b = math.isqrt(x.numerator), math.isqrt(x.denominator)
    assert a * a == x.numerator and b * b == x.denominator
    return F(a, b)


def reconstruct(v, z, single_levels, pair_level):
    assert single_levels == [1, 2]
    tags = tuple(range(1, len(v)))
    r = root(v, z)
    singleton = {str(i): [root(v, z, ((i, load),)) for load in single_levels] for i in tags}
    d, u = {}, {}
    for i in tags:
        x, y = [r - reading for reading in singleton[str(i)]]
        assert 0 < x < y < 2 * x
        d[i] = (2 * x - y) / (2 * (y - x))
        u[i] = sqrt_fraction(x * (1 + d[i]))
    pair_values = {}
    recovered_r = [[F(0) for _ in v] for _ in v]
    recovered_r[0][0] = r
    for i in tags:
        recovered_r[i][0] = recovered_r[0][i] = u[i]
        recovered_r[i][i] = d[i]
    for i, j in itertools.combinations(tags, 2):
        reading = root(v, z, ((i, pair_level), (j, pair_level)))
        pair_values[f"{i},{j}"] = reading
        drop = r - reading
        ai, aj = d[i] + F(1, pair_level), d[j] + F(1, pair_level)
        discriminant = (u[i] * u[j]) ** 2 - drop * (
            aj * u[i] ** 2 + ai * u[j] ** 2 - drop * ai * aj)
        numerator = u[i] * u[j]
        square_root = sqrt_fraction(discriminant)
        choices = [(numerator - square_root) / drop, (numerator + square_root) / drop]
        upper = min(d[j] * u[i] / u[j], d[i] * u[j] / u[i])
        physical = [c for c in choices if 0 < c <= upper]
        assert len(physical) == 1
        recovered_r[i][j] = recovered_r[j][i] = physical[0]
    recovered_v = inverse(recovered_r)
    for i in range(len(v)):
        recovered_v[i][i] -= z
    assert recovered_v == v
    return {
        "baseline": r, "singletons": singleton, "pairs": pair_values,
        "resolvent": recovered_r, "stiffness": recovered_v,
    }


def auxiliary(v, pin, spring, channel):
    n = len(v)
    vector = channel + [-1]
    expanded = [row + [0] for row in v] + [[0] * n + [pin]]
    return [[x + spring * vector[i] * vector[j] for j, x in enumerate(row)]
            for i, row in enumerate(expanded)]


def nonlinear(v, edge_set, target):
    seen = {0}
    frontier = {0}
    distance = 0
    while target not in seen:
        frontier = {b if a in frontier else a for a, b, _ in edge_set
                    if a in frontier or b in frontier} - seen
        seen.update(frontier)
        distance += 1
    column = eye(len(v))[0]
    for _ in range(distance):
        column = [sum(a * b for a, b in zip(row, column, strict=True)) for row in v]
    power = 2 * distance + 1
    coefficient = F((-1)**distance * column[target], math.factorial(power))
    exponent = 4 * power + 1
    return {
        "distance": distance,
        "root_to_tag_coefficient": coefficient,
        "root_to_tag_order": power,
        "quartic_time_order": exponent,
        "quartic_coefficient": -coefficient**4 * F(math.factorial(power) *
                                                    math.factorial(3 * power),
                                                    math.factorial(exponent)),
    }


def paired_stress(spring):
    s = spring["stress"]
    vertices = list(itertools.product(range(4), repeat=2))
    assert len(vertices) == s["hidden_sites"]
    steps = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    by_name = {}
    for name in ("rook", "shrikhande"):
        adjacency = [[] for _ in vertices]
        edges = []
        for i, j in itertools.combinations(range(len(vertices)), 2):
            a, b = vertices[i], vertices[j]
            joined = (sum(x == y for x, y in zip(a, b, strict=True)) == 1
                      if name == "rook" else
                      tuple((x - y) % 4 for x, y in zip(a, b, strict=True)) in steps)
            if joined:
                adjacency[i].append(j)
                adjacency[j].append(i)
                edges.append((i + 1, j + 1, s["hidden_edge_weight"]))
        assert len(edges) == 48 and all(len(v) == 6 for v in adjacency)
        assert all(len(set(adjacency[a]) & set(adjacency[b])) == 2
                   for a, b in itertools.combinations(range(len(vertices)), 2))
        cliques = sum(all(b in adjacency[a] for a, b in itertools.combinations(group, 2))
                      for group in itertools.combinations(range(len(vertices)), 4))
        rooted_edges = [(0, i + 1, s["root_edge_weight"]) for i in range(len(vertices))] + edges
        matrix = spring_stiffness(len(vertices) + 1, rooted_edges, s["pin"])
        targets = [vertices.index(tuple(coord)) + 1 for coord in s["pair"]]
        by_shift = {}
        for z in s["replica_shifts"]:
            base = root(matrix, z)
            common = (F(17 + z) - F(16, 2 + z)) ** -1
            assert base == common
            singleton = {
                str(i): [root(matrix, z, ((i, load),))
                         for load in spring["pin_shift_levels"]]
                for i in range(1, len(matrix))
            }
            hidden_h = F(8 + z)
            hh_inv_diag = hidden_h / (hidden_h ** 2 - 4)
            hh_inv_diag += 2 / ((hidden_h ** 2 - 4) * (hidden_h - 6))
            di = hh_inv_diag + base / (hidden_h - 6) ** 2
            ui = base / (hidden_h - 6)
            expected = [base - load * ui**2 / (1 + load * di)
                        for load in spring["pin_shift_levels"]]
            assert all(values == expected for values in singleton.values())
            by_shift[str(z)] = {
                "baseline": base,
                "singletons": singleton,
                "pair": root(matrix, z, tuple((i, spring["simultaneous_pin_level"])
                                               for i in targets)),
            }
        by_name[name] = {
            "four_cliques": cliques,
            "rooted_edge_count": len(rooted_edges),
            "by_shift": by_shift,
        }
    assert by_name["rook"]["four_cliques"] == 8
    assert by_name["shrikhande"]["four_cliques"] == 0
    for z in s["replica_shifts"]:
        l, r = (g["by_shift"][str(z)] for g in by_name.values())
        assert l["singletons"] == r["singletons"] and l["pair"] > r["pair"]
    return {
        "srg": [16, 6, 2, 2],
        "hidden_adjacency_polynomial": "(x-6)(x-2)^6(x+2)^9",
        "pair": s["pair"],
        "graphs": by_name,
    }


def spacetime(c):
    # Radial support has a constant plateau through r=1/16. Bound all
    # terms separately, including the 2/r coefficient of the Laplacian.
    shape = c["shape_control"]
    b_cap = 2 * 8**2 * shape["second_derivative_cap"]
    b_cap += 16 * shape["derivative_cap"] / F(c["plateau_radius"])
    assert b_cap == shape["box_cap"]
    width = 2 * F(c["bump_temporal_half_width"])
    volume_pi = width * F(4, 3) * F(c["bump_radius"])**3
    plateau_pi = 2 * F(c["plateau_temporal_half_width"])
    plateau_pi *= F(4, 3) * F(c["plateau_radius"])**3
    assert plateau_pi == F(1, 24576) and volume_pi / plateau_pi == 16
    inner = F(c["center_x"][0]) - F(c["bump_radius"])
    outer = F(c["center_x"][0]) + F(c["bump_radius"])
    distance_min = inner - F(c["source_spatial_radius"])
    assert c["center_t"] - width / 2 + distance_min > (
        c["probe_times"]["early"] + F(c["source_time_half_width"]))
    assert c["deformations"]["unobserved_center_t"] - width / 2 > c["tube_time"][1]
    # |J_1(x)| <= x/2 gives a tail m^2/(8 pi); use the frozen
    # conservative 1/6144 coefficient in the kernel bound.
    m2 = F(c["channels"]["massive_conformal"]["mass_squared"])
    free_lo, free_hi = F(4, 5) - m2 / 2, F(4, 3) + m2 / 2
    kernel_cap = width**2 / 2 + m2 / 6144
    threshold = c["thresholds"]
    positive = F(threshold["positive_scaled_wave_on_bump"])
    output = {}
    for name in ("weak", "strong"):
        eps = F(c["deformations"][name])
        potential = m2 * (2 * eps + eps**2)
        contraction = kernel_cap * potential
        wave_floor = free_lo - free_hi * contraction / (1 - contraction)
        mass_gap = positive**2 * potential * plateau_pi
        assert wave_floor > positive
        assert mass_gap > F(threshold[f"min_{name}_mass_signal"])
        output[name] = {
            "mass_kernel_contraction": contraction,
            "mass_wave_floor": wave_floor,
            "mass_gap_lower": mass_gap,
        }
    eps = F(c["deformations"]["weak"])
    lo = 2 * plateau_pi / outer**4
    hi = 2 * volume_pi / inner**4
    scale = eps * b_cap / 32
    remainder = volume_pi * F(4, 3)**2 * eps**2 * (b_cap + b_cap**2 / (32 * (1 - scale)))
    leading = eps * lo
    gap = leading - remainder
    assert remainder / leading < F(threshold["max_weak_curvature_remainder_to_linear_floor"])
    assert gap > F(threshold["min_weak_curvature_signal"])
    output["weak"].update(
        minimal_slope_lower=lo,
        minimal_slope_upper=hi,
        minimal_contraction=scale,
        minimal_remainder_upper=remainder,
        minimal_remainder_to_linear_floor=remainder / leading,
        minimal_gap_lower=gap,
    )
    output["strong"].update(minimal_contraction=F(c["deformations"]["strong"]) * b_cap / 32,
                            minimal_status="inconclusive under registered bound")
    return {
        "box_cap": b_cap,
        "volume_over_pi": volume_pi,
        "plateau_volume_over_pi": plateau_pi,
        "free_wave_bounds": [free_lo, free_hi],
        "mass_kernel_cap": kernel_cap,
        "minimal_to_massive_linearized_ratio_bounds": [
            1 / (outer**4 * m2 * free_hi**2),
            1 / (inner**4 * m2 * free_lo**2),
        ],
        "early_causal_blind": True,
        "future_causal_blind": True,
        "channels": output,
    }


def run(spec):
    f1 = (HERE / spec["f1_inputs_path"]).read_bytes()
    assert hashlib.sha256(f1).hexdigest() == spec["f1_inputs_sha256"]
    base = json.loads(f1)
    spring = spec["spring"]
    assert spring["reference_mass"] == spring["aux_mass"] == 1 and spring["observe"] == 0
    assert 1 + spring["stress"]["hidden_sites"] <= spec["compute"]["max_sites"]
    primary_z = spring["primary_z"]
    models = {}
    for name, edges in base["networks"].items():
        v = spring_stiffness(4, edges, base["units"]["common_pinning_frequency_squared"])
        auxiliary_v = [auxiliary(v, spring["aux_pin"], spring["aux_spring"], spring[channel])
                       for channel in ("collective_channel", "local_channel")]
        shifts = {}
        for z in spring["z"]:
            shifts[str(z)] = {
                "baseline": root(v, z),
                "root_pin": root(v, z, ((0, spring["root_pin_control"]),)),
                "uniform_mass": root(v, z, mass=spring["uniform_mass_multiplier"]),
                "pins": {str(i): root(v, z, ((i, 1),)) for i in range(1, len(v))},
                "masses": {str(i): root(v, z, ((i, z * spring["mass_increment"]),))
                           for i in range(1, len(v))},
                "aux_collective": response(auxiliary_v[0], z, ports=(0, 4)),
                "aux_local": response(auxiliary_v[1], z, ports=(0, 4)),
            }
        quartic = nonlinear(v, edges, spring["primary_tag"])
        assert quartic["quartic_time_order"] == spring["expected_quartic_time_orders"][name]
        models[name] = {
            "by_shift": shifts,
            "reconstruction": reconstruct(v, primary_z, spring["pin_shift_levels"],
                                          spring["simultaneous_pin_level"]),
            "quartic": quartic,
        }
    collisions = {}
    for z in spring["z"]:
        left, right = (m["by_shift"][str(z)] for m in models.values())
        assert all(left[k] == right[k] for k in ("baseline", "root_pin", "uniform_mass", "aux_collective"))
        assert left["aux_local"] != right["aux_local"]
        collisions[str(z)] = {kind: [[i, j] for i, j in itertools.product(left[kind], right[kind])
                                      if left[kind][i] == right[kind][j]]
                              for kind in ("pins", "masses")}
    tag = str(spring["primary_tag"])
    left, right = (m["by_shift"][str(primary_z)] for m in models.values())
    assert primary_z == 1 and left["pins"] == left["masses"] and right["pins"] == right["masses"]
    gap = abs(left["pins"][tag] - right["pins"][tag])
    time = spring["finite_time_cutoff"]
    # Each transfer has amplitude bounded by sqrt(2)<3/2. One positive
    # Taylor term of exp(time) suffices for an exact upper bound.
    tail_bound = F(3 * math.factorial(time), time**time)
    assert gap > tail_bound
    return {
        "method": "independent exact LDL solves and rational interval certificates",
        "spring": models,
        "collisions": collisions,
        "finite_time": {"untruncated_gap": gap, "sum_tail_upper": tail_bound,
                        "gap_lower": gap - tail_bound},
        "stress": paired_stress(spring),
        "continuum": spacetime(spec["continuum"]),
    }


def timeout(*_):
    raise TimeoutError("Frozen 30-second replay limit exceeded")


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
    print(json.dumps({key: result[key] for key in ("method", "collisions", "finite_time", "continuum")},
                     indent=2, default=str))


if __name__ == "__main__":
    main()
