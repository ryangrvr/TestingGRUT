"""Replay the frozen F2 spring experiments and analytic interval certificates."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import factorial, isqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent
FREEZE = "2db0de4c0c50aed3a7297964fc8ca9fedaa96991"


def inverse(matrix):
    size = len(matrix)
    augmented = [
        [F(a) for a in row] + [F(i == j) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for i in range(size):
        pivot = next(j for j in range(i, size) if augmented[j][i])
        augmented[i], augmented[pivot] = augmented[pivot], augmented[i]
        row = augmented[i]
        scale = row[i]
        row[i:] = [a / scale for a in row[i:]]
        for j, current in enumerate(augmented):
            if i != j and current[i]:
                factor = current[i]
                current[i:] = [a - factor * b for a, b in zip(current[i:], row[i:])]
    return [row[size:] for row in augmented]


def stiffness(size, edges, pin):
    matrix = [[pin if i == j else 0 for j in range(size)] for i in range(size)]
    for i, j, weight in edges:
        matrix[i][i] += weight
        matrix[j][j] += weight
        matrix[i][j] -= weight
        matrix[j][i] -= weight
    return matrix


def resolvent(stiffness_matrix, z, pins=None, masses=None):
    matrix = [row.copy() for row in stiffness_matrix]
    pins = {} if pins is None else pins
    masses = {} if masses is None else masses
    for i in range(len(matrix)):
        matrix[i][i] += z * masses.get(i, 1) + pins.get(i, 0)
    return inverse(matrix)


def rational_sqrt(number):
    root = F(isqrt(number.numerator), isqrt(number.denominator))
    if root * root != number:
        raise ValueError("Frozen rational fixture lacks its exact positive square root")
    return root


def recover(base, singleton, pairs, shifts, pair_shift, z):
    """Recover the stiffness using only the measured root scalars."""
    size = len(singleton) + 1
    green = [[F(0) for _ in range(size)] for _ in range(size)]
    green[0][0] = base
    d1, d2 = shifts
    for tag, (response1, response2) in singleton.items():
        a = d1 / (base - response1)
        b = d2 / (base - response2)
        slope = (b - a) / (d2 - d1)
        root_cross_squared = 1 / (a - d1 * slope)
        green[0][tag] = green[tag][0] = rational_sqrt(root_cross_squared)
        green[tag][tag] = slope * root_cross_squared

    for (i, j), response in pairs.items():
        loss = base - response
        u_i, u_j = green[0][i], green[0][j]
        v_i, v_j = green[i][i], green[j][j]
        alpha, beta = v_i + F(1, pair_shift), v_j + F(1, pair_shift)
        constant = beta * u_i**2 + alpha * u_j**2 - loss * alpha * beta
        radical = rational_sqrt((u_i * u_j)**2 - loss * constant)
        candidates = {(u_i * u_j + radical) / loss, (u_i * u_j - radical) / loss}
        ceiling = min(v_i * u_j / u_i, v_j * u_i / u_j)
        admissible = [x for x in candidates if 0 <= x <= ceiling]
        if len(admissible) != 1:
            raise ValueError("The positivity interval did not select a unique cross entry")
        green[i][j] = green[j][i] = admissible[0]
    matrix = inverse(green)
    for i in range(size):
        matrix[i][i] -= z
    return matrix


def auxiliary_green(matrix, z, vector, cfg):
    size = len(matrix)
    extended = [row.copy() + [0] for row in matrix] + [[0] * size + [cfg["aux_pin"]]]
    coupling = [*vector, -1]
    for i, j in product(range(size + 1), repeat=2):
        extended[i][j] += cfg["aux_spring"] * coupling[i] * coupling[j]
    green = resolvent(extended, z, masses={size: cfg["aux_mass"]})
    return [[green[i][j] for j in (0, size)] for i in (0, size)]


def quartic_response(matrix, tag):
    column = [int(i == tag) for i in range(len(matrix))]
    for k in range(len(matrix)):
        if column[0]:
            travel_power = 2 * k + 1
            travel_coefficient = (-1)**k * column[0]
            output_power = 4 * travel_power + 1
            coefficient = F(
                -travel_coefficient**4 * factorial(3 * travel_power),
                factorial(travel_power)**3 * factorial(output_power)
            )
            return {"travel_power": travel_power, "travel_derivative": travel_coefficient,
                    "cubic_time_power": output_power, "cubic_coefficient": coefficient}
        column = [sum(x * y for x, y in zip(row, column)) for row in matrix]
    raise ValueError("A disconnected tag fell outside the model")


def check_f1(definition, frozen):
    cfg = definition["spring"]
    names = ("cycle", "pendant_triangle")
    tag = cfg["primary_tag"]
    shifts = cfg["pin_shift_levels"]
    mass_shift = cfg["mass_increment"]
    matrices = {name: stiffness(4, frozen["networks"][name], 1) for name in names}
    rows = []
    for z in cfg["z"]:
        networks = {}
        for name, matrix in matrices.items():
            base = resolvent(matrix, z)[0][0]
            singleton = {i: [resolvent(matrix, z, pins={i: d})[0][0] for d in shifts]
                         for i in range(1, 4)}
            masses = {i: resolvent(matrix, z, masses={i: 1 + mass_shift})[0][0]
                      for i in range(1, 4)}
            row = {
                "base": base,
                "uniform_mass": resolvent(matrix, z, masses=dict.fromkeys(
                    range(4), cfg["uniform_mass_multiplier"]))[0][0],
                "root_pin": resolvent(matrix, z, pins={0: cfg["root_pin_control"]})[0][0],
                "single_pins": singleton,
                "single_masses": masses,
                "collective_channel": auxiliary_green(matrix, z, cfg["collective_channel"], cfg),
                "local_channel": auxiliary_green(matrix, z, cfg["local_channel"], cfg),
            }
            if z == cfg["primary_z"]:
                d = cfg["simultaneous_pin_level"]
                pairs = {(i, j): resolvent(matrix, z, pins={i: d, j: d})[0][0]
                         for i, j in combinations(range(1, 4), 2)}
                row["pair_pins"] = {f"{i},{j}": y for (i, j), y in pairs.items()}
                row["reconstructed_stiffness"] = recover(base, singleton, pairs, shifts, d, z)
            networks[name] = row

        cycle, triangle = (networks[name] for name in names)
        controls = {key: cycle[key] == triangle[key]
                    for key in ("base", "uniform_mass", "root_pin", "collective_channel")}
        contrasts = {
            "tagged_pin_gap": abs(cycle["single_pins"][tag][0] - triangle["single_pins"][tag][0]),
            "tagged_mass_gap": abs(cycle["single_masses"][tag] - triangle["single_masses"][tag]),
            "local_channel_distinguishes": cycle["local_channel"] != triangle["local_channel"],
            "single_pin_cross_hypothesis_collisions": [
                [i, j] for i, j in product(range(1, 4), repeat=2)
                if cycle["single_pins"][i][0] == triangle["single_pins"][j][0]
            ],
            "single_mass_cross_hypothesis_collisions": [
                [i, j] for i, j in product(range(1, 4), repeat=2)
                if cycle["single_masses"][i] == triangle["single_masses"][j]
            ],
        }
        rows.append({"z": z, "networks": networks, "equivalence_controls": controls,
                     "contrasts": contrasts})

    primary = next(r for r in rows if r["z"] == cfg["primary_z"])
    cutoff = cfg["finite_time_cutoff"]
    finite_tail_cap = 4 * F(2, 5)**cutoff
    quartic = {name: quartic_response(matrix, tag) for name, matrix in matrices.items()}
    checks = {
        "invariant_controls": all(all(r["equivalence_controls"].values()) for r in rows),
        "exact_primary_reconstruction": all(
            row["reconstructed_stiffness"] == matrices[name]
            for name, row in primary["networks"].items()
        ),
        "single_tag_and_independent_channel": all(
            r["contrasts"]["tagged_pin_gap"] > 0
            and r["contrasts"]["tagged_mass_gap"] > 0
            and r["contrasts"]["local_channel_distinguishes"] for r in rows
        ),
        "finite_time_gap": primary["contrasts"]["tagged_mass_gap"] > finite_tail_cap,
        "pin_mass_same_at_z_one": all(
            data["single_pins"][tag][0] == data["single_masses"][tag]
            for data in primary["networks"].values()
        ),
        "fixed_quartic_orders": all(
            r["cubic_time_power"] == cfg["expected_quartic_time_orders"][name]
            for name, r in quartic.items()
        ),
    }
    return {"rows": rows, "quartic": quartic, "finite_tail_cap": finite_tail_cap, "checks": checks}


def hidden_adjacency(name, labels):
    shrikhande_steps = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    matrix = []
    for a, b in labels:
        row = []
        for c, d in labels:
            if name == "rook":
                near = (a == c) != (b == d)
            else:
                near = ((c - a) % 4, (d - b) % 4) in shrikhande_steps
            row.append(int(near))
        matrix.append(row)
    return matrix


def has_four_clique(adjacency):
    return any(all(adjacency[i][j] for i, j in combinations(nodes, 2))
               for nodes in combinations(range(len(adjacency)), 4))


def check_stress(cfg, max_sites):
    stress = cfg["stress"]
    labels = list(product(range(4), repeat=2))
    names = ("rook", "shrikhande")
    n = len(labels)
    if n != stress["hidden_sites"] or n + 1 > max_sites:
        raise ValueError("Stress fixture changed the frozen size")
    adjacency = {name: hidden_adjacency(name, labels) for name in names}
    matrices = {}
    stats = {}
    for name, adj in adjacency.items():
        edges = [(0, i + 1, stress["root_edge_weight"]) for i in range(n)]
        edges.extend((i + 1, j + 1, stress["hidden_edge_weight"])
                     for i, j in combinations(range(n), 2) if adj[i][j])
        matrices[name] = stiffness(n + 1, edges, stress["pin"])
        stats[name] = {
            "n": n + 1, "springs": len(edges), "total_spring_strength": sum(w for _, _, w in edges),
            "has_hidden_four_clique": has_four_clique(adj),
            "degree_6": all(sum(row) == 6 for row in adj),
            "adjacency_identity": all(
                sum(adj[i][k] * adj[k][j] for k in range(n)) == 4 * (i == j) + 2
                for i, j in product(range(n), repeat=2)
            ),
        }

    pair = [labels.index(tuple(tag)) + 1 for tag in stress["pair"]]
    rows = []
    for z in stress["replica_shifts"]:
        data = {}
        for name, matrix in matrices.items():
            green = resolvent(matrix, z)
            singletons = [
                [green[0][0] - d * green[0][i]**2 / (1 + d * green[i][i])
                 for d in cfg["pin_shift_levels"]]
                for i in range(1, n + 1)
            ]
            joint = resolvent(matrix, z, pins=dict.fromkeys(pair, cfg["simultaneous_pin_level"]))
            data[name] = {"base": green[0][0], "single_pins": singletons, "pair_pin": joint[0][0],
                          "cross_entry": green[pair[0]][pair[1]]}
        left, right = (data[name] for name in names)
        rows.append({
            "z": z, "networks": data,
            "singleton_data_equal": left["base"] == right["base"]
            and left["single_pins"] == right["single_pins"],
            "pair_gap": abs(left["pair_pin"] - right["pair_pin"]),
        })
    return {"graphs": stats, "rows": rows, "checks": {
        "strongly_regular_identities": all(
            row["degree_6"] and row["adjacency_identity"] for row in stats.values()),
        "nonisomorphic_equal_resources": stats["rook"]["has_hidden_four_clique"]
        and not stats["shrikhande"]["has_hidden_four_clique"]
        and stats["rook"]["springs"] == stats["shrikhande"]["springs"]
        and stats["rook"]["total_spring_strength"] == stats["shrikhande"]["total_spring_strength"],
        "all_singletons_invariant": all(row["singleton_data_equal"] for row in rows),
        "frozen_double_load_distinguishes": all(row["pair_gap"] > 0 for row in rows),
    }}


def check_continuum(cfg):
    time_radius = F(cfg["bump_temporal_half_width"])
    r_bump = F(cfg["bump_radius"])
    r_device = F(cfg["source_spatial_radius"])
    t_half = F(cfg["source_time_half_width"])
    r_min = 1 - r_bump - r_device
    r_max = 1 + r_bump + r_device
    width = 2 * time_radius
    cap_shape = cfg["shape_control"]
    # Derivative bounds of the smoothed profile, including 2/r in the radial Laplacian.
    box_derivative_bound = (
        2 * cap_shape["second_derivative_cap"] / r_bump**2
        + 2 * cap_shape["derivative_cap"] / (r_bump * F(cfg["plateau_radius"]))
    )
    box_bound = cap_shape["box_cap"]
    volume = width * F(4, 3) * r_bump**3  # divided by pi
    volume_core = 2 * F(cfg["plateau_temporal_half_width"]) * F(4, 3) * F(cfg["plateau_radius"])**3
    mass_squared = F(cfg["channels"]["massive_conformal"]["mass_squared"])
    strong = F(cfg["deformations"]["strong"])
    weak = F(cfg["deformations"]["weak"])
    cutoffs = cfg["thresholds"]

    massive_kernel_bound = width**2 / 2 + mass_squared * width**4 / 24
    massless_kernel_bound = width**2 / 2
    massive_tail_bound = mass_squared * (2 * t_half) / 2  # wave times 4*pi
    wave_floor, wave_cap = 1 / r_max - massive_tail_bound, 1 / r_min + massive_tail_bound
    strong_potential = mass_squared * (2 * strong + strong**2)
    alpha = strong_potential * massive_kernel_bound
    distorted_floor = wave_floor - alpha / (1 - alpha) * wave_cap
    chosen_floor = F(cutoffs["positive_scaled_wave_on_bump"])
    mass_lower = {
        name: chosen_floor**2 * mass_squared * (2 * e + e**2) * volume_core
        for name, e in (("strong", strong), ("weak", weak))
    }
    slope_floor = 2 * volume_core / (1 + r_bump)**4
    alpha_weak = weak * box_bound * massless_kernel_bound
    slope_error = volume * (1 / r_min)**2 * weak**2 * (
        box_bound + box_bound**2 * massless_kernel_bound / (1 - alpha_weak)
    )
    error_ratio = slope_error / (weak * slope_floor)
    curvature_lower = weak * slope_floor - slope_error
    correction = (1 + r_bump) * massive_tail_bound
    slope_ratio = [
        mass_squared * (1 - r_bump)**2 * (1 - correction)**2,
        mass_squared * (1 + r_bump)**2 * (1 + correction)**2,
    ]
    arrival_floor = cfg["center_t"] - time_radius + r_min
    early_close = cfg["probe_times"]["early"] + t_half
    return {
        "box_derivative_bound": box_derivative_bound,
        "volume_over_pi": volume, "core_volume_over_pi": volume_core,
        "early_close": early_close, "earliest_bump_arrival": arrival_floor,
        "free_massive_wave_bound": [wave_floor, wave_cap],
        "massive_neumann_norm": alpha, "massive_distorted_wave_floor": distorted_floor,
        "massive_signal_lower": mass_lower,
        "massive_monotone_slope_magnitude_floor": chosen_floor**2 * 2 * mass_squared * volume_core,
        "weak_minimal_neumann_norm": alpha_weak,
        "weak_minimal_linear_floor": weak * slope_floor,
        "weak_minimal_remainder_bound": slope_error,
        "weak_minimal_remainder_ratio": error_ratio,
        "weak_minimal_signal_lower": curvature_lower,
        "strong_minimal_norm_certificate": strong * box_bound * massless_kernel_bound,
        "strong_minimal_status": "inconclusive by this bound",
        "linearized_signal_ratio_interval": slope_ratio,
        "checks": {
            "exact_aperture_and_timing": (time_radius + r_bump + r_device
                <= F(cfg["source_plateau_half_width"]) and early_close < arrival_floor),
            "profile_cap_verified": box_derivative_bound <= box_bound,
            "massive_positivity": 0 <= alpha < 1 and distorted_floor > chosen_floor,
            "strong_massive_lower_gate": mass_lower["strong"] > F(cutoffs["min_strong_mass_signal"]),
            "weak_massive_lower_gate": mass_lower["weak"] > F(cutoffs["min_weak_mass_signal"]),
            "weak_curvature_remainder_gate": error_ratio < F(
                cutoffs["max_weak_curvature_remainder_to_linear_floor"]),
            "weak_curvature_signal_gate": curvature_lower > F(cutoffs["min_weak_curvature_signal"]),
            "translated_bump_after_measurement": cfg["deformations"]["unobserved_center_t"]
            - time_radius > cfg["tube_time"][1],
        },
    }


def run(definition, frozen):
    spring = check_f1(definition, frozen)
    stress = check_stress(definition["spring"], definition["compute"]["max_sites"])
    continuum = check_continuum(definition["continuum"])
    passed = all(all(case["checks"].values()) for case in (spring, stress, continuum))
    return {"protocol": definition["protocol"], "freeze": FREEZE, "f1_pair": spring,
            "singleton_obstruction": stress, "continuum_certificates": continuum, "passed": passed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--f1-input", type=Path)
    args = parser.parse_args()
    definition = json.loads((HERE / "PARAMETERS.json").read_text())
    parent = args.f1_input or (HERE / definition["f1_inputs_path"])
    source = parent.read_bytes()
    if sha256(source).hexdigest() != definition["f1_inputs_sha256"]:
        raise ValueError("The frozen F1 spring input did not match its content hash")
    outcome = run(definition, json.loads(source))
    args.output.write_text(json.dumps(outcome, indent=2, default=str, allow_nan=False) + "\n")
    print("Pulse controls, exact reconstruction and continuum bounds:", outcome["passed"])
    if not outcome["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
