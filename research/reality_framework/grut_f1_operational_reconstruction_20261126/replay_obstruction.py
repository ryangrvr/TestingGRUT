"""Verify the frozen retarded-response obstruction using exact moments."""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import factorial
from pathlib import Path


HERE = Path(__file__).resolve().parent
FREEZE = "5ffdc3c8717d2a8331643780e98c14fd031566be"


def identity(size):
    return [[int(row == col) for col in range(size)] for row in range(size)]


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def multiply(left, right):
    columns = transpose(right)
    return [[sum(a * b for a, b in zip(row, column)) for column in columns] for row in left]


def stiffness(size, edges, pin):
    matrix = [[pin * value for value in row] for row in identity(size)]
    seen = set()
    for i, j, weight in edges:
        endpoints = frozenset((i, j))
        if endpoints in seen or not 0 <= i < j < size or weight <= 0:
            raise ValueError(f"Invalid spring {(i, j, weight)}")
        seen.add(endpoints)
        matrix[i][i] += weight
        matrix[j][j] += weight
        matrix[i][j] -= weight
        matrix[j][i] -= weight
    return matrix


def powers(matrix, count):
    current = identity(len(matrix))
    result = [current]
    for _ in range(count - 1):
        current = multiply(current, matrix)
        result.append(current)
    return result


def submatrix(matrix, indices):
    return [[matrix[i][j] for j in indices] for i in indices]


def row_elimination(matrix):
    """Return rank and the signed product of pivots."""
    reduced = [[Fraction(entry) for entry in row] for row in matrix]
    pivot_product = Fraction(1)
    rank = 0
    for col in range(len(reduced[0])):
        pivot = next((row for row in range(rank, len(reduced)) if reduced[row][col]), None)
        if pivot is None:
            continue
        if pivot != rank:
            reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
            pivot_product = -pivot_product
        top = reduced[rank]
        pivot_product *= top[col]
        for row in reduced[rank + 1:]:
            factor = row[col] / top[col]
            if factor:
                for k in range(col, len(row)):
                    row[k] -= factor * top[k]
        rank += 1
        if rank == len(reduced):
            break
    return rank, pivot_product


def controllability(moment_powers, observed):
    size = len(moment_powers[0])
    span = [
        [power[i][j] for power in moment_powers[:size] for j in observed]
        for i in range(size)
    ]
    rank, pivot_product = row_elimination(span)
    determinant = int(pivot_product) if rank == size and len(observed) == 1 else None
    return rank, determinant


def graph_invariants(size, edges):
    adjacency = [set() for _ in range(size)]
    for i, j, _ in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    frontier = [0]
    reached = {0}
    while frontier:
        neighbors = adjacency[frontier.pop()] - reached
        reached.update(neighbors)
        frontier.extend(neighbors)
    triangles = sum(
        j in adjacency[i] and k in adjacency[i] and k in adjacency[j]
        for i, j, k in combinations(range(size), 3)
    )
    return {
        "connected": len(reached) == size,
        "triangles": triangles,
        "edge_count": len(edges),
        "total_spring": sum(weight for _, _, weight in edges),
    }


def hidden_reflection(size, specification):
    reflection = identity(size)
    for i, row in enumerate(specification["numerator"], start=1):
        for j, number in enumerate(row, start=1):
            reflection[i][j] = Fraction(number, specification["denominator"])
    return reflection


def evaluate_case(definition, case):
    size = case["sites"]
    observed = case["observe"]
    pin = definition["units"]["common_pinning_frequency_squared"]
    edges = {name: data + case["append_edges"] for name, data in definition["networks"].items()}
    forms = {name: stiffness(size, data, pin) for name, data in edges.items()}
    topology = {name: graph_invariants(size, data) for name, data in edges.items()}
    matrices = {name: powers(form, 2 * size) for name, form in forms.items()}
    moments = {name: [submatrix(power, observed) for power in seq] for name, seq in matrices.items()}
    reachability = {name: controllability(seq, observed) for name, seq in matrices.items()}
    cycle, triangle = forms["cycle"], forms["pendant_triangle"]
    change_of_basis = hidden_reflection(size, definition["hidden_mixing"])
    adjoint = transpose(change_of_basis)
    q_unit = multiply(change_of_basis, adjoint)
    transformed = multiply(multiply(change_of_basis, cycle), adjoint)
    resource_keys = ("edge_count", "total_spring")
    checks = {
        "orthogonal": q_unit == identity(size),
        "fixes_measured_vectors": all(change_of_basis[i] == identity(size)[i] for i in observed),
        "fixes_constant_vector": all(sum(row) == 1 for row in change_of_basis),
        "conjugate_stiffness": transformed == triangle,
        "common_pinning": all(sum(row) == pin for form in forms.values() for row in form),
        "port_moments_equal": moments["cycle"] == moments["pendant_triangle"],
        "rank_as_registered": all(rank == case["rank"] for rank, _ in reachability.values()),
        "both_connected": all(graph["connected"] for graph in topology.values()),
        "different_triangle_counts": topology["cycle"]["triangles"] == 0
        and topology["pendant_triangle"]["triangles"] == 1,
        "resource_counts_match": all(topology["cycle"][key] == topology["pendant_triangle"][key]
                                     for key in resource_keys),
    }
    return {
        "name": case["name"],
        "sites": size,
        "observed": observed,
        "moment_count": 2 * size,
        "topology": topology,
        "rank": {name: rank for name, (rank, _) in reachability.items()},
        "all_modes_visible": all(rank == size for rank, _ in reachability.values()),
        "checks": checks,
    }, forms, reachability, moments


def separation(definition, forms):
    extra = definition["calibrated_extra_probe"]
    a, b = extra["observe"]
    derivatives = {name: -matrix[a][b] for name, matrix in forms.items()}
    norm = {name: max(sum(map(abs, row)) for row in matrix) for name, matrix in forms.items()}
    time = Fraction(extra["time"])
    bound = extra["stiffness_operator_norm_upper"]
    ratio = bound * time * time / (6 * 7)
    tail = bound**2 * time**5 / factorial(5) / (1 - ratio)
    difference = derivatives["pendant_triangle"] - derivatives["cycle"]
    leading = difference * time**3 / factorial(3)
    lower, upper = leading - 2 * tail, leading + 2 * tail
    return {
        "observed_added": extra["observe"],
        "third_derivative": derivatives,
        "third_derivative_gap": difference,
        "row_norms": norm,
        "time": str(time),
        "response_gap_interval": [str(lower), str(upper)],
        "checks": {
            "third_derivatives_as_registered": derivatives == extra["derivatives"],
            "valid_tail_bound": all(value <= bound for value in norm.values()) and ratio < 1,
            "gap_above_registered_floor": lower > Fraction(extra["required_response_gap_lower"]),
        },
    }


def run(definition):
    summary = []
    for case in definition["cases"]:
        if case["sites"] > definition["exact_compute"]["maximum_sites"]:
            raise ValueError("Site cap exceeded")
        result, forms, spaces, moments = evaluate_case(definition, case)
        summary.append(result)
        if case["name"] == "one_port":
            primary = result
            primary_forms = forms
            determinants = {name: determinant for name, (_, determinant) in spaces.items()}
            primary_moments = [matrix[0][0] for matrix in moments["cycle"]]

    extra = separation(definition, primary_forms)
    relabel = definition["hidden_relabel_of_cycle"]
    cycle = primary_forms["cycle"]
    relabeled = submatrix(cycle, relabel)
    relabeled_moments = [power[0][0] for power in powers(relabeled, len(primary_moments))]
    primary_checks = {
        "determinants_as_registered": all(abs(d) == definition["primary_krylov_determinant_magnitude"]
                                          for d in determinants.values()),
        "registered_edge_count_and_spring_sum": all(
            graph["edge_count"] == definition["primary_edge_count"]
            and graph["total_spring"] == definition["primary_spring_sum"]
            for graph in primary["topology"].values()
        ),
        "relabel_is_nontrivial_and_fixes_port": relabel[0] == 0
        and sorted(relabel) == list(range(len(cycle))) and cycle != relabeled,
        "relabel_port_moments_equal": relabeled_moments == primary_moments,
    }
    passed = all(primary_checks.values()) and all(extra["checks"].values())
    passed = passed and all(all(result["checks"].values()) for result in summary)
    return {
        "protocol": "GRUT-F1",
        "freeze": FREEZE,
        "cases": summary,
        "primary": {
            "port_moments": primary_moments,
            "krylov_determinants": determinants,
            "checks": primary_checks,
        },
        "calibrated_extra_probe": extra,
        "passed": passed,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    definition = json.loads((HERE / "WITNESSES.json").read_text())
    result = run(definition)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Verified {len(result['cases'])} pairs and calibrated extra probe: {result['passed']}")
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
