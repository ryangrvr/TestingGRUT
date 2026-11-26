#!/usr/bin/env python3
"""Independent Fraction replay of the frozen F1 resolvent witnesses."""
from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import signal


HERE = Path(__file__).resolve().parent


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def stiffness(n, edges, pinning):
    matrix = [[pinning * x for x in row] for row in identity(n)]
    for a, b, k in edges:
        matrix[a][a] += k
        matrix[b][b] += k
        matrix[a][b] -= k
        matrix[b][a] -= k
    return matrix


def shift(matrix, s):
    return [[x + s * int(i == j) for j, x in enumerate(row)]
            for i, row in enumerate(matrix)]


def exact_inverse_determinant(matrix):
    n = len(matrix)
    augmented = [[Fraction(x) for x in row + unit]
                 for row, unit in zip(matrix, identity(n), strict=True)]
    det = Fraction(1)
    for column in range(n):
        pivot_row = next(i for i in range(column, n) if augmented[i][column])
        if pivot_row != column:
            augmented[column], augmented[pivot_row] = augmented[pivot_row], augmented[column]
            det *= -1
        pivot = augmented[column][column]
        det *= pivot
        normalized = [x / pivot for x in augmented[column]]
        augmented[column] = normalized
        for i, values in enumerate(augmented):
            if i == column or values[column] == 0:
                continue
            multiplier = values[column]
            augmented[i] = [x - multiplier * y
                            for x, y in zip(values, normalized, strict=True)]
    return [row[n:] for row in augmented], det


def exact_rank(columns):
    array = [list(row) for row in zip(*columns, strict=True)]
    height, width = len(array), len(array[0])
    rank = 0
    for j in range(width):
        i = next((i for i in range(rank, height) if array[i][j]), None)
        if i is None:
            continue
        array[rank], array[i] = array[i], array[rank]
        pivot = array[rank][j]
        for i in range(rank + 1, height):
            if not array[i][j]:
                continue
            quotient = array[i][j] / pivot
            array[i][j:] = [x - quotient * y
                            for x, y in zip(array[i][j:], array[rank][j:], strict=True)]
        rank += 1
        if rank == height:
            break
    return rank


def multiply(a, b):
    columns = list(zip(*b, strict=True))
    return [[sum(x * y for x, y in zip(row, col, strict=True))
             for col in columns] for row in a]


def graph_data(n, edges):
    neighborhood = [set() for _ in range(n)]
    weighted_degrees = [0] * n
    for a, b, weight in edges:
        assert 0 <= a < n and 0 <= b < n and a != b and weight > 0
        assert b not in neighborhood[a]
        neighborhood[a].add(b)
        neighborhood[b].add(a)
        weighted_degrees[a] += weight
        weighted_degrees[b] += weight
    visited = {0}
    front = [0]
    while front:
        unseen = neighborhood[front.pop()] - visited
        visited.update(unseen)
        front.extend(unseen)
    triangle_count = sum(b in neighborhood[a] and c in neighborhood[b] and c in neighborhood[a]
                         for a, b, c in itertools.combinations(range(n), 3))
    return {
        "connected": len(visited) == n,
        "edges": len(edges),
        "weight_sum": sum(weight for _, _, weight in edges),
        "triangle_count": triangle_count,
        "degrees_sorted": sorted(map(len, neighborhood)),
        "weighted_degree_max": max(weighted_degrees),
    }


def mix_for_case(n, frozen):
    mix = identity(n)
    numerator = frozen["numerator"]
    denominator = frozen["denominator"]
    for i, row in enumerate(numerator, start=1):
        for j, element in enumerate(row, start=1):
            mix[i][j] = Fraction(element, denominator)
    return mix


def compare_resolvents(v, w, ports):
    n = len(v)
    vectors = ([], [])
    dets = ([], [])
    equal_port_samples = 0
    equal_spectra_samples = 0
    for s in range(1, 2 * n + 1):
        rv, dv = exact_inverse_determinant(shift(v, s))
        rw, dw = exact_inverse_determinant(shift(w, s))
        assert dv > 0 and dw > 0
        equal_port_samples += all(rv[a][b] == rw[a][b] for a in ports for b in ports)
        equal_spectra_samples += dv == dw
        for columns, inverse in zip(vectors, (rv, rw), strict=True):
            columns.extend([[inverse[i][port] for i in range(n)] for port in ports])
        dets[0].append(dv)
        dets[1].append(dw)
    return {
        "samples": 2 * n,
        "matrix_entries_compared": 2 * n * len(ports) ** 2,
        "equal_port_samples": equal_port_samples,
        "equal_spectral_determinant_samples": equal_spectra_samples,
        "rank": [exact_rank(columns) for columns in vectors],
        "determinants_at_shift_1": [str(det[0]) for det in dets],
    }, vectors, dets


def krylov_determinant_from_inverses(columns, dets, n):
    # The coefficient vectors of adj(V+sI)e_0 and the Krylov vectors
    # have a triangular basis change with unit absolute determinant.
    _, det_rational = exact_inverse_determinant(
        [list(row) for row in zip(*columns[:n], strict=True)]
    )
    vandermonde = math.prod(j - i for i, j in itertools.combinations(range(1, n + 1), 2))
    return abs(det_rational) * math.prod(dets[:n]) / vandermonde


def added_probe_bound(spec, matrices, graph_max_degree):
    first, last = spec["observe"]
    order = spec["derivative_order"]
    assert order == 3 and first != last
    derivatives = {name: -value[first][last] for name, value in matrices.items()}
    assert derivatives == spec["derivatives"]
    diff = abs(derivatives["cycle"] - derivatives["pendant_triangle"])
    t = Fraction(spec["time"])
    norm_upper = spec["stiffness_operator_norm_upper"]
    assert norm_upper >= 1 + 2 * graph_max_degree
    # For j >= 2 the term ratio decreases from the one at j=2.
    ratio = norm_upper * t**2 / (6 * 7)
    assert 0 <= ratio < 1
    tail = 2 * norm_upper**2 * t**5 / (math.factorial(5) * (1 - ratio))
    leading = diff * t**3 / math.factorial(order)
    lower = leading - tail
    required = Fraction(spec["required_response_gap_lower"])
    assert lower > required
    return {
        "derivatives": derivatives,
        "absolute_derivative_difference": diff,
        "leading_gap": str(leading),
        "tail_ratio_upper": str(ratio),
        "tail_upper": str(tail),
        "gap_lower": str(lower),
        "required_gap_lower": str(required),
        "excess_over_required": str(lower - required),
    }


def run(spec):
    spring_graphs = spec["networks"]
    pinning = spec["units"]["common_pinning_frequency_squared"]
    max_sites = spec["exact_compute"]["maximum_sites"]
    assert spec["units"]["mass"] == pinning == 1
    cases = []
    primary = None
    for case in spec["cases"]:
        n, ports = case["sites"], case["observe"]
        assert n <= max_sites
        edges = {name: values + case["append_edges"] for name, values in spring_graphs.items()}
        graphs = {name: graph_data(n, values) for name, values in edges.items()}
        v, w = matrices = tuple(stiffness(n, edge, pinning) for edge in edges.values())
        mix = mix_for_case(n, spec["hidden_mixing"])
        mt = list(map(list, zip(*mix, strict=True)))
        assert multiply(mt, mix) == identity(n)
        assert [sum(row) for row in mix] == [1] * n
        assert [mix[i] for i in ports] == [identity(n)[i] for i in ports]
        assert [sum(row) for row in v] == [pinning] * n
        assert [sum(row) for row in w] == [pinning] * n
        assert multiply(multiply(mt, v), mix) == w
        graph_left, graph_right = graphs.values()
        assert graph_left["connected"] and graph_right["connected"]
        assert graph_left["edges"] == graph_right["edges"]
        assert graph_left["weight_sum"] == graph_right["weight_sum"]
        assert graph_right["triangle_count"] - graph_left["triangle_count"] == 1
        comparison, columns, dets = compare_resolvents(v, w, ports)
        assert comparison["equal_port_samples"] == 2 * n
        assert comparison["equal_spectral_determinant_samples"] == 2 * n
        assert comparison["rank"] == [case["rank"]] * 2
        result = dict(name=case["name"], sites=n, ports=ports,
                      graphs=graphs, orthogonal_conjugacy=True, **comparison)
        if case["name"] == "one_port":
            primary = dict(zip(spring_graphs, matrices, strict=True))
            assert graph_left["edges"] == spec["primary_edge_count"]
            assert graph_left["weight_sum"] == spec["primary_spring_sum"]
            assert graph_left["degrees_sorted"] != graph_right["degrees_sorted"]
            determinants = [krylov_determinant_from_inverses(col, det, n)
                            for col, det in zip(columns, dets, strict=True)]
            assert determinants == [spec["primary_krylov_determinant_magnitude"]] * 2
            result["krylov_determinant_magnitudes"] = [str(d) for d in determinants]
        if case["name"] == "connected_dark_mode_control":
            dark = [[0], [0], [0], [0], [1], [-1]]
            assert multiply(v, dark) == multiply(w, dark) == [[2 * d[0]] for d in dark]
            result["dark_mode_eigenvalue"] = 2
        cases.append(result)

    assert primary is not None
    v, w = primary.values()
    perm = spec["hidden_relabel_of_cycle"]
    assert perm[0] == 0 and sorted(perm) == list(range(len(v)))
    permuted = [[v[i][j] for j in perm] for i in perm]
    relabel, _, _ = compare_resolvents(v, permuted, [0])
    assert relabel["equal_port_samples"] == 2 * len(perm)
    assert relabel["rank"] == [len(v)] * 2
    probed, _, _ = compare_resolvents(v, w, spec["calibrated_extra_probe"]["observe"])
    assert probed["equal_port_samples"] < probed["samples"]
    gap = added_probe_bound(
        spec["calibrated_extra_probe"], primary,
        max(g["weighted_degree_max"] for g in cases[0]["graphs"].values())
    )
    return {
        "method": "Fraction inverse sampling and resolvent-column visibility",
        "cases": cases,
        "primary_case_samples": sum(row["samples"] for row in cases),
        "primary_case_matrix_entry_comparisons": sum(row["matrix_entries_compared"] for row in cases),
        "relabel": {"permutation": perm, "isomorphic": True, **relabel},
        "added_probe": gap | probed,
    }


def timeout(*_):
    raise TimeoutError("Frozen independent replay wall limit exceeded")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    spec = json.loads((HERE / "WITNESSES.json").read_text())
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(spec["exact_compute"]["wall_seconds_per_program"])
    output = run(spec)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    signal.alarm(0)
    print(json.dumps({k: v for k, v in output.items() if k != "cases"}, indent=2))


if __name__ == "__main__":
    main()
