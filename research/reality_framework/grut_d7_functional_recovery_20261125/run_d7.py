#!/usr/bin/env python3
"""Evaluate the frozen D7 grid by exact ordered-path reduction."""

import argparse
import cmath
from collections import defaultdict
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
THETA = {"none": 0.0, "weak": math.pi / 12, "strong": math.pi / 4}
PHI = {"shallow": math.pi / 8, "deep": math.pi / 4}
COLUMNS = [
    "length", "assay", "severity", "port", "intervention", "policy", "source",
    "beta", "raw_record_correct", "latch_correct", "endpoint_bell",
    "left_receiver_bell", "right_receiver_bell", "recurrence_bell",
]


def angle(theta, phi, policy):
    if policy == "blind":
        return 0.0
    if policy == "full":
        return phi
    if policy != "calibrated":
        raise ValueError(policy)
    cosine = 0.0 if phi == math.pi / 4 else math.cos(2 * phi)
    return math.atan2(math.sin(2 * theta) * math.sin(2 * phi), cosine) / 2


def acquired_probabilities(theta, source):
    a = cmath.exp(-1j * (1 - 2 * source) * theta)
    v = (a - 1j * a.conjugate(), a.conjugate() - 1j * a)
    weights = tuple(abs(x) ** 2 for x in v)
    return tuple(w / math.fsum(weights) for w in weights)


def select_pointer(raw, passive_coin, mask_coin, source, port, intervention):
    if port == "direct":
        return source
    if port == "passive":
        return passive_coin
    if port != "main":
        raise ValueError(port)
    return {
        "intact": raw,
        "scrambled": raw ^ mask_coin,
        "inverted": raw ^ 1,
        "disabled": 0,
    }[intervention]


def endpoint(theta, phi, beta, source, port, intervention):
    pointer = acquired_probabilities(theta, source)
    latch = [0.0, 0.0]
    for raw, cp, cm in itertools.product((0, 1), repeat=3):
        decision = select_pointer(raw, cp, cm, source, port, intervention)
        latch[decision] += pointer[raw] / 4
    fidelity = math.fsum(
        p * math.cos((1 - 2 * source) * phi - (1 - 2 * bit) * beta) ** 2
        for bit, p in enumerate(latch)
    )
    return pointer[source], latch[source], fidelity


def clock_clearance(length, offset, first_ready):
    miss = math.sin(offset / 2) ** 2
    early = math.fsum(
        math.comb(length, j) * (1 - miss) ** j * miss ** (length - j)
        for j in range(first_ready)
    )
    return 1 - early, early


def direct_receiver(fidelity, clock):
    return 0.25 + clock * (fidelity - 0.25)


def programme_optima(theta, phi):
    states = []
    pointers = []
    for source in (0, 1):
        a = cmath.exp(-1j * (1 - 2 * source) * theta)
        vector = np.array([a, a.conjugate()]) / math.sqrt(2)
        states.append(np.outer(vector, vector.conj()))
        pointers.append(np.diag(acquired_probabilities(theta, source)))
    phase = cmath.exp(-2j * phi)

    def optimal(pair):
        weighted = (phase * pair[0] + phase.conjugate() * pair[1]) / 2
        return (1 + float(np.linalg.svd(weighted, compute_uv=False).sum())) / 2

    return optimal(states), optimal(pointers)


def calculate():
    spec = json.loads((HERE / "PROTOCOL.json").read_text())
    grid, gates, timing = spec["production_grid"], spec["gates"], spec["timing"]
    left, right = timing["window_offsets_from_pi"]
    first_ready = timing["first_ready_site"]
    clocks = {}
    for length in grid["length"]:
        ready_l, early_l = clock_clearance(length, left, first_ready)
        ready_r, early_r = clock_clearance(length, right, first_ready)
        clocks[length] = {
            "T": math.pi * (length + 1) / 4,
            "ready_left": ready_l,
            "ready_right": ready_r,
            "incomplete_left": early_l,
            "incomplete_right": early_r,
            "incomplete_bound": min(1.0, (2 ** length) * (0.01 ** (length - first_ready + 1))),
            "energy_mean": 0.0,
            "energy_variance": 4 * length / (length + 1) ** 2,
        }

    rows = []
    averages = defaultdict(list)
    for length, assay, severity, arm, policy, source in itertools.product(
        grid["length"], grid["assay"], grid["severity"], grid["arm"],
        grid["policy"], grid["source"]
    ):
        port, intervention = arm
        theta, phi = THETA[assay], PHI[severity]
        beta = angle(theta, phi, policy)
        raw, latch, fidelity = endpoint(theta, phi, beta, source, port, intervention)
        cl = clocks[length]
        row = [length, assay, severity, port, intervention, policy, source, beta, raw, latch,
               fidelity, direct_receiver(fidelity, cl["ready_left"]),
               direct_receiver(fidelity, cl["ready_right"]), 0.25]
        rows.append(row)
        averages[tuple(row[:6])].append(fidelity)

    checks = []
    errors = []
    benchmark = {}
    for assay, severity in itertools.product(grid["assay"], grid["severity"]):
        theta, phi = THETA[assay], PHI[severity]
        f_svd, f_pointer = programme_optima(theta, phi)
        f_exact = (1 + math.hypot(math.cos(2 * phi),
                                  math.sin(2 * theta) * math.sin(2 * phi))) / 2
        f_blind = (1 + abs(math.cos(2 * phi))) / 2
        errors += [abs(f_svd - f_exact), abs(f_pointer - f_exact), abs(f_svd - f_pointer)]
        benchmark[f"{assay}_{severity}"] = {
            "coherent_optimum": f_svd,
            "pointer_optimum": f_pointer,
            "closed_form_optimum": f_exact,
            "blind_optimum": f_blind,
        }

    def observed(length, assay, severity, policy="calibrated", arm=("main", "intact")):
        values = averages[(length, assay, severity, *arm, policy)]
        return math.fsum(values) / len(values)

    for length, assay, severity in itertools.product(
        grid["length"], grid["assay"], grid["severity"]
    ):
        b = benchmark[f"{assay}_{severity}"]
        errors += [
            abs(observed(length, assay, severity) - b["closed_form_optimum"]),
            abs(observed(length, assay, severity, "blind") - b["blind_optimum"]),
            abs(observed(length, assay, severity, "full", ("direct", "disabled")) - 1),
        ]

    weak_summary = {}
    for length in grid["length"]:
        soft = observed(length, "weak", "shallow")
        blind = observed(length, "weak", "shallow", "blind")
        full = observed(length, "weak", "shallow", "full")
        scramble = observed(length, "weak", "shallow", arm=("main", "scrambled"))
        inverted = observed(length, "weak", "shallow", arm=("main", "inverted"))
        errors += [
            abs(soft - (1 + math.sqrt(5 / 8)) / 2),
            abs(blind - (1 + 1 / math.sqrt(2)) / 2),
            abs(full - spec["exact_targets"]["weak_shallow_full"]),
            abs(observed(length, "weak", "deep") - spec["exact_targets"]["weak_deep_calibrated"]),
            abs(observed(length, "weak", "deep", "blind") - spec["exact_targets"]["weak_deep_blind"]),
        ]
        cl = clocks[length]
        ready = min(cl["ready_left"], cl["ready_right"])
        late_gain = ready * (soft - blind)
        checks += [
            (f"blind_advantage_L{length}", soft - blind >=
             gates["weak_shallow_calibrated_vs_optimal_blind_gain_min"]),
            (f"full_advantage_L{length}", soft - full >=
             gates["weak_shallow_calibrated_vs_full_gain_min"]),
            (f"scramble_reduction_L{length}", scramble <= gates["weak_shallow_same_angle_scramble_upper"]),
            (f"inversion_reduction_L{length}", inverted <= gates["weak_shallow_same_angle_inversion_upper"]),
        ]
        if length in gates["expected_ready_diagnostic_fail_length"]:
            checks.append((f"short_window_miss_L{length}", ready < gates["ready_diagnostic_threshold"]))
        if length in gates["expected_ready_diagnostic_pass_length"]:
            checks += [
                (f"padded_ready_L{length}", ready >= gates["ready_diagnostic_threshold"]),
                (f"tail_L{length}", cl["incomplete_bound"] <=
                 gates["padded_incompletion_bound_max"] and
                 max(cl["incomplete_left"], cl["incomplete_right"]) <= cl["incomplete_bound"]),
                (f"late_advantage_L{length}", late_gain >=
                 gates["weak_shallow_late_gain_min_for_length_at_least_16"]),
            ]
        weak_summary[str(length)] = dict(calibrated=soft, blind=blind, full=full, scrambled=scramble,
                                         inverted=inverted, gain_over_blind=soft-blind,
                                         min_window_gain=late_gain)

    errors.extend(abs(row[-1] - spec["exact_targets"]["blind_replacement"]) for row in rows)
    checks.append(("forecast_residual", max(errors) <= gates["analytic_residual_max"]))
    failed = [name for name, passed in checks if not passed]
    if failed:
        raise AssertionError(failed)

    return {
        "model": "ordered spin clock, Bell observable and independent programme trace norms",
        "columns": COLUMNS,
        "checks": dict(checks),
        "max_abs_analytic_residual": max(errors),
        "clock": clocks,
        "benchmarks": benchmark,
        "weak_shallow": weak_summary,
        "rows": rows,
    }


def serialize(result):
    metadata = json.dumps({k: v for k, v in result.items() if k != "rows"},
                          indent=2, allow_nan=False)
    rows = ",\n".join("    " + json.dumps(r, allow_nan=False) for r in result["rows"])
    return metadata[:-2] + ',\n  "rows": [\n' + rows + '\n  ]\n}\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--output", type=Path)
    output.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = calculate()
    serialized = serialize(result)
    if args.output:
        args.output.write_text(serialized)
    elif json.loads(args.verify.read_text()) != json.loads(serialized):
        raise SystemExit("Stored data differ from replay")
    print(f"{len(result['rows'])} reduced rows, residual {result['max_abs_analytic_residual']:.3g}")
    print(json.dumps(result["weak_shallow"], indent=2))


if __name__ == "__main__":
    main()
