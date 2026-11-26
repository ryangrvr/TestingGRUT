#!/usr/bin/env python3
"""Enumerate individual fuel molecules and independently propagate by uniformization."""

import argparse
from collections import defaultdict
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import poisson


HERE = Path(__file__).resolve().parent


def bond_count(x):
    return sum(x[i - 1] == 1 and x[i] == 1 for i in range(len(x)))


def replace_fuel(f, j):
    return f[:j] + (1 - f[j],) + f[j + 1:]


def transitions(state, cfg, spec):
    x, e, m, fuel = state
    yield (x, 1 - e, m, fuel), cfg["condition_flip_rate"], None, 0
    yield (x, e, 1 - m, fuel), cfg["record_thermal_flip_rate"], None, 0
    k = len(fuel)
    a = spec["a"][e if spec["actor"] == "condition" else m]
    j0, gap = cfg["J"], cfg["A"]
    for i in range(len(x)):
        j = (i + 1) % len(x)
        if x[i] == x[j]:
            continue
        moved = list(x)
        moved[i], moved[j] = moved[j], moved[i]
        moved = tuple(moved)
        d = bond_count(moved) - bond_count(x)
        barrier = cfg["patch_rates"][e] if i == 0 else cfg["background_rate"]
        yield (moved, e, m, fuel), barrier * math.exp(j0 * d / 2), None, 0
        if i == 0 and d != 0 and a:
            for u, f in enumerate(fuel):
                if (d > 0) == bool(f):
                    yield ((moved, e, m, replace_fuel(fuel, u)),
                           a / k * math.exp((j0 + gap) * d / 2), 0, d)
    sgn = 1 if m != e else -1
    if spec["s"]:
        for u, f in enumerate(fuel):
            if (sgn > 0) == bool(f):
                yield ((x, e, 1 - m, replace_fuel(fuel, u)),
                       spec["s"] / k * math.exp(gap * sgn / 2), 1, sgn)


def grow(q, t, initial, tail_tolerance, speed):
    """Positive uniformization standardizes probability tails for killed blocks."""
    embedded = np.eye(q.shape[0]) + q / speed
    expectation = speed * t
    cutoff = math.ceil(expectation)
    # The first Poisson-tail moment controls a linear reward accumulator.
    while max(poisson.sf(cutoff, expectation),
              expectation * poisson.sf(cutoff - 1, expectation)) > tail_tolerance:
        cutoff += 1
    weight = math.exp(-expectation)
    result = weight * initial
    step = initial.copy()
    for j in range(1, cutoff + 1):
        step = embedded @ step
        weight *= expectation / j
        result += weight * step
    return result, cutoff, float(poisson.sf(cutoff, expectation))


def marginals(states, p):
    groups = [defaultdict(float) for _ in range(4)]
    for state, mass in zip(states, p, strict=True):
        for index, value in enumerate(state):
            groups[index][value] += mass
    return groups


def entropy(xs):
    return -math.fsum(p * math.log(p) for p in xs if p > 0)


def relative(p, q):
    return math.fsum(x * math.log(x / y) for x, y in zip(p, q, strict=True) if x > 0)


def density_row(cfg, spec):
    size, k = (cfg["independent_grid"][c] for c in ("L", "K"))
    particles = size // 2
    gas = [x for x in itertools.product((0, 1), repeat=size) if sum(x) == particles]
    states = list(itertools.product(gas, (0, 1), (0, 1),
                                    tuple(itertools.product((0, 1), repeat=k))))
    rank = {state: i for i, state in enumerate(states)}
    n = len(states)
    aug = np.zeros((n + 2, n + 2))
    for j, state in enumerate(states):
        for successor, rate, channel, consumed in transitions(state, cfg, spec):
            aug[rank[successor], j] += rate
            aug[j, j] -= rate
            if channel is not None:
                aug[n + channel, j] += rate * consumed
    bonds = np.array([bond_count(x) for x, _, _, _ in states])
    fuel = np.array([sum(f) for _, _, _, f in states])
    keep = particles - bonds <= particles // 2
    init = (fuel == k).astype(float)
    if spec["start"] == "organized":
        init *= keep * np.exp(cfg["J"] * bonds)
    init /= init.sum()
    t, tol = cfg["duration"], cfg["thresholds"]["uniformization_tail"]
    speed = max(-aug.diagonal())
    outcome, terms, lower_tail = grow(aug, t, np.append(init, [0, 0]), tol, speed)
    p = outcome[:n]
    survival, _, killed_tail = grow(aug[:n, :n][np.ix_(keep, keep)], t, init[keep], tol, speed)
    energy = -cfg["J"] * bonds + cfg["A"] * fuel
    pi = np.exp(-energy)
    pi /= pi.sum()
    start_groups, end_groups, equilibrium = (marginals(states, v) for v in (init, p, pi))
    separate = math.fsum(relative(u.values(), e.values()) - relative(v.values(), e.values())
                         for u, v, e in zip(start_groups, end_groups, equilibrium, strict=True))
    correlation = sum(entropy(g.values()) for g in end_groups) - entropy(p)
    em = np.zeros((2, 2))
    for (_, e, m, _), value in zip(states, p, strict=True):
        em[e, m] += value
    past_information = None
    if spec["start"] == "organized":
        q = aug[:n, :n]
        mark = cfg["record_tag_time"]
        mid, _, _ = grow(q, mark, init, tol, speed)
        conditions = np.array([state[1] for state in states])
        pair = np.zeros((2, 2, 2))
        for old in (0, 1):
            propagated, _, _ = grow(q, t - mark, mid * (conditions == old), tol, speed)
            for (_, e, m, _), value in zip(states, propagated, strict=True):
                pair[old, e, m] += value
        past_information = (entropy(pair.sum(axis=2).ravel()) + entropy(pair.sum(axis=0).ravel())
                            - entropy(pair.sum(axis=(0, 2))) - entropy(pair.ravel()))
    return {
        "name": spec["name"], "labelled_dimension": n, "uniformization_terms": terms + 1,
        "probability_tail_bound": max(lower_tail, killed_tail),
        "first_exit_survival": float(survival.sum()), "end_organized": float(p[keep].sum()),
        "end_bond_fraction": float(p @ bonds / particles),
        "joint_reset": relative(init, pi) - relative(p, pi),
        "separate_reset": separate,
        "heat_to_solvent": float((init - p) @ energy),
        "final_total_correlation_nats": correlation,
        "fuel_used_repair": float(outcome[n]),
        "fuel_used_sensor": float(outcome[n + 1]),
        "fuel_remaining": float(p @ fuel),
        "fuel_distribution": np.bincount(fuel, weights=p, minlength=k + 1).tolist(),
        "current_record_information_nats": (entropy(em.sum(axis=0)) + entropy(em.sum(axis=1))
                                            - entropy(em.ravel())),
        "conditional_past_information_nats": past_information,
        "mass_residual": abs(float(p.sum()) - 1),
    }


def compare(rows, count, tolerance):
    expected = {r["name"]: r for r in count}
    if {r["name"] for r in rows} != set(expected) or len(expected) != len(rows):
        raise AssertionError("Reference and oracle grids differ.")
    differences = []
    for row in rows:
        target = expected[row["name"]]
        for metric in ("first_exit_survival", "end_organized", "end_bond_fraction",
                       "current_record_information_nats"):
            differences.append(abs(row[metric] - target[metric]))
        if row["conditional_past_information_nats"] is not None:
            differences.append(abs(row["conditional_past_information_nats"]
                                   - target["conditional_past_information_nats"]))
        for metric in ("joint_reset", "separate_reset", "heat_to_solvent",
                       "final_total_correlation_nats", "fuel_used_repair", "fuel_used_sensor",
                       "fuel_remaining"):
            differences.append(abs(row[metric] - target["ledger"][metric]))
        for a, b in zip(row["fuel_distribution"], target["ledger"]["fuel_distribution"], strict=True):
            differences.append(abs(a - b))
    error = max(differences)
    if error > tolerance:
        raise AssertionError(error)
    return {"scalar_comparisons": len(differences), "maximum_discrepancy": error}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--output", type=Path)
    output.add_argument("--verify", type=Path)
    args = parser.parse_args()
    cfg = json.loads((HERE / "PARAMETERS.json").read_text())
    rows = [density_row(cfg, spec) for spec in cfg["independent_grid"]["runs"]]
    primary = json.loads((HERE / "RESULTS.json").read_text())
    diagnostics = compare(rows, primary["independent_grid_count_solver"],
                          cfg["thresholds"]["numerical_abs"])
    result = {"method": "individual-fuel Markov chain with Poisson-uniformization",
              "rows": rows, "comparison": diagnostics}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads(args.verify.read_text()):
        raise SystemExit("Stored oracle output differs from replay.")
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
