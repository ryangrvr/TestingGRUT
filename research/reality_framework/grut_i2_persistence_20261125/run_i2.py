#!/usr/bin/env python3
"""Frozen aggregate persistence and finite-state thermodynamic frontier."""

import argparse
from dataclasses import dataclass
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply


HERE = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Policy:
    name: str
    a: tuple
    sensor: float
    actor: str = "record"


@dataclass(frozen=True)
class Gas:
    length: int
    bonds: np.ndarray
    longest: np.ndarray
    organized: np.ndarray
    hops: tuple


@lru_cache(maxsize=None)
def gas_states(length):
    particles = length // 2
    words = list(itertools.combinations(range(length), particles))
    masks = [sum(1 << j for j in word) for word in words]
    inverse = {mask: i for i, mask in enumerate(masks)}
    bits = np.array([[(mask >> j) & 1 for j in range(length)] for mask in masks])
    bonds = np.sum(bits * np.roll(bits, -1, axis=1), axis=1)
    longest = []
    for row in bits:
        doubled = "".join(map(str, np.tile(row, 2)))
        longest.append(max(map(len, doubled.split("0"))))
    hops = []
    for i, mask in enumerate(masks):
        for j in range(length):
            other = (j + 1) % length
            if bits[i, j] != bits[i, other]:
                target = inverse[mask ^ (1 << j) ^ (1 << other)]
                hops.append((i, target, bonds[target] - bonds[i], int(j == 0)))
    return Gas(length, bonds, np.array(longest), particles - bonds <= particles // 2,
               tuple(np.asarray(hops, dtype=int).T))


def entropy(probability):
    p = np.asarray(probability).ravel()
    p = p[p > 0]
    return -float(np.dot(p, np.log(p)))


def relative_entropy(probability, equilibrium):
    p, q = probability.ravel(), equilibrium.ravel()
    mask = p > 0
    return float(np.dot(p[mask], np.log(p[mask] / q[mask])))


def marginals(p):
    return [p.sum(axis=tuple(j for j in range(4) if j != i)) for i in range(4)]


def equilibrium(gas, fuel, attraction, fuel_gap):
    px = np.exp(attraction * gas.bonds)
    px /= px.sum()
    logdeg = np.log([math.comb(fuel, k) for k in range(fuel + 1)])
    pf = np.exp(logdeg - fuel_gap * np.arange(fuel + 1))
    pf /= pf.sum()
    pi = px[:, None, None, None] * np.ones((1, 2, 2, 1)) * pf[None, None, None, :] / 4
    return pi, logdeg


def construct_generator(cfg, gas, fuel, attraction, policy, pi):
    index = np.arange(pi.size).reshape(pi.shape)
    triplets = [[], [], []]
    reward = np.zeros((2, pi.size))

    def add(source, target, rate, ledger=None, spend=0):
        source, target, rate = np.broadcast_arrays(source, target, rate)
        positive = rate > 0
        src, dst, amp = source[positive], target[positive], rate[positive]
        for collection, field in zip(triplets, (dst, src, amp), strict=True):
            collection.append(field.ravel())
        if ledger is not None:
            reward[ledger] += np.bincount(src, weights=spend * amp, minlength=pi.size)

    add(index, index[:, ::-1, :, :], cfg["condition_flip_rate"])
    add(index, index[:, :, ::-1, :], cfg["record_thermal_flip_rate"])
    xs, xt, db, patch = gas.hops
    for marker in (0, 1):
        mask = patch == marker
        barrier = np.array(cfg["patch_rates"] if marker else [cfg["background_rate"]] * 2)
        rate = np.exp(attraction * db[mask] / 2)[:, None, None, None] * barrier[None, :, None, None]
        add(index[xs[mask]], index[xt[mask]], rate)

    if fuel > 0:
        rates = np.array(policy.a)
        gate = rates[None, None, :, None] if policy.actor == "record" else rates[None, :, None, None]
        fraction = np.arange(1, fuel + 1) / fuel
        slices = {1: (slice(1, None), slice(None, -1), fraction),
                  -1: (slice(None, -1), slice(1, None), fraction[::-1])}
        for d, (stay, leave, availability) in slices.items():
            take = (patch == 1) & (db == d)
            rate = (math.exp((attraction + cfg["A"]) * d / 2) * gate
                    * availability[None, None, None, :])
            add(index[xs[take], :, :, stay], index[xt[take], :, :, leave], rate,
                ledger=0, spend=d)
        for e, m in itertools.product((0, 1), repeat=2):
            d = 1 if e != m else -1
            stay, leave, availability = slices[d]
            rate = policy.sensor * math.exp(cfg["A"] * d / 2) * availability[None, :]
            add(index[:, e, m, stay], index[:, e, 1 - m, leave], rate, ledger=1, spend=d)

    dst, src, values = [np.concatenate(v) for v in triplets]
    off = coo_matrix((values, (dst, src)), shape=(pi.size, pi.size)).tocsr()
    edges = off.tocoo()
    reverse = np.asarray(off[edges.col, edges.row]).ravel()
    logpi = np.log(pi.ravel())
    detailed_balance = np.max(abs(np.log(edges.data) - np.log(reverse)
                                  - logpi[edges.row] + logpi[edges.col]))
    generator = off - diags(np.asarray(off.sum(axis=0)).ravel())
    return generator.tocsr(), reward, float(detailed_balance)


def evolve(generator, t, state):
    if t == 0:
        return state.copy()
    return expm_multiply(generator * t, state, traceA=t * float(generator.diagonal().sum()))


def initialize(pi, gas, fuel, start):
    canonical = pi.sum(axis=(1, 2, 3))
    if start == "organized":
        gas_input = canonical * gas.organized
    elif start == "uniform":
        gas_input = np.ones_like(canonical)
    else:
        raise ValueError(start)
    gas_input /= gas_input.sum()
    initial = np.zeros_like(pi)
    initial[:, :, :, fuel] = gas_input[:, None, None] / 4
    return initial


def propagation(cfg, generator, reward, initial, gas, start):
    n = initial.size
    meta = bmat([[generator, csr_matrix((n, 2))],
                 [csr_matrix(reward), csr_matrix((2, 2))]], format="csr")
    seed = np.concatenate((initial.ravel(), [0.0, 0.0]))
    horizon = cfg["duration"]
    info_past = None
    if start == "organized":
        at = cfg["record_tag_time"]
        mid = evolve(meta, at, seed)
        cond = np.zeros((n + 2, 2))
        by_e = cond[:n].reshape(*initial.shape, 2)
        tagged_mid = mid[:n].reshape(initial.shape)
        for e in (0, 1):
            by_e[:, e, :, :, e] = tagged_mid[:, e, :, :]
        cond[n:, 0] = mid[n:]
        end = evolve(meta, horizon - at, cond)
        triple = end[:n].reshape(*initial.shape, 2).sum(axis=(0, 3))
        info_past = (entropy(triple.sum(axis=1)) + entropy(triple.sum(axis=2))
                     - entropy(triple.sum(axis=(1, 2))) - entropy(triple))
        final = end.sum(axis=1)
    else:
        final = evolve(meta, horizon, seed)
    keep = np.broadcast_to(gas.organized[:, None, None, None], initial.shape).ravel()
    killed = generator[keep, :][:, keep]
    survival = float(evolve(killed, horizon, initial.ravel()[keep]).sum())
    return final[:n].reshape(initial.shape), final[n:], survival, info_past


def thermodynamics(initial, end, pi, logdeg, gas, attraction, fuel_gap, reward):
    init_m, end_m, eq_m = map(marginals, (initial, end, pi))
    init_free = relative_entropy(initial, pi)
    end_free = relative_entropy(end, pi)
    p, q = init_m[0], end_m[0]
    n_start = float(np.dot(init_m[3], np.arange(len(logdeg))))
    n_end = float(np.dot(end_m[3], np.arange(len(logdeg))))
    energy_initial = float(-attraction * np.dot(p, gas.bonds) + fuel_gap * n_start)
    energy_final = float(-attraction * np.dot(q, gas.bonds) + fuel_gap * n_end)
    heat = energy_initial - energy_final
    bath_multiplicity = float(np.dot(end_m[3] - init_m[3], logdeg))
    micro_entropy = entropy(end) - entropy(initial) + bath_multiplicity
    total_corr = sum(entropy(v) for v in end_m) - entropy(end)
    reset = [relative_entropy(p0, pe) - relative_entropy(pt, pe)
             for p0, pt, pe in zip(init_m, end_m, eq_m, strict=True)]
    joint = init_free - end_free
    return {
        "initial_free_energy": init_free,
        "remaining_free_energy": end_free,
        "shape_preparation": relative_entropy(init_m[0], eq_m[0]),
        "stock_preparation": relative_entropy(init_m[3], eq_m[3]),
        "joint_reset": joint,
        "separate_reset": sum(reset),
        "marginal_reset": dict(zip(("gas", "condition", "record", "fuel"), reset, strict=True)),
        "final_total_correlation_nats": total_corr,
        "heat_to_solvent": heat,
        "micro_entropy_change_nats": micro_entropy,
        "entropy_production_nats": heat + micro_entropy,
        "count_medium_entropy_nats": heat + bath_multiplicity,
        "fuel_used_repair": float(reward[0]),
        "fuel_used_sensor": float(reward[1]),
        "fuel_remaining": n_end,
        "fuel_distribution": end_m[3].tolist(),
        "fuel_budget_residual": abs(n_start - n_end - float(reward.sum())),
        "relative_free_energy_residual": abs(heat + micro_entropy - joint),
        "correlation_cost_residual": abs(sum(reset) - joint - total_corr),
    }


def evaluate(cfg, length, fuel, attraction, policy, start="organized"):
    gas = gas_states(length)
    pi, degeneracy = equilibrium(gas, fuel, attraction, cfg["A"])
    initial = initialize(pi, gas, fuel, start)
    gen, flux, balance = construct_generator(cfg, gas, fuel, attraction, policy, pi)
    final, reward, survival, lagged = propagation(cfg, gen, flux, initial, gas, start)
    thermo = thermodynamics(initial, final, pi, degeneracy, gas, attraction, cfg["A"], reward)
    gas_final = final.sum(axis=(1, 2, 3))
    joint = final.sum(axis=(0, 3))
    return {
        "name": policy.name, "actor": policy.actor, "a": list(policy.a), "s": policy.sensor,
        "length": length, "fuel": fuel, "J": attraction, "start": start, "dimension": pi.size,
        "first_exit_survival": survival,
        "end_organized": float(gas_final[gas.organized].sum()),
        "initial_bond_fraction": float(np.dot(initial.sum(axis=(1, 2, 3)), gas.bonds) / (length // 2)),
        "end_bond_fraction": float(np.dot(gas_final, gas.bonds) / (length // 2)),
        "end_longest_fraction": float(np.dot(gas_final, gas.longest) / (length // 2)),
        "current_record_information_nats": entropy(joint.sum(axis=0)) + entropy(joint.sum(axis=1)) - entropy(joint),
        "conditional_past_information_nats": lagged,
        "ledger": thermo,
        "numerics": {
            "normalization": abs(float(final.sum()) - 1),
            "minimum_probability": float(final.min()),
            "generator_conservation": float(np.max(abs(np.asarray(gen.sum(axis=0))))),
            "detailed_balance": balance,
        },
    }


def policy_grid(cfg):
    pairs = tuple(itertools.product(cfg["response_levels"], repeat=2))
    for actor, sensors in (("record", [0, *cfg["sensor_levels"]]), ("condition", [0])):
        for sensor, a in itertools.product(sensors, pairs):
            yield Policy(f"{actor}:s={sensor:g}:a={a[0]:g},{a[1]:g}", a, sensor, actor)


def named(raw):
    return Policy(raw["name"], tuple(raw["a"]), raw["s"], raw.get("actor", "record"))


def canonical_scaling(length, attraction):
    n = length // 2
    counts = [length * math.comb(n - 1, clusters - 1) ** 2 // clusters
              for clusters in range(1, n + 1)]
    weights_log = np.array([math.log(c) + attraction * (n - j)
                            for j, c in enumerate(counts, start=1)])
    weights = np.exp(weights_log - weights_log.max())
    weights /= weights.sum()
    return {
        "length": length, "J": attraction, "exact_count_check": sum(counts) == math.comb(length, n),
        "mean_cluster_fraction": float(np.dot(weights, np.arange(1, n + 1) / n)),
        "one_cluster_probability": float(weights[0]),
        "shape_class_probability": float(weights[:n // 2].sum()),
        "asymptotic_cluster_fraction": 1 / (1 + math.exp(attraction / 2)),
    }


def point(row, work):
    return {"work": row["ledger"][work], "survival": row["first_exit_survival"],
            "mixture": [[row["name"], 1.0]]}


def interpolate(left, right, w, fixed, value):
    v = (1 - w) * left["survival"] + w * right["survival"]
    c = (1 - w) * left["work"] + w * right["work"]
    return {"work": value if fixed == "work" else c,
            "survival": value if fixed == "survival" else v,
            "mixture": [[left["mixture"][0][0], 1 - w], [right["mixture"][0][0], w]]}


def optimum(points, limit, fixed):
    if fixed == "survival":
        feasible = [p for p in points if p["survival"] >= limit]
        key = lambda p: p["work"]
    else:
        feasible = [p for p in points if p["work"] <= limit]
        key = lambda p: -p["survival"]
    pure = min(feasible, key=key, default=None)
    for p, q in itertools.combinations(points, 2):
        left, right = sorted((p, q), key=lambda t: t[fixed])
        if left[fixed] < limit < right[fixed]:
            w = (limit - left[fixed]) / (right[fixed] - left[fixed])
            feasible.append(interpolate(left, right, w, fixed, limit))
    return {"pure": pure, "mixture": min(feasible, key=key, default=None)}


def envelopes(cfg, primary):
    record = [r for r in primary if r["actor"] == "record"]
    off = [r for r in record if r["s"] == 0]
    families = {
        "passive": [r for r in off if max(r["a"]) == 0],
        "sensor_off": off,
        "no_feedback": [r for r in record if r["s"] == 0 or r["a"][0] == r["a"][1]],
        "full_record": record,
        "direct_cue": [r for r in primary if r["actor"] == "condition"],
    }
    result = {}
    for work in ("joint_reset", "separate_reset"):
        result[work] = {}
        for group, rows in families.items():
            points = [point(r, work) for r in rows]
            result[work][group] = {
                "at_cost": {str(budget): optimum(points, budget, "work") for budget in cfg["reset_budgets"]},
                "at_reliability": {str(goal): optimum(points, goal, "survival") for goal in cfg["reliabilities"]},
            }
    return result


def scientific_status(cfg, primary, flat, formation, frontiers):
    lookup = {r["name"]: r for r in primary}
    passive = next(r for r in primary if r["actor"] == "record" and r["s"] == 0 and max(r["a"]) == 0)
    grouped = {r["name"]: r for r in formation}
    threshold = cfg["thresholds"]
    gains = []
    for budget in cfg["reset_budgets"]:
        by_work = {}
        for work, groups in frontiers.items():
            available = groups["full_record"]["at_cost"][str(budget)]["mixture"]
            control = groups["no_feedback"]["at_cost"][str(budget)]["mixture"]
            active_past = max((lookup[name]["conditional_past_information_nats"]
                               for name, w in available["mixture"] if w > 0 and
                               lookup[name]["s"] > 0 and len(set(lookup[name]["a"])) > 1),
                              default=0.0) if available is not None else 0.0
            by_work[work] = {
                "gain": (available["survival"] - control["survival"]) if available and control else None,
                "active_record_past_nats": active_past,
            }
        qualifies = all(x["gain"] is not None and x["gain"] >= threshold["memory_budget_survival_gain"] and
                        x["active_record_past_nats"] >= threshold["memory_conditional_past_information_nats"]
                        for x in by_work.values())
        gains.append({"budget": budget, "comparison": by_work, "passes": qualifies})
    binding = grouped["cohesion"]["end_bond_fraction"] - grouped["cohesion"]["initial_bond_fraction"]
    fuel = grouped["constant_fuel"]["end_bond_fraction"] - grouped["cohesion"]["end_bond_fraction"]
    passive_gain = passive["first_exit_survival"] - flat["first_exit_survival"]
    return {
        "passive_attraction_survival_gain": passive_gain,
        "passive_formation_gain": binding,
        "driven_formation_gain": fuel,
        "passive_pass": (passive_gain >= threshold["survival_advantage_passive_vs_flat"] and
                         binding >= threshold["passive_formation_bond_gain"]),
        "driven_formation_pass": fuel >= threshold["driven_formation_bond_gain"],
        "record_persistence_pass": any(r["passes"] for r in gains),
        "record_benefits_at_budgets": gains,
    }


def dominance_certificate(primary):
    """A consequence of the saved menu, separate from the frozen budget gate."""
    receptor = [r for r in primary if r["actor"] == "record"]
    blind = [r for r in receptor if r["s"] == 0]
    acquisitive = [r for r in receptor if r["s"] > 0]
    witness = max(blind, key=lambda r: r["first_exit_survival"])
    margin = min(witness["first_exit_survival"] - r["first_exit_survival"]
                 for r in acquisitive)
    saved = {kind: min(r["ledger"][kind] - witness["ledger"][kind] for r in acquisitive)
             for kind in ("joint_reset", "separate_reset")}
    return {
        "witness": witness["name"], "compared_settings": len(acquisitive),
        "minimum_survival_margin": margin, "minimum_work_savings": saved,
        "strict_dominance_in_both_reset_models": margin > 0 and min(saved.values()) > 0,
    }


def calculate(cfg):
    length, fuel, attraction = [cfg[k] for k in ("primary_length", "primary_fuel", "J")]
    primary = [evaluate(cfg, length, fuel, attraction, mode) for mode in policy_grid(cfg)]
    flat = evaluate(cfg, length, fuel, 0.0, Policy("flat_persistence", (0, 0), 0.0))
    formation = [evaluate(cfg, length, fuel, row["J"], named(row), start="uniform")
                 for row in cfg["formation_modes"]]
    scaling = [evaluate(cfg, size, fuel, attraction, named(mode))
               for size, mode in itertools.product(cfg["dynamic_scaling_lengths"], cfg["scaling_modes"])]
    stock = [evaluate(cfg, length, k, attraction, named(mode))
             for k, mode in itertools.product(cfg["fuel_controls"], cfg["scaling_modes"])]
    sizes = [canonical_scaling(size, j) for size, j in itertools.product(
        cfg["equilibrium_lengths"], cfg["equilibrium_attractions"])]
    test = cfg["independent_grid"]
    reference = [evaluate(cfg, test["L"], test["K"], attraction, named(mode), mode["start"])
                 for mode in test["runs"]]
    numerical_rows = [*primary, flat, *formation, *scaling, *stock, *reference]
    numerical = {
        k: max(row["numerics"][k] for row in numerical_rows)
        for k in ("normalization", "generator_conservation", "detailed_balance")
    }
    numerical["minimum_probability"] = min(r["numerics"]["minimum_probability"] for r in numerical_rows)
    numerical["minimum_entropy_production"] = min(r["ledger"]["entropy_production_nats"] for r in numerical_rows)
    for key in ("fuel_budget_residual", "relative_free_energy_residual", "correlation_cost_residual"):
        numerical[key] = max(r["ledger"][key] for r in numerical_rows)
    frontiers = envelopes(cfg, primary)
    tolerance = cfg["thresholds"]["numerical_abs"]
    passive = next(r for r in stock if r["fuel"] == 0 and r["name"] == "passive")
    numerical["no_fuel_survival_residual"] = max(abs(r["first_exit_survival"] - passive["first_exit_survival"])
                                                   for r in stock if r["fuel"] == 0)
    numerical["flat_uniform_bond_residual"] = abs(formation[0]["end_bond_fraction"] - formation[0]["initial_bond_fraction"])
    return {
        "primary": primary, "flat": flat, "formation": formation,
        "dynamic_size": scaling, "fuel_controls": stock, "equilibrium_scaling": sizes,
        "frontiers": frontiers, "independent_grid_count_solver": reference,
        "science": scientific_status(cfg, primary, flat, formation, frontiers),
        "post_run_menu_dominance": dominance_certificate(primary),
        "numerical_checks": numerical,
        "numerical_pass": (all(abs(value) <= tolerance for name, value in numerical.items()
                               if not name.startswith("minimum")) and
                           numerical["minimum_probability"] >= -tolerance and
                           numerical["minimum_entropy_production"] >= -tolerance and
                           all(r["exact_count_check"] for r in sizes)),
        "lab_rate_ratio": {
            "k8": math.exp(cfg["J"] + cfg["A"]) * 8,
            "k4": math.exp(cfg["J"] + cfg["A"]) * 4 / 5,
            "log_fill_contrast": math.log(10),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--output", type=Path)
    output.add_argument("--verify", type=Path)
    args = parser.parse_args()
    cfg = json.loads((HERE / "PARAMETERS.json").read_text())
    result = calculate(cfg)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads(args.verify.read_text()):
        raise SystemExit("Stored results differ from replay")
    print(json.dumps(result["numerical_checks"], indent=2))
    print({key: value for key, value in result["science"].items() if key != "record_benefits_at_budgets"})
    if not result["numerical_pass"]:
        raise SystemExit("Registered numerical check failed")


if __name__ == "__main__":
    main()
