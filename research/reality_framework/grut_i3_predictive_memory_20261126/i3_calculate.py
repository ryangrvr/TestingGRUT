"""Finite-state replay of the frozen screened assembly protocol."""
from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.special import expit, xlogy
from scipy.stats import binom


@dataclass(frozen=True)
class Rates:
    q: float
    mu: float
    nu: float
    eta: float
    t_s: float
    tau: float
    delta: float
    T_r: float
    a_s: float
    a_r: float
    J: float
    s: float
    c_h: float
    c_l: float

    @property
    def pi(self):
        return float(expit(self.a_r + self.J))


BITS = (0, 1)
TRIPLES = tuple(product(BITS, repeat=3))
SHAPE3 = (2, 2, 2)


def add_edge(generator, origin, destination, rate):
    source = np.ravel_multi_index(origin, SHAPE3)
    target = np.ravel_multi_index(destination, SHAPE3)
    generator[target, source] += rate
    generator[source, source] -= rate


def sensor_generator(p, enabled):
    generator = np.zeros((8, 8))
    for b, m, fuel in TRIPLES:
        add_edge(generator, (b, m, fuel), (1 - b, m, fuel), p.q)
        add_edge(generator, (b, m, fuel), (b, 1 - m, fuel), p.mu)
    if enabled:
        for b in BITS:
            depleted = (b, b, 0)
            charged = (b, 1 - b, 1)
            add_edge(generator, charged, depleted, p.s * math.exp(p.a_s / 2))
            add_edge(generator, depleted, charged, p.s * math.exp(-p.a_s / 2))
    return generator


def response_generator(p):
    generator = np.zeros((8, 8))
    for b, a, x in TRIPLES:
        add_edge(generator, (b, a, x), (1 - b, a, x), p.q)
        add_edge(generator, (b, a, x), (b, 1 - a, x), p.nu)
        rate_scale = p.c_h if b == a else p.c_l
        directional_rate = p.pi if x == 0 else 1 - p.pi
        add_edge(generator, (b, a, x), (b, a, 1 - x), rate_scale * directional_rate)
    return generator


def flip_kernel(rate, time):
    off_diagonal = -math.expm1(-2 * rate * time) / 2
    diagonal = 1 - off_diagonal
    return np.array([[diagonal, off_diagonal], [off_diagonal, diagonal]])


def linked_mass(distribution):
    return float(distribution[0, 0].sum() + distribution[1, 1].sum())


def entropy(distribution):
    return float(-xlogy(distribution, distribution).sum())


def marginal(distribution, axis):
    return distribution.sum(axis=tuple(i for i in range(distribution.ndim) if i != axis))


def check_balance(generator, energy, tol):
    np.testing.assert_allclose(generator.sum(axis=0), 0, atol=tol, rtol=0)
    gibbs_weight = np.exp(-np.array(energy))
    flux = generator * gibbs_weight[np.newaxis, :]
    np.testing.assert_allclose(flux, flux.T, atol=tol, rtol=0)


def check_distribution(distribution, tol):
    if float(distribution.min()) < -tol:
        raise AssertionError("Negative probability")
    if abs(float(distribution.sum()) - 1) > tol:
        raise AssertionError("Probability mass lost")


def acquire(p, armed, tol):
    generator = sensor_generator(p, armed)
    check_balance(generator, [p.a_s * fuel for _, _, fuel in TRIPLES], tol)
    initial = np.zeros(SHAPE3)
    initial[:, :, 1] = 0.25
    sensed = (expm(generator * p.t_s) @ initial.ravel()).reshape(SHAPE3)
    free_delay = p.tau - p.delta
    stored = np.einsum(
        "bi,mj,ijf->bmf", flip_kernel(p.q, free_delay), flip_kernel(p.mu, free_delay), sensed
    )
    arrivals = np.einsum("bi,imf->bmf", flip_kernel(p.q, p.delta), stored)
    for distribution in (sensed, stored, arrivals):
        check_distribution(distribution, tol)
    return sensed, stored, arrivals


def route(arrivals, arm, error):
    start = np.zeros((2,) * 6)  # B, M, F_s, A, C, X.
    for b, m, fuel in TRIPLES:
        for coin, mass in enumerate((1 - error, error)):
            if arm == "passive":
                stored, orientation = m, 0
            else:
                stored = 0
                orientation = m ^ coin ^ int(arm == "inverted")
            start[b, stored, fuel, orientation, coin, 0] += mass * arrivals[b, m, fuel]
    return start


def yields_and_propagator(p, tol):
    relative_rate = p.q + p.nu
    relative = np.array([[-relative_rate, relative_rate], [relative_rate, -relative_rate]])
    survival = expm((relative - np.diag([p.c_h, p.c_l])) * p.T_r).sum(axis=0)
    match, mismatch = p.pi * (1 - survival)

    generator = response_generator(p)
    check_balance(generator, [p.a_r * (1 - x) - p.J * x for _, _, x in TRIPLES], tol)
    propagation = expm(generator * p.T_r).reshape((2,) * 6)
    for initial_b, initial_a in product(BITS, repeat=2):
        end = propagation[..., initial_b, initial_a, 0]
        check_distribution(end, tol)
        target = match if initial_b == initial_a else mismatch
        np.testing.assert_allclose(end[:, :, 1].sum(), target, atol=tol, rtol=0)
    return float(match), float(mismatch), propagation


def work_and_entropy(final, p, error):
    # Full ordering is B, M, F_s, A, C, X, F_r. No reset access to B.
    joint = final.sum(axis=0)
    entropy_full = entropy(final)
    entropy_device = entropy(joint)
    entropy_separate = sum(entropy(marginal(joint, i)) for i in range(joint.ndim))
    entropy_start_device = math.log(2) + entropy(np.array([1 - error, error]))
    sensing_use = float(final[:, :, 0].sum())
    assembly_use = float(final[..., 0].sum())
    heat = p.a_s * sensing_use + (p.a_r + p.J) * assembly_use
    preparation = (
        np.logaddexp(0, p.a_s)
        + np.logaddexp(0, p.a_r)
        + np.logaddexp(0, p.J)
        + 2 * math.log(2)
        - entropy(np.array([1 - error, error]))
    )
    joint_reset = heat + entropy_device - entropy_start_device
    separate_reset = heat + entropy_separate - entropy_start_device
    return {
        "sensor_fuel_consumed": sensing_use,
        "repair_fuel_consumed": assembly_use,
        "heat": heat,
        "entropy_full": entropy_full,
        "entropy_device": entropy_device,
        "entropy_device_separate": entropy_separate,
        "entropy_production_full": heat + entropy_full - entropy_start_device - math.log(2),
        "source_device_mi": math.log(2) + entropy_device - entropy_full,
        "total_device_correlation": entropy_separate - entropy_device,
        "preparation_work": float(preparation),
        "joint_reset_work": joint_reset,
        "separate_reset_work": separate_reset,
    }


def reset_budget(p, gates):
    preparation_upper = float(
        np.logaddexp(0, p.a_s)
        + np.logaddexp(0, p.a_r)
        + np.logaddexp(0, p.J)
        + 2 * math.log(2)
    )
    total = preparation_upper + gates["common_work_allowance"]
    return {
        "preparation_upper": preparation_upper,
        "common_allowance": gates["common_work_allowance"],
        "total_upper": total,
        "cap_slack": gates["work_cap"] - total,
    }


def check_arm(final, stats, p, arm, error, target, tol):
    check_distribution(final, tol)
    for x in BITS:
        if np.any(final[..., x, x] != 0):
            raise AssertionError("Material fuel invariant violated")
    np.testing.assert_allclose(marginal(final, 0), [0.5, 0.5], atol=tol, rtol=0)
    np.testing.assert_allclose(
        stats["bond_probability"], stats["repair_fuel_consumed"], atol=tol, rtol=0
    )
    np.testing.assert_allclose(stats["bond_probability"], target, atol=tol, rtol=0)
    if stats["entropy_production_full"] < -tol or stats["source_device_mi"] < -tol:
        raise AssertionError("Thermodynamic or information inequality violated")
    if stats["joint_reset_work"] > stats["separate_reset_work"] + tol:
        raise AssertionError("Separate reset failed to charge discarded correlations")
    if stats["separate_reset_work"] > stats["preparation_work"] + tol:
        raise AssertionError("Nonnegative endpoint free energy was lost")
    eps = 1 - stats["acquisition_match"]
    switch_contrast = math.exp(-2 * p.q * p.tau - 2 * p.mu * (p.tau - p.delta))
    if arm == "passive":
        expected = 0.5
    else:
        sign = -1 if arm == "inverted" else 1
        expected = 0.5 + 0.5 * sign * (1 - 2 * eps) * switch_contrast * (1 - 2 * error)
    np.testing.assert_allclose(stats["actuator_match"], expected, atol=tol, rtol=0)
    if arm == "disabled":
        np.testing.assert_allclose(stats["acquisition_match"], 0.5, atol=tol, rtol=0)
        np.testing.assert_allclose(stats["sensor_fuel_consumed"], 0, atol=tol, rtol=0)


def run_row(name, p, arms, gates):
    tol = gates["probability_tol"]
    match, mismatch, propagation = yields_and_propagator(p, tol)
    acquisition = {active: acquire(p, active, tol) for active in (True, False)}
    metrics = {}
    for arm in arms:
        sensed, stored, arrivals = acquisition[arm != "disabled"]
        error = 0.5 if arm == "shuffled" else p.eta
        start = route(arrivals, arm, error)
        # Contract the active response indices B, A and X. Other carriers are disconnected.
        end = np.einsum("ijkbax,bmfacx->imfjck", propagation, start)
        final = np.zeros((2,) * 7)
        final[..., 0, 1] = end[..., 0]
        final[..., 1, 0] = end[..., 1]
        actuator = start.sum(axis=(1, 2, 4, 5))
        match_at_start = linked_mass(actuator)
        outcome = float(end[..., 1].sum())
        information = entropy(marginal(actuator, 0)) + entropy(marginal(actuator, 1)) - entropy(actuator)
        stats = {
            "acquisition_match": linked_mass(sensed),
            "stored_match": linked_mass(stored),
            "actuator_match": match_at_start,
            "action_source_mi": max(0.0, information),
            "bond_probability": outcome,
        }
        stats.update(work_and_entropy(final, p, error))
        check_arm(final, stats, p, arm, error, mismatch + (match - mismatch) * match_at_start, tol)
        metrics[arm] = stats

    c_bar = (p.c_h + p.c_l) / 2
    feedback_bound = -math.expm1(-p.pi * c_bar * p.T_r)
    nulls = {
        "fixed_orientation": (match + mismatch) / 2,
        "relaxed_open_loop": -p.pi * math.expm1(-c_bar * p.T_r),
        "relaxed_material_feedback": feedback_bound,
        "feedback_with_cue_bound": min(1.0, feedback_bound + gates["leakage_tv"]),
        "perfect_current_cue": match,
    }
    budget = reset_budget(p, gates)
    realized = metrics["predictive"]["bond_probability"]
    required_match = (gates["reliability"] - mismatch) / (match - mismatch)
    min_action_info = None
    if 0.5 <= required_match <= 1.0:
        min_action_info = math.log(2) - entropy(np.array([required_match, 1 - required_match]))
    decisions = {
        "record_reaches_target": realized >= gates["reliability"],
        "null_strictly_below_target_with_leakage": nulls["feedback_with_cue_bound"] < gates["reliability"],
        "work_passes_cap": budget["preparation_upper"] <= gates["worst_preparation_limit"]
        and budget["total_upper"] <= gates["work_cap"],
        "advantage_over_feedback_with_leakage": realized - nulls["feedback_with_cue_bound"],
        "required_locked_match": required_match,
        "required_action_information_nats": min_action_info,
    }
    return {
        "name": name,
        "parameters": vars(p),
        "arms": metrics,
        "yields": {"match": match, "mismatch": mismatch},
        "nulls": nulls,
        "resource_budget": budget,
        "decisions": decisions,
    }


def analytical_floors(p):
    forward = p.s * math.exp(p.a_s / 2)
    reverse = p.s * math.exp(-p.a_s / 2)
    sensor = math.exp(-(p.q + p.mu) * p.t_s) * (
        0.5 - 0.5 * forward / (forward + reverse) * math.expm1(-(forward + reverse) * p.t_s)
    )
    accuracy = 0.5 + (sensor - 0.5) * math.exp(
        -2 * p.q * p.tau - 2 * p.mu * (p.tau - p.delta)
    ) * (1 - 2 * p.eta)
    matched = -math.exp(-(p.q + p.nu) * p.T_r) * p.pi * math.expm1(-p.c_h * p.T_r)
    unmatched = -p.pi * math.expm1(-p.c_l * p.T_r)
    return {
        "acquisition": sensor,
        "decision": accuracy,
        "bond": unmatched + (matched - unmatched) * accuracy,
    }


def scale(primary, sites, reliability):
    predictive = primary["arms"]["predictive"]["bond_probability"]
    null = primary["nulls"]["relaxed_material_feedback"]
    leaky = primary["nulls"]["feedback_with_cue_bound"]
    results = []
    for count in sites:
        needed = math.ceil(reliability * count)
        results.append({
            "sites": count,
            "threshold_bonds": needed,
            "predictive_at_least_threshold": float(binom.sf(needed - 1, count, predictive)),
            "null_upper_at_least_threshold": float(binom.sf(needed - 1, count, null)),
            "leaky_null_upper_at_least_threshold": float(binom.sf(needed - 1, count, leaky)),
            "preparation_and_control_upper": count * primary["resource_budget"]["total_upper"],
            "fuel_capacity": 2 * count,
        })
    return results


def run(definition):
    main = Rates(**definition["primary"])
    cases = {}
    for i, q in enumerate(definition["sweep"]["q"]):
        for j, tau in enumerate(definition["sweep"]["tau"]):
            cases[f"grid_{i:02d}_{j:02d}"] = replace(main, q=q, tau=tau)
    for label, overrides in definition["stress"].items():
        cases[label] = replace(main, **overrides)
    if len(cases) != definition["production"]["parameter_rows"]:
        raise AssertionError("Registered row count changed")
    gates = definition["gates"]
    rows = [run_row(name, p, definition["arms"], gates) for name, p in cases.items()]
    primary = next(row for row in rows if row["parameters"] == vars(main))
    floors = analytical_floors(main)
    if floors["bond"] <= 0.71 or primary["nulls"]["feedback_with_cue_bound"] >= 0.67:
        raise AssertionError("Analytical feasibility certificate was not recovered")
    if floors["bond"] > primary["arms"]["predictive"]["bond_probability"]:
        raise AssertionError("Numerical result contradicts the analytical lower bound")
    # Independent sourced-telegraph and unreversed reaction limits.
    for q, tau in ((0, 1), (0.04, 1), (1, 4)):
        contrast = 0.5 + (0.5 - 0.13) * math.exp(-2 * q * tau)
        from_kernel = np.array([0.87, 0.13]) @ flip_kernel(q, tau)[:, 0]
        np.testing.assert_allclose(from_kernel, contrast, atol=gates["probability_tol"], rtol=0)
    zero_switch = replace(main, q=0, nu=0)
    h, low, _ = yields_and_propagator(zero_switch, gates["probability_tol"])
    np.testing.assert_allclose(
        [h, low], [-main.pi * math.expm1(-c * main.T_r) for c in (main.c_h, main.c_l)],
        atol=gates["probability_tol"], rtol=0,
    )
    return {
        "protocol": definition["protocol"],
        "preregistration": "b78f6199ec8a240e75e808cfd7204c915f5a6e54",
        "primary_row": primary["name"],
        "analytical_primary_floors": floors,
        "rows": rows,
        "scaling": scale(primary, definition["scaling_sites"], gates["reliability"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    definition = json.loads(Path(__file__).with_name("PARAMETERS.json").read_text())
    results = run(definition)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
