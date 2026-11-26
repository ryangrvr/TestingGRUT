#!/usr/bin/env python3
"""Independent full-state forward uniformization for frozen GRUT-I3."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import itertools
import json
import math
from pathlib import Path
import signal
import time

import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.special import gammainc, xlogy


HERE = Path(__file__).resolve().parent
STATES = tuple(itertools.product((0, 1), repeat=7))
INDEX = {state: i for i, state in enumerate(STATES)}
BITS = np.array(STATES)
B, M, FS, A, C, X, FR = range(7)
N = len(STATES)
LN_2 = math.log(2.0)
DEVICE_STATES, DEVICE_INVERSE = np.unique(BITS[:, 1:], axis=0, return_inverse=True)


@dataclass
class Audit:
    max_poisson_tail: float = 0.0
    max_mass_error: float = 0.0
    max_negative_probability: float = 0.0
    max_conservation_error: float = 0.0
    max_edge_log_ratio_error: float = 0.0
    max_energy_identity_error: float = 0.0
    limit_error: float = 0.0
    max_yield_mixture_error: float = 0.0
    max_prediction_formula_error: float = 0.0
    max_work_identity_error: float = 0.0
    min_entropy_production: float = math.inf
    propagations: int = 0

    def inspect(self, distribution):
        self.max_mass_error = max(self.max_mass_error, abs(float(distribution.sum()) - 1.0))
        self.max_negative_probability = max(self.max_negative_probability, -float(distribution.min()))
        self.max_conservation_error = max(
            self.max_conservation_error, float(distribution[BITS[:, X] + BITS[:, FR] != 1].sum())
        )

    def check_edges(self, rates, energies):
        for (start, end), rate in rates.items():
            reverse = rates[end, start]
            imbalance = math.log(rate / reverse) + energies[end] - energies[start]
            self.max_edge_log_ratio_error = max(self.max_edge_log_ratio_error, abs(imbalance))


def endpoint(state, changes):
    result = list(state)
    for site, value in changes:
        result[site] = value
    return INDEX[tuple(result)]


def choose_energy(config):
    return config["a_s"] * BITS[:, FS] + config["a_r"] * BITS[:, FR] - config["J"] * BITS[:, X]


def generator(cfg, stage, audit, sensor=True):
    options = {
        "acquisition": ((B, cfg["q"]), (M, cfg["mu"])),
        "storage": ((B, cfg["q"]), (M, cfg["mu"])),
        "pulses": ((B, cfg["q"]),),
        "catalysis": ((B, cfg["q"]), (A, cfg["nu"])),
    }[stage]
    rates = {}
    pi = 1 / (1 + math.exp(-cfg["a_r"] - cfg["J"]))
    for i, state in enumerate(STATES):
        for position, rate in options:
            if rate:
                j = endpoint(state, ((position, 1 - state[position]),))
                rates[i, j] = rate
        b, m, sensor_fuel, a, _, x, repair_fuel = state
        if stage == "acquisition" and sensor and cfg["s"]:
            if m != b and sensor_fuel:
                j = endpoint(state, ((M, b), (FS, 0)))
                rates[i, j] = cfg["s"] * math.exp(cfg["a_s"] / 2)
            if m == b and not sensor_fuel:
                j = endpoint(state, ((M, 1 - b), (FS, 1)))
                rates[i, j] = cfg["s"] * math.exp(-cfg["a_s"] / 2)
        if stage == "catalysis" and x + repair_fuel == 1:
            c = cfg["c_h"] if a == b else cfg["c_l"]
            rate = c * (1 - pi if x else pi)
            if rate:
                j = endpoint(state, ((X, 1 - x), (FR, 1 - repair_fuel)))
                rates[i, j] = rate
    audit.check_edges(rates, choose_energy(cfg))
    departures = np.zeros(N)
    rows, cols, vals = [], [], []
    for (start, end), rate in rates.items():
        rows.append(end)
        cols.append(start)
        vals.append(rate)
        departures[start] += rate
    rows.extend(range(N))
    cols.extend(range(N))
    vals.extend(-departures)
    return coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr(), max(departures)


def evolve(state, stage, duration, tail_limit, audit):
    if duration < 0:
        raise ValueError("Negative phase duration")
    matrix, uniform_rate = stage
    if not uniform_rate or not duration:
        audit.inspect(state)
        return state.copy()
    jump = eye(N, format="csr") + matrix / uniform_rate
    mean = duration * uniform_rate
    weight = math.exp(-mean)
    transient = state.copy()
    result = weight * transient
    k = 0
    while True:
        omitted = float(gammainc(k + 1, mean))
        if omitted <= tail_limit:
            break
        k += 1
        transient = jump @ transient
        weight *= mean / k
        result += weight * transient
    audit.propagations += 1
    audit.max_poisson_tail = max(audit.max_poisson_tail, omitted)
    audit.inspect(result)
    return result


def ready(eta):
    rho = np.zeros(N)
    for b, m, c in itertools.product((0, 1), repeat=3):
        rho[INDEX[b, m, 1, 0, c, 0, 1]] = (eta if c else 1 - eta) / 4
    return rho


def transfer(state, arm):
    if arm == "passive":
        return state
    result = np.zeros(N)
    for i, (b, m, fs, a, c, x, fr) in enumerate(STATES):
        control = 1 - c if arm == "inverted" else c
        result[INDEX[b, a, fs, m ^ control, c, x, fr]] += state[i]
    return result


def entropy(p):
    return float(-sum(xlogy(p, p)))


def h(p):
    return entropy(np.array([1 - p, p]))


def match(state, bit):
    return float(state[BITS[:, B] == BITS[:, bit]].sum())


def ledger(cfg, eta, end, audit):
    marginal = end @ BITS
    distribution = np.bincount(DEVICE_INVERSE, weights=end, minlength=len(DEVICE_STATES))
    entropy_device = entropy(distribution)
    entropy_full = entropy(end)
    entropy_device_separate = sum(h(p) for p in marginal[1:])
    initial_device = LN_2 + h(eta)
    source_device_mi = h(marginal[B]) + entropy_device - entropy_full
    preparation_work = (np.logaddexp(0.0, cfg["a_s"]) + np.logaddexp(0.0, cfg["a_r"])
                        + np.logaddexp(0.0, cfg["J"]) + 3 * LN_2 - initial_device)
    sensor_use, repair_use = 1 - marginal[FS], 1 - marginal[FR]
    heat = float(cfg["a_s"] + cfg["a_r"] - end @ choose_energy(cfg))
    expected_heat = cfg["a_s"] * sensor_use + cfg["a_r"] * repair_use + cfg["J"] * marginal[X]
    sigma = heat + entropy_full - initial_device - LN_2
    reset_joint = heat + entropy_device - initial_device
    reset_separate = heat + entropy_device_separate - initial_device
    audit.min_entropy_production = min(audit.min_entropy_production, sigma)
    audit.max_energy_identity_error = max(audit.max_energy_identity_error,
                                           abs(heat - expected_heat))
    audit.max_work_identity_error = max(audit.max_work_identity_error,
                                       abs(reset_joint - sigma - source_device_mi),
                                       abs(h(marginal[B]) - LN_2))
    return {
        "bond_probability": float(marginal[X]),
        "sensor_fuel_consumed": float(sensor_use),
        "repair_fuel_consumed": float(repair_use),
        "heat": heat,
        "entropy_full": entropy_full,
        "entropy_device": entropy_device,
        "entropy_device_separate": entropy_device_separate,
        "entropy_production_full": sigma,
        "source_device_mi": source_device_mi,
        "total_device_correlation": entropy_device_separate - entropy_device,
        "preparation_work": float(preparation_work),
        "joint_reset_work": reset_joint,
        "separate_reset_work": reset_separate,
    }


def row(name, cfg, arms, tail_limit, audit):
    stages = {tag: generator(cfg, tag, audit) for tag in
              ("acquisition", "storage", "pulses", "catalysis")}
    disabled = generator(cfg, "acquisition", audit, sensor=False)
    yields = {}
    for label, a in (("match", 0), ("mismatch", 1)):
        state = np.zeros(N)
        state[INDEX[0, 0, 1, a, 0, 0, 1]] = 1.0
        end = evolve(state, stages["catalysis"], cfg["T_r"], tail_limit, audit)
        yields[label] = float(end @ BITS[:, X])
    results = {}
    for arm in arms:
        eta = 0.5 if arm == "shuffled" else cfg["eta"]
        initial = ready(eta)
        acquired = evolve(initial, disabled if arm == "disabled" else stages["acquisition"],
                          cfg["t_s"], tail_limit, audit)
        stored = evolve(acquired, stages["storage"], cfg["tau"] - cfg["delta"], tail_limit, audit)
        pulse = evolve(stored, stages["pulses"], cfg["delta"], tail_limit, audit)
        actuator = transfer(pulse, arm)
        end = evolve(actuator, stages["catalysis"], cfg["T_r"], tail_limit, audit)
        result = {
            "acquisition_match": match(acquired, M),
            "stored_match": match(stored, M),
            "actuator_match": match(actuator, A),
            **ledger(cfg, eta, end, audit),
        }
        score = yields["mismatch"] + result["actuator_match"] * (yields["match"] - yields["mismatch"])
        audit.max_yield_mixture_error = max(audit.max_yield_mixture_error,
                                            abs(score - result["bond_probability"]))
        decay = math.exp(-2 * (cfg["q"] * cfg["tau"] + cfg["mu"] * (cfg["tau"] - cfg["delta"])))
        sign = -1 if arm == "inverted" else 1
        predicted = 0.5 if arm == "passive" else (
            0.5 + sign * (result["acquisition_match"] - 0.5) * decay * (1 - 2 * eta)
        )
        audit.max_prediction_formula_error = max(audit.max_prediction_formula_error,
                                                  abs(predicted - result["actuator_match"]))
        results[arm] = result
    return {"name": name, "parameters": cfg, "yields": yields, "arms": results}


def independent_limit(cfg, tail_limit, audit):
    fixed = dict(cfg, q=0.0, mu=0.0, nu=0.0, eta=0.0, tau=cfg["delta"])
    check = row("analytic_test_only", fixed, ("predictive",), tail_limit, audit)
    pos_rate = fixed["s"] * math.exp(fixed["a_s"] / 2)
    total_rate = pos_rate + fixed["s"] * math.exp(-fixed["a_s"] / 2)
    expected_sensor = 0.5 + 0.5 * pos_rate / total_rate * (-math.expm1(-total_rate * fixed["t_s"]))
    limit_pi = 1 / (1 + math.exp(-fixed["J"] - fixed["a_r"]))
    arm = check["arms"]["predictive"]
    audit.limit_error = max(
        abs(expected_sensor - arm["acquisition_match"]),
        abs(expected_sensor - arm["actuator_match"]),
        abs(limit_pi * (-math.expm1(-fixed["c_h"] * fixed["T_r"])) - check["yields"]["match"]),
        abs(limit_pi * (-math.expm1(-fixed["c_l"] * fixed["T_r"])) - check["yields"]["mismatch"]),
    )


def run(spec):
    audit = Audit()
    primary, sweep, arms, gates = (spec[k] for k in ("primary", "sweep", "arms", "gates"))
    tail_limit = gates["uniformization_tail"]
    rows = []
    for (i, q), (j, tau) in itertools.product(enumerate(sweep["q"]), enumerate(sweep["tau"])):
        rows.append(row(f"grid_{i:02d}_{j:02d}", dict(primary, q=q, tau=tau),
                        arms, tail_limit, audit))
    for name, changes in spec["stress"].items():
        rows.append(row(name, primary | changes, arms, tail_limit, audit))
    independent_limit(primary, tail_limit, audit)
    checks = vars(audit)
    tol = gates["probability_tol"]
    errors = [value for name, value in checks.items() if (
        name.endswith("_error") or name == "max_negative_probability"
    )]
    if (max(errors) > tol or checks["max_poisson_tail"] > tail_limit
            or checks["min_entropy_production"] < -tol
            or len(rows) != spec["production"]["parameter_rows"]):
        raise AssertionError(checks)
    return {"rows": rows, "method": "labelled-state sparse forward uniformization", "audit": checks}


def timeout(*_):
    raise TimeoutError("Recorded single-solver wall-clock budget exceeded")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    spec = json.loads((HERE / "PARAMETERS.json").read_text())
    limit = spec["gates"]["wall_seconds_per_solver"]
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(limit)
    start = time.monotonic()
    result = run(spec)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    signal.alarm(0)
    print(f"{len(result['rows'])} labelled-state rows in {time.monotonic() - start:.3f} s")
    print(json.dumps(result["audit"], indent=2))


if __name__ == "__main__":
    main()
