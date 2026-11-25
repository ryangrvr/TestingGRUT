#!/usr/bin/env python3
"""Replay the frozen GRUT-I1 conventional comparisons."""

import argparse
from dataclasses import dataclass
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import binom


HERE = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Rates:
    source_rate: float
    memory_rate: float
    hot_rate: float
    cold_rate: float
    cold_fraction: float


def entropy(x):
    return -sum(p * math.log(p) for p in (x, 1 - x) if p > 0)


def heat_record(x, initial, hot, cold, cfg):
    beta_h = math.log((1 - cfg["hot_excited_fraction"]) / cfg["hot_excited_fraction"])
    beta_c = math.log((1 - cfg["cold_excited_fraction"]) / cfg["cold_excited_fraction"])
    return {
        "excited_probability": x,
        "hot_heat_into_cells": hot,
        "cold_heat_into_cells": cold,
        "energy_residual": x - initial + hot + cold,
        "marginal_entropy_production_nats": entropy(x) - entropy(initial) + beta_h * hot + beta_c * cold,
    }


def exact_thermal(n, rates, initial, cfg):
    dt = cfg["episode_time"] / n
    qh, qc = cfg["hot_excited_fraction"], cfg["cold_excited_fraction"]
    ph, pc = rates.hot_rate * dt, rates.cold_rate * dt
    cold_count = round(n * rates.cold_fraction)
    hot_heat, cold_heat, x = 0.0, 0.0, initial
    for j in range(n):
        leaving = ph * (x - qh)
        hot_heat += leaving
        x -= leaving
        if j < cold_count:
            leaving = pc * (x - qc)
            cold_heat += leaving
            x -= leaving
    attenuation = (1 - ph) ** n * (1 - pc) ** cold_count
    return heat_record(x, initial, hot_heat, cold_heat, cfg) | {
        "attenuation": attenuation,
        "hot_cells_used": n if rates.hot_rate else 0,
        "cold_cells_used": cold_count if rates.cold_rate else 0,
        "flip_cells_used": n * (bool(rates.source_rate) + bool(rates.memory_rate)),
    }


def continuum_thermal(rates, initial, cfg):
    time = cfg["episode_time"]
    qh, qc = cfg["hot_excited_fraction"], cfg["cold_excited_fraction"]
    cold_time = time * rates.cold_fraction
    hot_heat, cold_heat, x = 0.0, 0.0, initial
    segments = ((cold_time, rates.cold_rate), (time - cold_time, 0.0))
    for duration, cold in segments:
        rate = rates.hot_rate + cold
        if rate == 0 or duration == 0:
            continue
        fixed = (rates.hot_rate * qh + cold * qc) / rate
        integral = fixed * duration + (x - fixed) * (-math.expm1(-rate * duration)) / rate
        hot_heat += rates.hot_rate * (integral - qh * duration)
        cold_heat += cold * (integral - qc * duration)
        x = fixed + (x - fixed) * math.exp(-rate * duration)
    return heat_record(x, initial, hot_heat, cold_heat, cfg) | {
        "attenuation": math.exp(-rates.hot_rate * time - rates.cold_rate * cold_time),
    }


def phase_factor(n, coupling, switch, time):
    dt = time / n
    p = switch * dt
    stochastic = np.array([[1 - p, p], [p, 1 - p]])
    phase = np.diag(np.exp(np.array([-2j, 2j]) * coupling * dt))
    return complex(np.linalg.matrix_power(stochastic @ phase, n)[:, 0].sum())


def phase_limit(coupling, switch, time):
    generator = np.array([[-switch - 2j * coupling, switch],
                          [switch, -switch + 2j * coupling]])
    return complex(expm(time * generator)[:, 0].sum())


def recovery(z, contrast, attenuation, phase_errors):
    coherence = math.sqrt(attenuation)
    correctable = math.hypot(z.real, contrast * z.imag)
    blind_coherence = abs(z.real)
    base = (1 + attenuation) / 4
    blind = base + coherence * blind_coherence / 2
    ideal = base + coherence * correctable / 2
    return {
        "fidelity_assisted": ideal,
        "fidelity_blind": blind,
        "responsiveness": ideal - blind,
        "phase_angle_magnitude": math.atan2(abs(contrast * z.imag), z.real) / 2,
        "jitter": [{
            "epsilon": epsilon,
            "assisted": max(blind, base + coherence * math.cos(2 * epsilon) * correctable / 2),
            "blind": blind,
            "responsiveness": max(0.0, coherence * (math.cos(2 * epsilon) * correctable - blind_coherence) / 2),
        } for epsilon in phase_errors],
    }


def channel_error_bound(n, coupling, rates, cfg):
    """CPTP contraction and first-order splitting bound, before cap at 2."""
    dt = cfg["episode_time"] / n
    accumulated = 0.0
    for cold, count in ((rates.cold_rate, round(n * rates.cold_fraction)),
                        (0.0, n - round(n * rates.cold_fraction))):
        flip = [rates.source_rate * dt, rates.memory_rate * dt]
        exchange = [rates.hot_rate * dt, cold * dt]
        exact_vs_linear = (
            sum(2 * p * p / (1 - 2 * p) for p in flip)
            + sum(p * p / (1 - p) for p in exchange)
        )
        norms = [2 * abs(coupling), 2 * rates.source_rate, 2 * rates.memory_rate,
                 2 * rates.hot_rate, 2 * cold]
        splitting = dt**2 * sum(a * b for a, b in itertools.combinations(norms, 2))
        accumulated += count * (exact_vs_linear + splitting)
    return min(2.0, accumulated)


def bridge(n, name, assay, d0, phase, coupling, rates, cfg):
    t = cfg["episode_time"]
    exact = exact_thermal(n, rates, cfg["payload_population_initial"]["bell"], cfg)
    limit = continuum_thermal(rates, cfg["payload_population_initial"]["bell"], cfg)
    z = phase_factor(n, coupling, rates.source_rate, t)
    z_limit = phase_limit(coupling, rates.source_rate, t)
    d = d0 * (1 - 2 * rates.memory_rate * t / n) ** n
    d_limit = d0 * math.exp(-2 * rates.memory_rate * t)
    measured = recovery(z, d, exact["attenuation"], cfg["phase_error_radians"])
    continuous = recovery(z_limit, d_limit, limit["attenuation"], cfg["phase_error_radians"])
    discrepancies = {
        "phase": abs(z - z_limit), "contrast": abs(d - d_limit),
        "population": abs(exact["excited_probability"] - limit["excited_probability"]),
        "assisted": abs(measured["fidelity_assisted"] - continuous["fidelity_assisted"]),
        "blind": abs(measured["fidelity_blind"] - continuous["fidelity_blind"]),
        "jitter": max(abs(e["assisted"] - l["assisted"])
                      for e, l in zip(measured["jitter"], continuous["jitter"], strict=True)),
    }
    return {
        "slices": n, "regime": name, "assay": assay, "phase": phase,
        "thermal": exact, "phase_factor": [z.real, z.imag],
        "memory_initial_contrast": d,
        "memory_current_sign_contrast": d * (1 - 2 * rates.source_rate * t / n) ** n,
        **measured,
        "continuum": {"thermal": limit, "phase_factor": [z_limit.real, z_limit.imag],
                      "memory_initial_contrast": d_limit, **continuous},
        "continuum_abs_errors": discrepancies,
        "channel_diamond_bound": channel_error_bound(n, coupling, rates, cfg),
    }


def reuse_thermal(n, rates, initial, cfg):
    # Closed three-qubit episode. The same hot and cold cells are reused.
    dt = cfg["episode_time"] / n
    qh, qc = cfg["hot_excited_fraction"], cfg["cold_excited_fraction"]
    basis = list(itertools.product((0, 1), repeat=3))
    index = {s: i for i, s in enumerate(basis)}

    def exchange(bath_index, rate):
        matrix = np.zeros((8, 8))
        for state in basis:
            if state[0] != state[bath_index]:
                swapped = list(state)
                swapped[0], swapped[bath_index] = swapped[bath_index], swapped[0]
                matrix[index[tuple(swapped)], index[state]] = 1
        return expm(-1j * math.asin(math.sqrt(rate * dt)) * matrix)

    step = exchange(2, rates.cold_rate) @ exchange(1, rates.hot_rate)
    energies = np.array(basis, dtype=float)
    diagonal = np.prod(np.where(energies, [initial, qh, qc], [1 - initial, 1 - qh, 1 - qc]), axis=1)
    unitary = np.linalg.matrix_power(step, n)
    final = abs(unitary) ** 2 @ diagonal
    x, h, c = map(float, final @ energies)
    fresh = exact_thermal(n, rates, initial, cfg)["excited_probability"]
    return {
        "slices": n, "initial_payload_excitation": initial, "hot_cells": 1, "cold_cells": 1,
        "remaining_payload_excitation": x, "fresh_payload_excitation": fresh,
        "reuse_fresh_abs_gap": abs(x - fresh), "hot_heat_into_cells": h - qh,
        "cold_heat_into_cells": c - qc, "energy_residual": x + h + c - initial - qh - qc,
        "unitarity_residual": float(np.max(abs(unitary.conj().T @ unitary - np.eye(8)))),
    }


def organization(cfg, rates):
    rows = []
    for name in cfg["organization_regimes"]:
        n = cfg["organization_slices"]
        thermal = exact_thermal(n, rates[name], cfg["payload_population_initial"]["organization"], cfg)
        x = thermal["excited_probability"]
        for count in cfg["organization_sizes"]:
            max_damage = math.floor(count * cfg["organization_damage_fraction_threshold"])
            rows.append({
                "regime": name, "sites": count, "mean_damage_fraction": x,
                "zero_damage_probability": (1 - x) ** count,
                "coarse_success_probability": float(binom.cdf(max_damage, count, x)),
                "cold_cells_per_site": thermal["cold_cells_used"],
                "cold_heat_per_site": thermal["cold_heat_into_cells"],
                "entropy_production_per_site_nats": thermal["marginal_entropy_production_nats"],
                "energy_residual_per_site": thermal["energy_residual"],
            })
    return rows


def collective(cfg):
    results = []
    for n in cfg["slices"]:
        p = cfg["collective_individual_rate"] * cfg["episode_time"] / n
        independent = (1 - p) ** n
        for match, bright_p in (("mean", 2 * p), ("peak", p)):
            bright = (1 - bright_p) ** n
            results.append({
                "slices": n, "match": match, "cold_cells_each_arm": 2 * n,
                "independent_survival": independent, "bright_survival": bright,
                "dark_survival": 1.0, "equal_site_defect_survival": (1 + bright) / 2,
                "independent_single_link_angle": math.asin(math.sqrt(p)),
                "collective_single_link_angle": math.asin(math.sqrt(bright_p)) / math.sqrt(2),
            })
    return results


def hopping(n, architecture, incident):
    graph = np.zeros((n, n))
    if architecture == "complete":
        graph.fill(incident / (n - 1))
        np.fill_diagonal(graph, 0)
    if architecture in ("ring", "ring_antipodal"):
        degree = 2 if architecture == "ring" else 3
        for i in range(n):
            graph[i, (i + 1) % n] = graph[i, (i - 1) % n] = incident / degree
            if degree == 3:
                graph[i, (i + n // 2) % n] = incident / degree
    if architecture == "dyadic":
        levels = n.bit_length() - 1
        norm = 1 - 2 ** (-levels)
        for i in range(n):
            for j in range(i + 1, n):
                level = (i ^ j).bit_length()
                graph[i, j] = graph[j, i] = (incident / norm) * 2 ** (1 - 2 * level)
    return graph


def spectral_average(values, amplitudes, tolerance):
    probability = np.zeros(amplitudes.shape[0])
    start = 0
    for j, eigenvalue in enumerate(values[1:], start=1):
        if eigenvalue - values[start] > tolerance:
            probability += abs(amplitudes[:, start:j].sum(axis=1)) ** 2
            start = j
    return probability + abs(amplitudes[:, start:].sum(axis=1)) ** 2


def geometry(cfg):
    options = cfg["networks"]
    j0, tolerance = options["incident_coupling"], options["eigenvalue_group_tolerance"]
    rows, residuals = [], []
    for n, kind in itertools.product(options["sizes"], options["classes"]):
        matrix = hopping(n, kind, j0)
        values, vectors = np.linalg.eigh(matrix)
        amplitudes = vectors * vectors[0]
        average = spectral_average(values, amplitudes, tolerance)
        residuals += [float(max(abs(matrix.sum(axis=1) - j0))), abs(float(average.sum()) - 1)]
        if kind == "ring":
            assert float(average[:2].sum()) <= 4 / n + cfg["thresholds"]["algebra_abs"]
        times = []
        for time in options["times"]:
            probability = abs(amplitudes @ np.exp(-1j * time * values)) ** 2
            residuals.append(abs(float(probability.sum()) - 1))
            if kind == "complete":
                analytic = 4 / n ** 2 * math.sin(j0 * n * time / (2 * (n - 1))) ** 2
                residuals += [abs(float(probability[n // 2]) - analytic),
                              abs(float(probability[0]) - (1 - (n - 1) * analytic))]
            times.append({
                "time": time, "local_trace_distance": float(probability[:2].sum()),
                "local_guessing_advantage": float(probability[:2].sum()) / 2,
                "opposite_receiver": float(probability[n // 2]),
            })
        rows.append({
            "architecture": kind, "sites": n, "edges": int(np.count_nonzero(matrix)) // 2,
            "incident_coupling": j0, "injected_mean_energy": options["bare_gap"],
            "injected_energy_variance": float(np.dot(matrix[0], matrix[0])),
            "local_cesaro_trace_distance": float(average[:2].sum()),
            "opposite_cesaro_probability": float(average[n // 2]),
            "samples": times,
        })
    return rows, max(residuals)


def calculate(cfg):
    rates = {name: Rates(**setting) for name, setting in cfg["regimes"].items()}
    bridge_rows = [bridge(n, regime, assay, d0, phase, math.pi * severity, rates[regime], cfg)
                   for n, regime, (assay, d0), (phase, severity) in itertools.product(
                       cfg["slices"], rates, cfg["assays"].items(), cfg["phase_over_pi"].items())]
    reuse = [reuse_thermal(n, rates[cfg["reuse_regime"]], x0, cfg)
             for n, x0 in itertools.product(cfg["slices"], cfg["payload_population_initial"].values())]
    counts = organization(cfg, rates)
    collective_rows = collective(cfg)
    graph, graph_error = geometry(cfg)
    gate = cfg["thresholds"]
    anchor = {"shallow": (1 + math.sqrt(5 / 8)) / 2, "deep": 0.75}
    anchor_errors = []
    assert all(0 <= r["thermal"]["excited_probability"] <= 1 for r in bridge_rows)
    for r in bridge_rows:
        if r["regime"] == "anchor":
            severity = r["phase"]
            target = anchor[severity]
            blind = (1 + 1 / math.sqrt(2)) / 2 if severity == "shallow" else 0.5
            optimum = {"none": blind, "weak": target, "strong": 1}[r["assay"]]
            anchor_errors.extend((abs(r["fidelity_assisted"] - optimum), abs(r["fidelity_blind"] - blind)))
    fine = max(error for r in bridge_rows if r["slices"] == max(cfg["slices"])
               for error in r["continuum_abs_errors"].values())
    coarse = max(error for r in bridge_rows if r["slices"] == min(cfg["slices"])
                 for error in r["continuum_abs_errors"].values())
    energy = max([abs(r["thermal"]["energy_residual"]) for r in bridge_rows]
                 + [abs(r["energy_residual"]) for r in reuse]
                 + [abs(r["energy_residual_per_site"]) for r in counts])
    entropy_min = min([r["thermal"]["marginal_entropy_production_nats"] for r in bridge_rows]
                      + [r["entropy_production_per_site_nats"] for r in counts])
    checks = {
        "d7_anchor": max(anchor_errors) <= gate["algebra_abs"],
        "geometry_analytic_controls": graph_error <= gate["algebra_abs"],
        "energy": energy <= gate["energy_abs"],
        "marginal_entropy": entropy_min >= gate["entropy_floor"],
        "fine_continuum": fine <= gate["fine_continuum_abs"] and fine < coarse,
        "reuse_failure_preserved": max(r["reuse_fresh_abs_gap"] for r in reuse) >=
        gate["reuse_min_max_population_gap"],
        "reuse_unitarity": max(r["unitarity_residual"] for r in reuse) <= gate["algebra_abs"],
        "channel_bounds_cover_observed_scores": all(
            max(r["continuum_abs_errors"][key] for key in ("assisted", "blind", "jitter"))
            <= r["channel_diamond_bound"] / 2 + gate["algebra_abs"] for r in bridge_rows),
        "collective_dark_and_no_prior_average_advantage": all(
            v["dark_survival"] == 1 and v["equal_site_defect_survival"] >= v["independent_survival"]
            for v in collective_rows),
    }
    shots = gate["laboratory_shots_per_setting_and_arm"]
    contrasts = gate["laboratory_confirmatory_contrasts"]
    sampling = math.sqrt(3 * math.log(2 * contrasts / gate["laboratory_familywise_alpha"]) / (4 * shots))
    return {
        "study": cfg["study"], "registered_checks": checks,
        "all_registered_checks_pass": all(checks.values()),
        "diagnostics": {
            "max_anchor_error": max(anchor_errors), "geometry_max_error": graph_error,
            "maximum_energy_residual": energy, "minimum_marginal_entropy_production": entropy_min,
            "maximum_continuum_error_by_slices": {
                str(n): max(e for r in bridge_rows if r["slices"] == n
                            for e in r["continuum_abs_errors"].values()) for n in cfg["slices"]},
            "hardware_familywise_sampling_halfwidth": sampling,
            "hardware_sampling_plus_systematic_halfwidth": sampling + gate["laboratory_systematic_contrast"],
        },
        "bridge": bridge_rows, "reuse": reuse, "organization": counts,
        "collective": collective_rows, "geometry": graph,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    outputs = parser.add_mutually_exclusive_group(required=True)
    outputs.add_argument("--output", type=Path)
    outputs.add_argument("--verify", type=Path)
    args = parser.parse_args()
    cfg = json.loads((HERE / "PARAMETERS.json").read_text())
    results = calculate(cfg)
    if args.output:
        args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    elif results != json.loads(args.verify.read_text()):
        raise SystemExit("Stored results differ from replay")
    print(json.dumps(results["diagnostics"], indent=2))
    print(f"{len(results['bridge'])} bridge rows, {len(results['geometry'])} graph rows")
    if not results["all_registered_checks_pass"]:
        raise SystemExit(results["registered_checks"])


if __name__ == "__main__":
    main()
