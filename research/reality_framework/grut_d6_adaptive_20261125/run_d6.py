#!/usr/bin/env python3
"""Evaluate the frozen D6 model using its exact one-head reduction."""

import argparse
import cmath
import itertools
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
NO_SIGNAL = "sealed"
ANGLES = {"sealed": 0, "weak": math.pi / 12, "strong": math.pi / 4}


def clock_distribution(length, phase):
    turns = phase / math.pi
    if turns.is_integer():
        # Use exact stroboscopic probabilities, including all recurrences.
        index = length if int(turns) % 2 else 0
        return [float(j == index) for j in range(length + 1)]
    arrive = math.sin(phase / 2) ** 2
    stay = math.cos(phase / 2) ** 2
    return [
        math.comb(length, j) * arrive ** j * stay ** (length - j)
        for j in range(length + 1)
    ]


def incomplete_log10(length, first_complete, phase):
    if phase == math.pi:
        return None
    log_arrive = math.log(math.sin(phase / 2) ** 2)
    log_stay = math.log(math.cos(phase / 2) ** 2)
    log_weights = [
        math.log(math.comb(length, j)) + j * log_arrive + (length - j) * log_stay
        for j in range(first_complete)
    ]
    peak = max(log_weights)
    return (peak + math.log(math.fsum(math.exp(x - peak) for x in log_weights))) / math.log(10)


def clock_data(n, offsets):
    first_complete = 3 * n + 4
    length = 2 * first_complete
    phase_a = 2 * math.asin(math.sqrt((3 * n + 7) / (2 * length)))
    stage = clock_distribution(length, phase_a)
    late = {}
    for offset in offsets:
        phase = math.pi + offset
        weights = clock_distribution(length, phase)
        late[str(offset)] = {
            "endpoint": weights[-1],
            "a_complete": math.fsum(weights[n + 3:]),
            "both_complete": math.fsum(weights[first_complete:]),
            "both_incomplete_log10": incomplete_log10(length, first_complete, phase),
        }
    return {
        "n": n,
        "length": length,
        "first_complete": first_complete,
        "phase_a": phase_a,
        "a_only": math.fsum(stage[n + 3:2 * n + 5]),
        "late": late,
        "late_failure_upper_bound": 0.01 * (0.04 ** first_complete),
        "clock_at_2T": clock_distribution(length, 2 * math.pi)[0],
        "clock_at_4T": clock_distribution(length, 4 * math.pi)[0],
    }


def probe_probabilities(strength, source):
    angle = (1 - 2 * source) * ANGLES[strength]
    zero = cmath.exp(-1j * angle)
    one = cmath.exp(1j * angle)
    # Probe starts in |+>. Apply exp(-i pi X/4) before copying it into M.
    rotated = ((zero - 1j * one) / 2, (one - 1j * zero) / 2)
    probabilities = tuple(abs(amplitude) ** 2 for amplitude in rotated)
    total = math.fsum(probabilities)
    return tuple(p / total for p in probabilities)


def ordered_outcomes(strength, port, intervention, source):
    for (dump, probability), passive_coin, mask_coin in itertools.product(
        enumerate(probe_probabilities(strength, source)), (0, 1), (0, 1)
    ):
        if intervention == "intact":
            memory = dump
        elif intervention == "scrambled":
            memory = dump ^ mask_coin
        elif intervention == "inverted":
            memory = dump ^ 1
        elif intervention == "disabled":
            memory = 0
        else:
            raise ValueError(intervention)

        if port == "main":
            decision = memory
        elif port == "passive":
            decision = passive_coin
        elif port == "direct":
            decision = source if strength != NO_SIGNAL else 0
        else:
            raise ValueError(port)

        yield dump, memory, decision, probability / 4


def constructed_words(n, decision, challenge):
    return (
        tuple([decision] * n),
        tuple(decision ^ challenge ^ (i % 2) for i in range(n)),
    )


def source_targets(n, source, challenge):
    even = source ^ challenge
    return (tuple([source] * n), tuple((even, 1 - even)[i % 2] for i in range(n)))


def witness_zero(strength, source, phase):
    # W starts along +X and precesses about Z under its own onsite term.
    sign = 1 if strength == NO_SIGNAL else 1 - 2 * source
    return (1 + sign * math.sin(phase / 2)) / 2


def conditional_row(n, strength, port, intervention, source, challenge, clock, offsets):
    target_a, target_b = source_targets(n, source, challenge)
    outcomes = []
    for dump, memory, decision, p in ordered_outcomes(strength, port, intervention, source):
        word_a, word_b = constructed_words(n, decision, challenge)
        outcomes.append((dump, memory, decision, p, word_a == target_a, word_b == target_b))

    def probability(indicator):
        return math.fsum(p * indicator(dump, memory, decision, a, b)
                         for dump, memory, decision, p, a, b in outcomes)

    terminal = clock["late"]["0"]
    witness_p0 = witness_zero(strength, source, math.pi)
    memory_zero = probability(lambda dump, memory, decision, a, b: memory == 0)
    decision_zero = probability(lambda dump, memory, decision, a, b: decision == 0)
    a_target = probability(lambda dump, memory, decision, a, b: a)
    b_target = probability(lambda dump, memory, decision, a, b: b)
    joint_target = probability(lambda dump, memory, decision, a, b: a and b)
    row = {
        "n": n,
        "strength": strength,
        "port": port,
        "intervention": intervention,
        "source": source,
        "challenge": challenge,
        "dump_correct": probability(lambda dump, memory, decision, a, b: dump == source),
        "memory_zero_at_cut": memory_zero,
        "endpoint_memory_zero": terminal["endpoint"] * memory_zero,
        "endpoint_memory_correct": terminal["endpoint"] * probability(
            lambda dump, memory, decision, a, b: memory == source
        ),
        "decision_zero": decision_zero,
        "endpoint": terminal["endpoint"],
        "endpoint_a_target": terminal["endpoint"] * a_target,
        "endpoint_b_target": terminal["endpoint"] * b_target,
        "endpoint_joint_target": terminal["endpoint"] * joint_target,
        "endpoint_witness_zero": witness_p0,
        "endpoint_forecast": terminal["endpoint"] * probability(
            lambda dump, memory, decision, a, b: witness_p0 if decision == 0 else 1 - witness_p0
        ),
        "endpoint_joint_with_forecast": terminal["endpoint"] * probability(
            lambda dump, memory, decision, a, b: (a and b) *
            (witness_p0 if decision == 0 else 1 - witness_p0)
        ),
        "a_only": clock["a_only"],
        "a_only_target": clock["a_only"] * a_target,
    }
    for label, offset in zip(("left", "right"), (offsets[0], offsets[-1]), strict=True):
        phase = math.pi + offset
        completed = clock["late"][str(offset)]["both_complete"]
        match_witness = witness_zero(strength, source, phase)
        if source:
            match_witness = 1 - match_witness
        row[f"late_{label}_joint_target"] = completed * joint_target
        row[f"late_{label}_joint_with_forecast"] = completed * joint_target * match_witness
    return row


def resources(clock):
    length = clock["length"]
    n = clock["n"]
    omega = 4 / (length + 1)
    kappa = omega / 4
    law = [
        (omega * (k - length / 2) + sign * kappa, weight / 2)
        for k, weight in ((k, math.comb(length, k) / 2 ** length) for k in range(length + 1))
        for sign in (-1, 1)
    ]
    mean = math.fsum(energy * weight for energy, weight in law)
    variance = math.fsum((energy - mean) ** 2 * weight for energy, weight in law)
    return {
        "n": n,
        "omega": omega,
        "kappa": kappa,
        "T": math.pi / omega,
        "link_norm_max": omega / 2 * math.sqrt((length // 2) * (length // 2 + 1)),
        "energy_law_norm": math.fsum(weight for energy, weight in law),
        "energy_mean": mean,
        "energy_variance": variance,
        "energy_variance_exact_formula": omega ** 2 * length / 4 + kappa ** 2,
        "energy_levels": len(law),
        "display_qutrits": 2 * n,
        "fair_coin_entropy_bits_conditional_on_b_and_z": 2,
        "full_dimension_log2": (length + 1) * math.log2(49) + math.log2(1536)
        + 2 * n * math.log2(3),
        "one_head_dimension_log2": math.log2((length + 1) * 48 * 1536)
        + 2 * n * math.log2(3),
        "work_during_episode": 0,
    }


def entropy_bit(p):
    return -math.fsum(x * math.log2(x) for x in (p, 1 - p) if x > 0)


def validate(protocol, rows, clocks, resource_rows):
    gates = protocol["numerical_gates"]
    forecasts = protocol["endpoint_forecasts"]
    tol = gates["endpoint_distribution_abs_error"]
    errors = []
    mismatches = []

    def match(actual, expected):
        errors.append(abs(actual - expected))

    def gate(condition, detail):
        if not condition:
            mismatches.append(detail)

    groups = {}
    for r in rows:
        n, strength, port, intervention, b = (r[k] for k in
                                             ("n", "strength", "port", "intervention", "source"))
        groups.setdefault((n, strength, port, intervention), []).append(r)
        p = forecasts["preintervention_record_correct_for_either_b"][strength]
        expected_m = {"intact": p, "scrambled": 0.5,
                      "inverted": 1 - p, "disabled": float(b == 0)}[intervention]
        expected_decision = {
            "main": expected_m, "passive": 0.5,
            "direct": float(strength != NO_SIGNAL or b == 0),
        }[port]
        expected_w0 = float(strength == NO_SIGNAL or b == 0)
        expected_d0 = expected_decision if b == 0 else 1 - expected_decision
        match(r["endpoint"], 1)
        match(r["dump_correct"], p)
        match(r["endpoint_memory_correct"], expected_m)
        match(r["endpoint_memory_zero"], expected_m if b == 0 else 1 - expected_m)
        match(r["endpoint_memory_zero"], r["memory_zero_at_cut"])
        match(r["decision_zero"], expected_d0)
        match(r["endpoint_witness_zero"], expected_w0)
        for metric in ("endpoint_a_target", "endpoint_b_target", "endpoint_joint_target"):
            match(r[metric], expected_decision)
        match(r["endpoint_forecast"], expected_decision if strength != NO_SIGNAL else expected_d0)
        expected_joint = expected_decision if strength != NO_SIGNAL or b == 0 else 0
        match(r["endpoint_joint_with_forecast"], expected_joint)
        if strength == "weak" and port == "main" and intervention == "intact":
            for side in ("left", "right"):
                gate(r[f"late_{side}_joint_with_forecast"]
                     >= gates["late_weak_main_intact_joint_with_forecast_min"],
                     [n, b, side, "weak_late_forecast"])
    balances = {}
    trace_distances = []
    for key, group in groups.items():
        n, strength, port, intervention = key
        rank = protocol["production_grid"]["intervention"].index(intervention)
        p = math.fsum(r["endpoint_joint_target"] for r in group) / len(group)
        p_forecast = math.fsum(r["endpoint_forecast"] for r in group) / len(group)
        balances[key] = p
        if port == "main":
            expected = forecasts["main_balanced_joint_structural_in_intervention_order"][strength][rank]
        elif port == "passive":
            expected = forecasts["passive_balanced_joint_structural_all"]
        else:
            expected = forecasts["direct_balanced_joint_structural_" +
                                 ("sealed" if strength == NO_SIGNAL else "open")]
        match(p, expected)
        if strength == NO_SIGNAL and port == "main" and intervention == "disabled":
            match(p_forecast, forecasts["silent_disabled_main_forecast_W_balanced"])
        if port == "main":
            by_b = {r["source"]: r["memory_zero_at_cut"] for r in group if r["challenge"] == 0}
            distance = abs(by_b[0] - by_b[1])
            forecast_key = f"memory_cut_trace_distance_in_intervention_order_{strength}_main"
            target = forecasts[forecast_key][rank] if strength != NO_SIGNAL else 0.0
            match(distance, target)
            if n == protocol["production_grid"]["n"][0]:
                trace_distances.append([strength, intervention, distance])

    stage_gate_rows = []
    for c in clocks:
        n = c["n"]
        active = balances[n, "weak", "main", "intact"]
        weak_forecast = min(r["endpoint_joint_with_forecast"] for r in
                            groups[n, "weak", "main", "intact"])
        gate(active >= gates["weak_main_intact_joint_structural_min"], [n, "weak_targets"])
        gate(weak_forecast >= gates["weak_main_intact_joint_with_forecast_min"], [n, "weak_forecast"])
        for control in (("passive", "intact"), ("main", "scrambled"), ("main", "disabled")):
            advantage = active - balances[(n, "weak", *control)]
            gate(advantage >= gates["weak_causal_advantage_over_each_of_passive_intact_scrambled_disabled_min"],
                 [n, *control, "weak_causal_advantage"])
        strong = balances[n, "strong", "main", "intact"]
        gate(strong >= gates["strong_main_intact_joint_structural_min"], [n, "strong_targets"])
        for strength, intervention in itertools.product(("weak", "strong"),
                                                        protocol["production_grid"]["intervention"]):
            gate(balances[n, strength, "direct", intervention] >=
                 gates["direct_open_joint_structural_min_in_all_interventions"],
                 [n, strength, intervention, "direct_targets"])
        for key, score in balances.items():
            if key[0] == n and key[1] == NO_SIGNAL:
                match(score, gates["sealed_balanced_structural_score"])
        gate(c["late_failure_upper_bound"] <= gates["late_incomplete_bound_max"], [n, "late_tail_bound"])
        for sample in c["late"].values():
            tail = sample["both_incomplete_log10"]
            if tail is not None:
                gate(tail <= math.log10(c["late_failure_upper_bound"]), [n, "late_tail"])
        above = c["a_only"] >= gates["a_only_probability_diagnostic_threshold"]
        if n in gates["a_only_expected_below_threshold_n"]:
            gate(not above, [n, "expected_short_clock_miss"])
        if n in gates["a_only_required_above_threshold_n"]:
            gate(above, [n, "required_long_clock_pass"])
        match(c["clock_at_2T"], 1)
        match(c["clock_at_4T"], 1)
        stage_gate_rows.append([n, c["a_only"], above])

    for r in resource_rows:
        match(r["energy_mean"], 0)
        match(r["energy_law_norm"], 1)
        match(r["energy_variance"], r["energy_variance_exact_formula"])
        gate(r["link_norm_max"] <= 1, [r["n"], "link_norm"])
    residual = max(errors)
    gate(residual <= tol, ["endpoint_and_resource_residual", residual])
    if mismatches:
        raise AssertionError(mismatches)
    return {
        "max_abs_check_residual": residual,
        "registered_checks_pass": True,
        "a_only_diagnostics_n_probability_threshold_pass": stage_gate_rows,
        "main_trace_distances_at_cut": trace_distances,
        "balanced_n_strength_port_intervention_joint": [
            [*key, value] for key, value in balances.items()
            if key[0] == protocol["production_grid"]["n"][0]
        ],
        "memory_information_bits": {
            s: 1 - entropy_bit(p)
            for s, p in forecasts["preintervention_record_correct_for_either_b"].items()
        },
        "largest_late_incomplete_upper_bound": max(c["late_failure_upper_bound"] for c in clocks),
        "min_late_weak_main_intact_joint_with_forecast": min(
            r[k] for r in rows if (r["strength"], r["port"], r["intervention"])
            == ("weak", "main", "intact")
            for k in ("late_left_joint_with_forecast", "late_right_joint_with_forecast")
        ),
    }


def calculate():
    protocol = json.loads((HERE / "PROTOCOL.json").read_text())
    grid = protocol["production_grid"]
    offsets = protocol["timings"]["late_omega_t_offsets_from_pi"]
    clocks = {n: clock_data(n, offsets) for n in grid["n"]}
    rows = [
        conditional_row(n, strength, port, intervention, source, challenge, clocks[n], offsets)
        for n, strength, port, intervention, source, challenge in itertools.product(
            grid["n"], grid["strength"], grid["port"], grid["intervention"],
            grid["source"], grid["stage_challenge"]
        )
    ]
    resource_rows = [resources(c) for c in clocks.values()]
    checks = validate(protocol, rows, list(clocks.values()), resource_rows)
    return {
        "study": "GRUT-D6",
        "solver": "analytic one-head gauge with explicit probe and witness probabilities",
        "checks": checks,
        "clock": list(clocks.values()),
        "resources": resource_rows,
        "columns": list(rows[0]),
        "conditional_rows": [list(r.values()) for r in rows],
    }


def serialize(results):
    metadata = {k: v for k, v in results.items() if k != "conditional_rows"}
    prefix = json.dumps(metadata, indent=2, allow_nan=False)
    row_block = ",\n".join("    " + json.dumps(row, allow_nan=False)
                           for row in results["conditional_rows"])
    return prefix[:-2] + ',\n  "conditional_rows": [\n' + row_block + '\n  ]\n}\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    destination = parser.add_mutually_exclusive_group(required=True)
    destination.add_argument("--output", type=Path)
    destination.add_argument("--verify", type=Path)
    args = parser.parse_args()
    results = calculate()
    serialized = serialize(results)
    if args.output is not None:
        args.output.write_text(serialized)
    elif json.loads(serialized) != json.loads(args.verify.read_text()):
        raise SystemExit(f"{args.verify} differs from the computed results")
    print(f"D6: {len(results['conditional_rows'])} conditional rows, all registered checks passed.")
    checks = results["checks"]
    print(f"Maximum forecast/resource residual: {checks['max_abs_check_residual']:.3g}")
    print(f"Weak late joint score, minimum: {checks['min_late_weak_main_intact_joint_with_forecast']:.12f}")
    print("A-only stage diagnostics [n, probability, threshold pass]:")
    print(checks["a_only_diagnostics_n_probability_threshold_pass"])


if __name__ == "__main__":
    main()
