"""Independent reachable physical-graph oracle for GRUT-D5.

No import of the builder's implementation. Local reversible head/tape hops
are implemented directly from MODEL.md, including inverses and pulse closure.
Small L=64/128 tests are implementation references, not the frozen large-L
external-success certificate. The program prints its own JSON review record.
"""
from __future__ import annotations

from collections import deque
import json
import math

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply


BLANK = 2
CONTROL_NAMES = (
    "intact", "G inverted", "G scrambled", "G disabled",
    "M inverted", "M scrambled", "M disabled",
)


def invert(value):
    return 1 - value if value != BLANK else BLANK


def swap_blank(value, color):
    if color == BLANK:
        return value
    if value == BLANK:
        return color
    return BLANK if value == color else value


def initial_data(n, h, g, coin=0, fault=0, control=0):
    # h,G,M,P,c,E_G,E_M,F,K; only data, not an external truth register.
    head = (h, g, BLANK, 0, coin, BLANK, BLANK, fault, control)
    tape = (BLANK,) * n
    return head, tape, tape, tape


def instruction(n, length, edge):
    quarter, half = n // 4, n // 2
    if edge == 0:
        return "assay", 0
    if 1 <= edge <= quarter:
        return "grow", edge - 1
    if edge == length // 4:
        return "fault1", 0
    start = length // 4 + 1
    if start <= edge < start + quarter:
        return "grow", quarter + edge - start
    if edge == length // 2:
        return "central", 1
    start = length // 2 + 1 + half
    if start <= edge < start + n:
        return "rebuild", edge - start
    if edge == 3 * length // 4:
        return "fault3", 0
    start = 3 * length // 4 + 1
    if start <= edge < start + n:
        return "maintain", n - 1 - (edge - start)
    return "idle", 0


def apply_control(head):
    mode = head[8]
    if mode in [1, 2, 4, 5]:
        if mode in [1, 4] or head[4] == 1:
            index = 1 if mode in [1, 2] else 2
            head[index] = invert(head[index])
    elif mode in [3, 6]:
        index, backup = (1, 5) if mode == 3 else (2, 6)
        head[index], head[backup] = head[backup], head[index]


def local_hop(data, n, length, edge, inverse=False, pilot_null=False):
    old_head, old_a, old_b, old_c = data
    kind, site = instruction(n, length, edge)
    if kind == "idle":
        return data
    head, a, b, c = list(old_head), list(old_a), list(old_b), list(old_c)

    def color():
        source = head[0] if pilot_null else head[2]
        return BLANK if source == BLANK or head[1] == BLANK else source ^ (head[1] & head[3])

    if kind == "assay":
        head[2] = swap_blank(head[2], head[0])
    elif kind == "grow":
        if inverse:
            head[3] ^= 1
        a[site] = swap_blank(a[site], color())
        if not inverse:
            head[3] ^= 1
    elif kind == "rebuild":
        if inverse:
            head[3] ^= 1
            b[site], a[site] = a[site], b[site]
            b[site] = swap_blank(b[site], color())
        else:
            b[site] = swap_blank(b[site], color())
            b[site], a[site] = a[site], b[site]
            head[3] ^= 1
    elif kind == "maintain":
        if inverse:
            c[site], a[site] = a[site], c[site]
            c[site] = swap_blank(c[site], color())
            head[3] ^= 1
        else:
            head[3] ^= 1
            c[site] = swap_blank(c[site], color())
            c[site], a[site] = a[site], c[site]
    elif kind == "fault1":
        if head[7] == 1:
            a[0] = invert(a[0])
    elif kind == "fault3":
        if head[7] == 3:
            a[0] = invert(a[0])
    elif kind == "central":
        if inverse:
            apply_control(head)
        if head[7] == 2:
            a[1] = invert(a[1])
        if not inverse:
            apply_control(head)
    else:
        raise AssertionError(kind)
    return tuple(head), tuple(a), tuple(b), tuple(c)


def flip_display(data, site):
    head, a, b, c = data
    a = list(a)
    a[site] = invert(a[site])
    return head, tuple(a), b, c


def ordered_history(data, n, length, pilot_null=False):
    path = [data]
    for edge in range(length):
        after = local_hop(path[-1], n, length, edge, pilot_null=pilot_null)
        if local_hop(after, n, length, edge, inverse=True, pilot_null=pilot_null) != path[-1]:
            raise AssertionError("Local gate failed its inverse")
        path.append(after)
    return path


def physical_graph(data, n, length, pulse_site=None, pilot_null=False):
    initial = (0, data)
    states, index, pending = [initial], {initial: 0}, deque([initial])
    while pending:
        position, current = pending.popleft()
        neighbors = []
        if position < length:
            neighbors.append((position + 1, local_hop(current, n, length, position,
                                                     pilot_null=pilot_null)))
        if position:
            neighbors.append((position - 1, local_hop(current, n, length, position - 1,
                                                     inverse=True, pilot_null=pilot_null)))
        if pulse_site is not None:
            neighbors.append((position, flip_display(current, pulse_site)))
        for neighbor in neighbors:
            if neighbor not in index:
                if len(states) >= 100000:
                    raise RuntimeError("Unexpectedly large physical pulse closure")
                index[neighbor] = len(states)
                states.append(neighbor)
                pending.append(neighbor)
    omega = 4 / (length + 1)
    row, column, values = [], [], []
    for source, (position, current) in enumerate(states):
        if position == length:
            continue
        target_state = (position + 1, local_hop(current, n, length, position,
                                               pilot_null=pilot_null))
        target = index[target_state]
        coupling = omega / 2 * math.sqrt((position + 1) * (length - position))
        row.extend([source, target])
        column.extend([target, source])
        values.extend([coupling, coupling])
    h = coo_matrix((values, (row, column)), shape=(len(states), len(states))).tocsr()
    permutation = None
    if pulse_site is not None:
        permutation = np.array([index[(position, flip_display(current, pulse_site))]
                                for position, current in states])
        if not np.array_equal(permutation[permutation], np.arange(len(states))):
            raise AssertionError("Pulse closure is not an involution")
    return states, index, h, permutation


def target_word(n, h, g):
    return tuple(h ^ (g & (site % 2)) for site in range(n))


def stage_band(n, length, trial):
    if trial == 0:
        return n // 4 + 1, length // 4
    if trial == 1:
        return length // 4 + n // 4 + 1, length // 2
    return length // 2 + n // 2 + n + 1, 3 * length // 4


def external_oracle(n, length, h, g, trial, pilot_null=False):
    probability, site = [(1/8, 0), (3/8, 1), (5/8, 0)][trial]
    data = initial_data(n, h, g)
    path = ordered_history(data, n, length, pilot_null)
    states, indices, ham, permutation = physical_graph(data, n, length, site, pilot_null)
    omega = 4 / (length + 1)
    pulse_time = 2 * math.asin(math.sqrt(probability)) / omega
    terminal = math.pi / omega
    initial = np.zeros(len(states), complex)
    initial[indices[(0, data)]] = 1
    before = expm_multiply(-1j * pulse_time * ham, initial, traceA=0.0)
    pulsed = np.zeros_like(before)
    pulsed[permutation] = before
    final = expm_multiply(-1j * (terminal - pulse_time) * ham, pulsed, traceA=0.0)

    # Independent PST endpoint identity: each position contributes its
    # binomial PROBABILITY as an amplitude, with one common final phase.
    weights = np.array([math.comb(length, j) * probability**j
                        * (1 - probability)**(length - j) for j in range(length + 1)])
    groups = {}
    for position in range(length + 1):
        current = flip_display(path[position], site)
        for edge in range(position, length):
            current = local_hop(current, n, length, edge, pilot_null=pilot_null)
        groups[current] = groups.get(current, 0.0) + weights[position]
    phase = (-1j)**length
    predicted = np.zeros(len(states), complex)
    for current, weight in groups.items():
        predicted[indices[(length, current)]] = phase * weight
    end_indices = [i for i, (position, _) in enumerate(states) if position == length]
    amplitude_error = float(np.max(abs(final[end_indices] - predicted[end_indices])))
    goal = target_word(n, h, g)
    archive, archive_site = (2, site) if trial < 2 else (3, site)

    def successful(current):
        return current[1] == goal and current[archive][archive_site] == invert(goal[archive_site])

    direct_success = float(sum(abs(final[i])**2 for i in end_indices if successful(states[i][1])))
    grouped_success = float(sum(weight**2 for current, weight in groups.items() if successful(current)))
    direct_end = float(np.sum(abs(final[end_indices])**2))
    grouped_end = float(sum(weight**2 for weight in groups.values()))
    lower, upper = stage_band(n, length, trial)
    band_mass = float(weights[lower:upper + 1].sum()) if lower <= upper else 0.0
    direct_stage = float(sum(abs(before[i])**2 for i, (position, _) in enumerate(states)
                             if lower <= position <= upper))
    before_energy = float(np.vdot(before, ham @ before).real)
    after_energy = float(np.vdot(pulsed, ham @ pulsed).real)
    variance_before = float(np.linalg.norm(ham @ before)**2 - before_energy**2)
    variance_after = float(np.linalg.norm(ham @ pulsed)**2 - after_energy**2)
    return {
        "n": n, "length": length, "h": h, "g": g, "trial": trial,
        "pilot_null": pilot_null, "physical_reachable_dimension": len(states),
        "hamiltonian_edges_max_norm": float(abs(ham.data).max()),
        "endpoint_amplitude_error": amplitude_error,
        "direct_head_endpoint_probability": direct_end,
        "grouped_head_endpoint_probability": grouped_end,
        "direct_target_and_archive_success": direct_success,
        "grouped_target_and_archive_success": grouped_success,
        "canonical_stage_band": [lower, upper],
        "binomial_stage_probability": band_mass,
        "direct_stage_probability": direct_stage,
        "mean_pulse_work": after_energy - before_energy,
        "energy_variance_before": variance_before,
        "energy_variance_after": variance_after,
        "normalization_error": float(abs(np.vdot(final, final).real - 1)),
        "complete_endpoint_groups": len(groups),
    }


def clock_energy_test(length):
    omega = 4 / (length + 1)
    edge = omega / 2 * np.sqrt(np.arange(1, length + 1) * np.arange(length, 0, -1))
    clock = np.diag(edge, 1) + np.diag(edge, -1)
    eigenvalues, vectors = eigh(clock)
    expected_energy = omega * (np.arange(length + 1) - length / 2)
    expected_law = np.array([math.comb(length, j) / 2**length for j in range(length + 1)])
    return {
        "length": length,
        "spectrum_error": float(np.max(abs(eigenvalues - expected_energy))),
        "initial_spectral_law_error": float(np.max(abs(abs(vectors[0])**2 - expected_law))),
        "mean_energy": float(expected_energy @ expected_law),
        "variance": float(expected_energy**2 @ expected_law),
        "expected_variance": 4*length/(length + 1)**2,
        "spectral_norm": float(abs(eigenvalues).max()),
    }


def internal_graph_oracle(n, length, coin, control, pilot_null=False):
    data = initial_data(n, 0, 1, coin, fault=2, control=control)
    path = ordered_history(data, n, length, pilot_null)
    states, indices, ham, _ = physical_graph(data, n, length, pilot_null=pilot_null)
    initial = np.zeros(len(states), complex)
    initial[indices[(0, data)]] = 1
    final = expm_multiply(-1j * math.pi * (length + 1) / 4 * ham, initial, traceA=0.0)
    target = target_word(n, 0, 1)
    endpoint = indices[(length, path[-1])]
    expected_success = int(path[-1][1] == target)
    success = float(sum(abs(final[index])**2 for index, (position, current) in enumerate(states)
                        if position == length and current[1] == target))
    return {"coin": coin, "control": CONTROL_NAMES[control], "pilot_null": pilot_null,
            "reachable_dimension": len(states), "expected_success": expected_success,
            "direct_success": success, "success_error": abs(success - expected_success),
            "endpoint_probability_error": abs(abs(final[endpoint])**2 - 1)}


def main():
    n, length = 8, 64
    internal = []
    observational_mismatches = 0
    for h in [0, 1]:
        for g in [0, 1]:
            target = target_word(n, h, g)
            for coin in [0, 1]:
                for fault in range(4):
                    initial = initial_data(n, h, g, coin, fault)
                    primary = ordered_history(initial, n, length)
                    null = ordered_history(initial, n, length, pilot_null=True)
                    observational_mismatches += sum(a != b for a, b in zip(primary, null))
                    assert primary[-1][1] == target
                    assert primary[length//4][0][2] == h
                    assert primary[length//2][0][2] == h
                    assert sum(x != BLANK for x in primary[length//4][1]) == n//4
                    assert sum(x != BLANK for x in primary[length//2][1]) == n//2
                    last_rebuild = length//2 + n//2 + n
                    assert BLANK in primary[last_rebuild][1]
                    assert primary[last_rebuild + 1][1] == target
                    if fault:
                        archive, cell = (2, 0) if fault == 1 else ((2, 1) if fault == 2 else (3, 0))
                        assert primary[-1][archive][cell] == invert(target[cell])
                for mode, name in enumerate(CONTROL_NAMES):
                    initial = initial_data(n, h, g, coin, fault=2, control=mode)
                    primary = ordered_history(initial, n, length)
                    null = ordered_history(initial, n, length, pilot_null=True)
                    intact_expected = mode == 0 or (mode in [2, 5] and coin == 0)
                    assert (primary[-1][1] == target) == intact_expected
                    if mode in [0, 4, 5, 6]:
                        assert null[-1][1] == target
                    internal.append({"h": h, "g": g, "coin": coin, "control": name,
                                     "main_target_success": int(primary[-1][1] == target),
                                     "null_target_success": int(null[-1][1] == target)})
    assert observational_mismatches == 0

    internal_graph_cases = [internal_graph_oracle(n, length, 1, mode, null)
                            for mode in range(7) for null in [False, True]]
    internal_graph_cases += [internal_graph_oracle(n, length, 0, mode, null)
                             for mode in [2, 5] for null in [False, True]]

    # Direct graph evolution includes arbitrary pulse closure. L=64 is an
    # intentionally small implementation test; the third pulse has broad
    # clock support across construction stages and need not meet large-L gates.
    external = []
    for h in [0, 1]:
        for g in [0, 1]:
            for trial in range(3):
                for null in [False, True]:
                    external.append(external_oracle(n, length, h, g, trial, null))
    for trial in range(3):
        external.append(external_oracle(n, 128, 0, 1, trial))

    worst_amplitude = max(case["endpoint_amplitude_error"] for case in external)
    worst_probability = max(abs(case["direct_target_and_archive_success"]
                                - case["grouped_target_and_archive_success"]) for case in external)
    worst_stage = max(abs(case["direct_stage_probability"]
                          - case["binomial_stage_probability"]) for case in external)
    worst_normalization = max(case["normalization_error"] for case in external)
    null_error = 0.0
    for i in range(0, 24, 2):
        null_error = max(null_error, abs(external[i]["direct_target_and_archive_success"]
                                          - external[i + 1]["direct_target_and_archive_success"]))
    energy_tests = [clock_energy_test(value) for value in [64, 128]]
    assert max(worst_amplitude, worst_probability, worst_stage, worst_normalization, null_error) < 1e-10
    assert max(case["success_error"] for case in internal_graph_cases) < 1e-10
    assert max(case["endpoint_probability_error"] for case in internal_graph_cases) < 1e-10
    assert max(test["spectrum_error"] for test in energy_tests) < 1e-10
    assert max(test["initial_spectral_law_error"] for test in energy_tests) < 1e-10
    print(json.dumps({
        "reviewer": "independent internal AI physical-graph implementation; no external RAI or human review",
        "small_oracle_scope": "N8, L64/128 tests implementation and clock tails, not frozen large-L success thresholds",
        "all_oracle_checks_pass": True,
        "internal_control_cases": internal,
        "internal_direct_graph_cases": internal_graph_cases,
        "intact_complete_history_null_mismatches": observational_mismatches,
        "external_physical_graph_cases": external,
        "worst_endpoint_amplitude_error": worst_amplitude,
        "worst_target_archive_probability_error": worst_probability,
        "worst_stage_probability_error": worst_stage,
        "worst_normalization_error": worst_normalization,
        "worst_null_external_success_error": null_error,
        "clock_energy_checks": energy_tests,
        "prospective_bound": {
            "independent_uniform_future_phase_and_rule_whole_word_upper": 0.25,
            "independent_uniform_future_rule_with_known_phase_whole_word_upper": 0.5,
            "known_phase_uniform_vs_alternating_per_cell_baseline": 0.75,
            "scope": "future target independent and unavailable before scored action; no anticipation inference from overlapping partial rewards",
        },
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
