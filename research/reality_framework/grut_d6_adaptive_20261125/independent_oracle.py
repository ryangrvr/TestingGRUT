"""Independent physical-transition oracle for frozen GRUT-D6.

No import of the builder's model. It constructs the reachable Hamiltonian
from local quantum transition matrices, including reverse hops and W's
stationary Hamiltonian, and propagates all 36 registered n=2 conditional
rows with both fair coins. Print JSON, or use --output ORACLE_RESULTS.json.
"""
from __future__ import annotations

import argparse
from collections import deque
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.stats import binom


BLANK = 2
THETA = {"sealed": 0.0, "weak": math.pi/12, "strong": math.pi/4}
ARMS = [
    ("weak", "main", "intact"), ("weak", "main", "scrambled"),
    ("weak", "main", "inverted"), ("weak", "main", "disabled"),
    ("weak", "passive", "intact"), ("weak", "direct", "disabled"),
    ("sealed", "main", "intact"), ("sealed", "direct", "disabled"),
    ("strong", "main", "intact"),
]


def swap_blank(value, color):
    return color if value == BLANK else (BLANK if value == color else value)


def initial_data(n, probe, witness):
    # P,M,A,Q,D,E,W, tape-A, tape-B. Conserved settings are fixed block labels.
    blank = (BLANK,) * n
    return (probe, 0, 0, 0, 0, 0, witness, blank, blank)


def local_transitions(data, link, parameters, inverse=False):
    n, strength, port, action, source, challenge, cp, cm = parameters
    value = list(data)
    if link == 0:
        angle = THETA[strength] * (-1)**source * (-1)**value[0]
        if inverse:
            angle = -angle
        return {data: np.exp(-1j * angle)}
    if link == 1:
        def remaining(v):
            v[1] ^= v[0]
            v[0], v[4] = v[4], v[0]
            if port == "passive":
                v[2] ^= cp
            elif port == "direct" and strength != "sealed":
                v[2] ^= source
            return v

        if inverse:
            if port == "passive":
                value[2] ^= cp
            elif port == "direct" and strength != "sealed":
                value[2] ^= source
            value[0], value[4] = value[4], value[0]
            value[1] ^= value[0]
            other = value.copy()
            other[0] ^= 1
            return {tuple(value): 1/math.sqrt(2), tuple(other): 1j/math.sqrt(2)}
        other = value.copy()
        other[0] ^= 1
        return {tuple(remaining(value)): 1/math.sqrt(2),
                tuple(remaining(other)): -1j/math.sqrt(2)}
    if link == 2:
        if action == "scrambled":
            value[1] ^= cm
        elif action == "inverted":
            value[1] ^= 1
        elif action == "disabled":
            value[1], value[5] = value[5], value[1]
        elif action != "intact":
            raise ValueError(action)
    elif 3 <= link <= n + 2:
        tape = list(value[7])
        color = value[1] if port == "main" else value[2]
        index = link - 3
        tape[index] = swap_blank(tape[index], color)
        value[7] = tuple(tape)
    elif link == 2*n + 3:
        value[3] ^= challenge
    elif 2*n + 4 <= link <= 3*n + 3:
        tape = list(value[8])
        if inverse:
            value[0] ^= 1
        decision = value[1] if port == "main" else value[2]
        color = decision ^ value[3] ^ value[0]
        index = link - (2*n + 4)
        tape[index] = swap_blank(tape[index], color)
        value[8] = tuple(tape)
        if not inverse:
            value[0] ^= 1
    return {tuple(value): 1.0 + 0j}


def physical_graph(parameters):
    n = parameters[0]
    length = 2*(3*n + 4)
    roots = [(0, initial_data(n, p, w)) for p in [0, 1] for w in [0, 1]]
    states = roots.copy()
    indices = {state: i for i, state in enumerate(states)}
    pending = deque(states)
    inverse_error = 0.0
    while pending:
        position, data = pending.popleft()
        neighbors = []
        if position < length:
            transitions = local_transitions(data, position, parameters)
            restored = {}
            for after, amplitude in transitions.items():
                neighbors.append((position + 1, after))
                for original, backward in local_transitions(after, position, parameters, inverse=True).items():
                    restored[original] = restored.get(original, 0j) + amplitude*backward
            inverse_error = max(inverse_error, abs(restored.get(data, 0)-1),
                                max((abs(a) for d, a in restored.items() if d != data), default=0.0))
        if position:
            neighbors.extend((position - 1, after) for after in
                             local_transitions(data, position - 1, parameters, inverse=True))
        for neighbor in neighbors:
            if neighbor not in indices:
                if len(states) > 10000:
                    raise RuntimeError("Unexpected physical graph growth")
                indices[neighbor] = len(states)
                states.append(neighbor)
                pending.append(neighbor)
    omega = 4/(length + 1)
    kappa = omega/4
    witness_sign = 1 if parameters[1] == "sealed" else (-1)**parameters[4]
    row, column, entries = [], [], []
    for source, (position, data) in enumerate(states):
        row.append(source)
        column.append(source)
        entries.append(kappa * witness_sign * (-1)**data[6])
        if position < length:
            coupling = omega/2 * math.sqrt((position + 1)*(length - position))
            for after, amplitude in local_transitions(data, position, parameters).items():
                target = indices[(position + 1, after)]
                row.extend([target, source])
                column.extend([source, target])
                entries.extend([coupling*amplitude, coupling*amplitude.conjugate()])
    ham = coo_matrix((entries, (row, column)), shape=(len(states), len(states))).tocsr()
    initial = np.zeros(len(states), complex)
    for root in roots:
        initial[indices[root]] = 0.5
    return states, indices, ham, initial, inverse_error


def quantum_prefix(parameters, links):
    n = parameters[0]
    vector = {initial_data(n, p, w): 0.5 for p in [0, 1] for w in [0, 1]}
    for link in range(links):
        after = {}
        for data, amplitude in vector.items():
            for target, local in local_transitions(data, link, parameters).items():
                after[target] = after.get(target, 0j) + local*amplitude
        vector = after
    return vector


def mobile_density(vector):
    rho = np.zeros((16, 16), complex)
    for first, amplitude in vector.items():
        i = 8*first[0] + 4*first[1] + 2*first[2] + first[3]
        for second, other in vector.items():
            if first[4:] == second[4:]:
                j = 8*second[0] + 4*second[1] + 2*second[2] + second[3]
                rho[i, j] += amplitude*other.conjugate()
    return rho


def memory_density(rho):
    shaped = rho.reshape(2, 2, 2, 2, 2, 2, 2, 2)
    result = np.zeros((2, 2), complex)
    for p in [0, 1]:
        for a in [0, 1]:
            for q in [0, 1]:
                result += shaped[p, :, a, q, p, :, a, q]
    return result


def trace_distance(first, second):
    return float(abs(np.linalg.eigvalsh(first-second)).sum()/2)


def score_state(vector, states, parameters, region):
    n, _, port, _, source, challenge, _, _ = parameters
    length = 2*(3*n + 4)
    target_a = (source,)*n
    target_b = tuple(source ^ challenge ^ (i % 2) for i in range(n))

    def selected(position):
        if region == "a_only":
            return n+3 <= position <= 2*n+4
        if region == "endpoint":
            return position == length
        if region == "complete":
            return position >= 3*n+4
        raise ValueError(region)

    result = {key: 0.0 for key in ["probability", "memory_zero", "memory_correct", "dump_correct",
                                  "A_success", "B_success", "both_success"]}
    pairs = {}
    for amplitude, (position, data) in zip(vector, states):
        if not selected(position):
            continue
        probability = float(abs(amplitude)**2)
        good_a, good_b = data[7] == target_a, data[8] == target_b
        result["probability"] += probability
        result["memory_zero"] += probability*(data[1] == 0)
        result["memory_correct"] += probability*(data[1] == source)
        result["dump_correct"] += probability*(data[4] == source)
        result["A_success"] += probability*good_a
        result["B_success"] += probability*good_b
        result["both_success"] += probability*good_a*good_b
        rest = (position, data[:6], data[7:])
        if rest not in pairs:
            pairs[rest] = [0j, 0j, data[1] if port == "main" else data[2], good_a and good_b]
        pairs[rest][data[6]] += amplitude
    witness_zero = forecast = joint = 0.0
    for a0, a1, decision, structural in pairs.values():
        mass = abs(a0)**2 + abs(a1)**2
        y_expectation = 2*(a0.conjugate()*a1).imag
        witness_zero += (mass+y_expectation)/2
        agreement = (mass+(-1)**decision*y_expectation)/2
        forecast += agreement
        joint += agreement*structural
    result.update({"W_zero": float(witness_zero), "forecast": float(forecast),
                   "both_with_forecast": float(joint)})
    return result


def energy_law_error(ham, initial, length):
    energies, vectors = eigh(ham.toarray())
    weights = abs(vectors.conj().T @ initial)**2
    omega, kappa = 4/(length+1), 1/(length+1)
    expected = {}
    for k in range(length+1):
        weight = math.comb(length, k)/2**length/2
        for sign in [-1, 1]:
            energy = omega*(k-length/2) + sign*kappa
            expected[energy] = expected.get(energy, 0.0) + weight
    actual = np.zeros(len(expected))
    targets = np.array(sorted(expected))
    spectral_error = 0.0
    for energy, probability in zip(energies, weights):
        nearest = int(abs(targets-energy).argmin())
        spectral_error = max(spectral_error, abs(targets[nearest]-energy))
        actual[nearest] += probability
    return max(spectral_error, float(np.max(abs(actual-np.array([expected[e] for e in targets])))))


def one_coin_run(parameters):
    n, strength, _, _, source, _, _, _ = parameters
    length = 2*(3*n + 4)
    omega = 4/(length+1)
    states, indices, ham, initial, inverse_error = physical_graph(parameters)
    at_a = 2*math.asin(math.sqrt((3*n+7)/(2*length)))/omega
    times = {"a_only": at_a, "endpoint": math.pi/omega,
             "late_minus": (math.pi-0.2)/omega, "late_plus": (math.pi+0.2)/omega}
    scored = {}
    normalization = 0.0
    for label, duration in times.items():
        evolved = expm_multiply(-1j*duration*ham, initial, traceA=0.0)
        region = label if label in ["a_only", "endpoint"] else "complete"
        scored[label] = score_state(evolved, states, parameters, region)
        normalization = max(normalization, abs(float(np.vdot(evolved, evolved).real)-1))
    prefix_before = quantum_prefix(parameters, 2)
    prefix_after = quantum_prefix(parameters, 3)
    before_correct = float(sum(abs(a)**2 for data, a in prefix_before.items() if data[1] == source))
    after_rho = mobile_density(prefix_after)
    after_memory = memory_density(after_rho)
    record_stability = abs(scored["endpoint"]["memory_zero"] - after_memory[0, 0].real)

    # Physical first reset versus full recurrence: W returns only at 4T.
    initial_minus = initial.copy()
    for i, (_, data) in enumerate(states):
        if data[6]:
            initial_minus[i] *= -1
    at_two = expm_multiply(-1j*(2*math.pi/omega)*ham, initial, traceA=0.0)
    at_four = expm_multiply(-1j*(4*math.pi/omega)*ham, initial, traceA=0.0)
    recurrence = {
        "two_T_reset_with_W_minus_fidelity": float(abs(np.vdot(initial_minus, at_two))**2),
        "two_T_full_initial_fidelity": float(abs(np.vdot(initial, at_two))**2),
        "four_T_full_initial_fidelity": float(abs(np.vdot(initial, at_four))**2),
    }
    return scored, after_rho, after_memory, {
        "dimension": len(states), "inverse_error": inverse_error,
        "normalization_error": normalization, "record_stability_error": record_stability,
        "preintervention_record_correct": before_correct,
        "energy_law_error": energy_law_error(ham, initial, length),
        "recurrence": recurrence,
    }, (states, ham, initial)


def run():
    n, length = 2, 20
    rows, mobile, memory = [], {}, {}
    worst = {name: 0.0 for name in ["inverse", "normalization", "energy_law", "record_stability",
                                   "endpoint", "dump_accuracy", "two_task_joint_law", "recurrence"]}
    graph_dimensions = set()
    sealed_cache = {}
    sealed_graph_error = 0.0
    for strength, port, action in ARMS:
        for source in [0, 1]:
            for challenge in [0, 1]:
                averages = None
                rho = np.zeros((16, 16), complex)
                mrho = np.zeros((2, 2), complex)
                before_accuracy = 0.0
                for cp in [0, 1]:
                    for cm in [0, 1]:
                        parameters = (n, strength, port, action, source, challenge, cp, cm)
                        scored, cut_rho, cut_m, diag, graph = one_coin_run(parameters)
                        if averages is None:
                            averages = {time: {key: 0.0 for key in values} for time, values in scored.items()}
                        for time, values in scored.items():
                            for key, value in values.items():
                                averages[time][key] += value/4
                        rho += cut_rho/4
                        mrho += cut_m/4
                        before_accuracy += diag["preintervention_record_correct"]/4
                        graph_dimensions.add(diag["dimension"])
                        for field, label in [("inverse_error", "inverse"), ("normalization_error", "normalization"),
                                             ("energy_law_error", "energy_law"), ("record_stability_error", "record_stability")]:
                            worst[label] = max(worst[label], diag[field])
                        expected_record = (1+math.sin(2*THETA[strength]))/2
                        worst["dump_accuracy"] = max(worst["dump_accuracy"],
                            abs(diag["preintervention_record_correct"]-expected_record),
                            abs(scored["endpoint"]["dump_correct"]-expected_record))
                        worst["endpoint"] = max(worst["endpoint"], abs(scored["endpoint"]["probability"]-1))
                        worst["two_task_joint_law"] = max(worst["two_task_joint_law"],
                            abs(scored["endpoint"]["A_success"]-scored["endpoint"]["both_success"]),
                            abs(scored["endpoint"]["B_success"]-scored["endpoint"]["both_success"]))
                        recur = diag["recurrence"]
                        worst["recurrence"] = max(worst["recurrence"],
                            abs(recur["two_T_reset_with_W_minus_fidelity"]-1),
                            abs(recur["two_T_full_initial_fidelity"]),
                            abs(recur["four_T_full_initial_fidelity"]-1))
                        if strength == "sealed":
                            key = (port, action, challenge, cp, cm)
                            states, ham, initial = graph
                            if source == 0:
                                sealed_cache[key] = (states, ham, initial)
                            else:
                                old_states, old_ham, old_initial = sealed_cache[key]
                                if states != old_states:
                                    raise AssertionError("Sealed source changed the accessible basis")
                                difference = ham-old_ham
                                sealed_graph_error = max(sealed_graph_error,
                                    max(abs(difference.data), default=0), float(np.max(abs(initial-old_initial))))
                key = (strength, port, action, source, challenge)
                mobile[key], memory[key] = rho, mrho
                endpoint = averages["endpoint"]
                row = {
                    "n": n, "strength": strength, "port": port, "intervention": action,
                    "source": source, "challenge": challenge,
                    "a_only_probability": averages["a_only"]["probability"],
                    "a_only_target": averages["a_only"]["A_success"],
                    "endpoint_probability": endpoint["probability"],
                    "endpoint_memory_zero": endpoint["memory_zero"],
                    "preintervention_dump_correct": endpoint["dump_correct"],
                    "preintervention_record_correct": before_accuracy,
                    "endpoint_A_success": endpoint["A_success"],
                    "endpoint_B_success": endpoint["B_success"],
                    "endpoint_both_success": endpoint["both_success"],
                    "endpoint_W_zero": endpoint["W_zero"],
                    "endpoint_forecast": endpoint["forecast"],
                    "endpoint_both_with_forecast": endpoint["both_with_forecast"],
                    "endpoint_memory_correct": endpoint["memory_correct"],
                    "late_minus": averages["late_minus"], "late_plus": averages["late_plus"],
                }
                rows.append(row)
    for row in rows:
        zero = (row["strength"], row["port"], row["intervention"], 0, row["challenge"])
        one = (*zero[:3], 1, zero[4])
        row["complete_response_cut_trace_distance"] = trace_distance(mobile[zero], mobile[one])
        row["memory_cut_trace_distance"] = trace_distance(memory[zero], memory[one])
    if max(worst.values()) > 1e-10 or sealed_graph_error > 1e-10:
        raise AssertionError((worst, sealed_graph_error))
    stage_diagnostics = []
    for size in [2, 64]:
        edges = 2*(3*size+4)
        probability = (3*size+7)/(2*edges)
        q = float(binom.cdf(2*size+4, edges, probability)-binom.cdf(size+2, edges, probability))
        stage_diagnostics.append({"n": size, "A_only_probability": q,
                                  "above_0_9": q > 0.9})
    return {
        "reviewer": "independent internal AI physical-transition implementation, not external RAI or human review",
        "physical_oracle": "complex local transitions and their adjoints; stationary witness term; direct sparse exponentials",
        "conditional_rows": len(rows), "fair_coin_runs": 4*len(rows),
        "reachable_graph_dimensions": sorted(graph_dimensions),
        "all_oracle_checks_pass": True, "worst_absolute_errors": worst,
        "sealed_source_graph_error": sealed_graph_error,
        "stage_diagnostics": stage_diagnostics,
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        print(json.dumps({"output": str(args.output), "conditional_rows": result["conditional_rows"],
                          "all_oracle_checks_pass": result["all_oracle_checks_pass"],
                          "worst_absolute_errors": result["worst_absolute_errors"],
                          "stage_diagnostics": result["stage_diagnostics"]}, sort_keys=True))
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
