"""Independent physical tensor-Hamiltonian and coding-bound audit for D4.

Physical matrices include qutrit blanks, W and every reservoir spin. The
ready-model reference is a separate three-physical-bit tensor construction,
not the builder's fermion-mode Gram calculation. Prints a reproducible JSON
review record; represents internal AI review, not external RAI execution.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.linalg import eigh

import model


DELTA = 4.0
KAPPA = 1 / 65536
LOWER = np.array([[0, 1], [0, 0]], complex)
NUMBER = LOWER.conj().T @ LOWER
MASKS = {"S": (1, 1), "M1": (1, 0), "M2": (0, 1)}
LOGICAL = {"S": (1, 1), "M1": (0, 1), "M2": (1, 0)}


def majority(bits):
    return int(sum(bits) >= 2)


def encode(bits):
    s, first, second = bits
    return majority(bits), s ^ first, s ^ second


def decode(logical, first, second):
    s = logical ^ (first & second)
    return s, s ^ first, s ^ second


def flat_core(s, first, second, writer):
    return ((3 * s + first) * 3 + second) * 2 + writer


def core_operators():
    blank = np.zeros((36, 36), complex)
    writer_number = np.zeros_like(blank)
    parity = [np.zeros_like(blank), np.zeros_like(blank)]
    lowering = [np.zeros_like(blank), np.zeros_like(blank)]
    flips = {site: np.zeros_like(blank) for site in MASKS}
    logical_z = np.zeros_like(blank)
    exchange = np.zeros_like(blank)
    for s in range(2):
        for first in range(3):
            for second in range(3):
                for w in range(2):
                    source = flat_core(s, first, second, w)
                    blank[source, source] = first == 2 and second == 2
                    writer_number[source, source] = w
                    ready = first < 2 and second < 2
                    logical = majority((s, first, second)) if ready else s
                    logical_z[source, source] = (-1)**logical
                    if ready:
                        r, x1, x2 = encode((s, first, second))
                        for channel, x in enumerate((x1, x2)):
                            parity[channel][source, source] = x
                            if x:
                                changed = [x1, x2]
                                changed[channel] = 0
                                ss, aa, bb = decode(r, *changed)
                                lowering[channel][flat_core(ss, aa, bb, w), source] = 1
                    for site in MASKS:
                        changed = [s, first, second]
                        index = {"S": 0, "M1": 1, "M2": 2}[site]
                        if changed[index] != 2:
                            changed[index] ^= 1
                        flips[site][flat_core(*changed, w), source] = 1
        source = flat_core(s, 2, 2, 0)
        target = flat_core(s, s, s, 1)
        exchange[source, target] = exchange[target, source] = 1
    return blank + writer_number, parity, lowering, exchange, flips, logical_z


def spin_operator(length, site, operator):
    result = np.ones((1, 1), complex)
    for j in range(length):
        result = np.kron(result, operator if j == site else np.eye(2))
    return result


def assemble(core_q, parity, lowering, writer_exchange, lengths, core_dimension):
    total_length = sum(lengths)
    bath_dimension = 2**total_length
    bath_identity = np.eye(bath_dimension)
    core_identity = np.eye(core_dimension)
    h0 = DELTA * np.kron(core_q + parity[0] + parity[1], bath_identity)
    bath_counts = []
    offset = 0
    for channel, length in enumerate(lengths):
        count = np.zeros((bath_dimension, bath_dimension), complex)
        for j in range(length):
            count += spin_operator(total_length, offset + j, NUMBER)
        bath_counts.append(count)
        h0 += DELTA * np.kron(core_identity, count)
        exchange = np.kron(lowering[channel], spin_operator(total_length, offset, LOWER.conj().T))
        h0 += exchange + exchange.conj().T
        for j in range(length - 1):
            forward = (spin_operator(total_length, offset + j, LOWER.conj().T)
                       @ spin_operator(total_length, offset + j + 1, LOWER))
            h0 += np.kron(core_identity, forward + forward.conj().T)
        offset += length
    writer_term = KAPPA * np.kron(writer_exchange, bath_identity)
    return h0 + writer_term, h0, bath_counts


def physical_hamiltonians(lengths):
    qw, parity, lowering, writer, flips, logical_z = core_operators()
    h, h0, counts = assemble(qw, parity, lowering, writer, lengths, 36)
    bath_identity = np.eye(2**sum(lengths))
    qi = [np.kron(parity[i], bath_identity) + np.kron(np.eye(36), counts[i])
          for i in range(2)]
    return h, h0, np.kron(qw, bath_identity), qi, {
        site: np.kron(flip, bath_identity) for site, flip in flips.items()
    }, np.kron(logical_z, bath_identity)


def ready_hamiltonian(lengths):
    parity = [np.zeros((8, 8), complex), np.zeros((8, 8), complex)]
    lowering = [np.zeros((8, 8), complex), np.zeros((8, 8), complex)]
    flips = {site: np.zeros((8, 8), complex) for site in MASKS}
    for source in range(8):
        bits = [(source >> j) & 1 for j in [2, 1, 0]]
        logical, x1, x2 = encode(bits)
        for channel, x in enumerate((x1, x2)):
            parity[channel][source, source] = x
            if x:
                changed = [x1, x2]
                changed[channel] = 0
                ss, aa, bb = decode(logical, *changed)
                lowering[channel][4*ss + 2*aa + bb, source] = 1
        for site in MASKS:
            changed = bits.copy()
            changed[{"S": 0, "M1": 1, "M2": 2}[site]] ^= 1
            flips[site][4*changed[0] + 2*changed[1] + changed[2], source] = 1
    zero = np.zeros((8, 8), complex)
    h, _, _ = assemble(zero, parity, lowering, zero, lengths, 8)
    bath_identity = np.eye(2**sum(lengths))
    return h, {site: np.kron(flip, bath_identity) for site, flip in flips.items()}


class DenseEvolution:
    def __init__(self, h):
        self.energy, self.vectors = eigh(h)

    def step(self, vector, elapsed):
        return self.vectors @ (np.exp(-1j * self.energy * elapsed)
                               * (self.vectors.conj().T @ vector))

    def sequence(self, initial, pulses, terminal, flips):
        state, previous = initial.copy(), 0.0
        for time, sites in pulses:
            state = self.step(state, time - previous)
            for site in sites:
                state = flips[site] @ state
            previous = time
        return self.step(state, terminal - previous)


def ready_joint(vector, lengths):
    L1, L2 = lengths
    joint = np.zeros((2, 2, 2, L1 + 1, L2 + 1))
    probabilities = abs(vector.reshape(8, -1))**2
    for code in range(8):
        logical, x1, x2 = encode([(code >> j) & 1 for j in [2, 1, 0]])
        for bath, probability in enumerate(probabilities[code]):
            first = (bath >> L2).bit_count()
            second = (bath & ((1 << L2) - 1)).bit_count()
            joint[logical, x1, x2, first, second] += probability
    return joint


def compare_padded(first, second):
    shape = tuple(max(a, b) for a, b in zip(first.shape, second.shape))
    a = np.zeros(shape)
    b = np.zeros(shape)
    a[tuple(slice(0, value) for value in first.shape)] = first
    b[tuple(slice(0, value) for value in second.shape)] = second
    return float(np.max(abs(a - b)))


def basis_state(lengths, s, first, second, writer=1, bath=0):
    vector = np.zeros(36 * 2**sum(lengths), complex)
    vector[flat_core(s, first, second, writer) * 2**sum(lengths) + bath] = 1
    return vector


def physical_scores(vector):
    probabilities = (abs(vector.reshape(2, 3, 3, 2, -1))**2).sum(axis=(3, 4))
    return {
        "historical_zero_codeword": float(probabilities[0, 0, 0]),
        "historical_one_codeword": float(probabilities[1, 1, 1]),
        "historical_zero_source": float(probabilities[0].sum()),
        "both_blank": float(probabilities[:, 2, 2].sum()),
        "not_ready": float(1 - probabilities[:, :2, :2].sum()),
    }


def main():
    errors = {key: 0.0 for key in [
        "ready_coordinate_bijection", "physical_fault_maps", "charge_commutators",
        "logical_conservation_between_faults", "global_color_symmetry",
        "writer_subspace_invariance", "writer_matrix", "ready_embedding",
        "ready_joint_versus_gram", "ready_normalization", "collision_state",
        "collision_energy_distribution", "wrong_majority_historical_codeword",
        "triple_fault_repair_activity", "undamaged_repair_activity",
        "three_state_recurrence_matrix", "three_state_recurrence_return",
    ]}

    def update(key, value):
        errors[key] = max(errors[key], float(value))

    for code in range(8):
        bits = tuple((code >> j) & 1 for j in [2, 1, 0])
        logical, x1, x2 = encode(bits)
        update("ready_coordinate_bijection", max(abs(a-b) for a, b in zip(bits, decode(logical, x1, x2))))
        for site in MASKS:
            changed = list(bits)
            changed[{"S": 0, "M1": 1, "M2": 2}[site]] ^= 1
            actual = encode(changed)
            d = LOGICAL[site]
            expected = (logical ^ ((d[0]*x1 + d[1]*x2) % 2),
                        x1 ^ MASKS[site][0], x2 ^ MASKS[site][1])
            update("physical_fault_maps", max(abs(a-b) for a, b in zip(actual, expected)))

    scenarios = {
        "S_four_times": ((0.0, ("S",)), (0.83, ("S",)), (1.91, ("S",)), (3.07, ("S",))),
        "M1_four_times": ((0.0, ("M1",)), (0.83, ("M1",)), (1.91, ("M1",)), (3.07, ("M1",))),
        "alternating_all_sites": ((0.0, ("S",)), (0.83, ("M1",)), (1.91, ("M2",)), (3.07, ("S",))),
        "alternating_memories": ((0.0, ("M1",)), (0.83, ("M2",)), (1.91, ("M1",)), (3.07, ("M2",))),
    }
    terminal = 4.51
    physical_cases = []
    perturbation_tests = []
    fault_work = {}
    for lengths in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        h, h0, qw, charges, flips, logical_z = physical_hamiltonians(lengths)
        for charge in [qw, *charges]:
            update("charge_commutators", np.max(abs(h @ charge - charge @ h)))
        update("logical_conservation_between_faults", np.max(abs(h @ logical_z - logical_z @ h)))
        global_flip = flips["S"] @ flips["M1"] @ flips["M2"]
        update("global_color_symmetry", np.max(abs(h @ global_flip - global_flip @ h)))
        ready_h, _ = ready_hamiltonian(lengths)
        bath_dimension = 2**sum(lengths)
        ready_indices = np.concatenate([
            np.arange(flat_core(s, a, b, 1) * bath_dimension,
                      (flat_core(s, a, b, 1) + 1) * bath_dimension)
            for s in range(2) for a in range(2) for b in range(2)
        ])
        update("ready_embedding", np.max(abs(h0[np.ix_(ready_indices, ready_indices)]
            - ready_h - DELTA * np.eye(len(ready_indices)))))
        for r in [0, 1]:
            writer = np.column_stack([basis_state(lengths, r, 2, 2, writer=0),
                                      basis_state(lengths, r, r, r)])
            writer_h = np.array([[DELTA, KAPPA], [KAPPA, DELTA]])
            update("writer_matrix", np.max(abs(writer.conj().T @ h @ writer - writer_h)))
            update("writer_subspace_invariance", np.max(abs(h @ writer - writer @ writer_h)))

        actual = DenseEvolution(h)
        reference = DenseEvolution(h0)
        # The writer subspace above is exactly invariant. At the symbolic
        # nominal phase pi/2 its ready state is exact; avoid evolving an
        # approximately degenerate full matrix over the huge Rabi time.
        encoded = basis_state(lengths, 0, 0, 0)
        for name, pulses in scenarios.items():
            on = actual.sequence(encoded, pulses, terminal, flips)
            off = reference.sequence(encoded, pulses, terminal, flips)
            distance = float(np.linalg.norm(on - off))
            bound = KAPPA * terminal
            if distance > bound + 1e-10:
                raise AssertionError("Writer Duhamel bound violated")
            perturbation_tests.append({"lengths": list(lengths), "scenario": name,
                "vector_distance_writer_on_vs_omitted": distance,
                "kappa_elapsed_time_upper": bound})

        # Also test actual writing-window residuals, with exact Rabi input.
        phase = math.pi / 2 + 1 / 512
        initial_window = (math.cos(phase) * basis_state(lengths, 0, 2, 2, writer=0)
                          - 1j * math.sin(phase) * encoded)
        on = actual.sequence(initial_window, scenarios["alternating_all_sites"], terminal, flips)
        off = reference.sequence(initial_window, scenarios["alternating_all_sites"], terminal, flips)
        assert np.linalg.norm(on - off) <= KAPPA * terminal + 1e-10

        # Historical ambiguity: identical full present states from two
        # distinct preparations/fault histories, including W and both baths.
        one_fault = flips["S"] @ encoded
        two_fault = flips["M1"] @ flips["M2"] @ basis_state(lengths, 1, 1, 1)
        update("collision_state", np.max(abs(one_fault - two_fault)))
        update("collision_energy_distribution", np.max(abs(
            abs(actual.vectors.conj().T @ one_fault)**2
            - abs(actual.vectors.conj().T @ two_fault)**2)))
        for sites in [("S", "M1"), ("S", "M2"), ("M1", "M2"), ("S", "M1", "M2")]:
            faulted = encoded.copy()
            for site in sites:
                faulted = flips[site] @ faulted
            result = actual.step(faulted, 1.37)
            update("wrong_majority_historical_codeword", physical_scores(result)["historical_zero_codeword"])
            if len(sites) == 3:
                shape = result.reshape(36, bath_dimension)
                update("triple_fault_repair_activity", float(abs(shape[:, 1:]).sum()))
        undamaged = actual.step(encoded, 1.37)
        update("undamaged_repair_activity", float(abs(undamaged.reshape(36, bath_dimension)[:, 1:]).sum()))
        if lengths == (1, 1):
            energy_before = float(np.vdot(encoded, h @ encoded).real)
            for site in MASKS:
                state = flips[site] @ encoded
                fault_work[site] = float(np.vdot(state, h @ state).real) - energy_before
        physical_cases.append({"lengths": list(lengths), "full_active_dimension": len(h)})

    # Independent dense three-bit physical reference for the Gram/Fock API.
    lengths = (3, 3)
    h, flips = ready_hamiltonian(lengths)
    ready = DenseEvolution(h)
    initial = np.zeros(len(h), complex)
    initial[0] = 1
    api_cases = dict(scenarios)
    api_cases.update({
        "single_S": ((0.0, ("S",)),),
        "single_M1": ((0.0, ("M1",)),),
        "single_M2": ((0.0, ("M2",)),),
        "simultaneous_pair": ((0.0, ("M1", "M2")),),
        "simultaneous_triple": ((0.0, ("S", "M1", "M2")),),
    })
    api_results = []
    for name, pulses in api_cases.items():
        state = ready.sequence(initial, pulses, terminal, flips)
        direct = ready_joint(state, lengths)
        supplied = model.ready_statistics(lengths, pulses, terminal)["joint"]
        discrepancy = compare_padded(direct, supplied)
        update("ready_joint_versus_gram", discrepancy)
        update("ready_normalization", max(abs(direct.sum() - 1), abs(supplied.sum() - 1)))
        api_results.append({"scenario": name, "absolute_joint_error": discrepancy,
                            "oracle_shape": list(direct.shape), "api_shape": list(supplied.shape)})

    # Actual H with one-spin baths: exact three-state recurrence.
    lengths = (1, 1)
    h, _, _, _, flips, _ = physical_hamiltonians(lengths)
    evolution = DenseEvolution(h)
    a = flips["M1"] @ basis_state(lengths, 0, 0, 0)
    b_state = basis_state(lengths, 0, 0, 0, bath=2)
    c = basis_state(lengths, 0, 2, 2, writer=0, bath=2)
    embedding = np.column_stack([a, b_state, c])
    expected_h = 2*DELTA*np.eye(3) + np.array([[0, 1, 0], [1, 0, KAPPA], [0, KAPPA, 0]])
    update("three_state_recurrence_matrix", np.max(abs(h @ embedding - embedding @ expected_h)))
    omega = math.sqrt(1 + KAPPA**2)
    period = 2*math.pi/omega
    returned = evolution.step(a, period)
    return_probability = float(abs(np.vdot(a, returned))**2)
    update("three_state_recurrence_return", abs(return_probability - 1))
    loaded_time = math.pi/(2*omega)
    loaded = evolution.step(b_state, loaded_time)
    vacuum = evolution.step(basis_state(lengths, 0, 0, 0), loaded_time)
    loaded_fault = float(abs(np.vdot(a, loaded))**2)
    vacuum_fault = float(abs(np.vdot(a, vacuum))**2)
    assert loaded_fault > 0.999999 and vacuum_fault < 1e-10

    if max(errors.values()) >= 1e-10:
        raise AssertionError(f"Independent physical reference mismatch: {errors}")
    scaling = []
    for d in [2, 3, 4, 8, 16, 64]:
        histories = 1 + 3*(d-1)
        scaling.append({"alphabet_size": d, "three_copy_dimension": d**3,
                        "single_fault_and_no_fault_histories_per_logical_label": histories,
                        "minimum_pure_environment_dimension_exact": histories,
                        "minimum_environment_entropy_capacity_bits": math.log2(histories),
                        "two_copy_one_unflagged_fault_minimax_upper": 0.5,
                        "three_copy_two_fault_ambiguity_max_labels": min(d, 3),
                        "three_copy_two_fault_minimax_upper": 1/min(d, 3)})
    output = {
        "reviewer": "separate internal AI implementation; not external RAI or human peer review",
        "physical_cases": physical_cases,
        "ready_gram_oracle_dimension": 8 * 2**6,
        "ready_gram_comparisons": api_results,
        "worst_absolute_errors": errors,
        "all_independent_physical_checks_pass": True,
        "writer_perturbation_tests": perturbation_tests,
        "single_fault_mean_work_from_codeword": fault_work,
        "short_chain_recurrence": {"period": period, "return_to_fault_probability": return_probability},
        "loaded_bath_counterexample": {"read_time": loaded_time,
            "loaded_bath_fault_probability": loaded_fault, "vacuum_bath_fault_probability": vacuum_fault,
            "scope": "different energy and preparation; resource counterexample, not an energy-matched causal control"},
        "d_ary_bounds": scaling,
        "scope": "small dense reference validates the actual H and ready comparator; large uniform writer-on claim requires the separate analytic proof",
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
