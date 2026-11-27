"""Independent exact audit via the signed clock-orbit representation."""

import json
from fractions import Fraction
from itertools import product
from math import isqrt
from pathlib import Path


def prefix(data, record, time, pulses):
    for step in range(1, time + 1):
        pulse = pulses.get(str(step))
        if pulse == "X_D":
            data ^= 1
        elif pulse == "CNOT_D_R":
            record ^= data
        elif pulse is not None:
            raise ValueError(f"Unrecognized registered pulse {pulse}")
    return data, record


def histories(parameters, strides, pulses):
    """Index basis vectors by the time-zero sector, then dress each point."""
    size = parameters["clock_size"]
    sectors = product(
        range(len(strides)),
        range(size),
        range(parameters["data_dimension"]),
        range(parameters["record_dimension"]),
    )
    all_histories = {}
    for seed in sectors:
        matter, offset, data, record = seed
        path = []
        for a in range(size):
            b = (offset + strides[matter] * a) % size
            current_data, current_record = data, record
            phase = 1
            if offset == parameters["writer_relative_sector"]:
                current_data, current_record = prefix(data, record, a, pulses)
            if offset == parameters["phase_relative_sector"] and a == parameters["phase_tick"]:
                phase = parameters["phase_sign"]
            path.append(((a, b, matter, current_data, current_record), phase))
        all_histories[seed] = path
    return all_histories


def signed_cycle_data(paths, wrap_sign=1):
    """Read the complete signed transition map and count invariant cycles."""
    successors = {}
    for path in paths.values():
        for step, (point, sign) in enumerate(path):
            next_step = (step + 1) % len(path)
            next_point, next_sign = path[next_step]
            if point in successors:
                raise ValueError("Different history sectors overlap")
            coefficient = sign * next_sign
            if next_step == 0:
                coefficient *= wrap_sign
            successors[point] = next_point, coefficient

    remaining = set(successors)
    cycles = []
    while remaining:
        start = min(remaining)
        position = start
        length = 0
        holonomy = 1
        while position in remaining:
            remaining.remove(position)
            position, sign = successors[position]
            holonomy *= sign
            length += 1
        if position != start:
            raise ValueError("Shift is not a permutation of the basis")
        cycles.append((length, holonomy))
    return successors, cycles


def sections(paths, clock, reading):
    """Return the support of each logical column at the designated reading."""
    return {
        seed: [
            (point[:clock] + point[clock + 1 :], phase)
            for point, phase in path
            if point[clock] == reading
        ]
        for seed, path in paths.items()
    }


def is_unitary(section):
    rows = [row for entries in section.values() for row, _ in entries]
    return (
        all(len(entries) == 1 for entries in section.values())
        and len(set(rows)) == len(section)
    )


def vector_in_section(section, seed_vector):
    result = {}
    for label, amplitude in seed_vector.items():
        (point, phase), = section[label]
        result[point] = amplitude * phase
    return result


def transport_vector(vector, source_section, destination_section):
    """Apply the change of coordinates without re-preparing the probe."""
    transported = {}
    for label, entries in source_section.items():
        (source, source_phase), = entries
        if source not in vector:
            continue
        (destination, destination_phase), = destination_section[label]
        transported[destination] = source_phase * destination_phase * vector[source]
    return transported


def seed_state(parameters, sector, relative_phase=1):
    data = parameters["initial_data"]
    record = parameters["initial_record"]
    return {(0, sector, data, record): 1, (1, sector, data, record): relative_phase}


def noisy_response(size, lag, error, raw, reported):
    if raw == reported:
        return 1 - error
    if (raw + lag) % size == reported:
        return error
    return Fraction(0)


def event_probability(paths, section, seed_vector, outcome_weight):
    """Transform a compressed diagonal effect into a chart before evaluating it."""
    state = vector_in_section(section, seed_vector)
    squared_norm = sum(amplitude**2 for amplitude in state.values())
    diagonal_effect = {}
    for seed, path in paths.items():
        (row, _), = section[seed]
        diagonal_effect[row] = Fraction(
            sum(outcome_weight(point) for point, _ in path), len(path)
        )
    return sum(
        Fraction(amplitude**2, squared_norm) * diagonal_effect[row]
        for row, amplitude in state.items()
    )


def record_given_reading(paths, section, seed_vector, response):
    denominator = event_probability(paths, section, seed_vector, response)
    joint = event_probability(
        paths, section, seed_vector, lambda point: response(point) * point[-1]
    )
    return denominator, joint, joint / denominator


def projector_probability(state, probe):
    overlap = sum(value * probe.get(row, 0) for row, value in state.items())
    state_norm = sum(value**2 for value in state.values())
    probe_norm = sum(value**2 for value in probe.values())
    return Fraction(overlap**2, state_norm * probe_norm)


def reduced_matter(state_at_a):
    """Trace out the B register in an ideal A coordinate slice."""
    entries = {}
    norm = sum(value**2 for value in state_at_a.values())
    for left, a_left in state_at_a.items():
        for right, a_right in state_at_a.items():
            if left[0] == right[0]:
                key = left[1:], right[1:]
                entries[key] = entries.get(key, Fraction(0)) + Fraction(a_left * a_right, norm)
    return {key: value for key, value in entries.items() if value}


def noisy_matter(paths, size, seed_vector, response):
    entries = {}
    for a in range(size):
        state = vector_in_section(sections(paths, 0, a), seed_vector)
        for key, value in reduced_matter(state).items():
            # Raw readings are uniform in the good clock. The likelihood is
            # therefore also the conditional weight of a given one report.
            entries[key] = entries.get(key, Fraction(0)) + response(a) * value
    return {key: value for key, value in entries.items() if value}


def matter_visibility(density):
    return sum(
        value
        for (left, right), value in density.items()
        if left[0] != right[0] and left[1:] == right[1:]
    )


def equal_state_distance(first, second):
    if first != second:
        raise ValueError("The two reduced states are not identical")
    return Fraction(0)


def exact_square_root(number):
    number = Fraction(number)
    root = Fraction(isqrt(number.numerator), isqrt(number.denominator))
    if root**2 != number:
        raise ValueError("The benchmark distance is not rational")
    return root


def clock_phase_minor(parameters, successors):
    controls = parameters["clock_interaction_minor"]
    eigenstate = controls["D_R_eigenstate"]
    if eigenstate != "++":
        raise ValueError("This audit uses the registered uniform matter eigenstate")
    data_records = tuple(
        product(range(parameters["data_dimension"]), range(parameters["record_dimension"]))
    )
    signs = []
    for a, b in product(controls["A_rows"], controls["B_columns"]):
        images = [successors[(a, b, controls["S"], d, r)] for d, r in data_records]
        if {point[3:] for point, _ in images} != set(data_records):
            raise ValueError("The uniform matter vector did not stay invariant")
        phases = {phase for _, phase in images}
        if len(phases) != 1:
            raise ValueError("The matter vector acquired a nonuniform phase")
        signs.append(phases.pop())
    return signs[0] * signs[3] - signs[1] * signs[2]


def calculate_observables(parameters):
    size = parameters["clock_size"]
    strides = parameters["clock_strides"]
    pulses = parameters["pulse_entries"]
    paths = histories(parameters, strides, pulses)
    successors, cycles = signed_cycle_data(paths)

    anchor_a = parameters["reference_chart_anchors"]["A"]
    anchor_b = parameters["reference_chart_anchors"]["B"]
    a_chart, b_chart = sections(paths, 0, anchor_a), sections(paths, 1, anchor_b)
    selected_a, selected_b = (parameters["selected_readings"][clock] for clock in ("A", "B"))
    seed = seed_state(parameters, parameters["main_relative_sector"])
    opposite_seed = seed_state(parameters, parameters["main_relative_sector"], -1)
    raw_a = lambda point: int(point[0] == selected_a)
    raw_b = lambda point: int(point[1] == selected_b)
    error_a, error_b = (Fraction(parameters["error_probabilities"][clock]) for clock in ("A", "B"))
    lag = parameters["reported_lag"]
    noisy_a = lambda point: noisy_response(size, lag, error_a, point[0], selected_a)
    noisy_b = lambda point: noisy_response(size, lag, error_b, point[1], selected_b)
    joint_response = lambda point: noisy_a(point) * noisy_b(point)

    joint_a = record_given_reading(paths, a_chart, seed, joint_response)
    joint_b = record_given_reading(paths, b_chart, seed, joint_response)
    no_write_pulses = {step: gate for step, gate in pulses.items() if gate != "CNOT_D_R"}
    no_write_paths = histories(parameters, strides, no_write_pulses)
    no_write_chart = sections(no_write_paths, 0, anchor_a)
    no_write = record_given_reading(no_write_paths, no_write_chart, seed, joint_response)

    a_sections = [sections(paths, 0, a) for a in range(size)]
    b_sections = [sections(paths, 1, b) for b in range(size)]
    ideal_a_record = []
    mixed_a_record = []
    for section in a_sections:
        state = vector_in_section(section, seed)
        ideal_a_record.append(
            str(sum(Fraction(amplitude**2, len(seed)) * row[-1] for row, amplitude in state.items()))
        )
        mixed_records = sum(entries[0][0][-1] for entries in section.values())
        mixed_a_record.append(str(Fraction(mixed_records, len(section))))

    selected_section = a_sections[selected_a]
    selected_state = vector_in_section(selected_section, seed)
    opposite_state = vector_in_section(selected_section, opposite_seed)
    both_clock_distance = exact_square_root(
        1 - projector_probability(selected_state, opposite_state)
    )
    reduced = reduced_matter(selected_state)
    other_reduced = reduced_matter(opposite_state)
    noisy_rows = noisy_matter(
        paths, size, seed, lambda a: noisy_response(size, lag, error_a, a, selected_a)
    )
    noisy_other_rows = noisy_matter(
        paths, size, opposite_seed, lambda a: noisy_response(size, lag, error_a, a, selected_a)
    )

    phase_seed = seed_state(parameters, parameters["phase_relative_sector"])
    phase_a = vector_in_section(a_chart, phase_seed)
    state_b = vector_in_section(b_chart, phase_seed)
    probe_b = transport_vector(phase_a, a_chart, b_chart)
    phase_b_wrong = {row: abs(amplitude) for row, amplitude in probe_b.items()}

    stopped = histories(parameters, parameters["controls"]["zero_stride"], pulses)
    stopped_section = sections(stopped, 1, anchor_b)
    _, twisted_cycles = signed_cycle_data(paths, parameters["controls"]["wrap_sign"])
    visibility_reference = vector_in_section(
        a_sections[parameters["calibration_visibility_tick"]], seed
    )

    return {
        "kinematic_dimension": len(successors),
        "physical_dimension": sum(sign == 1 for _, sign in cycles),
        "cyclic_power_identity": all(length == size and sign == 1 for length, sign in cycles),
        "A_chart_unitary": all(is_unitary(section) for section in a_sections),
        "B_chart_unitary": all(is_unitary(section) for section in b_sections),
        "clock_interaction_minor": clock_phase_minor(parameters, successors),
        "A_zero_visibility": str(matter_visibility(reduced_matter(visibility_reference))),
        "A_selected_visibility_after_B_loss": str(matter_visibility(reduced)),
        "B_readings_at_A_selected": sorted({row[0] for row in selected_state}),
        "record_by_ideal_A": ideal_a_record,
        "ideal_selected_A_record": str(record_given_reading(paths, a_chart, seed, raw_a)[2]),
        "ideal_selected_B_record": str(record_given_reading(paths, b_chart, seed, raw_b)[2]),
        "noisy_selected_A_record": str(record_given_reading(paths, a_chart, seed, noisy_a)[2]),
        "noisy_selected_B_record": str(record_given_reading(paths, b_chart, seed, noisy_b)[2]),
        "matched_joint_filter_A_chart": str(joint_a[0]),
        "matched_joint_filter_B_chart": str(joint_b[0]),
        "matched_joint_record_A_chart": str(joint_a[1]),
        "matched_joint_record_B_chart": str(joint_b[1]),
        "matched_conditional_A_chart": str(joint_a[2]),
        "matched_conditional_B_chart": str(joint_b[2]),
        "coherent_phase_probe_A_chart": str(projector_probability(phase_a, phase_a)),
        "coherent_phase_probe_B_chart": str(projector_probability(state_b, probe_b)),
        "coherent_phase_probe_B_wrong_sign": str(projector_probability(state_b, phase_b_wrong)),
        "no_write_joint_filter": str(no_write[0]),
        "no_write_matched_conditional": str(no_write[2]),
        "mixed_record_by_ideal_A": mixed_a_record,
        "orthogonal_seed_full_distinguishability": str(both_clock_distance),
        "missing_B_ideal_distinguishability": str(equal_state_distance(reduced, other_reduced)),
        "missing_B_noisy_distinguishability": str(equal_state_distance(noisy_rows, noisy_other_rows)),
        "zero_stride_B_chart_rank": sum(bool(entries) for entries in stopped_section.values()),
        "sign_flipped_wrap_invariant_dimension": sum(sign == 1 for _, sign in twisted_cycles),
    }


if __name__ == "__main__":
    data = json.loads(Path(__file__).with_name("PARAMETERS.json").read_text(encoding="utf-8"))
    print(json.dumps({"observables": calculate_observables(data)}, indent=2, sort_keys=True))
