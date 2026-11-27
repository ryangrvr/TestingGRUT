"""Orbit calculation for the frozen GRUT-F5 two-clock fixture."""
import json
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent


class Cycle:
    def __init__(self, params, *, write=True, strides=None, wrap_sign=1):
        self.n = params["clock_size"]
        self.strides = params["clock_strides"] if strides is None else strides
        self.writer_relative_sector = params["writer_relative_sector"]
        self.phase_relative_sector = params["phase_relative_sector"]
        self.phase_tick = params["phase_tick"]
        self.phase_sign = params["phase_sign"]
        self.wrap_sign = wrap_sign
        self.gates = sorted(
            (int(tick), gate)
            for tick, gate in params["pulse_entries"].items()
            if write or gate != "CNOT_D_R"
        )
        self.logical_basis = list(product(
            range(len(self.strides)), range(self.n),
            range(params["data_dimension"]), range(params["record_dimension"]),
        ))
        self.basis = [
            (a, b, s, data, record)
            for s, b, data, record in self.logical_basis
            for a in range(self.n)
        ]

    def prefix(self, a, rel, data, record, *, inverse):
        gates = reversed(self.gates) if inverse else self.gates
        if rel == self.writer_relative_sector:
            for tick, gate in gates:
                if tick <= a:
                    if gate == "X_D":
                        data ^= 1
                    elif gate == "CNOT_D_R":
                        record ^= data
                    else:
                        raise ValueError(f"Unimplemented gate {gate}")
        phase = self.phase_sign if (rel, a) == (self.phase_relative_sector, self.phase_tick) else 1
        return data, record, phase

    def step(self, row):
        a, b, s, data, record = row
        rel = (b - self.strides[s] * a) % self.n
        next_a = (a + 1) % self.n
        next_b = (b + self.strides[s]) % self.n
        data, record, before = self.prefix(a, rel, data, record, inverse=True)
        data, record, after = self.prefix(next_a, rel, data, record, inverse=False)
        sign = before * after * (self.wrap_sign if next_a == 0 else 1)
        return (next_a, next_b, s, data, record), sign

    def orbit(self, start):
        row, sign = start, 1
        path = []
        for _ in range(self.n):
            path.append((row, sign))
            row, factor = self.step(row)
            sign *= factor
        return path, row, sign

    def columns(self):
        return {
            key: self.orbit((0, key[1], key[0], key[2], key[3]))[0]
            for key in self.logical_basis
        }

    def cycle_residuals(self):
        return [
            (start, end, sign)
            for start in self.basis
            for _, end, sign in [self.orbit(start)]
            if end != start or sign != 1
        ]


def chart(columns, axis, tick):
    return {
        key: {
            row[:axis] + row[axis + 1:]: sign
            for row, sign in orbit
            if row[axis] == tick
        }
        for key, orbit in columns.items()
    }


def is_unitary(mapping, dimension):
    entries = [next(iter(image)) for image in mapping.values() if len(image) == 1]
    return len(entries) == len(set(entries)) == dimension


def move_state(mapping, source):
    target = defaultdict(int)
    for key, amplitude in source.items():
        for row, sign in mapping[key].items():
            target[row] += sign * amplitude
    return dict(target)


def pull_back_projector_vector(mapping, probe):
    return {
        key: phase * probe[row]
        for key, image in mapping.items()
        for row, phase in image.items()
        if row in probe
    }


def read_probability(underlying, reported, n, lag, error):
    if reported == underlying:
        return 1 - error
    if reported == (underlying + lag) % n:
        return error
    return Fraction(0)


def clock_effect(axis, reported, n, lag, error):
    return lambda row: read_probability(row[axis], reported, n, lag, error)


def record_effect(row):
    return row[4]


def compress_diagonal(columns, effect):
    return {
        key: Fraction(sum(effect(row) for row, _ in orbit), len(orbit))
        for key, orbit in columns.items()
    }


def mean(diagonal, state):
    norm2 = sum(amplitude**2 for amplitude in state.values())
    return sum(diagonal[k] * amplitude**2 for k, amplitude in state.items()) / norm2


def record_probability(columns, state, clock):
    record = compress_diagonal(columns, lambda row: record_effect(row) * clock(row))
    clock_filter = compress_diagonal(columns, clock)
    return mean(record, state) / mean(clock_filter, state)


def in_chart(diagonal, mapping):
    return {next(iter(mapping[key])): value for key, value in diagonal.items()}


def overlap_squared(left, right):
    dot = sum(amplitude * right.get(row, 0) for row, amplitude in left.items())
    return Fraction(dot**2, sum(v**2 for v in left.values()) * sum(v**2 for v in right.values()))


def visibility(vector):
    return Fraction(
        sum(
            amp * vector.get((b, 1 - s, data, record), 0)
            for (b, s, data, record), amp in vector.items()
        ),
        sum(amp**2 for amp in vector.values()),
    )


def discard_b(vector):
    norm2 = sum(amp**2 for amp in vector.values())
    density = defaultdict(Fraction)
    for left, amp in vector.items():
        for right, other in vector.items():
            if left[0] == right[0]:
                density[(left[1:], right[1:])] += Fraction(amp * other, norm2)
    return dict(density)


def noisy_a_without_b(columns, reported, n, lag, error, source):
    density = defaultdict(Fraction)
    for tick in (reported, (reported - lag) % n):
        weight = read_probability(tick, reported, n, lag, error)
        state_at_tick = move_state(chart(columns, 0, tick), source)
        for element, value in discard_b(state_at_tick).items():
            density[element] += weight * value
    return dict(density)


def clock_phase_minor(model, fixture):
    s = fixture["S"]
    assert fixture["D_R_eigenstate"] == "++"
    a0, a1 = fixture["A_rows"]
    b0, b1 = fixture["B_columns"]
    # X and CNOT both fix the two-qubit |++> state, so only the phase remains.
    _, w00 = model.step((a0, b0, s, 0, 0))
    _, w01 = model.step((a0, b1, s, 0, 0))
    _, w10 = model.step((a1, b0, s, 0, 0))
    _, w11 = model.step((a1, b1, s, 0, 0))
    return w00 * w11 - w01 * w10


def calculate(p):
    model = Cycle(p)
    columns = model.columns()
    n = model.n
    dim = len(columns)
    a, b = (p["selected_readings"][clock] for clock in ("A", "B"))
    a_anchor, b_anchor = (p["reference_chart_anchors"][clock] for clock in ("A", "B"))
    error_a, error_b = (Fraction(p["error_probabilities"][clock]) for clock in ("A", "B"))
    lag = p["reported_lag"]
    zero_error = Fraction(0)
    data, record = p["initial_data"], p["initial_record"]
    rel, phase_rel = p["main_relative_sector"], p["phase_relative_sector"]
    plus = {(0, rel, data, record): 1, (1, rel, data, record): 1}
    minus = {(0, rel, data, record): 1, (1, rel, data, record): -1}
    phase_seed = {(0, phase_rel, data, record): 1, (1, phase_rel, data, record): 1}
    a_chart, b_chart = chart(columns, 0, a_anchor), chart(columns, 1, b_anchor)
    a_slice = chart(columns, 0, a)
    a_state, b_state = move_state(a_chart, plus), move_state(b_chart, plus)
    a_phase, b_phase = move_state(a_chart, phase_seed), move_state(b_chart, phase_seed)
    # The declared phase probe has two even-phase components at A=2, r=1.
    # Pull that A-defined effect to the seed space before viewing it in B.
    a_probe = {
        ((phase_rel + stride * a_anchor) % n, s, data, record): 1
        for s, stride in enumerate(model.strides)
    }
    b_probe = move_state(b_chart, pull_back_projector_vector(a_chart, a_probe))
    a_noisy = clock_effect(0, a, n, lag, error_a)
    b_noisy = clock_effect(1, b, n, lag, error_b)
    joint_filter = lambda row: a_noisy(row) * b_noisy(row)
    joint_record = lambda row: joint_filter(row) * record_effect(row)
    joint_diagonal = compress_diagonal(columns, joint_filter)
    record_diagonal = compress_diagonal(columns, joint_record)
    filter_a = mean(in_chart(joint_diagonal, a_chart), a_state)
    filter_b = mean(in_chart(joint_diagonal, b_chart), b_state)
    record_a = mean(in_chart(record_diagonal, a_chart), a_state)
    record_b = mean(in_chart(record_diagonal, b_chart), b_state)
    a_without_b = discard_b(move_state(a_slice, plus))
    a_minus_without_b = discard_b(move_state(a_slice, minus))
    noisy_plus = noisy_a_without_b(columns, a, n, lag, error_a, plus)
    noisy_minus = noisy_a_without_b(columns, a, n, lag, error_a, minus)
    stopped = Cycle(p, strides=p["controls"]["zero_stride"])
    stopped_images = chart(stopped.columns(), 1, b_anchor).values()
    twisted = Cycle(p, wrap_sign=p["controls"]["wrap_sign"])
    no_writer = Cycle(p, write=False)
    no_writer_columns = no_writer.columns()
    no_writer_filter = mean(compress_diagonal(no_writer_columns, joint_filter), plus)
    no_writer_record = mean(compress_diagonal(no_writer_columns, joint_record), plus)

    history = []
    mixed_history = []
    for tick in range(n):
        perfect_tick = clock_effect(0, tick, n, lag, zero_error)
        history.append(str(record_probability(columns, plus, perfect_tick)))
        all_seeds = compress_diagonal(columns, lambda row: perfect_tick(row) * record_effect(row))
        mixed_history.append(str(sum(all_seeds.values()) * n / dim))

    zero_slice = chart(columns, 0, p["calibration_visibility_tick"])
    observables = {
        "kinematic_dimension": len(model.basis),
        "physical_dimension": sum(
            model.orbit((0, r, s, d, q))[2] == 1 for s, r, d, q in columns
        ),
        "cyclic_power_identity": not model.cycle_residuals(),
        "A_chart_unitary": all(is_unitary(chart(columns, 0, tick), dim) for tick in range(n)),
        "B_chart_unitary": all(is_unitary(chart(columns, 1, tick), dim) for tick in range(n)),
        "clock_interaction_minor": clock_phase_minor(model, p["clock_interaction_minor"]),
        "A_zero_visibility": str(visibility(move_state(zero_slice, plus))),
        "A_selected_visibility_after_B_loss": str(visibility(move_state(a_slice, plus))),
        "B_readings_at_A_selected": sorted({row[0] for row in move_state(a_slice, plus)}),
        "record_by_ideal_A": history,
        "ideal_selected_A_record": str(Fraction(history[a])),
        "ideal_selected_B_record": str(record_probability(
            columns, plus, clock_effect(1, b, n, lag, zero_error)
        )),
        "noisy_selected_A_record": str(record_probability(columns, plus, a_noisy)),
        "noisy_selected_B_record": str(record_probability(columns, plus, b_noisy)),
        "matched_joint_filter_A_chart": str(filter_a),
        "matched_joint_filter_B_chart": str(filter_b),
        "matched_joint_record_A_chart": str(record_a),
        "matched_joint_record_B_chart": str(record_b),
        "matched_conditional_A_chart": str(record_a / filter_a),
        "matched_conditional_B_chart": str(record_b / filter_b),
        "coherent_phase_probe_A_chart": str(overlap_squared(a_phase, a_probe)),
        "coherent_phase_probe_B_chart": str(overlap_squared(b_phase, b_probe)),
        "coherent_phase_probe_B_wrong_sign": str(
            overlap_squared(b_phase, {row: abs(amp) for row, amp in b_probe.items()})
        ),
        "no_write_joint_filter": str(no_writer_filter),
        "no_write_matched_conditional": str(no_writer_record / no_writer_filter),
        "mixed_record_by_ideal_A": mixed_history,
        "orthogonal_seed_full_distinguishability": "1" if overlap_squared(plus, minus) == 0 else None,
        "missing_B_ideal_distinguishability": "0" if a_without_b == a_minus_without_b else None,
        "missing_B_noisy_distinguishability": "0" if noisy_plus == noisy_minus else None,
        "zero_stride_B_chart_rank": sum(bool(image) for image in stopped_images),
        "sign_flipped_wrap_invariant_dimension": sum(
            twisted.orbit((0, r, s, d, q))[2] == 1 for s, r, d, q in columns
        ),
    }
    return observables


def main():
    p = json.loads((HERE / "PARAMETERS.json").read_text())
    observables = calculate(p)
    assert observables.keys() == p["expected"].keys()
    checks = {key: actual == p["expected"][key] for key, actual in observables.items()}
    result = {
        "observables": observables,
        "checks": checks,
        "passed": sum(checks.values()),
        "total": len(checks),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if all(checks.values()) else 1)


if __name__ == "__main__":
    main()
