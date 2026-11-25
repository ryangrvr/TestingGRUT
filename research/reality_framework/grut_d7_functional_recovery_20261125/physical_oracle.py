#!/usr/bin/env python3
"""Build the D7 physical transition matrix and propagate without the spin reduction."""

import argparse
from collections import defaultdict, deque
from dataclasses import dataclass
import itertools
import json
import math
from pathlib import Path
from typing import NamedTuple

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply


HERE = Path(__file__).resolve().parent
ROOT_TWO = math.sqrt(2)


class Data(NamedTuple):
    probe: int
    record: int
    bypass: int
    payload: int
    dump: int
    sink: int
    latch: int
    input: int
    receiver: int
    reference: int


@dataclass(frozen=True)
class Trial:
    length: int
    assay: str
    severity: str
    port: str
    intervention: str
    policy: str
    source: int
    passive_coin: int
    mask: int

    @property
    def theta(self):
        return {"none": 0.0, "weak": math.pi / 12, "strong": math.pi / 4}[self.assay]

    @property
    def phi(self):
        return {"shallow": math.pi / 8, "deep": math.pi / 4}[self.severity]

    @property
    def beta(self):
        if self.policy == "blind" or (self.policy == "calibrated" and self.assay == "none"):
            return 0.0
        if (self.policy, self.assay, self.severity) == ("calibrated", "weak", "shallow"):
            return math.atan(0.5) / 2
        return self.phi


def audit_geometry(length):
    sites = [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (2, 1),
             (1, 1), (1, 2), (2, 2)] + [(j - 6, 2) for j in range(9, length + 1)]
    assert len(set(sites)) == length + 1
    extra = {3: 5, 5: 1, 6: 5}
    for j in range(length):
        support = [sites[j], sites[j + 1]]
        assert sum(abs(y - x) for x, y in zip(*support, strict=True)) == 1
        if j in extra:
            support.append(sites[extra[j]])
        assert max(sum(abs(y - x) for x, y in zip(a, b, strict=True))
                   for a, b in itertools.combinations(support, 2)) <= 2


def transition(trial, link, data, reverse=False):
    if link == 0 or link == 5 or link == 6:
        if link == 0:
            phase = -trial.theta * (-1) ** (trial.source ^ data.probe)
        elif link == 5:
            phase = -trial.phi * (-1) ** (trial.source ^ data.payload)
        else:
            phase = trial.beta * (-1) ** (data.latch ^ data.payload)
        return [(data, np.exp(1j * (-phase if reverse else phase)))]

    if link == 1:
        extra = {"main": 0, "passive": trial.passive_coin, "direct": trial.source}[trial.port]

        def permute(v, back):
            if back:
                return v._replace(probe=v.dump, dump=v.probe, record=v.record ^ v.dump,
                                  bypass=v.bypass ^ extra)
            return v._replace(probe=v.dump, dump=v.probe, record=v.record ^ v.probe,
                              bypass=v.bypass ^ extra)

        if reverse:
            unpermuted = permute(data, True)
            return [(unpermuted, 1 / ROOT_TWO),
                    (unpermuted._replace(probe=1 ^ unpermuted.probe), 1j / ROOT_TWO)]
        return [(permute(data, False), 1 / ROOT_TWO),
                (permute(data._replace(probe=1 ^ data.probe), False), -1j / ROOT_TWO)]

    if link == 2:
        if trial.intervention == "scrambled":
            data = data._replace(record=data.record ^ trial.mask)
        elif trial.intervention == "inverted":
            data = data._replace(record=data.record ^ 1)
        elif trial.intervention == "disabled":
            data = data._replace(record=data.sink, sink=data.record)
    elif link == 3:
        select = data.record if trial.port == "main" else data.bypass
        data = data._replace(latch=data.latch ^ select)
    elif link == 4:
        data = data._replace(input=data.payload, payload=data.input)
    elif link == 7:
        data = data._replace(receiver=data.payload, payload=data.receiver)
    return [(data, 1.0)]


def reachable_hamiltonian(trial):
    roots = [(0, Data(probe, 0, 0, 0, 0, 0, 0, bit, 0, bit))
             for probe, bit in itertools.product((0, 1), repeat=2)]
    indices = {state: index for index, state in enumerate(roots)}
    queue = deque(roots)
    inverse_error = 0.0
    while queue:
        j, data = queue.popleft()
        neighbors = []
        if j < trial.length:
            forward = transition(trial, j, data)
            recovered = defaultdict(complex)
            for v, weight in forward:
                neighbors.append((j + 1, v))
                for old, undo in transition(trial, j, v, reverse=True):
                    recovered[old] += undo * weight
            recovered[data] -= 1
            inverse_error = max(inverse_error, max(map(abs, recovered.values())))
        if j > 0:
            neighbors.extend((j - 1, v) for v, _ in transition(trial, j - 1, data, reverse=True))
        for state in neighbors:
            if state not in indices:
                indices[state] = len(indices)
                queue.append(state)

    data_row, data_col, values = [], [], []
    for (j, registers), index in indices.items():
        if j >= trial.length:
            continue
        exchange = 2 * math.sqrt((j + 1) * (trial.length - j)) / (trial.length + 1)
        for out, entry in transition(trial, j, registers):
            target = indices[(j + 1, out)]
            value = exchange * entry
            data_row.extend((target, index))
            data_col.extend((index, target))
            values.extend((value, value.conjugate()))
    size = len(indices)
    matrix = coo_matrix((values, (data_row, data_col)), shape=(size, size), dtype=complex).tocsr()
    initial = np.zeros(size, dtype=complex)
    for root in roots:
        initial[indices[root]] = 0.5
    return tuple(indices), matrix, initial, inverse_error


def measure(vector, basis, source):
    bell_branches = defaultdict(complex)
    head_ready, dump_correct, latch_correct = 0.0, 0.0, 0.0
    for amplitude, (j, data) in zip(vector, basis, strict=True):
        probability = float(abs(amplitude) ** 2)
        head_ready += probability * (j >= 8)
        dump_correct += probability * (data.dump == source)
        latch_correct += probability * (data.latch == source)
        if data.receiver == data.reference:
            bell_branches[j, data[:-2]] += amplitude / ROOT_TWO
    fidelity = math.fsum(float(abs(a) ** 2) for a in bell_branches.values())
    return fidelity, head_ready, dump_correct, latch_correct


def energy_residual(hamiltonian, initial, length):
    eigenvalues, eigenvectors = eigh(hamiltonian.toarray())
    weights = abs(eigenvectors.conj().T @ initial) ** 2
    omega = 4 / (length + 1)
    indices = np.rint(eigenvalues / omega + length / 2).astype(int)
    if not np.all((0 <= indices) & (indices <= length)):
        raise AssertionError("Energy outside the registered spectrum")
    measured = np.bincount(indices, weights, minlength=length + 1)
    exact = [math.comb(length, k) / 2 ** length for k in range(length + 1)]
    errors = max(np.max(abs(eigenvalues - omega * (indices - length / 2))),
                 np.max(abs(measured - exact)))
    return float(errors)


def calculate(spec):
    rows = []
    worst = dict(inverse=0.0, energy=0.0, normalization=0.0, recurrence=0.0)
    dimensions = set()
    oracle = spec["oracle"]
    time_offsets = spec["timing"]["window_offsets_from_pi"]
    for length in oracle["length"]:
        audit_geometry(length)
    for length, source, severity, setting in itertools.product(
        oracle["length"], oracle["source"], oracle["severity"], oracle["rows"]
    ):
        assay, port, intervention, policy = setting
        accum = defaultdict(float)
        exact_time = math.pi * (length + 1) / 4
        times = dict(endpoint=exact_time, recurrence=2 * exact_time)
        for side, offset in zip(("left", "right"), time_offsets, strict=True):
            times[side] = (math.pi + offset) * (length + 1) / 4
        for cp, cm in itertools.product((0, 1), repeat=2):
            trial = Trial(length, assay, severity, port, intervention, policy, source, cp, cm)
            basis, hamiltonian, initial, inverse = reachable_hamiltonian(trial)
            dimensions.add(len(basis))
            worst["inverse"] = max(worst["inverse"], inverse)
            worst["energy"] = max(worst["energy"], energy_residual(hamiltonian, initial, length))
            for label, time in times.items():
                vector = expm_multiply(-1j * time * hamiltonian, initial, traceA=0.0)
                residual = abs(float(np.vdot(vector, vector).real) - 1)
                worst["normalization"] = max(worst["normalization"], residual)
                bell, ready, dump, latch = measure(vector, basis, source)
                accum[label + "_bell"] += bell / 4
                accum[label + "_ready"] += ready / 4
                if label == "endpoint":
                    accum["raw_record_correct"] += dump / 4
                    accum["latch_correct"] += latch / 4
                if label == "recurrence":
                    overlap = float(abs(np.vdot(initial, vector)) ** 2)
                    worst["recurrence"] = max(worst["recurrence"], abs(1.0 - overlap))
        rows.append(dict(length=length, assay=assay, severity=severity, port=port,
                         intervention=intervention, policy=policy, source=source, **accum))
    if max(worst.values()) > spec["gates"]["physical_oracle_residual_max"]:
        raise AssertionError(worst)
    return dict(rows=rows, coin_blocks=4 * len(rows), reachable_dimensions=sorted(dimensions),
                internal_residuals=worst)


KEY = ("length", "assay", "severity", "port", "intervention", "policy", "source")


def cross_check(result, reduced_path, spec):
    reduced = json.loads(reduced_path.read_text())
    lookup = {tuple(r[:len(KEY)]): dict(zip(reduced["columns"], r, strict=True)) for r in reduced["rows"]}
    row_keys = {tuple(row[k] for k in KEY) for row in result["rows"]}
    oracle = spec["oracle"]
    expected = {(length, assay, severity, port, action, policy, source)
                for length, severity, (assay, port, action, policy), source in itertools.product(
                    oracle["length"], oracle["severity"], oracle["rows"], oracle["source"])}
    if len(result["rows"]) != len(expected) or row_keys != expected:
        raise AssertionError("Missing or duplicate rows")
    errors = []
    for r in result["rows"]:
        forecast = lookup[tuple(r[k] for k in KEY)]
        clocks = reduced["clock"][str(r["length"])]
        for field in ("endpoint_bell", "recurrence_bell", "raw_record_correct", "latch_correct"):
            errors.append(abs(r[field] - forecast[field]))
        for side in ("left", "right"):
            errors.append(abs(r[side + "_bell"] - forecast[side + "_receiver_bell"]))
            errors.append(abs(r[side + "_ready"] - clocks["ready_" + side]))
        errors.extend([abs(r["endpoint_ready"] - 1), abs(r["recurrence_ready"])])
    if max(errors) > spec["gates"]["physical_oracle_residual_max"]:
        raise AssertionError(max(errors))
    return dict(scalars=len(errors), max_abs_error=max(errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--output", type=Path)
    output.add_argument("--verify", type=Path)
    args = parser.parse_args()
    parser_path = HERE / "PROTOCOL.json"
    spec = json.loads(parser_path.read_text())
    data = calculate(spec)
    data["cross_check"] = cross_check(data, HERE / "RESULTS.json", spec)
    if args.output:
        args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    elif json.loads(args.verify.read_text()) != data:
        raise SystemExit("Stored oracle differs from replay")
    print(f"{len(data['rows'])} physical rows in {data['coin_blocks']} coin blocks")
    print(dict(internal=data["internal_residuals"], cross_check=data["cross_check"]))


if __name__ == "__main__":
    main()
