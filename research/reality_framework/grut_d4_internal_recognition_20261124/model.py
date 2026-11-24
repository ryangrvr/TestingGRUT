"""Exact ready-sector dynamics and uniform bounds for the physical D4 model."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import isqrt, log10

import numpy as np


# Order is (comparison with M1, comparison with M2).
MASKS = {"S": (1, 1), "M1": (1, 0), "M2": (0, 1)}
MAJORITY_CHANGES = {"S": (1, 1), "M1": (0, 1), "M2": (1, 0)}
DELTA = 4.0
KAPPA = 1 / 65536


@lru_cache(maxsize=8)
def _creation_matrices(rank: int) -> np.ndarray:
    dimension = 1 << rank
    matrices = np.zeros((rank, dimension, dimension), dtype=complex)
    for mode in range(rank):
        bit = 1 << mode
        for state in range(dimension):
            if state & bit == 0:
                matrices[mode, state | bit, state] = (-1) ** (state & (bit - 1)).bit_count()
    return matrices


@lru_cache(maxsize=32)
def _endpoint_operators(length: int, times: tuple[float, ...]) -> tuple[dict, np.ndarray, dict]:
    """Endpoint observables in the Fock space of the exact time-mode span.

    Columns in the one-particle sine basis represent e^(+i h t)|endpoint>.
    Factor their overlap matrix via SVD, retaining its numerical rank. In the
    registered large calculation the distinct time columns are independent.
    """
    angles = np.arange(1, length + 2) * np.pi / (length + 2)
    mode_energy = DELTA + 2 * np.cos(angles)
    endpoint_weights = np.sqrt(2 / (length + 2)) * np.sin(angles)
    modes = endpoint_weights[:, None] * np.exp(1j * np.outer(mode_energy, times))
    _, singular, right = np.linalg.svd(modes, full_matrices=False)
    keep = singular > 5e-14
    coordinates = singular[keep, None] * right[keep, :]
    rank = int(keep.sum())
    creators = _creation_matrices(rank)
    dimension = 1 << rank
    identity = np.eye(dimension)
    operators = {}
    for index, time in enumerate(times):
        creation = np.tensordot(coordinates[:, index], creators, axes=1)
        annihilation = creation.conj().T
        occupation = creation @ annihilation
        operators[time] = {
            "X": creation + annihilation,
            "Z": identity - 2 * occupation,
            "occupation": occupation,
        }
    number = np.array([state.bit_count() for state in range(dimension)])
    reconstruction = float(np.max(abs(coordinates.conj().T @ coordinates - modes.conj().T @ modes)))
    return operators, number, {"rank": rank, "largest_mode_gram_error": reconstruction}


def _one_chain(length: int, channel: int, pulses, terminal_time: float):
    """Propagate both eigenvalues of the logical flip in the interaction picture."""
    times = tuple(sorted({terminal_time, *(time for time, _ in pulses)}))
    operators, number, diagnostics = _endpoint_operators(length, times)
    plus = np.zeros(len(number), dtype=complex)
    minus = np.zeros_like(plus)
    plus[0] = minus[0] = 1
    for time, sites in pulses:
        local = operators[time]
        for site in sites:
            if MAJORITY_CHANGES[site][channel]:
                minus = local["Z"] @ minus
            if MASKS[site][channel]:
                plus = local["X"] @ plus
                minus = local["X"] @ minus

    source = operators[terminal_time]["occupation"]
    count_limit = min(length, len(pulses))
    components = np.zeros((2, 2, count_limit + 1, len(number)), dtype=complex)
    for sign, state in enumerate((plus, minus)):
        for parity in (0, 1):
            for bath_count in range(count_limit + 1):
                component = np.where(number == bath_count + parity, state, 0)
                occupied = source @ component
                components[sign, parity, bath_count] = occupied if parity else component - occupied
    return components, diagnostics


def ready_statistics(
    lengths: tuple[int, int],
    pulses: tuple[tuple[float, tuple[str, ...]], ...],
    terminal_time: float,
) -> dict:
    """Return joint probabilities [relative majority, x1, x2, N_B1, N_B2].

    This exact reduction starts after writing in the ready sector. The writer
    exchange is omitted in this comparison only. Faults use the physical site
    maps, including logical changes on residual discrepancies. The coherent
    recombination below, rather than an incoherent average, reads history.
    """
    times = [time for time, _ in pulses]
    if min(lengths) < 1 or not np.isfinite(terminal_time) or terminal_time < 0:
        raise ValueError("Positive lengths and finite nonnegative read time required")
    if times != sorted(times) or any(not np.isfinite(t) or t < 0 or t > terminal_time for t in times):
        raise ValueError("Fault times must be ordered and at or before readout")
    if any(not sites or any(site not in MASKS for site in sites) for _, sites in pulses):
        raise ValueError("Each pulse must specify at least one physical site")

    first, first_diagnostics = _one_chain(lengths[0], 0, pulses, terminal_time)
    second, second_diagnostics = _one_chain(lengths[1], 1, pulses, terminal_time)
    joint = np.zeros((2, 2, 2, first.shape[2], second.shape[2]))
    for x1 in (0, 1):
        for x2 in (0, 1):
            for n1 in range(first.shape[2]):
                for n2 in range(second.shape[2]):
                    plus = np.outer(first[0, x1, n1], second[0, x2, n2])
                    minus = np.outer(first[1, x1, n1], second[1, x2, n2])
                    joint[0, x1, x2, n1, n2] = np.linalg.norm(plus + minus) ** 2 / 4
                    joint[1, x1, x2, n1, n2] = np.linalg.norm(plus - minus) ** 2 / 4
    return {"joint": joint, "mode_diagnostics": [first_diagnostics, second_diagnostics]}


def summarize(joint: np.ndarray) -> dict:
    codeword = float(joint[0, 0, 0].sum())
    site_probabilities = np.zeros(3)
    for majority, x1, x2 in np.ndindex(2, 2, 2):
        system = majority ^ (x1 & x2)
        colors = np.array([system, system ^ x1, system ^ x2])
        probability = joint[majority, x1, x2].sum()
        site_probabilities += (colors == 0) * probability
    flag_distribution = joint.sum(axis=(0, 3, 4))
    count_marginals = [joint.sum(axis=(0, 1, 2, 4)), joint.sum(axis=(0, 1, 2, 3))]
    mean_bath_counts = [float(np.arange(len(values)) @ values) for values in count_marginals]
    return {
        "historical_joint": codeword,
        "historical_system": float(site_probabilities[0]),
        "historical_M1": float(site_probabilities[1]),
        "historical_M2": float(site_probabilities[2]),
        "historical_both_memories": float(joint[0, 0, 0].sum() + joint[0, 1, 1].sum()),
        "ordinary_consensus": float(joint[:, 0, 0].sum()),
        "opposite_consensus": float(joint[1, 0, 0].sum()),
        "parity_probabilities": flag_distribution.tolist(),
        "mean_bath_counts": mean_bath_counts,
        "normalization_error": float(abs(joint.sum() - 1)),
    }


def _fraction_sqrt_above(value: Fraction, digits: int = 12) -> Fraction:
    scale = 10**digits
    floor_scaled = isqrt(value.numerator * scale**2 // value.denominator)
    return Fraction(floor_scaled + 1, scale)


def uniform_certificate(protocol: dict) -> dict:
    """Exact rational gates for the physical system including the writer."""
    cycles = protocol["repeated"]["max_faults"]
    wait_min, wait_max = protocol["repeated"]["gaps_and_terminal_wait_range"]
    length = min(protocol["correction_chain_spins"])
    assert length >= 4 * cycles * wait_max
    tail = 2 * Fraction(3, 4)**(2 * length)
    inverse_time = Fraction(100, 157)
    phase_error = Fraction(protocol["writer"]["certified_phase_radius"])
    writer_coupling = Fraction(protocol["units"]["kappa_over_g"])
    gates = protocol["strict_uniform_gates"]

    amplitudes = [Fraction(0)]
    comparisons = []
    for count in range(1, cycles + 1):
        harmonic = sum(Fraction(1, j) for j in range(1, count + 1))
        amplitude = inverse_time * harmonic / wait_min + count * tail
        accumulated = sum(amplitudes)
        amplitudes.append(amplitude)
        writer_error = phase_error + writer_coupling * count * wait_max
        dynamics_error = writer_error + 2 * accumulated
        source_error = Fraction(3, 2) * amplitude
        count_error = Fraction(3, 2) * accumulated
        joint_bad_comparison = _fraction_sqrt_above(source_error**2 + count_error**2)
        codeword_bound = 1 - (dynamics_error + source_error)**2
        count_bound = 1 - (dynamics_error + joint_bad_comparison)**2
        comparisons.append({
            "faults": count,
            "worst_one_source_amplitude_bound": float(amplitude),
            "writer_on_to_ready_comparison_vector_error_bound": float(writer_error),
            "historical_codeword_lower_bound": float(codeword_bound),
            "historical_codeword_and_counts_lower_bound": float(count_bound),
            "strict_codeword_gate_pass": codeword_bound > Fraction(
                str(gates["historical_codeword_after_each_of_up_to_four_separated_faults"])
            ),
            "strict_counts_gate_pass": count_bound > Fraction(
                str(gates["historical_codeword_and_exact_expected_bath_counts"])
            ),
        })
    write_lower = 1 - phase_error**2
    return {
        "write_probability_lower_bound": float(write_lower),
        "strict_write_gate_pass": write_lower > Fraction(str(gates["write_two_records"])),
        "log10_boundary_amplitude_bound": log10(2) + 2 * length * log10(0.75),
        "arithmetic": "Exact fractions for all comparisons. Printed decimal bounds are rounded displays.",
        "after_faults": comparisons,
    }
