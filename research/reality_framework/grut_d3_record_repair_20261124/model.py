"""Finite paths and the exact charge-sector reduction for GRUT-D3-01."""

from __future__ import annotations

from fractions import Fraction
from math import factorial

import numpy as np
from scipy.fft import dst
from scipy.sparse import coo_matrix, csr_matrix
from scipy.special import jv


def path_state(length: int, time: float) -> np.ndarray:
    """Evolve the origin of the open path with length + 1 sites and g = 1."""
    mode_angles = np.pi * np.arange(1, length + 2) / (length + 2)
    origin_modes = np.sqrt(2 / (length + 2)) * np.sin(mode_angles)
    phased_modes = origin_modes * np.exp(-2j * time * np.cos(mode_angles))
    return dst(phased_modes, type=1, norm="ortho")


def repair_adjacency(writer_length: int, repair_length: int) -> csr_matrix:
    """Joint grid in row-major (x, y) order, with the origin isolated."""
    repair_sites = repair_length + 1
    site_ids = np.arange((writer_length + 1) * repair_sites).reshape(-1, repair_sites)
    writer_edges = site_ids[:-1, :].ravel()[1:]
    repair_edges = site_ids[:, :-1].ravel()[1:]
    sources = np.concatenate((writer_edges, repair_edges))
    destinations = np.concatenate((writer_edges + repair_sites, repair_edges + 1))
    all_sources = np.concatenate((sources, destinations))
    all_destinations = np.concatenate((destinations, sources))
    return coo_matrix(
        (np.ones(len(all_sources)), (all_sources, all_destinations)),
        shape=(site_ids.size, site_ids.size),
    ).tocsr()


def state_after_fault(writer_length: int, repair_length: int, write_time: float) -> np.ndarray:
    state = np.zeros((writer_length + 1, repair_length + 1), dtype=complex)
    state[:, 0] = path_state(writer_length, write_time)
    return state.ravel()


def outcome_probabilities(
    vector: np.ndarray, writer_length: int, repair_length: int
) -> np.ndarray:
    """Rows are blank/recorded M. Columns are erased/restored S."""
    probabilities = abs(vector.reshape(writer_length + 1, repair_length + 1)) ** 2
    return np.array([
        [probabilities[0, 0], probabilities[0, 1:].sum()],
        [probabilities[1:, 0].sum(), probabilities[1:, 1:].sum()],
    ])


def window_certificate(protocol: dict) -> dict:
    """Check strict probability gates using rational upper bounds on amplitudes."""
    writer = protocol["writer"]
    repair = protocol["repair"]
    write_start, write_end = writer["window_gw"]
    repair_start, repair_end = repair["window_gu"]
    if not (
        writer["L"] >= 4 * (write_end + repair_end)
        and repair["L"] >= 4 * repair_end
        and 1 <= repair_start <= repair_end <= 16
    ):
        raise ValueError("The registered rational bounds require these path and time inequalities")

    writer_tail = 2 * Fraction(3, 4) ** (2 * writer["L"])
    repair_tail = 2 * Fraction(3, 4) ** (2 * repair["L"])
    inverse_time_coefficient = 2 / Fraction(314, 100)
    writer_amplitude = inverse_time_coefficient / write_start + writer_tail
    repair_amplitude = inverse_time_coefficient / repair_start + repair_tail
    # On [0,1], use unitarity. On [1,16], use the half-line bound and ln(16) < 3.
    coupling_error = writer_amplitude * (
        Fraction(8, 3) + Fraction(27, 2) * inverse_time_coefficient
        + Fraction(105, 2) * repair_tail
    )
    write_lower = 1 - writer_amplitude**2
    joint_lower = 1 - (writer_amplitude + repair_amplitude + coupling_error) ** 2
    advantage_lower = joint_lower / 2
    gates = protocol["analytic_gates"]
    thresholds = {
        "write": Fraction(str(gates["window_write_accuracy_strict_lower"])),
        "restore_and_retain": Fraction(str(gates["window_intact_restore_and_retain_strict_lower"])),
        "advantage": Fraction(str(gates["window_intact_minus_scrambled_joint_strict_lower"])),
    }
    lower_bounds = {"write": write_lower, "restore_and_retain": joint_lower, "advantage": advantage_lower}
    return {
        "strict_gates_pass": {name: value > thresholds[name] for name, value in lower_bounds.items()},
        "probability_lower_bounds": {name: float(value) for name, value in lower_bounds.items()},
        "writer_origin_amplitude_bound": float(writer_amplitude),
        "repair_origin_amplitude_bound": float(repair_amplitude),
        "finite_repair_boundary_amplitude_bound": float(repair_tail),
        "cut_corner_vector_error_bound": float(coupling_error),
        "inequalities": "Exact fractions with pi>3.14, e<3, ln(16)<3 and the norm estimates in PROOF.md.",
    }


def polynomial_propagation(
    adjacency: csr_matrix, vectors: np.ndarray, times: list[int]
) -> tuple[list[np.ndarray], Fraction]:
    """Evaluate a second propagation using Chebyshev polynomials and Bessel coefficients."""
    degree = 192
    first_omitted = degree + 1
    if max(abs(time) for time in times) > 16:
        raise ValueError("The polynomial tail bound assumes repair time at most 16")
    # For k >= 193, the Bessel series has absolute envelope
    # 32**k/k! times exp(1024/(k+1)), bounded by 3**6 times that term.
    tail_bound = (
        2 * 3**6 * Fraction(32**first_omitted, factorial(first_omitted))
        / (1 - Fraction(32, first_omitted + 1))
    )
    coefficients = (
        jv(np.arange(degree + 1)[None, :], 4 * np.asarray(times)[:, None])
        * (-1j) ** np.arange(degree + 1)[None, :]
    )
    sums = [coefficient[0] * vectors for coefficient in coefficients]
    previous = vectors
    current = (adjacency @ vectors) / 4
    for k in range(1, degree + 1):
        for result, coefficient in zip(sums, coefficients):
            result += (2 * coefficient[k]) * current
        if k < degree:
            previous, current = current, (adjacency @ current) / 2 - previous
    return sums, tail_bound
