"""Reversible local updates and exact endpoint reduction for GRUT-D5."""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from math import exp, isfinite, log, pi

import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom

BLANK = 2
MODES = ("intact", "G inverted", "G scrambled", "G disabled",
         "M inverted", "M scrambled", "M disabled")


@dataclass(frozen=True)
class Work:
    h: int
    gene: int
    coin: int
    fault: int
    mode: str
    a: bytes
    b: bytes
    c: bytes
    memory: int = BLANK
    phase: int = 0
    gene_backup: int = BLANK
    memory_backup: int = BLANK


@dataclass(frozen=True)
class Edge:
    link: int
    operation: str
    site: int = 0


def seed(size: int, h: int, gene: int, coin: int = 0, fault: int = 0, mode: str = "intact") -> Work:
    vacant = bytes([BLANK]) * size
    return Work(h, gene, coin, fault, mode, vacant, vacant, vacant)


def target(size: int, h: int, gene: int) -> bytes:
    return bytes(h ^ (gene & (site % 2)) for site in range(size))


def route(size: int, length: int) -> list[Edge]:
    if size < 8 or size % 4 or length < 8 * size or length % 8:
        raise ValueError("Route needs a size divisible by four and space for all passes")
    quarter = size // 4
    half = size // 2
    edges = [Edge(0, "assay")]
    edges.extend(Edge(1 + i, "grow", i) for i in range(quarter))
    edges.append(Edge(length // 4, "fault_one"))
    edges.extend(Edge(length // 4 + 1 + i - quarter, "grow", i)
                 for i in range(quarter, half))
    edges.append(Edge(length // 2, "control", 1))
    edges.extend(Edge(length // 2 + 1 + half + i, "rebuild", i) for i in range(size))
    edges.append(Edge(3 * length // 4, "fault_three"))
    edges.extend(Edge(3 * length // 4 + 1 + i, "maintain", size - 1 - i)
                 for i in range(size))
    assert edges[-1].link < length
    assert all(left.link < right.link for left, right in zip(edges, edges[1:]))
    return edges


def plateau_bands(size: int, length: int) -> tuple[tuple[int, int], ...]:
    """Inclusive clock-position intervals for the three pulse preparations."""
    return (
        (1 + size // 4, length // 4),
        (length // 4 + 1 + size // 4, length // 2),
        (length // 2 + 1 + size // 2 + size, 3 * length // 4),
    )


def change(tape: bytes, site: int, value: int) -> bytes:
    return tape[:site] + bytes([value]) + tape[site + 1:]


def invert(value: int) -> int:
    return BLANK if value == BLANK else value ^ 1


def flip_site(state: Work, site: int) -> Work:
    return replace(state, a=change(state.a, site, invert(state.a[site])))


def deposit(value: int, color: int) -> int:
    if color == BLANK:
        return value
    if value == color:
        return BLANK
    if value == BLANK:
        return color
    return value


def apply_control(state: Work) -> Work:
    if state.mode == "intact":
        return state
    field = "gene" if state.mode.startswith("G ") else "memory"
    value = getattr(state, field)
    if state.mode.endswith("disabled"):
        backup = field + "_backup"
        return replace(state, **{field: getattr(state, backup), backup: value})
    if state.mode.endswith("inverted") or state.coin == 1:
        return replace(state, **{field: invert(value)})
    return state


def forward(state: Work, edge: Edge, null: bool = False) -> Work:
    operation, site = edge.operation, edge.site
    if operation == "assay":
        return replace(state, memory=deposit(state.memory, state.h))
    if operation == "fault_one":
        return flip_site(state, site) if state.fault == 1 else state
    if operation == "control":
        damaged = flip_site(state, site) if state.fault == 2 else state
        return apply_control(damaged)
    if operation == "fault_three":
        return flip_site(state, site) if state.fault == 3 else state

    base = state.h if null else state.memory
    phase = state.phase ^ (operation == "maintain")
    color = base ^ (state.gene & phase) if BLANK not in (state.gene, base) else BLANK
    if operation == "grow":
        return replace(state, a=change(state.a, site, deposit(state.a[site], color)),
                       phase=state.phase ^ 1)
    if operation == "rebuild":
        return replace(state, a=change(state.a, site, deposit(state.b[site], color)),
                       b=change(state.b, site, state.a[site]), phase=state.phase ^ 1)
    if operation == "maintain":
        return replace(state, a=change(state.a, site, deposit(state.c[site], color)),
                       c=change(state.c, site, state.a[site]), phase=state.phase ^ 1)
    raise ValueError(operation)


def prefixes(state: Work, edges: list[Edge], null: bool = False) -> list[Work]:
    states = [state]
    for edge in edges:
        states.append(forward(states[-1], edge, null))
    return states


def archived_flip(state: Work, word: bytes, stage: int) -> bool:
    array, cell = ((state.b, 0), (state.b, 1), (state.c, 0))[stage - 1]
    return array[cell] == invert(word[cell])


def _encode_log(value: float) -> float | None:
    return float(value) if isfinite(value) else None


def pulse_result(initial: Work, length: int, stage: int, p: Fraction, null: bool = False) -> dict:
    """At T, amplitude for each complete endpoint state is a binomial sum.

    Keep sums in logarithms to retain rare failures when success rounds to one.
    Terms for distinct work states are orthogonal. The endpoint is part of the
    success event, never a condition on which to normalize the outputs.
    """
    size = len(initial.a)
    edges = route(size, length)
    states = prefixes(initial, edges, null)
    clock_log_weights = binom.logpmf(np.arange(length + 1), length, float(p))
    log_normalization = logsumexp(clock_log_weights)
    clock_log_weights -= log_normalization
    start, end = plateau_bands(size, length)[stage - 1]
    outside_band = logsumexp(np.concatenate((clock_log_weights[:start], clock_log_weights[end + 1:])))
    site = (0, 1, 0)[stage - 1]
    time = (length + 1) / 2 * np.arcsin(sqrt_fraction(p))

    groups: dict[Work, float] = {}
    left = 0
    for index, before in enumerate(states):
        right = edges[index].link if index < len(edges) else length
        weight = logsumexp(clock_log_weights[left:right + 1])
        after = flip_site(before, site)
        for edge in edges[index:]:
            after = forward(after, edge, null)
        groups[after] = np.logaddexp(groups.get(after, -np.inf), weight)
        left = right + 1

    word = target(size, initial.h, initial.gene)
    weights = list(groups.values())
    off_endpoint = logsumexp([log(2) + first + second
                              for i, first in enumerate(weights) for second in weights[i + 1:]])
    wrong_at_endpoint = [2 * weight for work, weight in groups.items()
                         if not (work.a == word and archived_flip(work, word, stage))]
    total_failure = logsumexp([off_endpoint, *wrong_at_endpoint])
    return {
        "pulse_time_in_inverse_J": float(time),
        "correct_stage_probability": 1 - exp(outside_band),
        "log_correct_stage_failure": _encode_log(outside_band),
        "joint_completion_and_archive_probability": 1 - exp(total_failure),
        "log_joint_completion_and_archive_failure": _encode_log(total_failure),
        "log_nonarrival_at_endpoint": _encode_log(off_endpoint),
        "log_incorrect_endpoint_score": _encode_log(logsumexp(wrong_at_endpoint)),
        "clock_log_normalization_roundoff": float(log_normalization),
        "orthogonal_endpoint_work_states": len(groups),
    }


def sqrt_fraction(value: Fraction) -> float:
    return float(value)**0.5


def window_certificate(size: int, length: int) -> dict:
    # Both boundaries of each band are at least L/10 from its binomial mean.
    phase_probability_radius = Fraction(1, 64)
    for (start, end), p in zip(plateau_bands(size, length), (Fraction(1, 8), Fraction(3, 8), Fraction(5, 8))):
        assert (p - phase_probability_radius) * length - start >= Fraction(length, 10)
        assert end - (p + phase_probability_radius) * length >= Fraction(length, 10)
    # Hoeffding: eps <= 2 exp(-L/50) < 2 (3/8)^floor(L/50).
    eps = 2 * Fraction(3, 8) ** (length // 50)
    assert eps < Fraction(1, 100000)
    assert 4 * eps < Fraction(1, 100000)
    return {
        "stage_failure_upper": float(eps),
        "completion_and_archive_failure_upper": float(4 * eps),
        "external_p_interval_radius": float(phase_probability_radius),
        "strict_stage_gate_pass": True,
        "strict_completion_gate_pass": True,
        "terminal_time_in_inverse_J": pi * (length + 1) / 4,
        "energy_mean_in_J": 0,
        "energy_variance_in_J_squared": 4 * length / (length + 1)**2,
        "local_exchange_upper_in_J": 1,
        "sector_norm_upper_in_J": 2 * length / (length + 1),
        "head_internal_states": 18144,
        "archive_qutrits": 2 * size,
    }
