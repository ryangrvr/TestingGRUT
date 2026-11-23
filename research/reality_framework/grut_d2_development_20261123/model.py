"""Exact reduction for the fixed GRUT-D2 spin Hamiltonian.

No fit, Lindblad bath, thermodynamic limit, or exponentially large truncation.
The one-particle wall+bath path has m+L sites. Onsite Delta*Q commutes
with hopping; its phase cancels in all reported number probabilities.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import math
import numpy as np


def size_card(n: int) -> dict:
    m = n - 1
    if m < 1:
        raise ValueError("n must be at least 2")
    s2 = m * (m + 1) * (2 * m + 1) // 6
    return {"n": n, "m": m, "S2": s2, "L": 64 * s2,
            "window_start": 8 * s2, "window_end": 16 * s2}


class WallPath:
    def __init__(self, m: int, L: int, g: float = 1.0):
        if min(m, L) < 1 or g <= 0:
            raise ValueError("positive m, L, g required")
        self.m, self.L, self.g = m, L, g
        self.N = m + L
        theta = np.arange(1, self.N + 1) * np.pi / (self.N + 1)
        self.energies = 2 * g * np.cos(theta)
        self.system_modes = np.sqrt(2 / (self.N + 1)) * np.sin(
            np.arange(1, m + 1)[:, None] * theta[None, :])

    def block(self, t: float) -> np.ndarray:
        """P_S exp(-i T t) P_S from the exact finite-path eigenmodes."""
        v = self.system_modes
        return (v * np.exp(-1j * self.energies * t)) @ v.T

    def rows_from_source(self, t: float, row_count: int) -> np.ndarray:
        """First rows of the first m propagator columns, O(row_count*m*N)."""
        theta = np.arange(1, self.N + 1) * np.pi / (self.N + 1)
        rows = np.sqrt(2 / (self.N + 1)) * np.sin(
            np.arange(1, min(row_count, self.N) + 1)[:, None] * theta[None, :])
        return (rows * np.exp(-1j * self.energies * t)) @ self.system_modes.T


def entropy(probabilities) -> float:
    p = np.asarray(probabilities, dtype=float)
    p = p[p > 0]
    return float(-p @ np.log2(p))


def order_probability(block: np.ndarray, p: float = 0.5) -> float:
    covariance = p * block @ block.conj().T
    return float(np.linalg.det(np.eye(len(block)) - covariance).real)


def number_distribution(eigenvalues) -> np.ndarray:
    dist = np.array([1.0])
    for value in eigenvalues:
        value = float(np.clip(value.real, 0.0, 1.0))
        dist = np.convolve(dist, [1 - value, value])
    return dist


def count_joint(block: np.ndarray, p: float = 0.5) -> np.ndarray:
    """Joint P(initial K, final bath number), enumerating source patterns.

    Conditional on each occupied initial set, wall number is the sum of
    Bernoulli variables from the restricted Slater covariance eigenvalues.
    This also gives an independent check of det(I-p*A*A^dagger).
    """
    m = len(block)
    joint = np.zeros((m + 1, m + 1))
    for k in range(m + 1):
        for chosen in combinations(range(m), k):
            weight = p ** k * (1 - p) ** (m - k)
            if not chosen:
                joint[0, 0] += weight
                continue
            a = block[:, chosen]
            ev = np.linalg.eigvalsh(a.conj().T @ a)
            distribution = number_distribution(ev)
            for walls, probability in enumerate(distribution):
                joint[k, k - walls] += weight * probability
    return joint


def count_information(joint: np.ndarray) -> dict:
    ref = joint.sum(axis=1)
    bath = joint.sum(axis=0)
    return {"success": float(np.trace(joint)),
            "initial_count_entropy_bits": entropy(ref),
            "mutual_information_bits": entropy(ref) + entropy(bath) - entropy(joint.ravel()),
            "normalization_error": float(abs(joint.sum() - 1.0))}


def bounds(m: int, L: int, start: int, end: int, p=Fraction(1, 2)) -> dict:
    """Uniform analytic certificate, exact rational inequalities.

    pi > 3.14 and e < 3 are conservative mathematical constants. We do
    not use rounded spectral data or floating Bessel values to certify.
    """
    s2 = sum(i * i for i in range(1, m + 1))
    K, x = 2 * (L + 1), 2 * end
    if not 4 * x < K:
        raise ValueError("family requires x/K < 1/4")
    eps = 2 * Fraction(3, 4) ** K
    pi_lower = Fraction(314, 100)

    def amplitude(t):
        return 2 * s2 / (pi_lower * t) + m * eps

    a = amplitude(start)
    b = amplitude(2 * start)
    order = 1 - p * a * a
    recovery = 1 - p * (b + 2 * a * a) ** 2 - a * a
    initial_kick_wall_upper = p * a * a
    shifted_count = recovery - initial_kick_wall_upper
    return {"method": "exact rational comparisons, pi>3.14, e<3",
            "finite_boundary_log10_upper": math.log10(2) + K * math.log10(0.75),
            "finite_boundary_first_member_exact_lt_2e_minus24":
                (2 * Fraction(32) ** 130 / math.factorial(130)
                 / (1 - Fraction(32, 131)) < Fraction(2, 10**24)),
            "order_lower_display": float(order),
            "record_count_lower_display": float(order),
            "kick_recovery_lower_display": float(recovery),
            "shifted_count_after_kick_lower_display": float(shifted_count),
            "all_time_order_gate_pass": order > Fraction(9968, 10000),
            "all_time_count_record_gate_pass": order > Fraction(9968, 10000),
            "kick_recovery_gate_pass": recovery > Fraction(9922, 10000)}


def kicked_expected_walls(path: WallPath, tau: float, s: float, p=0.5) -> float:
    """Exact normal-covariance result after endpoint physical X, no reset.

    The kick is the encoded X at wall k=m-1. A Gaussian vacuum determinant
    is not used after this change of excitation sector.
    """
    a, b = path.block(tau), path.block(s)
    q = np.array([2.0] * (path.m - 1) + [1.0])
    transferred = -path.block(tau + s) + (b * q) @ a
    occupied_before = p * float(np.linalg.norm(a[-1]) ** 2)
    return (p * float(np.linalg.norm(transferred) ** 2)
            + (1 - occupied_before) * float(np.linalg.norm(b[:, -1]) ** 2))


def bath_echo_block(block: np.ndarray) -> np.ndarray:
    """System block U(tau) D_B U(tau), with D_B identity on system.

    Full staggering D_j=(-1)^j reverses hopping. Applying it only to bath
    leaves the present system state and all bath-number readouts unchanged.
    """
    signs = (-1.0) ** np.arange(len(block))
    return np.diag(signs) + (block * (1 - signs)) @ block


def echo_energy_variance(path: WallPath, tau: float, delta=4.0, p=0.5) -> float:
    """Variance after the bath-only phase kick; mean work is zero.

    D_B T D_B = -T+F; F has support on system and at most first bath site.
    Use only local rows of U to evaluate Var(dGamma(h))=tr(h²C)-tr(hChC).
    The state remains number-conserving Gaussian, with C²=pC.
    """
    m = path.m
    q = min(path.N, m + 2)
    local_u = path.rows_from_source(tau, q)
    t = np.diag(path.g * np.ones(q - 1), 1)
    t += t.T
    signs = (-1.0)**np.arange(q)
    signs[:m] = 1.0
    f = signs[:, None] * t * signs[None, :] + t
    correction = f @ f - t @ f - f @ t
    t2_mean = p * (path.g**2 * (2*m - 1)
                   + np.trace(local_u.conj().T @ correction @ local_u).real)
    compressed = -t[:m, :m] + local_u.conj().T @ f @ local_u
    t_c_t_c = p*p * np.linalg.norm(compressed)**2
    return float(delta**2*m*p*(1-p) + t2_mean - t_c_t_c)


def subsystem_entropies(block: np.ndarray, p=0.5) -> dict:
    ev = np.clip(np.linalg.eigvalsh(block @ block.conj().T), 0, 1)
    s_wall = sum(entropy([p*v, 1-p*v]) for v in ev)
    s_bath = sum(entropy([p*(1-v), 1-p*(1-v)]) for v in ev)
    return {"physical_system_entropy_bits": 1+s_wall,
            "bath_entropy_bits": s_bath,
            "system_bath_mutual_information_bits": s_wall+s_bath-len(block),
            "active_global_entropy_bits": len(block)+1}


def initial_resources(m: int, L: int, delta=4.0, g=1.0, p=0.5) -> dict:
    mean = delta * m * p
    variance = delta**2 * m * p * (1 - p) + g**2 * (
        2 * (m - 1) * p * (1 - p) + p)
    return {"active_hilbert_log2_dimension": m + 1 + L,
            "reference_dimension": m + 1,
            "initial_active_entropy_bits": m + 1,
            "initial_bath_entropy_bits": 0,
            "prepared_bath_purity_bits": L,
            "mean_energy": mean,
            "energy_variance": variance,
            "disconnected_mean_energy": mean,
            "disconnected_energy_variance": variance - p * g**2,
            "energy_variance_mismatch": p * g**2,
            "disconnected_order_probability": (1 - p)**m,
            "maximally_mixed_bath_mean_energy": delta * (m * p + L / 2),
            "maximally_mixed_bath_active_entropy_bits": m + 1 + L}
