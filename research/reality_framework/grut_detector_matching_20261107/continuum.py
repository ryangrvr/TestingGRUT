"""Independent bosonized massless interval reference for frozen M0-01.

This module imports no lattice code or data. All histories retain ordered
Weyl phases. Numerical quadrature estimates and analytic mode tails are
reported separately.
"""
import itertools
import math

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm


def key(times):
    return ",".join(f"{t:g}" for t in times)


def tv(p, q):
    if set(p) != set(q):
        raise ValueError("Outcome dictionaries do not match")
    return 0.5 * sum(abs(p[k] - q[k]) for k in p)


def gaussian_histories(omega, h2, mean, times, coupling, phases):
    """Ordered vacuum expectation of the full bra/ket instrument expansion."""
    omega, h2 = np.asarray(omega), np.asarray(h2)
    n = len(times)
    cov = {(i, j): np.sum(h2 * np.exp(-1j * omega * (times[i] - times[j])))
           for i in range(n) for j in range(n)}
    # All coefficients are kept complex until the full word sum is formed.
    words = []
    for bra in itertools.product((0, 1), repeat=n):
        for ket in itertools.product((0, 1), repeat=n):
            ops = [(i, -coupling) for i in range(n) if bra[i]]
            ops += [(i, coupling) for i in reversed(range(n)) if ket[i]]
            z = 1j * mean * sum(lam for _, lam in ops)
            for a, (i, lam) in enumerate(ops):
                z -= 0.5 * lam * lam * cov[i, i]
                for j, other in ops[a + 1:]:
                    z -= lam * other * cov[i, j]
            z += 1j * sum(phases[i] * (ket[i] - bra[i]) for i in range(n))
            words.append((bra, ket, np.exp(z)))
    probabilities = {}
    imaginary_residual = 0.0
    for signs in itertools.product((1, -1), repeat=n):
        total = sum(value * math.prod(signs[i] ** (bra[i] + ket[i]) for i in range(n))
                    for bra, ket, value in words) / (4 ** n)
        imaginary_residual = max(imaginary_residual, abs(total.imag))
        probabilities["".join("+" if s == 1 else "-" for s in signs)] = float(total.real)
    return probabilities, float(imaginary_residual)


def interval_spectrum(card, modes, epsabs):
    length, sigma = card["length"], card["sigma"]
    norm = sigma * np.sqrt(2 * np.pi) * math.erf(length / (2 * np.sqrt(2) * sigma))
    def f(z):
        return np.exp(-0.5 * ((z - length / 2) / sigma) ** 2) / norm
    mu = 1 / np.sqrt(np.pi)
    mean, mean_err = quad(lambda z: f(z) * np.cosh(mu * (z - length / 2)) /
                         (4 * np.cosh(mu * length / 2)), 0, length,
                         epsabs=epsabs, epsrel=epsabs)
    integers = np.arange(1, modes + 1)
    momenta = np.pi * integers / length
    coefficients = np.zeros(modes)
    coefficient_errors = np.zeros(modes)
    # Symmetric smearing has exactly zero even sine coefficients.
    for index in range(0, modes, 2):
        val, err = quad(f, 0, length, weight="sin", wvar=momenta[index],
                        epsabs=epsabs, epsrel=epsabs, limit=150)
        coefficients[index] = np.sqrt(2 / length) * val
        coefficient_errors[index] = np.sqrt(2 / length) * err
    omega = np.sqrt(momenta ** 2 + mu ** 2)
    h2 = coefficients ** 2 / (2 * np.pi * omega)
    covariance_quad_estimate = np.sum((2 * abs(coefficients) * coefficient_errors +
                                      coefficient_errors ** 2) / (2 * np.pi * omega))
    b = 2 / norm
    tail_variance_bound = b * b * length * length / (2 * np.pi ** 4 * modes ** 2)
    metadata = {"modes": modes, "mean": mean, "mean_quad_error_estimate": mean_err,
                "variance": float(h2.sum()),
                "covariance_quad_error_estimate": float(covariance_quad_estimate),
                "tail_variance_bound": float(tail_variance_bound),
                "smearing_norm": float(norm), "integration_by_parts_B": float(b),
                "first_gap": float(omega[0])}
    return omega, h2, coefficients, coefficient_errors, metadata


def all_histories(card, omega, h2, mean):
    result = {}
    imaginary = 0.0
    for times in card["protocols"]:
        p, residual = gaussian_histories(omega, h2, mean, times, card["coupling"],
                                         [card["phase"]] * len(times))
        result[key(times)] = p
        imaginary = max(imaginary, residual)
    return result, imaginary


def direct_oscillator(card, cutoff):
    cfg = card["oscillator_verification"]
    annihilation = np.diag(np.sqrt(np.arange(1, cutoff + 1)), 1)
    f = cfg["mean"] * np.eye(cutoff + 1) + cfg["h"] * (annihilation + annihilation.T)
    w = expm(1j * card["coupling"] * f)
    energies = cfg["omega"] * np.arange(cutoff + 1)
    result = {}
    for signs in itertools.product((1, -1), repeat=len(cfg["times"])):
        state = np.zeros(cutoff + 1, complex)
        state[0] = 1
        previous = 0.0
        for time, phase, sign in zip(cfg["times"], cfg["phases"], signs):
            state *= np.exp(-1j * energies * (time - previous))
            state = (state + sign * np.exp(1j * phase) * (w @ state)) / 2
            previous = time
        result["".join("+" if s == 1 else "-" for s in signs)] = float(np.vdot(state, state).real)
    return result


def verify_reference(card):
    cfg = card["oscillator_verification"]
    gaussian, imag = gaussian_histories([cfg["omega"]], [cfg["h"] ** 2], cfg["mean"],
                                        cfg["times"], card["coupling"], cfg["phases"])
    low, high = [direct_oscillator(card, n) for n in cfg["cutoffs"]]
    return {"gaussian": gaussian, "direct_low": low, "direct_high": high,
            "tv_gaussian_vs_high": tv(gaussian, high), "tv_cutoff_change": tv(low, high),
            "imaginary_residual": imag}
