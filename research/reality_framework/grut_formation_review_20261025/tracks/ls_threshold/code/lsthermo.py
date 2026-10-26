#!/usr/bin/env python3
"""ls_threshold (Claude, own code): Bose-Einstein thermodynamics of a complex scalar by the Bessel series.

Method (different from both the Gauss-Legendre panels of b448f25/thermo.py and the adaptive quad of
gp_leftovers/equilibrium.py): for one bosonic species of mass m, chemical potential mu < m, temperature T,

  n   = (m^2 T / 2 pi^2)   sum_k  e^{k mu/T} K2(k m/T) / k
  P   = (m^2 T^2 / 2 pi^2) sum_k  e^{k mu/T} K2(k m/T) / k^2
  rho = (m^2 T / 2 pi^2)   sum_k  e^{k mu/T} [3 T K2(k m/T)/k^2 + m K1(k m/T)/k]
  s   = (rho + P - mu n)/T                                    [standard; e.g. Kolb & Turner ch. 3]

Everything is dimensionless here (m = 1, t = T/m, u = mu/m); callers multiply by m^3 or m^4.
kve(nu, x) = K_nu(x) e^x is used, so each term carries e^{-k(1-u)/t}; the series is truncated when that
factor is below e^{-38}. Bessel arrays depend only on t and are cached, so root solves over u are cheap.
"""
import numpy as np
from scipy.special import kve

TWO_PI2 = 2 * np.pi**2
_cache = {}


def bessel(t, umax):
    """k, K2e, K1e arrays (scaled) long enough for every u <= umax at temperature t."""
    d = max((1.0 - umax) / t, 1e-9)
    K = int(min(38.0 / d + 40, 300000))
    key = (t, K)
    if key not in _cache:
        if len(_cache) > 64:
            _cache.clear()
        k = np.arange(1, K + 1, dtype=float)
        _cache[key] = (k, kve(2, k / t), kve(1, k / t))
    return _cache[key]


def species(t, u, B):
    """Dimensionless n, P, rho, s for one species with chemical potential u, given Bessel arrays B."""
    k, k2, k1 = B
    e = np.exp(-k * (1.0 - u) / t)
    a = k2 * e
    n = t / TWO_PI2 * np.dot(a, 1 / k)
    P = t * t / TWO_PI2 * np.dot(a, 1 / (k * k))
    rho = t / TWO_PI2 * (3 * t * np.dot(a, 1 / (k * k)) + np.dot(k1 * e, 1 / k))
    return n, P, rho, (rho + P - u * n) / t


def log_n(t, u, B):
    """ln n for one species, stable for very negative u (factor e^{-(1-u)/t} taken out)."""
    k, k2, _ = B
    e = np.exp(-(k - 1) * (1.0 - u) / t)
    return -(1.0 - u) / t + np.log(t / TWO_PI2 * np.dot(k2 * e, 1 / k))


def net_charge(t, u, B):
    """n(u) - n(-u) for a particle/antiparticle pair at +-u, without cancellation."""
    k, k2, _ = B
    d = np.exp(-k * (1.0 - u) / t) * (-np.expm1(-2 * k * u / t))
    return t / TWO_PI2 * np.dot(k2 * d, 1 / k)


def pair(t, u, B):
    """Particle at +u and antiparticle at -u: dict of n+, n-, nq, rho, P, s (dimensionless)."""
    np_, Pp, rp, sp = species(t, u, B)
    nm, Pm, rm, sm = species(t, -u, B)
    return dict(npl=np_, nmi=nm, nq=net_charge(t, u, B), rho=rp + rm, P=Pp + Pm, s=sp + sm)
