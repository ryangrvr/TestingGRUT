#!/usr/bin/env python3
"""Independent spot check of the leading-order finite-(T, n) scheme stated in charge_phase NOTES.md §3
(verifier v_charge; written from the NOTES' verbal definition, not from the module code).

Scheme (v = 1):  F(phi; T, n) = U_T(phi) + f_gas(n_g; T, M(phi)) + E_cond(n - n_g; M(phi))
  U_T = lp/4 (phi^2-1)^2 + 1/2 (lp/4 + hD^2/12) T^2 phi^2
  M^2 = g^2 phi^2 + PiX T^2,  PiX = g^2/12 + lX/3
  f_gas = -P(mu_g) + mu_g n_g  (exact one-loop ideal Bose gas of X and Xbar, thermal part only), n(mu_g) = n_g
  n_g = min(n, n_max), n_max = n(mu = M)
  E_cond(q; M) = min_c [ q^2/(4c) + M^2 c + lX c^2 ]
n_local(T) = smallest n at which phi = 0 is a local minimum (dF/dphi^2 at phi=0 >= 0).
Computed two ways: (i) analytic derivative with partition response; (ii) finite difference of the full F.
Sensitivity variant B: partition frozen at its phi = 0 value (no dn_max/dM^2 response).
"""
import json, os, time
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

g, lp, lX, hD = 0.1, 0.01, 0.001, 0.01
PiX = g**2 / 12 + lX / 3
Pphi = lp / 4 + hD**2 / 12
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_thermal.json")


def _fm(E, mu, M, p, T):
    a = (p * p + M * M - mu * mu) / (E + mu)
    return 1.0 / np.expm1(a / T)


def _fp(E, mu, T):
    return 1.0 / np.expm1((E + mu) / T)


def Q(fun, T):
    # integrate over p in [0, 50 T]; split to help near p=0
    v1, _ = quad(fun, 0, 0.5 * T, limit=400, epsabs=0, epsrel=1e-11)
    v2, _ = quad(fun, 0.5 * T, 50 * T, limit=400, epsabs=0, epsrel=1e-11)
    return v1 + v2


def n_of(mu, M, T):
    f = lambda p: p * p * (_fm(np.sqrt(p * p + M * M), mu, M, p, T) - _fp(np.sqrt(p * p + M * M), mu, T))
    return Q(f, T) / (2 * np.pi**2)


def X2_of(mu, M, T):
    def f(p):
        E = np.sqrt(p * p + M * M)
        return p * p / (2 * E) * (_fm(E, mu, M, p, T) + _fp(E, mu, T))
    return Q(f, T) / (2 * np.pi**2)


def P_of(mu, M, T):
    def f(p):
        E = np.sqrt(p * p + M * M)
        return p**4 / E * (_fm(E, mu, M, p, T) + _fp(E, mu, T))
    return Q(f, T) / (6 * np.pi**2)


def mu_of_n(ng, M, T):
    if ng <= 0:
        return 0.0
    return brentq(lambda mu: n_of(mu, M, T) - ng, 0.0, M, xtol=1e-15, rtol=1e-13)


def cond(q, M):
    """returns (E_cond, c, mu_c)"""
    if q <= 0:
        return 0.0, 0.0, M
    c = brentq(lambda c: 4 * c * c * (M * M + 2 * lX * c) - q * q, 1e-300, 1e6, xtol=1e-16, rtol=1e-14)
    E = q * q / (4 * c) + M * M * c + lX * c * c
    return E, c, q / (2 * c)


def F_full(phi, T, n):
    M = np.sqrt(g**2 * phi**2 + PiX * T**2)
    UT = lp / 4 * (phi**2 - 1) ** 2 + 0.5 * Pphi * T**2 * phi**2
    nmax = n_of(M, M, T)
    ng = min(n, nmax)
    mug = mu_of_n(ng, M, T) if ng < nmax else M
    fgas = -P_of(mug, M, T) + mug * ng
    Ec = cond(n - ng, M)[0]
    return UT + fgas + Ec


def D_analytic(T, n, freeze_partition=False):
    """dF/d(phi^2) at phi=0."""
    M = np.sqrt(PiX) * T
    nmax = n_of(M, M, T)
    base = -lp / 2 + Pphi * T**2 / 2
    if n <= nmax:
        mug = mu_of_n(n, M, T)
        return base + g**2 * X2_of(mug, M, T)
    E, c, muc = cond(n - nmax, M)
    out = base + g**2 * (X2_of(M, M, T) + c)
    if not freeze_partition:
        # d n_max / d M^2 at fixed T (numerical, central difference in M^2)
        dM2 = 1e-6 * M * M
        dn = (n_of(np.sqrt(M * M + dM2), np.sqrt(M * M + dM2), T) - n_of(np.sqrt(M * M - dM2), np.sqrt(M * M - dM2), T)) / (2 * dM2)
        out += g**2 * (M - muc) * dn
    return out


def n_local(T, freeze=False):
    if T == 0:
        return np.sqrt(8 * lX) * 0.5**1.5
    if D_analytic(T, 0.0, freeze) >= 0:
        return 0.0
    ns = np.linspace(1e-6, 0.12, 121)
    Ds = np.array([D_analytic(T, n, freeze) for n in ns])
    idx = np.where(np.sign(Ds[:-1]) != np.sign(Ds[1:]))[0]
    roots = [brentq(lambda n: D_analytic(T, n, freeze), ns[i], ns[i + 1], xtol=1e-10) for i in idx]
    return roots[0] if roots else None, roots


def Tc0():
    return brentq(lambda T: D_analytic(T, 0.0), 1.0, 2.0, xtol=1e-12)


t0 = time.time()
Tc = Tc0()
Tc_leading = np.sqrt(lp / (lp / 4 + g**2 / 6 + hD**2 / 12))
# analytic one-step estimate with the -M T/(4pi) term: Pi_eff = Pi - g^2 sqrt(PiX)/(2 pi)
Tc_cubic = np.sqrt(lp / (lp / 4 + g**2 / 6 + hD**2 / 12 - g**2 * np.sqrt(PiX) / (2 * np.pi)))

# continuity check of zero-charge transition: phi_min at Tc - 0.0025
def phimin_n0(T):
    r = minimize_scalar(lambda p: F_full(p, T, 0.0), bounds=(0.0, 1.0), method="bounded", options=dict(xatol=1e-7))
    return r.x, r.fun - F_full(0.0, T, 0.0)

pm, dF = phimin_n0(Tc - 0.0025)

rows = []
for T in [0.5, 1.0, 1.3, 1.5, 1.545, 1.55, 1.555, 1.557]:
    nl, roots = n_local(T)
    nlB, rootsB = n_local(T, freeze=True)
    # finite-difference cross check of the root with the full F
    fd = None
    if nl and nl > 0:
        d = 2e-3
        D_fd = lambda n: (F_full(d, T, n) - F_full(0.0, T, n)) / d**2
        try:
            fd = brentq(D_fd, nl * 0.9, min(nl * 1.1, 0.2), xtol=1e-9)
        except ValueError:
            fd = "bracket failed"
    M = np.sqrt(PiX) * T
    rows.append(dict(T=T, n_local=nl, all_roots=roots, n_local_fd_full_F=fd, n_local_frozen_partition=nlB,
                     n_max_gas=n_of(M, M, T)))
    print(rows[-1], flush=True)

out = dict(Tc0=Tc, Tc_leading_T2=Tc_leading, Tc_with_cubic_estimate=Tc_cubic,
           phi_min_at_Tc0_minus_0p0025=pm, dF_at_phimin=dF, rows=rows, runtime_s=time.time() - t0)
with open(OUT, "w") as f:
    json.dump(out, f, indent=1, default=float)
print(json.dumps({k: out[k] for k in ["Tc0", "Tc_leading_T2", "Tc_with_cubic_estimate", "phi_min_at_Tc0_minus_0p0025", "runtime_s"]}, default=float))
