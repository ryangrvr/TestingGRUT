"""common.py -- constants and radiation-era background for module bias_threshold (numpy/scipy only).

Provenance tags: [D] derived/closed form, [Q] quoted (sibling module or literature, source given),
[P] preparation input outside the action, [E] estimate.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(HERE, "..", "..")

# ---------------------------------------------------------------- action (unchanged benchmark, BRIEF.md)
V = 1000.0                     # GeV
LAMBDA_PHI = 0.01
G_X = 0.1
M_PHI = math.sqrt(2 * LAMBDA_PHI) * V        # 141.42 GeV [D]
M_X = G_X * V                                # 100 GeV [D]
SIGMA = 2 * math.sqrt(2 * LAMBDA_PHI) / 3 * V ** 3   # 9.428090416e7 GeV^3 [D]
M_P = 2.435323e18              # reduced Planck mass [BRIEF]

# ---------------------------------------------------------------- calibration (retrospective)
T0 = 2.7255 * 8.617333262e-14  # GeV
RHO_OVER_S_MAX = 0.12 / 2.7440e8   # GeV; Omega h^2 <= 0.12 <-> rho/s <= 4.373e-10 GeV [Q minority_relic mr_common.py]

# ---------------------------------------------------------------- SM degrees of freedom
# Saikawa & Shirai arXiv:1803.01038 Appendix C fitting functions, as implemented in modules/cosmo_observable/code/core.py
# (copied, not imported, to keep this module self-contained).  The fit is used at all T (it gives 104.3-104.9 at
# 1e5-1e10 GeV, against the ideal-gas 106.75; a 0.8% effect on (g ratio)^(1/3)).  The hidden bath (xi = 5e-4)
# adds < 1e-11 to g and is neglected.
_A = [1.0, 1.11724e0, 3.12672e-1, -4.68049e-2, -2.65004e-2, -1.19760e-3, 1.82812e-4, 1.36436e-4, 8.55051e-5,
      1.22840e-5, 3.82259e-7, -6.87035e-9]
_B = [1.43382e-2, 1.37559e-2, 2.92108e-3, -5.38533e-4, -1.62496e-4, -2.87906e-5, -3.84278e-6, 2.78776e-6,
      7.40342e-7, 1.17210e-7, 3.72499e-9, -6.74107e-11]
_C = [1.0, 6.07869e-1, -1.54485e-1, -2.24034e-1, -2.82147e-2, 2.90620e-2, 6.86778e-3, -1.00005e-3, -1.69104e-4,
      1.06301e-5, 1.69528e-6, -9.33311e-8]
_D = [7.07388e1, 9.18011e1, 3.31892e1, -1.39779e0, -1.52558e0, -1.97857e-2, -1.60146e-1, 8.22615e-5, 2.02651e-2,
      -1.82134e-5, 7.83943e-5, 7.13518e-5]
_ME, _MMU, _MPI0, _MPIC = 511e-6, 0.1056, 0.135, 0.140
_M1, _M2, _M3, _M4 = 0.5, 0.77, 1.2, 2.0


def _poly(c, t):
    return sum(ci * t ** i for i, ci in enumerate(c))


def _frho(x): return math.exp(-1.04855 * x) * (1 + 1.03757 * x + 0.508630 * x ** 2 + 0.0893988 * x ** 3)
def _brho(x): return math.exp(-1.03149 * x) * (1 + 1.03317 * x + 0.398264 * x ** 2 + 0.0648056 * x ** 3)
def _fs(x): return math.exp(-1.04190 * x) * (1 + 1.03400 * x + 0.456426 * x ** 2 + 0.0595248 * x ** 3)
def _bs(x): return math.exp(-1.03365 * x) * (1 + 1.03397 * x + 0.342548 * x ** 2 + 0.0506182 * x ** 3)
def _Sfit(x): return 1 + 7.0 / 4.0 * math.exp(-1.0419 * x) * (1 + 1.034 * x + 0.456426 * x ** 2 + 0.0595249 * x ** 3)


def g_rho(T):
    if T >= 0.12:
        t = math.log(T)
        return _poly(_A, t) / _poly(_B, t)
    x = lambda m: m / T
    return (2.030 + 1.353 * _Sfit(x(_ME)) ** (4.0 / 3.0) + 3.495 * _frho(x(_ME)) + 3.446 * _frho(x(_MMU))
            + 1.05 * _brho(x(_MPI0)) + 2.08 * _brho(x(_MPIC)) + 4.165 * _brho(x(_M1)) + 30.55 * _brho(x(_M2))
            + 89.4 * _brho(x(_M3)) + 8209 * _brho(x(_M4)))


def g_s(T):
    if T >= 0.12:
        t = math.log(T)
        return g_rho(T) / (1.0 + _poly(_C, t) / _poly(_D, t))
    x = lambda m: m / T
    return (2.008 + 1.923 * _Sfit(x(_ME)) + 3.442 * _fs(x(_ME)) + 3.468 * _fs(x(_MMU))
            + 1.034 * _bs(x(_MPI0)) + 2.068 * _bs(x(_MPIC)) + 4.16 * _bs(x(_M1)) + 30.55 * _bs(x(_M2))
            + 90 * _bs(x(_M3)) + 6209 * _bs(x(_M4)))


def rho_rad(T):
    return math.pi ** 2 / 30 * g_rho(T) * T ** 4


def s_ent(T):
    return 2 * math.pi ** 2 / 45 * g_s(T) * T ** 3


def H_rad(T):
    """Hubble rate of the visible radiation bath (walls and hidden bath neglected) [D]."""
    return math.sqrt(rho_rad(T) / (3 * M_P ** 2))


def t_rad(T):
    """Radiation-era cosmic time t = 1/(2H) [D] (g changes ignored; adequate at the O(10%) level)."""
    return 1.0 / (2 * H_rad(T))


def T_of_H(H):
    return math.exp(brentq(lambda lt: math.log(H_rad(math.exp(lt))) - math.log(H), math.log(1e-9), math.log(1e19)))


def rho_wall(T, A, sQ=1.0):
    """Scaling network energy density rho_w = A_eff sigma / t, A_eff = sQ A [D]."""
    return sQ * A * SIGMA / t_rad(T)


def T_parity(A, sQ=1.0):
    """rho_w = rho_rad, no backreaction [D]."""
    return math.exp(brentq(lambda lt: math.log(rho_wall(math.exp(lt), A, sQ)) - math.log(rho_rad(math.exp(lt))),
                           math.log(1e-8), math.log(1e3)))


def T_relic_safe(A, f, kappa, sQ=1.0, Ymax=RHO_OVER_S_MAX):
    """Lowest decay temperature with f kappa rho_w(T)/s(T) <= Ymax [D]: rho_w/s grows as 1/T."""
    fn = lambda lt: math.log(f * kappa * rho_wall(math.exp(lt), A, sQ) / s_ent(math.exp(lt))) - math.log(Ymax)
    return math.exp(brentq(fn, math.log(1e-9), math.log(1e8)))


def K_of_T(T, T_f, L_f_phys, A, sQ=1.0):
    """Comoving coarsening factor between formation (T_f, network length V/Area = L_f_phys) and a scaling network at
    T: L(T) = t(T)/(sQ A) (physical), redshifted back to formation with a ~ g_s^(-1/3)/T [D]."""
    L_T = t_rad(T) / (sQ * A)
    a_ratio = (g_s(T) / g_s(T_f)) ** (1.0 / 3.0) * T / T_f      # a(T_f)/a(T)
    return L_T * a_ratio / L_f_phys


def load_json(rel):
    p = os.path.join(MOD, rel)
    return json.load(open(p)) if os.path.exists(p) else None
