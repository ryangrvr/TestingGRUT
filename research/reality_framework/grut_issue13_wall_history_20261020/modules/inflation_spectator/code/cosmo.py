"""Radiation-era bookkeeping for horizon re-entry (estimates; standard inputs quoted).

Visible sector only: the hidden scalars have no portal at mu_0 and stay cold in the histories considered
here, so they do not enter g_*. Instant reheating is an OUTSIDE-ACTION preparation assumption.
"""
import numpy as np
from scipy.optimize import brentq
from common import MP

T0 = 2.7255 * 8.617333262e-14          # GeV (T_CMB, BRIEF)
H0 = 1.4378e-42                         # GeV, h = 0.674 (Planck 2018 value; quoted, used only for e-fold counting)
GS0 = 3.909


def gstar(T):
    """(g_rho, g_s) piecewise-constant SM values (standard textbook steps; estimate)."""
    steps = [(173.0, 106.75), (80.0, 96.25), (4.2, 86.25), (1.2, 75.75), (0.15, 61.75), (0.1, 17.25),
             (5e-4, 10.75)]
    for Tmin, g in steps:
        if T >= Tmin:
            return g, g
    return 3.363, 3.909


def a_of_T(T):
    g, gs = gstar(T)
    return (T0 / T) * (GS0 / gs) ** (1.0 / 3.0)


def H_of_T(T):
    g, gs = gstar(T)
    return np.sqrt(np.pi ** 2 * g / 90.0) * T ** 2 / MP


def efolds_obs_to_reentry(T):
    """ln(a H at T / a0 H0): number of e-folds between horizon exit of today's Hubble scale and of the scale
    that re-enters at T (de Sitter, constant H)."""
    return np.log(a_of_T(T) * H_of_T(T) / H0)


def T_reheat_instant(Hinf):
    g = 106.75
    return (90.0 * Hinf ** 2 * MP ** 2 / (np.pi ** 2 * g)) ** 0.25   # 3 H^2 M_P^2 = (pi^2 g/30) T^4


def Nstar_instant(Hinf):
    """N_* for today's Hubble scale with instant reheating into the visible SM bath (estimate)."""
    Tr = T_reheat_instant(Hinf)
    return np.log(a_of_T(Tr) * Hinf / H0)


def T_reentry(k_over_aendHinf, Hinf):
    """Re-entry temperature of comoving k = (k/(a_end H_inf)) a_end H_inf, instant reheating, RD (estimate)."""
    Tr = T_reheat_instant(Hinf)
    target = np.log(k_over_aendHinf) + np.log(a_of_T(Tr) * Hinf)
    f = lambda lT: np.log(a_of_T(np.exp(lT)) * H_of_T(np.exp(lT))) - target
    return float(np.exp(brentq(f, np.log(1e-12), np.log(Tr) + 1e-9)))   # RD only: invalid below T_eq ~ 1e-9 GeV


def T_wall_dominates_hubble_patch(sigma):
    """T at which one wall of tension sigma per Hubble volume carries sigma/(M_P^2 H) = 1 of rho_tot."""
    Hd = sigma / MP ** 2
    f = lambda lT: np.log(H_of_T(np.exp(lT))) - np.log(Hd)
    return float(np.exp(brentq(f, np.log(1e-9), np.log(1e3))))
