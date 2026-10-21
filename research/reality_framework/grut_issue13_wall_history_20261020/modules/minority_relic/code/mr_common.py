"""Shared inputs for module minority_relic (issue #13 wall-history packet).

Unchanged Candidate 1 action. Inflation, H_inf, N_*, the start condition of the observable patch, the reheating
history (instant reheating at T_RH(H_inf) or a non-instant variant) are OUTSIDE-ACTION PREPARATION DATA [P].

Sibling code is imported READ-ONLY (no bytecode is written into sibling directories):
  * inflation_spectator/code: common.Spectator (couplings run to mu = H, xi included), fp1d.make_grid/odd_chain
    (SqRA Fokker-Planck chain), quoted as [Q-code].
  * cosmo_observable/code/core.py: Saikawa-Shirai g_*rho, g_*s fits (arXiv:1803.01038 App. C), T0, rho_c100,
    omega_m, Planck-like rho_Lambda (calibration only) [Q-code].
"""
import os
import sys
import json
import math
import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.dirname(HERE)
RESULTS = os.path.join(MOD, "results")
PACKET = os.path.dirname(os.path.dirname(MOD))
SPEC_CODE = os.path.join(PACKET, "modules", "inflation_spectator", "code")
COSMO_CODE = os.path.join(PACKET, "modules", "cosmo_observable", "code")
for p in (SPEC_CODE, COSMO_CODE):
    if p not in sys.path:
        sys.path.append(p)

import core as cosmo_core          # noqa: E402  (cosmo_observable)
import fp1d as spec_fp1d           # noqa: E402  (inflation_spectator)
from common import Spectator       # noqa: E402  (inflation_spectator/code/common.py; this module is mr_common)

# ------------------------------------------------------------------------------------------------------------
# Action / BRIEF inputs
# ------------------------------------------------------------------------------------------------------------
V = 1000.0                                   # GeV
LAM_PHI = 0.01
SIGMA = 2.0 * math.sqrt(2.0 * LAM_PHI) / 3.0 * V ** 3     # 9.428090416e7 GeV^3 (BRIEF)
M_P = 2.435323e18                            # GeV reduced Planck mass (BRIEF)
M_H = math.sqrt(2 * LAM_PHI) * V             # 141.42 GeV
M_X = 100.0                                  # GeV
DELTA_WALL = 2.0 / M_H                       # wall thickness (tanh width), GeV^-1

# ------------------------------------------------------------------------------------------------------------
# Verified packet numbers used here (quoted)
# ------------------------------------------------------------------------------------------------------------
PC_SITE = 0.3116077                          # simple-cubic site percolation (inflation_spectator, Wang et al.)
H_PERC = {"N50": 1792.0, "N60": 1616.6}      # GeV, v_inflation CORRECTED values (end of inflation)
H_PERC_TRAP = {"N50": 1424.0, "N60": 1348.0}  # GeV, after post-inflation trapping (v_inflation)
T_DOM_NET = 6.367e-6                         # GeV, cosmo_observable (network, backreaction) [Q]
T_DOM_ISO = 7.97e-6                          # GeV, sigma/(M_P^2 H) = 1 (critic I2) [Q]; recomputed below
RHO_OVER_S_MAX = 0.12 / 2.7440e8             # GeV: Omega h^2 <= 0.12 <-> rho/s <= 4.373e-10 GeV (v_relic strict)
Y_N_MAX_STRICT = 4.373e-12                   # at m = m_X = 100 GeV (v_relic) [Q]
T_BBN = 1.0e-3                               # GeV, onset of n/p freeze-out (criterion choice)
T_BBN_D = 7.0e-5                             # GeV, deuterium bottleneck (variant)
T_REC = cosmo_core.T0 * (1 + cosmo_core.Z_STAR)   # 0.2562 eV (calibration)

# ------------------------------------------------------------------------------------------------------------
# Background (standard LambdaCDM after reheating; hidden sector cold and subdominant) [Q-code + CALIB]
# ------------------------------------------------------------------------------------------------------------
T0 = cosmo_core.T0
RHO_L = cosmo_core.RHO_L_PLANCK              # closes a flat no-wall universe at h = 0.6736 (calibration)


def g_rho(T):
    return cosmo_core.gstar_rho(T)


def g_s(T):
    return cosmo_core.gstar_s(T)


def H_of_T(T):
    return math.sqrt((cosmo_core.rho_rad(T) + cosmo_core.rho_matter(T) + RHO_L) / (3 * M_P ** 2))


def s_of_T(T):
    return cosmo_core.s_ent(T)


H0 = H_of_T(T0)
LN_A0H0 = math.log(H0)                       # a0 = 1


def ln_aH(T):
    return math.log(cosmo_core.a_ratio(T)) + math.log(H_of_T(T))


# tabulate ln(aH) on the decelerating branch for inversion
_lnT = np.linspace(math.log(T0 * 2.2), math.log(1e14), 6000)
_lnaH = np.array([ln_aH(math.exp(x)) for x in _lnT])
assert np.all(np.diff(_lnaH) > 0), "ln(aH) must be monotonic on the tabulated branch"


def nprime_of_T(T):
    """n' = ln(a H / a0 H0) at temperature T: e-folds between the horizon exit of today's Hubble scale and the
    exit of the scale that re-enters at T."""
    return float(np.interp(math.log(T), _lnT, _lnaH) - LN_A0H0)


def T_of_nprime(nprime):
    """Re-entry temperature of a scale that exited n' e-folds after today's Hubble scale (vectorised)."""
    y = np.asarray(nprime) + LN_A0H0
    return np.exp(np.interp(y, _lnaH, _lnT, left=np.nan, right=np.nan))


def T_RH_instant(Hinf):
    """Instant reheating into the SM bath: 3 M_P^2 H^2 = (pi^2/30) g_rho(T) T^4 [P: outside-action assumption]."""
    T = (90.0 * Hinf ** 2 * M_P ** 2 / (math.pi ** 2 * 106.75)) ** 0.25
    for _ in range(30):
        T = (90.0 * Hinf ** 2 * M_P ** 2 / (math.pi ** 2 * g_rho(T))) ** 0.25
    return T


def N0_instant(Hinf):
    """N_0 = ln(a_end H_inf / a0 H0) for instant reheating: e-folds before the end of inflation at which today's
    Hubble scale exits."""
    Tr = T_RH_instant(Hinf)
    return math.log(cosmo_core.a_ratio(Tr)) + math.log(Hinf) - LN_A0H0


def T_dom_iso():
    """sigma/(M_P^2 H) = 1: wall energy of a sphere of radius 1/H equals the radiation energy inside it."""
    from scipy.optimize import brentq
    f = lambda lT: math.log(H_of_T(math.exp(lT))) - math.log(SIGMA / M_P ** 2)
    return math.exp(brentq(f, math.log(1e-8), math.log(1e-3)))


def dump(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
