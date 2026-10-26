#!/usr/bin/env python3
"""kin_thermo (Claude, own code): fixed inputs [Q], unit conventions and the visible background.

Inputs quoted from Codex INPUTS.json (225028b) light_sector/benchmark and from BRIEF_KINETIC.md.  Everything
dimensional is in GeV.  Hidden-sector thermodynamics is done in units m_X = 1 inside thermo.py; callers multiply
by MX**3 (densities) and MX**4 (energy densities, pressures).
"""
import numpy as np

# ---- light sector and benchmark [Q] ----
GX, LX, LP, V = 0.1, 1e-3, 0.03, 0.924089637375362
MX = GX * V                              # 0.0924089637 GeV
MS = np.sqrt(2 * LP) * V                 # 0.2263548088 GeV  (= sqrt(6) m_X)
RS = MS / MX                             # m_s/m_X = sqrt(6)
LAM3 = 2 * GX**2 * V                     # s X Xbar trilinear, 0.018482 GeV
GAMMA_S_VAC = LAM3**2 * np.sqrt(1 - 4 * MX**2 / MS**2) / (16 * np.pi * MS)   # 1.7333e-5 GeV
YQ, DM_PER_S, F_ALLOW = 4.8294009076988175e-9, 4.4e-10, 0.01
Q_REF = 200000
OMEGA_REF = 0.8611951308217225           # omega/m_X at Q_REF (archived; re-derived in objects.py)
E_OVER_QM_REF = 0.9233161353726057       # E/(Q m_X) at Q_REF (archived; re-derived in objects.py)
EPS = 0.14571744282991028                # |epsilon| of the spot point (INPUTS.json late.epsilon)
N_X, E_X = 2788.1767940228406, 38386.25063336385   # |comoving charge|/H_i^3, conformal energy/H_i^4
H_I = 133349619398.33566                 # GeV
MRED = 2.435e18                          # reduced Planck mass, GeV

# ---- onset (case-L density coordinate) [D], identical to Codex THERMODYNAMIC_CONTROL and kin_source ----
XC = np.sqrt(LP * V * V / (2 * GX * GX))         # critical amplitude |X| where phi turns on, 1.13177 GeV
N_C = 2 * np.sqrt(2 * LX) * XC**3                 # circular-orbit charge density at XC, 0.129665 GeV^3
RHO_C = 3 * LX * XC**4                            # 4.92221e-3 GeV^4
V0 = LP * V**4 / 4                                # 5.46912e-3 GeV^4 (= m_X^2 XC^2 / 2)
N_Q0 = EPS * N_C                                  # net charge density at onset, 0.0188945 GeV^3


def rho_onset(R, include_V0=True):
    """Hidden energy density at the onset coordinate for debris energy ratio R (old-proxy convention)."""
    return RHO_C * (1 + R) + (V0 if include_V0 else 0.0)


# ---- visible background ----
GSTAR_CONST = 106.75
_TS = np.log([1e-6, 5e-4, 1e-3, .1, .15, .3, 1., 5., 100., 200., 1e19])
_G = [3.36, 3.36, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]
_GS = [3.91, 3.91, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]


def gstar(T, table=False):
    """g* (energy).  Default: the constant 106.75 the task specifies.  table=True: the coarse SM table used by
    the old proxies (lsrelic/gp_leftovers), monotone, 10% accuracy."""
    return float(np.interp(np.log(T), _TS, _G)) if table else GSTAR_CONST


def gstar_s(T, table=False):
    return float(np.interp(np.log(T), _TS, _GS)) if table else GSTAR_CONST


def entropy_vis(T, table=False):
    return 2 * np.pi**2 / 45 * gstar_s(T, table) * T**3


def hubble(T, table=False):
    return np.sqrt(np.pi**2 * gstar(T, table) / 90) * T * T / MRED


def T_vis_of_nQ(nQ, ymult=1.0, table=False):
    """Visible temperature at which Y_Q*ymult*s_vis equals the hidden charge density nQ."""
    from scipy.optimize import brentq
    return brentq(lambda T: YQ * ymult * entropy_vis(T, table) - nQ, 1e-7, 1e6, xtol=1e-14, rtol=1e-14)


def lna_of_T(T, table=False):
    """ln a up to a constant (entropy conservation: a^3 g_s T^3 = const)."""
    return -np.log(gstar_s(T, table)**(1 / 3) * T)


def f_sym_from_pairs(pairs_per_net_charge):
    """Symmetric hidden mass relative to the dark-matter benchmark once the net yield is Y_Q."""
    return 2 * MX * YQ * pairs_per_net_charge / DM_PER_S
