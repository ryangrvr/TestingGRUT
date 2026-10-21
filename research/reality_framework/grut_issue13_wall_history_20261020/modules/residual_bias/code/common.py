"""Shared constants and helpers for module residual_bias (issue #13 wall-history packet).

Unchanged Candidate 1 action (BRIEF.md): U = (lambda_phi/4)(phi^2 - v^2)^2 + g_X^2 phi^2 |X|^2 + lambda_X |X|^4
+ 1/2 (m_D0^2 + h_D^2 phi^2) D^2 + (lambda_D/4) D^4.  Real basis X = (q1 + i q2)/sqrt2.
Nothing is added to the action.  Hidden-bath production, inflation and reheating are OUTSIDE-ACTION preparation
data [P].  Units GeV.  numpy/scipy only.
"""
import json
import math
import os

import numpy as np

# ---- action benchmark (BRIEF, quoted) ----
V = 1000.0
LPHI = 0.01
GX = 0.1
LX = 0.001
HD = 0.01
LD = 0.001
MD0 = 100.0
MP = 2.435323e18          # reduced Planck mass (BRIEF)
GSTAR = 106.75            # SM g* (BRIEF / task)

# ---- leading high-T Debye coefficients (thermal_potential NOTES section 3, re-derived there; quoted) ----
PI_PHI = LPHI / 4 + GX ** 2 / 6 + HD ** 2 / 12        # 0.004175
PI_Q = GX ** 2 / 12 + LX / 3                            # per real X component, 0.00116667
PI_D = HD ** 2 / 12 + LD / 4                            # 0.000258333
TC_LO = V * math.sqrt(LPHI / PI_PHI)                   # 1547.646 GeV (leading curvature zero)
TC = 1557.5                                            # GeV: Hartree/charge_phase T_c0 = 1.5575 v, inside the
                                                       # resummed range 1.5571-1.5641 v (thermal_potential) [Q]
TC_RANGE = (1557.14, 1564.15)                          # resummed range [Q thermal_potential]
EPS_G = 6.917e-4                                       # Ginzburg reduced temperature t_G [Q thermal_potential]

# ---- contact amplitudes (thermal_potential section 9, re-derived from U) ----
# channel list for the damping of a phi quantum at rest (a) by b -> c + d.  (name, b, c, d, |M|^2, S, multiplicity)
CHANNELS = [
    ("phi phi -> phi phi", "phi", "phi", "phi", (6 * LPHI) ** 2, 0.5, 1),
    ("phi q -> phi q", "q", "phi", "q", (2 * GX ** 2) ** 2, 1.0, 2),
    ("phi phi -> q q", "phi", "q", "q", (2 * GX ** 2) ** 2, 0.5, 2),
    ("phi D -> phi D", "D", "phi", "D", (2 * HD ** 2) ** 2, 1.0, 1),
    ("phi phi -> D D", "phi", "D", "D", (2 * HD ** 2) ** 2, 0.5, 1),
]

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")


def hubble(T_sm):
    """Radiation-era Hubble rate of the visible bath, H = sqrt(pi^2 g*/90) T^2 / M_P (GeV)."""
    return math.sqrt(math.pi ** 2 * GSTAR / 90.0) * T_sm ** 2 / MP


def masses(T_h, Tc=TC):
    """Leading-order quasi-particle masses at phi = 0 in a mu_N = 0 hidden bath of temperature T_h (GeV).
    m_phi^2 = Pi_phi (T^2 - Tc^2) (vanishes at the adopted T_c), m_q^2 = Pi_q T^2, m_D^2 = m_D0^2 + Pi_D T^2."""
    mphi = math.sqrt(max(PI_PHI * (T_h ** 2 - Tc ** 2), 0.0))
    mq = math.sqrt(PI_Q) * T_h
    mD = math.sqrt(MD0 ** 2 + PI_D * T_h ** 2)
    return {"phi": mphi, "q": mq, "D": mD}


def dump(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w") as fh:
        json.dump(obj, fh, indent=1, default=_conv)


def _conv(x):
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.bool_,)):
        return bool(x)
    raise TypeError(type(x))


def log10_m0(log10_r):
    """Second-order exit, Gaussian soft field at the decision time with mean phi_coh and rms sigma on the
    correlation scale: p(+) = P(phi > 0) = (1 + erf(r/sqrt2))/2, so m0 = 2p - 1 = erf(r/sqrt2), r = phi_coh/sigma.
    For small r, m0 = sqrt(2/pi) r.  Returns log10 m0 given log10 r (handles r far below 1e-300)."""
    if log10_r > -2:
        r = 10 ** log10_r
        return math.log10(math.erf(r / math.sqrt(2)))
    return math.log10(math.sqrt(2 / math.pi)) + log10_r
