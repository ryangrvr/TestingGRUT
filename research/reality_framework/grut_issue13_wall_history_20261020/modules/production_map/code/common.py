"""production_map / common: constants, Bose integrals, Hartree curvature at phi = 0, rates.

Module production_map of the issue #13 wall-history packet.  numpy/scipy only.  No random numbers here.

Labels used in the JSON outputs: [D] derived here, [C] computed here, [Q] quoted (BRIEF, sibling module,
verifier), [E] order-of-magnitude estimate, [P] outside-action preparation datum.
"""
import json
import math
import os

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
PKT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

# ------------------------------------------------------------------ action benchmark [Q BRIEF]
v = 1000.0
gX, lphi, lX, hD, lD = 0.1, 0.01, 0.001, 0.01, 0.001
gX2, hD2 = gX ** 2, hD ** 2
mD0 = 100.0
MP = 2.435323e18
mX_vac = gX * v

# ------------------------------------------------------------------ cosmology [Q]
gstar = 106.75
gh = 4
zeta3 = 1.2020569031595942
S_COEF = 2 * math.pi ** 2 * gstar / 45          # s_SM = S_COEF T^3
R_COEF = math.pi ** 2 * gstar / 30              # rho_SM = R_COEF T^4
H_COEF = math.sqrt(math.pi ** 2 * gstar / 90)   # H = H_COEF T^2 / M_P

# strict relic bound [Q v_relic VERIFY.md section 2]: every late quantum has m >= m_X
Y_MAX_STRICT = 4.373e-12
Y_MAX_THERMAL_COMPOSITION = 4.30e-12            # [Q relic_yield / v_relic], thermal-reference composition only
XI_MAX_OMEGA = 7.550e-4                          # [Q v_relic] Omega-based xi_max, all species populated
TC_RESUMMED = (1.5571397, 1.5641468)             # [Q thermal_potential, v_thermal CONFIRMED], units of v
TC_LEADING = 1.547646465                         # [Q BRIEF] units of v

# Pi coefficients: Pi_phi = 3 lphi I_phi + gX2 (I_q1 + I_q2) + hD2 I_D   [D, same as hidden_heating]
C_PHI, C_Q, C_D = 3 * lphi, gX2, hD2
C_EFF_EQUAL = (C_PHI + 2 * C_Q + C_D) / 4       # equal sharing among the 4 real dof


def H_of_T(T):
    return H_COEF * T ** 2 / MP


def T_of_H(H):
    return math.sqrt(H * MP / H_COEF)


# ------------------------------------------------------------------ Bose integrals (per real dof)
def _gl_nodes():
    edges = [0.0, 0.05, 0.5, 2.0, 6.0, 15.0, 35.0, 90.0]
    x, w = np.polynomial.legendre.leggauss(48)
    U, W = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        U.append(0.5 * (b - a) * x + 0.5 * (b + a))
        W.append(0.5 * (b - a) * w)
    return np.concatenate(U), np.concatenate(W)


_U, _W = _gl_nodes()


def bose(mT, z):
    """Return (I/T^2, n/T^3, rho/T^4) for one real Bose dof of mass m (mT = m/T) and fugacity z = exp(mu/T) <= 1
    (mu <= 0 <= m, so no condensate).  f = z/(exp(E/T) - z);  I = int d^3k f/((2 pi)^3 E)."""
    k = _U
    E = np.sqrt(k * k + mT * mT)
    f = z / np.expm1(E) if z == 1.0 else z / (np.exp(E) - z)
    pref = 1.0 / (2 * math.pi ** 2)
    I = pref * np.sum(_W * k * k * f / E)
    n = pref * np.sum(_W * k * k * f)
    rho = pref * np.sum(_W * k * k * E * f)
    return I, n, rho


# ------------------------------------------------------------------ Hartree at phi = 0
def hartree_phi0(T, z, mphi_T=0.0, D_in_bath=True, iters=60):
    """Self-consistent (Hartree / one-loop gap) masses at phi = 0 for a kinetically thermal bath with common
    fugacity z for phi, q1, q2 (and D if D_in_bath).  phi is evaluated at the critical point m_phi = mphi_T*T
    (default 0, the restoration boundary).  Returns dict with Pi_phi (GeV^2), masses, n, rho of the bath.

    Gap equations [D]:
      m_q^2 = gX2 I_phi + 4 lX I_q            (per real component, both components equal)
      m_D^2 = mD0^2 + hD2 I_phi + 3 lD I_D
      Pi_phi = 3 lphi I_phi + 2 gX2 I_q + hD2 I_D      ;  restoration iff Pi_phi > lphi v^2
    """
    Iphi, nphi, rphi = bose(mphi_T, z)
    mq2, mD2 = 0.0, mD0 ** 2
    for _ in range(iters):
        Iq, nq, rq = bose(math.sqrt(mq2) / T, z)
        ID, nD, rD = bose(math.sqrt(mD2) / T, z) if D_in_bath else (0.0, 0.0, 0.0)
        mq2_new = gX2 * Iphi * T ** 2 + 4 * lX * Iq * T ** 2
        mD2_new = mD0 ** 2 + hD2 * Iphi * T ** 2 + 3 * lD * ID * T ** 2
        if abs(mq2_new - mq2) < 1e-13 * (mq2 + 1e-300) and abs(mD2_new - mD2) < 1e-13 * mD2:
            mq2, mD2 = mq2_new, mD2_new
            break
        mq2, mD2 = mq2_new, mD2_new
    Iq, nq, rq = bose(math.sqrt(mq2) / T, z)
    ID, nD, rD = bose(math.sqrt(mD2) / T, z) if D_in_bath else (0.0, 0.0, 0.0)
    Pi = (C_PHI * Iphi + 2 * C_Q * Iq + C_D * ID) * T ** 2
    n = (nphi + 2 * nq + nD) * T ** 3
    rho = (rphi + 2 * rq + rD) * T ** 4
    return dict(Pi=Pi, mq=math.sqrt(mq2), mD=math.sqrt(mD2), n=n, rho=rho,
                n_phi=nphi * T ** 3, Iphi=Iphi * T ** 2, Iq=Iq * T ** 2, ID=ID * T ** 2)


def Tc_of_z(z, D_in_bath=True):
    """Kinetic temperature at which Pi_phi(T, z) = lphi v^2 (Hartree, m_phi = 0 at the boundary)."""
    target = lphi * v ** 2
    f = lambda T: hartree_phi0(T, z, D_in_bath=D_in_bath)["Pi"] - target
    # leading massless MB/Bose estimate for bracketing
    Li2 = sum(z ** j / j ** 2 for j in range(1, 400)) if z > 0.5 else z * (1 + z / 4 + z * z / 9)
    T0 = TC_LEADING * v * math.sqrt((math.pi ** 2 / 6) / max(Li2, 1e-300))
    return brentq(f, 0.7 * T0, 1.5 * T0, xtol=1e-12 * T0, rtol=1e-13)


# ------------------------------------------------------------------ elastic rates [Q thermal_potential section 9]
# Gamma_a = sum_b n_b sum_f |M|^2 S_f /(16 pi <s>), <s> = 2<E>^2, n_b per real dof.
COUP_PHI = 18 * lphi ** 2 + 12 * gX2 ** 2 + 6 * hD2 ** 2     # = 3.0e-3
COUP_X = 24 * lX ** 2 + 6 * gX2 ** 2                          # = 6.24e-4


def gamma_el(n_dof, E_mean, coup=COUP_PHI):
    return n_dof * coup / (16 * math.pi * 2 * E_mean ** 2)


def _default(o):
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return float(o)


def dump(name, obj):
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, name), "w") as fh:
        json.dump(obj, fh, indent=1, default=_default)


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)
