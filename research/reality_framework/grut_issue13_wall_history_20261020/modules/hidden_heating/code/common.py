"""Shared constants and helpers for module hidden_heating (issue #13 wall-history packet).

Units: GeV, natural units. Reduced Planck mass M_P. Action conventions of
GRUT_CANDIDATE_1.md section 2 (metric -,+,+,+). Internally the amplitude code uses
mostly-minus (+,-,-,-) for Dirac algebra; every physical output is metric independent.
"""
import json
import os
import numpy as np
from scipy.special import zeta

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")
os.makedirs(RESULTS, exist_ok=True)

# ---------------- benchmark (BRIEF.md, GRUT_CANDIDATE_1.md) ----------------
M_P = 2.435323e18          # reduced Planck mass [GeV]
V_BENCH = 1000.0           # v [GeV]
G_X = 0.1
LAM_PHI = 0.01
LAM_X = 0.001
H_D = 0.01
LAM_D = 0.001
M_D0 = G_X * V_BENCH       # 100 GeV
MU0 = G_X * V_BENCH        # matching scale mu_0 = g_X v
TC_OVER_V = 1.547646465    # leading high-T curvature zero (BRIEF); task: use 1.5476 v
TC_BAND = 0.30             # +-30 % sensitivity band requested

# couplings in the quantum-packet normalisation (REPORT.md section 1)
C_PORT = G_X**2            # c = g_X^2   (phi^2 R^2 / 2 with R^2 = q1^2 + q2^2)
D_PORT = H_D**2            # d = h_D^2
ALPHA0 = 0.0               # |X|^2 D^2 portal at mu_0

# hidden-sector thermal/induced-mass coefficients: Pi_phi = sum_i c_i <chi_i^2>
C_I = {"phi": 3 * LAM_PHI, "q1": G_X**2, "q2": G_X**2, "D": H_D**2}
C_AVG = sum(C_I.values()) / 4.0   # equal sharing among the 4 real dof
G_H = 4                            # hidden bosonic dof (phi, q1, q2, D)

# ---------------- Standard Model bath ----------------
G_STAR_SM = 106.75   # nu_R not thermalised (Dirac Yukawa ~1e-12); stated assumption
N_HIGGS_REAL = 4
N_WEYL = 45
N_VECTOR = 12


def c_H(gstar=G_STAR_SM):
    """H = c_H T^2 / M_P in radiation domination."""
    return np.sqrt(np.pi**2 * gstar / 90.0)


def rho_rad(T, g):
    return np.pi**2 * g * T**4 / 30.0


def T_of_rho(rho, g):
    return (30.0 * rho / (np.pi**2 * g)) ** 0.25


# thermal integrals J_n = int_0^inf E^n f(E) dE / T^(n+1)
J3 = {"B": np.pi**4 / 15.0, "F": 7.0 * np.pi**4 / 120.0, "MB": 6.0}
J4 = {"B": 24.0 * zeta(5), "F": (15.0 / 16.0) * 24.0 * zeta(5), "MB": 24.0}


def tc(v=V_BENCH, factor=1.0):
    return TC_OVER_V * v * factor


def dump(name, obj):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=False, default=_conv)
    return path


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


def beta_xi_coeffs(lp=LAM_PHI, lx=LAM_X, c=C_PORT, d=D_PORT, ld=LAM_D, al=ALPHA0):
    """One-loop nonminimal-coupling running, 16 pi^2 d xi_a / dln mu = lambda_aacd (xi_cd - delta_cd/6),
    convention L ⊃ -(1/2) xi_phi R phi^2 - xi_X R |X|^2 - (1/2) xi_D R D^2 (mostly-plus, conformal xi = 1/6).
    Returns matrix A with 16pi^2 beta(xi) = A @ (xi - 1/6) for xi = (xi_phi, xi_X, xi_D)."""
    # lambda_{aa cc} entries in the (phi, q1, q2, D) real basis, summed over q1,q2 for X
    A = np.array([
        [6 * lp, 4 * c, 2 * d],                 # phi row: 6lp, 2c (q1)+2c (q2), 2d
        [2 * c, 6 * lx + 2 * lx, 2 * al],       # q1 row: 2c, 6lx (q1q1) + 2lx (q2q2), 2alpha
        [2 * d, 4 * al, 6 * ld],                # D row: 2d, 2alpha+2alpha, 6 ld
    ])
    return A


def xi_running(mu, mu0=MU0, xi0=(0.0, 0.0, 0.0)):
    """Leading-log xi(mu) from xi(mu0) = 0 (couplings frozen; their running is O(1e-3) relative)."""
    A = beta_xi_coeffs()
    x0 = np.array(xi0)
    b = A @ (x0 - 1.0 / 6.0) / (16 * np.pi**2)
    return x0 + b * np.log(mu / mu0)
