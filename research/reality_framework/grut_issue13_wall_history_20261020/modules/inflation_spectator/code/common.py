"""Shared definitions for module inflation_spectator.

Unchanged Candidate 1 scalar sector (GRUT_CANDIDATE_1.md sections 2-3; BRIEF.md).
Units: GeV. Metric (-,+,+,+); de Sitter R = 12 H^2.

Nonminimal convention (stated here, zero at mu_0 by matching):
    L  >=  - (1/2) xi_phi R phi^2 - xi_X R |X|^2 - (1/2) xi_D R D^2
so that the potential gains +6 xi_phi H^2 phi^2 and the Hubble-induced mass^2 of phi is 12 xi_phi H^2.
With X = (q1 + i q2)/sqrt2 this is (1/2) xi_X R (q1^2 + q2^2): every real component has the same
"(1/2) xi R z^2" normalization. (GRUT_CANDIDATE_1.md writes "xi_phi R phi^2" without the 1/2; its
coefficient is xi_phi/2 in this convention. Both vanish at mu_0.)

Inflation, H_inf, N_* and the inflaton are OUTSIDE-ACTION preparation data. Nothing here adds an
operator to the action.
"""
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

V = 1000.0                       # GeV, benchmark v
MU0 = 100.0                      # GeV, matching scale mu_0 = g_X v
MP = 2.435323e18                 # GeV, reduced Planck mass (BRIEF)
K16 = 16.0 * np.pi ** 2
SIGMA_TREE = 2.0 * np.sqrt(2.0 * 0.01) / 3.0 * V ** 3   # GeV^3, planar tree wall tension (BRIEF)
PC_SITE_SC = 0.3116077            # 3D simple-cubic site percolation threshold, Wang et al. arXiv:1302.0421 (abstract)
DIFF = 1.0 / (8.0 * np.pi ** 2)   # stochastic diffusion constant per e-fold in units H = 1
BETA_FP = 8.0 * np.pi ** 2 / 3.0  # P_eq ~ exp(-BETA_FP * V/H^4)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")

# parameter vector at mu_0 (MSbar), normalization of grut_candidate_1_quantum_20261018/RG_DERIVATION.md:
# V0 = C + mp phi^2/2 + mx R^2/2 + md D^2/2 + lp phi^4/4 + lx R^4/4 + ld D^4/4 + c phi^2 R^2/2 + d phi^2 D^2/2
#      + alpha R^2 D^2/2,  R^2 = q1^2 + q2^2 = 2|X|^2.
NAMES = ["lp", "lx", "ld", "c", "d", "alpha", "mp", "mx", "md", "xp", "xx", "xd"]
P0 = np.array([0.01, 0.001, 0.001, 0.01, 1e-4, 0.0,
               -0.01 * V ** 2, 0.0, (0.1 * V) ** 2,
               0.0, 0.0, 0.0])


def B_all(p):
    """16 pi^2 beta for every parameter. Quartic and mass betas: rg_beta_functions.txt (quoted).
    xi betas: derived in NOTES.md section 1 from B_{xi_ab} = sum_cd Lambda_abcd (xi_cd - delta_cd/6)."""
    lp, lx, ld, c, d, a, mp, mx, md, xp, xx, xd = p
    s = 1.0 / 6.0
    return np.array([
        18 * lp ** 2 + 4 * c ** 2 + 2 * d ** 2,
        20 * lx ** 2 + 2 * c ** 2 + 2 * a ** 2,
        18 * ld ** 2 + 2 * d ** 2 + 4 * a ** 2,
        6 * lp * c + 8 * lx * c + 8 * c ** 2 + 2 * d * a,
        6 * lp * d + 6 * ld * d + 4 * c * a + 8 * d ** 2,
        8 * lx * a + 6 * ld * a + 2 * c * d + 8 * a ** 2,
        6 * mp * lp + 4 * mx * c + 2 * md * d,
        2 * mp * c + 8 * mx * lx + 2 * md * a,
        2 * mp * d + 4 * mx * a + 6 * md * ld,
        6 * lp * (xp - s) + 4 * c * (xx - s) + 2 * d * (xd - s),
        2 * c * (xp - s) + 8 * lx * (xx - s) + 2 * a * (xd - s),
        2 * d * (xp - s) + 4 * a * (xx - s) + 6 * ld * (xd - s),
    ])


class Running:
    """One-loop running of all pure-scalar parameters including the xi's, t = ln(mu/mu_0)."""

    def __init__(self, tmin=-3.0, tmax=np.log(1e15 / MU0)):
        f = lambda t, y: B_all(y) / K16
        self.up = solve_ivp(f, (0.0, tmax), P0, method="DOP853", rtol=1e-12, atol=1e-16, dense_output=True)
        self.dn = solve_ivp(f, (0.0, tmin), P0, method="DOP853", rtol=1e-12, atol=1e-16, dense_output=True)

    def __call__(self, mu):
        t = np.log(mu / MU0)
        y = self.up.sol(t) if t >= 0 else self.dn.sol(t)
        return dict(zip(NAMES, y))


_RUN = None


def couplings(mu):
    global _RUN
    if _RUN is None:
        _RUN = Running()
    return _RUN(mu)


class Spectator:
    """Dimensionless stochastic potential w = V/H^4 for given H (GeV) and RG scale mu = kmu*H.
    Field variables y = phi/H, r = rho/H with rho^2 = 2|X|^2. D is set to zero (see NOTES)."""

    def __init__(self, H, kmu=1.0, xi_on=True):
        self.H = H
        p = couplings(kmu * H)
        self.p = p
        self.lp, self.lx, self.c = p["lp"], p["lx"], p["c"]
        self.xp = p["xp"] if xi_on else 0.0
        self.xx = p["xx"] if xi_on else 0.0
        self.a2 = p["mp"] / H ** 2 + 12.0 * self.xp          # phi mass^2 coefficient at origin / H^2
        self.b2 = p["mx"] / H ** 2 + 12.0 * self.xx          # X (per real comp.) mass^2 at phi=0 / H^2
        self.vH = V / H

    def w(self, y, r=0.0):
        return (0.5 * self.a2 * y ** 2 + 0.25 * self.lp * y ** 4
                + 0.5 * self.c * y ** 2 * r ** 2 + 0.25 * self.lx * r ** 4 + 0.5 * self.b2 * r ** 2)

    def s(self, y, r=0.0):
        return BETA_FP * self.w(y, r)

    def ymin(self):
        return np.sqrt(-self.a2 / self.lp) if self.a2 < 0 else 0.0

    def barrier_B(self):
        """Hawking-Moss / stochastic exponent 8 pi^2 DeltaV/(3H^4) along X = D = 0."""
        if self.a2 >= 0:
            return 0.0
        return BETA_FP * self.a2 ** 2 / (4.0 * self.lp)


def dump(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=False, default=float)
