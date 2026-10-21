"""Shared definitions for the thermal_potential module.

Units: v = 1 (the matching unit of Candidate 1; v = 1000 GeV).  All masses, temperatures and fields are in units of v.
Fields: phi (real), X = (q1 + i q2)/sqrt(2), D (real).  Background (phi, X=0, D=0) unless stated otherwise.

Potential in the normalisation of grut_candidate_1_quantum_20261018 (REPORT.md section 1), R^2 = q1^2 + q2^2:
  V0 = C + mp/2 phi^2 + mx/2 R^2 + md/2 D^2 + lp/4 phi^4 + lx/4 R^4 + ld/4 D^4 + c/2 phi^2 R^2 + d/2 phi^2 D^2 + al/2 R^2 D^2
with c = g_X^2, d = h_D^2.  At mu_0 = g_X v = 0.1 v the unchanged benchmark gives
  (lp, lx, c, d, ld, al) = (0.01, 0.001, 0.01, 1e-4, 0.001, 0), (mp, mx, md) = (-0.01, 0, 0.01), C = 0.0025.

Only numpy/scipy.  Deterministic.
"""
import numpy as np
from scipy import integrate, special
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

PI = np.pi
V_GEV = 1000.0                      # v in GeV
MP_RED_GEV = 2.435323e18            # reduced Planck mass (BRIEF)
MP_RED = MP_RED_GEV / V_GEV         # in units of v
MU0 = 0.1                           # mu_0 = g_X v

BENCH = dict(lp=0.01, lx=0.001, ld=0.001, c=0.01, d=1e-4, al=0.0,
             mp=-0.01, mx=0.0, md=0.01, C=0.0025)
KEYS = ["lp", "lx", "ld", "c", "d", "al", "mp", "mx", "md", "C"]

# ----------------------------------------------------------------------------------------------------
# One-loop beta functions (B = 16 pi^2 beta), transcribed from grut_candidate_1_quantum_20261018/rg_beta_functions.txt
# ----------------------------------------------------------------------------------------------------

def betas(p):
    lp, lx, ld, c, d, al = p["lp"], p["lx"], p["ld"], p["c"], p["d"], p["al"]
    mp, mx, md = p["mp"], p["mx"], p["md"]
    B = dict(
        lp=18 * lp**2 + 4 * c**2 + 2 * d**2,
        lx=20 * lx**2 + 2 * c**2 + 2 * al**2,
        ld=18 * ld**2 + 2 * d**2 + 4 * al**2,
        c=6 * lp * c + 8 * lx * c + 8 * c**2 + 2 * d * al,
        d=6 * lp * d + 6 * ld * d + 4 * c * al + 8 * d**2,
        al=8 * lx * al + 6 * ld * al + 2 * c * d + 8 * al**2,
        mp=6 * mp * lp + 4 * mx * c + 2 * md * d,
        mx=2 * mp * c + 8 * mx * lx + 2 * md * al,
        md=2 * mp * d + 4 * mx * al + 6 * md * ld,
        C=0.5 * mp**2 + mx**2 + 0.5 * md**2,
    )
    return {k: B[k] / (16 * PI**2) for k in KEYS}


def run_couplings(mu, p0=BENCH, mu0=MU0):
    """One-loop RG evolution of all ten coefficients from mu0 to mu (t = ln mu/mu0)."""
    if abs(np.log(mu / mu0)) < 1e-15:
        return dict(p0)
    y0 = np.array([p0[k] for k in KEYS])

    def rhs(t, y):
        b = betas(dict(zip(KEYS, y)))
        return np.array([b[k] for k in KEYS])

    sol = solve_ivp(rhs, (0.0, np.log(mu / mu0)), y0, method="DOP853", rtol=1e-13, atol=1e-18)
    return dict(zip(KEYS, sol.y[:, -1]))


# ----------------------------------------------------------------------------------------------------
# Bosonic thermal function J_B(y2) = int_0^inf dx x^2 ln(1 - exp(-sqrt(x^2 + y2)))  (real part for y2 < 0)
# ----------------------------------------------------------------------------------------------------
LN_AB = 2 * np.log(4 * PI) - 2 * np.euler_gamma + 1.5      # ln a_b, a_b = 16 pi^2 exp(3/2 - 2 gamma_E)
_L = np.arange(1, 41)
_SER = (-2 * PI**3.5 * (-1.0) ** _L * special.zeta(2 * _L + 1) / special.factorial(_L + 2)
        * special.gamma(_L + 0.5) / (4 * PI**2) ** (_L + 2))   # coefficient of y2^(l+2)
JB0 = -PI**4 / 45


def jb_nonanalytic(y2):
    """-pi^4/45 + pi^2 y2/12 - (pi/6) Re y2^{3/2} - (y2^2/32) (ln|y2| - ln a_b): the truncated high-T expansion."""
    y2 = np.asarray(y2, dtype=float)
    a = np.abs(y2)
    lg = np.where(a > 0, np.log(np.where(a > 0, a, 1.0)), 0.0)
    return JB0 + PI**2 / 12 * y2 - PI / 6 * np.maximum(y2, 0.0) ** 1.5 - y2**2 / 32 * (lg - LN_AB)


def djb_nonanalytic(y2):
    y2 = np.asarray(y2, dtype=float)
    a = np.abs(y2)
    lg = np.where(a > 0, np.log(np.where(a > 0, a, 1.0)), 0.0)
    return PI**2 / 12 - PI / 4 * np.sqrt(np.maximum(y2, 0.0)) - y2 / 16 * (lg - LN_AB) - y2 / 32


def _series_tail(y2, nterms=40):
    """sum_{l=1}^{nterms} c_l y2^(l+2) = y2^3 P(y2) (Horner)."""
    y2 = np.asarray(y2, dtype=float)
    return y2**3 * np.polyval(_SER[:nterms][::-1], y2)


def _dseries_tail(y2, nterms=40):
    y2 = np.asarray(y2, dtype=float)
    return y2**2 * np.polyval((_SER[:nterms] * (_L[:nterms] + 2))[::-1], y2)


def _jb_bessel(y2, kmax=60):
    y = np.sqrt(y2)
    k = np.arange(1, kmax + 1)[:, None]
    return -np.sum(y2[None, :] * special.kv(2, k * y[None, :]) / k**2, axis=0)


def _djb_bessel(y2, kmax=60):
    y = np.sqrt(y2)
    k = np.arange(1, kmax + 1)[:, None]
    return np.sum(y[None, :] * special.kv(1, k * y[None, :]) / (2 * k), axis=0)


def jb_quad(y2):
    """Direct numerical quadrature of the defining integral (real part for y2 < 0, valid for |y| < 2 pi)."""
    y2 = float(y2)
    opts = dict(epsabs=1e-15, epsrel=1e-14, limit=500)
    if y2 >= 0:
        f = lambda x: x * x * np.log(-np.expm1(-np.sqrt(x * x + y2)))
        return integrate.quad(f, 0, np.inf, **opts)[0]
    a = np.sqrt(-y2)
    if a >= 2 * PI:
        raise ValueError("Re J_B quadrature implemented for |y| < 2 pi only")
    f1 = lambda x: x * x * np.log(2 * abs(np.sin(0.5 * np.sqrt(max(a * a - x * x, 0.0)))))
    f2 = lambda x: x * x * np.log(-np.expm1(-np.sqrt(max(x * x - a * a, 0.0))))
    return integrate.quad(f1, 0, a, **opts)[0] + integrate.quad(f2, a, np.inf, **opts)[0]


class JBSeries:
    """'Exact' J_B through its convergent representations: high-T series (|y2| <= 10, 40 terms, convergent for
    |y| < 2 pi) and the Bessel sum (y2 > 10).  Real part for y2 < 0.  For y2 < -10 falls back to quadrature."""
    name = "series+bessel"

    def __call__(self, y2):
        y2 = np.atleast_1d(np.asarray(y2, dtype=float))
        out = np.empty_like(y2)
        m1 = np.abs(y2) <= 10
        out[m1] = jb_nonanalytic(y2[m1]) + _series_tail(y2[m1])
        m2 = y2 > 10
        if m2.any():
            out[m2] = _jb_bessel(y2[m2])
        m3 = y2 < -10
        if m3.any():
            out[m3] = np.array([jb_quad(t) for t in y2[m3]])
        return out

    def d(self, y2):
        y2 = np.atleast_1d(np.asarray(y2, dtype=float))
        out = np.empty_like(y2)
        m1 = np.abs(y2) <= 10
        out[m1] = djb_nonanalytic(y2[m1]) + _dseries_tail(y2[m1])
        m2 = y2 > 10
        if m2.any():
            out[m2] = _djb_bessel(y2[m2])
        m3 = y2 < -10
        if m3.any():
            h = 1e-5
            out[m3] = np.array([(jb_quad(t + h) - jb_quad(t - h)) / (2 * h) for t in y2[m3]])
        return out


class JBTable:
    """Tabulated quadrature: R(y2) = J_quad(y2) - jb_nonanalytic(y2) is analytic in y2; it is computed by direct
    quadrature on a uniform grid in [-25, 25] and interpolated with a cubic spline.  The known non-analytic terms
    (cubic, log) are added back exactly.  Bessel sum for y2 > 25."""
    name = "quadrature table"

    def __init__(self, ymin=-25.0, ymax=25.0, step=0.01):
        self.grid = np.round(np.arange(ymin, ymax + 0.5 * step, step), 10)
        self.jq = np.array([jb_quad(t) for t in self.grid])
        self.R = self.jq - jb_nonanalytic(self.grid)
        self.spl = CubicSpline(self.grid, self.R)
        self.dspl = self.spl.derivative()
        self.ymin, self.ymax = ymin, ymax

    def __call__(self, y2):
        y2 = np.atleast_1d(np.asarray(y2, dtype=float))
        if (y2 < self.ymin).any():
            raise ValueError("y2 below table")
        out = np.empty_like(y2)
        m = y2 <= self.ymax
        out[m] = jb_nonanalytic(y2[m]) + self.spl(y2[m])
        if (~m).any():
            out[~m] = _jb_bessel(y2[~m])
        return out

    def d(self, y2):
        y2 = np.atleast_1d(np.asarray(y2, dtype=float))
        out = np.empty_like(y2)
        m = y2 <= self.ymax
        out[m] = djb_nonanalytic(y2[m]) + self.dspl(y2[m])
        if (~m).any():
            out[~m] = _djb_bessel(y2[~m])
        return out


# ----------------------------------------------------------------------------------------------------
# Thermal (Debye) masses: Pi_ab = (T^2/24) d_a d_b Tr M^2  (leading high-T, all four real scalars light)
# ----------------------------------------------------------------------------------------------------

def debye_coeffs(p):
    """Pi_i / T^2 for phi, each real X mode (q1, q2), and D."""
    return np.array([(3 * p["lp"] + 2 * p["c"] + p["d"]) / 12,
                     (p["c"] + 4 * p["lx"] + p["al"]) / 12,
                     (p["d"] + 3 * p["ld"] + 2 * p["al"]) / 12])


N_DOF = np.array([1.0, 2.0, 1.0])     # phi, X (two real modes), D


# ----------------------------------------------------------------------------------------------------
# Effective potential along the background (phi, 0, 0) as a function of s = phi^2
# ----------------------------------------------------------------------------------------------------
SCHEMES = ["oneloop", "AE", "Parwani", "DR"]
SCHEME_DOC = {
    "oneloop": "tree + CW(m^2) + exact J_B(m^2/T^2); no resummation (IR-unreliable where X is massless at phi=0)",
    "AE": "Arnold-Espinosa: oneloop - (T/12pi) sum n [Re(m^2+Pi)^{3/2} - Re(m^2)^{3/2}] for all scalars",
    "Parwani": "Parwani: tree + CW(m^2+Pi) + exact J_B((m^2+Pi)/T^2)",
    "DR": ("dimensional-reduction split: AE for X and D; for phi only the hard (nonzero-Matsubara) part of its own loop, "
           "i.e. J_B minus its zero-mode term -(pi/6)Re(y^2)^{3/2}, no phi daisy.  = mean-field potential of the soft phi "
           "zero mode in the 3D effective theory; the soft phi loop is left to the nonperturbative 3D (Ising) treatment"),
}


class Veff:
    def __init__(self, p, mu, scheme, jb):
        assert scheme in SCHEMES
        self.p, self.mu, self.scheme, self.jb = dict(p), mu, scheme, jb
        self.pic = debye_coeffs(p)
        self.ds = np.array([3 * p["lp"], p["c"], p["d"]])

    def masses(self, s):
        p = self.p
        s = np.asarray(s, dtype=float)
        return [p["mp"] + 3 * p["lp"] * s, p["mx"] + p["c"] * s, p["md"] + p["d"] * s]

    @staticmethod
    def _cw(M2, mu):
        a = np.abs(M2)
        lg = np.log(np.where(a > 0, a, 1.0) / mu**2)
        return np.where(a > 0, M2**2 * (lg - 1.5), 0.0) / (64 * PI**2)

    @staticmethod
    def _dcw(M2, mu):
        a = np.abs(M2)
        lg = np.log(np.where(a > 0, a, 1.0) / mu**2)
        return np.where(a > 0, M2 * (lg - 1.0), 0.0) / (32 * PI**2)

    def W(self, s, T):
        """V_eff(phi = sqrt(s), T) with the field-independent constants (C, -pi^2 T^4/90 per dof) dropped."""
        p, mu, jb = self.p, self.mu, self.jb
        s = np.atleast_1d(np.asarray(s, dtype=float))
        out = p["mp"] / 2 * s + p["lp"] / 4 * s**2
        m2 = self.masses(s)
        Pi = self.pic * T**2
        for i in range(3):
            n = N_DOF[i]
            if self.scheme == "Parwani":
                M2 = m2[i] + Pi[i]
            else:
                M2 = m2[i]
            out = out + n * self._cw(M2, mu)
            if T > 0:
                th = jb(M2 / T**2) - JB0
                if self.scheme == "DR" and i == 0:
                    th = th + PI / 6 * np.maximum(M2 / T**2, 0.0) ** 1.5     # remove phi zero mode
                out = out + n * T**4 / (2 * PI**2) * th
                if self.scheme == "AE" or (self.scheme == "DR" and i > 0):
                    out = out - n * T / (12 * PI) * (np.maximum(m2[i] + Pi[i], 0.0) ** 1.5 - np.maximum(m2[i], 0.0) ** 1.5)
        return out

    def dW(self, s, T):
        """dW/ds."""
        p, mu, jb = self.p, self.mu, self.jb
        s = np.atleast_1d(np.asarray(s, dtype=float))
        out = p["mp"] / 2 + p["lp"] / 2 * s
        m2 = self.masses(s)
        Pi = self.pic * T**2
        for i in range(3):
            n, ds = N_DOF[i], self.ds[i]
            M2 = m2[i] + Pi[i] if self.scheme == "Parwani" else m2[i]
            g = self._dcw(M2, mu)
            if T > 0:
                dth = T**2 / (2 * PI**2) * jb.d(M2 / T**2)
                if self.scheme == "DR" and i == 0:
                    dth = dth + T**2 / (2 * PI**2) * PI / 4 * np.sqrt(np.maximum(M2 / T**2, 0.0))
                g = g + dth
                if self.scheme == "AE" or (self.scheme == "DR" and i > 0):
                    g = g - T / (8 * PI) * (np.sqrt(np.maximum(m2[i] + Pi[i], 0.0)) - np.sqrt(np.maximum(m2[i], 0.0)))
            out = out + n * ds * g
        return out

    def V_phi(self, phi, T):
        """V as an explicit function of phi (for the parity check)."""
        return self.W(np.asarray(phi) ** 2, T)


# ----------------------------------------------------------------------------------------------------
# Full-field-space effective potential (phi, q1, D) with q2 = 0 (U(1)_X rotation), 4x4 Hessian eigenvalues
# ----------------------------------------------------------------------------------------------------

def hessian_full(p, phi, q1, D):
    phi, q1, D = np.broadcast_arrays(np.asarray(phi, float), np.asarray(q1, float), np.asarray(D, float))
    lp, lx, ld, c, d, al = p["lp"], p["lx"], p["ld"], p["c"], p["d"], p["al"]
    mp, mx, md = p["mp"], p["mx"], p["md"]
    R2 = q1**2
    H = np.zeros(phi.shape + (4, 4))
    H[..., 0, 0] = mp + 3 * lp * phi**2 + c * R2 + d * D**2
    H[..., 1, 1] = mx + c * phi**2 + lx * (R2 + 2 * q1**2) + al * D**2
    H[..., 2, 2] = mx + c * phi**2 + lx * R2 + al * D**2
    H[..., 3, 3] = md + d * phi**2 + 3 * ld * D**2 + al * R2
    H[..., 0, 1] = H[..., 1, 0] = 2 * c * phi * q1
    H[..., 0, 3] = H[..., 3, 0] = 2 * d * phi * D
    H[..., 1, 3] = H[..., 3, 1] = 2 * al * q1 * D
    return H


def tree_full(p, phi, q1, D):
    R2 = q1**2
    return (p["mp"] / 2 * phi**2 + p["mx"] / 2 * R2 + p["md"] / 2 * D**2 + p["lp"] / 4 * phi**4 + p["lx"] / 4 * R2**2
            + p["ld"] / 4 * D**4 + p["c"] / 2 * phi**2 * R2 + p["d"] / 2 * phi**2 * D**2 + p["al"] / 2 * R2 * D**2)


def V_full(p, mu, scheme, jb, phi, q1, D, T):
    """Full one-loop V_eff at a general background (phi, q1, 0, D); constants dropped.  scheme in oneloop/AE/Parwani."""
    phi, q1, D = np.broadcast_arrays(np.asarray(phi, float), np.asarray(q1, float), np.asarray(D, float))
    H = hessian_full(p, phi, q1, D)
    pic = debye_coeffs(p)
    Pi = np.diag([pic[0], pic[1], pic[1], pic[2]]) * T**2
    ev = np.linalg.eigvalsh(H)
    evP = np.linalg.eigvalsh(H + Pi)
    out = tree_full(p, phi, q1, D)
    use = evP if scheme == "Parwani" else ev
    out = out + np.sum(Veff._cw(use, mu), axis=-1)
    if T > 0:
        sh = use.shape
        out = out + T**4 / (2 * PI**2) * np.sum((jb((use / T**2).ravel()) - JB0).reshape(sh), axis=-1)
        if scheme == "AE":
            out = out - T / (12 * PI) * np.sum(np.maximum(evP, 0) ** 1.5 - np.maximum(ev, 0) ** 1.5, axis=-1)
    return out


def hubble(T_sm, gstar=106.75):
    """H in units of v for visible temperature T_sm (units of v): H = sqrt(pi^2 g*/90) T^2 / M_P."""
    return np.sqrt(PI**2 * gstar / 90) * T_sm**2 / MP_RED


def to_builtin(o):
    if isinstance(o, dict):
        return {str(k): to_builtin(x) for k, x in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_builtin(x) for x in o]
    if isinstance(o, np.ndarray):
        return [to_builtin(x) for x in o.tolist()]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o
