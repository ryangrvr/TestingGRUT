"""st_physics: shared tools (own code).

Conformal X variables of a Codex background row, quartic-orbit construction from the invariants (E, L), and the
relative Floquet monodromy of the isolated complex quartic orbit.

The Floquet part is adapted from our own tracks/rm_review/code/quartic_floquet.py (same construction: the orbit
is rebuilt from E and L, radial period and precession by quadrature between the turning points, monodromy from
apoapsis, precession removed by a rotation).  DOP853 replaces Radau for speed; a Radau control is kept.

Conventions (checked against Codex's reheating.py / timing_common.responses):
  bg columns: 0 inflaton, 1 its tau-velocity, 2..7 = (Re x, Im x, Re y, Im y, Re s, Im s), 8..13 their
  tau-velocities, 14 radiation, 15 ln a, 16 H/H_i.  V = lam_X |x|^4 + ... with complex x (|D x|^2 kinetic term),
  so the real canonical radius is rho = sqrt2 |x| and the conformal field chi = sqrt2 a (Re x, Im x) obeys
  chi'' = -lam |chi|^2 chi with lam = lam_X = 1e-3 in a radiation era (a'' = 0, R = 0), ' = d/d eta, d eta = d tau / a.
  chi' = sqrt2 a^2 (vel + H x).
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

LAM = 1e-3
J4 = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])


def conformal_X(bg):
    a = np.exp(bg[:, 15])
    x = bg[:, 2] + 1j * bg[:, 3]
    vel = bg[:, 8] + 1j * bg[:, 9]
    chi = np.sqrt(2) * a * x
    chip = np.sqrt(2) * a ** 2 * (vel + bg[:, 16] * x)
    c = np.stack([chi.real, chi.imag], 1)
    cp = np.stack([chip.real, chip.imag], 1)
    E = 0.5 * np.sum(cp ** 2, 1) + LAM / 4 * np.sum(c ** 2, 1) ** 2
    L = c[:, 0] * cp[:, 1] - c[:, 1] * cp[:, 0]
    return a, c, cp, E, L


class Orbit:
    """Isolated quartic orbit with invariants E, L (conformal units, lam = 1e-3)."""

    def __init__(self, E, L, lam=LAM):
        self.E, self.L, self.lam = float(E), float(L), lam
        pr2 = self.pr2
        # turning points: roots of 2E - L^2/r^2 - lam r^4/2
        rc = (L * L / lam) ** (1 / 6)  # radius where d/dr of effective potential vanishes
        self.rmin = brentq(pr2, 1e-9 * rc, rc)
        hi = rc
        while pr2(hi) > 0:
            hi *= 2
        self.rmax = brentq(pr2, rc, hi)
        c, h = (self.rmin + self.rmax) / 2, (self.rmax - self.rmin) / 2

        def integrand(s, f):
            r = c + h * np.sin(s)
            return f(r) * h * np.cos(s) / np.sqrt(max(pr2(r), 1e-300))

        self.T = 2 * quad(integrand, -np.pi / 2, np.pi / 2, args=(lambda r: 1.0,), epsabs=0, epsrel=1e-13, limit=400)[0]
        self.dth = 2 * quad(integrand, -np.pi / 2, np.pi / 2, args=(lambda r: L / r ** 2,), epsabs=0, epsrel=1e-13,
                            limit=400)[0]
        # circular orbit at the same charge
        self.r_circ = (abs(L) / np.sqrt(lam)) ** (1 / 3)
        self.E_circ = 0.75 * lam ** (1 / 3) * abs(L) ** (4 / 3)
        self.eps = abs(L) / ((4 * E / (3 * lam ** (1 / 3))) ** 0.75)  # L / L_max(E); E_circ/E = eps^(4/3)

    def pr2(self, r):
        return 2 * self.E - self.L ** 2 / r ** 2 - self.lam * r ** 4 / 2

    def _rhs(self, t, y, k):
        q, v, M = y[0:2], y[2:4], y[4:].reshape(4, 4)
        rr = q @ q
        hess = k * k * np.eye(2) + self.lam * (rr * np.eye(2) + 2 * np.outer(q, q))
        dM = np.vstack((M[2:], -hess @ M[:2]))
        return np.r_[v, -self.lam * rr * q, dM.ravel()]

    def monodromy(self, k, method='DOP853', rtol=1e-12, atol=1e-12):
        y0 = np.r_[self.rmax, 0.0, 0.0, self.L / self.rmax, np.eye(4).ravel()]
        s = solve_ivp(self._rhs, (0, self.T), y0, args=(k,), method=method, rtol=rtol, atol=atol)
        yT = s.y[:, -1]
        ang = np.arctan2(yT[1], yT[0])
        Rm = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
        Mon = np.kron(np.eye(2), Rm) @ yT[4:].reshape(4, 4)
        return Mon, ang

    def mu(self, k, **kw):
        Mon, _ = self.monodromy(k, **kw)
        m = float(abs(np.linalg.eigvals(Mon)).max())
        return np.log(m) / self.T, m, float(np.max(abs(Mon.T @ J4 @ Mon - J4)))

    def summary(self):
        return dict(E=self.E, L=self.L, rmin=self.rmin, rmax=self.rmax, rmax_over_rmin=self.rmax / self.rmin,
                    radial_period_eta=self.T, precession=self.dth, r_circ=self.r_circ, E_circ=self.E_circ,
                    E_circ_over_E=self.E_circ / self.E, eps_charge_fraction=self.eps,
                    amplitude_scale_E14=(4 * self.E / self.lam) ** 0.25)


def eta_of(a, C):
    """Radiation-era conformal time up to a constant: a^2 = 2C(tau - tau_*) => d eta = d tau / a = da / C."""
    return a / C
