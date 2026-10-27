#!/usr/bin/env python3
"""st_response skeptic (Claude, own code): a third, differently formulated linear phase-seed response.

Purpose: adversarial check of track st_response's confirmations of Codex b73de64 (grut_source_timing_20261026).
Nothing is imported from Codex's modes.py / response.c / timing_common.py, nor from st_response's st_tangent.py.
Only the model constants and the preparation z0 come from tracks/rp_reheating/code/rp.py ([Q] inputs).

Formulation (different from both Codex's conformal canonical Ward-projected C engine and st_response's conformal
Hamiltonian with an analytic Wirtinger Hessian):
  * cosmic time tau = H_i t, PHYSICAL complex fields z = (X, Y, Sigma), second-order form;
  * background geometry SIMPLIFIED: H = sqrt(rho/3), R = rho - 3p, i.e. the nonminimal factor f = 1 - 2 xi eta2 |X|^2
    and the xi-dependent eta2 terms in R are dropped (they are O(eta2 xi |X|^2) ~ 1e-11; the size of the resulting
    background difference is measured against Codex's archived a, H);
  * the Hessian is NEVER written down: dF = d/deps F(z + eps dz) is obtained by forward-mode dual-number
    differentiation of our own force F = dV/dz^* (checked against a numerical gradient of V);
  * e = 0 (the Coulomb/longitudinal term, an O(e^2/k^2) effect of ~1e-8 of the envelope at e = 1e-7, is omitted);
  * perturbation equations: dz'' = -3 H dz' - (k^2/a^2) dz - dF[z](dz) - xi R (dX, 0, 0), no metric perturbation;
  * seed: dz(0) = dz/d delta from a 4-point (O(h^4)) difference of rp.prepare(b, delta), dz'(0) = 0.
Observables: A = Re(X^* dX)/|X|^2 ;  N_X = 2 a^3 Im(X^* X') (canonical sign, = q0 p1 - q1 p0) ;
  dN_X = 2 a^3 Im(dX^* X' + X^* dX').
"""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

TRACKS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TRACKS / 'rp_reheating' / 'code'))
import rp  # noqa: E402  (model constants and preparation only)

PP = rp.P
LX, LY, LS, KAP, MM, HH, MU = PP['lx'], PP['ly'], PP['ls'], PP['kap'], PP['M'], PP['h'], PP['mu']
XI, ETA2 = rp.XI, rp.ETA2


class D:
    """Forward-mode dual number over the reals: value v (complex scalar), tangents d (complex array, one per k)."""
    __slots__ = ('v', 'd')

    def __init__(self, v, d):
        self.v, self.d = v, d

    def __add__(s, o):
        return D(s.v + o.v, s.d + o.d) if isinstance(o, D) else D(s.v + o, s.d)
    __radd__ = __add__

    def __sub__(s, o):
        return D(s.v - o.v, s.d - o.d) if isinstance(o, D) else D(s.v - o, s.d)

    def __mul__(s, o):
        return D(s.v * o.v, s.v * o.d + s.d * o.v) if isinstance(o, D) else D(s.v * o, s.d * o)
    __rmul__ = __mul__

    def conj(s):
        return D(np.conj(s.v), np.conj(s.d))


def force(x, y, s):
    """F_i = dV/dz_i^* (own hand derivation from V below); works on complex numbers or D."""
    cx = x.conj() if isinstance(x, D) else np.conj(x)
    cy = y.conj() if isinstance(y, D) else np.conj(y)
    cs = s.conj() if isinstance(s, D) else np.conj(s)
    Fx = 2 * LX * (x * cx) * x + 3 * KAP * y * cx * cx
    Fy = MM * MM * y + 2 * LY * (y * cy) * y + KAP * x * x * x + 2 * MU * s * cy
    Fs = 2 * LS * ((s * cs) - HH * HH / 2) * s + MU * y * y
    return Fx, Fy, Fs


def potential(x, y, s):
    x2, y2, s2 = abs(x) ** 2, abs(y) ** 2, abs(s) ** 2
    return (LX * x2 * x2 + MM * MM * y2 + LY * y2 * y2 + LS * (s2 - HH * HH / 2) ** 2
            + 2 * np.real(KAP * np.conj(y) * x ** 3 + MU * np.conj(s) * y * y))


class Model:
    """State: [ln a, i, j, rho_r, X, Y, S (re, im), X', Y', S' (re, im)] (16) + per k [dz (6), dz' (6)]."""

    def __init__(self, b, gam, ks, simplified=True):
        self.b, self.gam = float(b), float(gam)
        self.ks = np.atleast_1d(np.asarray(ks, float)); self.nk = len(self.ks); self.k2 = self.ks ** 2
        self.simplified = simplified

    def geo(self, y):
        z = y[4:10:2] + 1j * y[5:10:2]; v = y[10:16:2] + 1j * y[11:16:2]
        V = potential(*z); K = float(np.sum(np.abs(v) ** 2))
        rhoI = 0.5 * y[2] ** 2 + 0.5 * self.b ** 2 * y[1] ** 2
        rho = rhoI + y[3] + ETA2 * (K + V)
        H = np.sqrt(rho / 3.0)
        R = 2 * self.b ** 2 * y[1] ** 2 - y[2] ** 2 + ETA2 * (4 * V - 2 * K)
        return z, v, H, R

    def rhs(self, t, y):
        z, v, H, R = self.geo(y)
        out = np.empty_like(y)
        out[0] = H; out[1] = y[2]; out[2] = -(3 * H + self.gam) * y[2] - self.b ** 2 * y[1]
        out[3] = self.gam * y[2] ** 2 - 4 * H * y[3]
        T = y[16:].reshape(self.nk, 12)
        dz = T[:, 0:6:2] + 1j * T[:, 1:6:2]; dv = T[:, 6:12:2] + 1j * T[:, 7:12:2]
        Fv = force(D(z[0], dz[:, 0]), D(z[1], dz[:, 1]), D(z[2], dz[:, 2]))
        acc = np.array([-3 * H * v[i] - Fv[i].v for i in range(3)]); acc[0] -= XI * R * z[0]
        out[4:10:2], out[5:10:2] = v.real, v.imag
        out[10:16:2], out[11:16:2] = acc.real, acc.imag
        a2 = np.exp(2 * y[0])
        dacc = np.stack([-3 * H * dv[:, i] - (self.k2 / a2) * dz[:, i] - Fv[i].d for i in range(3)], axis=1)
        dacc[:, 0] -= XI * R * dz[:, 0]
        o = out[16:].reshape(self.nk, 12)
        o[:, 0:6:2], o[:, 1:6:2] = dv.real, dv.imag
        o[:, 6:12:2], o[:, 7:12:2] = dacc.real, dacc.imag
        return out


def seed(b):
    """dz/d delta by a 4-point central difference (O(h^4)) of rp.prepare; returns z0 and dz (complex 3-vectors)."""
    z0, info = rp.prepare(b)
    h = 1e-4
    zs = {s: np.asarray(rp.prepare(b, delta=s * h)[0]) for s in (-2, -1, 1, 2)}
    dz = (8 * (zs[1] - zs[-1]) - (zs[2] - zs[-2])) / (12 * h)
    return np.asarray(z0), dz, info


def initial(model):
    z0, dz, info = seed(model.b)
    V = potential(*z0)
    y = np.zeros(16 + 12 * model.nk)
    y[1] = np.sqrt(2 * (3.0 - ETA2 * V)) / model.b          # H(0) = 1 with the simplified constraint
    y[4:10:2], y[5:10:2] = z0.real, z0.imag
    T = y[16:].reshape(model.nk, 12)
    T[:, 0:6:2], T[:, 1:6:2] = dz.real, dz.imag
    return y, z0, dz, info


def observables(model, Y):
    """Y: (nt, nvar). Returns a, H, N_X, |X|, A (nt, nk), dN (nt, nk)."""
    a = np.exp(Y[:, 0])
    X = Y[:, 4] + 1j * Y[:, 5]; vX = Y[:, 10] + 1j * Y[:, 11]
    T = Y[:, 16:].reshape(len(Y), model.nk, 12)
    dX = T[:, :, 0] + 1j * T[:, :, 1]; dvX = T[:, :, 6] + 1j * T[:, :, 7]
    A = np.real(np.conj(X)[:, None] * dX) / (np.abs(X) ** 2)[:, None]
    N = 2 * a ** 3 * np.imag(np.conj(X) * vX)
    dN = 2 * (a ** 3)[:, None] * np.imag(np.conj(dX) * vX[:, None] + np.conj(X)[:, None] * dvX)
    rhoI = 0.5 * Y[:, 2] ** 2 + 0.5 * model.b ** 2 * Y[:, 1] ** 2
    H = np.array([model.geo(Y[n])[2] for n in range(len(Y))])
    return dict(a=a, H=H, NX=N, absX=np.abs(X), A=A, dN=dN, rhoI=rhoI)


def selftest():
    rng = np.random.default_rng(7)
    worst_F = worst_J = 0.0
    for _ in range(5):
        z = (rng.normal(size=3) + 1j * rng.normal(size=3)) * np.array([60, 3, 25])
        F = np.array(force(*z))
        # F = dV/dz^* = (dV/dRe + i dV/dIm)/2, by a 4-point numerical gradient of V
        g = np.empty(3, complex)
        for i in range(3):
            def Vat(dr, di):
                zz = z.copy(); zz[i] += dr + 1j * di; return potential(*zz)
            h = 1e-3 * max(1, abs(z[i]))
            gr = (8 * (Vat(h, 0) - Vat(-h, 0)) - (Vat(2 * h, 0) - Vat(-2 * h, 0))) / (12 * h)
            gi = (8 * (Vat(0, h) - Vat(0, -h)) - (Vat(0, 2 * h) - Vat(0, -2 * h))) / (12 * h)
            g[i] = 0.5 * (gr + 1j * gi)
        worst_F = max(worst_F, float(np.max(np.abs(F - g)) / np.max(np.abs(F))))
        # dual-number directional derivative vs central difference of F
        d = (rng.normal(size=3) + 1j * rng.normal(size=3))
        Fd = force(*[D(z[i], np.array([d[i]])) for i in range(3)])
        J = np.array([Fd[i].d[0] for i in range(3)])
        h = 1e-4
        Jfd = (8 * (np.array(force(*(z + h * d))) - np.array(force(*(z - h * d))))
               - (np.array(force(*(z + 2 * h * d))) - np.array(force(*(z - 2 * h * d))))) / (12 * h)
        worst_J = max(worst_J, float(np.max(np.abs(J - Jfd)) / np.max(np.abs(J))))
    return dict(force_vs_numerical_gradient_of_V=worst_F, dual_jacobian_vs_fd=worst_J)


if __name__ == '__main__':
    import json
    print(json.dumps(selftest(), indent=1))
