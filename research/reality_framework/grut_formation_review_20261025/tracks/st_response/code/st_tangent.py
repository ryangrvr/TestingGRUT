#!/usr/bin/env python3
"""st_response track (Claude, own code): independent linear phase-seed response on the quiet reheating histories.

Review of Codex b73de64 (grut_source_timing_20261026), issue #13.  Nothing from Codex's modes.py, response.c or
timing_common.py is imported for the dynamics.  The background is OUR rp_reheating engine's equations (rp.py /
rpreheat.c, mode 0, oscillating inflaton, H reconstructed from the Friedmann constraint at every evaluation) transcribed
to Python here so that background and tangent are integrated together; the transcription is checked against the rp C
engine (rp.run) separately in check_background().

Background state y (own layout, tau = H_i t):
  0 ln a | 1 i = I/M_r | 2 j = di/dtau | 3..8 X, Y, Sigma (re, im) | 9..14 their tau-velocities | 15 rho_r | 16 H_aux
  H from 3 f H^2 + 3 H fdot = rho (positive root), f = 1 - 2 xi eta2 |X|^2;  H_aux integrates Raychaudhuri (check only).
Tangent per comoving k (conformal canonical variables, own derivation in NOTES.md section 1):
  chi = a r, r = sqrt2 (Re, Im) of (X, Y, Sigma); q = delta chi, p canonical;  g = T chi, g' = T chi', T = diag(J, 3J, 6J)
  H2 = p^2/2 + q^T (k^2 + M) q / 2 + (e^2 / 2 k^2) (g.p - g'.q)^2,
  M = a^2 [V_rr + xi R P_X] - (a^2 R / 6) 1,   V_rr = analytic Hessian of rp.potential in r (own Wirtinger derivation)
  dq/dtau = (p + c g j)/a,  dp/dtau = (-(k^2 + M) q + c g' j)/a,  c = e^2/k^2,  j = g.p - g'.q.
  The Hessian is used directly (no Ward-identity projection): an independent numerical route from Codex's C engine.
Seed: q0 = d chi/d delta at a = 1 from a central difference of rp.prepare(b, delta = +-1e-5) (Sigma phase -pi/2 + delta,
  X real positive, Y and radii re-solved); p0 = q0 (zero physical velocity perturbation: chi' = a^2 (rdot + H r) = r at
  a = H = 1 and zero velocities; the seed current g.p0 - g'.q0 = 0 because g' = g at tau = 0, so q' = p0).
Observables: A = chi_X . q_X / |chi_X|^2 (fractional |X| response per radian);
  N_X = chi_0 chi_1' - chi_1 chi_0' (= q0 p1 - q1 p0 for the background; the NEGATIVE of rp's NX = -2 a^3 Im(X* Xdot));
  dN_X = dq0 chi1' + chi0 dp1 - dq1 chi0' - chi1 dp0 (canonical p, gauge covariant); charge fraction = dN_X / |<N_X>_{64..80}|.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
TRACKS = HERE.parents[1]
sys.path.insert(0, str(TRACKS / 'rp_reheating' / 'code'))
import rp  # noqa: E402  (Claude's own reheating engine)

P, XI, ETA2 = rp.P, rp.XI, rp.ETA2
E_GAUGE = 1e-7          # [Q] tracks/hw_gates/code/spotmodel.py SPOT['e'] (heavy-window spot point)
LX, LY, LS, KAP, MM, HH, MU = P['lx'], P['ly'], P['ls'], P['kap'], P['M'], P['h'], P['mu']
CHARGE = (1.0, 3.0, 6.0)
SQ2 = np.sqrt(2.0)


def pot_force(x, y, s):
    """Scalar copy of rp.potential (V and F_i = dV/dz_i^*); checked against rp.potential in selftest()."""
    x2, y2, s2 = (x * x.conjugate()).real, (y * y.conjugate()).real, (s * s.conjugate()).real
    u = s2 - HH * HH / 2
    V = LX * x2 * x2 + MM * MM * y2 + LY * y2 * y2 + LS * u * u + 2 * (KAP * y.conjugate() * x ** 3 + MU * s.conjugate() * y * y).real
    Fx = 2 * LX * x2 * x + 3 * KAP * y * x.conjugate() ** 2
    Fy = MM * MM * y + 2 * LY * y2 * y + KAP * x ** 3 + 2 * MU * s * y.conjugate()
    Fs = 2 * LS * u * s + MU * y * y
    return V, (Fx, Fy, Fs)


def hessian_r(x, y, s, xiR):
    """6x6 Hessian of V + xiR |X|^2 in r = sqrt2 (Re, Im).  A_ij = d^2V/dz_i^* dz_j, B_ij = d^2V/dz_i^* dz_j^*;
    block (i, j) = [[Re(A+B), -Im(A-B)], [Im(A+B), Re(A-B)]] (own derivation, NOTES.md 1.2)."""
    A = np.zeros((3, 3), complex); B = np.zeros((3, 3), complex)
    A[0, 0] = 4 * LX * abs(x) ** 2 + xiR
    B[0, 0] = 2 * LX * x * x + 6 * KAP * y * x.conjugate()
    A[0, 1] = 3 * KAP * x.conjugate() ** 2; A[1, 0] = 3 * KAP * x ** 2
    A[1, 1] = MM * MM + 4 * LY * abs(y) ** 2
    B[1, 1] = 2 * LY * y * y + 2 * MU * s
    A[1, 2] = 2 * MU * y.conjugate(); A[2, 1] = 2 * MU * y
    A[2, 2] = 4 * LS * abs(s) ** 2 - LS * HH * HH
    B[2, 2] = 2 * LS * s * s
    Hm = np.empty((6, 6))
    Hm[0::2, 0::2] = (A + B).real; Hm[0::2, 1::2] = -(A - B).imag
    Hm[1::2, 0::2] = (A + B).imag; Hm[1::2, 1::2] = (A - B).real
    return Hm


class Model:
    def __init__(self, b, gam, ks, e=E_GAUGE, xi=XI, eta2=ETA2):
        self.b, self.gam, self.xi, self.eta2, self.e = float(b), float(gam), xi, eta2, e
        self.ks = np.atleast_1d(np.asarray(ks, float)); self.nk = len(self.ks)
        self.k2 = self.ks ** 2
        self.c = (e * e / self.k2) if e else np.zeros(self.nk)

    def geometry(self, y):
        x, yy, s = complex(y[3], y[4]), complex(y[5], y[6]), complex(y[7], y[8])
        vx, vy, vs = complex(y[9], y[10]), complex(y[11], y[12]), complex(y[13], y[14])
        V, F = pot_force(x, yy, s)
        K = abs(vx) ** 2 + abs(vy) ** 2 + abs(vs) ** 2
        x2 = abs(x) ** 2
        xi, eta2 = self.xi, self.eta2
        f = 1 - 2 * xi * eta2 * x2
        fd = -4 * xi * eta2 * (x.conjugate() * vx).real
        rhoI = 0.5 * y[2] ** 2 + 0.5 * self.b ** 2 * y[1] ** 2
        pI = 0.5 * y[2] ** 2 - 0.5 * self.b ** 2 * y[1] ** 2
        rho = rhoI + eta2 * (K + V) + y[15]
        H = (-3 * fd + np.sqrt(9 * fd * fd + 12 * f * rho)) / (6 * f)
        R = (rhoI - 3 * pI + eta2 * (4 * V - 2 * K + 12 * xi * abs(vx) ** 2 - 12 * xi * (x.conjugate() * F[0]).real)) \
            / (f + 12 * xi * xi * eta2 * x2)
        return dict(z=(x, yy, s), v=(vx, vy, vs), V=V, F=F, K=K, f=f, fd=fd, H=H, R=R, rho=rho, rhoI=rhoI)

    def bg_rhs(self, y, g):
        d = np.zeros(17)
        H, b, gam, xi, eta2 = g['H'], self.b, self.gam, self.xi, self.eta2
        d[0] = H; d[1] = y[2]; d[2] = -(3 * H + gam) * y[2] - b * b * y[1]
        d[15] = gam * y[2] ** 2 - 4 * H * y[15]
        acc = [-3 * H * g['v'][k] - g['F'][k] for k in range(3)]
        acc[0] -= xi * g['R'] * g['z'][0]
        for k in range(3):
            d[3 + 2 * k], d[4 + 2 * k] = g['v'][k].real, g['v'][k].imag
            d[9 + 2 * k], d[10 + 2 * k] = acc[k].real, acc[k].imag
        x, vx = g['z'][0], g['v'][0]
        fdd = -4 * xi * eta2 * (abs(vx) ** 2 + (x.conjugate() * acc[0]).real)
        d[16] = -(y[2] ** 2 + 2 * eta2 * g['K'] + 4 * y[15] / 3 + fdd - H * g['fd']) / (2 * g['f'])
        return d

    def chis(self, y, g):
        a = np.exp(y[0]); H = g['H']
        chi = np.empty(6); chip = np.empty(6)
        for k in range(3):
            z, v = g['z'][k], g['v'][k]
            chi[2 * k], chi[2 * k + 1] = SQ2 * a * z.real, SQ2 * a * z.imag
            w = SQ2 * a * a * (v + H * z)
            chip[2 * k], chip[2 * k + 1] = w.real, w.imag
        return a, chi, chip

    def gauge(self, chi):
        gv = np.empty(6)
        for i, Q in enumerate(CHARGE):
            gv[2 * i] = -Q * chi[2 * i + 1]; gv[2 * i + 1] = Q * chi[2 * i]
        return gv

    def Mmat(self, y, g):
        a = np.exp(y[0]); a2 = a * a
        Hm = hessian_r(*g['z'], self.xi * g['R'])
        return a2 * Hm - (a2 * g['R'] / 6) * np.eye(6)

    def rhs(self, t, y):
        g = self.geometry(y)
        out = np.empty_like(y)
        out[:17] = self.bg_rhs(y, g)
        a, chi, chip = self.chis(y, g)
        M = self.Mmat(y, g); gv = self.gauge(chi); gp = self.gauge(chip)
        T = y[17:].reshape(self.nk, 12); q, p = T[:, :6], T[:, 6:]
        j = p @ gv - q @ gp
        cj = (self.c * j)[:, None]
        dq = (p + cj * gv) / a
        dp = (-self.k2[:, None] * q - q @ M + cj * gp) / a      # M symmetric
        o = out[17:].reshape(self.nk, 12); o[:, :6] = dq; o[:, 6:] = dp
        return out


def initial(b, step=1e-5):
    z0, info = rp.prepare(b)
    zp, ip = rp.prepare(b, delta=step); zm, im = rp.prepare(b, delta=-step)
    dz = (np.asarray(zp) - np.asarray(zm)) / (2 * step)
    q0 = np.empty(6); q0[0::2] = SQ2 * dz.real; q0[1::2] = SQ2 * dz.imag
    y = np.zeros(17)
    V, _ = rp.potential(z0)
    f = 1 - 2 * XI * ETA2 * abs(z0[0]) ** 2
    need = 3 * f - ETA2 * V
    y[1] = np.sqrt(2 * need) / b
    for k in range(3):
        y[3 + 2 * k], y[4 + 2 * k] = z0[k].real, z0[k].imag
    y[16] = 1.0
    return y, q0, dict(prep=info, prep_plus_residual=ip['residual'], prep_minus_residual=im['residual'], z0=z0, dz=dz)


def observables(model, Y):
    """Y: (nt, 17 + 12 nk).  Returns a, H, NX (background, q0p1-q1p0 sign), amplitude (nt, nk), dNX (nt, nk), extras."""
    nt = Y.shape[0]
    a = np.empty(nt); H = np.empty(nt); NX = np.empty(nt); R = np.empty(nt); Haux = np.empty(nt); absX = np.empty(nt)
    amp = np.empty((nt, model.nk)); dN = np.empty((nt, model.nk)); gauss = np.empty(nt)
    for n in range(nt):
        y = Y[n]; g = model.geometry(y); a_, chi, chip = model.chis(y, g)
        T = y[17:].reshape(model.nk, 12); q, p = T[:, :6], T[:, 6:]
        a[n], H[n], R[n], Haux[n], absX[n] = a_, g['H'], g['R'], y[16], abs(g['z'][0])
        NX[n] = chi[0] * chip[1] - chi[1] * chip[0]
        amp[n] = (q[:, 0] * chi[0] + q[:, 1] * chi[1]) / (chi[0] ** 2 + chi[1] ** 2)
        dN[n] = q[:, 0] * chip[1] + chi[0] * p[:, 1] - q[:, 1] * chip[0] - chi[1] * p[:, 0]
        nq = [chi[2 * i] * chip[2 * i + 1] - chi[2 * i + 1] * chip[2 * i] for i in range(3)]
        gauss[n] = (nq[0] + 3 * nq[1] + 6 * nq[2]) / max(1.0, abs(nq[0]) + 3 * abs(nq[1]) + 6 * abs(nq[2]))
    return dict(a=a, H=H, NX=NX, R=R, Haux=Haux, absX=absX, amp=amp, dN=dN, gauss=gauss)


def run(b, gam, ks, times, rtol=1e-10, atol=1e-13, method='DOP853', max_step=0.05, e=E_GAUGE, y0=None, t0=0.0):
    model = Model(b, gam, ks, e=e)
    ybg, q0, info = initial(b)
    if y0 is None:
        y0 = np.r_[ybg, np.tile(np.r_[q0, q0], model.nk)]
    t = time.time()
    sol = solve_ivp(model.rhs, (t0, times[-1]), y0, method=method, t_eval=times, rtol=rtol, atol=atol,
                    max_step=max_step, dense_output=False)
    if not sol.success:
        raise RuntimeError(sol.message)
    return model, sol, info, time.time() - t


def selftest():
    rng = np.random.default_rng(1)
    worst_V = worst_F = worst_H = 0.0
    for _ in range(5):
        z = rng.normal(size=3) * np.array([60, 3, 25]) + 1j * rng.normal(size=3) * np.array([60, 3, 25])
        V, F = pot_force(*z); Vr, Fr = rp.potential(z)
        worst_V = max(worst_V, abs(V - Vr) / max(1, abs(Vr)))
        worst_F = max(worst_F, float(np.max(np.abs(np.array(F) - Fr)) / max(1, np.max(np.abs(Fr)))))
        # Hessian vs central differences of the gradient grad_r V = sqrt2 (Re F, Im F) built from rp.potential
        xiR = -7.3
        Hm = hessian_r(*z, xiR)

        def grad(r):
            zz = (r[0::2] + 1j * r[1::2]) / SQ2
            _, FF = rp.potential(zz)
            FF = np.array(FF, complex); FF[0] = FF[0] + xiR * zz[0]
            out = np.empty(6); out[0::2] = SQ2 * FF.real; out[1::2] = SQ2 * FF.imag
            return out
        r0 = np.empty(6); r0[0::2] = SQ2 * z.real; r0[1::2] = SQ2 * z.imag
        Hfd = np.empty((6, 6))
        for i in range(6):
            h = 1e-3 * max(1.0, abs(r0[i]))
            e1 = np.zeros(6); e1[i] = h
            d1 = (grad(r0 + e1) - grad(r0 - e1)) / (2 * h)
            e2 = np.zeros(6); e2[i] = 2 * h
            d2 = (grad(r0 + e2) - grad(r0 - e2)) / (4 * h)
            Hfd[:, i] = (4 * d1 - d2) / 3                     # Richardson: exact for the cubic gradient up to roundoff
        worst_H = max(worst_H, float(np.max(np.abs(Hm - Hfd)) / np.max(np.abs(Hm))))
    return dict(potential_rel=worst_V, force_rel=worst_F, hessian_vs_richardson_fd_rel=worst_H)


if __name__ == '__main__':
    print(json.dumps(selftest(), indent=1))
