#!/usr/bin/env python3
"""mf_chain SKEPTIC (key mf_chain_skeptic; Claude, multi-field inflaton pass, GitHub issue #13 context).

INDEPENDENT re-implementation of the homogeneous M1r release on an external (bench-format) background with the holder on a
spectator s.  Nothing of the track's engine (rpmfi.c) is used: the equations are re-derived from the M1r potential
(BRIEF_DECISIVE.md), the integrator is scipy's DOP853 (not the engine's hand-coded Dormand-Prince 5(4) with a step cap),
the background is interpolated with not-a-knot cubic splines (H, R, m_s^2; ln a = exact spline integral of H) instead of
cubic Hermite with finite-difference slopes, the preparation is my own Newton/least-squares solve of the stationary point,
and N_YL, N_YH are read with my own forced-point subtraction (finite-difference forced velocity).

Model [Q BRIEF_DECISIVE, lambda = 0.5]:
  V = lx|X|^4 + M^2|Y|^2 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + 2 Re[kap conj(Y) X^3 + u conj(S) Y^2],  + xi R |X|^2
  u = mu + G c,  c = s/s_*,  G = 10i - mu  (mu_eff(slice) = 10i at c = 1: the M1r calibration)
  z'' + 3H z' + dV/dconj(z) = 0 (complex fields, L = |z'|^2 - V);  c'' + (3H + Gs) c' + m_s^2(t) c + (H_i/s_*)^2 dV/dc = 0.
Units H_i = 1, a = 1 at the slice (t = 0).  Every history and parameter is a model choice; nothing follows from closure.
Labels: [C] given the model.
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

LX, LY, LS = 1e-3, 0.5, 0.5
KAP = 0.0034586847211042523
MM = 30.0
HH = 40.0
MU = 961402368973.181 / 134284273491.30815          # 7.159456...
XI = -5.0 / 6.0
G_DEF = 10j - MU
QCH = np.array([1.0, 3.0, 6.0])


def Fvec(x, y, s, u, R):
    """dV/dconj(z) for (X, Y, Sigma), including xi R X."""
    x2 = (x * x.conjugate()).real; y2 = (y * y.conjugate()).real; s2 = (s * s.conjugate()).real
    Fx = 2 * LX * x2 * x + 3 * KAP * y * x.conjugate() ** 2 + XI * R * x
    Fy = MM * MM * y + 2 * LY * y2 * y + KAP * x ** 3 + 2 * u.conjugate() * s * y.conjugate()
    Fs = 2 * LS * (s2 - HH * HH / 2) * s + u * y * y
    return Fx, Fy, Fs


def Vpot(x, y, s, u, R=0.0):
    x2 = abs(x) ** 2; y2 = abs(y) ** 2; s2 = abs(s) ** 2
    return (LX * x2 * x2 + MM ** 2 * y2 + LY * y2 * y2 + LS * (s2 - HH ** 2 / 2) ** 2
            + 2 * (KAP * y.conjugate() * x ** 3 + u * s.conjugate() * y * y).real + XI * R * x2)


_PREP_CACHE = {}


def prep_cont(c0, G=G_DEF, R=12.0):
    """prepared point at u = mu + G c0 by continuation in c0 from the M1r point (c0 = 1), steps of <= 0.02."""
    key = (round(c0, 12), complex(G), R)
    if key in _PREP_CACHE:
        return _PREP_CACHE[key]
    n = max(1, int(math.ceil(abs(c0 - 1.0) / 0.02)))
    guess = (90.8, -6.888, -28.576j)
    for c in np.linspace(1.0, c0, n + 1):
        z, gres = prep(MU + G * c, R, guess=guess)
        guess = z
    _PREP_CACHE[key] = (z, gres)
    return z, gres


def prep(u, R=12.0, guess=(90.8, -6.888, -28.576j)):
    """stationary point of V + xi R|X|^2 at fixed u, gauge X real > 0 (my own least-squares on the analytic gradient)."""
    def g(r):
        x = r[0] + 0j; y = r[1] + 1j * r[2]; s = r[3] + 1j * r[4]
        Fx, Fy, Fs = Fvec(x, y, s, u, R)
        return [2 * Fx.real, 2 * Fy.real, 2 * Fy.imag, 2 * Fs.real, 2 * Fs.imag, 2 * Fx.imag]
    x0, y0, s0 = guess
    r0 = [abs(x0), complex(y0).real, complex(y0).imag, complex(s0).real, complex(s0).imag]
    sol = least_squares(lambda r: g(r), r0, xtol=1e-15, ftol=1e-15, gtol=1e-15, method='lm', max_nfev=20000)
    r = sol.x
    z = (r[0] + 0j, r[1] + 1j * r[2], r[3] + 1j * r[4])
    return z, float(np.linalg.norm(g(r)))


class PW:
    """two-piece not-a-knot cubic spline joined at t_j (the bench's switch from the resolved homogeneous stage to the
    sparsely sampled fluid): avoids ringing of one global spline across the change of sampling (Delta t 4e-4 -> 0.06)."""
    def __init__(self, t, f, tj):
        i = int(np.searchsorted(t, tj, side='right'))          # rows t <= tj
        self.tj = t[i - 1]
        self.a = CubicSpline(t[:i], f[:i])
        self.b = CubicSpline(t[i - 1:], f[i - 1:]) if len(t) - i + 1 >= 2 else None
        if self.b is not None and len(t) - i + 1 < 4:
            from scipy.interpolate import interp1d
            self.b = interp1d(t[i - 1:], f[i - 1:])
        self.ia = self.a.antiderivative()
        self.ib = self.b.antiderivative() if hasattr(self.b, 'antiderivative') else None
        self.Ia_j = float(self.ia(self.tj))

    def __call__(self, t):
        return float(self.a(t)) if (t <= self.tj or self.b is None) else float(self.b(t))

    def integral(self, t):
        if t <= self.tj or self.ib is None:
            return float(self.ia(t))
        return self.Ia_j + float(self.ib(t) - self.ib(self.tj))


class Background:
    """bench-format history -> two-piece splines (t >= 0); exact radiation (a^2 H fixed, R = 0) beyond the last row."""
    def __init__(self, bench, R_mode='raw', ms_scale=1.0, ms_override=None, t_join=None):
        h = bench['history']
        t = np.array(h['t'], float); N = np.array(h['N'], float); H = np.array(h['H'], float); R = np.array(h['R'], float)
        ms = np.array(h['m_s'], float) * ms_scale
        N0 = float(np.interp(0.0, t, N)); H0 = float(np.interp(0.0, t, H)); R0 = float(np.interp(0.0, t, R))
        m0 = float(np.interp(0.0, t, ms))
        sel = t > 1e-12
        t = np.r_[0.0, t[sel]]; N = np.r_[N0, N[sel]] - N0; H = np.r_[H0, H[sel]]; R = np.r_[R0, R[sel]]; ms = np.r_[m0, ms[sel]]
        keep = np.r_[True, np.diff(t) > 1e-13]
        t, N, H, R, ms = t[keep], N[keep], H[keep], R[keep], ms[keep]
        if t_join is None:
            ro = h.get('R_osc_sig')
            t_join = float(max(tt for tt, r in zip(h['t'], ro) if r is None)) if ro is not None else float(t[-1])
        self.t, self.N, self.Hs, self.Rs = t, N, H, R
        self.sH = PW(t, H, t_join); self.sR = PW(t, R, t_join)
        m2 = ms * ms if ms_override is None else ms_override(t, ms)
        self.sM2 = PW(t, m2, t_join)
        self.te = float(t[-1]); self.Ne = self.sH.integral(self.te); self.He = float(H[-1])
        self.lnaN_dev = float(np.max(np.abs(np.array([self.sH.integral(x) for x in t[::7]]) - N[::7])))
        self.m2e = float(m2[-1])

    def at(self, t):
        if t <= self.te:
            return self.sH.integral(t), self.sH(t), self.sR(t), self.sM2(t)
        x = 1 + 2 * self.He * (t - self.te)
        return self.Ne + 0.5 * math.log(x), self.He / x, 0.0, self.m2e


def run(bench_path, T=100.0, dt_out=0.005, c0=1.0, Gs=None, sstar_GeV=None, R_mode='raw', ms_scale=1.0, rtol=1e-10, atol=1e-12,
        hom_cap=4e-4, late_cap=0.02, ms_override=None, hold_s=False, tag=None, G=G_DEF):
    bench = json.load(open(bench_path))
    bg = Background(bench, R_mode, ms_scale, ms_override)
    HiG = float(bench['params'].get('H_i_GeV', bench.get('H_i_GeV')))
    sst = float(bench['params'].get('s_star_GeV')) if sstar_GeV is None else sstar_GeV
    gs = float(bench.get('spectator', {}).get('Gamma_s_assumed', 0.0) or 0.0) if Gs is None else Gs
    BR = (HiG / sst) ** 2
    z0, gres = prep_cont(c0, G)
    # end of the resolved (homogeneous) part: where R_osc columns start
    ro = bench['history'].get('R_osc_sig')
    th = float(max(tt for tt, r in zip(bench['history']['t'], ro) if r is None)) if ro is not None else 0.6

    def rhs(t, y):
        la, H, R, m2 = bg.at(t)
        x = y[0] + 1j * y[3]; yy = y[1] + 1j * y[4]; s = y[2] + 1j * y[5]
        vx = y[6] + 1j * y[9]; vy = y[7] + 1j * y[10]; vs = y[8] + 1j * y[11]
        c, cd = y[12], y[13]
        u = MU + G * c
        Fx, Fy, Fs = Fvec(x, yy, s, u, R)
        ax = -3 * H * vx - Fx; ay = -3 * H * vy - Fy; as_ = -3 * H * vs - Fs
        dVdc = 2 * (G * s.conjugate() * yy * yy).real
        cdd = 0.0 if hold_s else (-(3 * H + gs) * cd - m2 * c - BR * dVdc)
        return [vx.real, vy.real, vs.real, vx.imag, vy.imag, vs.imag, ax.real, ay.real, as_.real, ax.imag, ay.imag, as_.imag,
                0.0 if hold_s else cd, cdd]

    y0 = [z0[0].real, z0[1].real, z0[2].real, z0[0].imag, z0[1].imag, z0[2].imag] + [0.0] * 6 + [c0, 0.0]
    tg = np.arange(0.0, T + 1e-9, dt_out)
    t0 = time.time()
    s1 = solve_ivp(rhs, (0.0, th), y0, method='DOP853', rtol=rtol, atol=atol, max_step=hom_cap, t_eval=tg[tg <= th], dense_output=True)
    y1 = s1.sol(th)
    s2 = solve_ivp(rhs, (th, T), y1, method='DOP853', rtol=rtol, atol=atol, max_step=late_cap, t_eval=tg[tg > th])
    Y = np.concatenate([s1.y, s2.y], axis=1); tt = np.concatenate([s1.t, s2.t])
    secs = time.time() - t0
    D = diagnostics(tt, Y, bg, G)
    late = tt >= 0.75 * T
    out = dict(tag=tag, bench=str(bench_path), T=T, c0=c0, Gamma_s=gs, s_star_GeV=sst, BR=BR, ms_scale=ms_scale, hold_s=hold_s,
               prep=dict(z=[[c.real, c.imag] for c in z0], amps=[abs(c) for c in z0], grad=gres), t_hom=th, seconds=secs,
               nfev=int(s1.nfev + s2.nfev), lnaN_dev=bg.lnaN_dev, rtol=rtol,
               late_mean={k: float(np.mean(v[late])) for k, v in D.items() if k not in ('t',)},
               at_T={k: float(v[-1]) for k, v in D.items()},
               gauss_max=float(np.max(np.abs(D['gauss']))),
               timeline=[{k: float('%.6g' % D[k][i]) for k in ('t', 'a', 'NX', 'EX', 'NYL', 'NYH', 'Nsig', 'c')}
                         for i in np.unique(np.r_[np.arange(0, min(len(tt), 400), 20), np.geomspace(1, len(tt) - 1, 50).astype(int)])])
    return out, (tt, Y, D)


def diagnostics(tt, Y, bg, G):
    la = np.array([bg.at(t)[0] for t in tt]); a = np.exp(la); a3 = a ** 3
    x = Y[0] + 1j * Y[3]; yy = Y[1] + 1j * Y[4]; s = Y[2] + 1j * Y[5]
    vx = Y[6] + 1j * Y[9]; vy = Y[7] + 1j * Y[10]; vs = Y[8] + 1j * Y[11]
    c = Y[12]
    n = [-2 * np.imag(np.conj(z) * v) * a3 for z, v in ((x, vx), (yy, vy), (s, vs))]
    EX = a ** 4 * (np.abs(vx) ** 2 + LX * np.abs(x) ** 4)
    Nc = 2 * math.sqrt(2 * LX) * (EX / (3 * LX)) ** 0.75
    eps = n[0] / Nc
    u = MU + G * c
    # forced point of the linear Y equation  A Y + B conj(Y) = -J,  A = M^2, B = 2 conj(u) Sigma, J = kap X^3
    A = MM * MM; B = 2 * np.conj(u) * s; J = KAP * x ** 3
    Dt = A * A - np.abs(B) ** 2
    yp = (-A * J + B * np.conj(J)) / Dt
    # one fixed-point pass with the Y quartic (|yp| << 1 late)
    Ae = A + 2 * LY * np.abs(yp) ** 2
    yp = (-Ae * J + B * np.conj(J)) / (Ae * Ae - np.abs(B) ** 2)
    ypd = np.gradient(yp, tt)
    th = np.angle(B)
    rot = np.exp(-0.5j * th)
    thd = np.gradient(np.unwrap(th), tt)
    dy = (yy - yp) * rot
    dv = (vy - ypd - 0.5j * thd * (yy - yp)) * rot
    A2 = A + 4 * LY * np.abs(yp) ** 2
    bb = np.abs(B + 2 * LY * yp * yp)
    wH = np.sqrt(A2 + bb); wL = np.sqrt(np.maximum(A2 - bb, 1e-300))
    NYH = a3 * (dv.real ** 2 + wH ** 2 * dy.real ** 2) / wH
    NYL = a3 * (dv.imag ** 2 + wL ** 2 * dy.imag ** 2) / wL
    sr = np.abs(s); sd = np.real(np.conj(s) * vs) / sr
    Nsig = a3 * (sd ** 2 + LS * (sr ** 2 - HH ** 2 / 2) ** 2) / (math.sqrt(2 * LS) * HH)
    scale = np.maximum(np.abs(n[0]) + 3 * np.abs(n[1]) + 6 * np.abs(n[2]), 1.0)
    gauss = (n[0] + 3 * n[1] + 6 * n[2]) / scale
    return dict(t=tt, a=a, NX=n[0], NY=n[1], NS=n[2], EX=EX, eps=eps, NYL=NYL, NYH=NYH, Nsig=Nsig, c=c, gauss=gauss,
                absX=np.abs(x), mL=wL)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('bench'); ap.add_argument('tag'); ap.add_argument('--T', type=float, default=100.0)
    ap.add_argument('--c0', type=float, default=1.0); ap.add_argument('--ms_scale', type=float, default=1.0)
    ap.add_argument('--rtol', type=float, default=1e-10); ap.add_argument('--hold_s', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    o, _ = run(a.bench, T=a.T, c0=a.c0, ms_scale=a.ms_scale, rtol=a.rtol, hold_s=a.hold_s, tag=a.tag)
    lm = o['late_mean']
    print('%s: prep %s grad %.1e | late(%.0f-%.0f): NX %.6g EX %.6g eps %.5f NYL %.5g NYH %.4g Nsig %.4g | gauss %.1e | %.1f s, %d fev' % (
        a.tag, ['%.6g' % x for x in o['prep']['amps']], o['prep']['grad'], 0.75 * a.T, a.T, lm['NX'], lm['EX'], lm['eps'], lm['NYL'], lm['NYH'],
        lm['Nsig'], o['gauss_max'], o['seconds'], o['nfev']), flush=True)
    if a.out:
        json.dump(o, open(a.out, 'w'), indent=1)


def prep_hessian(c0=1.0, R=12.0, G=G_DEF, h=1e-5):
    """eigenvalues of the 6x6 real Hessian of V + xi R|X|^2 at the prepared point (kinetic |z'|^2 -> mass^2 = eig/2)."""
    u = MU + G * c0
    z, gres = prep_cont(c0, G, R)
    r0 = np.array([z[0].real, z[0].imag, z[1].real, z[1].imag, z[2].real, z[2].imag])

    def V(r):
        return Vpot(r[0] + 1j * r[1], r[2] + 1j * r[3], r[4] + 1j * r[5], u, R)
    Hm = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            e_i = np.zeros(6); e_j = np.zeros(6); e_i[i] = h; e_j[j] = h
            Hm[i, j] = (V(r0 + e_i + e_j) - V(r0 + e_i - e_j) - V(r0 - e_i + e_j) + V(r0 - e_i - e_j)) / (4 * h * h)
    ev = np.linalg.eigvalsh(0.5 * (Hm + Hm.T)) / 2.0
    return dict(c0=c0, u=[u.real, u.imag], amps=[abs(x) for x in z], grad=gres, mass2=sorted(ev.tolist()))
