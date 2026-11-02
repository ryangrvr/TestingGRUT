#!/usr/bin/env python3
"""fx_holder_skeptic (Claude, adversarial reviewer of track fx_holder; own code written for this review).

Shared model for the skeptic's independent recomputation.  Nothing here imports or runs the track's fx_common.py,
rptri.c or tri_point.py, nor any Codex code.  Pins are the M1 GeV values [Q d4747cb], re-entered here.

Model: V = lx|X|^4 + ly|Y|^4 + M^2|Y|^2 + ls(|S|^2 - h^2/2)^2 + 2Re[kap Y* X^3 + u S* Y^2] + mc |X|^2
       (+ optional phase-blind counterterm  cport (|u|^2 - mu^2)(|S|^2 + 2|Y|^2), see sk_evolve.py)
with u = mu + (holder), mc = xiX R.  Wirtinger derivatives derived by hand for this review:
  dV/dX* = 2 lx|X|^2 X + 3 kap Y conj(X)^2 + mc X
  dV/dY* = (M^2 + 2 ly|Y|^2) Y + kap X^3 + 2 conj(u) S conj(Y)
  dV/dS* = 2 ls(|S|^2 - h^2/2) S + u Y^2
Units: H_i = 1.
"""
import numpy as np
import mpmath as mp

H_GEV = 134284273491.30815
MR_GEV = 2.435e18
M = 4028528204739.244 / H_GEV
HV = 5371370939652.326 / H_GEV
MU = 961402368973.181 / H_GEV
KAP = 0.0034586847211042523
LX = 1e-3
XIX = -5.0 / 6.0
LOCK = (10j - MU) / 12.0
YQ = 5.07941423e-9
GSTAR = 106.75


def par(lam=1.0, cport=0.0):
    return dict(lx=LX, ly=lam, ls=lam, kap=KAP, M=M, h=HV, mu=MU, cport=cport)


def grad_c(z, u, p, mc):
    """complex gradient (dV/dX*, dV/dY*, dV/dS*) incl. optional phase-blind counterterm."""
    x, y, s = z
    x2, y2, s2 = abs(x) ** 2, abs(y) ** 2, abs(s) ** 2
    du = p['cport'] * (abs(u) ** 2 - p['mu'] ** 2)
    gx = 2 * p['lx'] * x2 * x + 3 * p['kap'] * y * np.conj(x) ** 2 + mc * x
    gy = (p['M'] ** 2 + 2 * p['ly'] * y2 + 2 * du) * y + p['kap'] * x ** 3 + 2 * np.conj(u) * s * np.conj(y)
    gs = (2 * p['ls'] * (s2 - p['h'] ** 2 / 2) + du) * s + u * y * y
    return np.array([gx, gy, gs])


def V_c(z, u, p, mc):
    x, y, s = z
    x2, y2, s2 = abs(x) ** 2, abs(y) ** 2, abs(s) ** 2
    du = p['cport'] * (abs(u) ** 2 - p['mu'] ** 2)
    return float(p['lx'] * x2 ** 2 + p['ly'] * y2 ** 2 + p['M'] ** 2 * y2 + p['ls'] * (s2 - p['h'] ** 2 / 2) ** 2
                 + 2 * np.real(p['kap'] * np.conj(y) * x ** 3 + u * np.conj(s) * y * y) + mc * x2 + du * (s2 + 2 * y2))


# ---------------------------------------------------------------- mpmath versions (40 digits)
def mp_V(r, u, p, mc):
    x = mp.mpc(r[0], r[1]) / mp.sqrt(2); y = mp.mpc(r[2], r[3]) / mp.sqrt(2); s = mp.mpc(r[4], r[5]) / mp.sqrt(2)
    x2, y2, s2 = abs(x) ** 2, abs(y) ** 2, abs(s) ** 2
    lx, ly, ls, kap, MM, h = [mp.mpf(p[k]) for k in ('lx', 'ly', 'ls', 'kap', 'M', 'h')]
    du = mp.mpf(p['cport']) * (abs(u) ** 2 - mp.mpf(p['mu']) ** 2)
    return (lx * x2 ** 2 + ly * y2 ** 2 + MM ** 2 * y2 + ls * (s2 - h ** 2 / 2) ** 2
            + 2 * mp.re(kap * mp.conj(y) * x ** 3 + u * mp.conj(s) * y * y) + mc * x2 + du * (s2 + 2 * y2))


def mp_grad_c(z, u, p, mc):
    x, y, s = z
    lx, ly, ls, kap, MM, h = [mp.mpf(p[k]) for k in ('lx', 'ly', 'ls', 'kap', 'M', 'h')]
    du = mp.mpf(p['cport']) * (abs(u) ** 2 - mp.mpf(p['mu']) ** 2)
    gx = 2 * lx * abs(x) ** 2 * x + 3 * kap * y * mp.conj(x) ** 2 + mc * x
    gy = (MM ** 2 + 2 * ly * abs(y) ** 2 + 2 * du) * y + kap * x ** 3 + 2 * mp.conj(u) * s * mp.conj(y)
    gs = (2 * ls * (abs(s) ** 2 - h ** 2 / 2) + du) * s + u * y * y
    return gx, gy, gs


def mp_stationary(u, p, mc, guess):
    """true stationary point with X real (gauge), unknowns (x, Re Y, Im Y, Re S, Im S); 5 equations
    Re gX, Re gY, Im gY, Re gS, Im gS (Im gX follows from the U(1) identity)."""
    u = mp.mpc(u)
    def eqs(x, yr, yi, sr, si):
        gx, gy, gs = mp_grad_c((mp.mpc(x, 0), mp.mpc(yr, yi), mp.mpc(sr, si)), u, p, mc)
        return [mp.re(gx), mp.re(gy), mp.im(gy), mp.re(gs), mp.im(gs)]
    g = [mp.mpf(float(guess[0].real)), mp.mpf(float(guess[1].real)), mp.mpf(float(guess[1].imag)),
         mp.mpf(float(guess[2].real)), mp.mpf(float(guess[2].imag))]
    sol = mp.findroot(eqs, g, tol=mp.mpf(10) ** (-30), maxsteps=200)
    z = (mp.mpc(sol[0], 0), mp.mpc(sol[1], sol[2]), mp.mpc(sol[3], sol[4]))
    gx, gy, gs = mp_grad_c(z, u, p, mc)
    res = max(abs(gx), abs(gy), abs(gs))
    return z, res


def mp_held(u, p, mc, guess):
    """held-phase preparation (the 'original' of the M1 chain): X = x > 0, S = -i s (phase held), Y free; equations
    Re gX = 0, gY = 0, and the radial S equation Re(gS * conj(-i)) = 0.  Written for this review from the description
    'conditional minimum at held phases' [Q REVIEW_REHEATING]."""
    u = mp.mpc(u)
    d = mp.mpc(0, -1)
    def eqs(x, yr, yi, s):
        gx, gy, gs = mp_grad_c((mp.mpc(x, 0), mp.mpc(yr, yi), d * s), u, p, mc)
        return [mp.re(gx), mp.re(gy), mp.im(gy), mp.re(gs * mp.conj(d))]
    sol = mp.findroot(eqs, [mp.mpf(v) for v in guess], tol=mp.mpf(10) ** (-30), maxsteps=200)
    return (mp.mpc(sol[0], 0), mp.mpc(sol[1], sol[2]), d * sol[3])


def mp_hessian(z, u, p, mc, dps=30):
    """6x6 Hessian in canonical real fields r = sqrt2 (Re, Im) by mpmath numerical differentiation."""
    r0 = [mp.sqrt(2) * c for zz in z for c in (mp.re(zz), mp.im(zz))]
    Hm = mp.matrix(6, 6)
    for i in range(6):
        for j in range(i, 6):
            def f(a, b, i=i, j=j):
                rr = list(r0); rr[i] += a; rr[j] += b
                return mp_V(rr, mp.mpc(u), p, mc)
            if i == j:
                def g1(a, i=i):
                    rr = list(r0); rr[i] += a
                    return mp_V(rr, mp.mpc(u), p, mc)
                Hm[i, i] = mp.diff(g1, 0, 2)
            else:
                Hm[i, j] = Hm[j, i] = mp.diff(f, (0, 0), (1, 1))
    E, Q = mp.eigsy(Hm)
    return Hm, sorted([E[k] for k in range(6)]), Q


def phases(z):
    x, y, s = [complex(c) for c in z]
    th = np.angle(x)
    return float(np.angle(s * np.exp(-6j * th))), float(np.angle(y * np.exp(-3j * th)))
