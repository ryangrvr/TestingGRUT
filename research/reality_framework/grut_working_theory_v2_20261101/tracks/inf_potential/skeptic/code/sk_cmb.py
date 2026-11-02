#!/usr/bin/env python3
"""inf_potential SKEPTIC (Claude, own independent code): recompute the CMB numbers of the three benchmarks.

Independent of inf_common.py: background in COSMIC time t (not e-folds), scalar modes as the flat-gauge field
perturbation Q (Sasaki-Mukhanov variable), Q'' + 3H Q' + [k^2/a^2 + V'' - a^-3 d/dt(a^3 phidot^2/H)] Q = 0,
R = H Q/phidot; tensors u'' + 3H u' + k^2/a^2 u = 0 with P_t = 2 * 4 * k^3/(2 pi^2) |u|^2 (u = v_t/a, BD).
Bunch-Davies start at k/aH = 250, read-out at k/aH = 2e-3.  Units M_P = 1; V = V0 v(phi); V0 is fixed by A_s at the pivot.
N_* (pivot exit before eps_H = 1) is recomputed from our own reheating model: matter-like inflaton (cycle averaged)
with constant width Gamma, radiation g* = 106.75, entropy conservation to T0 = 2.7255 K, g_s0 = 3.909
(and 3.938 as a variant), and also instant reheating.  The pivot is k_* = 0.05 Mpc^-1.
Only the benchmark SHAPE parameters are taken from the track's bench_*.json (read-only).
Output: ../results/sk_cmb.json
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]                     # tracks/inf_potential
OUT = HERE.parent / 'results'
MP = 2.435e18
GEV_MPC = 1.0 / (1.97327e-16 / 3.08568e22)  # 1 GeV in Mpc^-1
KP = 0.05
T0 = 2.7255 * 8.617333e-14
AS = 2.1e-9


class Pot:
    def __init__(self, name, par):
        self.name = name
        if name in ('T1',):
            self.M = math.sqrt(6 * par['alpha']); self.k = 'T'
        elif name in ('E1',):
            self.M = math.sqrt(1.5 * par['alpha']); self.k = 'E'
        else:
            self.k = 'HT'; self.vv = par['vev']; self.p = par['p']; self.b = par.get('beta', 0.0)
            xt = -(self.b / self.p) ** (1 / (self.p - 1)) if self.b > 0 else 0.0
            self.ft = 1 - self.b * xt - xt ** self.p
            self.xmin = brentq(lambda x: 1 - self.b * x - x ** self.p, 0.5, 1.0)

    def V(self, f):
        if self.k == 'T':
            return np.tanh(f / self.M) ** 2
        if self.k == 'E':
            return (1 - np.exp(-f / self.M)) ** 2
        x = f / self.vv
        return (1 - self.b * x - x ** self.p) ** 2 / self.ft ** 2

    def dV(self, f, h=None):
        if self.k == 'T':
            y = f / self.M; return 2 * np.tanh(y) / np.cosh(y) ** 2 / self.M
        if self.k == 'E':
            e = np.exp(-f / self.M); return 2 * (1 - e) * e / self.M
        x = f / self.vv; g = 1 - self.b * x - x ** self.p; gp = (-self.b - self.p * x ** (self.p - 1)) / self.vv
        return 2 * g * gp / self.ft ** 2

    def d2V(self, f):
        if self.k == 'T':
            y = f / self.M; c = np.cosh(y); return 2 * (1 - 2 * np.sinh(y) ** 2) / c ** 4 / self.M ** 2
        if self.k == 'E':
            e = np.exp(-f / self.M); return 2 * (2 * e * e - e) / self.M ** 2
        x = f / self.vv; g = 1 - self.b * x - x ** self.p; gp = (-self.b - self.p * x ** (self.p - 1)) / self.vv
        gpp = -self.p * (self.p - 1) * x ** (self.p - 2) / self.vv ** 2
        return 2 * (gp * gp + g * gpp) / self.ft ** 2

    def phimin(self):
        return 0.0 if self.k in ('T', 'E') else self.xmin * self.vv


def bg_rhs(t, y, P):
    f, fd, lna = y
    H = math.sqrt((0.5 * fd * fd + P.V(f)) / 3)
    return [fd, -3 * H * fd - P.dV(f), H]


def background(P, f0):
    """slow-roll start at f0, cosmic time; returns dense solution and t_end (eps_H = 1)."""
    H0 = math.sqrt(P.V(f0) / 3); fd0 = -P.dV(f0) / (3 * H0)
    ev = lambda t, y, P: 0.5 * y[1] ** 2 / ((0.5 * y[1] ** 2 + P.V(y[0])) / 3) - 1.0   # eps_H - 1
    ev.terminal = True; ev.direction = 1
    s = solve_ivp(bg_rhs, (0, 1e9), [f0, fd0, 0.0], args=(P,), method='DOP853', rtol=1e-12, atol=1e-14,
                  dense_output=True, events=ev)
    return s, s.t_events[0][0]


def f_start(P, Ntot=70.0):
    """field value Ntot slow-roll e-folds before the end (rough, then the exact background is used)."""
    if P.k == 'T' or P.k == 'E':
        sgn = 1.0
        g = lambda f: np.sqrt(2) * abs(P.dV(f) / P.V(f)) - 1
        fe = brentq(g, 1e-6 * P.M, 30 * P.M)
        fs = np.linspace(fe, 60 * P.M, 200001)
        integ = P.V(fs) / P.dV(fs)
        N = np.concatenate([[0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(fs))])
        return float(np.interp(Ntot, N, fs))
    # hilltop: field rolls from near the top (x_top < 0) toward x_min ~ 1
    xt = -(P.b / P.p) ** (1 / (P.p - 1)) if P.b > 0 else 0.0
    fe = brentq(lambda f: np.sqrt(2) * abs(P.dV(f) / P.V(f)) - 1, (xt + 1e-6) * P.vv + 1e-9, P.xmin * P.vv * 0.999)
    fs = np.linspace(fe, xt * P.vv + 1e-7, 400001)
    integ = P.V(fs) / P.dV(fs)
    N = np.concatenate([[0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(fs))])
    N = np.abs(N)
    return float(np.interp(Ntot, N, fs))


def mode(P, sol, k, t_lo, t_hi, kind):
    """integrate one mode from k/aH = 250 to 2e-3; kind 's' (Q) or 't' (u)."""
    def aH(t):
        f, fd, lna = sol.sol(t); H = math.sqrt((0.5 * fd * fd + P.V(f)) / 3); return math.exp(lna) * H
    ta = brentq(lambda t: math.log(k / aH(t)) - math.log(250.0), t_lo, t_hi)
    tb = brentq(lambda t: math.log(k / aH(t)) - math.log(2e-3), t_lo, t_hi)
    f, fd, lna = sol.sol(ta); a = math.exp(lna); H = math.sqrt((0.5 * fd * fd + P.V(f)) / 3)
    Q0 = 1 / (a * math.sqrt(2 * k)); Qd0 = -(1j * k / a + H) * Q0

    def rhs(t, y):
        f, fd, lna = y[0], y[1], y[2]
        H = math.sqrt((0.5 * fd * fd + P.V(f)) / 3); a = math.exp(lna)
        fdd = -3 * H * fd - P.dV(f)
        eps = 0.5 * fd * fd / H ** 2
        if kind == 's':
            meff = P.d2V(f) - (3 * fd * fd + 2 * fd * fdd / H + fd * fd * eps)   # a^-3 d/dt(a^3 fd^2/H)
        else:
            meff = 0.0
        Q = y[3] + 1j * y[4]; Qd = y[5] + 1j * y[6]
        Qdd = -3 * H * Qd - (k * k / (a * a) + meff) * Q
        return [fd, fdd, H, Qd.real, Qd.imag, Qdd.real, Qdd.imag]
    y0 = [f, fd, lna, Q0.real, Q0.imag, Qd0.real, Qd0.imag]
    s = solve_ivp(rhs, (ta, tb), y0, method='DOP853', rtol=1e-10, atol=1e-30)
    f, fd, lna = s.y[0, -1], s.y[1, -1], s.y[2, -1]
    H = math.sqrt((0.5 * fd * fd + P.V(f)) / 3)
    Q = s.y[3, -1] + 1j * s.y[4, -1]
    if kind == 's':
        return k ** 3 / (2 * math.pi ** 2) * abs(H * Q / fd) ** 2
    return 8 * k ** 3 / (2 * math.pi ** 2) * abs(Q) ** 2


def reheat_dN(gam_over_Hend, Hend_GeV, gs0=3.909, grh=106.75):
    """a_end/a_0 for: matter-like inflaton (rho_phi, cycle averaged) decaying at constant Gamma into radiation, started at
    eps_H = 1 with rho = 3 H_end^2 (all in the oscillating field, a deliberately simple model), integrated until
    rho_phi/rho < 1e-8, then entropy conservation.  Returns ln(a_end/a_0)."""
    if gam_over_Hend is None:
        rho = 3 * (Hend_GeV) ** 2 * MP ** 2
        T = (rho * 30 / (math.pi ** 2 * grh)) ** 0.25
        return math.log(T0 * (gs0 / grh) ** (1 / 3) / T)
    G = gam_over_Hend

    def rhs(N, y):
        rp, rr = math.exp(y[0]), math.exp(y[1])
        H = math.sqrt((rp + rr) / 3)
        return [-3 - G / H, -4 + G * rp / (H * rr)]
    ev = lambda N, y: y[0] - y[1] - math.log(1e-8); ev.terminal = True
    s = solve_ivp(rhs, (0, 60), [math.log(3.0), math.log(3e-30)], method='LSODA', rtol=1e-10, atol=1e-12, events=ev)
    N = s.t[-1]; rr = math.exp(s.y[1, -1])                 # in units H_end^2 M_P^2
    rho = rr * Hend_GeV ** 2 * MP ** 2
    T = (rho * 30 / (math.pi ** 2 * grh)) ** 0.25
    return math.log(T0 * (gs0 / grh) ** (1 / 3) / T) - N   # ln(a_end/a_0) = ln(a_th/a_0) - N


def analyse(name, Hi=1.29e11, Nstar_track=None, gam_Hi=0.5):
    d = json.load(open(TRK / 'results' / ('bench_%s.json' % name)))
    rec = d['by_Hi'][0]
    P = Pot(name, rec['params'])
    f0 = f_start(P, 72.0)
    sol, tend = background(P, f0)
    # e-folds before end, H, eps on a grid
    ts = np.linspace(0, tend, 400001)
    Y = sol.sol(ts); lna = Y[2]; Nb = lna[-1] - lna
    Hs = np.sqrt((0.5 * Y[1] ** 2 + P.V(Y[0])) / 3)
    t_of_N = lambda N: float(np.interp(-N, -Nb, ts))
    H_of_N = lambda N: float(np.interp(-N, -Nb, Hs))
    tS1 = t_of_N(1.0)
    # V0 from H(S1) = H_i (the track's calibration convention) -> check A_s afterwards
    V0 = (Hi / MP / H_of_N(1.0)) ** 2
    Hend = math.sqrt(V0) * Hs[-1] * MP
    res = dict(benchmark=name, params=rec['params'], V0_MP4_skeptic=V0, V0_MP4_track=rec['V0_MP4'], Hend_GeV=Hend,
               Hend_over_Hi=Hend / Hi, phi_S1=float(sol.sol(tS1)[0]), phi_end=float(Y[0, -1]))
    # N_* by our own matching
    out_N = {}
    for lab, G in (('instant', None), ('Gamma0.5Hi', gam_Hi * Hi / Hend)):
        lnae = reheat_dN(G, Hend)
        # N_*: k_* = a_* H_*  with a_* = a_end e^{-N}
        g = lambda N: math.log(KP / GEV_MPC) - (lnae - N + math.log(math.sqrt(V0) * H_of_N(N) * MP))
        out_N[lab] = brentq(g, 40, 65)
    lnae = reheat_dN(gam_Hi * Hi / Hend, Hend, gs0=3.938)
    out_N['Gamma0.5Hi_gs3.938'] = brentq(lambda N: math.log(KP / GEV_MPC) - (lnae - N + math.log(math.sqrt(V0) * H_of_N(N) * MP)), 40, 65)
    res['Nstar_skeptic'] = out_N
    res['Nstar_track'] = rec['Nstar']
    Ns = out_N['Gamma0.5Hi']
    # spectra around the pivot (k in code units, a(t=0) = 1): k_* = a H at Ns before the end
    lnk = []
    Ps = []; Pt = []
    for dN in (-0.6, -0.3, 0.0, 0.3, 0.6):
        N = Ns + dN
        aH = math.exp(lna[-1] - N) * H_of_N(N)
        k = aH
        Ps.append(mode(P, sol, k, 0, tend, 's')); Pt.append(mode(P, sol, k, 0, tend, 't'))
        lnk.append(math.log(k))
    lnk = np.array(lnk); lP = np.log(Ps)
    c = np.polyfit(lnk - lnk[2], lP, 3)
    ns = 1 + c[2]; als = 2 * c[1]
    As = V0 * Ps[2]
    r = Pt[2] / Ps[2]
    nt = float(np.polyfit(lnk - lnk[2], np.log(Pt), 2)[1])
    # slow-roll cross-check at the pivot
    fN = float(np.interp(-Ns, -Nb, Y[0]))
    eV = 0.5 * (P.dV(fN) / P.V(fN)) ** 2; etaV = P.d2V(fN) / P.V(fN)
    # recalibrated V0 so that A_s = 2.1e-9 exactly, and the H_i that implies at S1
    V0c = AS / Ps[2]
    res.update(Nstar_used=Ns, As_at_trackcal=As, ns=ns, alpha_s=als, r=r, n_t=nt, ns_SR=1 - 6 * eV + 2 * etaV, r_SR=16 * eV,
               H_S1_for_As2p1=math.sqrt(V0c) * H_of_N(1.0) * MP, Hstar_over_HS1=H_of_N(Ns) / H_of_N(1.0),
               m_I_over_Hi=math.sqrt(V0 * P.d2V(P.phimin())) * MP / Hi)
    return res


def main(names):
    try:
        allr = json.load(open(OUT / 'sk_cmb.json'))
    except Exception:
        allr = dict(method=__doc__, runs={})
    for n in names:
        r = analyse(n)
        allr['runs'][n] = r
        print(json.dumps({k: (v if not isinstance(v, float) else float('%.6g' % v)) for k, v in r.items()}), flush=True)
        json.dump(allr, open(OUT / 'sk_cmb.json', 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1:] or ['HT6t', 'T1', 'E1'])
