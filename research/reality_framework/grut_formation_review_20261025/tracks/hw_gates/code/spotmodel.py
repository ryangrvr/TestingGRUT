#!/usr/bin/env python3
"""hw_gates (Claude, own code): shared model of the ls_scaling spot point, H_i = 1 units.

Spot point (ls_scaling B_base_Hscaled, yield 1.00 Y_Q): H = 1.3335e11 GeV, M_Y = 30 H, h = 40 H, mu = 7.1595 H
(r = sqrt2 mu h/M_Y^2 = 0.45), kappa = 3.4587e-3, lambda_X = 1e-3, lambda_Y = lambda_Sigma = 1, e = 1e-7, c_H = 10.
Potential during inflation (physical fields):
  V = -cH|X|^2 + lx|X|^4 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + M^2|Y|^2 + 2Re(kap Y* X^3) + 2Re(mu S* Y^2).
Charges X:Y:Sigma = 1:3:6. Gauge-invariant phases (z = R e^{i theta}): psi1 = theta_Y - 3 theta_X,
psi2 = theta_Sigma - 2 theta_Y, delta = theta_Sigma - 6 theta_X = psi2 + 2 psi1 (the relative phase that sources
the X charge in the Y-eliminated EFT: V_eff ⊃ 2|alpha6| X^6 |Sigma|... cos(delta - pi)).

Preparations at supplied delta (Sigma = -s e^{i delta}, X = x > 0):
  codex   : aligned radial minimum, Y left real (= exact minimum of a pure mu-rotation, mu -> mu e^{i delta});
  relaxed : conditional minimum over (x, Re Y, Im Y, s) at fixed arg X, arg Sigma (ls_scaling's variant; = exact
            minimum when the holding term acts on X^6 Sigma* only, i.e. does not enter the Y equation);
  kaprot  : aligned moduli with Y -> -y e^{i delta/2}, Sigma -> -s e^{i delta} (= exact minimum of a kappa-rotation,
            kappa -> kappa e^{i delta/2}, written in the gauge where X is real).
Hessian: sympy-exact 6x6 in canonical real components (x1,x2,y1,y2,s1,s2), lambdified once with coupling symbols.
"""
import os, subprocess
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.optimize import minimize, fsolve

H_SPOT = 133350496019.93945           # GeV [Q, own ls_scaling track]
MRED, MPL, GSTAR = 2.435e18, 1.2209e19, 106.75
YQ = 4.8294009076988175e-9
MX_LIGHT = 0.1 * 0.924089637
Q = np.array([1.0, 3.0, 6.0])
SPOT = dict(lx=1e-3, ly=1.0, ls=1.0, kap=0.0034586847211042523, M=30.0, h=40.0, mu=7.159456159513793, e=1e-7, cH=10.0)
SCR = Path(os.environ.get('HWG_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                          'scratchpad/review/hw_gates'))
BIN = SCR / 'build' / 'hwgrel'
SRC = Path(__file__).resolve().parent / 'hwgrel.c'
COLS = ['t', 'a', 'H', 'NX', 'NY', 'NS', 'EX', 'Ncirc', 'eps', 'NYH', 'NYL', 'JYcodex', 'Nsig', 'gauss', 'ledger',
        'a4rho', 'BoverA', 'EYexc', 'Esig', 'thX', 'thS', 'thG', 'absS_over_vev', 'absX', 'absY', 'F2']


def params(**over):
    p = dict(SPOT); p.update(over); return p


# ------------------------------------------------------------------ potential and preparations
def V(z, p, cH=None):
    cH = p['cH'] if cH is None else cH
    x, y, s = z
    ax2, ay2, as2 = abs(x)**2, abs(y)**2, abs(s)**2
    return float(-cH * ax2 + p['lx'] * ax2**2 + p['ly'] * ay2**2 + p['ls'] * (as2 - p['h']**2 / 2)**2
                 + p['M']**2 * ay2 + 2 * np.real(p['kap'] * np.conj(y) * x**3 + p['mu'] * np.conj(s) * y * y))


def grad(z, p, cH=None):
    """dV/dz_i^* (complex gradient; the force in the field equations)."""
    cH = p['cH'] if cH is None else cH
    x, y, s = z
    ax2, ay2, as2 = abs(x)**2, abs(y)**2, abs(s)**2
    return np.array([-cH * x + 2 * p['lx'] * ax2 * x + 3 * p['kap'] * y * np.conj(x)**2,
                     p['M']**2 * y + 2 * p['ly'] * ay2 * y + p['kap'] * x**3 + 2 * p['mu'] * s * np.conj(y),
                     2 * p['ls'] * (as2 - p['h']**2 / 2) * s + p['mu'] * y * y])


def aligned(p, cH=None):
    """Aligned radial minimum X = x, Y = -y, Sigma = -s (x, y, s > 0): multistart in log variables + Newton polish."""
    cH = p['cH'] if cH is None else cH
    f = lambda q: V(np.array([np.exp(q[0]), -np.exp(q[1]), -np.exp(q[2])]), p, cH)
    xr = np.sqrt(cH / (2 * p['lx']))
    best = None
    for lx0 in (np.log(xr) - 0.3, np.log(xr), np.log(xr) + 0.5):
        for ly0 in (-3.0, 0.0, 1.5):
            for ls0 in (np.log(p['h'] / np.sqrt(2)), np.log(p['h'] / np.sqrt(2)) + 0.5):
                r = minimize(f, [lx0, ly0, ls0], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-18, maxiter=30000))
                if best is None or r.fun < best.fun:
                    best = r
    q = np.exp(best.x)
    g = lambda q: [np.real(grad(np.array([q[0], -q[1], -q[2]]), p, cH)[i] * np.array([1, -1, -1])[i]) for i in range(3)]
    q = fsolve(g, q, xtol=1e-15)
    return q


def prep(kind, delta, p, cH=None):
    cH = p['cH'] if cH is None else cH
    x, y, s = aligned(p, cH)
    ph = np.exp(1j * delta)
    if kind == 'codex':
        return np.array([x, -y, -s * ph]), dict(x=x, y=y, s=s)
    if kind == 'kaprot':
        return np.array([x, -y * np.exp(1j * delta / 2), -s * ph]), dict(x=x, y=y, s=s)
    if kind == 'relaxed':
        def f(q):
            return V(np.array([np.exp(q[0]), q[1] + 1j * q[2], -np.exp(q[3]) * ph]), p, cH)
        best = None
        for jit in (0.0, 0.4, -0.4):
            q0 = [np.log(x), -y, jit * y, np.log(s)]
            r = minimize(f, q0, method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-18, maxiter=60000))
            r = minimize(f, r.x, method='BFGS', options=dict(gtol=1e-13))
            if best is None or r.fun < best.fun:
                best = r
        q = best.x
        z = np.array([np.exp(q[0]), q[1] + 1j * q[2], -np.exp(q[3]) * ph])
        g = grad(z, p, cH)
        torque = float(np.imag(np.conj(z[0]) * g[0]))
        return z, dict(x=float(abs(z[0])), Y=[float(z[1].real), float(z[1].imag)], s=float(abs(z[2])), V=best.fun,
                       torque_X=torque, cond_grad=float(abs(np.real(np.conj(z[0]) * g[0]) / abs(z[0])) + abs(g[1])
                                                        + abs(np.real(np.conj(z[2]) * g[2]) / abs(z[2]))))
    raise ValueError(kind)


# ------------------------------------------------------------------ Hessian (canonical real components)
_SYM = {}


def _build():
    if _SYM:
        return _SYM
    q = sp.symbols('x1 x2 y1 y2 s1 s2', real=True)
    lx, ly, ls, kap, M, h, mu, cH = sp.symbols('lx ly ls kap M h mu cH', real=True)
    kpr, kpi, mpr, mpi, c7r, c7i, cS, cY = sp.symbols('kpr kpi mpr mpi c7r c7i cS cY', real=True)
    X = (q[0] + sp.I * q[1]) / sp.sqrt(2); Y = (q[2] + sp.I * q[3]) / sp.sqrt(2); S = (q[4] + sp.I * q[5]) / sp.sqrt(2)
    aX = sp.expand(X * sp.conjugate(X)); aY = sp.expand(Y * sp.conjugate(Y)); aS = sp.expand(S * sp.conjugate(S))
    kapc = kap + kpr + sp.I * kpi; muc = mu + mpr + sp.I * mpi; c7 = c7r + sp.I * c7i
    inter = kapc * sp.conjugate(Y) * X**3 + muc * sp.conjugate(S) * Y**2 + c7 * X**6 * sp.conjugate(S)
    Vs = (-cH * aX + lx * aX**2 + ly * aY**2 + ls * (aS - h**2 / 2)**2 + M**2 * aY + cS * aS + cY * aY
          + 2 * sp.re(sp.expand(inter)))
    Vs = sp.expand(Vs)
    g = [sp.diff(Vs, v) for v in q]
    Hs = [[sp.diff(gi, v) for v in q] for gi in g]
    pars = (lx, ly, ls, kap, M, h, mu, cH, kpr, kpi, mpr, mpi, c7r, c7i, cS, cY)
    _SYM['V'] = sp.lambdify([q, pars], Vs, 'numpy')
    _SYM['g'] = sp.lambdify([q, pars], g, 'numpy')
    _SYM['H'] = sp.lambdify([q, pars], Hs, 'numpy')
    return _SYM


def _pars(p, extra):
    e = dict(kpr=0, kpi=0, mpr=0, mpi=0, c7r=0, c7i=0, cS=0, cY=0); e.update(extra or {})
    return (p['lx'], p['ly'], p['ls'], p['kap'], p['M'], p['h'], p['mu'], p['cH'],
            e['kpr'], e['kpi'], e['mpr'], e['mpi'], e['c7r'], e['c7i'], e['cS'], e['cY'])


def canon(z):
    return np.sqrt(2) * np.ravel([[c.real, c.imag] for c in z])


def hessian(z, p, extra=None):
    S = _build()
    return np.array(S['H'](canon(z), _pars(p, extra)), float)


def gradient_canon(z, p, extra=None):
    S = _build()
    return np.array(S['g'](canon(z), _pars(p, extra)), float)


def spectrum(z, p, extra=None):
    """Full 6x6 spectrum with the exact gauge zero mode (dz_j = i q_j z_j) projected out; phase content of each mode:
    d theta_j per unit canonical amplitude, and d delta = d thS - 6 d thX, d psi1, d psi2; f_eff = 1/|d delta|."""
    Hs = hessian(z, p, extra)
    gz = 1j * Q * z
    g = canon(gz); g /= np.linalg.norm(g)
    P = np.eye(6) - np.outer(g, g)
    w, U = np.linalg.eigh(P @ Hs @ P)
    gi = int(np.argmax(np.abs(U.T @ g)))
    rows = []
    for k in range(6):
        v = U[:, k]
        dz = (v[0::2] + 1j * v[1::2]) / np.sqrt(2)
        dth = np.imag(dz / z)            # z = R e^{i th}: d th = Im(dz/z)
        dlnR = np.real(dz / z)
        rows.append(dict(m2=float(w[k]), gauge=(k == gi), dtheta=[float(t) for t in dth],
                         dlnR=[float(t) for t in dlnR], ddelta=float(dth[2] - 6 * dth[0]),
                         dpsi1=float(dth[1] - 3 * dth[0]), dpsi2=float(dth[2] - 2 * dth[1]),
                         angular_weight=float(np.sum((np.abs(z) * dth)**2) / (np.sum((np.abs(z) * dth)**2) + np.sum((np.abs(z) * dlnR)**2)))))
    phys = [r for r in rows if not r['gauge']]
    ang = [r for r in phys if r['angular_weight'] > 0.5]
    light = min(ang, key=lambda r: r['m2']) if ang else None
    return dict(modes=rows, light_phase=(dict(light, mass_over_H=float(np.sqrt(max(light['m2'], 0))),
                                             f_eff=float(1 / abs(light['ddelta'])) if light['ddelta'] else None)
                                        if light else None),
                hessian_unprojected_eigs=np.linalg.eigvalsh(Hs).tolist())


def higgs_norm(z):
    """F_g^2 = 2 sum q_i^2 |z_i|^2 (H units); m_A = e F_g."""
    return float(2 * np.sum(Q**2 * np.abs(z)**2))


# ------------------------------------------------------------------ release driver (hwgrel.c)
def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-o', str(BIN), str(SRC), '-lm'], check=True)


def run(z0, p, a_end=60.0, rtol=1e-10, atol=1e-13, nsamp=60, win0=30.0, extra=None):
    """extra: dict(pmode, Dl, muP, kapP, c7, cX). pmode 0 none, 1 exp(-(a-1)/Dl), 2 H-tracking (A-terms a^-2,
    masses a^-4), 3 constant (never switched off)."""
    build()
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']), g(p['M']), g(p['h']), g(p['mu']),
            g(0.5 * a_end**2), g(rtol), g(atol), str(int(nsamp)), g(win0)] + [g(v) for c in z0 for v in (c.real, c.imag)]
    env = dict(os.environ)
    ex = dict(pmode=0, Dl=1.0, muP=0j, kapP=0j, c7=0j, cX=0.0); ex.update(extra or {})
    env.update(HWG_PMODE=str(int(ex['pmode'])), HWG_DL=g(ex['Dl']), HWG_MUPR=g(np.real(ex['muP'])), HWG_MUPI=g(np.imag(ex['muP'])),
               HWG_KPR=g(np.real(ex['kapP'])), HWG_KPI=g(np.imag(ex['kapP'])), HWG_C7R=g(np.real(ex['c7'])),
               HWG_C7I=g(np.imag(ex['c7'])), HWG_CX=g(ex['cX']))
    out = subprocess.run(args, capture_output=True, text=True, env=env)
    S, W, END = [], None, None
    for line in out.stdout.split('\n'):
        if line.startswith('S '):
            S.append(dict(zip(COLS, map(float, line.split()[1:]))))
        elif line.startswith('W '):
            v = line.split(); W = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_ledger=float(v[4]),
                                        minS_over_vev=float(v[5]))
        elif line.startswith('FAIL'):
            raise RuntimeError(line + out.stderr)
    if W is None:
        raise RuntimeError('no window: ' + out.stderr + out.stdout[:500])
    return dict(samples=S, window=W, end=END)


# ------------------------------------------------------------------ ledger (same as ls_scaling / Codex b448f25) [D]
def masses_GeV(p, H=H_SPOT):
    M, h, mu, e = p['M'] * H, p['h'] * H, p['mu'] * H, p['e']
    mL2 = M**2 - np.sqrt(2) * mu * h
    return dict(mYH=np.sqrt(M**2 + np.sqrt(2) * mu * h), mYL=np.sqrt(mL2) if mL2 > 0 else float('nan'),
                msig=np.sqrt(2 * p['ls']) * h, mA=6 * e * h, r=np.sqrt(2) * mu * h / M**2)


def ledger(p, w, H=H_SPOT):
    m = masses_GeV(p, H); kap, e = p['kap'], p['e']
    GL = 3 * kap**2 * m['mYL'] / (256 * np.pi**3)
    thr = (1 - 4 * MX_LIGHT**2 / m['mA']**2)**1.5
    GA = e**2 * m['mA'] / (48 * np.pi) * thr
    aL = np.sqrt(np.pi * H / (2 * GL)); aA = np.sqrt(np.pi * H / (2 * GA))
    dep_L = m['mYL'] / H * aL; dep_A = m['mA'] / H * aA
    EX = w['EX']
    RY = (w['NYH'] + w['NYL']) * dep_L / EX; RS = 2 * w['Nsig'] * dep_A / EX; RC = w['NYH'] * dep_A / EX
    return dict(R_Y=RY, R_Sigma=RS, R_cascade=RC, R_lower=RY + RS, R_central=RY + RS + RC,
                deposit_per_A_in_EX=dep_A / EX, deposit_per_YL_in_EX=dep_L / EX, Gamma_A_GeV=GA, Gamma_L_GeV=GL,
                a_dec_A=aA, a_dec_L=aL, masses=m)


RL_TABLE = [(0.0677, 29.756595100551028), (0.1, 49.83550533984897), (0.15, 91.80895700634682),
            (0.19, 135.89463084582223), (0.240201575085664, 172.49014985200023)]   # our sp_critic case-L [C]


def R_L(eps):
    e = abs(eps)
    if e < RL_TABLE[0][0] or e > RL_TABLE[-1][0]:
        return None
    xs = [np.log(a) for a, _ in RL_TABLE]; ys = [np.log(b) for _, b in RL_TABLE]
    return float(np.exp(np.interp(np.log(e), xs, ys)))


def entropy_RH(H=H_SPOT):
    T = (90 * MRED**2 * H**2 / (np.pi**2 * GSTAR))**0.25
    return 2 * np.pi**2 / 45 * GSTAR * T**3, T


def yield_over_YQ(NX, H=H_SPOT):
    s, _ = entropy_RH(H)
    return abs(NX) * H**3 / s / YQ


def alpha6(p, H=H_SPOT):
    M, h, mu, kap = p['M'] * H, p['h'] * H, p['mu'] * H, p['kap']
    Sig = h / np.sqrt(2); D = M**4 - 4 * (mu * Sig)**2
    a6 = mu * kap**2 * Sig / D
    return dict(alpha6_GeVm2=a6, c6=a6 * MRED**2, beta6_GeVm2=-M**2 * kap**2 / D, r=2 * mu * Sig / M**2)
