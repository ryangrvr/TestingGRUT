#!/usr/bin/env python3
"""rp_reheating track (Claude, own code): driver for rpreheat.c (own C engine) at the heavy-mediator spot point.

Spot point [Q REVIEW_HEAVYWINDOW.md / hw_critic]: H_i = 1.3334962e11 GeV, M_Y = 30 H_i, h = 40 H_i, mu = 7.15945616 H_i,
kappa = 3.458684721e-3, lambda_X = 1e-3, lambda_Y = lambda_Sigma = 1, c_H = 10 (xi_X = -5/6 gives xi R = -10 H^2 in de
Sitter), Y_Q = 4.8294009077e-9.  Preparation: conditional minimum over (|X|, Y, |Sigma|) at held phases arg X = 0,
arg Sigma = -pi/2 with the actual initial curvature mass xi R(0) (own Newton solve from an analytic starting guess),
zero velocities, a = 1, H = H_i; the inflaton amplitude follows from the Friedmann constraint.
Yield per visible entropy: Y_X = N_X H_i^3 / (a^3 s), s = 2 pi^2 g* T^3 / 45, T = (30 rho_r / pi^2 g*)^(1/4), g* = 106.75
(the convention of ls_scaling, hw_verify and Codex).  Dilution factor relative to the instant release: D = (Rr/3)^(3/4),
Rr = rho_r a^4 at the end in H_i units (hw_verify's bookkeeping [D]).
"""
import json, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import root

H0 = 133349619398.33566          # GeV [Q INPUTS.json of the parent packet; = ls_scaling's normalised scale]
MRED, GSTAR = 2.435e18, 106.75
ETA2 = (H0 / MRED)**2
XI = -5.0 / 6.0
YQ = 4.8294009076988175e-9
P = dict(lx=1e-3, ly=1.0, ls=1.0, kap=0.0034586847211042523, M=30.0, h=40.0, mu=7.159456159513794, cH=10.0)
NX_INSTANT = 2788.1767940228406  # |N_X| of the instant release at this point [Q parent BACKGROUND late; hw_verify 2788.1]
SCR = Path(os.environ.get('RP_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                          'scratchpad/review/rp_reheating'))
BIN = SCR / 'build' / 'rpreheat'
SRC = Path(__file__).resolve().parent / 'rpreheat.c'
COLS = ['tau', 'a', 'H', 'R', 'NX', 'NY', 'NS', 'Gauss', 'EX', 'eps', 'NYH', 'NYL', 'Nsig', 'absX', 'absY', 'S_over_vev',
        'rho_r', 'rho_I', 'radfrac', 'Hcheck', 'Wcheck', 'w_I', 'm2curv', 'KV', 'Haux', 'rho']


def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-std=c11', '-Wall', '-o', str(BIN), str(SRC), '-lm'], check=True)


def potential(z, p=P):
    x, y, s = z
    x2, y2, s2 = abs(x)**2, abs(y)**2, abs(s)**2
    V = p['lx'] * x2**2 + p['M']**2 * y2 + p['ly'] * y2**2 + p['ls'] * (s2 - p['h']**2 / 2)**2 \
        + 2 * np.real(p['kap'] * np.conj(y) * x**3 + p['mu'] * np.conj(s) * y * y)
    F = np.array([2 * p['lx'] * x2 * x + 3 * p['kap'] * y * np.conj(x)**2,
                  p['M']**2 * y + 2 * p['ly'] * y2 * y + p['kap'] * x**3 + 2 * p['mu'] * s * np.conj(y),
                  2 * p['ls'] * (s2 - p['h']**2 / 2) * s + p['mu'] * y * y])
    return float(V), F


def ricci0(z, b, xi=XI, eta2=ETA2, mode=0):
    """R/H_i^2 at a zero-velocity state with H = 1: inflaton (or fluid) density from the Friedmann constraint."""
    V, F = potential(z)
    x2 = abs(z[0])**2
    f = 1 - 2 * xi * eta2 * x2
    need = 3 * f - eta2 * V                       # rho_I such that 3 f H^2 = rho at H = 1
    trI = 4 * need if mode == 0 else (need if mode in (1, 2) else 0.0)    # oscillator at rest: rho - 3p = 4 U
    return (trI + eta2 * (4 * V - 12 * xi * np.real(np.conj(z[0]) * F[0]))) / (f + 12 * xi**2 * eta2 * x2), need


def prepare(b=2.0, delta=0.0, xi=XI, eta2=ETA2, mode=0, cH_override=None):
    """Conditional stationary point at held phases: X = x > 0, Sigma = -i e^{i delta} s, Y free; curvature mass from
    the actual R(0) (modes 0, 1) or -cH rho_I/3 (mode 2).  Starting guess analytic: x^2 = cH/(2 lx), Y = -kap x^3/M^2,
    s = h/sqrt2 (not Codex's numbers)."""
    direc = -1j * np.exp(1j * delta)

    def state(q):
        return np.array([q[0], q[1] + 1j * q[2], direc * q[3]], complex)

    def m2c(z):
        if mode == 2:
            _, need = ricci0(z, b, xi, eta2, mode)
            return -(P['cH'] if cH_override is None else cH_override) * need / 3
        R, _ = ricci0(z, b, xi, eta2, mode)
        return xi * R

    def eq(q):
        z = state(q); V, F = potential(z); m = m2c(z)
        return [np.real(F[0] + m * z[0]), F[1].real, F[1].imag, np.real(F[2] * np.conj(direc))]
    x0 = np.sqrt(P['cH'] / (2 * P['lx']))
    q0 = [x0, -P['kap'] * x0**3 / P['M']**2, 0.0, P['h'] / np.sqrt(2)]
    sol = root(eq, q0, tol=1e-14)
    z = state(sol.x); res = float(np.linalg.norm(eq(sol.x)))
    R, need = ricci0(z, b, xi, eta2, mode)
    return z, dict(q=sol.x.tolist(), residual=res, R0=float(R), m2curv0=float(m2c(z)), rhoI0=float(need), success=bool(sol.success))


def run(tag, mode=0, b=2.0, gam=0.1, end=None, xi=XI, eta2=ETA2, rtol=1e-10, atol=1e-13, nsamp=2000, z0=None, delta=0.0,
        cH=None, save=True):
    build()
    if z0 is None:
        z0, info = prepare(b, delta, xi, eta2, mode)
    else:
        info = None
    end = end if end is not None else 25.0 / gam
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), str(mode), g(b), g(gam), g(xi), g(eta2), g(P['cH'] if cH is None else cH), g(P['lx']), g(P['ly']), g(P['ls']),
            g(P['kap']), g(P['M']), g(P['h']), g(P['mu']), g(end), g(rtol), g(atol), str(int(nsamp))] + \
           [g(v) for c in z0 for v in (c.real, c.imag)]
    t0 = time.time()
    out = subprocess.run(args, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(tag + ': ' + out.stderr)
    S, W, Y, END = [], {}, None, None
    for line in out.stdout.split('\n'):
        if line.startswith('S '):
            S.append([float(v) for v in line.split()[1:]])
        elif line.startswith('W') or line.startswith('A3040'):
            v = line.split(); W[v[0]] = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('Y '):
            Y = [float(v) for v in line.split()[1:]]
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_H_consistency=float(v[4]),
                                        max_hidden_work=float(v[5]))
    S = np.array(S)
    if save:
        (SCR / 'runs').mkdir(parents=True, exist_ok=True)
        np.save(SCR / 'runs' / (tag + '.npy'), S)
    last = dict(zip(COLS, S[-1]))
    a, rr = last['a'], last['rho_r']
    rho_rad = rr * MRED**2 * H0**2
    T = (30 * rho_rad / (np.pi**2 * GSTAR))**0.25
    s = 2 * np.pi**2 / 45 * GSTAR * T**3
    w = W.get('W80', last)
    nx = w['NX']
    yld = abs(nx) * H0**3 / (a**3 * s) / YQ
    Rr = rr * a**4
    D = (Rr / 3) ** 0.75
    res = dict(tag=tag, mode=mode, b=b, gamma=gam, end=end, xi=xi, eta2=eta2, rtol=rtol, delta=delta, prep=info,
               z0=[[c.real, c.imag] for c in z0], final_a=a, final_H=last['H'], final_radfrac=last['radfrac'], Rr_end=Rr,
               dilution=D, Tvis_GeV=T,
               late=dict(NX=w['NX'], EX=w['EX'], eps=w['eps'], NYH=w['NYH'], NYL=w['NYL'], Nsig=w['Nsig'], NY=w['NY']),
               late90=({k: W['W90'][k] for k in ('NX', 'EX', 'eps', 'NYH', 'NYL', 'Nsig')} if 'W90' in W else None),
               a3040=({k: W['A3040'][k] for k in ('NX', 'EX', 'eps', 'NYH', 'NYL', 'Nsig')} if 'A3040' in W else None),
               yield_over_YQ=yld, charge_ratio_to_instant=abs(nx) / NX_INSTANT, yield_decomposed=abs(nx) / NX_INSTANT / D,
               half_window_charge_change=(abs(W['W90']['NX'] / nx - 1) if 'W90' in W else None),
               checks=END, final_state=Y, seconds=time.time() - t0)
    return res, S


if __name__ == '__main__':
    build()
    z, info = prepare()
    print('prepared', z, info)
