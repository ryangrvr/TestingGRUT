#!/usr/bin/env python3
"""ls_scaling track (Claude, own code): couplings-general release, ledger and threshold helpers.

Release: scalrel.c (parametrised copy of our sp_yledger yrelease.c; cosmic time, physical fields, DP5(4)).
Minimum: own multi-start radial minimisation of
    V = -cH|X|^2 + lx|X|^4 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + M^2|Y|^2 + 2Re(kap Y* X^3) + 2Re(mu S* Y^2)
with X > 0, Y < 0, Sigma = -|Sigma| e^{i delta} (the Codex 373d37e / b448f25 preparation), zero physical velocity,
instant radiation release at a = 1.  All quantities in H_i units unless a name ends in _GeV.

Ledger [D]: same disjoint decay ledger as Codex b448f25 and our ls_ledger track, written in closed form:
    R_Y    = (N_YH + N_YL)/E_X * (m_L/H) <a>_L,   <a>_L = sqrt(pi H/(2 Gamma_L)),  Gamma_L = 3 kap^2 m_L/(256 pi^3)
           = 64.47 kap^-1 sqrt(m_L/H) N_Y/E_X
    R_Sig  = 2 N_sig/E_X * (m_A/H) <a>_A,         Gamma_A = e^2 m_A/(48 pi) (1-4m_X^2/m_A^2)^{3/2},  m_A = 6 e h
           = 24 pi sqrt(h/(e H)) N_sig/E_X        (threshold factor ~1)
    R_casc = N_YH/E_X * (m_A/H) <a>_A             (one cascade vector per Y_H; central ledger only)
Threshold: case-L R_L(eps) at normalised yield, log-log interpolation of our own sp_critic rmax_relic tables
(eps = 0.0677..0.2402) [C, earlier track]; values outside the table are flagged, never extrapolated silently.
"""
import os, subprocess, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, fsolve

warnings.filterwarnings('ignore')

MRED = 2.435e18
GSTAR = 106.75
YQ = 4.8294009076988175e-9
MXL = 0.092408963737536          # light-target m_X = g_X v [Q]
STATED = dict(lx=1e-3, ly=1.0, ls=1.0, kap=1e-3, M_GeV=1e12, h_GeV=1e11, mu_GeV=5e11, e=1e-7)
H_CAND = 434867798206.6179
SCR = Path(os.environ.get('LSS_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                          'scratchpad/review/ls_scaling/build'))
BIN = SCR / 'scalrel'
SRC = Path(__file__).resolve().parent / 'scalrel.c'
COLS = ['t', 'a', 'H', 'NX', 'NY', 'NS', 'EX', 'Ncirc', 'eps', 'NYH', 'NYL', 'JYcodex', 'Nsig', 'gauss', 'ledger',
        'a4rho', 'BoverA', 'EYexc', 'Esig', 'phi_frac', 'EYfull', 'Neff']
# our sp_critic case-L thresholds at normalised yield (rmax_relic_Llow.json, rmax_relic_L.json) [C]
RL_TABLE = [(0.0677, 29.756595100551028), (0.1, 49.83550533984897), (0.15, 91.80895700634682),
            (0.19, 135.89463084582223), (0.240201575085664, 172.49014985200023)]


def build():
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        SCR.mkdir(parents=True, exist_ok=True)
        subprocess.run(['gcc', '-O3', '-o', str(BIN), str(SRC), '-lm'], check=True)


def pars(cp, H):
    return dict(lx=cp['lx'], ly=cp['ly'], ls=cp['ls'], kap=cp['kap'], M=cp['M_GeV'] / H, h=cp['h_GeV'] / H,
                mu=cp['mu_GeV'] / H)


def V_radial(R, p, cH):
    x, y, s = R
    return (-cH * x**2 + p['lx'] * x**4 + p['ly'] * y**4 + p['ls'] * (s**2 - p['h']**2 / 2)**2 + p['M']**2 * y**2
            - 2 * p['kap'] * y * x**3 - 2 * p['mu'] * s * y**2)


def grad_radial(R, p, cH):
    x, y, s = R
    return np.array([-2 * cH * x + 4 * p['lx'] * x**3 - 6 * p['kap'] * y * x**2,
                     4 * p['ly'] * y**3 + 2 * p['M']**2 * y - 2 * p['kap'] * x**3 - 4 * p['mu'] * s * y,
                     4 * p['ls'] * s * (s**2 - p['h']**2 / 2) - 2 * p['mu'] * y**2])


def hess_radial(R, p, cH, d=1e-6):
    H = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = d * max(abs(R[j]), 1e-8)
        H[:, j] = (grad_radial(R + e, p, cH) - grad_radial(R - e, p, cH)) / (2 * e[j])
    return 0.5 * (H + H.T)


def minima(cp, H, cH):
    """All distinct radial local minima found by multi-start (log variables) + Newton polish; sorted by |X|."""
    p = pars(cp, H)
    x_ref = np.sqrt(cH / (2 * p['lx']))
    found = []
    for lxs in (np.log(x_ref) - 0.5, np.log(x_ref), np.log(x_ref) + 0.7, np.log(x_ref) + 1.5):
        for lys in (-6.0, -3.0, 0.0, 2.0):
            for lss in (np.log(p['h'] / np.sqrt(2)), np.log(p['h'] / np.sqrt(2)) + 1.0):
                f = lambda q: V_radial(np.exp(q), p, cH)
                r = minimize(f, [lxs, lys, lss], method='Nelder-Mead',
                             options=dict(xatol=1e-12, fatol=1e-18, maxiter=40000))
                r = minimize(f, r.x, method='BFGS', options=dict(gtol=1e-14))
                if not np.all(np.isfinite(r.x)):
                    continue
                R0 = np.exp(r.x)
                R1 = fsolve(lambda R: grad_radial(R, p, cH), R0, xtol=1e-15)
                if np.all(np.isfinite(R1)) and np.all(R1 > 0) and \
                        np.linalg.norm(grad_radial(R1, p, cH)) < np.linalg.norm(grad_radial(R0, p, cH)):
                    R0 = R1
                ev = np.linalg.eigvalsh(hess_radial(R0, p, cH))
                if ev.min() <= 0:
                    continue
                if not any(np.allclose(R0, q['R'], rtol=1e-6) for q in found):
                    found.append(dict(R=R0, V=V_radial(R0, p, cH), grad=float(np.linalg.norm(grad_radial(R0, p, cH))),
                                      min_hess=float(ev.min())))
    found.sort(key=lambda q: q['R'][0])
    for q in found:
        x, y, s = q['R']
        q.update(x=x, y=y, s=s, x_over_ref=x / x_ref, lamY_y2_over_M2=p['ly'] * y * y / p['M']**2,
                 kapx2_over_M2=p['kap'] * x * x / p['M']**2)
    return found


def initial(cp, H, cH, delta, which='small'):
    ms = minima(cp, H, cH)
    m = ms[0] if which == 'small' else min(ms, key=lambda q: q['V'])
    return np.array([m['x'], -m['y'], -m['s'] * np.exp(1j * delta)]), m, ms


def run(cp, H, z0, a_end=300.0, rtol=1e-11, atol=1e-14, nsamp=80):
    build()
    p = pars(cp, H)
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']), g(p['M']), g(p['h']), g(p['mu']),
            g(0.5 * a_end**2), g(rtol), g(atol), str(int(nsamp))] + [g(v) for c in z0 for v in (c.real, c.imag)]
    out = subprocess.run(args, capture_output=True, text=True).stdout.split('\n')
    S, W, F, END = [], {}, None, None
    for line in out:
        if line.startswith('S '):
            S.append(dict(zip(COLS, map(float, line.split()[1:]))))
        elif line.startswith('W '):
            v = line.split()
            W[f'{float(v[1]):g}-{float(v[2]):g}'] = dict(span=float(v[3]), **dict(zip(COLS, map(float, v[4:]))))
        elif line.startswith('F '):
            F = list(map(float, line.split()[1:]))
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_ledger=float(v[4]))
        elif line.startswith('FAIL'):
            raise RuntimeError(line)
    return dict(samples=S, windows=W, final=F, end=END)


def masses(cp):
    M, h, mu, e = cp['M_GeV'], cp['h_GeV'], cp['mu_GeV'], cp['e']
    mH = np.sqrt(M**2 + np.sqrt(2) * mu * h); mL2 = M**2 - np.sqrt(2) * mu * h
    mL = np.sqrt(mL2) if mL2 > 0 else float('nan')
    return dict(mYH=mH, mYL=mL, msig=np.sqrt(2 * cp['ls']) * h, mA=6 * e * h, r=np.sqrt(2) * mu * h / M**2)


def ledger(cp, H, NYH, NYL, Nsig, EX):
    m = masses(cp); kap, e = cp['kap'], cp['e']
    GL = 3 * kap**2 * m['mYL'] / (256 * np.pi**3)
    thr = (1 - 4 * MXL**2 / m['mA']**2)**1.5 if m['mA'] > 2 * MXL else float('nan')
    GA = e**2 * m['mA'] / (48 * np.pi) * thr
    aL = np.sqrt(np.pi * H / (2 * GL)); aA = np.sqrt(np.pi * H / (2 * GA))
    dep_L = m['mYL'] / H * aL          # conformal deposit per Y quantum (all Y counted at m_L, as Codex/ls_ledger)
    dep_A = m['mA'] / H * aA           # per vector
    RY = (NYH + NYL) * dep_L / EX
    RS = 2 * Nsig * dep_A / EX
    RC = NYH * dep_A / EX
    flags = []
    if not m['msig'] > 2 * m['mA']:
        flags.append('sigma->AA closed')
    if m['msig'] > 2 * m['mYL']:
        flags.append('sigma->Y_L Y_L open (vector channel kept: conservative-high for R_Sigma)')
    return dict(R_Y=RY, R_Sigma=RS, R_cascade=RC, R_lower=RY + RS, R_central=RY + RS + RC,
                deposit_per_YL=dep_L, deposit_per_A=dep_A, a_dec_L=aL, a_dec_A=aA,
                closed_form_R_Y=64.47 / kap * np.sqrt(m['mYL'] / H) * (NYH + NYL) / EX,
                closed_form_R_Sigma=24 * np.pi * np.sqrt(cp['h_GeV'] / (e * H)) * Nsig / EX / np.sqrt(thr),
                masses=m, flags=flags)


def R_L(eps):
    """Case-L threshold at normalised yield (our sp_critic tables), log-log interpolation; None outside table."""
    e = abs(eps)
    xs = [np.log(a) for a, _ in RL_TABLE]; ys = [np.log(b) for _, b in RL_TABLE]
    if e < RL_TABLE[0][0] or e > RL_TABLE[-1][0]:
        return None
    return float(np.exp(np.interp(np.log(e), xs, ys)))


def yield_over_YQ(NX, H):
    T = (90 * MRED**2 * H**2 / (np.pi**2 * GSTAR))**0.25
    s = 2 * np.pi**2 / 45 * GSTAR * T**3
    return abs(NX) * H**3 / s / YQ


def c6(cp):
    M, h, mu, kap = cp['M_GeV'], cp['h_GeV'], cp['mu_GeV'], cp['kap']
    Sig = h / np.sqrt(2)
    D = M**4 - 4 * (mu * Sig)**2
    a6 = mu * kap**2 * Sig / D
    return dict(alpha6_GeVm2=a6, c6=a6 * MRED**2, c6_ceiling=1.486e13, beta6_GeVm2=-M**2 * kap**2 / D)


def kappa_max(cp):
    """Tree-level large-field boundedness along the aligned direction (Sigma subleading at large field) [D]:
    lx > (3/2) kap^{4/3} (2 ly)^{-1/3}  <=>  kap < (lx/1.1906)^{3/4} ly^{1/4}."""
    return (cp['lx'] / (1.5 * 2**(-1 / 3)))**0.75 * cp['ly']**0.25


def V_full(z, p, cH):
    x, y, s = z
    ax2, ay2, as2 = abs(x)**2, abs(y)**2, abs(s)**2
    return float(-cH * ax2 + p['lx'] * ax2**2 + p['ly'] * ay2**2 + p['ls'] * (as2 - p['h']**2 / 2)**2
                 + p['M']**2 * ay2 + 2 * np.real(p['kap'] * np.conj(y) * x**3 + p['mu'] * np.conj(s) * y * y))


def relaxed_initial(cp, H, cH, delta, start=None):
    """Y-relaxed supplied-phase preparation (our variant, not Codex's): conditional minimum over |X|, complex Y
    and |Sigma| at fixed arg X = 0 and Sigma = -|Sigma| e^{i delta}.  A heavy Y (M_Y > 3H/2) tracks this
    conditional minimum during inflation; Codex's prescription instead keeps the delta = 0 Y."""
    p = pars(cp, H)
    if start is None:
        z0, m, ms = initial(cp, H, cH, delta)
        start = [np.log(m['x']), -m['y'], 0.0, np.log(m['s'])]
    ph = np.exp(1j * delta)

    def f(q):
        x = np.exp(q[0]); y = q[1] + 1j * q[2]; s = np.exp(q[3])
        return V_full(np.array([x, y, -s * ph]), p, cH)
    best = None
    for jitter in (0.0, 0.3, -0.3):
        q0 = np.array(start, float); q0[2] += jitter * abs(q0[1])
        r = minimize(f, q0, method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-18, maxiter=60000))
        r = minimize(f, r.x, method='BFGS', options=dict(gtol=1e-13))
        if best is None or r.fun < best.fun:
            best = r
    q = best.x
    z = np.array([np.exp(q[0]), q[1] + 1j * q[2], -np.exp(q[3]) * ph])
    # conditional gradient (radial X, full Y, radial Sigma)
    x, y, s = z
    gX = 2 * p['lx'] * abs(x)**2 * x + 3 * p['kap'] * y * np.conj(x)**2 - cH * x
    gY = p['M']**2 * y + 2 * p['ly'] * abs(y)**2 * y + p['kap'] * x**3 + 2 * p['mu'] * s * np.conj(y)
    gS = 2 * p['ls'] * (abs(s)**2 - p['h']**2 / 2) * s + p['mu'] * y * y
    cond = abs(np.real(gX * np.conj(x)) / abs(x)) + abs(gY) + abs(np.real(gS * np.conj(s)) / abs(s))
    torque_X = float(np.imag(np.conj(x) * gX))
    return z, dict(V=best.fun, cond_grad=float(cond), torque_X=torque_X, x=float(abs(x)), Y=[float(y.real), float(y.imag)],
                   s=float(abs(s)))
