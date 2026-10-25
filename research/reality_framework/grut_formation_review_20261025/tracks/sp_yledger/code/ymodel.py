#!/usr/bin/env python3
"""sp_yledger track (Claude, own code): parameters, inflationary minimum and a wrapper around yrelease.c.

Units H_i = 1. Couplings unchanged: M_Y=1e12, h=1e11, mu=5e11 GeV, kappa=1e-3, lambda_X=1e-3,
lambda_Y=lambda_Sigma=1. Preparation as stated in Codex 373d37e CLAUDE_HANDOFF.md: full aligned inflationary
radial minimum with -c_H H^2 |X|^2, X>0, Y<0, Sigma = -|Sigma| exp(i delta), zero physical velocity.
The minimum is found here by our own multi-start radial minimisation plus a full complex-gradient check.
"""
import os, subprocess, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

warnings.filterwarnings('ignore')
MRED = 2.435e18
GSTAR = 106.75
YQ = 4.8294009076988175e-9
MY, HV, MU = 1e12, 1e11, 5e11
LX, LY, LS, KAP, E_GAUGE = 1e-3, 1.0, 1.0, 1e-3, 1e-7
MXL = 0.092408963737536        # light-target m_X = g_X v [Q]
BDM = 4.4e-10                  # calibrated DM per visible entropy, GeV [Q]
SCR = Path(os.environ.get('SPY_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/sp_yledger/build'))
BIN = SCR / 'yrelease'
SRC = Path(__file__).resolve().parent / 'yrelease.c'
COLS = ['t', 'a', 'H', 'NX', 'NY', 'NS', 'EX', 'Ncirc', 'eps', 'NYH', 'NYL', 'JYcodex', 'Nsig', 'gauss', 'ledger',
        'a4rho', 'BoverA', 'EYexc', 'Esig', 'phi_frac', 'EYfull', 'Neff']


def build():
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        SCR.mkdir(parents=True, exist_ok=True)
        subprocess.run(['gcc', '-O3', '-o', str(BIN), str(SRC), '-lm'], check=True)


def params(H):
    return dict(M=MY / H, h=HV / H, mu=MU / H)


def V_radial(R, H, cH, c2=1.0):
    p = params(H)
    x, y, s = R
    return (-cH * x**2 + LX * x**4 + LY * y**4 + LS * (s**2 - p['h']**2 / 2)**2 + p['M']**2 * y**2
            - 2 * KAP * y * x**3 - 2 * p['mu'] * s * y**2 * c2)


def minimum(H, cH):
    """Aligned radial minimum (|X|,|Y|,|Sigma|) in H units: multi-start Nelder-Mead/BFGS in log variables,
    then a Newton polish of the radial gradient, accepted only if it is finite and lowers |grad V|."""
    from scipy.optimize import fsolve
    p = params(H)

    def grad(R):
        x, y, s = R
        return np.array([-2 * cH * x + 4 * LX * x**3 - 6 * KAP * y * x**2,
                         4 * LY * y**3 + 2 * p['M']**2 * y - 2 * KAP * x**3 - 4 * p['mu'] * s * y,
                         4 * LS * s * (s**2 - p['h']**2 / 2) - 2 * p['mu'] * y**2])
    best = None
    for lx in (1.0, 2.5, 3.5):
        for ly in (-2.0, 0.0, 1.5):
            for lsg in (np.log(p['h'] / np.sqrt(2)), 0.0):
                f = lambda q: V_radial(np.exp(q), H, cH)
                r = minimize(f, [lx, ly, lsg], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-16, maxiter=20000))
                r = minimize(f, r.x, method='BFGS', options=dict(gtol=1e-13))
                if np.all(np.isfinite(r.x)) and np.isfinite(r.fun) and (best is None or r.fun < best.fun):
                    best = r
    R0 = np.exp(best.x)
    R1 = fsolve(grad, R0, xtol=1e-15)
    if np.all(np.isfinite(R1)) and np.all(R1 > 0) and np.linalg.norm(grad(R1)) < np.linalg.norm(grad(R0)):
        R0 = R1
    x, y, s = R0
    X, Y, S = x + 0j, -y + 0j, -s + 0j
    F = np.array([2 * LX * abs(X)**2 * X + 3 * KAP * Y * np.conj(X)**2 - cH * X,
                  p['M']**2 * Y + 2 * LY * abs(Y)**2 * Y + KAP * X**3 + 2 * p['mu'] * S * np.conj(Y),
                  2 * LS * (abs(S)**2 - p['h']**2 / 2) * S + p['mu'] * Y * Y])
    return dict(x=x, y=y, s=s, V=V_radial(R0, H, cH), grad_norm=float(np.linalg.norm(F)))


def initial(H, cH, delta):
    m = minimum(H, cH)
    return np.array([m['x'], -m['y'], -m['s'] * np.exp(1j * delta)]), m


def run(H, z0, t_end, bg=0, gphi=0.0, cH=1.0, rtol=1e-11, atol=1e-14, nsamp=400):
    build()
    p = params(H)
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), str(int(bg)), g(gphi), g(cH), g(p['M']), g(p['h']), g(p['mu']), g(t_end),
            g(rtol), g(atol), str(int(nsamp))] + [g(v) for c in z0 for v in (c.real, c.imag)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout.split('\n')
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


def entropy_inst(H):
    """Instant visible reheating at H_i (GeV): s_i and T_R."""
    T = (90 * MRED**2 * H**2 / (np.pi**2 * GSTAR))**0.25
    return 2 * np.pi**2 / 45 * GSTAR * T**3, T


def TR_of_Gamma(G_GeV):
    return (90 / (np.pi**2 * GSTAR))**0.25 * np.sqrt(G_GeV * MRED)


def Gamma_of_TR(TR):
    return TR**2 / MRED * np.sqrt(np.pi**2 * GSTAR / 90)
