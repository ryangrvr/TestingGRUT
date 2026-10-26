#!/usr/bin/env python3
"""ls_critic (Claude): semi-independent check of ls_scaling's Y-relaxed preparation numbers at candidates 10 / 17.

ls_scaling (scal.relaxed_initial: Nelder-Mead + BFGS on V; scalrel.c release; internal C projection; R_L from a
log-interpolated table) reported R_lower 148.9 / 192.8 with the Y-relaxed preparation.  Here:
  * the conditional minimum (arg X = 0, Sigma = -|Sigma| e^{i delta}; free |X|, complex Y, |Sigma|) is found by
    solving the analytic stationarity equations with scipy root (hybr) from the aligned radial minimum, and its
    4x4 Hessian is checked;
  * the release uses our sp_yledger yrelease.c via ymodel.run (separate binary build) to a = 1000, and the late state
    is projected with ls_ledger's Python projection (ls_project.project) as well as the C window means;
  * the decay ledger is re-written here in closed form (widths and <a> = sqrt(pi H/(2 Gamma))).
Usage: python3 -B relaxed_check.py <out.json>     (~3 min)
"""
import json, os, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.optimize import root

sys.dont_write_bytecode = True
T = Path(__file__).resolve().parents[2]
os.environ.setdefault('SPY_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                      'scratchpad/review/ls_critic/build')
sys.path.insert(0, str(T / 'sp_yledger' / 'code'))
sys.path.insert(0, str(T / 'ls_ledger' / 'code'))
import ymodel as Y          # noqa: E402
import ls_project as P      # noqa: E402

H = 434867798206.6179
CH = 10.0
MY, HV, MU, KAP, LX, LY, LS, E = 1e12, 1e11, 5e11, 1e-3, 1e-3, 1.0, 1.0, 1e-7
MXL = 0.092408963737536
CANDS = {10: np.pi / 6, 17: 5 * np.pi / 24}


def grads(q, d):
    M, h, mu = MY / H, HV / H, MU / H
    x, yr, yi, s = q
    X = x + 0j; Yf = yr + 1j * yi; S = -s * np.exp(1j * d)
    gX = 2 * LX * abs(X)**2 * X + 3 * KAP * Yf * np.conj(X)**2 - CH * X
    gY = M * M * Yf + 2 * LY * abs(Yf)**2 * Yf + KAP * X**3 + 2 * mu * S * np.conj(Yf)
    gS = 2 * LS * (abs(S)**2 - h * h / 2) * S + mu * Yf * Yf
    return np.array([2 * gX.real, 2 * gY.real, 2 * gY.imag, 2 * np.real(gS * (-np.exp(-1j * d)))]), gX


def relaxed(d):
    m = Y.minimum(H, CH)
    q0 = np.array([m['x'], -m['y'], 0.0, m['s']])
    sol = root(lambda q: grads(q, d)[0], q0, method='hybr', options=dict(xtol=1e-14))
    q = sol.x
    J = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4); e[j] = 1e-6 * max(abs(q[j]), 1.0)
        J[:, j] = (grads(q + e, d)[0] - grads(q - e, d)[0]) / (2 * e[j])
    ev = np.linalg.eigvalsh(0.5 * (J + J.T))
    g, gX = grads(q, d)
    z0 = np.array([q[0] + 0j, q[1] + 1j * q[2], -q[3] * np.exp(1j * d)])
    return z0, dict(x=q[0], Y=[q[1], q[2]], s=q[3], grad_norm=float(np.linalg.norm(g)), hess_eigs=ev.tolist(),
                    torque_X=float(np.imag(np.conj(z0[0]) * gX)), success=bool(sol.success))


def ledger(NYH, NYL, Nsig, EX):
    mL = np.sqrt(MY**2 - np.sqrt(2) * MU * HV); mA = 6 * E * HV
    GL = 3 * KAP**2 * mL / (256 * np.pi**3)
    GA = E**2 * mA / (48 * np.pi) * (1 - 4 * MXL**2 / mA**2)**1.5
    dL = mL / H * np.sqrt(np.pi * H / (2 * GL)); dA = mA / H * np.sqrt(np.pi * H / (2 * GA))
    RY = (NYH + NYL) * dL / EX; RS = 2 * Nsig * dA / EX; RC = NYH * dA / EX
    return dict(R_Y=RY, R_Sigma=RS, R_cascade=RC, R_lower=RY + RS, R_central=RY + RS + RC)


def yield_over_YQ(NX):
    Tr = (90 * 2.435e18**2 * H**2 / (np.pi**2 * 106.75))**0.25
    return abs(NX) * H**3 / (2 * np.pi**2 / 45 * 106.75 * Tr**3) / 4.8294009076988175e-9


def job(cand):
    d = CANDS[cand]
    z0, info = relaxed(d)
    tend = 0.5 * 1000.0**2
    r = Y.run(H, z0, t_end=tend, rtol=1e-11, atol=1e-14, nsamp=60)
    w = r['windows']['600-1000']
    F = r['final']
    z = np.array([F[0] + 1j * F[1], F[2] + 1j * F[3], F[4] + 1j * F[5]])
    zd = np.array([F[6] + 1j * F[7], F[8] + 1j * F[9], F[10] + 1j * F[11]])
    sn = P.project(tend, z, zd, H)
    out = dict(candidate=cand, delta=d, prep=info,
               window=dict(eps=w['eps'], NX=w['NX'], EX=w['EX'], NYH=w['NYH'], NYL=w['NYL'], Nsig=w['Nsig']),
               snapshot={k: float(sn[k]) for k in ('NYH', 'NYL', 'Nsig', 'EX', 'NX')})
    out['ledger_window'] = ledger(w['NYH'], w['NYL'], w['Nsig'], w['EX'])
    out['ledger_snapshot'] = ledger(sn['NYH'], sn['NYL'], sn['Nsig'], sn['EX'])
    out['yield_over_YQ'] = yield_over_YQ(w['NX'])
    out['controls'] = r['end']
    return out


if __name__ == '__main__':
    Y.build()
    with Pool(2) as p:
        rows = p.map(job, list(CANDS))
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        print(json.dumps({k: r[k] for k in ('candidate', 'prep', 'window', 'ledger_window', 'ledger_snapshot',
                                            'yield_over_YQ')}, default=float))
