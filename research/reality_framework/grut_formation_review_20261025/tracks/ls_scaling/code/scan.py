#!/usr/bin/env python3
"""ls_scaling (Claude, own code): coupling-direction spot checks of epsilon, yield, debris R, R_L and c6.

Part A (stated regime): one-at-a-time coupling changes about the stated point at the candidates' H_i and c_H = 10,
  Codex's preparation (aligned radial minimum, then Sigma -> -|Sigma| e^{i delta}), 6 phases, a_end = 200
  (window a = 100-200; the candidate-10 control shows R within 0.03% of a = 1000), plus the Y-relaxed preparation
  at the base point.
Part B (heavy-mediator direction): M_Y/H = m, r = sqrt2 mu h/M_Y^2, h/H, kappa set by the EFT edge fraction
  eta/eta_edge (eta = alpha6 X0^2/lambda_X, eta_edge = r/12 at phi0 = pi/2), Y-relaxed preparation,
  a_end = 40 (window a = 30-40; numbers are frozen from a ~ 10, control at a_end = 100 in the log).
Usage: python3 -B scan.py <out.json> [A|B|AB]     (A ~20 min, B ~10 min on 4 cores)
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scal as S  # noqa: E402

H = S.H_CAND
CH = 10.0
LAM = 1e-3


def heavy(m=30.0, r=0.45, hH=40.0, frac=0.5, ly=1.0, Hs=H, e=1e-7):
    kap = m * LAM * np.sqrt(frac * 4 * (1 - r * r) / (12 * CH))
    M = m * Hs; h = hH * Hs; mu = r * M * M / (np.sqrt(2) * h)
    return dict(S.STATED, kap=kap, M_GeV=M, h_GeV=h, mu_GeV=mu, ly=ly, e=e)


def stated(**f):
    cp = dict(S.STATED)
    for k, v in f.items():
        cp[k] = cp[k] * v
    return cp


def job(spec):
    name, cp, Hs, deltas, prep, aend, win = spec
    out = []
    t0 = time.time()
    if prep == 'codex':
        ms = S.minima(cp, Hs, CH)
    for d in deltas:
        if prep == 'codex':
            m = ms[0]
            z0 = np.array([m['x'], -m['y'], -m['s'] * np.exp(1j * d)]); info = dict(x=m['x'], y=m['y'], s=m['s'],
                                                                                   n_minima=len(ms))
        else:
            z0, info = S.relaxed_initial(cp, Hs, CH, d)
        try:
            rr = S.run(cp, Hs, z0, a_end=aend, rtol=1e-10, atol=1e-13, nsamp=40)
        except Exception as exc:
            out.append(dict(name=name, delta=d, prep=prep, status=repr(exc))); continue
        w = rr['windows'][win]
        L = S.ledger(cp, Hs, w['NYH'], w['NYL'], w['Nsig'], w['EX'])
        eps = w['eps']
        RL = S.R_L(eps)
        out.append(dict(name=name, delta=d, prep=prep, H_GeV=Hs, eps=eps, NX=w['NX'], EX=w['EX'], NYH=w['NYH'],
                        NYL=w['NYL'], Nsig=w['Nsig'], R_Y=L['R_Y'], R_Sigma=L['R_Sigma'], R_cascade=L['R_cascade'],
                        R_lower=L['R_lower'], R_central=L['R_central'], R_L=RL,
                        RL_note=None if RL is not None else ('eps above table: R_L > 172' if abs(eps) > 0.24
                                                             else 'eps below table: R_L < 30'),
                        central_over_RL=(L['R_central'] / RL) if RL else None,
                        yield_over_YQ=S.yield_over_YQ(w['NX'], Hs), c6=S.c6(cp)['c6'], kappa=cp['kap'],
                        kappa_max=S.kappa_max(cp), masses=L['masses'], flags=L['flags'], prep_info=info,
                        controls=rr['end'],
                        couplings_over_H=dict(M=cp['M_GeV'] / Hs, h=cp['h_GeV'] / Hs, mu=cp['mu_GeV'] / Hs)))
    print(name, prep, 'done %.0fs' % (time.time() - t0), flush=True)
    return out


def specs_A():
    D6 = [(k + 0.5) * np.pi / 6 for k in range(6)]
    A = [('A_base', stated()), ('A_M_x2', stated(M_GeV=2)), ('A_M_x4', stated(M_GeV=4)),
         ('A_kap_x0.5', stated(kap=0.5)), ('A_kap_x3', stated(kap=3)),
         ('A_mu_x0.5', stated(mu_GeV=0.5)), ('A_mu_x2', stated(mu_GeV=2)),
         ('A_h_x3_mu_x1/3', stated(h_GeV=3, mu_GeV=1 / 3)), ('A_h_x3', stated(h_GeV=3)),
         ('A_lamY_x0.1', stated(ly=0.1)), ('A_lamY_x10', stated(ly=10))]
    sp = [(n, cp, H, D6, 'codex', 200.0, '100-300') for n, cp in A]
    sp.append(('A_base', stated(), H, D6, 'relaxed', 200.0, '100-300'))
    return sp


def specs_B():
    P2, P4 = np.pi / 2, np.pi / 4
    sp = []
    for m in (10.0, 15.0, 20.0, 30.0, 40.0):
        sp.append((f'B_m{m:g}', heavy(m=m), H, [P2, P4], 'relaxed', 40.0, '30-100'))
    for hH in (1.0, 3.0, 10.0, 20.0):
        sp.append((f'B_h{hH:g}', heavy(hH=hH), H, [P2], 'relaxed', 40.0, '30-100'))
    for fr in (0.25, 0.75, 0.9):
        sp.append((f'B_frac{fr:g}', heavy(frac=fr), H, [P2], 'relaxed', 40.0, '30-100'))
    for r in (0.0707, 0.2, 0.7):
        sp.append((f'B_r{r:g}', heavy(r=r), H, [P2], 'relaxed', 40.0, '30-100'))
    for ly in (0.1, 10.0):
        sp.append((f'B_lamY{ly:g}', heavy(ly=ly), H, [P2], 'relaxed', 40.0, '30-100'))
    sp.append(('B_phase', heavy(), H, [(k + 0.5) * np.pi / 6 for k in range(6)], 'relaxed', 40.0, '30-100'))
    sp.append(('B_base_codexprep', heavy(), H, [P2], 'codex', 40.0, '30-100'))
    sp.append(('B_base_a100', heavy(), H, [P2], 'relaxed', 100.0, '30-100'))
    s = 5.889032093211876**(-2 / 3)   # rescale every scale by s so the raw yield is ~ Y_Q (scaling symmetry check)
    sp.append(('B_base_Hscaled', heavy(Hs=s * H), s * H, [P2], 'relaxed', 40.0, '30-100'))
    return sp


if __name__ == '__main__':
    out = Path(sys.argv[1]); which = sys.argv[2] if len(sys.argv) > 2 else 'AB'
    sp = (specs_A() if 'A' in which else []) + (specs_B() if 'B' in which else [])
    sp.sort(key=lambda s: -(s[5]**2) * (s[1]['M_GeV'] / s[2]) * len(s[3]))   # longest first
    rows = []
    with Pool(4) as p:
        for res in p.imap_unordered(job, sp):
            rows += res
            out.write_text(json.dumps(dict(method=__doc__, H_GeV=H, cH=CH, rows=rows), indent=1, default=float))
    print('rows', len(rows))
