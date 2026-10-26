#!/usr/bin/env python3
"""hw_verify (Claude, own code): the exact scaling symmetry (M_Y, h, mu, H_i) -> s (M_Y, h, mu, H_i).

Derivation [D]: with the fields measured in H_i units and time in 1/H_i, the release equations contain only
lx, ly, ls, kap, c_H, M/H, h/H, mu/H; the ledger deposits contain m_L/H and H/Gamma_L (Gamma_L ∝ m_L) and
m_A/H and H/Gamma_A (Gamma_A ∝ m_A up to the (1-4m_X^2/m_A^2)^{3/2} factor); so eps, N_i, E_X and every R are
invariant.  Yield = |N_X| H^3/s(T_R) with s ∝ (H M_red)^{3/2} gives yield ∝ H^{3/2}.  c6 = alpha6 M_red^2 with
alpha6 = mu kap^2 Sigma/(M^4 - 4 mu^2 Sigma^2) ∝ H^{-2}.  Numerical check: the GeV couplings are scaled and
re-divided by H (a non-trivial floating-point path), for s = 0.1, 0.30665 (yield-normalising), 3 and 10.
Usage: python3 -B scaling.py ../results/scaling.json    (~1 min)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hw

H = hw.H_CAND


def job(s):
    Hs = s * H
    cp = hw.heavy(Hs=Hs)
    z0, info = hw.prep_relaxed(cp, Hs, hw.CH, np.pi / 2)
    r = hw.run(cp, Hs, z0, a_end=100.0, windows=(1, 10, 30, 100))
    w = r['windows']['30-100']
    o = hw.summarize(cp, Hs, w, prep_info=info, extra=dict(s=s, couplings_GeV={k: cp[k] for k in ('M_GeV', 'h_GeV', 'mu_GeV', 'kap')}))
    o['masses']['threshold_factor_A'] = (1 - 4 * hw.MX_LIGHT**2 / o['masses']['mA']**2)**1.5
    return o


if __name__ == '__main__':
    ss = [1.0, 0.1, 5.889090163799063**(-2 / 3), 3.0, 10.0]
    with Pool(4) as p:
        rows = p.map(job, ss)
    base = rows[0]
    checks = []
    for r in rows[1:]:
        s = r['s']
        checks.append(dict(s=s, H_GeV=r['H_GeV'], d_eps=r['eps'] - base['eps'], d_NX_rel=r['NX'] / base['NX'] - 1,
                           d_EX_rel=r['EX'] / base['EX'] - 1, d_NYL_rel=r['NYL'] / base['NYL'] - 1,
                           d_NYH_rel=r['NYH'] / base['NYH'] - 1, d_Rcentral_rel=r['R_central'] / base['R_central'] - 1,
                           yield_ratio=r['yield_over_YQ'] / base['yield_over_YQ'], s_pow_1p5=s**1.5,
                           c6_ratio=r['c6'] / base['c6'], s_pow_m2=s**-2.0, yield_over_YQ=r['yield_over_YQ'], c6=r['c6'],
                           R_L_normalised=r['R_L']))
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows, checks=checks), indent=1, default=float))
    for c in checks:
        print(c)
