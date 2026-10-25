"""sp_point: epoch at which the X charge is set, versus the kappa-resonance estimate a_res = 3 omega_X H/M_Y (own code).

omega_X = 1.198 sqrt(lambda_X) (E_X/lambda_X)^(1/4) is the conformal frequency of the quartic X oscillation through
the origin (real orbit; 2 pi / (4 x 1.311029), with 1.311029 = int_0^1 du/sqrt(1-u^4)); [D] for a real orbit, [E] here; a25/a75: first a at which the running maximum of
|N_X|/N_circ reaches 25%/75% of its overall maximum. Also the Sigma vacuum-mass crossing omega_X/(sqrt(2 lambda_S) h/H).
Usage: python3 formation_epoch.py ../results/formation_epoch.json
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import spmodel as m


def job(spec):
    H, d = spec
    z = m.initial_z(H, d)
    r = m.run(H, z, a_end=1000, cpar=0.025, windows=(1, 600, 1000), ns=1500)
    h = r['history']
    a = np.array([x['a'] for x in h]); NX = np.array([x['NX'] for x in h]); EX = np.array([x['EX'] for x in h])
    Nc = 2 * np.sqrt(2 * m.LX) * (EX / (3 * m.LX))**0.75
    cm = np.maximum.accumulate(np.abs(NX) / Nc)
    a25 = float(a[np.argmax(cm >= 0.25 * cm[-1])]); a75 = float(a[np.argmax(cm >= 0.75 * cm[-1])])
    w = r['windows']['600-1000']
    wX = 1.198 * np.sqrt(m.LX) * (w['EX'] / m.LX)**0.25
    M, hh, mu = m.pars(H)
    return dict(H=H, delta=d, eps=w['epsX'], a25=a25, a75=a75, omega_X=wX, a_res_kappa=3 * wX / M,
                a_sigma=wX / (np.sqrt(2 * m.LS) * hh))


def main(out):
    specs = [(H, d) for H in [2e12, 3.3e12, 5.5e12, 8e12, m.H_STAR, 1.5e13, 2.2e13, 3.2e13, 4.7e13]
             for d in [np.pi / 2, np.pi / 4]]
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(job, specs))
    for r in rows:
        print('H=%.3e d=%.3f eps=%+.3f a25=%6.1f a75=%6.1f a_res=%6.1f a_sigma=%6.1f' %
              (r['H'], r['delta'], r['eps'], r['a25'], r['a75'], r['a_res_kappa'], r['a_sigma']))
    ratio = [r['a75'] / r['a_res_kappa'] for r in rows]
    json.dump(dict(rows=rows, a75_over_a_res=dict(min=min(ratio), max=max(ratio), median=float(np.median(ratio)))),
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
