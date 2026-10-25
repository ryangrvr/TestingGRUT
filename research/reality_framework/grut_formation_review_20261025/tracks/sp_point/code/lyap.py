"""sp_point: finite-time sensitivity / Lyapunov growth of the release versus H_inf (own code, exact tangent map).

For each H (c_H = 1, Codex's supplied preparation at delta) integrate the variational equations along the
source-phase direction to a = 1000 and record
  S = |d eps_X / d delta| (window 600-1000 mean), ln|T(a)| (phase-space tangent norm) history, and the epoch
  a_90 by which ln|T| has completed 90% of its growth up to its running maximum.
Mechanism diagnostics (clearly not the stated model): rescale one of M_Y, mu, h at fixed H.

Usage: python3 lyap.py ../results/lyap.json
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import spmodel as m


def job(spec):
    H, delta, scale, aend = spec
    m.SCALE.update(scale)
    z = m.initial_z(H, delta)
    t0 = time.time()
    r = m.run(H, z, a_end=aend, cpar=0.025, windows=(1, 30, 100, 300, 600, 1000), ns=240,
              tangent=m.delta_tangent(H, delta))
    w = r['windows']['600-1000']
    hist = r['history']
    a = np.array([h['a'] for h in hist]); lt = np.array([h['lnT'] for h in hist])
    lt0 = lt[0]
    runmax = np.maximum.accumulate(lt - lt0)
    tot = runmax[-1]
    a90 = float(a[np.argmax(runmax >= 0.9 * tot)]) if tot > 0 else 1.0
    a50 = float(a[np.argmax(runmax >= 0.5 * tot)]) if tot > 0 else 1.0
    Mh, hh, muh = m.pars(H)
    out = dict(H=H, delta=delta, scale=dict(scale), eps=w['epsX'], deps=w['deps'], S=abs(w['deps']),
               lnS=float(np.log(abs(w['deps']))), YX_over_YQ=w['YX_over_YQ'], JY=w['JY'], EX=w['EX'],
               lnT_growth=float(tot), lnT_final=float(lt[-1] - lt0), a50=a50, a90=a90,
               a_Y=1.0 / Mh, M_hat=Mh, mu_hat=muh, h_hat=hh, runtime=time.time() - t0,
               eps_300_600=r['windows']['300-600']['epsX'],
               hist_a=a[::4].tolist(), hist_lnT=(lt[::4] - lt0).tolist(),
               hist_deps=[h['deps'] for h in hist][::4], hist_eps=[h['epsX'] for h in hist][::4])
    m.SCALE.update(M=1.0, h=1.0, mu=1.0)
    return out


def main(out):
    Hs = np.geomspace(3e11, 4.7e13, 41)
    specs = [(float(H), np.pi / 2, {}, 1000.) for H in Hs]
    specs += [(float(H), np.pi / 4, {}, 1000.) for H in Hs[::2]]
    specs += [(m.H_STAR, np.pi / 2, {}, 1000.)]
    for f in [0.5, 0.71, 1.41, 2.0]:
        specs.append((m.H_STAR, np.pi / 2, dict(M=f), 1000.))
        specs.append((4.7e13, np.pi / 2, dict(M=f), 1000.))
    for f in [0.5, 2.0]:
        specs.append((m.H_STAR, np.pi / 2, dict(mu=f), 1000.))
        specs.append((m.H_STAR, np.pi / 2, dict(h=f), 1000.))
    specs.sort(key=lambda s: s[0] / s[2].get('M', 1.0))  # cheapest (high H_hat) first is irrelevant; keep order stable
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as pool:
        res = list(pool.map(job, specs))
    for r in res:
        print('H=%.4e d=%.3f sc=%s eps=%+.4f lnS=%6.2f growth=%6.2f a50=%7.1f a90=%7.1f a_Y=%6.1f t=%.0fs' %
              (r['H'], r['delta'], r['scale'], r['eps'], r['lnS'], r['lnT_growth'], r['a50'], r['a90'], r['a_Y'],
               r['runtime']), flush=True)
    json.dump(dict(rows=res, runtime=time.time() - t0,
                   scope='Exact discrete tangent of the Yoshida-8 map along the supplied source phase; homogeneous k=0.'),
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
