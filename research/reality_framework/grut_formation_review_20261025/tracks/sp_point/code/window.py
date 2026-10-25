"""sp_point: extent of the productive supplied-phase window around H* = 1.1211167e13 GeV (own code).

Lines (with exact tangent along the source phase, so each point carries its own chaos indicator ln|T| growth):
  line_H     : delta = pi/2, H/H* = geomspace(0.8, 1.25, 301)
  line_delta : H = H*, delta in (0.01, pi-0.01), 315 points
  fine_delta : H = H*, delta = pi/2 + [-0.02, 0.02], 201 points (no tangent)
Map (with tangent, focused on the band): H/H* = linspace(0.96, 1.04, 33) x delta = linspace(1.30, 1.85, 23).
Local ensembles (Codex's metric, sigma = 0.03 rad, N = 32): eps_net(H) on H/H* = linspace(0.9, 1.1, 21), delta = pi/2,
  and eps_net(delta0) at H* for delta0 = k pi/12, k = 1..11.

Usage: python3 window.py ../results/window.json [part ...]
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import spmodel as m

HS = m.H_STAR
C = 0.025


def job(spec):
    H, delta, tang = spec
    z = m.initial_z(H, delta)
    r = m.run(H, z, a_end=1000, cpar=C, windows=(1, 300, 600, 1000), ns=8,
              tangent=m.delta_tangent(H, delta) if tang else None)
    w = r['windows']['600-1000']
    hist = r['history']
    lnT = [h['lnT'] for h in hist]
    return dict(H=H, delta=delta, eps=w['epsX'], NX=w['NX'], Ncirc=w['NcircX'], JY=w['JY'],
                YX_over_YQ=w['YX_over_YQ'], eps300=r['windows']['300-600']['epsX'],
                deps=(w['deps'] if tang and not r['controls']['tangent_rescaled'] else None),
                lnT_growth=(max(lnT) - lnT[0]) if tang else None, steps=r['controls']['steps'],
                gauss=r['controls']['max_gauss_rel'])


def ens_stats(rows):
    nx = np.array([r['NX'] for r in rows]); nm = np.array([r['Ncirc'] for r in rows])
    mean = nx.mean() / nm.mean()
    sem = np.std(nx - mean * nm, ddof=1) / (np.sqrt(len(nx)) * nm.mean())
    return dict(eps_net=float(mean), sem=float(sem), frac_neg=float(np.mean(nx < 0)),
                pairs=float(1.5 * np.mean([r['JY'] for r in rows]) / abs(nx.mean())))


def main(out, parts):
    res = json.load(open(out)) if os.path.exists(out) else {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as pool:
        if 'lines' in parts:
            specs = [(HS * float(f), np.pi / 2, True) for f in np.geomspace(0.8, 1.25, 301)]
            res['line_H'] = list(pool.map(job, specs, chunksize=4))
            print('line_H done', time.time() - t0, flush=True)
            specs = [(HS, float(d), True) for d in np.linspace(0.01, np.pi - 0.01, 315)]
            res['line_delta'] = list(pool.map(job, specs, chunksize=4))
            print('line_delta done', time.time() - t0, flush=True)
            specs = [(HS, float(np.pi / 2 + d), False) for d in np.linspace(-0.02, 0.02, 201)]
            res['fine_delta'] = list(pool.map(job, specs, chunksize=4))
            print('fine_delta done', time.time() - t0, flush=True)
            json.dump(res, open(out, 'w'))
        if 'ens' in parts:
            rows = []
            for f in np.linspace(0.9, 1.1, 21):
                rng = np.random.default_rng(int(9000 + round(f * 1000)))
                d = rng.normal(0, 0.03, 32)
                rr = list(pool.map(job, [(HS * float(f), float(np.pi / 2 + x), False) for x in d]))
                rows.append(dict(H_factor=float(f), **ens_stats(rr), member_eps=[r['eps'] for r in rr]))
                print('ensH', f, rows[-1]['eps_net'], rows[-1]['sem'], flush=True)
            res['ens_H'] = rows
            rows = []
            for k in range(1, 12):
                rng = np.random.default_rng(9500 + k)
                d = rng.normal(0, 0.03, 32)
                rr = list(pool.map(job, [(HS, float(k * np.pi / 12 + x), False) for x in d]))
                rows.append(dict(delta0=k * np.pi / 12, **ens_stats(rr), member_eps=[r['eps'] for r in rr]))
                print('ensD', k, rows[-1]['eps_net'], rows[-1]['sem'], flush=True)
            res['ens_delta'] = rows
            json.dump(res, open(out, 'w'))
        if 'map' in parts:
            specs = [(HS * float(f), float(d), True) for f in np.linspace(0.96, 1.04, 33)
                     for d in np.linspace(1.30, 1.85, 23)]
            res['map'] = list(pool.map(job, specs, chunksize=8))
            print('map done', time.time() - t0, flush=True)
            json.dump(res, open(out, 'w'))
    res['scope'] = ('Supplied aligned-radii preparation, c_H=1, homogeneous k=0, Yoshida-8 c=%g, window 600-1000. '
                    'Not an inflationary distribution or a spatial calculation.' % C)
    json.dump(res, open(out, 'w'))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:] or ['lines', 'ens', 'map'])
