#!/usr/bin/env python3
"""sp_phaseprep: transfer from a long-wavelength source-phase shift to the local yield (own code).

Needed for the isocurvature estimate S = T * dpsi_L with T = d ln|Y_X| / d psi_L:
 - local T: finite difference of the single release at delta = pi/2 +- 1e-7 (aligned radii, instant release);
 - smoothed T: net yield of an ensemble with short-scale Gaussian phase noise sigma_S (0.03 rad, Codex's tolerance
   level; and 0.29 rad, the free light mode's spread at x_c = 1/40), evaluated at delta0 = pi/2 + {-0.04..0.04} with
   common random numbers (the same standard-normal draws at every delta0), so that the slope is not swamped by
   sampling noise.

Usage: python3 transfer.py ../results/transfer.json
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ppi  # noqa: E402

XA, YA, SA = ppi.XA, ppi.YA, ppi.SA


def z_of(delta):
    return np.array([XA, -YA, -SA * np.exp(1j * delta)], complex)


def one(delta):
    r = ppi.run(z_of(delta))
    return dict(NX=r['NX'], Ncirc=r['Ncirc'], eps=r['eps'], JY=r['JY'])


def main(out):
    ppi.build()
    t0 = time.time()
    res = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        h = 1e-7
        r = list(pool.map(one, [np.pi / 2 - h, np.pi / 2, np.pi / 2 + h]))
        dNX = (r[2]['NX'] - r[0]['NX']) / (2 * h)
        deps = (r[2]['eps'] - r[0]['eps']) / (2 * h)
        res['local'] = dict(h=h, NX=r[1]['NX'], eps=r[1]['eps'], dNX_ddelta=dNX, deps_ddelta=deps,
                            T_local_lnNX=dNX / r[1]['NX'], T_local_lneps=deps / r[1]['eps'])
        print('local', res['local'], flush=True)
        offs = [-0.04, -0.02, 0.0, 0.02, 0.04]
        for sig, n, seed in ((0.03, 128, 777), (0.29, 128, 778)):
            xi = np.random.default_rng(seed).normal(size=n)
            rows = {}
            for o in offs:
                ds = np.pi / 2 + o + sig * xi
                rr = list(pool.map(one, ds))
                nx = np.array([q['NX'] for q in rr]); nc = np.array([q['Ncirc'] for q in rr])
                rows[o] = dict(NX_mean=float(nx.mean()), eps_net=float(nx.mean() / nc.mean()),
                               NX_sem=float(nx.std(ddof=1) / np.sqrt(n)), NX=nx.tolist())
                print(sig, o, rows[o]['eps_net'], rows[o]['NX_mean'], time.time() - t0, flush=True)
            # slope with common random numbers: per-member differences
            NXm = np.array([rows[o]['NX'] for o in offs])          # (5, n)
            A = np.vstack([np.ones(len(offs)), offs]).T
            coef = np.linalg.lstsq(A, NXm, rcond=None)[0]          # per-member linear fits (2, n)
            slope = coef[1].mean(); slope_sem = coef[1].std(ddof=1) / np.sqrt(n)
            c = rows[0.0]['NX_mean']
            res[f'smoothed_sigma{sig}'] = dict(
                sigma_S=sig, N=n, offsets=offs, eps_net=[rows[o]['eps_net'] for o in offs],
                NX_mean=[rows[o]['NX_mean'] for o in offs], dNXmean_ddelta0=float(slope),
                dNXmean_ddelta0_sem=float(slope_sem), T_smoothed=float(slope / c), T_smoothed_sem=float(slope_sem / abs(c)))
            print(res[f'smoothed_sigma{sig}'], flush=True)
    Path(out).write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
