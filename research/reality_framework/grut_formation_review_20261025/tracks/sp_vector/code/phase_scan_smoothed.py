#!/usr/bin/env python3
"""Track sp_vector (Claude, own code): large-scale response of the noise-averaged late Sigma phase.

For background source phases psi0 = pi/2 + {-0.6, -0.3, 0, 0.3, 0.6} and small-scale Gaussian noise of 0.489 rad
(the psi2 power above k_* at the end of inflation, from vmodel.py), N = 128 members each (common random numbers),
compute <theta_S(inf)>, <N_S>, <eps_X> and fit the slope d<theta_S(inf)>/dpsi0 with its standard error.
Usage: VREL=/path/vrel python3 phase_scan_smoothed.py ../results/phase_scan_smoothed.json
"""
import json
import sys
from pathlib import Path

import numpy as np

from phase_transfer import prun, D0


def main(out):
    xi = np.random.default_rng(2026102507).standard_normal(128)
    offs = [-0.6, -0.3, 0.0, 0.3, 0.6]
    pts = []
    for o in offs:
        rows = [r for r in prun([(D0 + o + 0.489 * x,) for x in xi]) if 'fail' not in r]
        th = np.array([r['thS_inf'] for r in rows]); ns = np.array([r['NS_mean'] for r in rows])
        ep = np.array([r['eps_ratio'] for r in rows]); nx = np.array([r['NX_mean'] for r in rows])
        pts.append(dict(offset=o, n=len(rows), thS_mean=float(th.mean()), thS_sem=float(th.std(ddof=1) / np.sqrt(len(th))),
                        NS_mean=float(ns.mean()), eps_mean=float(ep.mean()), eps_sem=float(ep.std(ddof=1) / np.sqrt(len(ep))),
                        NX_mean=float(nx.mean())))
        print(pts[-1], flush=True)
    x = np.array(offs); y = np.array([p['thS_mean'] for p in pts]); w = 1 / np.array([p['thS_sem'] for p in pts])**2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    b = cov @ (A.T @ (w * y))
    ye = np.array([p['eps_mean'] for p in pts]); we = 1 / np.array([p['eps_sem'] for p in pts])**2
    cove = np.linalg.inv(A.T @ (A * we[:, None])); be = cove @ (A.T @ (we * ye))
    res = dict(points=pts, slope_thS=float(b[1]), slope_thS_se=float(np.sqrt(cov[1, 1])),
               slope_eps=float(be[1]), slope_eps_se=float(np.sqrt(cove[1, 1])),
               thS_per_eps_nominal=-2310.74 / -0.240661,
               note='members are correlated across offsets only through common noise draws; SEs treat points as independent')
    print(json.dumps({k: v for k, v in res.items() if k != 'points'}, indent=1))
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
