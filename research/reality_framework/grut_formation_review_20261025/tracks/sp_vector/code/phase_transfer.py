#!/usr/bin/env python3
"""Track sp_vector (Claude, own code): separate-universe transfer from the source phase to the late eaten phase.

Runs our vrel release (code/vrel.c, compiled to $VREL) for:
  (1) linear response: symmetric differences of theta_S(inf), N_S, eps_X in the supplied Sigma phase delta and in |Sigma|;
  (2) Gaussian ensembles of the Sigma phase about pi/2 (as in Codex's ensembles) for several widths sigma;
  (3) the smoothed large-scale response d<theta_S(inf)>/d psi0 under small-scale noise of 0.489 rad
      (common random numbers at psi0 = pi/2 +- dpsi).
theta_S(inf) = theta_S(a_end) - <N_S>/(h^2 a_end) is the late unwrapped Sigma phase; the gauge-invariant late Goldstone
angle is theta_S(inf)/6 (charge 6), compared with its inflationary value theta_G0 (vrel output).
Usage: VREL=/path/vrel python3 phase_transfer.py ../results/phase_transfer.json
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

VREL = os.environ.get('VREL', 'vrel')
HI = '1.1211167304108865e13'
D0 = np.pi / 2
RTOL = '1e-10'


def run(delta, sfac=1.0, aend=1000, aw0=600, rtol=RTOL):
    p = subprocess.run([VREL, repr(float(delta)), repr(float(sfac)), str(aend), str(aw0), rtol, HI],
                       capture_output=True, text=True)
    return json.loads(p.stdout)


def prun(args, workers=4):
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(lambda a: run(*a), args))


def stats(rows, key):
    v = np.array([r[key] for r in rows if 'fail' not in r])
    return dict(mean=float(v.mean()), std=float(v.std(ddof=1)), sem=float(v.std(ddof=1) / np.sqrt(len(v))), n=len(v))


def main(out):
    res = dict(method=__doc__, rtol=RTOL)
    # (1) linear response
    lin = []
    for dd in [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]:
        rp, rm = prun([(D0 + dd,), (D0 - dd,)])
        lin.append(dict(dd=dd, dthS_inf=(rp['thS_inf'] - rm['thS_inf']) / (2 * dd),
                        dNS=(rp['NS_mean'] - rm['NS_mean']) / (2 * dd),
                        deps=(rp['eps_ratio'] - rm['eps_ratio']) / (2 * dd),
                        minS=[rp['minS_over_vac'], rm['minS_over_vac']]))
        print('lin', lin[-1], flush=True)
    res['linear_phase'] = lin
    linr = []
    for ds in [1e-9, 1e-7, 1e-5]:
        rp, rm = prun([(D0, 1 + ds), (D0, 1 - ds)])
        linr.append(dict(ds=ds, dthS_inf=(rp['thS_inf'] - rm['thS_inf']) / (2 * ds),
                         dNS=(rp['NS_mean'] - rm['NS_mean']) / (2 * ds),
                         deps=(rp['eps_ratio'] - rm['eps_ratio']) / (2 * ds)))
        print('linr', linr[-1], flush=True)
    res['linear_radial_Sigma'] = linr
    # (2) ensembles
    ens = {}
    rng = np.random.default_rng(2026102505)
    for sig in [1e-6, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.26, 0.489]:
        deltas = D0 + sig * rng.standard_normal(128)
        rows = prun([(d,) for d in deltas])
        good = [r for r in rows if 'fail' not in r]
        thG = np.array([r['thS_inf'] / 6 for r in good])
        e = dict(sigma=sig, n=len(good), fails=len(rows) - len(good),
                 eps=stats(good, 'eps_ratio'), NX=stats(good, 'NX_mean'), NS=stats(good, 'NS_mean'),
                 thS_inf=stats(good, 'thS_inf'), thG_late_std=float(thG.std(ddof=1)),
                 minS_min=float(min(r['minS_over_vac'] for r in good)),
                 frac_minS_below_0p01=float(np.mean([r['minS_over_vac'] < 0.01 for r in good])),
                 members=[dict(delta=r['delta'], eps=r['eps_ratio'], NX=r['NX_mean'], NS=r['NS_mean'],
                               thS_inf=r['thS_inf'], minS=r['minS_over_vac']) for r in good])
        ens[str(sig)] = e
        print('ens', sig, e['eps'], e['NS'], e['thS_inf'], e['minS_min'], flush=True)
    res['ensembles'] = ens
    # (3) smoothed large-scale response with common random numbers
    sm = []
    xi = np.random.default_rng(2026102506).standard_normal(256)
    for dpsi in [0.05, 0.15]:
        rp = prun([(D0 + dpsi + 0.489 * x,) for x in xi])
        rm = prun([(D0 - dpsi + 0.489 * x,) for x in xi])
        dif = np.array([(a['thS_inf'] - b['thS_inf']) for a, b in zip(rp, rm) if 'fail' not in a and 'fail' not in b])
        difN = np.array([(a['NS_mean'] - b['NS_mean']) for a, b in zip(rp, rm) if 'fail' not in a and 'fail' not in b])
        difE = np.array([(a['eps_ratio'] - b['eps_ratio']) for a, b in zip(rp, rm) if 'fail' not in a and 'fail' not in b])
        sm.append(dict(dpsi=dpsi, n=len(dif), dthS_dpsi=float(dif.mean() / (2 * dpsi)),
                       dthS_dpsi_sem=float(dif.std(ddof=1) / np.sqrt(len(dif)) / (2 * dpsi)),
                       dNS_dpsi=float(difN.mean() / (2 * dpsi)), dNS_dpsi_sem=float(difN.std(ddof=1) / np.sqrt(len(dif)) / (2 * dpsi)),
                       deps_dpsi=float(difE.mean() / (2 * dpsi)), deps_dpsi_sem=float(difE.std(ddof=1) / np.sqrt(len(dif)) / (2 * dpsi)),
                       mean_thS_plus=float(np.mean([a['thS_inf'] for a in rp])), mean_thS_minus=float(np.mean([b['thS_inf'] for b in rm])),
                       mean_eps_plus=float(np.mean([a['eps_ratio'] for a in rp])), mean_eps_minus=float(np.mean([b['eps_ratio'] for b in rm]))))
        print('smooth', sm[-1], flush=True)
    res['smoothed_response'] = sm
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
