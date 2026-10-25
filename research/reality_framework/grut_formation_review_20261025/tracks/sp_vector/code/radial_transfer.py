#!/usr/bin/env python3
"""Track sp_vector (Claude, own code): radial inflationary fluctuations and the late eaten phase.

Radial eigenmodes of the inflationary Hessian (gp_matching model, canonical real components rotated to the field
directions). For each radial eigenmode u (unit canonical vector over |X|, |Y|, |Sigma| directions), perturb the
initial moduli by |z_i| -> |z_i| + c u_i / sqrt2 at fixed phases (delta = pi/2), and compute:
  - the linear response d theta_S(inf)/dc and d eps/dc (symmetric differences, c = 1e-7 and 1e-6);
  - the BD amplitude per ln k at k_* at the end of inflation, (H/2pi) exp(-(3/2 - nu) N_*) for m^2 < 9/4 H^2;
  - the resulting late Goldstone angle delta theta_G(k_*) = |d theta_S/dc|/6 * delta c(k_*), against the
    0.226 rad per ln k criterion from vector_ledger.py;
  - an ensemble (N = 128) with Hubble-scale radial noise only (canonical sigma = H/2pi on each light radial mode),
    source phase fixed at pi/2: is the release decorrelated without any phase noise?
Usage: VREL=/path/vrel2 python3 radial_transfer.py ../results/radial_transfer.json
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'gp_matching' / 'code'))
import model as gpm

VREL = os.environ.get('VREL', 'vrel2')
HI = 1.1211167304108865e13
D0 = np.pi / 2
MODS = np.array([23.826064298969072, 1.8937619219396156, 0.43086809346988353])


def run(fx, fy, fs, delta=D0, rtol='1e-10'):
    p = subprocess.run([VREL, repr(float(delta)), repr(float(fs)), '1000', '600', rtol, repr(HI), '-', repr(float(fx)),
                        repr(float(fy))], capture_output=True, text=True)
    return json.loads(p.stdout)


def factors(u, c):
    d = c * u / np.sqrt(2)
    return 1 + d / MODS


def main(out):
    bg = json.loads((HERE.parent / 'results' / 'background.json').read_text())
    m = gpm.build(c_H=1.0, Hinf=HI)
    # canonical minimum point in the aligned gauge: X>0, Y<0, Sigma<0 (real); radial directions = +|z| directions
    z = np.array([np.sqrt(2) * MODS[0], 0, -np.sqrt(2) * MODS[1], 0, -np.sqrt(2) * MODS[2], 0])
    Hs = m['hess'](z)
    sgn = np.array([1, -1, -1])
    rad_idx = [0, 2, 4]
    Hrr = Hs[np.ix_(rad_idx, rad_idx)] * np.outer(sgn, sgn)   # Hessian along +|z_i| directions
    w, v = np.linalg.eigh(Hrr)
    N_star = bg['efolds_kstar_to_end']
    modes = []
    for j in range(3):
        u = v[:, j] * np.sign(v[np.argmax(abs(v[:, j])), j])
        nu2 = 9 / 4 - w[j]
        decay = np.exp(-(1.5 - np.sqrt(nu2)) * N_star) if nu2 > 0 else np.exp(-1.5 * N_star) * np.sqrt(1.5 / np.sqrt(w[j]))
        amp_k = (1 / (2 * np.pi)) * decay
        resp = {}
        for c in (1e-7, 1e-6):
            fp = factors(u, c); fm = factors(u, -c)
            with ThreadPoolExecutor(2) as ex:
                rp, rm = ex.map(lambda f: run(*f), [tuple(fp), tuple(fm)])
            resp[str(c)] = dict(dthS_dc=(rp['thS_inf'] - rm['thS_inf']) / (2 * c),
                                deps_dc=(rp['eps_ratio'] - rm['eps_ratio']) / (2 * c))
        dth = resp['1e-07']['dthS_dc']
        modes.append(dict(m2_over_H2=float(w[j]), u_canonical_XYS=[float(x) for x in u], response=resp,
                          dc_per_lnk_kstar=float(amp_k), dthetaG_kstar=float(abs(dth) / 6 * amp_k),
                          ratio_to_criterion=float(abs(dth) / 6 * amp_k / 0.2257)))
        print(modes[-1], flush=True)
    # radial-noise-only ensemble on the light radial modes
    rng = np.random.default_rng(2026102508)
    light = [j for j in range(3) if w[j] < 9 / 4]
    args = []
    for n in range(128):
        dz = np.zeros(3)
        for j in light:
            dz += rng.standard_normal() * (1 / (2 * np.pi)) * np.array(modes[j]['u_canonical_XYS'])
        f = 1 + dz / np.sqrt(2) / MODS
        args.append(tuple(f))
    with ThreadPoolExecutor(4) as ex:
        rows = list(ex.map(lambda f: run(*f), args))
    good = [r for r in rows if 'fail' not in r]
    arr = lambda k: np.array([r[k] for r in good])
    ens = dict(n=len(good), sigma_canonical=1 / (2 * np.pi), light_modes=light,
               eps_mean=float(arr('eps_ratio').mean()), eps_std=float(arr('eps_ratio').std(ddof=1)),
               NX_mean=float(arr('NX_mean').mean()), NS_std=float(arr('NS_mean').std(ddof=1)),
               thS_inf_std=float(arr('thS_inf').std(ddof=1)), minS_min=float(arr('minS_over_vac').min()))
    print(ens, flush=True)
    Path(out).write_text(json.dumps(dict(modes=modes, radial_noise_ensemble=ens,
                                         members=[dict(f=[r['xfac'], r['yfac'], r['sfac']], eps=r['eps_ratio'],
                                                       NS=r['NS_mean'], thS_inf=r['thS_inf']) for r in good]),
                                    indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
