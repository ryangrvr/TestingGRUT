#!/usr/bin/env python3
"""Track sp_critic (Claude, own driver; reuses our sp_vector release code vrel.c and its P0 transfer).

Cross-track question: does the best hypothetical phase preparation of sp_phaseprep (DRT-type term with c_S = 1,
light phase mass 1.153 H, f_eff = 0.635 H, residual source-phase spread 0.026 rad at x_c = 1/40) also clear
sp_vector's release-induced Goldstone piece B1 at the GMR peak scale k_*?

1. Noise-averaged response d<theta_S(inf)>/dpsi0 under phase-only Gaussian noise sigma_S = 0.026 rad, at
   psi0 = pi/2 + {-0.04, -0.02, 0, 0.02, 0.04}, N = 256 common draws per offset (vrel, rtol 1e-10, a_end 1000).
2. B1 vectors: N_A = int dlnk [dN_A^{P0}/dlnk] (slope/6)^2 P_psi(k) / (H/(2 pi F_vac))^2 (sp_vector conversion),
   with the PPI light-phase power P_psi(k) = [H/(2 pi f_eff)]^2 C_nu^2 (k/a_end H)^{3 - 2 nu} at the end of inflation.
   Phase only; the PPI radial modes (1.40, 1.45 H) are not included.
Usage: python3 -B b1_ppi_slope.py <sp_vector_dir> <out.json>    (~8 min on 4 cores)
"""
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.special import gamma as G

SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/sp_critic/bin')
HI = '1.1211167304108865e13'
D0 = math.pi / 2
SIG_S = 0.026
OFFS = [-0.04, -0.02, 0.0, 0.02, 0.04]
NMEM = 256
M_LIGHT, F_EFF = 1.152856997633605, 0.6349408014962153     # sp_phaseprep ppi.json, tuned_cS1.0
NX_PPI = 0.2416 / 0.24066 * 238.65                            # net charge at eps_net = -0.2416 (sp_phaseprep)
ALLOW = 0.0049296413205678555


def build(vdir):
    SCR.mkdir(parents=True, exist_ok=True)
    exe = SCR / 'vrel'
    subprocess.run(['gcc', '-O3', '-o', str(exe), str(Path(vdir) / 'code' / 'vrel.c'), '-lm'], check=True)
    return str(exe)


def run(exe, delta):
    p = subprocess.run([exe, repr(float(delta)), '1.0', '1000', '600', '1e-10', HI], capture_output=True, text=True)
    return json.loads(p.stdout)


def main(vdir, out):
    exe = build(vdir)
    xi = np.random.default_rng(2026102531).standard_normal(NMEM)
    pts = []
    with ThreadPoolExecutor(4) as ex:
        for o in OFFS:
            rows = [r for r in ex.map(lambda d: run(exe, d), D0 + o + SIG_S * xi) if 'fail' not in r]
            th = np.array([r['thS_inf'] for r in rows]); ep = np.array([r['eps_ratio'] for r in rows])
            nx = np.array([r['NX_mean'] for r in rows])
            pts.append(dict(offset=o, n=len(rows), thS_mean=float(th.mean()), thS_sem=float(th.std(ddof=1) / math.sqrt(len(th))),
                            thS_std=float(th.std(ddof=1)), eps_mean=float(ep.mean()),
                            eps_sem=float(ep.std(ddof=1) / math.sqrt(len(ep))), NX_mean=float(nx.mean()),
                            thS=th.tolist(), eps=ep.tolist()))
            print({k: v for k, v in pts[-1].items() if k not in ('thS', 'eps')}, flush=True)
    x = np.array(OFFS)
    # weighted LSQ (points treated as independent) and a paired-difference estimate using the common draws
    def fit(key, sem):
        y = np.array([p[key] for p in pts]); w = 1 / np.array([p[sem] for p in pts])**2
        A = np.vstack([np.ones_like(x), x]).T
        cov = np.linalg.inv(A.T @ (A * w[:, None])); b = cov @ (A.T @ (w * y))
        return float(b[1]), float(math.sqrt(cov[1, 1]))
    sl_th, se_th = fit('thS_mean', 'thS_sem')
    sl_ep, se_ep = fit('eps_mean', 'eps_sem')
    TH = np.array([p['thS'] for p in pts]) if all(p['n'] == NMEM for p in pts) else None
    paired = None
    if TH is not None:
        per = np.array([np.polyfit(x, TH[:, j], 1)[0] for j in range(NMEM)])
        paired = dict(slope=float(per.mean()), se=float(per.std(ddof=1) / math.sqrt(NMEM)))
    # ---- B1 conversion (sp_vector P0 transfer)
    pm = json.loads((Path(vdir) / 'results' / 'proca_modes.json').read_text())
    modes = pm['profiles']['P0']['modes']
    ks = np.array([m['k'] for m in modes]); dN0 = ks**3 * np.array([m['n_k'] for m in modes]) / (2 * math.pi**2)
    PG = (1 / (2 * math.pi * 6e11 / float(HI)))**2
    nu = math.sqrt(2.25 - M_LIGHT**2)
    Cnu = 2**(nu - 1.5) * G(nu) / G(1.5)
    amp = 1 / (2 * math.pi * F_EFF)
    Ppsi = (amp * Cnu)**2 * ks**(3 - 2 * nu)
    kstar = pm['background']['k_star_over_Hi']
    res = dict(method=__doc__, sigma_S=SIG_S, offsets=OFFS, N=NMEM, points=pts,
               slope_thS=sl_th, slope_thS_se=se_th, slope_thS_paired=paired, slope_eps=sl_ep, slope_eps_se=se_ep,
               nu=nu, C_nu=Cnu, dpsi_per_lnk_at_kstar=float(amp * Cnu * kstar**(1.5 - nu)), NX_PPI=NX_PPI, B1=[])
    lk = np.log(ks)
    for lab, sl in (('central', abs(sl_th)), ('minus_2se', max(abs(sl_th) - 2 * se_th, 0.0)),
                    ('plus_2se', abs(sl_th) + 2 * se_th)):
        for kmax in (1 / 3000.0, 1 / 1000.0):
            sel = ks <= kmax
            NA = float(np.trapezoid(dN0[sel] * (sl / 6)**2 * Ppsi[sel] / PG, lk[sel]))
            res['B1'].append(dict(which=lab, slope=sl, kmax=kmax, dthetaG_kstar=sl / 6 * res['dpsi_per_lnk_at_kstar'],
                                  N_A=NA, pairs_per_net=NA / NX_PPI, ratio_to_allowance=NA / NX_PPI / ALLOW))
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k not in ('points', 'method')}, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
