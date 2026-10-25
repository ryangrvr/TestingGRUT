"""sp_point: phase-noise ensembles at the leading point (own code).

A. Reproduce Codex's sigma=0.03, N=128, seed 2026102600 source-phase ensemble member by member.
B. Wider supplied source-phase spreads sigma = 0.1, 0.26, 0.5 rad (N=128 each).
   0.26 rad is the per-e-fold de Sitter kick of the light source mode, H/(2 pi sqrt2 |Sigma|) at this point;
   0.5 rad is the order of its equilibrium spread (estimate, see NOTES).
C. De Sitter (Starobinsky-Yokoyama) Gaussian draws of all six field components around the aligned inflationary
   minimum (gauge zero mode removed), N=256, and the same draws shifted so Sigma is centred on Codex's supplied
   Sigma = -|Sigma| e^{i pi/2}. Estimates only: Gaussian approximation, zero velocities.

Usage: python3 ensembles.py ../results/ensembles.json
"""
import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import spmodel as m

H = m.H_STAR
C = 0.0125
PKT = ('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/sp_point/pkt/'
       'research/reality_framework/grut_inflation_scan_20261025/ENSEMBLE_MEMBERS.csv')


def one(z):
    r = m.run(H, np.asarray(z, complex), a_end=1000, cpar=C, windows=(1, 300, 600, 1000), ns=4)
    w = r['windows']['600-1000']
    return dict(NX=w['NX'], Ncirc=w['NcircX'], epsX=w['epsX'], JY=w['JY'], NY=w['NY'],
                eps300=r['windows']['300-600']['epsX'], gauss=r['controls']['max_gauss_rel'])


def stats(rows):
    nx = np.array([r['NX'] for r in rows]); nm = np.array([r['Ncirc'] for r in rows])
    eps = nx / nm
    mean = nx.mean() / nm.mean()
    sem = np.std(nx - mean * nm, ddof=1) / (np.sqrt(len(nx)) * nm.mean())
    jy = np.array([r['JY'] for r in rows])
    return dict(N=len(rows), eps_net=float(mean), eps_net_sem=float(sem), eps_mean=float(np.mean([r['epsX'] for r in rows])),
                local_eps_rms=float(np.sqrt(np.mean(eps**2))), local_eps_std=float(eps.std(ddof=1)),
                frac_abs_ge_019=float(np.mean(np.abs(eps) >= 0.19)), frac_negative=float(np.mean(nx < 0)),
                YX_over_YQ_mean=float(nx.mean() * m.yield_factor(H) / m.YQ),
                pairs_per_net_charge=float(1.5 * jy.mean() / abs(nx.mean())))


def sdy_draws(n, seed, centre_sigma_on_supplied=False):
    """Gaussian de Sitter equilibrium P ~ exp(-8 pi^2 V/(3 H^4)) around the aligned inflationary minimum."""
    x, y, s = m.minimum(H, 1.0)
    z0 = np.array([x, -y, -s], complex)
    def Vr(u):
        z = z0 + (u[0::2] + 1j * u[1::2])
        return m.V(z, H, 1.0)
    e = 1e-4
    Hs = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            di = np.zeros(6); dj = np.zeros(6); di[i] = e; dj[j] = e
            Hs[i, j] = (Vr(di + dj) - Vr(di - dj) - Vr(-di + dj) + Vr(-di - dj)) / (4 * e * e)
    # canonical phi = sqrt2 (Re z, Im z): Hess_phi = Hess_u / 2
    Hp = Hs / 2
    w, U = np.linalg.eigh(Hp)
    keep = w > 1e-8 * w.max()
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n):
        g = rng.normal(size=keep.sum())
        phi = U[:, keep] @ (g * np.sqrt(3 / (8 * np.pi**2) / w[keep]))
        u = phi / np.sqrt(2)
        z = z0 + (u[0::2] + 1j * u[1::2])
        if centre_sigma_on_supplied:
            z[2] = z[2] - z0[2] + (-s * np.exp(1j * np.pi / 2))
        draws.append(z)
    sig = np.sqrt(3 / (8 * np.pi**2) / w[keep]) / np.sqrt(2)
    return draws, dict(eigen_mass2_over_H2=w.tolist(), per_mode_rms_field_over_H=sig.tolist())


def main(out, parts=('A', 'B', 'C')):
    import os
    res = json.load(open(out))['results'] if os.path.exists(out) else {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as pool:
        if 'A' in parts:
            rng = np.random.default_rng(2026102600)
            d = rng.normal(0, 0.03, 128)
            zs = [m.initial_z(H, np.pi / 2 + float(x)) for x in d]
            rows = list(pool.map(one, zs))
            codex = [r for r in csv.DictReader(open(PKT)) if r['ensemble_id'] == 'candidate_0_sigma003']
            ce = np.array([float(r['epsX']) for r in codex])
            me = np.array([r['epsX'] for r in rows])
            res['A_sigma0.03_seed2026102600'] = dict(stats=stats(rows), member_eps=me.tolist(),
                                                     codex_member_eps=ce.tolist(),
                                                     max_member_abs_diff=float(np.max(np.abs(me - ce))),
                                                     median_member_abs_diff=float(np.median(np.abs(me - ce))))
            print('A', res['A_sigma0.03_seed2026102600']['stats'], 'member diff max', np.max(np.abs(me - ce)),
                  time.time() - t0, flush=True)
        if 'B' in parts:
            for sig, seed in [(0.1, 7101), (0.26, 7102), (0.5, 7103)]:
                rng = np.random.default_rng(seed)
                d = rng.normal(0, sig, 128)
                rows = list(pool.map(one, [m.initial_z(H, np.pi / 2 + float(x)) for x in d]))
                res['B_sigma%g' % sig] = dict(stats=stats(rows), member_delta=d.tolist(),
                                              member_eps=[r['epsX'] for r in rows])
                print('B', sig, res['B_sigma%g' % sig]['stats'], time.time() - t0, flush=True)
        if 'C' in parts:
            for lab, flag in [('C_dS_aligned', False), ('C_dS_supplied_sigma_centre', True)]:
                draws, info = sdy_draws(256, 7201, flag)
                rows = list(pool.map(one, draws))
                sI = np.array([z[2].imag for z in draws])
                eps = np.array([r['epsX'] for r in rows])
                res[lab] = dict(stats=stats(rows), info=info, member_eps=eps.tolist(), member_ImSigma=sI.tolist(),
                                corr_eps_ImSigma=float(np.corrcoef(eps, sI)[0, 1]),
                                mean_eps_ImSigma_neg=float(eps[sI < 0].mean()),
                                mean_eps_ImSigma_pos=float(eps[sI > 0].mean()))
                print(lab, res[lab]['stats'], 'corr', res[lab]['corr_eps_ImSigma'], time.time() - t0, flush=True)
        if 'D' in parts:
            # which fluctuations wash the supplied charge out: Sigma-only vs X,Y-only de Sitter draws (same draws as C)
            draws, info = sdy_draws(256, 7201, True)
            x, y, s = m.minimum(H, 1.0)
            zsup = np.array([x, -y, -s * np.exp(1j * np.pi / 2)], complex)
            for lab, keep in [('D_sigma_only', [2]), ('D_XY_only', [0, 1])]:
                zz = []
                for z in draws:
                    q = zsup.copy()
                    for k in keep:
                        q[k] = z[k]
                    zz.append(q)
                rows = list(pool.map(one, zz[:128]))
                res[lab] = dict(stats=stats(rows), member_eps=[r['epsX'] for r in rows])
                print(lab, res[lab]['stats'], time.time() - t0, flush=True)
    json.dump(dict(H_GeV=H, integrator='Yoshida-8, c=%g' % C, results=res,
                   scope='Homogeneous k=0 releases; supplied or Gaussian de Sitter draws; not a spatial calculation.'),
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1], tuple(sys.argv[2:]) or ('A', 'B', 'C'))
