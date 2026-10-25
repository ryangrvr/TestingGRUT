#!/usr/bin/env python3
"""sp_yledger (Claude, own code): does a non-instant (inflaton matter-dominated) reheating dilute the 66-fold
yield while keeping the charge window?  Item (c).

Background (yrelease.c bg 1): inflaton matter era with perturbative decay Gamma_phi set by a stated T_R,
T_R = (90/(pi^2 g*))^(1/4) sqrt(Gamma_phi M_red), g* = 106.75; X Hubble mass removed at a = 1 as in Codex's
preparation.  bg 2: the Hubble mass persists as -c_H (rho_phi/3)|X|^2 until the inflaton decays.
Yield: Y_X = N_X H_i^3 / (s a^3) with s a^3 from the late radiation after the inflaton has decayed.
Usage: python3 md_release.py <out.json>     (~20 min on 4 cores)
"""
import json, sys, time
from multiprocessing import Pool
import numpy as np
from ymodel import minimum, run, entropy_inst, Gamma_of_TR, MRED, GSTAR, YQ

HLEAD = 1.1211167304108865e13
TR_STATED = entropy_inst(HLEAD)[1] / 65.98227759048056   # = 4.28e13 GeV: 1/66 dilution if N_X were unchanged


def one(args):
    H, delta, TR, bg, tend, amp = args
    x, y, s = amp
    z0 = np.array([x, -y, -s * np.exp(1j * delta)])
    G = Gamma_of_TR(TR) / H
    t0 = time.time()
    try:
        r = run(H, z0, t_end=tend, bg=bg, gphi=G, rtol=1e-10, atol=1e-13, nsamp=10)
    except Exception as exc:
        return dict(H_GeV=H, delta=delta, T_R_GeV=TR, bg=bg, status='failed %r' % exc)
    F = r['final']
    Ta = (30 * F[14] * MRED**2 * H**2 / (np.pi**2 * GSTAR))**0.25
    S = 2 * np.pi**2 / 45 * GSTAR * Ta**3
    out = dict(H_GeV=H, delta=delta, T_R_GeV=TR, bg=bg, Gamma_phi_over_H=G,
               T_R_over_instant=TR / entropy_inst(H)[1], runtime_s=time.time() - t0, controls=r['end'])
    for k in ('300-600', '600-1000'):
        w = r['windows'].get(k)
        if w:
            out[k] = dict(eps=w['eps'], NX=w['NX'], EX=w['EX'], NY_modes=w['NYH'] + w['NYL'], NYH=w['NYH'],
                          Nsig=w['Nsig'], YX_over_YQ=w['NX'] * H**3 / S / YQ)
    return out


def main():
    jobs = []
    # (1) the leading RD point under the stated T_R, both Hubble-mass prescriptions, 11 phases
    for bg in (1, 2):
        for j in range(1, 12):
            jobs.append((HLEAD, j * np.pi / 12, TR_STATED, bg, 1.5e6))
    # (2) coarse H scan at the stated T_R and at T_R = 1.5e14 GeV (bg 1), 11 phases each
    for TR in (TR_STATED, 1.5e14):
        for H in (3e12, 6e12, 2e13, 3e13, 4.7e13):
            for j in range(1, 12):
                jobs.append((H, j * np.pi / 12, TR, 1, 1.5e6))
    amps = {}
    for H in sorted({j[0] for j in jobs}):
        m = minimum(H, 1.0); amps[H] = (m['x'], m['y'], m['s'])
    jobs = [j + (amps[j[0]],) for j in jobs]
    with Pool(4) as p:
        rows = p.map(one, jobs)
    failed = [r for r in rows if 'status' in r]
    rows = [r for r in rows if 'status' not in r]
    best = max(rows, key=lambda r: abs(r['600-1000']['eps']))
    ok = [r for r in rows if abs(r['600-1000']['eps']) >= 0.19 and abs(r['600-1000']['YX_over_YQ']) >= 1]
    res = dict(T_R_stated_GeV=TR_STATED, n_runs=len(rows), failed=failed, minima_over_H={f'{k:.6g}': v for k, v in amps.items()}, rows=rows,
               max_abs_eps_row=best, rows_with_eps_ge_019_and_yield_ge_YQ=ok,
               scope='Homogeneous k=0 release; coarse grid (H x 11 phases); a finite scan is not an exclusion.')
    json.dump(res, open(sys.argv[1], 'w'), indent=1, default=float)
    for r in rows:
        w = r['600-1000']
        print(f"bg{r['bg']} H={r['H_GeV']:.3g} TR={r['T_R_GeV']:.3g} d={r['delta']:.3f} eps={w['eps']:+.4f} "
              f"Y/YQ={w['YX_over_YQ']:+.4g} NY={w['NY_modes']:.3g} Nsig={w['Nsig']:.3g}")
    print('max |eps|', best['600-1000']['eps'], best['H_GeV'], best['delta'], best['T_R_GeV'], 'passing rows', len(ok))


if __name__ == '__main__':
    main()
