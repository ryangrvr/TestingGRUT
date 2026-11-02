#!/usr/bin/env python3
"""Track inf_higgs: danger exponent with the benchmarks' actual plasma timing, the R-even portal, and controls.

  python3 -B ih_final.py plasmaD    # D(M_t) with the S1 hard-gas induced mass m_ind^2(tau) of ih_plasma.py in m_c^2
  python3 -B ih_final.py portal     # mode runs with (lam_IH/2) phi_I^2 |H|^2 at lam_min (m_H^2 = lam phi_I^2/2)
  python3 -B ih_final.py controls   # Hartree at the critical M_t; Nb0 = 3 vs 5 start; n = 400 modes; rtol 1e-10
Writes ../results/final_plasmaD.json, final_portal.json, final_controls.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np

import ih_common as C
import ih_modes as M
import ih_run as IR

MT_EXT = (172.57, 172.86, 173.15, 173.45, 173.75, 174.05, 174.35, 174.6)


def crit(Ds, Mts):
    return IR.crit_Mt(Ds, list(Mts))


def cmd_plasmaD():
    full = C.SP / 'plasma_full.json'                # full induced-mass curves (scratchpad); thinned copy in results
    P = json.load(open(full if full.exists() else C.RES / 'plasma.json'))
    out = {}
    for name in C.BENCH:
        H = C.load_bench(name)['H_GeV']
        xi = IR.smtab(172.57).xi(H, 0.1)
        R = IR.run_cached(name, 'MD', xi, 0)
        Nob = IR.N_obs(name, 'MD', 0)
        tau = R['rows'][:, 0]
        o = {}
        grav = [IR.evalD(R, xi, H, Mt, Ns=(60, Nob))['best'][60]['D'] for Mt in MT_EXT]
        o['gravitational_only'] = dict(Mt=MT_EXT, D60=grav)
        rows = P[name]['rows']
        sel = list(rows)
        for k in sel:
            c = rows[k]['curve']
            mp2 = np.interp(tau, c['t'], c['m2'], left=0.0, right=c['m2'][-1])
            Ds, tcs = [], []
            for Mt in MT_EXT:
                d = IR.evalD(R, xi, H, Mt, Ns=(60, Nob), mp2=mp2)
                Ds.append(d['best'][60]['D']); tcs.append(d['best'][60].get('tau_c'))
            o[k] = dict(D60=Ds, tau_c=tcs, Mt_crit=crit(Ds, MT_EXT), tau_m2_1em3=rows[k]['tau_m2_1em3'])
            print(name, k, 'D60', [round(x, 1) if np.isfinite(x) else x for x in Ds], 'crit', o[k]['Mt_crit'], flush=True)
        o['Mt_crit_gravitational'] = crit(grav, MT_EXT)
        out[name] = o
        C.jdump(out, C.RES / 'final_plasmaD.json')


def cmd_portal():
    P = json.load(open(C.RES / 'plasma.json'))
    out = {}
    for name, hist in (('HT6t', 'MD'), ('T1', 'MD')):   # homogeneous: keeps <phi_I^2> (the w = 0 fluid would drop it)
        B = C.load_bench(name)
        H = B['H_GeV']
        xi = IR.smtab(172.57).xi(H, 0.1)
        lam = P[name]['lambda_IH']['lam_min']
        o = {}
        for f in (1.0, 0.1):
            t0 = time.time()
            R = IR.run_cached(name, hist, xi, 0, lamIH=lam * f)
            Ds = [IR.evalD(R, xi, H, Mt, Ns=(60,))['best'][60]['D'] for Mt in MT_EXT]
            rows = R['rows']
            te = R['res']['tau_endinf']
            post = rows[:, 0] > te
            o[f'lam={f}*lam_min'] = dict(lamIH=lam * f, D60=Ds, sigma2_endinf=R['res']['sigma2_L_endinf'],
                                       sigma2_tau0=R['res']['sigma2_tau0'],
                                       min_m2x_over_H2_post=float(np.min(R['m2xbar'][post] / rows[post, 2] ** 2)),
                                       seconds=time.time() - t0)
            print(name, f, o[f'lam={f}*lam_min'], flush=True)
        out[name] = o
        C.jdump(out, C.RES / 'final_portal.json')


def cmd_controls():
    out = {}
    name, hist = 'HT6t', 'MD'
    B = C.load_bench(name)
    H = B['H_GeV']
    xi = IR.smtab(172.57).xi(H, 0.1)
    base = IR.run_cached(name, hist, xi, 0)
    Mts = (172.86, 172.92, 173.0)
    ref = {f'{m}': IR.evalD(base, xi, H, m, Ns=(60,))['best'][60]['D'] for m in Mts}
    out['reference'] = ref
    vari = {}
    t0 = time.time()
    R = IR.run_cached(name, hist, xi, 0, sm_mt=172.92)
    vari['Hartree at M_t = 172.92'] = {f'{m}': IR.evalD(R, xi, H, m, Ns=(60,))['best'][60]['D'] for m in Mts}
    R = IR.run_cached(name, hist, xi, 0, Nb0=3.0)
    vari['start Nb0 = 3'] = {f'{m}': IR.evalD(R, xi, H, m, Ns=(60,))['best'][60]['D'] for m in Mts}
    R = IR.run_cached(name, hist, xi, 0, n=400)
    vari['n = 400 modes'] = {f'{m}': IR.evalD(R, xi, H, m, Ns=(60,))['best'][60]['D'] for m in Mts}
    R = IR.run_cached(name, hist, xi, 0, Kmax=15.0)
    vari['Kmax = 15'] = {f'{m}': IR.evalD(R, xi, H, m, Ns=(60,))['best'][60]['D'] for m in Mts}
    # rtol: separate run (not cached under the same key)
    bg = IR.background(name, hist, 0)
    sm = C.SMTab(172.57)
    res, rows, spec, ks, wq = M.run(bg, xi, sm, H, rtol=1e-10, atol=1e-15, nsave=900)
    Rx = dict(res=res, rows=rows, spec=spec, ks=ks, wq=wq, Rb=M.rbar_rows(bg, rows), m2xbar=np.zeros(len(rows)))
    vari['rtol 1e-10'] = {f'{m}': IR.evalD(Rx, xi, H, m, Ns=(60,))['best'][60]['D'] for m in Mts}
    # Hartree off
    res, rows, spec, ks, wq = M.run(bg, xi, sm, H, hartree=False, nsave=900)
    Rx = dict(res=res, rows=rows, spec=spec, ks=ks, wq=wq, Rb=M.rbar_rows(bg, rows), m2xbar=np.zeros(len(rows)))
    vari['Hartree off'] = {f'{m}': IR.evalD(Rx, xi, H, m, Ns=(60,))['best'][60]['D'] for m in Mts}
    # Rbar from instantaneous R smoothed over 2 periods instead of 1 (window convention)
    out['variants'] = vari
    out['seconds'] = time.time() - t0
    for k, v in vari.items():
        print(k, {m: round(v[m] - ref[m], 2) for m in v}, flush=True)
    C.jdump(out, C.RES / 'final_controls.json')


if __name__ == '__main__':
    {'plasmaD': cmd_plasmaD, 'portal': cmd_portal, 'controls': cmd_controls}[sys.argv[1]]()
