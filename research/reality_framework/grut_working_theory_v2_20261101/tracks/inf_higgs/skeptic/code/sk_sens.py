#!/usr/bin/env python3
"""Skeptic of inf_higgs: sensitivities the track did not vary.
  python3 -B sk_sens.py Hi <H_GeV>   # HT6t MD background in H_i units, SM at a different H_i (xi(H_i) re-run):
                                     # the slow-decay histories were never re-calibrated (track caveat)
  python3 -B sk_sens.py alphas       # alpha_s = 0.1180 -/+ 0.0009 on the HT6t MD and G05 runs -> M_t,crit
  python3 -B sk_sens.py tauth        # saturation of D(tau_th) by tau = 80 (MD) and the N-normalisation conventions
"""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np

import sk_common as C
import sk_bench as B

MTS = [172.28, 172.57, 172.70, 172.80, 172.86, 172.90, 172.95, 173.0, 173.05, 173.15]


def save(k, o):
    p = C.RES / 'sk_sens.json'
    try:
        d = json.load(open(p))
    except Exception:
        d = {}
    d[k] = o
    C.jdump(d, p)


def cmd_Hi(H):
    ev = B.run('HT6t', 'MD', 0, H_GeV=H)
    g = B.Dgrid(ev, Ns=(60.0,), H_GeV=H, Mts=MTS)
    Ds = [g[m][60.0]['D'] for m in MTS]
    o = dict(H_GeV=H, xi_H=ev['xi'], D60=dict(zip(map(str, MTS), Ds)), Mt_crit_60=C.crit(MTS, Ds))
    save(f'HT6t MD H_i={H:.3g}', o)
    print('H_i %.3g xi %.4f' % (H, ev['xi']), {m: round(d, 1) for m, d in zip(MTS, Ds)}, 'crit', o['Mt_crit_60'], flush=True)


def cmd_alphas():
    for hist in ('MD', 'G05'):
        ev = B.load(B.key('HT6t', hist, 0))
        ev['Rbar'] = C.period_mean(ev, 1)
        for als in (0.1171, 0.1180, 0.1189):
            Ds = [C.danger(ev, C.sm_in(m, 125.20, als), ev['Rbar'], Ns=(60.0,))[60.0]['D'] for m in MTS]
            o = dict(alpha_s=als, D60=dict(zip(map(str, MTS), Ds)), Mt_crit_60=C.crit(MTS, Ds))
            save(f'HT6t {hist} alpha_s={als}', o)
            print(hist, als, {m: round(d, 1) for m, d in zip(MTS, Ds)}, 'crit', o['Mt_crit_60'], flush=True)


def cmd_tauth():
    ev = B.load(B.key('HT6t', 'MD', 0))
    ev['Rbar'] = C.period_mean(ev, 1)
    o = {}
    for Mt in (172.86, 172.95, 173.15):
        smi = C.sm_in(Mt)
        row = {}
        for tt in (10, 20, 30, 40, 60, 80):
            row[str(tt)] = C.danger(ev, smi, ev['Rbar'], Ns=(60.0,), tau_th=tt)[60.0]['D']
        o[str(Mt)] = row
        print(Mt, {k: round(v, 1) for k, v in row.items()}, flush=True)
    # region-count conventions: (4 pi/3) and the comoving radius of the last-scattering sphere (~3.2/(a0 H0))
    o['conventions'] = dict(four_pi_over_3=math.log10(4 * math.pi / 3), r_LS_over_hubble_cubed=3 * math.log10(3.2),
                            total=math.log10(4 * math.pi / 3) + 3 * math.log10(3.2))
    save('HT6t MD D(tau_th)', o)


if __name__ == '__main__':
    c = sys.argv[1]
    if c == 'Hi':
        cmd_Hi(float(sys.argv[2]))
    elif c == 'alphas':
        cmd_alphas()
    elif c == 'tauth':
        cmd_tauth()
