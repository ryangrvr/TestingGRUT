#!/usr/bin/env python3
"""Track mf_higgs: alpha_s = 0.1180 -+ 0.0009 (N5) on the pieces not covered by mh_run bench / mh_inflation sens:
D with the N -> L H plasma on EXP and PAULI0, and the instantaneous-HM inflationary bracket cap.  -> results/alphas.json"""
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
import mh_common as C
import mh_run as RN
import mh_inflation as INF

out = {}
for hist in ('EXP', 'PAULI0'):
    bg = RN.bgc('HFa', hist)
    R = RN.run_cached('HFa', hist, RN.xi_at(bg.Hi_GeV))
    N = R['res']['N_obs']
    o = {}
    for als in (0.1171, 0.1189):
        rows = {}
        for Mt in (172.28, 172.57, 172.86, 173.45, 174.05, 174.6):
            rows[f'{Mt}'] = {'grav': RN.D_of(R, Mt, N, als=als)['D'], 'plasma': RN.D_of(R, Mt, N, als=als, plasma=1.0)['D'],
                             'plasma_c1e-4': RN.D_of(R, Mt, N, als=als, plasma=1e-4)['D']}
        o[f'{als}'] = rows
        print(hist, als, {k: (round(v['grav'], 1), v['plasma'], v['plasma_c1e-4']) for k, v in rows.items()}, flush=True)
    out[hist] = o
d = C.load_bench('HFa')
W = INF.waterfall_table('HFa')
N0 = C.N_obs('HFa', d['reheating']['K']) + W['N_end']
caps = {}
for als in (0.1171, 0.1180, 0.1189):
    caps[f'{als}'] = INF.Mt_cap(lambda m: INF.log10N_inst(INF.SMfull(round(m, 3), als=als), 'HFa', d['H_i_GeV'], 0.1, N0, False)[0],
                                lo=172.0, hi=175.5)
    print('inst cap', als, caps[f'{als}'], flush=True)
out['inst_HM_cap_tree'] = caps
C.jdump(out, C.RES / 'alphas.json')
