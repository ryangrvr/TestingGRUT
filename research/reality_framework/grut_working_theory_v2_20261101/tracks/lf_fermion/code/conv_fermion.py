"""Convergence of the line fermion runs: fermion grid (N_f 8192 -> 16384), time step (cfl 0.5 -> 0.25), realisations
(n_r 4 -> 8, new seed).  Claude, lf_fermion (issue #13 context). -> ../results/conv_fermion.json"""
import json, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
import lfcommon as LC
v = 367046.6168416639
D = LC.SP / 'lfq' / 'REC_P0'


def load(tag, line):
    d = json.loads((D / ('%s_l%d_chi.json' % (tag, line))).read_text())
    t = np.array([r['t'] for r in d['meas']])
    out = {}
    for nm in ('N2', 'N3', 'N4'):
        e = np.array([r[nm]['e_proj_cut'] for r in d['meas']])
        out[nm] = (t, (e - e[t < 0.15].mean()) / v ** 2, e[t < 0.15].mean() / v ** 2)
    return out


res = {}
for name, (tag, line) in {'nf16k': ('nf16k', 1), 'cfl25': ('cfl25', 1), 'nr8': ('nr8', 0)}.items():
    a, b = load('base', line), load(tag, line)
    rows = {}
    for nm in ('N2', 'N3', 'N4'):
        ts = [0.3, 0.5, 0.6, 0.7, 0.8]
        ea = np.interp(ts, a[nm][0], a[nm][1])
        eb = np.interp(ts, b[nm][0], b[nm][1])
        rows[nm] = dict(t=ts, base=list(ea), test=list(eb), rel=list(eb / ea - 1), baseline_offset=[a[nm][2], b[nm][2]])
        print(name, 'line', line, nm, ' '.join('t%.1f %+.3f' % (t, r) for t, r in zip(ts, eb / ea - 1)))
    res[name] = dict(line=line, species=rows)
Path(__file__).resolve().parents[1].joinpath('results', 'conv_fermion.json').write_text(json.dumps(res, indent=1))
