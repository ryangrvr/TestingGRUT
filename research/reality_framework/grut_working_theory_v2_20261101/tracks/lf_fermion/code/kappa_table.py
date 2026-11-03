"""Calibrated closure table kappa_j(t) for the drain (model choice), from calib_seq_<RUN>.json:
   t <= 0.45: kappa_hom,j = geometric mean of the sequential fits over windows ending 0.25-0.45 (pre-fragmentation);
   t >= 0.55: kappa_frag,j = geometric mean over windows ending 0.55-1.10 (fragmented stage), held constant after 1.1
   (extrapolation [A]); linear in between.  Variants: 'mid' (mean fit), 'lo' (fits to Q - sem), 'hi' (Q + sem, capped
   at 1e3).  Claude, lf_fermion (issue #13 context). Usage: python3 -B kappa_table.py RUN -> ../results/kappa_<RUN>.json"""
import json, math, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
import lfcommon as LC

SPI = {'N1': 0, 'N2': 1, 'N3': 2, 'N4': 3}


def table(run, t_sw=(0.45, 0.55)):
    d = json.loads((LC.TRACK / 'results' / ('calib_seq_%s.json' % run)).read_text())
    tw = np.array(d['t_window_end'])
    out = {}
    for var, key in (('mid', 'kappa_inhom'), ('lo', 'kappa_inhom_minus_sem'), ('hi', 'kappa_inhom_plus_sem')):
        row = {}
        for nm, r in d['species'].items():
            k = np.minimum(np.array(r[key]), 1e3)
            early = (tw >= 0.249) & (tw <= t_sw[0] + 1e-3)
            late = (tw >= t_sw[1] - 1e-3) & (tw <= 1.101)
            ke = float(np.exp(np.mean(np.log(np.maximum(k[early], 1.0)))))
            kl = float(np.exp(np.mean(np.log(np.maximum(k[late], 1.0)))))
            row[nm] = dict(k_hom=ke, k_frag=kl, t_sw=list(t_sw))
        out[var] = row
    return out


def par_lines(tab_var):
    lines = []
    for nm, r in tab_var.items():
        j = SPI[nm]
        t0, t1 = r.get('t_sw', (0.45, 0.55))
        lines += ['kappa %d %.4g %.6g' % (j, t0, r['k_hom']), 'kappa %d %.4g %.6g' % (j, t1, r['k_frag'])]
    return lines


if __name__ == '__main__':
    run = sys.argv[1] if len(sys.argv) > 1 else 'REC_P0'
    sw = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (0.45, 0.55)
    T = table(run, sw)
    for v, r in T.items():
        print(v, {nm: (round(x['k_hom'], 2), round(x['k_frag'], 1)) for nm, x in r.items()})
    LC.dump(LC.TRACK / 'results' / ('kappa_%s.json' % run), dict(run=run, label='C/A', note=__doc__, t_switch=list(sw), table=T))
