#!/usr/bin/env python3
"""dq_heat (Claude, own code): S/Q and e_crit if the GMR modes left the horizon with the inflaton at s = phi_I/phi_I(0) > 1.

Only the superhorizon part (k <= 1) of the vector number is rescaled by (F_inf(1)/F_inf(s))^2 (finf_family.json); the
release-era part (k >= 1), the mediators and the cascade are unchanged.  S/Q from the copied Hartree solver (heat_table).
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
from heat_table import SQ_of_W, W_THR, SIG

RES = Path(__file__).resolve().parent.parent / 'results'
ht = json.load(open(RES / 'heat_table.json'))['budgets']['portals_vanish_at_MY_benchmark']
fam = json.load(open(RES / 'finf_family.json'))['rows']
out = []
for r in fam:
    f = r['N_A_factor_vs_s1']
    rows = []
    for x in ht['rows']:
        Rv = f * x['R_vec_k_le_1'] + x['R_vec_k_ge_1_release_era'] + x['R_vec_pieceB_SU']
        Rt = Rv + x['R_Y_hom'] + x['R_Y_k'] + x['R_cascade']
        rows.append((x['e'], Rt, ht['W0'] * (1 + Rt)))
    le = np.log([a for a, _, _ in rows]); lW = np.log([c for _, _, c in rows])
    g = lambda y: np.interp(y, le, lW) - math.log(W_THR['bulk_Q_inf'])
    ec = math.exp(brentq(g, le[0], le[-1])) if g(le[0]) > 0 > g(le[-1]) else (None if g(le[0]) <= 0 else float('inf'))
    r6 = [t for t in rows if t[0] == 1e-6][0]
    s6 = SQ_of_W(r6[2])['S_per_Q']
    out.append(dict(s=r['s'], F_inf=r['F_inf'], N_A_factor=f, R_total_e1e6=r6[1], W_e1e6=r6[2], S_per_Q_e1e6=s6,
                    margin_bulk_e1e6=s6 / SIG['bulk_Q_inf'], e_crit_bulk=ec))
    print('s=%.2f F=%.1f factor=%.3f R_tot(1e-6)=%.3f S/Q(1e-6)=%.3f (x%.3f) e_crit(bulk)=%s' % (r['s'], r['F_inf'], f, r6[1], s6, s6 / SIG['bulk_Q_inf'], ec))
# the s at which e = 1e-6 sits exactly at the bulk threshold (interpolate the factor)
fs = [o['N_A_factor'] for o in out]; m = [o['margin_bulk_e1e6'] for o in out]
s_star = float(np.interp(1.0, m[::-1], [o['s'] for o in out][::-1]))
json.dump(dict(method=__doc__, rows=out, s_at_threshold_e1e6=s_star), open(RES / 'finf_sq.json', 'w'), indent=1, default=float)
print('s at which e = 1e-6 reaches S/Q = 3.4716:', s_star)
