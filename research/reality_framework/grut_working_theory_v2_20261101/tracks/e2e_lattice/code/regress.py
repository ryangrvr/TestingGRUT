"""Regression R1: the lattice code with no fluctuations (all sites identical) against mf_chain's homogeneous stage 1
(hist_build.stage1, read-only) for the three closures, and the time-step convergence of the homogeneous solution.
Usage: python3 -B regress.py -> ../results/regress_homog.json
"""
import json, math, sys
import numpy as np
import latcommon as LC

out = {}
d, P = LC.model('HFa')
v = P['v']
for clos in ('track', 'pauli', 'hubble'):
    fr = {'track': 'track', 'pauli': 'pauli', 'hubble': 'hubble'}[clos]
    ref = LC.HB.stage1(P, fric=fr, Gth=0.0, pert_free=True, t_sw=0.6, dt_rec=2e-4)
    row = dict()
    for dtv in (2e-4, 1e-4, 5e-5):
        tag = 'reg_%s_%g' % (clos, dtv)
        LC.make_run(tag, 'HFa', closure=clos, N=4, L=0.05, dt=dtv, fluct=False, t0=-0.02, t_end=0.6,
                    ndiag=max(1, int(round(2e-4 / dtv))), snaps=[])
        LC.run(tag, budget=300, threads=1)
        D = LC.load_diag(tag)
        t = D['t']
        rec = {}
        for tt in (0.1, 0.19, 0.3, 0.45, 0.6):
            rec['t%.2f' % tt] = dict(chi_v=float(np.interp(tt, t, D['chi_mean'])), H=float(np.interp(tt, t, D['H'])),
                                     frad=float(np.interp(tt, t, D['fN1'] + D['fN2'] + D['fN3'] + D['fN4'] + D['fSM'])),
                                     fchi=float(np.interp(tt, t, D['fchi'])))
        j = int(np.argmax(D['eps'] >= 1.0))
        rec['t_eps1'] = float(t[j]); rec['H_eps1'] = float(D['H'][j]); rec['R_min'] = float(D['R_mean'][t > 0].min())
        x = D['chi_mean']
        tp = [(float(t[k]), float(x[k])) for k in range(1, len(t) - 1) if t[k] > 0.2 and x[k] < x[k - 1] and x[k] <= x[k + 1]]
        rec['inner_turning_points'] = tp[:6]
        rec['budget_max'] = float(np.abs(D['budget']).max())
        row['lat_dt%g' % dtv] = rec
    tr = ref['t']
    rr = {}
    frad = ref['fr']
    for tt in (0.1, 0.19, 0.3, 0.45, 0.6):
        rr['t%.2f' % tt] = dict(chi_v=float(np.interp(tt, tr, ref['x'] / v)), H=float(np.interp(tt, tr, ref['H'])),
                                frad=float(np.interp(tt, tr, frad)), fchi=float(np.interp(tt, tr, ref['fc'])))
    j = int(np.argmax(ref['eps'] >= 1.0))
    rr['t_eps1'] = float(tr[j]); rr['H_eps1'] = float(ref['H'][j]); rr['R_min'] = float(ref['R'].min())
    x = ref['x'] / v
    rr['inner_turning_points'] = [(float(tr[k]), float(x[k])) for k in range(1, len(tr) - 1)
                                  if tr[k] > 0.2 and x[k] < x[k - 1] and x[k] <= x[k + 1]][:6]
    row['mf_chain_stage1'] = rr
    out[clos] = row
    print(clos, json.dumps({k: (v_['t0.30'], v_['t_eps1'], v_['inner_turning_points'][:3]) for k, v_ in row.items()}), flush=True)
LC.dump(LC.RES / 'regress_homog.json', dict(label='C: lattice code (no fluctuations) vs mf_chain hist_build.stage1', cases=out))
