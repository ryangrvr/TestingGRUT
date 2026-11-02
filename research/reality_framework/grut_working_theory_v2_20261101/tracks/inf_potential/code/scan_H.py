#!/usr/bin/env python3
"""inf_potential (Claude, own code): reach of the families in the slice Hubble rate H_i = 1e11 - 1e12 GeV (history
Gamma_I = 0.5 H_i).  For each H_i: the calibrated shape parameter, N_*, n_s, r, m_I/H_i, and which Y windows hold
(Y resonance needs m_I < 5 or > 30 H_i [Q brief]; the fundamental inflaton -> Y_L Y_L channel opens at m_I > 44.50 H_i [D]).
Output ../results/scan_H.json.  Labels [C]."""
import json, math, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inf_common as ic, calibrate as cal
RES = Path(__file__).resolve().parent.parent / 'results'
out = dict(method=__doc__, rows=[])
t0 = time.time()
for Hi in (1.0e11, 1.29e11, 1.5e11, 2.0e11, 3.0e11, 5.0e11, 1.0e12):
    for name in ('T1', 'HT6t'):
        if name == 'T1':
            m, r = cal.calibrate_alpha('T', Hi, ('gamma', 0.5)); par = dict(alpha=m.p['alpha'])
        else:
            m, r = cal.calibrate_hilltop_ns(Hi, ('gamma', 0.5), 6, x0=(0.3469 * (Hi / 1.29e11) ** (2 / 3), 1.553e-5 * (Hi / 1.29e11) ** 0.5))
            par = dict(vev=m.vev, beta=m.beta)
        hs = r['h_slice']; mI = math.sqrt(m.m2_min()) / hs
        row = dict(H_i_GeV=Hi, benchmark=name, params=par, Nstar=r['Nstar'], ns=r['ns'], r=r['r'], alpha_s=r['alpha_s'],
                   m_I_over_Hi=mI, m_I_GeV=mI * Hi, Y_resonance_avoided=(mI < 5 or mI > 30), YY_fundamental_closed=(mI < 44.50))
        out['rows'].append(row)
        print('%-5s H_i=%.2e %s N*=%.2f ns=%.4f r=%.2e mI/Hi=%.1f (mI=%.2e GeV) Yres-avoided=%s YY-closed=%s [%.0fs]' % (
            name, Hi, {k: '%.4g' % v for k, v in par.items()}, r['Nstar'], r['ns'], r['r'], mI, mI * Hi, row['Y_resonance_avoided'],
            row['YY_fundamental_closed'], time.time() - t0), flush=True)
        json.dump(out, open(RES / 'scan_H.json', 'w'), indent=1, default=float)
