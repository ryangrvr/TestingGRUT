"""Orbit-cycle means of the late diagnostics (own N48 run): windows of one zero-mode period in conformal time."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indep3d as I
R = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/en_indep/runs'
D = I.derived(np.load(f'{R}/bm6_h004/records.npy'))
s, tau = D['s'], D['tau']
# period from successive maxima of outside_mean for tau>=400 (records every 0.25 tau there)
sel = np.flatnonzero(tau >= 300); v = D['outside_mean'][sel]
pk = sel[[i for i in range(1, len(v) - 1) if v[i] > v[i - 1] and v[i] >= v[i + 1]]]
P = float(np.mean(np.diff(s[pk])))
out = dict(period_s=P, peak_taus=[float(x) for x in tau[pk]], windows=[])
for k in range(len(pk) - 1):
    w = (s >= s[pk[k]]) & (s < s[pk[k + 1]])
    # time-weighted mean in s
    row = dict(tau_start=float(tau[pk[k]]), tau_end=float(tau[pk[k + 1]]))
    for key in ('outside_mean', 'phi_plus_mix_fraction', 'Ephi_fraction', 'mean_charge_fraction', 'charge_contrast', 'Xvar_fraction'):
        ss, vv = s[w], D[key][w]
        row[key] = dict(mean=float(np.trapezoid(vv, ss) / (ss[-1] - ss[0])), min=float(vv.min()), max=float(vv.max()))
    out['windows'].append(row)
    print(json.dumps({k: (v if not isinstance(v, dict) else round(v['mean'], 5)) for k, v in row.items()}))
print('period_s', P)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
