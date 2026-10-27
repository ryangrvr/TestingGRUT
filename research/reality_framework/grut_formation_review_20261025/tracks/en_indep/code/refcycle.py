"""Time dependence of the reference projections and of the tau-800 diagnostics over the last orbit cycle (own fields)."""
import sys, os, json, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indep3d as I
from refproj import refproj
R = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/en_indep/runs'
out = {}
files = [f'{R}/bm6_h004/snap_tau{t:.2f}.npz' for t in (60, 70, 100, 150, 200, 400)] + sorted(glob.glob(f'{R}/lastcycle/snap_tau*.npz'))
for f in files:
    z = np.load(f); q, p = z['q'], z['p']; tau = float(z['tau'])
    r = refproj(q, p)
    mq = q.mean(axis=(1, 2, 3)); mp = p.mean(axis=(1, 2, 3))
    out[f'{tau:.2f}'] = dict(massless=r['massless']['nminus_over_net'], hartree=r['isotropic_Hartree']['nminus_over_net'],
                             hartree_mass2=r['isotropic_Hartree']['mass2'], mean_field_radius=float(np.hypot(mq[0], mq[1])))
    print(f'{tau:.2f}', json.dumps(out[f'{tau:.2f}']), flush=True)
# compare restart run with main run at tau 800
a = I.derived(np.load(f'{R}/bm6_h004/records.npy')); b = I.derived(np.load(f'{R}/lastcycle/records.npy'))
out['restart_vs_main_at_800_outside_mean_rel'] = float(abs(a['outside_mean'][-1] / b['outside_mean'][-1] - 1))
v = np.array([[out[k]['massless'], out[k]['hartree']] for k in out if k.startswith(('75', '76', '77', '78', '79', '80'))])
out['last_cycle_752_800'] = dict(massless_min=float(v[:, 0].min()), massless_max=float(v[:, 0].max()), massless_mean=float(v[:, 0].mean()),
                                 hartree_min=float(v[:, 1].min()), hartree_max=float(v[:, 1].max()), hartree_mean=float(v[:, 1].mean()),
                                 ratio_min=float((v[:, 0] / v[:, 1]).min()), ratio_max=float((v[:, 0] / v[:, 1]).max()))
print(json.dumps(out['last_cycle_752_800']), out['restart_vs_main_at_800_outside_mean_rel'])
json.dump(out, open(sys.argv[1], 'w'), indent=1)
