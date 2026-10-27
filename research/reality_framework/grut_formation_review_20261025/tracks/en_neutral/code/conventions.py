"""en_neutral item 1: convention bookkeeping from Codex's PUBLISHED curves (git 7d9c93b FIELD_CURVES.npz,
LINEAR_PHI_CONTROL.json), read as data.  Own arithmetic only.
usage: python3 -I conventions.py <published_dir> <out_json>
"""
import sys, json
from pathlib import Path
import numpy as np

pub, out = Path(sys.argv[1]), Path(sys.argv[2])
d = np.load(pub / 'FIELD_CURVES.npz')
lin = json.load(open(pub / 'LINEAR_PHI_CONTROL.json'))
lt = np.array([r['tau'] for r in lin['rows']]); lv = np.array([r['phi_plus_interaction_raw_fraction'] for r in lin['rows']])
res = dict(label='own arithmetic on Codex published curves (7d9c93b)', controls={})
for tag in ('N32', 'N48', 'halfstep48', 'N64', 'seed2_48', 'cut10_48', 'box16_64'):
    t = d[tag + '_tau']; f = d[tag + '_phi_plus_interaction_energy_fraction']; e = d[tag + '_Ephi_fraction']
    res['controls'][tag] = dict(phi_plus_int_tau20=float(f[0]), Ephi_tau20=float(e[0]), phi_plus_int_mean_tau20_40=float(f[t <= 40].mean()),
                                phi_plus_int_tau800=float(f[-1]), Ephi_tau800=float(e[-1]), phi_plus_int_mean_tau700_800=float(f[t >= 700].mean()),
                                seed_share_of_final=float(f[0] / f[-1]), net_increase_over_seed=float(f[-1] - f[0]),
                                final_over_seed=float(f[-1] / f[0]))
n48 = d['N48_tau']; f48 = d['N48_phi_plus_interaction_energy_fraction']
res['N48_nonlinear_vs_linear_same_grid'] = {str(x): dict(nonlinear=float(f48[np.argmin(abs(n48 - x))]), linear=float(lv[np.argmin(abs(lt - x))]),
                                                         ratio_lin_over_nl=float(lv[np.argmin(abs(lt - x))] / f48[np.argmin(abs(n48 - x))]))
                                            for x in (20, 40, 60, 70, 80, 100, 134, 200, 300, 400, 600, 800)}
res['linear_tau20_vs_nonlinear_tau20_abs_diff'] = float(abs(lv[0] - f48[0]))
res['range_final_phi_plus_int'] = [min(c['phi_plus_int_tau800'] for c in res['controls'].values()), max(c['phi_plus_int_tau800'] for c in res['controls'].values())]
res['range_final_Ephi'] = [min(c['Ephi_tau800'] for c in res['controls'].values()), max(c['Ephi_tau800'] for c in res['controls'].values())]
res['range_net_increase'] = [min(c['net_increase_over_seed'] for c in res['controls'].values()), max(c['net_increase_over_seed'] for c in res['controls'].values())]
out.write_text(json.dumps(res, indent=1) + '\n')
for k, v in res['controls'].items(): print(k, {a: round(b, 6) for a, b in v.items()})
print(res['range_final_phi_plus_int'], res['range_final_Ephi'], res['range_net_increase'], res['linear_tau20_vs_nonlinear_tau20_abs_diff'])
