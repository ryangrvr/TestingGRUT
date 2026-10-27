"""Skeptic (en_indep_skeptic): checks on Codex's PUBLISHED curves only (FIELD_CURVES.npz, LINEAR_PHI_CONTROL.json,
FIELD_RESULTS.json, RUNS.json from the 7d9c93b archive):
  (a) Codex N48 vs Codex halfstep48 diagnostic differences by tau segment (do step-size differences grow?);
  (b) control crossing times interpolated on the published dtau = 0.5 control curves vs the reported thresholds;
  (c) linear fixed-orbit phi control vs nonlinear N48 phi+interaction on the same tau grid;
  (d) exchange-work ledger residual by tau window from a Codex-engine N48 trace (en_rerun regeneration, read-only).
usage: published_checks.py PUB_PACKET_DIR EN_RERUN_N48_STATE_TSV OUT.json
"""
import sys, json
import numpy as np
pub, tsv, outp = sys.argv[1:4]
z = np.load(pub + '/FIELD_CURVES.npz'); fr = json.load(open(pub + '/FIELD_RESULTS.json'))['cases']
out = {}
# (a)
t48, th = z['N48_tau'], z['halfstep48_tau']; i48 = {round(t, 6): i for i, t in enumerate(t48)}
pairs = [(i48[round(t, 6)], j) for j, t in enumerate(th) if round(t, 6) in i48]
a = np.array([x[0] for x in pairs]); b = np.array([x[1] for x in pairs]); tt = th[b]
segs = [(20, 40), (40, 70), (70, 150), (150, 300), (300, 500), (500, 650), (650, 800.1)]
out['N48_vs_halfstep48_max_rel_by_segment'] = {
    k: {f'{s0}-{s1:g}': float((np.abs(z['N48_' + k][a] - z['halfstep48_' + k][b]) / np.abs(z['halfstep48_' + k][b]))[(tt >= s0) & (tt < s1)].max()) for s0, s1 in segs}
    for k in ['Ephi_fraction', 'phi_plus_interaction_energy_fraction', 'energy_outside_spatial_mean_fraction',
              'spatial_mean_charge_fraction', 'charge_contrast_RMS', 'X_field_fluctuation_variance_fraction']}
# (b)
cr = {}
for t in ['N32', 'N48', 'halfstep48', 'N64', 'seed2_48', 'cut10_48', 'box16_64']:
    tau = z[t + '_tau']; row = {}
    for key, lev, rep in (('charge_contrast_RMS', 1., 'first_charge_contrast_one'), ('X_field_fluctuation_variance_fraction', .5, 'first_X_variance_half')):
        v = z[t + '_' + key]; i = np.flatnonzero(v >= lev)[0]
        row[key] = dict(curve_dtau=float(tau[i] - tau[i - 1]), first_sample=float(tau[i]),
                        interpolated=float(tau[i - 1] + (lev - v[i - 1]) * (tau[i] - tau[i - 1]) / (v[i] - v[i - 1])),
                        reported_dense=fr[t]['thresholds'][rep])
    cr[t] = row
out['control_crossings'] = cr
# (c)
lc = json.load(open(pub + '/LINEAR_PHI_CONTROL.json'))
lt = np.array([r['tau'] for r in lc['rows']]); lv = np.array([r['phi_plus_interaction_raw_fraction'] for r in lc['rows']])
nt = z['N48_tau'][:-1]; nv = z['N48_phi_plus_interaction_energy_fraction'][:-1]
cmp = {}
for T in [20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300, 400, 600, 800]:
    i = int(np.argmin(abs(lt - T))); n = float(np.interp(lt[i], nt, nv)); cmp[str(T)] = dict(linear=float(lv[i]), nonlinear=n, ratio=float(lv[i] / n))
w = lt <= 100.01
out['linear_vs_nonlinear_phi'] = dict(rows=cmp, max_abs_log_ratio_tau_le_100=float(np.max(np.abs(np.log(lv[w] / np.interp(lt[w], nt, nv))))))
# (d)
t = np.genfromtxt(tsv, names=True, delimiter='\t')
rho = t['rho']; EX = rho - (t['Kphi'] + t['Gphi'] + t['Uphi'])
I = np.r_[0., np.cumsum(.5 * (t['work_phi_to_X'][1:] + t['work_phi_to_X'][:-1]) * np.diff(t['eta']))]
r = np.abs((EX - EX[0]) - I) / rho[0]; tau = t['tau']
I2 = np.r_[0., np.cumsum(.5 * (t['work_phi_to_X'][2::2] + t['work_phi_to_X'][:-2:2]) * np.diff(t['eta'][::2]))]
r2 = np.abs((EX[::2] - EX[0]) - I2) / rho[0]
out['work_ledger_en_rerun_N48'] = dict(overall=float(r.max()), tau_at_max=float(tau[np.argmax(r)]), tau_le_150=float(r[tau <= 150].max()),
                                       tau_ge_150=float(r[tau >= 150].max()), half_resolution_overall=float(r2.max()))
json.dump(out, open(outp, 'w'), indent=1)
print(json.dumps(out['work_ledger_en_rerun_N48']), out['linear_vs_nonlinear_phi']['max_abs_log_ratio_tau_le_100'])
for k, v in cr.items(): print(k, {kk: (round(vv['interpolated'], 3), round(vv['reported_dense'], 3)) for kk, vv in v.items()})
