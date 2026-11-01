"""Aggregate the skeptic's runs into results JSONs and print the tables used in SKEPTIC.md:
sk_scan.json (fine eps scan), sk_linear.json (linear growth), sk_epst.json (tilt sensitivity), sk_slow.json
(eps <= 1e-6), sk_budget.json (energy/action budget of selected runs), sk_track_ratio.json (2:3 test on the track's
own per-run files)."""
import glob, json, os, numpy as np
from sk_common import *

def load_all():
    out = []
    for fn in glob.glob(os.path.join(RUNS, '*.json')):
        if fn.endswith('.ckpt.json'): continue
        d = json.load(open(fn)); out.append(d)
    return out

R = load_all()
def sel(**kw):
    rows = []
    for d in R:
        j = d['job']
        if all(abs(j[k] - v) < 1e-12 if isinstance(v, float) else j[k] == v for k, v in kw.items()):
            rows.append(d)
    return sorted(rows, key=lambda d: (d['job']['orbit'], -d['job']['eps'], d['job']['phase']))

def brief(d):
    j = d['job']; J0 = d['J0_per_q']
    return dict(orbit=j['orbit'], eps=j['eps'], phase=j['phase'], epst=j['epst'], growth=j['growth'], rm=j['rm'], mend=j['mend'],
                nminus=d['nminus'], Nphi=d['Nphi'], loss=1 - d['nminus'] / J0, ratio=(J0 - d['nminus']) / d['Nphi'] if d['Nphi'] > 1e-9 else None,
                combo_over_2J0=(2 * d['nminus'] + 3 * d['Nphi']) / (2 * J0), hext=d['max_abs_hext'], qrel=d['max_qrel'], steps=d['steps'],
                windows={k: (round(v['nminus'], 6), round(v['Nphi'], 6)) for k, v in d['windows'].items()}, runtime=d['runtime'],
                chi_end=d['trace'][-1]['chi_over_M'])

# ---- fine eps scan (exp growth, epst 0.002, mend 40)
scan = [brief(d) for d in sel(growth=0, epst=0.002, mend=40.0, mstart=0.5, nper=32, method=8)]
print('== fine eps scan (exp growth, epst 0.002, 40 Mapo): loss 1 - n-/J0 per phase')
for orb in ('tau20', 'late'):
    byeps = {}
    for r in scan:
        if r['orbit'] == orb: byeps.setdefault(r['eps'], []).append(r)
    for e in sorted(byeps, reverse=True):
        rs = sorted(byeps[e], key=lambda r: r['phase'])
        print(f"  {orb:5s} eps {e:.1e}: " + '  '.join(f"ph{r['phase']:.2f} {r['loss']:+.4f}" for r in rs) +
              f"   | captured (loss>5%): {sum(r['loss'] > 0.05 for r in rs)}/{len(rs)}; ratio of lossy runs: " +
              ', '.join(f"{r['ratio']:.3f}" for r in rs if r['loss'] > 0.05))
save_json(os.path.join(RESULTS, 'sk_scan.json'), scan)

# ---- linear growth
lin = [brief(d) for d in sel(growth=1)]
print('== linear growth M = M0 (1 + kappa s), kappa matched to eps*omega_r at rm*Mperi, 10 Mapo')
for r in lin:
    print(f"  {r['orbit']:5s} eps_eff {r['eps']:.0e} rm {r['rm']:g} ph {r['phase']:.1f}: n- {r['nminus']:.4f} Nphi {r['Nphi']:.4f} loss {r['loss']:+.4f}")
save_json(os.path.join(RESULTS, 'sk_linear.json'), lin)

# ---- epst
ept = [brief(d) for d in R if abs(d['job']['epst'] - 0.002) > 1e-9]
ept.sort(key=lambda r: (r['orbit'], r['epst'], r['phase']))
print('== tilt sensitivity (exp growth; tau20 at eps 1e-5, late at eps 1e-4)')
for r in ept:
    print(f"  {r['orbit']:5s} epst {r['epst']:.4f} eps {r['eps']:.0e} ph {r['phase']:.3f}: n- {r['nminus']:.4f} Nphi {r['Nphi']:.4f} loss {r['loss']:+.4f} ratio {r['ratio'] if r['ratio'] is None else round(r['ratio'],3)} chi_end {r['chi_end']:+.4f}")
save_json(os.path.join(RESULTS, 'sk_epst.json'), ept)

# ---- slow runs
slow = [brief(d) for d in R if d['job']['eps'] <= 1.0001e-6 or d['job']['mend'] < 40 and d['job']['growth'] == 0]
slow.sort(key=lambda r: (r['orbit'], -r['eps'], r['phase']))
print('== slow runs (eps <= 1e-6 or reduced M_end)')
for r in slow:
    print(f"  {r['orbit']:5s} eps {r['eps']:.0e} ph {r['phase']:.6f} mend {r['mend']:g}: n- {r['nminus']:.5f} Nphi {r['Nphi']:.5f} loss {r['loss']:.4f} ratio {r['ratio']:.4f} windows {r['windows']} hext {r['hext']:.1e} q {r['qrel']:.1e} t {r['runtime']:.0f}s")
save_json(os.path.join(RESULTS, 'sk_slow.json'), slow)

# ---- budget of selected runs (final window)
print('== budget at the final window (per q): E_X/(m_X q) = 1 + 2 n-, E_phi/(m_phi q) = N_phi, E_int/(m_X q), quartic/E_X, work')
bud = []
for d in sel(growth=0, epst=0.002, mstart=0.5, nper=32, method=8):
    j = d['job']
    if j['phase'] not in (0.0,) or j['eps'] not in (1e-3, 1e-4, 1e-5, 1e-6): continue
    k = f"{j['mend']:g}"; w = d['windows'][k]; J0 = d['J0_per_q']
    row = dict(orbit=j['orbit'], eps=j['eps'], mend=j['mend'], EX_over_mXq=w['EX_over_mXq'], one_plus_2n=1 + 2 * w['nminus'],
               Ephi_lin=w['Ephi_lin_over_mphiq'], Ephi_exact=w['Ephi_exact_over_mphiq'], Eint=w['Eint_over_mXq'], quart=w['Equart_over_EX'],
               dchi_rms=w['dchi_rms_over_M'], nminus=w['nminus'], nminus_actual=w['nminus_actual'], Nphi=w['Nphi'], J0=J0,
               action_in=J0, action_out_X=w['nminus'], action_out_phi=w['Nphi'], deficit_X=J0 - w['nminus'], ratio=(J0 - w['nminus']) / w['Nphi'],
               W_over_E=w['W'] / w['E'], hext=w['hext_rel'], qrel=w['qrel'])
    bud.append(row)
    print(f"  {row['orbit']:5s} eps {row['eps']:.0e} @ {row['mend']:g} Mapo: EX/(mXq) {row['EX_over_mXq']:.6f} vs 1+2n- {row['one_plus_2n']:.6f}; Ephi_lin {row['Ephi_lin']:.5f} exact {row['Ephi_exact']:.5f}; Eint {row['Eint']:+.1e}; quartic/EX {row['quart']:.1e}; n-(gM) {row['nminus']:.5f} n-(g|chi|) {row['nminus_actual']:.5f}; dJ_r {row['deficit_X']:.4f} = {row['ratio']:.3f} x Nphi; W/E {row['W_over_E']:.4f}; Hext {row['hext']:.0e} q {row['qrel']:.0e}")
save_json(os.path.join(RESULTS, 'sk_budget.json'), bud)

# ---- 2:3 test on the track's own per-run files
J0 = {'tau20': 1.7637252793969966, 'late': 0.48049184685032387}
tr = {}
for orb in ('tau20', 'late'):
    for eps in ('1e-02', '1e-03', '3e-04', '1e-04', '3e-05', '1e-05', '1e-06'):
        fs = glob.glob(os.path.join(TRACK_RUNS, f'{orb}_m3_et0.002_b0_phi0qs_eps{eps}_ph*_rt1e-11.json'))
        if not fs: continue
        rat, comb = [], []
        for f in fs:
            d = json.load(open(f)); dJ = J0[orb] - d['nminus']
            comb.append((2 * d['nminus'] + 3 * d['Nphi']) / (2 * J0[orb]))
            if dJ / J0[orb] > 0.05: rat.append(dJ / d['Nphi'])
        tr[f'{orb}_eps{eps}'] = dict(n=len(fs), n_lossy=len(rat), ratio_mean=float(np.mean(rat)) if rat else None, ratio_sd=float(np.std(rat)) if rat else None,
                                  ratio_min=float(np.min(rat)) if rat else None, ratio_max=float(np.max(rat)) if rat else None,
                                  combo_over_2J0_mean=float(np.mean(comb)), combo_over_2J0_sd=float(np.std(comb)))
print('== track per-run files: (J0 - n-)/Nphi over lossy runs (loss > 5%), and (2n- + 3Nphi)/(2J0)')
for k, v in tr.items():
    print(f"  {k}: n {v['n']} lossy {v['n_lossy']} ratio {v['ratio_mean']} sd {v['ratio_sd']} [{v['ratio_min']}, {v['ratio_max']}]  combo/2J0 {v['combo_over_2J0_mean']:.5f} +- {v['combo_over_2J0_sd']:.5f}")
save_json(os.path.join(RESULTS, 'sk_track_ratio.json'), tr)
