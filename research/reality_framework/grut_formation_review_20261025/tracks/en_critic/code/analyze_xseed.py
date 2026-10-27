"""Summarize X-seed-amplitude control runs against the published N48 curve and the published same-seed linear control.

usage: python3 -I analyze_xseed.py <packet_dir> <runs_dir> <out_json>
"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import trace_diag as D


def summarize(d, lin_t, lin_v):
    out = dict(energy_drift=float(np.max(abs(d['rho'] / d['rho'][0] - 1))),
               first_charge_contrast_one=D.first(d['tau'], d['contrast'], 1.),
               first_X_variance_half=D.first(d['tau'], d['fluct'], .5),
               first_phimix_1pct=D.first(d['tau'], d['phimix'], .01),
               first_phimix_10pct=D.first(d['tau'], d['phimix'], .1))
    i = int(np.argmax(d['phimix'])); out['max_phimix'] = [float(d['phimix'][i]), float(d['tau'][i])]
    out['at_800'] = {k: float(d[k][-1]) for k in ['phimix', 'Ephi', 'outside', 'Q0Q', 'fluct', 'contrast']}
    out['window_700_800'] = {k: D.window(d, k, 700, 800) for k in ['phimix', 'Ephi', 'outside', 'Q0Q']}
    rows = []
    for tau in [20, 60, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800]:
        i = int(np.argmin(abs(d['tau'] - tau))); j = int(np.argmin(abs(lin_t - tau)))
        rows.append(dict(tau=tau, phimix=float(d['phimix'][i]), linear_control=float(lin_v[j]),
                         contrast=float(d['contrast'][i]), outside=float(d['outside'][i])))
    out['samples'] = rows
    # tau where the nonlinear curve first falls 10% below the linear control (after tau 30)
    lin_on = np.interp(d['tau'], lin_t, lin_v); m = (d['tau'] > 30) & (d['phimix'] < .9 * lin_on)
    out['first_10pct_below_linear_control'] = float(d['tau'][m][0]) if m.any() else None
    return out


def main():
    pkt, runs, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    L = json.loads((pkt / 'LINEAR_PHI_CONTROL.json').read_text())['rows']
    lt = np.array([r['tau'] for r in L]); lv = np.array([r['phi_plus_interaction_raw_fraction'] for r in L])
    c = np.load(pkt / 'FIELD_CURVES.npz')
    pub = dict(tau=c['N48_tau'], eta=c['N48_eta'], phimix=c['N48_phi_plus_interaction_energy_fraction'],
               Ephi=c['N48_Ephi_fraction'], outside=c['N48_energy_outside_spatial_mean_fraction'],
               Q0Q=c['N48_spatial_mean_charge_fraction'], contrast=c['N48_charge_contrast_RMS'],
               fluct=c['N48_X_field_fluctuation_variance_fraction'], rho=c['N48_rho'])
    res = dict(method='Codex early3d engine (unmodified, compiled from git 7d9c93b) on own initial data: published '
                      'N48/L12/seed1 tau-20 field (f=1 file SHA-256 equals published), X nonzero Fourier modes '
                      '(chi1,chi2,pi1,pi2) times f, X zero mode and phi seed unchanged; threads 2, hmax .004, '
                      'sample .05 tau. Reproduction-type engine, new control.',
               published_N48_f1=summarize(pub, lt, lv), runs={})
    for info in sorted(runs.glob('*/info.json')):
        meta = json.loads(info.read_text())
        if meta['end_tau'] < 800: continue
        d = D.diag(D.load(info.parent / 'state.tsv'))
        res['runs'][f"f={meta['f']:g}"] = dict(f=meta['f'], seconds=meta['seconds'],
                                                initial_X_fluctuation_variance=meta['initial_X_fluctuation_variance'],
                                                **summarize(d, lt, lv))
    out.write_text(json.dumps(res, indent=2) + '\n')
    for k, v in [('published f=1', res['published_N48_f1'])] + list(res['runs'].items()):
        print(k, 'contrast1', v['first_charge_contrast_one'], 'half', v['first_X_variance_half'],
              'max', v['max_phimix'], '800', round(v['at_800']['phimix'], 5),
              'mean700-800', round(v['window_700_800']['phimix']['eta_mean'], 5),
              'below-linear', v['first_10pct_below_linear_control'], 'drift', v['energy_drift'])


if __name__ == '__main__':
    main()
