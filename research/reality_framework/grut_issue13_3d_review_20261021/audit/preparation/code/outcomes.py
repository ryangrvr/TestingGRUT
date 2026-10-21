"""Post-formation bias growth and outcome classes for Codex's 12 runs plus this audit's reruns.

Usage: python3 -I outcomes.py <codex_raw_dir> <scratch_runs_dir> <results_dir>
Reads csv/topology series (Codex's: read only). Writes results/outcomes.json.
Bias law (Claude bias_threshold): m(K) = erf(y_f K^gamma), gamma = 1.41, x10 = 1.54 (network < 10%).
K = comoving V/A growth = A_basin(max after formation) / A_basin(a=2.2), basin = Codex's basin interface estimator.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
from scipy.special import erf, erfinv

RAW, RUNS, RES = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
GAMMA, X10 = 1.41, 1.54
CODEX = ['main_s101', 'main_s202', 'main_s303', 'cold_s101', 'kappa4_s101', 'dt_half_s101', 'grid256_s101',
         'grid384_s101', 'large_s101', 'gaussian_s101', 'gaussian_s202', 'gaussian_large_s101']


def y_of(m):
    return float(erfinv(np.clip(m, -1 + 1e-12, 1 - 1e-12)))


def one(name, base):
    d = np.genfromtxt(base / f'{name}.csv', delimiter=',', names=True)
    top = [json.loads(l) for l in open(base / f'{name}.topology.jsonl')]
    fp, fm = d['fplus_05'], d['fminus_05']
    form = np.flatnonzero((fp >= 0.1) & (fm >= 0.1))
    sat = np.flatnonzero(fp + fm >= 0.5)
    i0 = int(form[0]) if len(form) else (int(sat[0]) if len(sat) else None)
    if i0 == 0:  # preparation already in one broken vacuum (cold control): no tachyonic exit
        i0 = None
    rec = {'name': name, 'formed_10_10': bool(len(form)), 'a_form_sample': float(d['a'][form[0]]) if len(form) else None,
           'end': {'fplus': float(fp[-1]), 'fminus': float(fm[-1]), 'm': float((fp[-1] - fm[-1]) / max(fp[-1] + fm[-1], 1e-300))}}
    if i0 is None:
        rec['class_observed'] = 'no tachyonic exit (single preparation vacuum, wall-free)'
        rec['end_topology'] = {'minus_occupied': [json.loads(l) for l in open(base / f'{name}.topology.jsonl')][-1]['minus']['occupied_fraction']}
        return rec
    m0 = float((fp[i0] - fm[i0]) / (fp[i0] + fm[i0]))
    rec['first_saturated_or_formed'] = {'a': float(d['a'][i0]), 'm': m0, 'y': y_of(m0)}
    post = [t for t in top if t['a'] >= d['a'][i0] - 1e-9 and t['basin_link_area_comoving_over_V'] > 0]
    end_t = top[-1]
    if post:
        jmax = int(np.argmax([t['basin_link_area_comoving_over_V'] for t in post]))
        tmax = post[jmax]
        jj = int(np.argmin(np.abs(d['a'] - tmax['a'])))
        mmax = float((fp[jj] - fm[jj]) / (fp[jj] + fm[jj]))
        Amax, Aend = tmax['basin_link_area_comoving_over_V'], end_t['basin_link_area_comoving_over_V']
        K = Amax / Aend if Aend > 0 else float('inf')
        yf, ye = y_of(mmax), y_of(rec['end']['m'])
        Kd = 1.48  # dense unbiased reference (Codex gaussian_s101 basin K, a=1.6->2.2)
        rec['bias_law_dense_Kref'] = {'K_ref': Kd, 'takeover_pred': bool(abs(y_of(m0)) * Kd ** GAMMA >= X10),
                                      'm_f_used': m0, 'threshold_m_f': float(erf(X10 * Kd ** -GAMMA))}
        rec['bias_law_sparse_Kref'] = {'K_ref': 1.0, 'takeover_pred': bool(abs(y_of(m0)) >= X10), 'threshold_m_f': float(erf(X10))}
        rec['bias_law'] = {'a_at_max_basin_area': tmax['a'], 'A_max': Amax, 'A_end': Aend, 'K_basin': K, 'm_f': mmax,
                           'y_f': yf, 'y_end': ye,
                           'gamma_eff': float(np.log(ye / yf) / np.log(K)) if (np.isfinite(K) and K > 1.02 and yf * ye > 0) else None,
                           'K10_needed': float((X10 / abs(yf)) ** (1 / GAMMA)) if yf != 0 else float('inf'),
                           'm_end_pred_gamma1.41': float(erf(yf * K ** GAMMA)) if np.isfinite(K) else 1.0 * np.sign(yf),
                           'takeover_pred_gamma1.41': bool(np.isfinite(K) and abs(yf) * K ** GAMMA >= X10) or not np.isfinite(K)}
    gone = (min(fp[-1], fm[-1]) < 0.01) or (end_t['basin_link_area_comoving_over_V'] < 0.1 * max(
        [t['basin_link_area_comoving_over_V'] for t in top]))
    rec['end_topology'] = {'plus_winding': end_t['plus']['any_winding'], 'minus_winding': end_t['minus']['any_winding'],
                           'basin_area': end_t['basin_link_area_comoving_over_V']}
    if gone:
        rec['class_observed'] = 'one-sign (minority removed by a=2.2)' + ('' if len(form) else '; 10%/10% formation never reached')
    elif abs(rec['end']['m']) >= 0.1:
        rec['class_observed'] = 'biased two-sign network'
    else:
        rec['class_observed'] = 'near-unbiased two-sign network'
    return rec


def main():
    out = {'codex_runs': {}, 'audit_reruns': {}}
    for n in CODEX:
        out['codex_runs'][n] = one(n, RAW)
    for p in sorted(RUNS.glob('int_*.csv')) + sorted((RUNS / 'fb').glob('fb_*.csv')) + sorted(RUNS.glob('repro_*.csv')):
        n = p.stem
        if (p.parent / f'{n}.topology.jsonl').exists() and (p.parent / f'{n}.meta.json').exists():
            out['audit_reruns'][n] = one(n, p.parent)
    for grp in out.values():
        for n, r in grp.items():
            b = r.get('bias_law', {})
            print('%-22s %-58s m_end=% .3f K=%s m_f=%s gam_eff=%s K10=%s takeover_pred=%s' % (
                n, r['class_observed'], r['end']['m'], ('%.2f' % b['K_basin']) if b else '-', ('%.3f' % b['m_f']) if b else '-',
                ('%.2f' % b['gamma_eff']) if b.get('gamma_eff') else '-', ('%.2f' % b['K10_needed']) if b else '-',
                b.get('takeover_pred_gamma1.41', '-')))
    (RES / 'outcomes.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
