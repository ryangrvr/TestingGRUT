"""Paired forced-bias analysis: biased run vs its unbiased partner with identical fluctuations.

Usage: python3 -I paired_boxes.py <codex_raw_dir> <fb_dir> <out_json>
Pairs: L=51.2 fb_mu{mu} vs fb_mu0 (this audit); L=76.8 fb_N192_mu0.02 vs Codex gaussian_s101;
L=115.2 fb_N288_mu0.02 vs Codex gaussian_large_s101. Response y_b - y_0 isolates the systematic part.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
from scipy.special import erf, erfinv

RAW, FB, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])


def series(base, name):
    d = np.genfromtxt(base / f'{name}.csv', delimiter=',', names=True)
    top = [json.loads(l) for l in open(base / f'{name}.topology.jsonl')]
    return d, top


def m_at(d, a):
    i = int(np.argmin(np.abs(d['a'] - a)))
    fp, fm = d['fplus_05'][i], d['fminus_05'][i]
    return float((fp - fm) / max(fp + fm, 1e-300))


def basin(top, a):
    t = min(top, key=lambda r: abs(r['a'] - a))
    return t['basin_link_area_comoving_over_V']


def y(m):
    return float(erfinv(np.clip(m, -1 + 1e-12, 1 - 1e-12)))


pairs = [(51.2, f'fb_mu{mu:g}', FB, 'fb_mu0', FB) for mu in (0.02, 0.04, 0.07, 0.1)]
pairs += [(76.8, 'fb_N192_mu0.02', FB, 'gaussian_s101', RAW), (115.2, 'fb_N288_mu0.02', FB, 'gaussian_large_s101', RAW)]
out = []
for L, b, bdir, u, udir in pairs:
    if not (bdir / f'{b}.csv').exists():
        continue
    db, tb = series(bdir, b)
    du, tu = series(udir, u)
    if db['a'][-1] < 2.199:
        continue
    Kref = basin(tu, 1.6) / basin(tu, 2.2)
    rec = {'L': L, 'biased': b, 'unbiased_partner': u, 'K_ref_basin_1.6_to_2.2': Kref}
    for a in (1.6, 1.72, 1.84, 1.96, 2.08, 2.2):
        mb, mu_ = m_at(db, a), m_at(du, a)
        rec[f'a{a}'] = {'m_biased': mb, 'm_unbiased': mu_, 'y_response': y(mb) - y(mu_), 'f_minus_biased': float(
            db['fminus_05'][int(np.argmin(np.abs(db['a'] - a)))]), 'basin_biased': basin(tb, a), 'basin_unbiased': basin(tu, a)}
    r0, r1 = rec['a1.6']['y_response'], rec['a2.2']['y_response']
    rec['response_ratio_2.2_over_1.6'] = r1 / r0
    rec['gamma_eff_vs_Kref'] = float(np.log(r1 / r0) / np.log(Kref))
    yf = y(rec['a1.6']['m_biased'])
    rec['law_pred_m_end'] = float(erf(yf * Kref ** 1.41))
    rec['law_takeover_threshold_m_f'] = float(erf(1.54 * Kref ** -1.41))
    rec['observed_m_end'] = rec['a2.2']['m_biased']
    rec['observed_minority_end'] = rec['a2.2']['f_minus_biased']
    out.append(rec)
    print('L=%5.1f %-16s Kref=%.3f m_f=%.3f (partner %.3f) m_end=%.3f (partner %.3f) f-_end=%.4f resp_ratio=%.2f gamma_eff=%.2f law_m_end=%.3f thr=%.2f' % (
        L, b, Kref, rec['a1.6']['m_biased'], rec['a1.6']['m_unbiased'], r1 and rec['a2.2']['m_biased'], rec['a2.2']['m_unbiased'],
        rec['observed_minority_end'], rec['response_ratio_2.2_over_1.6'], rec['gamma_eff_vs_Kref'], rec['law_pred_m_end'],
        rec['law_takeover_threshold_m_f']))
OUT.write_text(json.dumps(out, indent=1))
