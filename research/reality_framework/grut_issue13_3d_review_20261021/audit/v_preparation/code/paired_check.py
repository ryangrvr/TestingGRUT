"""Recompute the paired forced-bias gamma_eff (biased fb run minus unbiased partner) and test its sensitivity to
the start sample a0 in {1.57, 1.60, 1.72, 1.84}. gamma_eff = ln[(y_b-y_0)(2.2)/(y_b-y_0)(a0)] / ln[A_0(a0)/A_0(2.2)],
y = erfinv(m05), A_0 = unbiased partner's basin link area. Usage: python3 -I paired_check.py <codex_raw> <audit_fb_dir> <out>"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.special import erfinv
raw, fb, outp = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
def load(d, n):
    c = np.genfromtxt(d / f'{n}.csv', delimiter=',', names=True); t = [json.loads(l) for l in open(d / f'{n}.topology.jsonl')]
    return c, t
def m(c, a):
    i = int(np.argmin(abs(c['a'] - a))); return float((c['fplus_05'][i] - c['fminus_05'][i]) / (c['fplus_05'][i] + c['fminus_05'][i])), float(c['a'][i])
def A(t, a):
    r = min(t, key=lambda r: abs(r['a'] - a)); return r['basin_link_area_comoving_over_V'], r['a']
res = {}
for L, b, bd, u, ud in ((51.2, 'fb_mu0.02', fb, 'fb_mu0', fb), (76.8, 'fb_N192_mu0.02', fb, 'gaussian_s101', raw), (115.2, 'fb_N288_mu0.02', fb, 'gaussian_large_s101', raw)):
    cb, tb = load(bd, b); cu, tu = load(ud, u)
    out = {'topology_a_samples': [round(r['a'], 3) for r in tu]}
    yend = erfinv(m(cb, 2.2)[0]) - erfinv(m(cu, 2.2)[0])
    for a0 in (1.57, 1.6, 1.72, 1.84):
        y0 = erfinv(m(cb, a0)[0]) - erfinv(m(cu, a0)[0]); A0, aA = A(tu, a0); A1, _ = A(tu, 2.2)
        K = A0 / A1
        out[f'a0={a0}'] = {'csv_a': m(cb, a0)[1], 'topology_a': aA, 'm_b': m(cb, a0)[0], 'm_0': m(cu, a0)[0], 'y_resp': float(y0), 'K_ref': K,
                           'resp_ratio': float(yend / y0), 'gamma_eff': float(np.log(yend / y0) / np.log(K)) if K > 1 else None}
    res[str(L)] = out
    print(L, {k: (v if k == 'topology_a_samples' else {kk: round(vv, 3) if isinstance(vv, float) else vv for kk, vv in v.items()}) for k, v in out.items()})
outp.write_text(json.dumps(res, indent=1))
