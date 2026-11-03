"""Skeptic S2: split lf_fermion's output-level lambda_PB into the common-mode (homogeneous-control) closure deficit and
the locality part. Read-only on results/cmp_REC_P0.json. Claude, lf_fermion skeptic (issue #13 context).
lambda_out = (P0 - Q)/(P0 - GM) (lf_fermion). r_hom = Q_hom/closure_hom (same line-mean history, P0 = GM there).
lambda_loc = (r_hom P0 - Q)/(r_hom P0 - r_hom GM): both closures first corrected by the common-mode factor [model choice].
Also with Q(N1) divided by 1.15 (sk_measure_bias: N1 Q biased high by 14-23% after t = 0.4).
Usage: python3 -B sk_lambda_split.py -> ../results/sk_lambda_split.json"""
import json
from pathlib import Path
import numpy as np
TR = Path(__file__).resolve().parents[2]
c = json.load(open(TR / 'results' / 'cmp_REC_P0.json'))
tg = np.array(c['t'])
out = {}
for nm, d in c['species'].items():
    q = np.array(d['Q_mean']); p0 = np.array(d['tracer_mean']['0.0']); gm = np.array(d['tracer_mean']['1.0'])
    h = c['homog'][nm]; r = np.array(h['Q_hom']) / np.array(h['closure_hom'])
    row = {}
    for tt in (0.6, 0.7, 0.9, 1.1):
        i = int(np.argmin(abs(tg - tt)))
        lam = (p0[i] - q[i]) / (p0[i] - gm[i])
        lam_loc = (r[i] * p0[i] - q[i]) / (r[i] * (p0[i] - gm[i]))
        qc = q[i] / 1.15 if nm == 'N1' else q[i]
        row[tt] = dict(Q_over_P0=float(q[i] / p0[i]), r_hom=float(r[i]), extra_factor_over_hom=float(q[i] / p0[i] / r[i]),
                       lambda_out=float(lam), lambda_loc=float(lam_loc),
                       lambda_loc_N1corr=float((r[i] * p0[i] - qc) / (r[i] * (p0[i] - gm[i]))))
    out[nm] = row
    print(nm, ' '.join('t%.1f: Q/P0 %.2f r_hom %.2f extra %.2f lam_out %.2f lam_loc %.2f (N1corr %.2f)' % (
        tt, v['Q_over_P0'], v['r_hom'], v['extra_factor_over_hom'], v['lambda_out'], v['lambda_loc'], v['lambda_loc_N1corr'])
        for tt, v in row.items() if tt in (0.7, 1.1)))
(TR / 'skeptic' / 'results' / 'sk_lambda_split.json').write_text(json.dumps(dict(label='C', note=__doc__, species=out), indent=1))
