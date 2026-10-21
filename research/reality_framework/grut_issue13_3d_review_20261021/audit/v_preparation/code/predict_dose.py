"""Pre-registered predictions (written BEFORE the runs) for s101 with the X perturbation scaled by s (S unchanged).
Quasi-linear: S0 scales exactly with s, sigma_f unchanged. Usage: python3 -I predict_dose.py <raw> <out_json>"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import leapfrog_ql as lq
raw, outp = Path(sys.argv[1]), Path(sys.argv[2])
t = np.loadtxt(raw / 'modes_L76.8_s101.txt'); n = t[:, :3].astype(int); c = t[:, 3:]
res = {'written_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
       'rule': ('Class by the audit track\'s empirical kappa=2 sparse thresholds: S0/sigma_f(rms 0.5) <= 0.13 -> two-sign network '
                'retained at a=2.2 (both vacua > 1%, sign of bias = sign S0); >= 0.34 -> minority removed by a=2.2; between: '
                'untested, verifier predicts network retained if < 0.25. Formation m_f (first f+ + f- >= 0.5) predicted '
                'as m_lin x [1.4, 2.5] (the audit\'s sparse under-prediction band).'), 'variants': {}}
for s in (0.1, 0.3):
    o, k2 = lq.evolve(n, c, 192, 0.4, 2.0, xscale=(s, s))
    tot = np.sqrt(o['sigma_f'] ** 2 + o['S0'] ** 2); j = int(np.argmax(tot >= 0.5))
    oo, _ = lq.evolve(n, c, 192, 0.4, 2.0, T=o['t'][j], xscale=(s, s))
    st = lq.field_sign_stats(n, oo['_final_ch'], oo['S0'][-1], oo['a'][-1], 192, sub=2)
    r = o['S0'][j] / o['sigma_f'][j]; m = st['mean_sign_total']
    cls = 'network retained' if abs(r) < 0.25 else ('minority removed' if abs(r) >= 0.34 else 'untested')
    res['variants'][f'Xscale{s:g}'] = {'a_rms05': float(o['a'][j]), 'S0': float(o['S0'][j]), 'sigma_f': float(o['sigma_f'][j]),
        'S0_over_sigma_f': float(r), 'm_lin_rms05': m, 'm_f_pred_band': [1.4 * m, 2.5 * m], 'predicted_class': cls,
        'S0_pred_a1.54': float(o['S0'][int(round(0.54 * 40 / 0.06))]), 'S0_pred_a1.57': float(o['S0'][int(round(0.57 * 40 / 0.06))])}
    print(s, res['variants'][f'Xscale{s:g}'])
outp.write_text(json.dumps(res, indent=1))
