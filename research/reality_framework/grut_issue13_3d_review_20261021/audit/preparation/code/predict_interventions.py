"""Pre-registered quasi-linear predictions for the s101 interventions (written BEFORE the runs).

Usage: python3 -I predict_interventions.py <codex_raw_dir> <out_json>
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sparse_diagnostics as sd  # noqa: E402

RAW, OUT = Path(sys.argv[1]), Path(sys.argv[2])
base = np.loadtxt(RAW / 'modes_L76.8_s101.txt')
res = {'written_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'rule': (
    'Outcome classes from the linear field at total rms 0.5 (m_lin = <sign S_lin> incl. sourced zero mode): '
    '|m_lin| >= 0.45 -> one-sign selection before formation (no 10%/10% formation); '
    '0.1 <= |m_lin| < 0.45 -> biased two-sign network, sign = sign(m_lin); |m_lin| < 0.1 -> near-unbiased two-sign. '
    'Threshold 0.45 chosen from the 4 Codex sparse charged runs (s101 0.52 one-sign; kappa4 0.21 two-sign), so it is '
    'a calibrated classifier, not a derived constant.'), 'variants': {}}
for name in ['baseline', 'flipX', 'noX', 'gapfill', 'noU', 'noW']:
    t = base.copy()
    if name == 'noX':
        t[:, 3:7] = 0
    elif name == 'flipX':
        t[:, 3:7] *= -1
    elif name == 'noW':
        t[:, 5:7] = 0
    elif name == 'noU':
        t[:, 3:5] = 0
    elif name == 'gapfill':
        t = np.vstack([t, [1, 0, 0, 0, 0, 0, 0, 0.125, 0.125]])
    n, c = t[:, :3].astype(int), t[:, 3:]
    r = sd.analyze_arrays(name, 'modified s101', n, c, 192, 0.4, 2.0, None)
    p = r['predictions']['rms_0.5']
    m = p['mean_sign_with_zero_mode']
    cls = 'one-sign' if abs(m) >= 0.45 else ('biased two-sign' if abs(m) >= 0.1 else 'near-unbiased two-sign')
    res['variants'][name] = {'predicted_class': cls, 'predicted_sign': int(np.sign(m)) if abs(m) >= 0.1 else 0,
                             'm_lin_rms0.5': m, 'S0_over_sigma_f_rms0.5': p['S0_over_sigma_f'], 'a_rms0.5': p['a'],
                             'predictions': r['predictions'], 'C_uS': r['C_uS'], 'C_wS': r['C_wS'], 'k_min': r['k_min'],
                             'n_triads_ordered': r['n_triads_ordered']}
    print(name, cls, 'm_lin=%.3f S0/sig=%.3f a=%.3f' % (m, p['S0_over_sigma_f'], p['a']), flush=True)
OUT.write_text(json.dumps(res, indent=1))
