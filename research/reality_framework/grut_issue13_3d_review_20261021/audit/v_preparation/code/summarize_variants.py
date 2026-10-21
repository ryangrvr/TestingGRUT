"""Summarize the verifier's Codex-binary variant runs; compare flipX with the track's record; copy small records.
Usage: python3 -I summarize_variants.py <runs_dir> <track_rerun_series_dir> <out_json>"""
import json, shutil, sys
from pathlib import Path
import numpy as np
runs, track, outp = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
pred = json.loads((outp.parent / 'dose_predictions_preregistered.json').read_text())['variants']
dest = outp.parent / 'variant_series'; dest.mkdir(exist_ok=True)
res = {}
for name in ('flipX', 'Xscale0.1', 'Xscale0.3'):
    p = runs / f'v_{name}.csv'
    if not (runs / f'v_{name}.meta.json').exists():
        res[name] = 'not finished'; continue
    for ext in ('.csv', '.topology.jsonl', '.meta.json'):
        shutil.copy(runs / f'v_{name}{ext}', dest / f'v_{name}{ext}')
    shutil.copy(runs / f'modes_v_{name}.txt', dest / f'modes_v_{name}.txt')
    d = np.genfromtxt(p, delimiter=',', names=True); t = [json.loads(l) for l in open(runs / f'v_{name}.topology.jsonl')]
    fp, fm = d['fplus_05'], d['fminus_05']; s = np.flatnonzero(fp + fm >= 0.5); i0 = int(s[0]) if len(s) else None
    r = {'end': {'a': float(d['a'][-1]), 'fplus': float(fp[-1]), 'fminus': float(fm[-1]), 'genuine_plus': float(d['genuine_plus'][-1]),
                 'genuine_minus': float(d['genuine_minus'][-1]), 'm': float((fp[-1] - fm[-1]) / (fp[-1] + fm[-1])),
                 'plus_winding': t[-1]['plus']['any_winding'], 'minus_winding': t[-1]['minus']['any_winding'],
                 'basin_area': t[-1]['basin_link_area_comoving_over_V']},
         'formed_10_10': bool(np.any((fp >= 0.1) & (fm >= 0.1))),
         'first_saturation': {'a': float(d['a'][i0]), 'm': float((fp[i0] - fm[i0]) / (fp[i0] + fm[i0]))} if i0 is not None else None,
         'mean_S_a1.54': float(d['mean_S'][int(np.argmin(abs(d['a'] - 1.54)))]), 'mean_S_a1.57': float(d['mean_S'][int(np.argmin(abs(d['a'] - 1.57)))]),
         'minority_min_after_1.7': float(min(np.min(fp[d['a'] > 1.7]), np.min(fm[d['a'] > 1.7])))}
    r['series_m'] = [[round(float(a), 3), round(float((x - y) / max(x + y, 1e-300)), 5), round(float(x), 5), round(float(y), 5)] for a, x, y in zip(d['a'], fp, fm) if a >= 1.6]
    if name in pred:
        r['preregistered'] = pred[name]
        r['mean_S_rel_dev_a1.54'] = r['mean_S_a1.54'] / pred[name]['S0_pred_a1.54'] - 1
        r['mean_S_rel_dev_a1.57'] = r['mean_S_a1.57'] / pred[name]['S0_pred_a1.57'] - 1
    if name == 'flipX':
        b = np.genfromtxt(track / 'int_flipX.csv', delimiter=',', names=True)
        r['max_abs_diff_vs_track_int_flipX'] = float(max(np.max(np.abs(d[c] - b[c])) for c in d.dtype.names))
        r['topology_identical_to_track'] = t == [json.loads(l) for l in open(track / 'int_flipX.topology.jsonl')]
        r['mode_table_identical_to_track'] = (runs / 'modes_v_flipX.txt').read_bytes() == (track / 'modes_L76.8_s101_flipX.txt').read_bytes()
    res[name] = r
    print(name, json.dumps({k: v for k, v in r.items() if k not in ('series_m', 'preregistered')}))
outp.write_text(json.dumps(res, indent=1))
