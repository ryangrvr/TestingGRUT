"""Does rigid core translation account for the old canonical relative-shift amplification?
Uses the old code's own R8 local translation norms (amplitude+S only). Usage: python -I rigid_norm_estimate.py <raw_dir>"""
import json, sys, numpy as np
R = sys.argv[1]
for h in [0.5, 1/3, 0.25]:
    tr = json.load(open(f'{R}/large_modes_h{h:.6f}.json'))['trace']
    c0 = tr[0]['late_core_shape']; c1 = tr[-1]['late_core_shape']
    G0 = [c['translation_norm'] for c in c0]
    num = sum(c['fitted_displacement']**2*c['translation_norm'] for c in c1)
    raw_end = sum(c['raw'] for c in c1); raw0 = sum(c['raw'] for c in c0)
    print('h=%.4f initial core translation norms %s initial fitted disp %s | rigid amplitude+S amplification sqrt(sum X^2 G_end / sum G0)=%.3f ; raw local ratio sqrt(raw_end/raw0)=%.3f ; reported canonical amp=%.3f' % (
        h, np.round(G0, 3), np.round([c['fitted_displacement'] for c in c0], 4), np.sqrt(num/sum(G0)), np.sqrt(raw_end/raw0), tr[-1]['norm'][2]/tr[0]['norm'][2]))
