"""st_rerun_skeptic: definitional sensitivity of the 'unit modulation' threshold (own code, own engine outputs).

The fractional response dA = d|X|/|X| peaks at the X orbit's radial minima (periapsis), where |X| is ~11-12% of the
orbit's maximum.  Compare its first crossings with an orbit-scale measure  d|X| / max_{|s-t|<=8}|X|(s),
using |X|(t) from our own C background (bg_*.npz, dtau 0.01) and dA from our tangent engine (every other 0.005 sample).
Usage: python3 -B orbit_scale.py <scratch_dir> <out_json>
"""
import sys, json
from pathlib import Path
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
L, outp = Path(sys.argv[1]), Path(sys.argv[2]); out = {}
for tag in ('quiet3', 'quiet5'):
    m = np.load(L / f'{tag}_main_800_r1e-10_d0.005.npz'); bg = np.load(L / f'bg_{tag}.npz')
    t = m['t'][::2]; amp = np.abs(m['amp'][::2]); chg = np.abs(m['charge'][::2]); assert np.allclose(t, bg['tau'])
    X = bg['absX']; W = 800
    Xmax = sliding_window_view(np.pad(X, (W, W), mode='edge'), 2 * W + 1).max(axis=1)
    orb = amp * (X / Xmax)[:, None]
    f = lambda y, l: (float(t[np.flatnonzero(y >= l)[0]]) if (y >= l).any() else None)
    row = dict(min_X_over_orbitmax_100_150=float((X / Xmax)[(t > 100) & (t < 150)].min()),
               min_X_over_orbitmax_700_800=float((X / Xmax)[t > 700].min()))
    for s in (.03, 1e-6):
        e = s * amp.max(axis=1); o = s * orb.max(axis=1); c = s * chg.max(axis=1)
        row[f'seed_{s:g}'] = dict(fractional_unity=f(e, 1.), orbit_scale_unity=f(o, 1.), orbit_scale_0p3=f(o, .3),
                                  orbit_scale_0p1=f(o, .1), charge_unity=f(c, 1.), orbit_scale_max_by_800=float(o.max()))
    out[tag] = row
    print(tag, json.dumps(row))
outp.write_text(json.dumps(out, indent=1) + '\n')
