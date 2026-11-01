"""INDEPENDENT analysis (Claude, pt_response): damped-pole fits of the normal-channel retarded response.

G_N(t) ~ sum_{j=1..n} c_j exp(-i w_j t - G_j t), t in [0, 6], complex c_j, real w_j, G_j >= 0 (lab convention of
response.py, whose transform uses e^{-i w t}; a free field sin(w t)/w has poles at +-w).
Fitted to (a) our low-noise measured traces (Hutchinson M=8 shells 1-4, exact shell 1) and Codex's published traces
(shells 5-15, where probe noise is <= 14%), and (b) the homogeneous mean-Hessian (Floquet-HFB) model traces, which
have no intrinsic damping: their fitted widths are the baseline from Floquet sidebands / finite window.
Usage: python3 -I fit_poles.py <probes dir> <hfb series npz> <packet dir> <out json>
"""
import sys, json, glob, os
import numpy as np
from scipy.optimize import least_squares

D = sys.argv[1].rstrip('/') + '/'; S = np.load(sys.argv[2]); P = sys.argv[3].rstrip('/') + '/'; OUT = sys.argv[4]
t = np.arange(151) * .04

def basis(x, n):
    return np.array([np.exp(-1j * x[2 * j] * t - abs(x[2 * j + 1]) * t) for j in range(n)]).T

def resid(x, G, n):
    B = basis(x, n); c, *_ = np.linalg.lstsq(B, G, rcond=None); r = B @ c - G
    return np.r_[r.real, r.imag]

def fit(G, n=2):
    """Variable projection: amplitudes linear, (w_j, Gamma_j) nonlinear; starts from transform peaks."""
    w = np.linspace(-15, 15, 1501); taper = np.cos(np.pi * t / 12) ** 2
    rho = -2 * (np.exp(-1j * w[:, None] * t[None, :]) @ (G * taper * .04)).imag
    wp = w[w > 0][np.argmax(rho[w > 0])]; wn = w[w < 0][np.argmax(-rho[w < 0])]
    starts = [[wp, .2, wn, .2], [wp, .02, wn, .02], [wp, .6, wn, .6]]
    if n == 4: starts = [s + [wp + 2., .3, wn - 2., .3] for s in starts] + [s + [wp - 1.5, .3, wn + 1.5, .3] for s in starts]
    best = None
    for x0 in starts:
        r = least_squares(resid, x0, args=(G, n), max_nfev=3000)
        if best is None or r.cost < best.cost: best = r
    x = best.x; B = basis(x, n); c, *_ = np.linalg.lstsq(B, G, rcond=None)
    res = np.linalg.norm(B @ c - G) / np.linalg.norm(G)
    poles = sorted([dict(w=float(x[2 * j]), Gamma=float(abs(x[2 * j + 1])), amp=float(abs(c[j]))) for j in range(n)], key=lambda p: p['w'])
    return dict(rel_residual=float(res), poles=poles)

traces = {}
for f in sorted(glob.glob(D + 'MEAS_hutch_*.npz')):
    d = np.load(f); s0 = float(d['s0'])
    for k in range(4): traces[(f'{s0:g}', k + 1, 'ours_hutch')] = d['normal'][:, k]
for f in sorted(glob.glob(D + 'MEAS_exact1_*.npz')):
    d = np.load(f); traces[(f"{float(d['s0']):g}", 1, 'ours_exact')] = d['normal'][:, 0]
for f in sorted(glob.glob(P + 'RESPONSE_s*_seed*.npz')):
    d = np.load(f); s0 = float(d['s0']); tag = os.path.basename(f)[9:-4]
    for k in range(4, 15): traces[(f'{s0:g}', k + 1, 'codex_' + tag)] = d['normal'][:, k]
res = []
for (s0, sh, src), G in sorted(traces.items()):
    row = dict(s0=s0, shell=sh, source=src)
    for n in (2, 4):
        row[f'meas_{n}pole'] = fit(G, n)
        key = f'{s0}_meanH_{sh}_N'
        if key in S.files: row[f'meanH_{n}pole'] = fit(S[key], n)
    res.append(row)
    m2 = row['meas_2pole']; h2 = row.get('meanH_2pole')
    print(s0, sh, src, '| meas 2p res %.3f' % m2['rel_residual'], ' '.join('(w %.2f G %.2f)' % (p['w'], p['Gamma']) for p in m2['poles']),
          '| 4p res %.3f' % row['meas_4pole']['rel_residual'],
          '| meanH 2p res %.3f' % h2['rel_residual'] if h2 else '', ' '.join('(w %.2f G %.2f)' % (p['w'], p['Gamma']) for p in h2['poles']) if h2 else '', flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
