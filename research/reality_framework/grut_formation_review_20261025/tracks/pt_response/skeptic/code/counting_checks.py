"""SKEPTIC (pt_response_skeptic), INDEPENDENT checks on published Codex arrays (5f5513b), no Codex code.

1. S = A+ + A- per net charge over all 18 WINDOW files, from two routes:
   (a) S/Q = spectral_charge_ratio + 2 * total_negative (the pt_response route),
   (b) per-shell positive + negative wave actions + coherent part recomputed from FIELD_TRACES (own Hann FFT).
2. Share of the wave negative action in shells 1-4 (k <= 2.1) and wave charge (A+ - A-) by shell.
3. Published response arrays: conjugate/normal RMS ratio per shell, seed 7 vs 11 differences.
4. Radial action J_r of the coherent mode from FIELD_TRACES by a DIFFERENT method from pt_response:
   cubic-spline resampling (ds = 0.002) and period segmentation at r MINIMA, loop integral (1/2pi) oint p_r dr,
   compared with the per-maximum segmentation and with our own Hann lab statistic in the same window.
Usage: python3 -I counting_checks.py <packet dir> <out json>
"""
import sys, json, glob, os
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import argrelmin, argrelmax

P = sys.argv[1].rstrip('/') + '/'; OUT = sys.argv[2]
res = {}

# ---------- 1/2: windows
F = np.load(P + 'FIELD_TRACES.npz'); names = list(F['names'])
wins = sorted(glob.glob(P + 'WINDOW_*.json'))
rows = []
for w in wins:
    d = json.load(open(w)); tag = d['tag']; lo, hi = d['window_s']
    sh = [x for x in d['shells'] if x['shell'] > 0]
    wm = sum(x['negative_frequency_action_per_net'] for x in sh); wp = sum(x['positive_frequency_action_per_net'] for x in sh)
    # own coherent lab statistic from the trace (Hann, uniform samples within window)
    s = F[f'{tag}_s']; v = F[f'{tag}_values']
    sel = (s >= lo - 1e-9) & (s <= hi + 1e-9); ss = s[sel]
    psi = v[sel, names.index('q0')] + 1j * v[sel, names.index('q1')]
    nQ = float(np.mean(v[sel, names.index('nQ')]))
    n = len(ss); h = np.hanning(n); ds = ss[1] - ss[0]
    X = np.fft.fft(h * psi); Pw = np.abs(X) ** 2 / (n * np.sum(h * h)); om = 2 * np.pi * np.fft.fftfreq(n, ds)
    cm = float(np.sum((-om * Pw)[om < 0]) / nQ); cp = float(np.sum((om * Pw)[om > 0]) / nQ)
    S_a = d['spectral_charge_ratio'] + 2 * d['total_negative_frequency_action']
    S_b = wm + wp + cm + cp
    low = sum(x['negative_frequency_action_per_net'] for x in sh if x['k_center'] <= 2.1 + 1e-9)
    rows.append(dict(tag=tag, window=[lo, hi], A_minus_total=d['total_negative_frequency_action'],
                     coherent_A_minus_codex=d['coherent_negative_frequency_action'], coherent_A_minus_ours=cm, coherent_A_plus_ours=cp,
                     wave_A_minus=wm, wave_A_plus=wp, wave_charge_over_Q=wp - wm, coherent_charge_over_Q_spectral=cp - cm,
                     S_over_Q_route_a=S_a, S_over_Q_route_b=S_b, charge_identity=(wp + cp) - (wm + cm),
                     wave_A_minus_share_k_le_2p1=low / wm, wave_A_minus_share_k_le_1p6=sum(x['negative_frequency_action_per_net'] for x in sh if x['k_center'] <= 1.6) / wm,
                     n_samples=n))
res['windows'] = rows
Sa = [r['S_over_Q_route_a'] for r in rows]; Sb = [r['S_over_Q_route_b'] for r in rows]
res['S_range_route_a'] = [min(Sa), max(Sa)]; res['S_range_route_b'] = [min(Sb), max(Sb)]
res['coherent_A_minus_max_abs_diff_vs_codex'] = max(abs(r['coherent_A_minus_ours'] - r['coherent_A_minus_codex']) for r in rows)
print('S/Q route a %.4f-%.4f  route b %.4f-%.4f  (n=%d windows)' % (min(Sa), max(Sa), min(Sb), max(Sb), len(rows)))
for r in rows:
    print('  %-22s %5.0f-%-3.0f A- %.4f  coh %.4f/%.4f  waves A- %.4f A+ %.4f (Qw %.3f)  share k<=2.1 %.3f  S %.4f/%.4f' % (
        r['tag'], *r['window'], r['A_minus_total'], r['coherent_A_minus_codex'], r['coherent_A_minus_ours'], r['wave_A_minus'], r['wave_A_plus'],
        r['wave_charge_over_Q'], r['wave_A_minus_share_k_le_2p1'], r['S_over_Q_route_a'], r['S_over_Q_route_b']))

# ---------- 3: response arrays
resp = {}
for f in sorted(glob.glob(P + 'RESPONSE_*.npz')):
    d = np.load(f); tag = os.path.basename(f)[9:-4]
    nr = d['normal']; an = d['anomalous']
    resp[tag] = dict(ratio=[float(np.linalg.norm(an[:, k]) / np.linalg.norm(nr[:, k])) for k in range(nr.shape[1])], shells=d['shell'].tolist())
res['response_ratio_by_shell'] = resp
d7 = np.load(P + 'RESPONSE_s20_seed7.npz'); d11 = np.load(P + 'RESPONSE_s20_seed11.npz')
res['seed7_vs_11_normal_reldiff'] = [float(np.linalg.norm(d7['normal'][:, k] - d11['normal'][:, k]) / np.linalg.norm(d7['normal'][:, k])) for k in range(15)]
res['seed7_vs_11_anomalous_reldiff'] = [float(np.linalg.norm(d7['anomalous'][:, k] - d11['anomalous'][:, k]) / np.linalg.norm(d7['anomalous'][:, k])) for k in range(15)]
for k, v in resp.items():
    print(k, np.round(v['ratio'], 3))
print('seed diff normal', np.round(res['seed7_vs_11_normal_reldiff'], 3))

# ---------- 4: radial action
s = F['base48_s']; v = F['base48_values']
q0, q1, p0, p1 = [v[:, names.index(x)] for x in ('q0', 'q1', 'p0', 'p1')]
nQ = float(v[0, names.index('nQ')])

def jr_window(lo, hi, seg='min'):
    sel = (s >= lo - 1.5) & (s <= hi + 1.5)
    sp = {k: CubicSpline(s[sel], a[sel]) for k, a in (('q0', q0), ('q1', q1), ('p0', p0), ('p1', p1))}
    sf = np.arange(max(lo - 1.0, s[0]), min(hi + 1.0, s[-1]) + 1e-9, .002)
    a, b, pa, pb = [sp[k](sf) for k in ('q0', 'q1', 'p0', 'p1')]
    r = np.hypot(a, b); pr = (a * pa + b * pb) / r
    ext = argrelmin(r)[0] if seg == 'min' else argrelmax(r)[0]
    J = []; mids = []
    for i, j in zip(ext[:-1], ext[1:]):
        rr = r[i:j + 1]; pp = pr[i:j + 1]
        # loop integral (1/2pi) oint p_r dr by trapezoid in r (parametric), robust to sampling
        J.append(np.sum(.5 * (pp[1:] + pp[:-1]) * np.diff(rr)) / (2 * np.pi)); mids.append(.5 * (sf[i] + sf[j]))
    J = np.array(J); mids = np.array(mids)
    inwin = (mids >= lo) & (mids <= hi)
    wts = np.sin(np.pi * (mids[inwin] - lo) / (hi - lo)) ** 2
    return dict(Jr_mean_per_Q=float(J[inwin].mean() / nQ), Jr_hann2_per_Q=float(np.sum(wts * J[inwin]) / np.sum(wts) / nQ),
                periods=int(inwin.sum()), first_last_per_Q=[float(J[inwin][0] / nQ), float(J[inwin][-1] / nQ)])

jr = {}
for lo, hi in ((14., 24.), (24., 34.), (0., 6.)):
    a_ = jr_window(lo, hi, 'min'); b_ = jr_window(lo, hi, 'max')
    lab = [r for r in rows if r['tag'] == 'base48' and r['window'] == [lo, hi]][0]['coherent_A_minus_codex']
    jr[f'{lo:g}_{hi:g}'] = dict(seg_min=a_, seg_max=b_, codex_lab_coherent=lab,
                               excess_hann_min=a_['Jr_hann2_per_Q'] / lab - 1, excess_hann_max=b_['Jr_hann2_per_Q'] / lab - 1)
    print('J_r %g-%g: min-seg %.4f (Hann %.4f)  max-seg %.4f (Hann %.4f)  lab %.4f  excess %.2f%% / %.2f%%  periods %d' % (
        lo, hi, a_['Jr_mean_per_Q'], a_['Jr_hann2_per_Q'], b_['Jr_mean_per_Q'], b_['Jr_hann2_per_Q'], lab,
        100 * jr[f'{lo:g}_{hi:g}']['excess_hann_min'], 100 * jr[f'{lo:g}_{hi:g}']['excess_hann_max'], a_['periods']))
res['radial_action'] = jr
json.dump(res, open(OUT, 'w'), indent=1)
