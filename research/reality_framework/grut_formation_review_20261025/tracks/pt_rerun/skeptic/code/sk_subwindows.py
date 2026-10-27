"""Skeptic of pt_rerun: is the late decline of the k<=8 statistic really 'cascade into k>8'? (own code, Claude)

Rerun IR12_s1_w1 psi.npy, 5-unit Hann windows (M=126) stepped by 2.5 in s from 14 to 34; for each window the k<=8
coherent and wave negative-frequency action (explicit FFT, own code) plus the Hann^2-weighted equal-time k>8 tail
(trace columns tail_H and tail_0) as published by pt_production. Reports the drift of (k<=8) and (k<=8 + tail).
Usage: python3 -B sk_subwindows.py <rerun dir IR12_s1_w1> <out.json>
"""
import sys, json
import numpy as np

rundir, out = sys.argv[1:3]
meta = np.load(f'{rundir}/modes.npz'); psi = np.load(f'{rundir}/psi.npy', mmap_mode='r'); s = meta['s']; Q = float(meta['nQ'])
kc = meta['kcont']; zero = kc == 0
tr = np.load(f'{rundir}/trace.npz'); names = [str(x) for x in tr['names']]; V = tr['values']
tH = V[:, names.index('tail_H')]; t0 = V[:, names.index('tail_0')]
rows = []
for lo in np.arange(14., 29.01, 2.5):
    hi = lo + 5.; ix = np.flatnonzero((s >= lo - 1e-9) & (s <= hi + 1e-9)); M = len(ix); ds = s[1] - s[0]
    h = np.hanning(M); w = 2 * np.pi * np.fft.fftfreq(M, ds); a = np.asarray(psi[ix], dtype=np.complex128)
    P = abs(np.fft.fft(h[:, None] * a, axis=0)) ** 2 / (M * np.sum(h * h))
    coh = float(np.sum((-w * P[:, zero].sum(1))[w < 0]) / Q); wav = float(np.sum((-w * P[:, ~zero].sum(1))[w < 0]) / Q)
    wq = h * h / np.sum(h * h)
    rows.append(dict(window=[float(lo), float(hi)], coherent=coh, waves=wav, k_le_8=coh + wav, tail_H=float(wq @ tH[ix]), tail_0=float(wq @ t0[ix]),
                     with_tail=coh + wav + float(wq @ tH[ix])))
    print(rows[-1], flush=True)
a = np.array([r['k_le_8'] for r in rows]); b = np.array([r['with_tail'] for r in rows])
res = dict(rows=rows, k_le_8_first_to_last_percent=float(100 * (a[-1] / a[0] - 1)), with_tail_first_to_last_percent=float(100 * (b[-1] / b[0] - 1)),
           k_le_8_range=[float(a.min()), float(a.max())], with_tail_range=[float(b.min()), float(b.max())])
json.dump(res, open(out, 'w'), indent=1)
print({k: v for k, v in res.items() if k != 'rows'})
