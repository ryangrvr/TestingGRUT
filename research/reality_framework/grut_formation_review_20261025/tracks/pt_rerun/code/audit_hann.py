"""Audit (c): normalisation, leakage and resolution of the Hann-windowed negative-frequency action.

Own code (Claude, pt_rerun). The estimator is re-implemented from MATHEMATICS.md:
  P_k(w_n) = |sum_j h_j psi_k(s_j) e^{-i w_n s_j}|^2 / (M sum_j h_j^2),  A_- = sum_{k, w_n<0} |w_n| P_k(w_n).
Parts: (i) synthetic free tones psi = A e^{+ims} + B e^{-ims} -> expect m|B|^2;
(ii) the exact homogeneous starting orbit (lab-frequency Fourier series 1.737788 per coherent charge);
(iii) the rerun IR12_s1_w1 psi.npy with other window lengths and tapers.
Usage: python3 -B audit_hann.py <early_packet_dir> <rerun workdir/IR12_s1_w1> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp

early, rundir, out = sys.argv[1:4]
LX = 1e-3
res = {}


def taper(M, kind):
    if kind == 'hann': return np.hanning(M)
    if kind == 'rect': return np.ones(M)
    if kind == 'blackman': return np.blackman(M)
    raise ValueError(kind)


def action(series, ds, kind='hann'):
    """series: (M, n) complex. Returns (A_minus, A_plus, total power)."""
    M = series.shape[0]; h = taper(M, kind); w = 2 * np.pi * np.fft.fftfreq(M, ds)
    X = np.fft.fft(h[:, None] * series, axis=0); P = np.sum(abs(X) ** 2, axis=1) / (M * np.sum(h * h))
    return float(np.sum((-w * P)[w < 0])), float(np.sum((w * P)[w > 0])), float(P.sum())


# ---------------------------------------------------------------- (i) synthetic free tones
rng = np.random.default_rng(20261027)
ds = .04; s = 14 + ds * np.arange(251)  # identical sampling to the packet's 14-24 window (M=251)
free = []
for m in [0.2, 0.4, 0.63, 1.0, 1.5, 2.0, 3.0, 4.5, 7.0, 9.2]:
    n = 400
    A = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2); B = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2) * .7
    psi = A[None, :] * np.exp(1j * m * s)[:, None] + B[None, :] * np.exp(-1j * m * s)[:, None]
    am, ap, _ = action(psi, ds)
    exact_m = m * np.sum(abs(B) ** 2); exact_p = m * np.sum(abs(A) ** 2)
    # single tone, B only, random phase offsets of the frequency relative to the DFT grid
    errs = []
    for off in np.linspace(0, 1, 9):
        mm = m + off * 2 * np.pi / (251 * ds)
        one = np.exp(-1j * mm * s)[:, None]
        errs.append(action(one, ds)[0] / mm - 1)
    free.append(dict(m=m, n_minus_estimator=am, n_minus_exact=exact_m, rel_error=am / exact_m - 1,
                     n_plus_rel_error=ap / exact_p - 1, single_tone_rel_error_range=[float(min(errs)), float(max(errs))]))
res['synthetic_free'] = dict(window_s=10.0, M=251, ds=ds, resolution=2 * np.pi / (251 * ds), rows=free,
                             note='n_- = m|B|^2 per mode; tones at |w| >~ 2 resolution widths (1.25) are recovered to <~1e-3 (pair statistics) ;'
                                  ' below that, leakage across w=0 biases the split')

# ---------------------------------------------------------------- (ii) homogeneous orbit through the windowed estimator
d = np.load(f'{early}/MATCH_L12.npz'); y = np.r_[d['light_mean_q'][:2], d['light_mean_p'][:2]]
Q = y[0] * y[3] - y[1] * y[2]


def rhs(t, v): return np.r_[v[2:], -LX * np.sum(v[:2] ** 2) * v[:2]]


T = 60.; tt = np.arange(0, T + 1e-9, ds)
sol = solve_ivp(rhs, (0, T), y, t_eval=tt, rtol=2e-13, atol=1e-11, method='DOP853')
psi0 = sol.y[0] + 1j * sol.y[1]
rows = []
for start in np.arange(0, 40, 0.37):
    i0 = int(round(start / ds)); seg = psi0[i0:i0 + 251, None]
    rows.append(action(seg, ds)[0] / Q)
rows = np.array(rows)
rows5 = np.array([action(psi0[int(round(st / ds)):int(round(st / ds)) + 126, None], ds)[0] / Q for st in np.arange(0, 40, 0.37)])
rows20 = np.array([action(psi0[int(round(st / ds)):int(round(st / ds)) + 501, None], ds)[0] / Q for st in np.arange(0, 35, 0.37)])
rect = np.array([action(psi0[int(round(st / ds)):int(round(st / ds)) + 251, None], ds, 'rect')[0] / Q for st in np.arange(0, 40, 0.37)])
res['homogeneous_orbit'] = dict(exact_lab_negative_action_per_coherent_charge=1.7377883531120057,
                                hann10_mean=float(rows.mean()), hann10_range=[float(rows.min()), float(rows.max())],
                                hann5_range=[float(rows5.min()), float(rows5.max())], hann20_range=[float(rows20.min()), float(rows20.max())],
                                rect10_range=[float(rect.min()), float(rect.max())],
                                note='window placed at 109 start phases; exact value from COHERENT_MATCHING.json (Fourier series over one radial period)')

# ---------------------------------------------------------------- (iii) rerun data: window length and taper sensitivity
meta = np.load(f'{rundir}/modes.npz'); psi = np.load(f'{rundir}/psi.npy', mmap_mode='r'); t = meta['s']; Qn = float(meta['nQ'])
kc = meta['kcont']; L = float(meta['L']); sid = np.floor(kc / (2 * np.pi / L) + .5).astype(int); zero = np.flatnonzero(sid == 0); waves = np.flatnonzero(sid > 0)


def window_action(lo, hi, kind='hann'):
    ix = np.flatnonzero((t >= lo - 1e-9) & (t <= hi + 1e-9))
    a = np.asarray(psi[ix], dtype=np.complex128)
    c = action(a[:, zero], ds, kind); wv = action(a[:, waves], ds, kind)
    return dict(window=[lo, hi], taper=kind, coherent=c[0] / Qn, waves=wv[0] / Qn, total=(c[0] + wv[0]) / Qn,
                charge_ratio=(c[1] + wv[1] - c[0] - wv[0]) / Qn)


sens = []
for lo, hi in [(14, 24), (24, 34)]:
    for kind in ['hann', 'blackman', 'rect']: sens.append(window_action(lo, hi, kind))
for lo, hi in [(14, 19), (19, 24), (24, 29), (29, 34), (14, 34), (19, 29)]:
    sens.append(window_action(lo, hi, 'hann'))
res['rerun_window_sensitivity'] = sens

# power and action near w=0 (within two Hann bins) in the packet windows
near = []
for lo, hi in [(14, 24), (24, 34)]:
    ix = np.flatnonzero((t >= lo - 1e-9) & (t <= hi + 1e-9)); M = len(ix); h = np.hanning(M); w = 2 * np.pi * np.fft.fftfreq(M, ds)
    a = np.asarray(psi[ix], dtype=np.complex128)
    X = np.fft.fft(h[:, None] * a, axis=0); P = abs(X) ** 2 / (M * np.sum(h * h))
    Pw = P[:, waves].sum(1); Pc = P[:, zero].sum(1); band = abs(w) <= 2 * 2 * np.pi / (M * ds) + 1e-9
    near.append(dict(window=[lo, hi], resolution=2 * np.pi / (M * ds),
                     wave_power_fraction_abs_w_le_2bins=float(Pw[band].sum() / Pw.sum()),
                     wave_negative_action_from_w_in_minus2bins_to_0=float(np.sum((-w * Pw)[band & (w < 0)]) / Qn),
                     wave_positive_action_from_0_to_2bins=float(np.sum((w * Pw)[band & (w > 0)]) / Qn),
                     coherent_power_fraction_abs_w_le_2bins=float(Pc[band].sum() / Pc.sum()),
                     wave_negative_action_per_frequency=[[float(a_), float(b_)] for a_, b_ in zip(np.fft.fftshift(w), np.fft.fftshift(np.where(w < 0, -w * Pw, 0) / Qn)) if -14 < a_ < 0]))
res['near_zero_frequency'] = near
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != 'near_zero_frequency'}, indent=1))
for n in near: print({k: v for k, v in n.items() if k != 'wave_negative_action_per_frequency'})
