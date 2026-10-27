"""Skeptic of pt_rerun: the rerun IR12_s1_w1 budget by an explicit (non-FFT) DFT and a zero-padded grid (own code, Claude).

* X_n = sum_j h_j psi_k(s_j) exp(-i w_n s_j) by an explicit complex matrix product (absolute times s_j, so the phase
  convention differs from the FFT's), P = |X|^2/(M sum h^2), A_- = sum_{w_n<0} |w_n| P, per net charge. Shells are not
  used; the coherent zero mode is identified by |k| = 0.
* The same statistic on a 4x zero-padded frequency grid (sum divided by 4) tests the discreteness of the |w| weight.
* Q0 statistics from the trace (does the zero-mode charge approach zero, where the radial-orbit picture degenerates?).
Usage: python3 -B sk_spectrum.py <rerun dir IR12_s1_w1> <published WINDOW json 14_24> <WINDOW json 24_34> <FIELD_TRACES.npz> <out.json>
"""
import sys, json
import numpy as np

rundir, w1, w2, traces, out = sys.argv[1:6]
meta = np.load(f'{rundir}/modes.npz'); psi = np.load(f'{rundir}/psi.npy', mmap_mode='r'); s = meta['s']; Q = float(meta['nQ'])
kc = meta['kcont']; zero = np.flatnonzero(kc == 0); waves = np.flatnonzero(kc > 0)
res = dict(n_modes=int(len(kc)), n_zero=int(len(zero)), nQ=Q)
for (lo, hi), wj in (((14., 24.), w1), ((24., 34.), w2)):
    ix = np.flatnonzero((s >= lo - 1e-9) & (s <= hi + 1e-9)); M = len(ix); ds = s[1] - s[0]; h = np.hanning(M)
    a = np.asarray(psi[ix], dtype=np.complex128)  # (M, modes)
    pub = json.load(open(wj))
    row = dict(window=[lo, hi], M=M)
    for pad in (1, 4):
        Mf = M * pad; w = 2 * np.pi * np.fft.fftfreq(Mf, ds)
        Emat = np.exp(-1j * np.outer(w, s[ix]))  # (Mf, M) explicit DFT with absolute time
        X = Emat @ (h[:, None] * a)
        P = abs(X) ** 2 / (M * np.sum(h * h)) / pad
        coh = float(np.sum((-w * P[:, zero].sum(1))[w < 0]) / Q); wav = float(np.sum((-w * P[:, waves].sum(1))[w < 0]) / Q)
        chg = float((np.sum(w * P.sum(1))) / Q)
        row[f'pad{pad}'] = dict(coherent=coh, waves=wav, total=coh + wav, spectral_charge_ratio=chg)
    row['published'] = dict(coherent=pub['coherent_negative_frequency_action'], waves=pub['wave_negative_frequency_action'],
                            total=pub['total_negative_frequency_action'], spectral_charge_ratio=pub['spectral_charge_ratio'])
    row['pad1_minus_published_total'] = row['pad1']['total'] - pub['total_negative_frequency_action']
    row['pad4_minus_pad1_total'] = row['pad4']['total'] - row['pad1']['total']
    res[f'{lo:g}_{hi:g}'] = row
    print(json.dumps(row), flush=True)
Z = np.load(traces); names = [str(x) for x in Z['names']]; c = {n: i for i, n in enumerate(names)}
q0 = {}
for tag in ['base48', 'IR12_s1_w1', 'IR12_s2_w1', 'IR12_s1_w3', 'IR12_fine_s1_w1', 'IR24_s1_w1', 'IR24_s2_w1']:
    V = Z[tag + '_values']; ss = Z[tag + '_s']; r = V[:, c['Q0_over_Q']]; sel = ss >= 14 - 1e-9
    q0[tag] = dict(min_s14_34=float(r[sel].min()), max_s14_34=float(r[sel].max()))
res['Q0_over_Q_extremes'] = q0
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(q0, indent=1))
