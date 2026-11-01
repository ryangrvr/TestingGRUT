"""INDEPENDENT code (Claude, pt_response): coherent-mode counting in two bases over the SAME interval.

From our recorded zero-mode history (record_mean_hessian.py, Delta s = 0.01, s 20-26.8):
  * lab-frequency negative action of the coherent mode (Codex's statistic, Hann window over whole radial periods),
  * instantaneous-orbit radial action J_r = (1/2pi) * integral p_r dr over each radial period (= the k=0 radial-branch
    quasiparticle action about the circular condensate of the same charge; what slow matching preserves),
both per measured net charge Q.  Also the same two numbers from Codex's published base48 trace (Delta s=0.04) over
s 14-24 and 24-34, for comparison with Codex's coherent window values 0.6075 / 0.4760.
Usage: python3 -I coherent_counting.py <hessian npz> <packet dir> <out json>
"""
import sys, json
import numpy as np
from scipy.signal import argrelmax

hz = np.load(sys.argv[1]); P = sys.argv[2].rstrip('/') + '/'; OUT = sys.argv[3]

def analyse(s, q1, q2, p1, p2, Q, lo, hi):
    sel = (s >= lo - 1e-9) & (s <= hi + 1e-9); s = s[sel]; q1, q2, p1, p2 = q1[sel], q2[sel], p1[sel], p2[sel]
    psi = q1 + 1j * q2; ds = s[1] - s[0]; n = len(s)
    h = np.hanning(n); X = np.fft.fft(h * psi); Pw = abs(X) ** 2 / (n * np.sum(h * h)); w = 2 * np.pi * np.fft.fftfreq(n, ds)
    Am = float(np.sum((-w * Pw)[w < 0]) / Q); Ap = float(np.sum((w * Pw)[w > 0]) / Q)
    r = np.sqrt(q1 ** 2 + q2 ** 2); pr = (q1 * p1 + q2 * p2) / r; qc = q1 * p2 - q2 * p1
    mx = argrelmax(r)[0]; J = []; qq = []
    for a, b in zip(mx[:-1], mx[1:]):
        J.append(np.trapezoid(pr[a:b + 1] ** 2, s[a:b + 1]) / (2 * np.pi)); qq.append(np.trapezoid(qc[a:b + 1], s[a:b + 1]) / (s[b] - s[a]))
    J = np.array(J); qq = np.array(qq)
    # Hann-weighted mean of per-period J_r (weights at period mid-times) for a like-for-like comparison
    mids = np.array([(s[a] + s[b]) / 2 for a, b in zip(mx[:-1], mx[1:])])
    wts = np.sin(np.pi * (mids - s[0]) / (s[-1] - s[0])) ** 2
    return dict(window=[float(s[0]), float(s[-1])], lab_A_minus_per_Q=Am, lab_A_plus_per_Q=Ap, lab_charge_check=Ap - Am,
                coherent_charge_over_Q=float(np.mean(qc) / Q),
                Jr_per_Q_periods=(J / Q).tolist(), Jr_per_Q_mean=float(np.mean(J) / Q), Jr_per_Q_hannweighted=float(np.sum(wts * J) / np.sum(wts) / Q),
                Jr_per_coherent_charge_mean=float(np.mean(J / qq)), radial_periods=len(J))

res = {}
m = hz['mean']; Q = float(hz['Q'][0])
res['ours_0.01_s20.09_26.40'] = analyse(hz['s'], m[:, 0], m[:, 1], m[:, 3], m[:, 4], Q, 20.09, 26.40)
F = np.load(P + 'FIELD_TRACES.npz'); names = list(F['names']); v = F['base48_values']; s = F['base48_s']
Qp = float(v[0, names.index('nQ')])
ix = [names.index(x) for x in ['q0', 'q1', 'p0', 'p1']]
for lo, hi in [(14., 24.), (24., 34.), (0., 6.)]:
    res[f'codex_trace_{lo:g}_{hi:g}'] = analyse(s, *[v[:, i] for i in ix], Qp, lo, hi)
W = {k: json.load(open(P + f'WINDOW_base48_{k}.json'))['coherent_negative_frequency_action'] for k in ['14_24', '24_34', '0_6']}
res['codex_window_coherent_A_minus'] = W
for k, r in res.items():
    if 'lab_A_minus_per_Q' in r:
        print(k, 'lab A- %.4f' % r['lab_A_minus_per_Q'], 'Jr(mean) %.4f' % r['Jr_per_Q_mean'], 'Jr(Hann) %.4f' % r['Jr_per_Q_hannweighted'],
              'ratio %.4f' % (r['Jr_per_Q_hannweighted'] / r['lab_A_minus_per_Q']), 'periods', r['radial_periods'], 'charge chk %.4f' % r['lab_charge_check'], 'Qc/Q %.4f' % r['coherent_charge_over_Q'])
print('codex windows', W)
json.dump(res, open(OUT, 'w'), indent=1)
