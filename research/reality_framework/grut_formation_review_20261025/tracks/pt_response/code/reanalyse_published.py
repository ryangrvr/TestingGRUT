"""INDEPENDENT analysis code (Claude, track pt_response) on Codex's PUBLISHED arrays (5f5513b).

Reads RESPONSE_s*.npz / .json and RESPONSE_FREE_CHECK.json from the extracted packet and recomputes,
with our own code (no Codex import):
  * anomalous/normal RMS ratio per shell (Codex: anomalous_over_normal_RMS),
  * the canonical f-sum Z+ + Z- with our own transform (and the exact t->0 slope),
  * the free control: our own Yoshida-4 integration of a single FD4 oscillator with the same step rule,
  * decomposition algebra: delta psi = G_N p + G_A p*  (checked symbolically on random matrices).
Usage: python3 -I reanalyse_published.py <packet dir> <out json>
"""
import sys, json
import numpy as np

P = sys.argv[1].rstrip('/') + '/'
OUT = sys.argv[2]
cases = [('s20_seed7', 20.0, 7), ('s20.44_seed7', 20.44, 7), ('s20_seed11', 20.0, 11)]
res = dict(source='Codex 5f5513b published RESPONSE_*.npz/.json; our own analysis code', cases={})

# 1. decomposition algebra: for real 2x2 R (chi_a response to pi_b kick) and complex kick p = p1 + i p2,
#    delta psi = (R00 p1 + R01 p2) + i (R10 p1 + R11 p2) == N p + A conj(p)
rng = np.random.default_rng(0)
err = 0.
for _ in range(100):
    R = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))  # complex allowed (shell traces)
    p1, p2 = rng.normal(size=2)
    p = p1 + 1j * p2
    dpsi = (R[0, 0] * p1 + R[0, 1] * p2) + 1j * (R[1, 0] * p1 + R[1, 1] * p2)
    N = .5 * (R[0, 0] + R[1, 1] + 1j * (R[1, 0] - R[0, 1]))
    A = .5 * (R[0, 0] - R[1, 1] + 1j * (R[1, 0] + R[0, 1]))
    err = max(err, abs(dpsi - (N * p + A * np.conj(p))))
res['decomposition_identity_max_error'] = float(err)
res['decomposition_note'] = ('With real kicks p1,p2 (as in response.py: each tangent kicks one real momentum), '
                             'normal = coefficient of p, anomalous = coefficient of conj(p) for p = d(psi)/ds kick, '
                             'psi = chi1 + i chi2. Holds for complex shell-trace entries too.')

for tag, s0, seed in cases:
    d = np.load(P + f'RESPONSE_{tag}.npz', allow_pickle=False)
    j = json.load(open(P + f'RESPONSE_{tag}.json'))
    t = d['delay']; Nn = d['normal']; An = d['anomalous']; w = d['omega']; rho_pub = d['rho']
    rows = []
    for k, sh in enumerate(d['shell']):
        ratio = np.linalg.norm(An[:, k]) / np.linalg.norm(Nn[:, k])
        # our own transform: same taper & trapezoid weights as stated in response.py docstring of the rows,
        # but evaluated with an independent implementation (explicit sine transform for real/complex parts)
        T = t[-1]; dt = t[1] - t[0]
        taper = np.cos(np.pi * t / (2 * T)) ** 2
        trap = np.full(len(t), dt); trap[[0, -1]] *= .5
        f = Nn[:, k] * taper * trap
        # G(w) = sum_t e^{-i w t} f(t); rho = -2 Im G
        Gw = np.array([np.sum(np.exp(-1j * ww * t) * f) for ww in w])
        rho = -2 * Gw.imag
        dw = w[1] - w[0]
        Zp = np.sum((w * rho)[w > 0]) * dw / (2 * np.pi); Zm = np.sum((w * rho)[w < 0]) * dw / (2 * np.pi)
        # alternative transform-free canonical check: exact slope G'(0)=1 (Taylor through t=0,0.04,0.08)
        slope3 = (-3 * 0 + 4 * Nn[1, k] - Nn[2, k]) / (2 * dt)  # one-sided 2nd-order (G(0)=0)
        pub = j['rows'][k]
        rows.append(dict(shell=int(sh), k_center=pub['k_center'], ratio_ours=float(ratio),
                         ratio_codex=pub['anomalous_over_normal_RMS'], Zsum_ours=float(Zp + Zm),
                         Zsum_codex=pub['Z_sum'], Zp_ours=float(Zp), Zm_ours=float(Zm),
                         rho_max_abs_diff=float(np.max(abs(rho - rho_pub[:, k]))),
                         slope_t0_second_order=float(slope3.real),
                         imag_normal_rms_over_real=float(np.linalg.norm(Nn[:, k].imag) / np.linalg.norm(Nn[:, k].real)),
                         normal_late_over_early_envelope=float(np.max(abs(Nn[100:, k])) / np.max(abs(Nn[:50, k]))),
                         phi_to_X_rms_over_normal=float(np.linalg.norm(d['phi_to_X'][:, k, :]) / np.linalg.norm(Nn[:, k])),
                         X_to_phi_rms_over_normal=float(np.linalg.norm(d['X_to_phi'][:, k, :]) / np.linalg.norm(Nn[:, k]))))
    res['cases'][tag] = dict(s0=s0, seed=seed, rows=rows,
                             max_ratio_abs_diff=max(abs(r['ratio_ours'] - r['ratio_codex']) for r in rows),
                             max_Zsum_abs_diff=max(abs(r['Zsum_ours'] - r['Zsum_codex']) for r in rows),
                             max_Zsum_dev_from_1=max(abs(r['Zsum_ours'] - 1) for r in rows))

ir = [r['ratio_ours'] for c in res['cases'].values() for r in c['rows'][:3]]
res['IR_shells_1to3_ratio_range_ours'] = [min(ir), max(ir)]
res['max_Zsum_deviation_all'] = max(c['max_Zsum_dev_from_1'] for c in res['cases'].values())
# probe-noise estimate: two probes at the same phase s0=20
a = np.load(P + 'RESPONSE_s20_seed7.npz'); b = np.load(P + 'RESPONSE_s20_seed11.npz')
res['probe_difference_s20'] = [dict(shell=int(sh),
    normal_rel_diff=float(np.linalg.norm(a['normal'][:, k] - b['normal'][:, k]) / np.linalg.norm(a['normal'][:, k])),
    anomalous_rel_diff=float(np.linalg.norm(a['anomalous'][:, k] - b['anomalous'][:, k]) / np.linalg.norm(a['anomalous'][:, k])))
    for k, sh in enumerate(a['shell'])]

# 2. free control, own code: FD4 symbol and Yoshida-4 for x'' = -w^2 x, x(0)=0, x'(0)=1, h<=0.004, samples every 0.04
N, L = 8, 5.
u = 2 * np.pi / N; dx = L / N
w2 = 4 * np.sin(u / 2) ** 2 * (1 + np.sin(u / 2) ** 2 / 3) / dx ** 2   # one axis, mode n=1
om = np.sqrt(w2)
w1 = 1 / (2 - 2 ** (1 / 3)); w0 = 1 - 2 * w1
x, v = 0., 1.; times = np.arange(151) * .04; err = 0.; prev = 0.
for tt in times:
    dur = tt - prev
    if dur > 0:
        n = int(np.ceil(dur / .004 - 1e-12)); h = dur / n
        for _ in range(n):
            for c, dd in [(w1 / 2, w1), ((w1 + w0) / 2, w0), ((w0 + w1) / 2, w1)]:
                v -= c * h * w2 * x; x += dd * h * v
            v -= w1 / 2 * h * w2 * x
    prev = tt
    err = max(err, abs(x - np.sin(om * tt) / om))
fc = json.load(open(P + 'RESPONSE_FREE_CHECK.json'))
res['free_control'] = dict(omega_ours=float(om), omega_codex=fc['exact_frequency'], max_err_ours=float(err),
                           max_err_codex=fc['maximum_response_absolute_error'],
                           codex_cross=fc['maximum_cross_response'], codex_anomalous=fc['maximum_anomalous_response'])
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != 'cases'}, indent=1)[:3000])
for tag, c in res['cases'].items():
    print(tag, 'max|ratio diff|', c['max_ratio_abs_diff'], 'max|Zsum diff|', c['max_Zsum_abs_diff'], 'max|Z-1|', c['max_Zsum_dev_from_1'])
    for r in c['rows'][:6]:
        print('  ', r['shell'], round(r['ratio_ours'], 4), 'Zsum', round(r['Zsum_ours'], 5), 'slope', round(r['slope_t0_second_order'], 5),
              'ImN/ReN', round(r['imag_normal_rms_over_real'], 3), 'late/early', round(r['normal_late_over_early_envelope'], 3),
              'phiX', round(r['phi_to_X_rms_over_normal'], 4))
