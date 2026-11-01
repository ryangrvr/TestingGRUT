"""SKEPTIC (pt_response_skeptic), INDEPENDENT code: UV check of pt_response's 'meanH reproduces poles to <= 2.3% and the
ratio to 5-15%' against Codex's published single-probe responses (shells 5-15, where the probe noise is 1.5-14%).
Same driven homogeneous model as closures_phase.py (X block, RK4 h = .0025, every FD4 kappa^2 with |k| <= 8 in the
shell, mode-weighted), closures meanH and the Goldstone-repaired meanH_Gt; s0 = 20 (mean of seeds 7 and 11) and 20.44.
Usage: python3 -I closures_uv.py <hessian npz> <packet dir> <out json>
"""
import sys, json
import numpy as np
from scipy.interpolate import CubicSpline
hz = np.load(sys.argv[1]); P = sys.argv[2].rstrip('/') + '/'; OUT = sys.argv[3]
S = hz['s']; m = hz['mean']; q = m[:, :2]; p = m[:, 3:5]
H0 = hz['Hbar'][:, :2, :2].copy()
J = np.array([[0., -1.], [1., 0.]]); u = q @ J.T; upp = CubicSpline(S, p, axis=0)(S, 1) @ J.T
uu = np.sum(u * u, 1); tt = u / np.sqrt(uu)[:, None]
rho0 = upp + np.einsum('tij,tj->ti', H0, u); Dt0 = np.sum(rho0 * tt, 1) / np.sqrt(uu)
closures = dict(meanH=H0, meanH_Gt=H0 - Dt0[:, None, None] * np.einsum('ti,tj->tij', tt, tt))
sp = {k: CubicSpline(S, v.reshape(len(S), 4), axis=0) for k, v in closures.items()}
N, L = 48, 12.
mm = np.fft.fftfreq(N) * N
kc2 = ((2 * np.pi / L) ** 2 * (mm[:, None, None] ** 2 + mm[None, :, None] ** 2 + mm[None, None, :] ** 2)).ravel()
f = lambda x: 4 * np.sin(x / 2) ** 2 * (1 + np.sin(x / 2) ** 2 / 3) / (L / N) ** 2
x = 2 * np.pi * mm / N
kl = (f(x)[:, None, None] + f(x)[None, :, None] + f(x)[None, None, :]).ravel()
keep = kc2 <= 64 + 1e-12
sid = np.floor(np.sqrt(kc2) / (2 * np.pi / L) + .5).astype(int)
HS = .0025; delays = np.arange(151) * .04
def respond(name, s0, k2s):
    n = int(round(6. / HS)); Hs = sp[name](s0 + np.arange(2 * n + 1) * HS / 2).reshape(-1, 2, 2)
    K = len(k2s); Q = np.zeros((K, 2, 2)); Pm = np.repeat(np.eye(2)[None], K, 0)
    acc = lambda Q, Hm: -(k2s[:, None, None] * Q + np.einsum('ij,kjb->kib', Hm, Q))
    out = [Q.copy()]
    for i in range(n):
        a1q, a1p = Pm, acc(Q, Hs[2 * i]); a2q, a2p = Pm + .5 * HS * a1p, acc(Q + .5 * HS * a1q, Hs[2 * i + 1])
        a3q, a3p = Pm + .5 * HS * a2p, acc(Q + .5 * HS * a2q, Hs[2 * i + 1]); a4q, a4p = Pm + HS * a3p, acc(Q + HS * a3q, Hs[2 * i + 2])
        Q = Q + HS / 6 * (a1q + 2 * a2q + 2 * a3q + a4q); Pm = Pm + HS / 6 * (a1p + 2 * a2p + 2 * a3p + a4p)
        if (i + 1) % 16 == 0: out.append(Q.copy())
    R = np.array(out)
    return .5 * (R[..., 0, 0] + R[..., 1, 1] + 1j * (R[..., 1, 0] - R[..., 0, 1])), .5 * (R[..., 0, 0] - R[..., 1, 1] + 1j * (R[..., 1, 0] + R[..., 0, 1]))
w = np.linspace(-20, 20, 2001); taper = np.cos(np.pi * delays / 12) ** 2; trap = np.full(151, .04); trap[[0, -1]] *= .5
TR = np.exp(-1j * w[:, None] * delays[None, :])
def peaks(g):
    rho = -2 * (TR @ (g * trap * taper)).imag; pos, neg = w > 0, w < 0
    return float(w[pos][np.argmax(rho[pos])]), float(w[neg][np.argmax(-rho[neg])])
meas = {20.: [np.load(P + 'RESPONSE_s20_seed7.npz'), np.load(P + 'RESPONSE_s20_seed11.npz')], 20.44: [np.load(P + 'RESPONSE_s20.44_seed7.npz')]}
rows = []
for s0, ds in meas.items():
    for n in range(1, 16):
        k2s, cnt = np.unique(np.round(kl[keep & (sid == n)], 10), return_counts=True)
        mN_meas = np.mean([d['normal'][:, n - 1] for d in ds], 0); mA_meas = np.mean([d['anomalous'][:, n - 1] for d in ds], 0)
        row = dict(s0=s0, shell=n, k=float(n * 2 * np.pi / L), ratio_meas=float(np.linalg.norm(mA_meas) / np.linalg.norm(mN_meas)), peaks_meas=peaks(mN_meas))
        for name in closures:
            gN, gA = respond(name, s0, k2s)
            mN = (gN * cnt).sum(1) / cnt.sum(); mA = (gA * cnt).sum(1) / cnt.sum()
            pk = peaks(mN)
            row[name] = dict(ratio=float(np.linalg.norm(mA) / np.linalg.norm(mN)), peaks=pk,
                             pole_err=[abs(pk[0] / row['peaks_meas'][0] - 1), abs(pk[1] / row['peaks_meas'][1] - 1)],
                             errN=float(np.linalg.norm(mN - mN_meas) / np.linalg.norm(mN_meas)))
        rows.append(row)
        print('s0 %.2f sh%2d k %.2f meas r %.3f pk %+.2f %+.2f | meanH r %.3f pk %+.2f %+.2f eN %.2f | Gt r %.3f pk %+.2f %+.2f eN %.2f' % (
            s0, n, row['k'], row['ratio_meas'], *row['peaks_meas'], row['meanH']['ratio'], *row['meanH']['peaks'], row['meanH']['errN'],
            row['meanH_Gt']['ratio'], *row['meanH_Gt']['peaks'], row['meanH_Gt']['errN']))
summ = {}
for name in closures:
    uv = [r for r in rows if r['shell'] >= 5]
    summ[name] = dict(uv_pole_err_max=float(max(max(r[name]['pole_err']) for r in uv)),
                      uv_pole_err_max_k_ge_4p7=float(max(max(r[name]['pole_err']) for r in uv if r['k'] >= 4.7)),
                      uv_ratio_rel_err_range=[float(min(r[name]['ratio'] / r['ratio_meas'] - 1 for r in uv)), float(max(r[name]['ratio'] / r['ratio_meas'] - 1 for r in uv))],
                      uv_ratio_rel_err_median_abs=float(np.median([abs(r[name]['ratio'] / r['ratio_meas'] - 1) for r in uv])),
                      ir_pole_err_range=[float(min(min(r[name]['pole_err']) for r in rows if r['shell'] <= 4)), float(max(max(r[name]['pole_err']) for r in rows if r['shell'] <= 4))],
                      uv_errN_range=[float(min(r[name]['errN'] for r in uv)), float(max(r[name]['errN'] for r in uv))])
print(json.dumps(summ, indent=1))
json.dump(dict(method=__doc__, rows=rows, summary=summ), open(OUT, 'w'), indent=1)
