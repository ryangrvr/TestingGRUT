"""SKEPTIC (pt_response_skeptic), INDEPENDENT code.  Does any homogeneous closure reproduce the IR response, and its
orbital-phase dependence?  pt_response tested bare / isotropic Hartree / meanH and concluded "every homogeneous closure
fails in the IR".  Two of those closures are not Goldstone-consistent with the measured orbit:
  * bare  : the measured orbit is held by more than the bare force (tangential mass defect -Sigma, tachyonic);
  * meanH : the full Gaussian (HFB) mean Hessian includes the anomalous pair density <dpsi^2>; like static HFB it is
            gapped: the phase mode i psi_c is not a solution (Hohenberg-Martin / Griffin 'conserving vs gapless').
We measure the Goldstone residual of each closure on the recorded orbit, rho(s) = u'' + H u with u = J qbar (the
global-phase rotation of the coherent mode, the k = 0 Goldstone direction), and build two Goldstone-repaired closures:
  * meanH_Gt : meanH - Delta_t(s) t t^T              (tangential shift only; Delta_t = t.rho/|u|, t = u/|u|)
  * meanH_G2 : meanH + Delta H(s), Delta H u = -rho   (symmetric rank-2, also removes the radial residual)
For each closure, phase s0 in {20, 20.12, 20.24, 20.44, 20.68} and shell 1-4 (every FD4 kappa^2 of the shell, mode
weighted), we integrate d'' + kappa^2 d + H(s) d = 0 (X block, 2x2; RK4 h = 0.0025, own code) for unit momentum kicks
and form normal / conjugate exactly as response.py.  Compared with our measured low-noise traces (measured_traces.py):
L2 relative error of the normal and conjugate channels, delays 0-2 and 0-6, RMS ratio, and the dominant spectral
peaks (Codex's transform: cos^2 taper, trapezoid, omega grid +-20 step 0.02).
Usage: python3 -I closures_phase.py <hessian npz> <meas_traces npz> <out json>
"""
import sys, json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import argrelmax, argrelmin

hz = np.load(sys.argv[1]); me = np.load(sys.argv[2]); OUT = sys.argv[3]
LX = 1e-3
S = hz['s']; m = hz['mean']; q = m[:, :2]; p = m[:, 3:5]
H0 = hz['Hbar'][:, :2, :2].copy(); Hb = hz['Hcoh'][:, :2, :2].copy()
r = np.hypot(q[:, 0], q[:, 1])
J = np.array([[0., -1.], [1., 0.]])
u = q @ J.T                      # J qbar
pp = CubicSpline(S, p, axis=0)(S, 1)   # qbar'' = pbar'
upp = pp @ J.T
def resid(H):
    return upp + np.einsum('tij,tj->ti', H, u)
uu = np.sum(u * u, 1); tt = u / np.sqrt(uu)[:, None]; rr_ = q / r[:, None]
out = dict(method=__doc__)
mx = argrelmax(r)[0]; mn = argrelmin(r)[0]
inner = (S > 20.1) & (S < 26.7)
gold = {}
for name, H in (('bare', Hb), ('meanH', H0)):
    rho = resid(H)
    Dt = np.sum(rho * tt, 1) / np.sqrt(uu); Dr = np.sum(rho * rr_, 1) / np.sqrt(uu)
    gold[name] = dict(tangential_defect_mean=float(Dt[inner].mean()), tangential_defect_apo=float(Dt[mx[(S[mx] > 20.1) & (S[mx] < 26.7)]].mean()),
                      tangential_defect_peri=float(Dt[mn[(S[mn] > 20.1) & (S[mn] < 26.7)]].mean()),
                      radial_defect_mean=float(Dr[inner].mean()), radial_defect_rms=float(np.sqrt(np.mean(Dr[inner] ** 2))))
out['goldstone_residual'] = gold
print('Goldstone residual (tangential mass defect, positive = gapped, negative = tachyonic):')
for k, v in gold.items(): print('  ', k, {a: round(b, 3) for a, b in v.items()})
rho0 = resid(H0)
Dt0 = np.sum(rho0 * tt, 1) / np.sqrt(uu)
H_Gt = H0 - Dt0[:, None, None] * np.einsum('ti,tj->tij', tt, tt)
v = -rho0
H_G2 = H0 + (np.einsum('ti,tj->tij', v, u) + np.einsum('ti,tj->tij', u, v)) / uu[:, None, None] \
          - (np.sum(u * v, 1) / uu ** 2)[:, None, None] * np.einsum('ti,tj->tij', u, u)
assert np.max(np.abs(np.einsum('tij,tj->ti', H_G2, u) + upp)) < 1e-6 * np.max(np.abs(upp))
closures = dict(meanH=H0, meanH_Gt=H_Gt, meanH_G2=H_G2, bare=Hb)
sp = {k: CubicSpline(S, v_.reshape(len(S), 4), axis=0) for k, v_ in closures.items()}

# FD4 kappa^2 of each shell (N48, L12)
N, L = 48, 12.
mm = np.fft.fftfreq(N) * N
kc2 = (2 * np.pi / L) ** 2 * (mm[:, None, None] ** 2 + mm[None, :, None] ** 2 + mm[None, None, :] ** 2)
f = lambda x: 4 * np.sin(x / 2) ** 2 * (1 + np.sin(x / 2) ** 2 / 3) / (L / N) ** 2
x = 2 * np.pi * mm / N
kl = (f(x)[:, None, None] + f(x)[None, :, None] + f(x)[None, None, :]).ravel()
sid = np.floor(np.sqrt(kc2.ravel()) / (2 * np.pi / L) + .5).astype(int)
shk = {n: np.unique(np.round(kl[sid == n], 10), return_counts=True) for n in (1, 2, 3, 4)}
HS = .0025; delays = np.arange(151) * .04

def respond(name, s0, k2s):
    n = int(round(6. / HS)); Hs = sp[name](s0 + np.arange(2 * n + 1) * HS / 2).reshape(-1, 2, 2)
    K = len(k2s); Q = np.zeros((K, 2, 2)); P = np.repeat(np.eye(2)[None], K, 0)
    def acc(Q, Hm): return -(k2s[:, None, None] * Q + np.einsum('ij,kjb->kib', Hm, Q))
    out_ = [Q.copy()]
    for i in range(n):
        a1q, a1p = P, acc(Q, Hs[2 * i]); a2q, a2p = P + .5 * HS * a1p, acc(Q + .5 * HS * a1q, Hs[2 * i + 1])
        a3q, a3p = P + .5 * HS * a2p, acc(Q + .5 * HS * a2q, Hs[2 * i + 1]); a4q, a4p = P + HS * a3p, acc(Q + HS * a3q, Hs[2 * i + 2])
        Q = Q + HS / 6 * (a1q + 2 * a2q + 2 * a3q + a4q); P = P + HS / 6 * (a1p + 2 * a2p + 2 * a3p + a4p)
        if (i + 1) % 16 == 0: out_.append(Q.copy())
    R = np.array(out_)  # [t, k, a, b]
    gN = .5 * (R[..., 0, 0] + R[..., 1, 1] + 1j * (R[..., 1, 0] - R[..., 0, 1])); gA = .5 * (R[..., 0, 0] - R[..., 1, 1] + 1j * (R[..., 1, 0] + R[..., 0, 1]))
    return gN, gA

w = np.linspace(-20, 20, 2001); taper = np.cos(np.pi * delays / 12) ** 2; trap = np.full(151, .04); trap[[0, -1]] *= .5
TR = np.exp(-1j * w[:, None] * delays[None, :])
def peaks(g):
    rho = -2 * (TR @ (g * trap * taper)).imag
    pos, neg = w > 0, w < 0
    return float(w[pos][np.argmax(rho[pos])]), float(w[neg][np.argmax(-rho[neg])])
TN, TA = me['tr_normal'], me['tr_anom']; phases = [float(x) for x in me['phases']]
rows = []
for i, s0 in enumerate(phases):
    for name in closures:
        for n in (1, 2, 3, 4):
            k2s, cnt = shk[n]
            gN, gA = respond(name, s0, k2s)
            mN = (gN * cnt).sum(1) / cnt.sum(); mA = (gA * cnt).sum(1) / cnt.sum()
            row = dict(s0=s0, closure=name, shell=n)
            for lab, sel in (('0_2', delays < 2.0001), ('0_6', delays < 7)):
                row[f'errN_{lab}'] = float(np.linalg.norm(mN[sel] - TN[i, n - 1, sel]) / np.linalg.norm(TN[i, n - 1, sel]))
                row[f'errA_{lab}'] = float(np.linalg.norm(mA[sel] - TA[i, n - 1, sel]) / np.linalg.norm(TA[i, n - 1, sel]))
            row['ratio_model'] = float(np.linalg.norm(mA) / np.linalg.norm(mN)); row['ratio_meas'] = float(np.linalg.norm(TA[i, n - 1]) / np.linalg.norm(TN[i, n - 1]))
            row['peaks_model'] = peaks(mN); row['peaks_meas'] = peaks(TN[i, n - 1])
            rows.append(row)
    print('s0', s0)
    for name in closures:
        rr = [x for x in rows if x['s0'] == s0 and x['closure'] == name]
        print('  %-9s' % name, ' | '.join('sh%d eN %.2f/%.2f eA %.2f r %.2f(%.2f) pk %+.2f %+.2f (%+.2f %+.2f)' % (
            x['shell'], x['errN_0_2'], x['errN_0_6'], x['errA_0_6'], x['ratio_model'], x['ratio_meas'], *x['peaks_model'], *x['peaks_meas']) for x in rr))
out['rows'] = rows
# summary: per closure, median errors over phases, shells 1-3
summ = {}
for name in closures:
    rr = [x for x in rows if x['closure'] == name and x['shell'] <= 3]
    rm = np.array([x['ratio_model'] for x in rr]); rme = np.array([x['ratio_meas'] for x in rr])
    summ[name] = dict(errN_0_2_median=float(np.median([x['errN_0_2'] for x in rr])), errN_0_6_median=float(np.median([x['errN_0_6'] for x in rr])),
                      errN_0_6_range=[float(min(x['errN_0_6'] for x in rr)), float(max(x['errN_0_6'] for x in rr))],
                      errA_0_6_median=float(np.median([x['errA_0_6'] for x in rr])),
                      ratio_model_range=[float(rm.min()), float(rm.max())], ratio_meas_range=[float(rme.min()), float(rme.max())],
                      ratio_phase_correlation=float(np.corrcoef(rm, rme)[0, 1]))
out['summary_shells_1_3'] = summ
for k, v_ in summ.items(): print(k, {a: (np.round(b, 3).tolist() if isinstance(b, list) else round(b, 3)) for a, b in v_.items()})
json.dump(out, open(OUT, 'w'), indent=1)
