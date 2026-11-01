"""INDEPENDENT code (Claude, pt_response): Floquet-Bogoliubov analysis of the eccentric medium.

Idealisation: one measured radial period [sa, sb] of the background (s ~ 20-23) is repeated periodically in the frame
co-rotating with the coherent mode (the orbit precesses by dth = theta(sb) - theta(sa) per period).
For each FD4 kappa^2 (shell means), the lab 4x4 system d'' + kappa^2 d + H(s) d = 0 (X block; phi mixing neglected,
<= 4% of the normal response except one noisy shell) is integrated over one period (RK4, h=0.0025), giving
Phi(s, sa).  The co-rotating monodromy M = R4(-th_b) Phi(sb,sa) R4(th_a) gives Floquet multipliers mu (|mu|=1:
stable quasiparticle; |mu|>1: parametric instability, rate ln|mu|/T).  For each stable positive-norm Floquet mode
(Wronskian W = i sum(q* p - p* q) > 0, which equals the classical action of the real k,-k mode pair), its lab
trajectory over NP periods is built exactly from mu^n R4(n dth) Phi, and the lab-frequency spectra of
delta psi_k = q1 + i q2 and delta psi_{-k} = conj(q1) + i conj(q2) give A-/W and A+/W per unit quasiparticle action
(Codex's Hann |omega|-weighted statistic per Floquet quasiparticle).
Closures: 'meanH' (spatially averaged Hessian, Gaussian/HFB) and 'bare' (coherent mode only).
Usage: python3 -I floquet.py <hessian npz> <packet dir> <out json>      (real background)
       python3 -I floquet.py --test                                    (circular & free controls)
"""
import sys, json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import argrelmax

H = 0.0025

def R4(a):
    c, s = np.cos(a), np.sin(a); R = np.array([[c, -s], [s, c]]); Z = np.zeros((2, 2))
    return np.block([[R, Z], [Z, R]])

def one_period(Hfun, sa, sb, k2):
    n = int(round((sb - sa) / H)); h = (sb - sa) / n
    Hm = Hfun(sa + np.arange(2 * n + 1) * h / 2)
    Y = np.eye(4); traj = [Y.copy()]
    def f(Y, Hk):
        A = np.zeros((4, 4)); A[:2, 2:] = np.eye(2); A[2:, :2] = -(k2 * np.eye(2) + Hk); return A @ Y
    for i in range(n):
        a1 = f(Y, Hm[2 * i]); a2 = f(Y + .5 * h * a1, Hm[2 * i + 1]); a3 = f(Y + .5 * h * a2, Hm[2 * i + 1]); a4 = f(Y + h * a3, Hm[2 * i + 2])
        Y = Y + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4); traj.append(Y.copy())
    return np.array(traj), h

def floquet_modes(Hfun, sa, sb, tha, thb, k2, NP=64):
    T = sb - sa; dth = thb - tha
    traj, h = one_period(Hfun, sa, sb, k2)
    M = R4(-thb) @ traj[-1] @ R4(tha)
    mu, V = np.linalg.eig(M)
    growth = float(np.max(np.log(abs(mu))) / T)
    row = dict(k2=float(k2), max_growth_rate=growth, multipliers_abs=sorted(abs(mu).tolist()), modes=[])
    if growth > 1e-6: return row
    sub = traj[::4][:-1]
    for j in range(4):
        y0 = R4(tha) @ V[:, j]; q = y0[:2]; p = y0[2:]
        W = float(np.real(1j * (np.vdot(q, p) - np.vdot(p, q))))
        if W <= 0: continue
        Ys = np.einsum('tij,j->ti', sub, y0)
        Yl = np.concatenate([(mu[j] ** nper) * (R4(nper * dth) @ Ys.T).T for nper in range(NP)])
        tt = np.arange(len(Yl)) * 4 * h
        psik = Yl[:, 0] + 1j * Yl[:, 1]; psimk = np.conj(Yl[:, 0]) + 1j * np.conj(Yl[:, 1])
        nn = len(tt); hw = np.hanning(nn); w = 2 * np.pi * np.fft.fftfreq(nn, tt[1] - tt[0])
        Ptot = (abs(np.fft.fft(hw * psik)) ** 2 + abs(np.fft.fft(hw * psimk)) ** 2) / (nn * np.sum(hw * hw))
        Am = float(np.sum((-w * Ptot)[w < 0])); Ap = float(np.sum((w * Ptot)[w > 0]))
        ipk = np.argsort(Ptot)[::-1][:3]
        row['modes'].append(dict(W=W, quasi_phase=float(np.angle(mu[j])), A_minus_over_W=Am / W, A_plus_over_W=Ap / W,
                                 charge_over_W=(Ap - Am) / W, top_lab_lines=[(float(w[i]), float(Ptot[i] / Ptot.max())) for i in ipk]))
    row['modes'].sort(key=lambda d: -d['A_minus_over_W'])
    return row

def test():
    out = {}
    # free massive: H = m^2 I -> pure particle (A-/W = 0) and antiparticle (A-/W = 1)
    m2 = 16.; Hf = lambda s: np.broadcast_to(m2 * np.eye(2), (len(s), 2, 2))
    out['free'] = floquet_modes(Hf, 0., .9, 0., 0., 2.7)
    # circular HFB condensate (C1 parameters): H_lab = R(w0 s) diag(Hrr, Htt) R(-w0 s); compare with circular r-
    Hrr, Htt, w0 = 24.21490972, 16.275851, np.sqrt(16.275851)
    def Hc(s):
        c, sn = np.cos(w0 * s), np.sin(w0 * s)
        R = np.array([[c, -sn], [sn, c]]).transpose(2, 0, 1); D = np.diag([Hrr, Htt])
        return R @ D @ np.transpose(R, (0, 2, 1))
    out['circular_C1'] = [floquet_modes(Hc, 0., .9, 0., w0 * .9, k2) for k2 in (0.457, 2.70, 57.79)]
    # eccentric harmonic background without fluctuations: psi_c = a e^{iws} + b e^{-iws}, H = lam|psi|^2 I + 2 lam psi psi^T
    lam, w, a, b = 1e-3, 3.75, 62.4, 37.0
    def He(s):
        z = a * np.exp(1j * w * s) + b * np.exp(-1j * w * s); v = np.array([z.real, z.imag]).T
        return lam * (abs(z) ** 2)[:, None, None] * np.eye(2) + 2 * lam * v[:, :, None] * v[:, None, :] + 10. * np.eye(2)
    out['static_anomalous_demo'] = [floquet_modes(He, 0., 2 * np.pi / w, 0., 0., k2) for k2 in (57.79,)]
    print(json.dumps(out, indent=1)[:6000])
    return out

def main():
    hz = np.load(sys.argv[1]); P = sys.argv[2].rstrip('/') + '/'; OUT = sys.argv[3]
    S = hz['s']; m = hz['mean']
    th = np.unwrap(np.arctan2(m[:, 1], m[:, 0])); r = np.hypot(m[:, 0], m[:, 1])
    mx = argrelmax(r)[0]
    shells = json.load(open(P + 'WINDOW_base48_14_24.json'))['shells']
    k2s = [sh['kfd2_mean'] for sh in shells if sh['shell'] > 0]
    res = dict(method=__doc__, periods={})
    for (ia, ib) in [(mx[0], mx[1]), (mx[2], mx[3])]:
        sa, sb = S[ia], S[ib]; key = f'{sa:.2f}_{sb:.2f}'
        res['periods'][key] = dict(T=float(sb - sa), precession=float(th[ib] - th[ia]), closures={})
        for closure, arr in [('meanH', hz['Hbar']), ('bare', hz['Hcoh'])]:
            sp = CubicSpline(S, arr.reshape(len(S), 9), axis=0)
            Hfun = lambda s, sp=sp: sp(s).reshape(-1, 3, 3)[:, :2, :2]
            rows = [floquet_modes(Hfun, sa, sb, th[ia], th[ib], k2) for k2 in k2s]
            res['periods'][key]['closures'][closure] = rows
            print(key, closure)
            for r_ in rows:
                print('   k2 %.2f growth %.3f' % (r_['k2'], r_['max_growth_rate']), ' | '.join('A-/W %.3f A+/W %.3f' % (d['A_minus_over_W'], d['A_plus_over_W']) for d in r_['modes']))
    inv = {}
    for key, pr in res['periods'].items():
        rows = pr['closures']['meanH']; Jrad = Jph = Am = 0.; ok = True; per = []; conds = []
        for sh, r_ in zip([x for x in shells if x['shell'] > 0], rows):
            y = np.array([sh['negative_frequency_action_per_net'], sh['positive_frequency_action_per_net']])
            if len(r_['modes']) != 2: ok = False; per.append(dict(shell=sh['shell'], unstable=True, A_minus=float(y[0]))); continue
            Mx = np.array([[r_['modes'][0]['A_minus_over_W'], r_['modes'][1]['A_minus_over_W']], [r_['modes'][0]['A_plus_over_W'], r_['modes'][1]['A_plus_over_W']]])
            J = np.linalg.solve(Mx, y); Jrad += J[0]; Jph += J[1]; Am += y[0]; conds.append(float(np.linalg.cond(Mx)))
            per.append(dict(shell=sh['shell'], J_antiparticle_like=float(J[0]), J_particle_like=float(J[1]), A_minus=float(y[0]), cond=conds[-1]))
        inv[key] = dict(all_stable=ok, waves_A_minus_lab_stable_shells=Am, waves_J_antiparticle_like=Jrad, waves_J_particle_like=Jph,
                        max_condition_number=max(conds), per_shell=per)
        print('inversion', key, 'all stable', ok, 'lab A- %.4f  J_anti-like %.4f  J_part-like %.4f  max cond %.1f' % (Am, Jrad, Jph, max(conds)))
    res['counting_inversion_window_14_24'] = inv
    json.dump(res, open(OUT, 'w'), indent=1)

if __name__ == '__main__':
    if sys.argv[1] == '--test': test()
    else: main()
