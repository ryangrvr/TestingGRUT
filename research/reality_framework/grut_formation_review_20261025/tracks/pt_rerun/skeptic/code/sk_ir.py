"""Skeptic of pt_rerun, audit (a): continuum handoff RMS and the obstruction structure by a different method (own code).

pt_rerun recomputed transfer matrices at new k with Codex's integrator. Here no integrator is used at all:
* the light transfer rows T(k) (4x12: chi1, chi2, pi1, pi2 from the 12 initial phase-space variables) are recovered at
  every published L24 shell from Codex's factor F = T S(k) by T = F S^-1 (S is the stated initial-state factor);
* A = T_q U and B = T_p U (smooth, analytic in k^2) are interpolated with cubic splines in k^2 over the 780 L24 shells;
  the interpolation error is measured by leave-every-other-shell-out;
* the singular state factor is applied exactly: per mode Var_x(k) = sum_a A_xa^2/(2 W_a) + B_xa^2 W_a/2, with
  W_a = sqrt(k^2+lam_a) (stable) or sqrt(k^2+sigma^2) (unstable), so the 1/W_a squeeze at the obstructions is exact;
* the continuum integral (1/2pi^2) int_0^8 k^2 dk Var(k) is done on panels with breaks at k0a, k0b and 1.8 and the
  substitution k = k0 + u^2 on the stable side of each obstruction.
Usage: python3 -B sk_ir.py <PREPARATION_L24.npz> <PREPARATION_L12.npz> <out.json>
"""
import sys, json
import numpy as np
from scipy.interpolate import CubicSpline

p24, p12, out = sys.argv[1:4]
d = np.load(p24); d12 = np.load(p12)
lam = d['initial_mass_eigenvalues']; U = d['initial_mass_eigenvectors']; ks = d['ks']; F1 = d['factor_width1']; F3 = d['factor_width3']
r0 = np.linalg.norm(d['light_mean_q'][:2]); p0 = np.linalg.norm(d['light_mean_p'][:2])
k0 = np.sort(np.sqrt(-lam[lam < 0])); k0a, k0b = k0
res = dict(eigenvalues=lam.tolist(), obstruction_k=k0.tolist())
# overlap of the unstable eigen-directions with the charged X fields (rows 0,1 of the 6 real fields)
neg = np.flatnonzero(lam < 0)
res['unstable_direction_X_content'] = [float(U[0, a] ** 2 + U[1, a] ** 2) for a in neg]


def W(k, sigma):
    w2 = k * k + lam
    return np.sqrt(np.where(w2 > 0, w2, k * k + sigma * sigma))


def S_of(k, sigma):
    w = W(k, sigma); S = np.zeros((12, 12)); S[:6, :6] = U / np.sqrt(2 * w)[None, :]; S[6:, 6:] = U * np.sqrt(w / 2)[None, :]
    return S


# recover A, B on the L24 shells
A = np.zeros((len(ks), 4, 6)); B = np.zeros((len(ks), 4, 6))
for i, k in enumerate(ks):
    T = F1[i][[0, 1, 3, 4], :12] @ np.linalg.inv(S_of(k, 1.))
    A[i] = T[:, :6] @ U; B[i] = T[:, 6:] @ U
# consistency: rebuild width-3 factor from the width-1-recovered T
err3 = 0.
for i, k in enumerate(ks):
    T = np.c_[A[i] @ U.T, B[i] @ U.T]
    err3 = max(err3, np.max(abs(T @ S_of(k, 3.) - F3[i][[0, 1, 3, 4], :12])) / np.max(abs(F3[i])))
res['width3_rebuilt_from_width1_transfer_max_rel'] = float(err3)

x = ks ** 2
splA = CubicSpline(x, A.reshape(len(ks), -1), axis=0); splB = CubicSpline(x, B.reshape(len(ks), -1), axis=0)
# leave-every-other-out interpolation error, relative to the local max entry
ev, od = np.arange(0, len(ks), 2), np.arange(1, len(ks) - 1, 2)
sA = CubicSpline(x[ev], A[ev].reshape(len(ev), -1), axis=0); sB = CubicSpline(x[ev], B[ev].reshape(len(ev), -1), axis=0)
eA = np.max(abs(sA(x[od]) - A[od].reshape(len(od), -1)), axis=1) / np.max(abs(A[od].reshape(len(od), -1)), axis=1)
eB = np.max(abs(sB(x[od]) - B[od].reshape(len(od), -1)), axis=1) / np.max(abs(B[od].reshape(len(od), -1)), axis=1)
res['spline_leave_half_out_max_rel'] = dict(A=float(eA.max()), B=float(eB.max()), A_k_lt_3=float(eA[ks[od] < 3].max()), B_k_lt_3=float(eB[ks[od] < 3].max()))
# check against the L12 shells (all of them are L24 shells, 4n): sanity of the transfer recovery
d12k = d12['ks']; F12 = d12['factor_width1']; worst = 0.
for i, k in enumerate(d12k):
    a = splA(k * k).reshape(4, 6); b = splB(k * k).reshape(4, 6); T = np.c_[a @ U.T, b @ U.T]
    worst = max(worst, np.max(abs(T @ S_of(k, 1.) - F12[i][[0, 1, 3, 4], :12])) / np.max(abs(F12[i])))
res['spline_vs_L12_factors_max_rel'] = float(worst)


def var_mode(k, sigma, rows):
    a = splA(k * k).reshape(-1, 4, 6)[:, rows]; b = splB(k * k).reshape(-1, 4, 6)[:, rows]; w = np.array([W(kk, sigma) for kk in np.atleast_1d(k)])
    return np.sum(a ** 2 / (2 * w[:, None, :]) + b ** 2 * w[:, None, :] / 2, axis=(1, 2))


def gl(a, b, n):
    x_, w_ = np.polynomial.legendre.leggauss(n); return (a + b) / 2 + (b - a) / 2 * x_, (b - a) / 2 * w_


def panels(n):
    P = []
    P.append(('0_k0a', *gl(1e-12, k0a, n)))
    u, w = gl(0, np.sqrt(k0b - k0a), n); P.append(('k0a_k0b', k0a + u * u, 2 * u * w))
    u, w = gl(0, np.sqrt(1.8 - k0b), n); P.append(('k0b_1.8', k0b + u * u, 2 * u * w))
    P.append(('1.8_3', *gl(1.8, 3., n))); P.append(('3_8', *gl(3., 8., 2 * n)))
    return P


cont = {}
for sigma in (1., 3.):
    for n in (32, 64, 128):
        tot = {'q': 0., 'p': 0.}; parts = {}
        for name, kk, ww in panels(n):
            vq = var_mode(kk, sigma, [0, 1]); vp = var_mode(kk, sigma, [2, 3])
            iq = np.sum(ww * kk ** 2 * vq) / (2 * np.pi ** 2); ip = np.sum(ww * kk ** 2 * vp) / (2 * np.pi ** 2)
            parts[name] = dict(q=iq, p=ip); tot['q'] += iq; tot['p'] += ip
        cont[f'sigma{sigma:g}_n{n}'] = dict(q_RMS_over_mean=float(np.sqrt(tot['q']) / r0), p_RMS_over_mean=float(np.sqrt(tot['p']) / p0),
                                            q_share={k_: float(v['q'] / tot['q']) for k_, v in parts.items()},
                                            p_share={k_: float(v['p'] / tot['p']) for k_, v in parts.items()},
                                            q_share_k_lt_1p8=float(sum(v['q'] for k_, v in parts.items() if k_ in ('0_k0a', 'k0a_k0b', 'k0b_1.8')) / tot['q']))
res['continuum'] = cont

# per-mode variance structure at the obstructions: left limit finite; right side ~ c (k-k0)^-1/2
struct = {}
for name, kz in (('k0a', k0a), ('k0b', k0b)):
    left = float(var_mode(np.array([kz - 1e-9]), 1., [0, 1])[0])
    eps = np.array([1e-6, 1e-4, 1e-3, 1e-2, 3e-2])
    right = var_mode(kz + eps, 1., [0, 1])
    struct[name] = dict(left_limit_per_mode=left, right_per_mode=right.tolist(), eps=eps.tolist(),
                        right_times_sqrt_eps=(right * np.sqrt(eps)).tolist(), ratio_right_over_left=(right / left).tolist())
# what the boxes see: ratio across each obstruction between the nearest box shells
for L, dd in (('L12', d12), ('L24', d)):
    kk = dd['ks']; FF = dd['factor_width1']; v = np.einsum('kia,kia->k', FF[:, [0, 1], :], FF[:, [0, 1], :])
    for name, kz in (('k0a', k0a), ('k0b', k0b)):
        lo = np.flatnonzero(kk < kz)[-1]; hi = np.flatnonzero(kk > kz)[0]
        struct[f'{L}_{name}_nearest_shells'] = dict(k_below=float(kk[lo]), k_above=float(kk[hi]), ratio=float(v[hi] / v[lo]))
res['obstruction_structure'] = struct
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
