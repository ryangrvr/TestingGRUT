"""Skeptic of pt_rerun: how much does the tau=20 handoff depend on the near-obstruction state rule? (own code, Claude)

Same transfer recovery and k^2-spline as sk_ir.py (no integrator). Continuum RMS for alternative pure Gaussian rules on
k < 1.8 only (the old seeded band k >= 1.8 is kept exactly as published):
  codex_s1, codex_s3, codex_s100, codex_s300 : W = sqrt(k^2+lam) if positive else sqrt(k^2+sigma^2)
  clamp1    : W = sqrt(max(k^2+lam, 1))  (no squeeze below W = 1; removes the 1/sqrt divergence)
  absval    : W = sqrt(|k^2+lam|)        (inverted-oscillator magnitude; divergent on both sides, integrable)
Also counts the Fourier modes (lattice vectors) below each obstruction in L12 and L24.
Usage: python3 -B sk_ir_alt.py <PREPARATION_L24.npz> <out.json>
"""
import sys, json
import numpy as np
from scipy.interpolate import CubicSpline

p24, out = sys.argv[1:3]
d = np.load(p24)
lam = d['initial_mass_eigenvalues']; U = d['initial_mass_eigenvectors']; ks = d['ks']; F1 = d['factor_width1']
r0 = np.linalg.norm(d['light_mean_q'][:2]); p0 = np.linalg.norm(d['light_mean_p'][:2])
k0a, k0b = np.sort(np.sqrt(-lam[lam < 0]))


def W_codex(k, sigma):
    w2 = k * k + lam; return np.sqrt(np.where(w2 > 0, w2, k * k + sigma * sigma))


def S_of(w):
    S = np.zeros((12, 12)); S[:6, :6] = U / np.sqrt(2 * w)[None, :]; S[6:, 6:] = U * np.sqrt(w / 2)[None, :]; return S


A = np.zeros((len(ks), 4, 6)); B = np.zeros((len(ks), 4, 6))
for i, k in enumerate(ks):
    T = F1[i][[0, 1, 3, 4], :12] @ np.linalg.inv(S_of(W_codex(k, 1.))); A[i] = T[:, :6] @ U; B[i] = T[:, 6:] @ U
splA = CubicSpline(ks ** 2, A.reshape(len(ks), -1), axis=0); splB = CubicSpline(ks ** 2, B.reshape(len(ks), -1), axis=0)

rules = {
    'codex_s1': lambda k: W_codex(k, 1.), 'codex_s3': lambda k: W_codex(k, 3.),
    'codex_s100': lambda k: W_codex(k, 100.), 'codex_s300': lambda k: W_codex(k, 300.),
    'clamp1': lambda k: (np.sqrt(np.maximum(k * k + lam, 1.)) if k < 1.8 else W_codex(k, 1.)),
    'absval': lambda k: (np.sqrt(np.abs(k * k + lam)) if k < 1.8 else W_codex(k, 1.)),
}


def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n); return (a + b) / 2 + (b - a) / 2 * x, (b - a) / 2 * w


def nodes(n):
    P = []
    # square-root substitutions on both sides of each obstruction (covers absval's two-sided divergence)
    for a, b in ((0., k0a), (k0a, k0b), (k0b, 1.8)):
        m = (a + b) / 2
        if a > 0:
            u, w = gl(0, np.sqrt(m - a), n); P.append((a + u * u, 2 * u * w))
        else:
            P.append(gl(1e-12, m, n))
        u, w = gl(0, np.sqrt(b - m), n); P.append((b - u * u, 2 * u * w)) if b != 1.8 else P.append(gl(m, b, n))
    P.append(gl(1.8, 3., n)); P.append(gl(3., 8., 2 * n))
    return np.concatenate([p[0] for p in P]), np.concatenate([p[1] for p in P])


res = dict(obstruction_k=[float(k0a), float(k0b)])
for name, rule in rules.items():
    row = {}
    for n in (64, 128):
        kk, ww = nodes(n); vq = np.zeros(len(kk)); vp = np.zeros(len(kk))
        for i, k in enumerate(kk):
            a = splA(k * k).reshape(4, 6); b = splB(k * k).reshape(4, 6); w = rule(k)
            v = a ** 2 / (2 * w[None, :]) + b ** 2 * w[None, :] / 2; vq[i] = v[:2].sum(); vp[i] = v[2:].sum()
        fac = ww * kk ** 2 / (2 * np.pi ** 2)
        row[f'n{n}'] = dict(q_RMS_over_mean=float(np.sqrt(fac @ vq) / r0), p_RMS_over_mean=float(np.sqrt(fac @ vp) / p0),
                            q_share_k_lt_1p8=float((fac * vq)[kk < 1.8].sum() / (fac @ vq)))
    res[name] = row
    print(name, row, flush=True)

# Fourier-mode counts below each obstruction
counts = {}
for L in (12., 24.):
    n = int(np.ceil(1.8 * L / (2 * np.pi))) + 1; a = np.arange(-n, n + 1)
    kk = 2 * np.pi / L * np.sqrt((a[:, None, None] ** 2 + a[None, :, None] ** 2 + a[None, None, :] ** 2).ravel())
    kk = kk[kk > 0]
    counts[f'L{L:g}'] = dict(modes_below_k0b=int(np.sum(kk < k0b)), modes_below_k0a=int(np.sum(kk < k0a)),
                             shells_below_k0b=int(len(np.unique(np.round(kk[kk < k0b], 9)))), shells_below_k0a=int(len(np.unique(np.round(kk[kk < k0a], 9)))))
res['mode_counts'] = counts
print(counts)
json.dump(res, open(out, 'w'), indent=1)
