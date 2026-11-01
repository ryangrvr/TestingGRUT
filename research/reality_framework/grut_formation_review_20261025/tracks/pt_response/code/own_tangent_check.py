"""INDEPENDENT integrator (Claude, pt_response; numpy only, no Codex import) for the background + tangent.

Checks the probe-study engine output: starting from the regenerated base48 field at s = 20, evolves the nonlinear
X/phi field and two tangents (pi_1 and pi_2 kicked with cos(k.x) of the first half-space shell-1 mode) with our own
FD4 (spectral application of the FD4 symbol) and a DRIFT-first Yoshida-4 composition (Codex's is kick-first), h=0.004,
and compares the shell-1..4 Fourier projections of the tangent chi_1, chi_2 with PROJ_exact1_s20.npz (tangents 0,1).
Usage: python3 -I own_tangent_check.py <state_s20.npy> <PROJ_exact1_s20.npz> <out json>
"""
import sys, json
import numpy as np
from numpy.fft import fftn, ifftn

z = np.load(sys.argv[1]); proj = np.load(sys.argv[2]); OUT = sys.argv[3]
N, L = 48, 12.; LX, G2, LP = .001, .01, .03
mm = np.fft.fftfreq(N) * N
u = 2 * np.pi * mm / N; f1 = 4 * np.sin(u / 2) ** 2 * (1 + np.sin(u / 2) ** 2 / 3) / (L / N) ** 2
sym = f1[:, None, None] + f1[None, :, None] + f1[None, None, :]
lap = lambda f: np.real(ifftn(-sym * fftn(f, axes=(-3, -2, -1)), axes=(-3, -2, -1)))
allm = proj['allmodes']
# first half-space shell-1 mode used by probe_study (meta[0])
ix0 = int(eval(str(proj['meta'][0]))[0])
i, j, k = np.unravel_index(ix0, (N, N, N)); kv = 2 * np.pi / L * np.array([mm[i], mm[j], mm[k]])
x = np.arange(N) * L / N; X, Y, Z = np.meshgrid(x, x, x, indexing='ij'); pat = np.cos(kv[0] * X + kv[1] * Y + kv[2] * Z)
q = z[:3].copy(); p = z[3:].copy()
dq = np.zeros((2, 3, N, N, N)); dp = np.zeros((2, 3, N, N, N)); dp[0, 0] = pat; dp[1, 1] = pat

def force(q):
    c1, c2, ph = q; r2 = c1 * c1 + c2 * c2; F = np.empty_like(q)
    F[0] = lap(c1) - (LX * r2 + G2 * ph * ph) * c1; F[1] = lap(c2) - (LX * r2 + G2 * ph * ph) * c2
    F[2] = lap(ph) - (LP * ph * ph + G2 * r2) * ph
    return F

def dforce(q, d):
    c1, c2, ph = q; r2 = c1 * c1 + c2 * c2; f = LX * r2 + G2 * ph * ph; dot = c1 * d[0] + c2 * d[1]; F = np.empty_like(d)
    F[0] = lap(d[0]) - f * d[0] - 2 * LX * c1 * dot - 2 * G2 * c1 * ph * d[2]
    F[1] = lap(d[1]) - f * d[1] - 2 * LX * c2 * dot - 2 * G2 * c2 * ph * d[2]
    F[2] = lap(d[2]) - (3 * LP * ph * ph + G2 * r2) * d[2] - 2 * G2 * ph * dot
    return F

w1 = 1 / (2 - 2 ** (1 / 3)); w0 = 1 - 2 * w1
cs = [w1 / 2, (w1 + w0) / 2, (w0 + w1) / 2, w1 / 2]; ds = [w1, w0, w1]
h = .004; out = []; t = 0.
cmp_samples = [10, 25, 50]
for nstep in range(1, 501):
    # drift-first composition: D(c) K(d) D(c) K(d) D(c) K(d) D(c)
    for a in range(4):
        q += cs[a] * h * p; dq += cs[a] * h * dp
        if a < 3:
            p += ds[a] * h * force(q)
            for b in range(2): dp[b] += ds[a] * h * dforce(q, dq[b])
    if nstep % 10 == 0:
        sidx = nstep // 10
        if sidx in cmp_samples:
            mine = np.array([fftn(dq[b, :2], axes=(-3, -2, -1), norm='forward').reshape(2, -1)[:, allm] for b in range(2)])
            theirs = proj['proj'][sidx, 0:2]
            rel = float(np.linalg.norm(mine - theirs) / np.linalg.norm(theirs))
            out.append(dict(delay=sidx * .04, rel_diff=rel)); print(out[-1], flush=True)
json.dump(dict(check=out, note='own numpy drift-first Yoshida vs Codex kick-first fieldstep tangent; O(h^4) differences expected'), open(OUT, 'w'), indent=1)
