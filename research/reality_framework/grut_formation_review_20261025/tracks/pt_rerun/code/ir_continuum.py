"""Continuum (infinite-volume) estimate of the tau=20 handoff RMS for the stated IR completion.

Uses Codex's transfer machinery (prepare.M.integrate / neutral, prepare.factor; reproduction) on a fine k grid with
dense nodes near the obstruction momenta 1.414 and 1.683, and our own quadrature
sigma^2 = (1/2pi^2) int k^2 dk [F F^T]_{xx} (own code). The per-mode stable-side variance diverges like 1/Omega,
Omega = sqrt(k^2 - k0^2), which is integrable; the substitution k = k0 + u^2 removes the singularity.
Usage (from the packet code dir): python3 -B ir_continuum.py <packet_dir> <early_dir> <out.json>
"""
import sys, os, json
pk, early, out = sys.argv[1:4]
NS = float(sys.argv[4]) if len(sys.argv) > 4 else 1.  # node-count multiplier for a convergence check
sys.path.insert(0, os.path.join(pk, 'code'))
import numpy as np
import prepare as P  # Codex code (reproduction machinery)
from pt_common import E, ROOT

M = E.M; exe = M.build()
ne = exe.parent / 'neutral_pt'
import subprocess
subprocess.run(['gcc', '-O3', '-std=c11', '-DNEUTRAL_SECTOR', '-I', str(exe.parent), str(E.T.MODE / 'code/modes.c'), '-lm', '-o', str(ne)], check=True)  # as in prepare.py
old = np.load(os.path.join(early, 'MATCH_L12.npz')); lam = old['initial_mass_eigenvalues']; U = old['initial_mass_eigenvectors']
v0, _ = M.prepared(3.); phi_mass = P.NE.mass(v0, 3.)
k0 = np.sort(np.sqrt(-lam[lam < 0]))  # 1.414, 1.683
mq = old['light_mean_q']; mp = old['light_mean_p']; r0 = np.linalg.norm(mq[:2]); p0 = np.linalg.norm(mp[:2])

# piecewise grid: [0.02, k0a), (k0a, k0b), (k0b, 8]; near each k0 from above use k = k0 + u^2
def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n); return (a + b) / 2 + (b - a) / 2 * x, (b - a) / 2 * w
# below first obstruction and between: smooth on the unstable side (finite width) -> plain Gauss-Legendre
nodes = []; weights = []
x, w = gl(0., k0[0], int(48 * NS)); nodes.append(x); weights.append(w)
# (k0a, k0b): stable dir1 near k0a from above -> k = k0a + u^2
ua = np.sqrt(k0[1] - k0[0]); x, w = gl(0., ua, int(48 * NS)); nodes.append(k0[0] + x ** 2); weights.append(w * 2 * x)
# (k0b, 3): k = k0b + u^2
ub = np.sqrt(3. - k0[1]); x, w = gl(0., ub, int(64 * NS)); nodes.append(k0[1] + x ** 2); weights.append(w * 2 * x)
x, w = gl(3., 8., int(64 * NS)); nodes.append(x); weights.append(w)
ks = np.concatenate(nodes); wq = np.concatenate(weights)
order = np.argsort(ks); ks = ks[order]; wq = wq[order]
Fw = {1.: np.zeros((len(ks), 6, 14)), 3.: np.zeros((len(ks), 6, 14))}
for start in range(0, len(ks), 96):
    kk = ks[start:start + 96]
    bg, T, c = M.integrate(exe, f'ptrr_cont_{start}', 3., .5, kk, np.array([0., 20.]), rtol=2e-12, e=0.)
    pb, PT, pc = P.NE.run(ne, dict(tag=f'ptrr_cont_{start}', beta=3., Gamma_over_H=.5), kk, np.array([0., 20.]), rtol=2e-12, suffix='_pt')
    for width in (1., 3.):
        Fw[width][start:start + 96] = P.factor(T[-1], PT[-1], kk, lam, U, phi_mass, width)
    print('batch', start, len(kk), flush=True)
res = dict(obstruction_k=k0.tolist(), nodes=len(ks), node_multiplier=NS)
for width in (1., 3.):
    F = Fw[width]; C = np.einsum('kia,kja->kij', F, F)
    q = C[:, 0, 0] + C[:, 1, 1]; p = C[:, 3, 3] + C[:, 4, 4]
    fac = ks ** 2 / (2 * np.pi ** 2)
    Q2 = np.sum(wq * fac * q); P2 = np.sum(wq * fac * p)
    lo = ks < 1.8
    Q2old = np.sum((wq * fac * q)[~lo]); P2old = np.sum((wq * fac * p)[~lo])
    near = (ks > k0[1]) & (ks < 1.8)
    res[f'width{width:g}'] = dict(q_RMS_over_mean=float(np.sqrt(Q2) / r0), p_RMS_over_mean=float(np.sqrt(P2) / p0),
                                  q_RMS_band_ge_1p8=float(np.sqrt(Q2old) / r0), q_variance_share_k_lt_1p8=float(1 - Q2old / Q2),
                                  q_variance_share_k0b_to_1p8=float(np.sum((wq * fac * q)[near]) / Q2),
                                  p_variance_share_k_lt_1p8=float(1 - P2old / P2))
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
