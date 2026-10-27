"""Skeptic (en_indep_skeptic): turn a Codex-engine snapshot (early3d.cpp format) into a GRUT3D1 initial file
for a restart, optionally adding a perturbation. Own code (format read from early3d.cpp / en_common.write_initial).

usage: snap2init.py SNAP OUT.bin KIND EPS SEED
  KIND = none | ir | white | mid
    ir    : Gaussian random field with Fourier support 0 < |k| <= 2.5 (in H_i units), q only
    mid   : support 2.5 < |k| <= 6
    white : independent N(0,1) per site (all lattice modes), q only (same type as en_indep's twin)
  Each of the three q fields gets ||dq_c|| = EPS * ||q_c - mean(q_c)||.
Prints JSON with the restart clock (tau0, a0) computed from the snapshot's elapsed eta.
"""
import sys, json, struct
import numpy as np
AM, C, LX, G2, LP = 8.869218738806, 1.9857962750634526, 1e-3, 1e-2, 3e-2
snap, out, kind, eps, seed = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5])
with open(snap, 'rb') as f:
    N, L, s, bits = struct.unpack('<iddi', f.read(24)); assert bits == 64
    a = np.fromfile(f, dtype='<f8').reshape(6, N, N, N)
q, p = a[:3].copy(), a[3:].copy()
rng = np.random.default_rng(seed)
if kind != 'none':
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    K = np.sqrt(k1[:, None, None] ** 2 + k1[None, :, None] ** 2 + k1[None, None, :] ** 2)
    for c in range(3):
        w = rng.normal(size=(N, N, N))
        if kind in ('ir', 'mid'):
            lo, hi = (0.0, 2.5) if kind == 'ir' else (2.5, 6.0)
            mask = (K > lo) & (K <= hi)
            w = np.real(np.fft.ifftn(np.fft.fftn(w) * mask))
        w -= w.mean()
        fl = q[c] - q[c].mean()
        q[c] += eps * np.linalg.norm(fl) / np.linalg.norm(w) * w
with open(out, 'wb') as f:
    f.write(b'GRUT3D1\n'); f.write(struct.pack('<id4d', N, L, 0., G2, LX, LP))
    q.astype('<f8').tofile(f); p.astype('<f8').tofile(f)
tau0 = 20 + AM * s + 0.5 * C * s * s; a0 = AM + C * s
print(json.dumps(dict(snapshot=snap, s=s, tau0=tau0, a0=a0, C=C, kind=kind, eps=eps, seed=seed)))
