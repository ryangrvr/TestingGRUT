"""Skeptic (en_indep_skeptic): independent recomputation of Codex's fixed-homogeneous-orbit linear phi control (claim 3b).
Own code and own method, nothing imported from Codex:
  * homogeneous reduced X orbit  Xbar'' = -lambda_X |Xbar|^2 Xbar  from the spatial mean of the tau=20 field;
  * linear phi with Omega_k^2(s) = K_FD4^2 + g^2 |Xbar(s)|^2, evolved as SECOND MOMENTS per shell
    (A=<phi^2>, B=<phi pi>, D=<pi^2>:  A'=2B, B'=D-Omega^2 A, D'=-2 Omega^2 B) -- Codex uses 2x2 transfer matrices + DOP853;
  * fixed-step classical RK4 in conformal time s (own stepper), step h, two step sizes for convergence;
  * shells = unique FD4 symbol values, moments from a full complex FFT (Codex: rfft with half-space weights).
Energy convention copied from LINEAR_PHI_CONTROL.json scope: (K_phi + G_phi + g^2/2 |Xbar|^2 phi^2) summed over all phi
modes (incl. the k=0 mode), no lambda_phi quartic, no vacuum subtraction, divided by the total active energy at tau=20.
usage: linear_phi_moments.py INITIAL.bin OUT.json
"""
import sys, json, struct
import numpy as np
LX, G2, LP = 1e-3, 1e-2, 3e-2
AM, C = 8.869218738806, 1.9857962750634526
path, outp = sys.argv[1], sys.argv[2]
with open(path, 'rb') as f:
    assert f.read(8) == b'GRUT3D1\n'
    N, L, nu, g2, lx, lp = struct.unpack('<id4d', f.read(44))
    a = np.fromfile(f, dtype='<f8').reshape(6, N, N, N)
q, p = a[:3], a[3:]
dx = L / N; th = 2*np.pi*np.fft.fftfreq(N); s2 = np.sin(th/2)**2; f1 = 4*s2*(1+s2/3)/dx**2
K2 = f1[:, None, None] + f1[None, :, None] + f1[None, None, :]
# own total energy at tau=20 (spectral FD4 gradient energy)
def energy(q, p):
    qh = np.fft.fftn(q, axes=(1, 2, 3))
    grad = 0.5*np.sum(K2*np.abs(qh)**2)/N**6
    r = q[0]**2 + q[1]**2
    return float(0.5*np.mean(np.sum(p*p, axis=0)) + grad + np.mean(LX/4*r*r + G2/2*q[2]**2*r + LP/4*q[2]**4))
rho0 = energy(q, p)
ph = np.fft.fftn(q[2])/N**3; pp = np.fft.fftn(p[2])/N**3
use = (np.abs(ph) + np.abs(pp)) > 1e-12
vals, inv = np.unique(np.round(K2[use], 10), return_inverse=True)
A0 = np.bincount(inv, weights=np.abs(ph[use])**2); D0 = np.bincount(inv, weights=np.abs(pp[use])**2)
B0 = np.bincount(inv, weights=np.real(np.conj(ph[use])*pp[use]))
xb = q[:2].mean(axis=(1, 2, 3)); vb = p[:2].mean(axis=(1, 2, 3))
def rhs(y):
    x, v = y[0:2], y[2:4]; A, B, D = y[4:4+n], y[4+n:4+2*n], y[4+2*n:]
    r = x @ x; W = vals + G2*r
    return np.concatenate([v, -LX*r*x, 2*B, D - W*A, -2*W*B])
n = len(vals)
def run(h):
    taus = np.arange(20.0, 800.0 + 1e-9, 0.5)
    ss = 2*(taus-20)/(AM + np.sqrt(AM*AM + 2*C*(taus-20)))
    y = np.concatenate([xb, vb, A0, B0, D0]); s = 0.0; out = []
    for st in ss:
        m = int(np.ceil((st - s)/h - 1e-12))
        if m > 0:
            hh = (st - s)/m
            for _ in range(m):
                k1 = rhs(y); k2 = rhs(y + .5*hh*k1); k3 = rhs(y + .5*hh*k2); k4 = rhs(y + hh*k3)
                y = y + hh/6*(k1 + 2*k2 + 2*k3 + k4)
        s = st
        x = y[0:2]; r = x @ x; A, D = y[4:4+n], y[4+2*n:]
        out.append(0.5*np.sum(D + (vals + G2*r)*A)/rho0)
    return taus, np.array(out), y
t, e1, y1 = run(0.004)
_, e2, y2 = run(0.002)
res = dict(method=__doc__, nshells=int(n), rho_tau20=rho0, phi_zero_mode_included=bool(use[0, 0, 0]),
           fraction_tau20=float(e1[0]), fraction_tau800_h004=float(e1[-1]), fraction_tau800_h002=float(e2[-1]),
           first_tau_fraction_ge_1=float(t[np.flatnonzero(e2 >= 1)[0]]) if np.any(e2 >= 1) else None,
           max_fraction=float(e2.max()),
           rows={f'{T:g}': float(v) for T, v in zip(t, e2) if abs(T % 50) < 1e-9 or T in (20.5, 60.5, 702.5)})
json.dump(res, open(outp, 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ('method', 'rows')}))
print({k: round(v, 6) for k, v in res['rows'].items()})
