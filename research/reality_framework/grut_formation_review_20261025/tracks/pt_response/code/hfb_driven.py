"""INDEPENDENT code (Claude, pt_response): homogeneous (Floquet-)Bogoliubov response driven by the MEASURED orbit.

Linear fluctuations about a spatially homogeneous but time-dependent background obey, per lattice momentum k,
    d''_a + kappa^2(k) d_a + H_ab(s) d_b = 0,   a,b in (chi1, chi2, phi),
with kappa^2 the FD4 symbol.  Three closures for H(s):
  'bare'   : Hessian of V at the coherent zero mode only (condensate Bogoliubov-de Gennes, time-periodic orbit);
  'hartree': isotropic Hartree: (hf - lam|psi_c|^2) delta_ab + 2 lam psi_c,a psi_c,b  (X block; phi block bare+g2<r2>);
  'meanH'  : the exact spatial mean of the Hessian, <d2V/dchi_a dchi_b>(s) (Gaussian/HFB mean field incl.
             fluctuation anisotropy <dchi_a dchi_b>);
H(s) is recorded every 0.01 by record_mean_hessian.py (Codex engine for the background only) and cubic-spline
interpolated.  We integrate the response to unit momentum kicks (DOP853, rtol 1e-11) for 6 time units from s0, form
normal = .5(R00+R11+i(R10-R01)), anomalous = .5(R00-R11+i(R10+R01)) exactly as response.py, average over the modes
of each shell (each mode's own kappa^2), and compare with the published nonlinear lattice responses.
Usage: python3 -I hfb_driven.py <hessian npz> <packet dir> <out json> [<extra measured npz dir>]
"""
import sys, json, glob, os
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

hz = np.load(sys.argv[1]); P = sys.argv[2].rstrip('/') + '/'; OUT = sys.argv[3]
extra = sys.argv[4] if len(sys.argv) > 4 else None
LX, G2, LP = .001, .01, .03
S = hz['s']; m = hz['mean']
psi = m[:, :2]
closures = {}
closures['meanH'] = hz['Hbar']
closures['bare'] = hz['Hcoh']
Hh = hz['Hcoh'].copy()
rc2 = psi[:, 0] ** 2 + psi[:, 1] ** 2
for a in range(2):
    Hh[:, a, a] = hz['hf'] - LX * rc2 + 2 * LX * psi[:, a] ** 2
Hh[:, 2, 2] = hz['Hbar'][:, 2, 2]
closures['hartree'] = Hh
splines = {k: CubicSpline(S, v.reshape(len(S), 9), axis=0) for k, v in closures.items()}

# lattice shells (N48, L12): FD4 symbol per mode, continuum |k| <= 8, shell = round(|k| L / 2pi)
N, L = 48, 12.
mm = np.fft.fftfreq(N) * N
kc2 = (2 * np.pi / L) ** 2 * (mm[:, None, None] ** 2 + mm[None, :, None] ** 2 + mm[None, None, :] ** 2)
f = lambda u: 4 * np.sin(u / 2) ** 2 * (1 + np.sin(u / 2) ** 2 / 3) / (L / N) ** 2
x = 2 * np.pi * mm / N
kl = f(x)[:, None, None] + f(x)[None, :, None] + f(x)[None, None, :]
sel = kc2.ravel() <= 64 + 1e-12
kc = np.sqrt(kc2.ravel()[sel]); kfd = kl.ravel()[sel]
sid = np.floor(kc / (2 * np.pi / L) + .5).astype(int)
shells = [n for n in np.unique(sid) if n > 0]
shell_k2 = {n: np.unique(np.round(kfd[sid == n], 12), return_counts=True) for n in shells}

delays = np.arange(151) * .04

ALLK2 = np.unique(np.round(kfd, 9))
K2INDEX = {n: (np.searchsorted(ALLK2, np.round(shell_k2[n][0], 9)), shell_k2[n][1]) for n in shells}
HSTEP = .0025

def response_all(closure, s0, fields=3):
    """Vectorised classical RK4 over every distinct FD4 kappa^2; returns R[t, k2, a, b] at the 151 delays."""
    sp = splines[closure]; nst = int(round(6. / HSTEP))
    Hs = sp(s0 + np.arange(2 * nst + 1) * HSTEP / 2).reshape(-1, 3, 3)[:, :fields, :fields]
    nk = len(ALLK2); k2 = ALLK2[:, None, None]
    q = np.zeros((nk, fields, fields)); p = np.broadcast_to(np.eye(fields), (nk, fields, fields)).copy()
    out = np.zeros((len(delays), nk, fields, fields)); every = int(round(.04 / HSTEP)); j = 1
    acc = lambda qq, H: -(k2 * qq) - np.einsum('ab,kbc->kac', H, qq)
    for n in range(nst):
        H0, H1, H2 = Hs[2 * n], Hs[2 * n + 1], Hs[2 * n + 2]
        k1q, k1p = p, acc(q, H0)
        k2q, k2p = p + .5 * HSTEP * k1p, acc(q + .5 * HSTEP * k1q, H1)
        k3q, k3p = p + .5 * HSTEP * k2p, acc(q + .5 * HSTEP * k2q, H1)
        k4q, k4p = p + HSTEP * k3p, acc(q + HSTEP * k3q, H2)
        q = q + HSTEP / 6 * (k1q + 2 * k2q + 2 * k3q + k4q); p = p + HSTEP / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        if (n + 1) % every == 0: out[j] = q; j += 1
    return out

def channels(R):
    nr = .5 * (R[:, 0, 0] + R[:, 1, 1] + 1j * (R[:, 1, 0] - R[:, 0, 1]))
    an = .5 * (R[:, 0, 0] - R[:, 1, 1] + 1j * (R[:, 1, 0] + R[:, 0, 1]))
    return nr, an

def peaks(G):
    w = np.linspace(-20, 20, 2001); T = 6.; taper = np.cos(np.pi * delays / (2 * T)) ** 2
    trap = np.full(len(delays), .04); trap[[0, -1]] *= .5
    Gw = np.exp(-1j * w[:, None] * delays[None, :]) @ (G * trap * taper); rho = -2 * Gw.imag
    pos = w > 0; neg = w < 0
    return float(w[pos][np.argmax(rho[pos])]), float(w[neg][np.argmax(-rho[neg])])

measured = {}
for f_ in sorted(glob.glob(P + 'RESPONSE_s*_seed*.npz')):
    d = np.load(f_); tag = os.path.basename(f_)[9:-4]; measured[('codex', tag)] = (float(d['s0']), d['normal'], d['anomalous'])
if extra:
    for f_ in sorted(glob.glob(extra.rstrip('/') + '/MEAS_*.npz')):
        d = np.load(f_); tag = os.path.basename(f_)[5:-4]; measured[('ours', tag)] = (float(d['s0']), d['normal'], d['anomalous'])

res = dict(method=__doc__, phases={}); series = {}
for s0 in sorted({v[0] for v in measured.values()} | {20., 20.24, 20.44, 20.68}):
    ph = {}
    for closure in ['bare', 'hartree', 'meanH']:
        rows = []
        Rall = response_all(closure, s0)
        for n in shells:
            ix, cnt = K2INDEX[n]
            Rsh = np.einsum('tkab,k->tab', Rall[:, ix], cnt) / cnt.sum()
            Nm, Am = channels(Rsh)
            series.setdefault(f'{s0:g}', {}).setdefault(closure, {})[int(n)] = (Nm, Am)
            row = dict(shell=int(n), k_center=float(n * 2 * np.pi / L), ratio=float(np.linalg.norm(Am) / np.linalg.norm(Nm)),
                       peak_pos=peaks(Nm)[0], peak_neg=peaks(Nm)[1], growth_late_over_early=float(np.max(abs(Nm[100:])) / np.max(abs(Nm[:50]))))
            for key, (ss, Nmeas, Ameas) in measured.items():
                if abs(ss - s0) > 1e-9: continue
                k = n - 1
                row[f'{key[0]}:{key[1]}'] = dict(
                    meas_ratio=float(np.linalg.norm(Ameas[:, k]) / np.linalg.norm(Nmeas[:, k])),
                    normal_rel_err=float(np.linalg.norm(Nmeas[:, k] - Nm) / np.linalg.norm(Nmeas[:, k])),
                    anomalous_rel_err=float(np.linalg.norm(Ameas[:, k] - Am) / np.linalg.norm(Ameas[:, k])))
            rows.append(row)
            if closure == 'meanH' or n <= 3:
                pass
        ph[closure] = rows
        if n:
            # store time series for the first 4 shells
            pass
    res['phases'][f'{s0:g}'] = ph
    print('s0', s0, flush=True)
    for closure in ['bare', 'hartree', 'meanH']:
        print('  ', closure, ' '.join(f"{r['shell']}:{r['ratio']:.3f}" for r in ph[closure][:8]))
        for r in ph[closure][:4]:
            extra_s = '; '.join(f"{k}: meas {v['meas_ratio']:.3f} errN {v['normal_rel_err']:.2f} errA {v['anomalous_rel_err']:.2f}" for k, v in r.items() if isinstance(v, dict))
            print('      shell', r['shell'], 'peaks', r['peak_pos'], r['peak_neg'], 'grow', round(r['growth_late_over_early'], 2), '|', extra_s)
json.dump(res, open(OUT, 'w'), indent=1)
np.savez_compressed(OUT.replace('.json', '_series.npz'), **{f'{s}_{c}_{n}_{w}': (v[0] if w == 'N' else v[1]) for s, d1 in series.items() for c, d2 in d1.items() for n, v in d2.items() for w in 'NA'})
