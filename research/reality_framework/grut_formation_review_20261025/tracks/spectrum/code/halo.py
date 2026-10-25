#!/usr/bin/env python3
"""Direction and momentum content of the slow 'halo' between the core and r=240 in a snapshot.
The positive-frequency (particle) field P(r) = sum_n alpha_n s_n(r) is formed from box-mode amplitudes of the
windowed exterior field; local Fourier content at +k (e^{+ikr}) is outgoing, at -k incoming.
Usage: halo.py SNAP R1 R2 OUTJSON"""
import json, os, sys
import numpy as np
from scipy.fft import dst
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec_lib as S

snap = S.load_snap(sys.argv[1]); R1, R2 = float(sys.argv[2]), float(sys.argv[3])
N, h = snap['N'], snap['h']; L = N * h
r = np.arange(1, N) * h; n = np.arange(1, N); k = n * np.pi / L
w = np.sqrt(1 + S.K2_6(k, h))
win = S.ramp(r, R1, 15.) * (1 - S.ramp(r, R2 - 15., 15.))
T = lambda a: dst(a, type=1, norm='ortho')
U = (snap['x'] + 1j * snap['y'])[1:N] * win; V = (snap['vx'] + 1j * snap['vy'])[1:N] * win
u = T(U.real) + 1j * T(U.imag); v = T(V.real) + 1j * T(V.imag)
a = 0.5 * (u + 1j * v / w); b = 0.5 * (u - 1j * v / w)
NX = 4 * np.pi * h * w * abs(a)**2; NXb = 4 * np.pi * h * w * abs(b)**2
# number-weighted field: phi = IDST(sqrt(4 pi h w) a), so sum|phi|^2 = NX
ca = np.sqrt(4 * np.pi * h * w) * a; phi = T(ca.real) + 1j * T(ca.imag)
m = (r > R1 - 5) & (r < R2 + 5); seg = phi[m]
F = np.abs(np.fft.fft(seg, 8 * len(seg)))**2; kk = 2 * np.pi * np.fft.fftfreq(8 * len(seg), h)
out_frac = float(F[kk > 0].sum() / F.sum())
kb = [0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.2, 12]
res = dict(snapshot_time=snap['t'], R1=R1, R2=R2, NX_in_window=float(NX.sum()), NXbar_in_window=float(NXb.sum()),
           Q_window_direct=float(4 * np.pi * h * np.sum(-np.imag(np.conj(U) * V))),
           outgoing_fraction_of_particle_number=out_frac,
           X_k_bins=kb, X_by_k=[float(NX[(k >= lo) & (k < hi)].sum()) for lo, hi in zip(kb[:-1], kb[1:])],
           outgoing_by_k=[float(F[(kk >= lo) & (kk < hi)].sum() / max(F[(np.abs(kk) >= lo) & (np.abs(kk) < hi)].sum(), 1e-300)) for lo, hi in zip(kb[:-1], kb[1:])])
json.dump(res, open(sys.argv[4], 'w'), indent=1)
print(json.dumps(res, indent=1))
