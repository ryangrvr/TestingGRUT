#!/usr/bin/env python3
"""Phase-space occupation of the extracted s-wave quanta (validity of the classical-to-quanta map).
f(r,k) = dN / (dr dk / pi) for the s-wave (l=0) radial problem, from a Gabor analysis of the
particle (alpha) and antiparticle (beta) fields of a snapshot outside Rw.  Physical occupation = f / g_X^2.
Usage: occupation.py SNAP Rw OUTJSON"""
import json, sys, os
import numpy as np
from scipy.fft import dst
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec_lib as S

snap = S.load_snap(sys.argv[1]); Rw = float(sys.argv[2])
N, h = snap['N'], snap['h']; L = N * h
r = np.arange(1, N) * h; n = np.arange(1, N); k = n * np.pi / L
w = np.sqrt(1 + S.K2_6(k, h)); win = S.ramp(r, Rw, 30.)
T = lambda a: dst(a, type=1, norm='ortho')
U = (snap['x'] + 1j * snap['y'])[1:N] * win; V = (snap['vx'] + 1j * snap['vy'])[1:N] * win
u = T(U.real) + 1j * T(U.imag); v = T(V.real) + 1j * T(V.imag)
res = {}
for name, amp in [('X', 0.5 * (u + 1j * v / w)), ('Xbar', 0.5 * (u - 1j * v / w))]:
    c = np.sqrt(4 * np.pi * h * w) * amp
    phi = T(c.real) + 1j * T(c.imag)           # sum |phi|^2 = N
    Ntot = float(np.sum(np.abs(phi)**2))
    Wd = 400.0; nw = int(Wd / h); step = nw // 2
    g = np.hanning(nw); g /= np.sqrt(np.sum(g**2))
    best = (0, None, None); cells = []
    kk = 2 * np.pi * np.fft.fftfreq(nw, h)
    for s0 in range(0, len(phi) - nw, step):
        seg = phi[s0:s0 + nw] * g
        F = np.abs(np.fft.fft(seg))**2           # sum_k F = sum |seg|^2 ~ number in this window (Hann-normalised)
        # F sums to the number in the window.  A |k|-bin (one +k and one -k FFT bin, width 2pi/Wd) holds
        # Wd*(2pi/Wd)/pi = 2 s-wave states, so the occupation is f = (F(k) + F(-k)) / 2.
        f = 0.5 * (F + F[(-np.arange(nw)) % nw])
        i = np.argmax(f)
        if f[i] > best[0]: best = (float(f[i]), float(r[s0 + nw // 2]), float(abs(kk[i])))
    res[name] = dict(N=Ntot, max_occupation_dimensionless=best[0], at_r=best[1], at_k=best[2],
                     max_occupation_physical_gX_0p1=best[0] / 0.01)
res['note'] = ('Gabor estimate with 400-unit Hann windows; physical occupation = f/g_X^2 (g_X=0.1 quoted from bb211d1). '
               'Quantum (vacuum) corrections to the coherent-state map are relative order 1/f_phys per occupied cell.')
json.dump(res, open(sys.argv[3], 'w'), indent=1)
print(json.dumps(res, indent=1))
