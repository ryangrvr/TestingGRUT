#!/usr/bin/env python3
"""INDEPENDENT estimator (own code) of the absorber sensitivity of the conditional antiparticle
fraction in 9636950, on the regenerated bb211d1 probe traces (N0125: eta=.03; N0125_abs: eta=.015).

Two own estimators, both conditional linear-mode diagnostics (not particle abundances):
 (A) outgoing projection, own implementation: a(w) = [F + (D + F/r)/(i k)]/2, k = sign(w) sqrt(w^2-1),
     number density |k| |a|^2, Hann or Tukey window, |k| > kmin.
 (B) net measured energy flux per frequency divided by |w| (no in/out split): fraction of
     sum_{w<0} max(Eflux,0)/|w| over sum max(Eflux,0)/|w|.
Usage: python3 indep_antiparticle.py TRACE_DIR OUT.json
"""
import json, sys
from pathlib import Path
import numpy as np


def window(n, kind):
    if kind == 'hann':
        return np.hanning(n)
    a = float(kind.split(':')[1]); w = np.ones(n); m = int(a * (n - 1) / 2)
    ramp = .5 * (1 - np.cos(np.pi * np.arange(m) / m)); w[:m] = ramp; w[n - m:] = ramp[::-1]
    return w


def estimate(trace, r, kind, start=0., kmin=.01):
    tr = trace[trace[:, 0] >= start]; t = tr[:, 0]; dt = t[1] - t[0]; n = len(t)
    w = window(n, kind)
    psi = tr[:, 1] + 1j * tr[:, 2]; dpsi = tr[:, 4] + 1j * tr[:, 5]; pt = tr[:, 7] + 1j * tr[:, 8]
    F = np.fft.fft(w * psi); D = np.fft.fft(w * dpsi); P = np.fft.fft(w * pt)
    om = -2 * np.pi * np.fft.fftfreq(n, dt)
    k2 = om ** 2 - 1; ok = k2 > kmin ** 2
    k = np.sign(om[ok]) * np.sqrt(k2[ok])
    a = .5 * (F[ok] + (D[ok] + F[ok] / r) / (1j * k))
    num = np.abs(k) * np.abs(a) ** 2
    neg = om[ok] < 0
    fracA = num[neg].sum() / num.sum()
    eflux = -np.real(np.conj(P) * D)        # outward energy flux density (up to positive constants)
    nb = np.maximum(eflux[ok], 0) / np.abs(om[ok])
    fracB = nb[neg].sum() / nb.sum()
    return float(fracA), float(fracB)


def main(tdir, out):
    tdir = Path(tdir); rows = {}
    for name in ['N0125', 'N0125_abs', 'N0125_box', 'N025', 'N00625']:
        f = tdir / f'{name}_trace.npz'
        if not f.exists():
            continue
        trace = np.load(f)['trace']
        h = {'N0125': .125, 'N0125_abs': .125, 'N0125_box': .125, 'N025': .25, 'N00625': .0625}[name]
        r = round(120 / h) * h
        rows[name] = {}
        for kind in ['hann', 'tukey:0.05', 'tukey:0.2']:
            for start in [0., 6000.]:
                rows[name][f'{kind}|start{int(start)}'] = estimate(trace, r, kind, start)
    sens = {}
    if 'N0125' in rows and 'N0125_abs' in rows:
        for key in rows['N0125']:
            a, b = rows['N0125'][key], rows['N0125_abs'][key]
            sens[key] = dict(A_relative_shift=b[0] / a[0] - 1, B_relative_shift=b[1] / a[1] - 1)
    json.dump(dict(label='INDEPENDENT (own spectral code) conditional antiparticle-fraction estimators', fractions=rows,
                   absorber_halving_sensitivity=sens), open(out, 'w'), indent=1)
    print(json.dumps(rows, indent=1)); print(json.dumps(sens, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
