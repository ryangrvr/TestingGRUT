#!/usr/bin/env python3
"""Controlled 3-D NR onset only. No soliton identification or target charge."""
import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np
from scipy import fft, ndimage

ROOT = Path(__file__).resolve().parents[1]
L = 4 * np.pi
N0 = 1.8218345666575245e-13
M = .1
ATTR = .009
K0 = np.sqrt(ATTR * N0 / M)
Q_UNIT = N0 / K0**3


def initial(n, seed):
    # Identical finite continuous Fourier series at every resolution.
    rng = np.random.Generator(np.random.PCG64(seed))
    modes = [(i, j, k) for i in range(-3, 4) for j in range(-3, 4)
             for k in range(-3, 4) if .5 <= .5*np.linalg.norm([i, j, k]) <= 1.5]
    c = rng.normal(size=len(modes)) + 1j*rng.normal(size=len(modes))
    c /= np.linalg.norm(c)
    f = np.zeros((n, n, n), dtype=np.complex128)
    for ijk, a in zip(modes, c):
        f[tuple(x % n for x in ijk)] = a*n**3
    psi = 1 + .001*fft.ifftn(f, workers=2)
    psi /= np.sqrt(np.mean(np.abs(psi)**2))
    return psi


def periodic_components(rho, threshold, dx):
    labels, count = ndimage.label(rho >= threshold)
    parent = np.arange(count+1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for axis in range(3):
        a, b = np.take(labels, 0, axis), np.take(labels, -1, axis)
        mask = (a > 0) & (b > 0)
        for p, q in np.unique(np.stack([a[mask], b[mask]], axis=1), axis=0):
            parent[find(int(q))] = find(int(p))
    lookup = np.array([find(i) for i in range(count+1)])
    labels = lookup[labels]
    ids = np.unique(labels)
    out = []
    for i in ids[ids > 0]:
        mask = labels == i
        out.append(dict(charge=float(np.sum(rho[mask])*dx**3*Q_UNIT),
                        volume_tilde=float(np.sum(mask)*dx**3),
                        peak=float(np.max(rho[mask])), cells=int(np.sum(mask))))
    return sorted(out, key=lambda x: -x['charge'])


def run(n=64, seed=17, max_dt=.02, end=12., output=None):
    start = time.monotonic()
    tag = f'n{n}_seed{seed}_dt{max_dt:g}'
    output = output or ROOT/'data'/tag
    output.parent.mkdir(exist_ok=True, parents=True)
    psi = initial(n, seed)
    dx = L/n
    k = 2*np.pi*fft.fftfreq(n, d=dx)
    k2 = k[:, None, None]**2+k[None, :, None]**2+k[None, None, :]**2
    # Radial upper third of the axis Nyquist frequency, not corner Nyquist.
    high = k2 > ((2/3)*np.pi/dx)**2
    coarse_index = (fft.fftfreq(32)*32).astype(int) % n
    snapshots = {}
    components = {}
    diagnostics = []
    t = 0.
    steps = 0
    min_dt = max_dt
    save_times = [0., 6., 7., 7.5, 8., 8.5, 9., 10., 12.]
    events = sorted(set(save_times + list(np.round(np.arange(.1, end+.001, .1), 9))))
    previous_safe = None
    reason = 'end_time'

    def record():
        rho = np.abs(psi)**2
        f = fft.fftn(psi, workers=2)
        power = np.abs(f)**2/n**6
        gradient = np.sum(k2*power)
        h = (gradient-.5*np.mean(rho**2))*L**3
        row = dict(tau=float(t), steps=steps, mean_density=float(np.mean(rho)),
                   peak_density=float(np.max(rho)), hamiltonian=float(h),
                   high_k_norm_fraction=float(np.sum(power[high])/np.sum(power)),
                   high_k_gradient_fraction=float(np.sum((k2*power)[high])/max(gradient, 1e-300)))
        diagnostics.append(row)
        return row, rho, f

    def save_state(rho, f):
        key = f'tau_{t:.9g}'
        fc = f[np.ix_(coarse_index, coarse_index, coarse_index)]*(32/n)**3
        snapshots[key+'_psi32'] = fft.ifftn(fc).astype(np.complex64)
        snapshots[key+'_density_midplane'] = rho[:, :, n//2].astype(np.float32)
        components[key] = {str(th): periodic_components(rho, th, dx) for th in (2, 4, 8)}

    def is_unsafe(row):
        return row['high_k_norm_fraction'] > 1e-6 or row['high_k_gradient_fraction'] > 1e-3

    for target in events:
        if target > end+1e-10:
            continue
        while t < target-1e-11:
            dt = min(max_dt, .04/np.max(np.abs(psi)**2), target-t)
            min_dt = min(min_dt, dt)
            psi *= np.exp(.5j*dt*np.abs(psi)**2)
            f = fft.fftn(psi, workers=2)
            f *= np.exp(-1j*dt*k2)
            psi = fft.ifftn(f, workers=2)
            psi *= np.exp(.5j*dt*np.abs(psi)**2)
            t += dt
            steps += 1
            # Collapse accelerates within a .1 output interval. Inspect every
            # five steps once rho>8 rather than stepping blindly to that output.
            if steps % 5 == 0 and np.max(np.abs(psi)**2) > 8:
                row, rho, f = record()
                if is_unsafe(row):
                    reason = 'spectral_stop'
                    save_state(rho, f)
                    print(tag, 'STOP', json.dumps(row), flush=True)
                    break
                previous_safe = row
        if reason == 'spectral_stop':
            break
        t = target
        row, rho, f = record()
        unsafe = is_unsafe(row)
        if any(abs(t-s) < 1e-9 for s in save_times) or unsafe:
            save_state(rho, f)
        if unsafe:
            reason = 'spectral_stop'
            print(tag, 'STOP', json.dumps(row), flush=True)
            break
        previous_safe = row
        if round(t*10) % 10 == 0:
            print(tag, 'tau', round(t, 2), 'peak', round(row['peak_density'], 5), flush=True)
    result = dict(n=n, seed=seed, maximum_dt=max_dt, minimum_dt=min_dt,
                  box_length=L, initial_charge=N0*(L/K0)**3,
                  physical_charge_per_dimensionless_integral=Q_UNIT,
                  reason=reason, final=diagnostics[-1], last_safe=previous_safe,
                  max_norm_drift=max(abs(r['mean_density']-1) for r in diagnostics),
                  max_relative_hamiltonian_drift=max(abs(r['hamiltonian']/diagnostics[0]['hamiltonian']-1)
                                                      for r in diagnostics),
                  elapsed_s=time.monotonic()-start, component_diagnostics=components)
    Path(str(output)+'.json').write_text(json.dumps(result, indent=2)+'\n')
    with Path(str(output)+'.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=diagnostics[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(diagnostics)
    np.savez_compressed(str(output)+'.npz', **snapshots)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=64)
    parser.add_argument('--seed', type=int, default=17)
    parser.add_argument('--max-dt', type=float, default=.02)
    parser.add_argument('--end', type=float, default=12.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run(args.n, args.seed, args.max_dt, args.end, args.output)
    print(json.dumps({k: v for k, v in result.items() if k != 'component_diagnostics'}, indent=2))
