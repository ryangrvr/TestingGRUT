#!/usr/bin/env python3
"""Compare mode evolution on the actual H55 spatial mean to the old fluid."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline

from kinetic_replay import load_source
from trace_audit import minimum_energy

CASES = {
    'coarse_q': (720, 12, 45000),
    'baseline': (1440, 12, 45000),
    'fine_q': (2880, 12, 45000),
    'half_step': (1440, 24, 45000),
    'wide_q': (2160, 12, 67500),
}
SOURCE = Path(__file__).with_suffix('.c')


def compile_engine(work):
    binary = work / 'mean_modes'
    subprocess.run(['gcc', '-std=c11', '-Wall', '-Wextra', '-O3', '-march=native',
                    '-fopenmp', '-o', str(binary), str(SOURCE), '-lm'], check=True)
    return binary


def evolve(binary, directory, name, t, a, H, theta, theta_dot, obs, masses, settings, threads):
    nq, substeps, qmax = settings
    assert obs[0] == 0 and obs[-1] == len(t)-1 and np.all(np.diff(obs) > 0)
    dt = (t[-1] - t[0]) / ((len(t)-1) * substeps)
    nstep = (len(t)-1) * substeps
    assert np.max(abs(np.diff(t) / dt / substeps - 1)) < 1e-10
    positions = t[0] + dt * (np.arange(nstep)[:, None] + [0.5-math.sqrt(3)/6, 0.5+math.sqrt(3)/6])
    log_a = CubicHermiteSpline(t, np.log(a), H)
    phases = CubicHermiteSpline(t, theta, theta_dot, axis=0)
    drive = np.concatenate([np.exp(-log_a(positions))[..., None], phases(positions, 1) / 2], axis=2)
    observation = np.column_stack([1/a[obs], H[obs], theta_dot[obs]/2, phases(t[obs], 2)/2])
    infile = directory / f'{name}.bin'
    outfile = directory / f'{name}.txt'
    with infile.open('wb') as f:
        np.array([nstep, nq, len(obs)], dtype='<i4').tofile(f)
        np.array([dt, qmax, *masses], dtype='<f8').tofile(f)
        np.array(obs * substeps, dtype='<i4').tofile(f)
        np.array(observation, dtype='<f8').tofile(f)
        np.array(drive, dtype='<f8').tofile(f)
    subprocess.run([str(binary), str(infile), str(outfile)], check=True,
                   env=dict(os.environ, OMP_NUM_THREADS=str(threads)), timeout=1800)
    values = np.loadtxt(outfile).reshape(4, 2, nq, 4 + 2 * (len(obs)-1))
    assert np.all(np.isfinite(values))
    assert np.max(values[:, :, :, 3]) < 1e-8
    assert np.min(values[:, :, :, 4:]) >= 0 and np.max(values[:, :, :, 4:]) <= 1 + 1e-14
    return values


def controls(binary, directory, threads):
    t = np.linspace(0, 1.8, 61)
    masses = np.array([0.7, 1.2, 0.9, 0.4])
    hubble = 0.07
    a = np.exp(hubble * t)
    amplitude = np.array([1.7, 0.9])
    frequency = 2.3
    angles = (1 - np.cos(frequency * t[:, None])) * amplitude
    velocities = frequency * np.sin(frequency * t[:, None]) * amplitude
    obs = np.array([0, len(t)-1])
    zero = np.zeros_like(angles)
    static = evolve(binary, directory, 'static', t, np.ones_like(a), np.zeros_like(a),
                    zero, zero, obs, masses, (8, 2, 6), threads)
    assert np.max(static[..., 4:]) < 1e-23

    def reference(mass, parent, hq):
        phase = CubicHermiteSpline(t, angles[:, parent], velocities[:, parent])

        def endpoint(time):
            k = hq * math.exp(-hubble*time)
            z = k + float(phase(time, 1))/2
            w = math.hypot(mass, z)
            zp = -hubble*k + float(phase(time, 2))/2
            v = np.array([-mass/w, -mass*zp/(2*w**3), -z/w])
            return v/np.linalg.norm(v)

        def rhs(time, state):
            angle = phase(time)
            b = [mass * math.cos(angle), -mass * math.sin(angle), hq * math.exp(-hubble*time)]
            return 2 * np.cross(b, state)

        value = solve_ivp(rhs, (t[0], t[-1]), endpoint(t[0]), rtol=3e-13, atol=1e-14,
                          method='DOP853', max_step=0.02).y[:, -1]
        angle = phase(t[-1])
        rotated = np.array([math.cos(angle)*value[0] - math.sin(angle)*value[1],
                            math.sin(angle)*value[0] + math.cos(angle)*value[1], value[2]])
        return np.sum((rotated - endpoint(t[-1]))**2)/4

    targets = np.array([[[reference(mass, int(j > 0), h * (i + 0.5) * 0.75)
                          for i in range(8)] for h in (-1, 1)] for j, mass in enumerate(masses)])
    errors = {}
    for substeps in (1, 2, 4):
        actual = evolve(binary, directory, f'ode_control_{substeps}', t, a, np.full_like(a, hubble),
                        angles, velocities, obs, masses, (8, substeps, 6), threads)
        errors[str(substeps)] = float(np.max(abs(actual[..., 5] - targets)))
    assert errors['4'] < 2e-8 and errors['4'] < errors['1'] / 30
    return dict(maximum_static_occupation=float(np.max(static[..., 4:])),
                maximum_ODE_occupation_errors=errors)


def analyze(values, times, scale_factor, native, masses, settings):
    nq, substeps, qmax = settings
    q = values[0, 0, :, 2]
    edges = np.linspace(0, qmax, nq+1)
    weights = np.diff(edges**3) / (6 * math.pi**2)
    rows = []
    for index, (time, a) in enumerate(zip(times, scale_factor, strict=True)):
        entries = []
        for j, mass in enumerate(masses):
            p = q/a
            energy = np.hypot(p, mass)
            record = dict(species=f'N{j+1}', legacy_rho=float(native[index, j]))
            for order, name in enumerate(('instantaneous', 'dressed')):
                occup = np.sum(values[j, :, :, 4 + 2 * index + order], axis=0)
                n = (weights @ occup) / a**3
                rho = (weights @ (energy * occup)) / a**3
                p_density = (weights @ (p*p * occup / energy)) / (3 * a**3)
                ratio = rho / minimum_energy(n, mass)
                assert ratio >= 1 - 3/nq**2
                high = q >= 0.9 * qmax
                tail = (weights[high] @ (energy[high] * occup[high])) / a**3 / rho
                record[name] = dict(n=float(n), rho=float(rho), pressure=float(p_density),
                                    rho_over_minimum=float(ratio), high_tail_fraction=float(tail),
                                    rho_over_legacy=float(rho/native[index, j]))
            entries.append(record)
        rows.append(dict(t=float(time), a=float(a), species=entries))
    return dict(nq=nq, substeps=substeps, qmax=qmax, samples=rows,
                norm_error=float(np.max(values[..., 3])), maximum_occupation=float(np.max(values[..., 4:])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('tape', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = compile_engine(work)
    check = controls(binary, work, args.threads)
    print('controls', check, flush=True)

    meta, d, last, _ = load_source(args.tape)
    sl = slice(None, last+1)
    v = meta['P']['v']
    t, a, H = (d[k][sl] for k in ('t', 'a', 'H'))
    phases, velocities = [], []
    for parent, f in [('sigma', meta['channels'][0]['f']), ('chi', meta['channels'][1]['f'])]:
        field = d[f'{parent}_field'][sl]
        phases.append((field-field[0]) * v/f)
        velocities.append(d[f'{parent}_dot'][sl] * v/f)
    theta, theta_dot = np.array(phases).T, np.array(velocities).T
    obs = np.array([0, *(np.argmin(abs(t-x)) for x in (0.3, 0.5, 0.75, 1.5)), last])
    masses = np.array([s['M'] for s in meta['channels']])
    native = np.column_stack([d[f'rho_{j}'][obs[1:]] * v*v for j in range(1, 5)])
    results = {}
    for name in args.cases:
        values = evolve(binary, work, name, t, a, H, theta, theta_dot, obs, masses, CASES[name], args.threads)
        results[name] = analyze(values, t[obs[1:]], a[obs[1:]], native, masses, CASES[name])
        terminal = results[name]['samples'][-1]['species']
        print(name, 'dressed energies over native', [x['dressed']['rho_over_legacy'] for x in terminal], flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(scope='collisionless homogeneous projection of an actual lattice history, driven externally',
                   source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), controls=check,
                   max_a_theta_dot_over_two=list(np.max(abs(a[:, None]*theta_dot/2), axis=0)),
                   cases=results, accepted_cosmological_point=False)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
