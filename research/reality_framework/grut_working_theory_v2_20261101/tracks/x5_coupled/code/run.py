#!/usr/bin/env python3
"""Replay finite-regulator coupled H55 histories and collect their diagnostics."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np

TRACK = Path(__file__).resolve().parents[1]
TRACKS = TRACK.parent
SOURCE = Path(__file__).with_name('coupled.c')
ARCHIVE = TRACKS / 'x2_strip/results/lower_H_lattice.json'
DT = 0.0000328125
CASES = {
    'reference': (720, 45000, DT, True, False, True),
    'half_step': (720, 45000, DT / 2, True, False, True),
    'double_bins': (1440, 45000, DT, True, False, True),
    'wide_cutoff': (1080, 67500, DT, True, False, True),
    'scalar_control': (720, 45000, DT, True, False, False),
    'decay_control': (720, 45000, DT, True, True, True),
    'unsubtracted_control': (720, 45000, DT, False, False, True),
}
BASE_CASES = tuple(CASES)
CASES.update({
    'four_bins': (2880, 45000, DT, True, False, True),
    'eight_bins': (5760, 45000, DT, True, False, True),
    'wide_four_bins': (4320, 67500, DT, True, False, True),
    'four_bins_half_step': (2880, 45000, DT / 2, True, False, True),
})
MODEL_KEYS = 'mp lambda v g2 sc c1 m2 mu2 cb offset gs2 gamma_s'.split()
SAVE_TIMES = (0.0, 0.1, 0.2, 0.225, 0.24, 0.25, 0.275, 0.3, 0.4, 0.5,
              0.6, 0.75, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 2.9)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compile_engine(work):
    binary = work / 'coupled'
    subprocess.run(['gcc', '-std=c11', '-O3', '-march=native', '-Wall', '-Wextra', '-Werror',
                    '-fopenmp', '-o', str(binary), str(SOURCE), '-lm'], check=True)
    return binary


def write_input(path, model):
    header = [int(model[k]) for k in ('expand', 'subtract', 'dressed_start', 'nq')]
    header += [len(model['species']), model['nsteps'], model['stride']]
    rows = [header, [model[k] for k in ('t0', 't1', 'qmax')],
            [model[k] for k in MODEL_KEYS], [model[k] for k in ('a0', 'H0', 'rho_sm', 'rho_dr')],
            model['fields'], model['velocities'], *model['species']]
    path.write_text(''.join(' '.join(format(v, '.17g') for v in row) + '\n' for row in rows))


def read_columns(path):
    with path.open() as f:
        keys = f.readline().lstrip('# ').split()
        data = np.loadtxt(f, ndmin=2)
    if data.shape[1] != len(keys) or not np.all(np.isfinite(data)):
        raise ValueError(f'Bad history or spectrum in {path}')
    return dict(zip(keys, data.T, strict=True))


def run_engine(binary, work, name, model, threads=1):
    folder = work / name
    folder.mkdir(parents=True, exist_ok=True)
    infile, history, spectra = (folder / item for item in ('input.txt', 'history.txt', 'spectra.txt'))
    write_input(infile, model)
    start = time.monotonic()
    proc = subprocess.run([str(binary), str(infile), str(history), str(spectra)], check=False,
                          capture_output=True, text=True, timeout=1800,
                          env=dict(os.environ, OMP_NUM_THREADS=str(threads)))
    if proc.returncode != 0:
        raise RuntimeError(f'{name}: {proc.stderr.strip()}')
    d = read_columns(history)
    if abs(d['t'][-1] - model['t1']) > 1e-12:
        raise ValueError(f'{name}: truncated evolution')
    return d, read_columns(spectra) if model['species'] else {}, {
        'input_sha256': digest(infile), 'history_sha256': digest(history), 'spectra_sha256': digest(spectra),
        'runtime_seconds': time.monotonic() - start,
    }


def h55_model(meta, settings):
    nq, qmax, dt_goal, subtract, decay, coupled = settings
    model = meta['run']['P']
    start = meta['run']['start']
    t1 = meta['termination']['last_diag_t']
    nsteps = math.ceil((t1 - start['t0']) / dt_goal)
    step = (t1 - start['t0']) / nsteps
    species = [(s['w'], s['M'], s['f'], s['Grest'] if decay else 0.0)
               for s in meta['run']['channels']] if coupled else []
    return dict(
        expand=True, subtract=subtract, dressed_start=True, nq=nq, nsteps=nsteps,
        stride=max(1, round(0.003 / step)), t0=start['t0'], t1=t1, qmax=qmax,
        mp=model['mP'], **{'lambda': model['lam']}, v=model['v'], g2=model['g']**2,
        sc=model['sc'], c1=model['c1'], m2=model['m2'], mu2=model['mu2'], cb=model['cb'],
        offset=model['Vc'], gs2=model['g_s']**2, gamma_s=1.3e-6 * model['H_GeV'] / 1.3e11,
        fields=[start['s'], start['x'], model['s_star']],
        velocities=[start['sd'], start['xd'], 0.0], a0=start['a'], H0=start['H'],
        rho_sm=0.0, rho_dr=0.0, species=species,
    )


def archived_history():
    path = ARCHIVE.with_name('lower_H_diag.txt.gz')
    with gzip.open(path, 'rt') as f:
        labels = f.readline().lstrip('# ').split()
        matrix = np.loadtxt(f)
    return dict(zip(labels, matrix.T, strict=True))


def average_last(d, key, duration=0.10):
    t = d['t']
    left = t[-1] - duration
    end = t > left
    return float(np.trapezoid(np.r_[np.interp(left, t, d[key]), d[key][end]], np.r_[left, t[end]]) / duration)


def summarize(case, model, data, spec, info, old):
    a, rho = data['a'], data['rho']
    t, h = data['t'], data['H']
    bare_fermion = sum(data[f'rho_{i}'] for i in range(1, 5))
    dressed_fermion = sum(data[f'dressed_rho_{i}'] for i in range(1, 5))
    net = bare_fermion - data['scalar_ct'] - data['gravity_ct']
    final = {key: float(value[-1]) for key, value in data.items()}
    turns = np.flatnonzero((data['chi'][1:] > model['v']) &
                           (data['dchi'][:-1] > 0) & (data['dchi'][1:] <= 0)) + 1
    i_turn = int(turns[0]) if len(turns) else None
    records = []
    for i, (parent, mass, f, width) in enumerate(model['species'], 1):
        band = spec['species'] == i
        p = spec['q'][band] / a[-1]
        weights = spec['w'][band] / a[-1]**3
        e = np.hypot(p, mass)
        f0, f1 = spec['f_free'][band], spec['f_dressed'][band]
        tail = spec['q'][band] > 0.9 * model['qmax']
        assert abs(weights @ (f0 * e) / final[f'rho_{i}'] - 1) < 5e-10
        assert abs(weights @ f0 / final[f'n_{i}'] - 1) < 5e-10
        records.append(dict(
            species=i, mass=mass, f=f, vacuum_width=width,
            free_energy=final[f'rho_{i}'], free_number=final[f'n_{i}'],
            dressed_energy=final[f'dressed_rho_{i}'], dressed_number=final[f'dressed_n_{i}'],
            mean_free_energy=average_last(data, f'rho_{i}'),
            mean_free_number=average_last(data, f'n_{i}'),
            mean_dressed_energy=average_last(data, f'dressed_rho_{i}'),
            mean_dressed_number=average_last(data, f'dressed_n_{i}'),
            free_energy_over_legacy=final[f'rho_{i}'] / (old[f'rN{i}'][-1] * model['v']**2),
            tail_free_energy_fraction=float(weights[tail] @ (f0[tail]*e[tail]) / final[f'rho_{i}']),
            tail_dressed_energy_fraction=float(weights[tail] @ (f1[tail]*e[tail]) / final[f'dressed_rho_{i}']),
            free_occupation_range=[float(min(f0)), float(max(f0))],
            dressed_occupation_range=[float(min(f1)), float(max(f1))],
            helicity_number_asymmetry=float((weights @ (f0 * spec['h'][band])) / (weights @ f0)),
        ))
    indices = sorted({0, len(t)-1, *(int(np.argmin(abs(t - when))) for when in SAVE_TIMES)})
    history = [dict(
        t=float(t[i]), a=float(a[i]), H=float(h[i]), sigma_over_v=float(data['sigma'][i]/model['v']),
        chi_over_v=float(data['chi'][i]/model['v']), s=float(data['s'][i]),
        rho=float(rho[i]), f_fermion=float(bare_fermion[i]/rho[i]),
        f_fermion_subtracted=float(net[i]/rho[i]), f_s=float(data['rho_s'][i]/rho[i]),
        f_sm=float(data['rho_sm'][i]/rho[i]), f_dr=float(data['rho_dr'][i]/rho[i]),
        residual=float(data['residual'][i]),
    ) for i in indices]
    return dict(
        case=case, model=model, terminal=final, species=records, history=history,
        first_turn_t=float(t[i_turn]) if i_turn is not None else None,
        max_energy_budget_error=float(max(abs(data['residual']))),
        max_norm=float(max(data['norm_max'])), occupation_range=[float(min(data['f_min'])), float(max(data['f_max']))],
        terminal_fermion_fraction=float(bare_fermion[-1]/rho[-1]),
        terminal_subtracted_fermion_fraction=float(net[-1]/rho[-1]),
        terminal_dressed_fermion_fraction=float(dressed_fermion[-1]/rho[-1]),
        terminal_projector_energy_gap_fraction=float((dressed_fermion[-1]-bare_fermion[-1])/rho[-1]),
        terminal_subtracted_to_dressed_gap_fraction=float((dressed_fermion[-1]-net[-1])/rho[-1]),
        t_peak_free_fraction=float(t[np.argmax(bare_fermion/rho)]),
        peak_free_fermion_fraction=float(max(bare_fermion/rho)),
        output_records=len(t), scientific_acceptance=False, **info,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=BASE_CASES)
    parser.add_argument('--workers', type=int, default=5)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = compile_engine(work)
    archived = json.loads(ARCHIVE.read_text())
    old = archived_history()

    def run(name):
        model = h55_model(archived, CASES[name])
        data, spec, info = run_engine(binary, work, name, model, args.threads)
        result = summarize(name, model, data, spec, info, old)
        (work / name / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        print(name, f'fermion={result["terminal_fermion_fraction"]:.7g}',
              f'budget={result["max_energy_budget_error"]:.3g}', f'{info["runtime_seconds"]:.1f}s',
              flush=True)
        return result

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        outcomes = list(pool.map(run, args.cases))
    output = dict(parent_archive_sha256=digest(ARCHIVE), engine_sha256=digest(SOURCE),
                  cases=outcomes, scope='finite-regulator homogeneous H55 with mutual backreaction',
                  scientific_acceptance=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
