#!/usr/bin/env python3
"""Evolve a finite H55 scalar lattice with local Majorana continua."""
import argparse
import hashlib
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import trapezoid

TRACK = Path(__file__).resolve().parents[1]
TRACKS = TRACK.parent
sys.path.insert(0, str(TRACKS / 'x5_coupled/code'))
from run import ARCHIVE, DT, h55_model, read_columns, write_input

END = 0.55
LATE = 2.9980937
# Grid size, radial bins, comoving cutoff, target step and endpoint.
CASES = {
    'homogeneous': (1, 360, 45000, DT, END),
    'reference': (3, 360, 45000, DT, END),
    'half_step': (3, 360, 45000, DT/2, END),
    'radial': (3, 720, 45000, DT, END),
    'cutoff': (3, 540, 67500, DT, END),
    'cutoff_half_step': (3, 540, 67500, DT/2, END),
    'spatial': (5, 360, 45000, DT, END),
    'early_boundary': (3, 360, 45000, DT, .35),
    'late': (3, 360, 45000, DT, LATE),
    'late_homogeneous': (1, 2880, 45000, DT, LATE),
    'scalar_coarse': (3, 360, 45000, DT, END),
    'scalar_fine': (17, 360, 45000, DT, END),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project(field, n, L):
    size = len(field)
    index = np.rint(n * np.fft.fftfreq(n)).astype(int)
    k = (2*np.pi/L) * np.rint(size*np.fft.fftfreq(size))
    k2 = k[:, None, None]**2 + k[None, :, None]**2 + k[None, None, :]**2
    select = np.ix_(index, index, index)
    full = np.fft.fftn(field) / size**3
    retained = full[select]
    data = np.fft.ifftn(retained).real * n**3
    norm = abs(full)**2
    grad = k2 * norm
    fractions = dict(variance=float(norm[select].sum()/norm.sum()),
                     gradient=float(grad[select].sum()/grad.sum()))
    return data, fractions


def compile_engine(work):
    binary = work / 'production'
    subprocess.run(['gcc', '-std=c11', '-O3', '-march=native', '-Wall', '-Wextra', '-Werror', '-fopenmp',
                    str(TRACK/'code/production.c'), '-lm', '-o', str(binary)], check=True)
    return binary


def summarize(data, spectra, model):
    times, rho, a = (data[k] for k in ('t','rho','a'))
    mass = [sp[1] for sp in model['species']]
    order = range(1, len(mass)+1)
    powers = sum((data[f'work_{j}'] for j in order), np.zeros_like(times))
    data['N_fraction'] = sum(data[f'rho_{j}'] for j in range(1, 5)) / rho
    data['N_count_comoving'] = a**3 * sum(data[f'n_{j}'] for j in range(1, 5))
    data['signed_work_over_density'] = powers/rho
    stress_jump = sum(data[f'p_free_{j}']-data[f'p_{j}'] for j in range(1, 5)) / rho
    data['dephased_pressure_jump'] = stress_jump
    data['scalar_subtraction_fraction'] = data['scalar_ct']/rho
    # Use local field and fermion stresses without allocating a subtraction twice.
    qenergy = a**3 * sum(data[f'rho_{j}'] for j in range(1, 5))
    qpressure = sum(data[f'p_{j}'] for j in range(1, 5))
    external = (a[1:]**3-a[:-1]**3) * 0.5 * (qpressure[1:]+qpressure[:-1])
    reciprocity = qenergy[1:]-qenergy[0]+np.cumsum(external)
    work = powers * a**3
    work_prefix = np.cumsum(np.diff(times) * (work[1:]+work[:-1]) / 2)
    err = np.max(abs(reciprocity-work_prefix))/(rho[0]*a[0]**3)
    idx = sorted({0, len(times)-1, int(np.argmax(data['N_fraction'])),
                  *(int(np.argmin(abs(times-t))) for t in (.12,.16,.20,.22,.23,.24,.25,.27,.30,.35,.40,.45,.50,1.,2.))})
    keys = ('t','a','H','rho','P','rho_s','rho_dr','G','sigma','chi','s','chi2','rms_sigma','rms_chi',
            'phase_sigma','phase_chi','N_fraction','N_count_comoving','signed_work_over_density',
            'dephased_pressure_jump','scalar_subtraction_fraction','drive_1','drive_2','drive_3','drive_4')
    details = []
    for i, M in zip(order, mass, strict=True):
        s = spectra['species']==i
        occupation = spectra['f_free'][s]
        weights = spectra['w'][s]
        E = np.hypot(spectra['q'][s]/a[-1], M)
        moments = np.array([weights @ occupation, (weights*E) @ occupation]) / a[-1]**3
        expected = [data[f'n_{i}'][-1], data[f'rho_{i}'][-1]]
        contribution = weights*E*occupation
        cutoff50 = float(contribution[spectra['q'][s] >= 0.5*model['qmax']].sum()/contribution.sum())
        details.append(dict(species=i, number=float(moments[0]), energy=float(moments[1]),
                            native_moment_relative_error=float(max(abs(moments/expected-1))),
                            terminal_energy_fraction_above_half_cutoff=cutoff50))
    return dict(
        terminal={k:float(v[-1]) for k,v in data.items()},
        samples=[{k:float(data[k][i]) for k in keys} for i in idx],
        max_total_energy_ledger_error=float(max(abs(data['residual']))),
        max_reciprocal_work_ledger_error_at_output_cadence=float(err),
        max_reciprocal_work_ledger_error=float(max(abs(data['Q_ledger_error']))),
        integrated_comoving_scalar_work=float(data['Q_scalar_work'][-1]),
        sampled_quadrature_scalar_work=float(trapezoid(work,times)),
        max_quantum_work_over_density=float(max(abs(data['signed_work_over_density']))),
        max_mass_phase_deviation={k:float(max(data[k])) for k in ('phase_sigma','phase_chi')},
        occupation_range=[float(min(data['nmin'])),float(max(data['nmax']))],
        max_bloch_norm_squared_error=float(max(data['norm_error'])),
        moment_projection=details)


def run_case(name, work, meta, fields, binary, threads):
    n, nq, qmax, dt, end = CASES[name]
    coupled = not name.startswith('scalar_')
    model = h55_model(meta, (nq, qmax, dt, True, False, coupled))
    model['t1'] = end
    model['nsteps'] = math.ceil((end-model['t0'])/dt)
    actual_step = (end-model['t0'])/model['nsteps']
    model['stride'] = max(1, round(.0004/actual_step))
    folder = work / name
    folder.mkdir(parents=True, exist_ok=True)
    inp, seed, history, spectra = (folder/f'{f}.txt' for f in ('input','seed','history','spectra'))
    columns, fractions = [], {}
    for field in ('sigma','chi','spectator_a','dsigma','dchi','dspectator_a'):
        values, fractions[field] = project(fields[field], n, meta['run']['L'])
        columns.append(values.ravel())
    write_input(inp, model)
    np.savetxt(seed,np.column_stack(columns),fmt='%.17g',
               header=f"{n} {meta['run']['L']:.17g} {model['H0']:.17g}",comments='')
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_MAX_ACTIVE_LEVELS='1')
    print(f'{name} starting grid {n}^3 and {nq} radial bins',flush=True)
    process = subprocess.run([str(binary),str(inp),str(seed),str(history),str(spectra)],
                             env=env,capture_output=True,text=True)
    if process.returncode:
        if name != 'cutoff' or 'Midpoint scalar/fermion update did not converge' not in process.stderr:
            raise RuntimeError(f'{name}: {process.stderr.strip()}')
        result = dict(grid=n,radial_bins=nq,comoving_cutoff=qmax,timestep=actual_step,
                      completed=False,last_recorded_time=float(read_columns(history)['t'][-1]),
                      reason=process.stderr.strip(),input_sha256=sha(inp),seed_sha256=sha(seed))
        (folder/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
        print(f'{name} retained midpoint failure',flush=True)
        return
    data = read_columns(history)
    if not math.isclose(data['t'][-1],end,rel_tol=0,abs_tol=1e-13):
        raise ValueError(f'{name} did not reach endpoint')
    result = summarize(data, read_columns(spectra) if coupled else {}, model)
    result.update(grid=n, radial_bins=nq, comoving_cutoff=qmax, timestep=actual_step, completed=True,
                  seed_fractions=fractions, input_sha256=sha(inp), seed_sha256=sha(seed))
    (folder/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'{name} complete. N fraction {result["terminal"]["N_fraction"]:.6f}, '
          f'energy error {result["max_total_energy_ledger_error"]:.3g}',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('--only',nargs='+',choices=CASES,default=list(CASES))
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--workers',type=int,default=1)
    args = parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0, str(TRACKS / 'x10_handoff/code'))
    from h55_seed import archive, reconstruct

    meta = json.loads(ARCHIVE.read_text())
    run, _, _ = archive()
    fields = reconstruct(run)
    binary = compile_engine(args.work)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(run_case,name,args.work,meta,fields,binary,args.threads) for name in args.only]
        for job in jobs:
            job.result()


if __name__ == '__main__':
    main()
