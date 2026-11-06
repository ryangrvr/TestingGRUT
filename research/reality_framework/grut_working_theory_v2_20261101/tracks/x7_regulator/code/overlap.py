#!/usr/bin/env python3
"""Replay the archived H55 state and measure projection gaps past its old boundary."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq

TRACKS = Path(__file__).resolve().parents[2]
X5 = TRACKS / 'x5_coupled'
sys.path.insert(0, str(X5 / 'code'))
from run import read_columns, write_input

DURATION = 0.12
WINDOW = 0.02
TOLERANCE = 1e-6
GAPS = ('dressed_gap', 'pressure_gap', 'scalar_work')
CASES = {'coherent': (1, False), 'coherent_half': (2, False),
         'vacuum_decay': (1, True), 'vacuum_decay_half': (2, True)}


def initial_data():
    results = X5 / 'results'
    refined = json.loads((results / 'h55_refined.json').read_text())
    case = next(row for row in refined['cases'] if row['case'] == 'eight_bins')
    audit = json.loads((results / 'diagnosis.json').read_text())
    reference = next(row for row in audit['handoffs'] if row['case'] == 'eight_bins')
    state = results / 'h55_terminal_state.npz'
    sha = hashlib.sha256(state.read_bytes()).hexdigest()
    if sha != audit['published_terminal_state_sha256']:
        raise ValueError('Checkpoint digest differs from X5 record')
    with np.load(state) as f:
        data = dict(f)
    model = case['model']
    n = model['nq']
    if data['q'].shape != (n,) or data['w'].shape != (n,) or data['bloch'].shape != (4, 2, n, 3):
        raise ValueError('Checkpoint array dimensions differ from model')
    if any(not np.all(np.isfinite(a)) for a in data.values()):
        raise ValueError('Nonfinite checkpoint')
    x2 = json.loads((TRACKS / 'x2_strip/results/lower_H_lattice.json').read_text())
    widths = [row['Grest'] for row in x2['run']['channels']]
    return model, case['terminal'], reference, data, sha, widths, audit


def frozen_tail(model, terminal, data, widths):
    """Exact worst-phase envelope in a fixed free basis with vacuum amplitude damping."""
    a = terminal['a']
    k = np.array([-1., 1.])[:,None]*data['q']/a
    mass = np.array([row[1] for row in model['species']])[:,None,None]
    energy = np.hypot(k,mass)
    direction = np.stack((mass/energy,np.zeros_like(energy),k/energy),axis=-1)
    r = data['bloch']
    parallel = np.sum(direction*r,axis=-1)
    occupation = (1+parallel)/2
    transverse = r-direction*parallel[...,None]
    z_amplitude = np.hypot(transverse[...,2],np.cross(direction,transverse)[...,2])
    weight = data['w']/a**3
    pressure_mode = weight*abs(k)*z_amplitude/(6*terminal['rho'])
    gas_mode = weight*energy*occupation
    rate = np.array(widths)[:,None,None]*mass/energy
    parent0 = gas_mode.sum()

    def surviving(t):
        return np.sum(gas_mode*np.exp(-rate*t))/parent0

    def envelope(t):
        return np.sum(pressure_mode*np.exp(-rate*t/2))

    def root(f):
        upper = 1.
        while f(upper)>0:
            upper *= 2
        return brentq(f,0,upper)

    half = root(lambda t: surviving(t)-0.5)
    pressure_time = root(lambda t: envelope(t)-TOLERANCE)
    p_free = np.sum(weight*k*k*occupation/(3*energy),axis=(1,2))
    p_full = np.sum(weight*k*(r[...,2]+k/energy)/6,axis=(1,2))
    return dict(
        free_pressure_gap=float(np.sum(p_free-p_full)/terminal['rho']),
        pressure_envelope_initial=float(envelope(0)),
        pressure_envelope_at_half_parent_energy=float(envelope(half)),
        pressure_envelope_at_live_window_end=float(envelope(DURATION)),
        time_to_half_parent_energy=float(half),
        time_to_phase_independent_pressure_tolerance=float(pressure_time),
        parent_fraction_at_pressure_tolerance=float(surviving(pressure_time)),
        hubble_time_initial=1/terminal['H'],
        actual_widths=list(widths),
        max_gamma_over_H_initial=float(np.max(rate)/terminal['H']),
        energy_half_before_phase_independent_pressure_gate=bool(half<pressure_time),
    )


def moment_check(model,terminal,data):
    """Reconstruct the archive in free-basis moments and test phase reflection."""
    a,H,rho = (terminal[k] for k in ('a','H','rho'))
    k = np.array([-1.,1.])[:,None]*data['q']/a
    m = np.array([row[1] for row in model['species']])[:,None,None]
    fscale = np.array([row[2] for row in model['species']])[:,None,None]
    v = np.array([terminal[('dsigma','dchi','ds')[row[0]]] for row in model['species']])[:,None,None]
    E = np.hypot(k,m)
    rx,ry,rz = np.moveaxis(data['bloch'],-1,0)
    parallel = (m*rx+k*rz)/E
    u = (k*rx-m*rz)/E
    f = (1+parallel)/2
    w = data['w']/a**3
    shift = v/(2*fscale)
    basis_rate = H*m*k/E**2
    h_perp = -m*shift/E
    occupation_derivative = h_perp*ry+basis_rate*u/2
    pressure_gas = np.sum(w*k*k*f/(3*E),axis=(1,2))
    pressure_gap = np.sum(w*k*m*u/(6*E),axis=(1,2))
    pressure = pressure_gas-pressure_gap
    work = np.sum(-w*m*ry*shift,axis=(1,2))
    derivative = np.sum(w*(E*occupation_derivative-H*k*k*f/E),axis=(1,2))
    balance = work-3*H*pressure
    rho_F = np.sum(w*E*f,axis=(1,2))
    columns = dict(energy=rho_F,pressure=pressure,work=work)
    want = dict(energy=[terminal[f'rho_{i}'] for i in range(1,5)],
                pressure=[terminal[f'p_{i}'] for i in range(1,5)],
                work=[terminal[f'J_{i}']*v[i-1,0,0] for i in range(1,5)])
    return dict(
        archive_reconstruction_errors={key:float(np.max(abs(value-want[key]))/rho)
                                       for key,value in columns.items()},
        max_continuity_identity_error=float(np.max(abs(derivative-balance))/rho),
        max_bloch_ball_excess=float(np.max(u*u+ry*ry-4*f*(1-f))),
        initial_pressure_gap_by_species_over_rho=(pressure_gap/rho).tolist(),
        initial_work_by_species_over_rho=(work/rho).tolist(),
        phase_reflected_work_over_rho=float(-work.sum()/rho),
        phase_reflection_derivative_separation_over_rho=float(2*abs(work.sum())/rho),
    )


def window_results(t, score):
    """A sampled passing window must contain no violating sample."""
    good = score<TOLERANCE
    bounds = np.flatnonzero(np.diff(np.r_[False,good,False]))
    starts = [(t[lo],t[hi-1]) for lo,hi in bounds.reshape(-1,2)]
    best = max((hi-lo for lo,hi in starts),default=0.)
    # Overlapping failing windows are also summarized by their lowest maximum gap.
    maxima = []
    for lo in range(len(t)):
        hi = int(np.searchsorted(t,t[lo]+WINDOW-1e-13))
        if hi<len(t):
            maxima.append((float(max(score[lo:hi+1])),float(t[lo])))
    min_max, best_start = min(maxima)
    return dict(sampled_passing_points=int(good.sum()), longest_passing_window=float(best),
                best_window_start=best_start, best_window_max_gap=min_max, gate_pass=bool(best>=WINDOW-1e-13))


def measure(d, projection, parents):
    rho = d['rho']
    raw_p = np.array([d[f'p_{i}'] for i in range(1,5)]).T
    work = np.column_stack([d[f'J_{i}']*d[('dsigma','dchi','ds')[j]]
                            for i,j in enumerate(parents,1)])/rho[:,None]
    coherent_pressure = (projection[:,:4]-raw_p)/rho[:,None]
    raw_e = sum(d[f'rho_{i}'] for i in range(1,5))
    dressed_e = sum(d[f'dressed_rho_{i}'] for i in range(1,5))
    pressure_gap = coherent_pressure.sum(axis=1)
    scalar_work = work.sum(axis=1)
    return dict(t=d['t'], dressed_gap=(dressed_e-raw_e+d['scalar_ct']+d['gravity_ct'])/rho,
                pressure_gap=pressure_gap, scalar_work=scalar_work,
                scalar_work_matched=projection[:,4],
                pressure_gap_matched=pressure_gap+projection[:,5],
                free_energy_gap_matched=d['scalar_ct']/rho,
                derivative_jump=-scalar_work-3*d['H']*pressure_gap,
                max_species_abs_work=np.max(abs(work),axis=1),
                max_species_abs_pressure_gap=np.max(abs(coherent_pressure),axis=1),
                energy_residual=d['residual'], fermion_fraction=raw_e/rho,
                radiation_sm_fraction=d['rho_sm']/rho, radiation_dark_fraction=d['rho_dr']/rho)


def summarize(d, trace):
    score = np.max(abs(np.column_stack([trace[key] for key in GAPS])),axis=1)
    t = trace['t']-trace['t'][0]
    best = int(np.argmin(score))
    matched = np.max(abs(np.column_stack([trace[key] for key in
                     ('free_energy_gap_matched','pressure_gap_matched','scalar_work_matched')])),axis=1)
    return dict(
        samples=len(t), sample_spacing=float(max(np.diff(t))),
        initial={key:float(value[0]) for key,value in trace.items()},
        final={key:float(value[-1]) for key,value in trace.items()},
        min_abs_gaps={key:float(min(abs(trace[key]))) for key in GAPS},
        max_abs_gaps={key:float(max(abs(trace[key]))) for key in GAPS},
        best_simultaneous_gap=float(score[best]), time_at_best_gap=float(t[best]),
        budget_max=float(max(abs(d['residual']))), pauli_min=float(min(d['f_min'])),
        pauli_max=float(max(d['f_max'])), bloch_norm_max=float(max(d['norm_max'])),
        max_iterations=int(max(d['max_iterations'])), **window_results(t,score),
        canonical_scalar_free_projection=dict(best_simultaneous_gap=float(min(matched)),
                                               **window_results(t,matched)),
    )


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    model,terminal,reference,state,sha,widths,audit=initial_data()
    checkpoint=work/'checkpoint.bin'
    with checkpoint.open('wb') as f:
        for key in ('q','w','bloch'):
            np.asarray(state[key],dtype='=f8',order='C').tofile(f)
    source=Path(__file__).with_name('overlap.c')
    engine=work/'overlap'
    subprocess.run(['gcc','-std=c11','-O3','-march=native','-fopenmp','-Wall','-Wextra','-Werror',
                    str(source),'-lm','-o',str(engine)],check=True)
    base_steps=math.ceil(DURATION*model['nsteps']/(model['t1']-model['t0']))

    def run(name):
        multiplier,decay=CASES[name]
        species=[[*row[:3],width if decay else 0.] for row,width in zip(model['species'],widths,strict=True)]
        case=dict(model,species=species,t0=terminal['t'],t1=terminal['t']+DURATION,
                  nsteps=multiplier*base_steps,stride=4*multiplier,a0=terminal['a'],H0=terminal['H'],
                  fields=[terminal[k] for k in ('sigma','chi','s')],
                  velocities=[terminal[k] for k in ('dsigma','dchi','ds')],
                  rho_sm=terminal['rho_sm'],rho_dr=terminal['rho_dr'])
        infile,history,projection=[work/f'{name}.{ext}' for ext in ('input','history','projection')]
        write_input(infile,case)
        subprocess.run([str(engine),str(infile),str(checkpoint),str(history),str(projection)],check=True,
                       env=dict(os.environ,OMP_NUM_THREADS='1'))
        d=read_columns(history)
        if abs(d['t'][-1]-case['t1'])>1e-12:
            raise ValueError(f'Truncated replay: {name}')
        trace=measure(d,np.loadtxt(projection,ndmin=2),[row[0] for row in species])
        np.savez_compressed(work/f'{name}.npz',**trace)
        result=summarize(d,trace)
        result['dt']=DURATION/case['nsteps']
        result['restart_max_relative_difference']=max(
            abs(d[k][0]-terminal[k])/max(1,abs(terminal[k]))
            for k in ('a','H','rho','pressure','sigma','chi','s','dsigma','dchi','ds',
                      'rho_1','rho_2','rho_3','rho_4','rho_dr','scalar_ct','gravity_ct'))
        print(name,'budget',result['budget_max'],'best window gap',result['best_window_max_gap'],flush=True)
        return name,(result,trace)

    with ThreadPoolExecutor(max_workers=4) as pool:
        histories=dict(pool.map(run,CASES))
    differences={}
    for key in ('coherent','vacuum_decay'):
        d=histories[key][1]
        ref=histories[f'{key}_half'][1]
        differences[key]={field:float(max(abs(d[field]-ref[field]))) for field in
                          (*GAPS,'fermion_fraction','radiation_sm_fraction')}
    result=dict(
        scope='archived homogeneous H55 overlap at the X5 momentum grid and matching convention',
        duration=DURATION, window=WINDOW, tolerance=TOLERANCE, checkpoint_sha256=sha,
        checkpoint_case='eight_bins', spectator_decay_retained=True,
        species=[dict(parent=row[0],mass=row[1],phase_scale=row[2],width=w)
                 for row,w in zip(model['species'],widths,strict=True)],
        inherited_terminal_refinement_pass=audit['numerical_terminal_pass'],
        reference_gaps=dict(dressed_gap=reference['replacing_subtracted_stress_with_dressed_energy_gap'],
                            pressure_gap=reference['total_pressure_change_if_dephased_over_rho'],
                            scalar_work=reference['total_coherent_fermion_power_over_rho']),
        cases={name:case for name,(case,_) in histories.items()},
        timestep_max_abs_differences=differences,
        frozen_tail=frozen_tail(model,terminal,state,widths), scientific_acceptance=False,
        archive_moment_check=moment_check(model,terminal,state),
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
