#!/usr/bin/env python3
"""Measure numerical convergence and the stress change at a particle handoff."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from run import ARCHIVE, digest, read_columns

MOMENTS = ('free_energy', 'free_number', 'dressed_energy', 'dressed_number')
PAIRS = (
    ('reference', 'half_step'), ('reference', 'double_bins'), ('reference', 'wide_cutoff'),
    ('double_bins', 'four_bins'), ('four_bins', 'eight_bins'),
    ('four_bins', 'wide_four_bins'), ('four_bins', 'four_bins_half_step'),
)
TARGET = 0.01


def compare(runs):
    comparisons = []
    for a, b in PAIRS:
        first = runs[a]['species']
        second = runs[b]['species']
        errors = {key: [y[key] / x[key] - 1 for x, y in zip(first, second, strict=True)]
                  for key in MOMENTS + tuple('mean_' + key for key in MOMENTS)}
        terminal = max(abs(x) for key in MOMENTS for x in errors[key])
        window = max(abs(x) for key in MOMENTS for x in errors['mean_' + key])
        comparisons.append(dict(base=a, control=b, fractional_changes=errors,
                                terminal_maximum=terminal, final_window_maximum=window,
                                terminal_pass=terminal < TARGET, final_window_pass=window < TARGET))
    return comparisons


def spectrum_data(run, path):
    d = read_columns(path)
    nq = run['model']['nq']
    state = np.stack([d[f'r_{part}'] for part in ('x', 'y', 'z')], axis=-1)
    state = state.reshape(4, 2, nq, 3)
    q = d['q'][:nq]
    w = d['w'][:nq]
    expected = np.repeat(np.arange(1, 5), 2 * nq)
    assert np.array_equal(d['species'], expected)
    assert np.array_equal(d['h'].reshape(4, 2, nq)[:, 0], np.full((4, nq), -1))
    assert np.array_equal(d['h'].reshape(4, 2, nq)[:, 1], np.full((4, nq), 1))
    return q, w, state


def handoff(run, q, w, state, widths):
    """Dephasing along the free-mass axis leaves energy but removes work."""
    terminal = run['terminal']
    a, rho = terminal['a'], terminal['rho']
    k = np.stack([-q, q]) / a
    weights = w / a**3
    rows = []
    for j, (parent, mass, scale, _) in enumerate(run['model']['species'], 1):
        r = state[j-1]
        energy = np.hypot(k, mass)
        parallel = (mass*r[...,0]+k*r[...,2])/energy
        f = (1+parallel)/2
        free_energy = float(np.sum((f*energy) @ weights))
        number = float(np.sum(f @ weights))
        coherent_pressure = float(np.sum((k*(r[...,2]+k/energy)/6) @ weights))
        diagonal_pressure = float(np.sum((k*k*f/(3*energy)) @ weights))
        force = -float(np.sum((mass*r[...,1]/(2*scale)) @ weights))
        assert abs(free_energy / terminal[f'rho_{j}'] - 1) < 5e-10
        assert abs(number / terminal[f'n_{j}'] - 1) < 5e-10
        assert abs(coherent_pressure / terminal[f'p_{j}'] - 1) < 5e-10
        assert abs(force - terminal[f'J_{j}']) < 5e-10 * max(1, abs(force))
        vel = terminal[('dsigma','dchi')[parent]]
        kf = (3 * math.pi**2 * number)**(1/3)
        minimum = (kf*math.hypot(kf,mass)*(2*kf*kf+mass*mass)
                   - mass**4*math.asinh(kf/mass))/(8*math.pi**2)
        assert free_energy >= minimum
        rows.append(dict(species=j, free_number=number, free_energy=free_energy,
                         energy_over_Pauli_minimum=free_energy/minimum,
                         coherent_pressure=coherent_pressure, diagonal_pressure=diagonal_pressure,
                         pressure_change_if_dephased_over_rho=(diagonal_pressure-coherent_pressure)/rho,
                         coherent_fermion_power_over_rho=vel*force/rho,
                         absolute_phase_rate_over_mass=abs(vel/(2*scale*mass)),
                         vacuum_coherence_damping_time_lower_bound=2/widths[j-1],
                         occupation_range=[float(np.min(f)),float(np.max(f))]))
    return dict(
        case=run['case'], t=terminal['t'], rho=rho, species=rows,
        total_pressure_change_if_dephased_over_rho=sum(s['pressure_change_if_dephased_over_rho'] for s in rows),
        total_coherent_fermion_power_over_rho=sum(s['coherent_fermion_power_over_rho'] for s in rows),
        replacing_free_with_dressed_energy_gap=run['terminal_projector_energy_gap_fraction'],
        replacing_subtracted_stress_with_dressed_energy_gap=run['terminal_subtracted_to_dressed_gap_fraction'],
    )


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline',type=Path)
    parser.add_argument('refinements',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--state',type=Path)
    args=parser.parse_args()
    runs={}
    for work in (args.baseline,args.refinements):
        for path in sorted(work.glob('*/summary.json')):
            row=json.loads(path.read_text())
            runs[row['case']]=row
    widths=[s['Grest'] for s in json.loads(ARCHIVE.read_text())['run']['channels']]
    ledger=[]
    for name in ('reference','four_bins','eight_bins','wide_four_bins','four_bins_half_step','decay_control'):
        path=(args.baseline if name in ('reference','decay_control') else args.refinements)/name/'spectra.txt'
        assert digest(path)==runs[name]['spectra_sha256']
        q,w,state=spectrum_data(runs[name],path)
        ledger.append(handoff(runs[name],q,w,state,widths))
        if name=='eight_bins' and args.state is not None:
            np.savez_compressed(args.state,q=q,w=w,bloch=state)
    c=compare(runs)
    scientific_pass=False
    payload=dict(
        comparisons=c, numerical_terminal_pass=all(r['terminal_pass'] for r in c[-3:]),
        numerical_window_pass=all(r['final_window_pass'] for r in c[-3:]),
        handoffs=ledger, published_terminal_state_case='eight_bins' if args.state else None,
        published_terminal_state_sha256=digest(args.state) if args.state else None,
        scientific_acceptance=scientific_pass,
        scope='coherent-to-kinetic stress diagnosis at the fixed archived H55 time',
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print('Fine terminal and time-averaged gates:',
          payload['numerical_terminal_pass'],payload['numerical_window_pass'],flush=True)


if __name__=='__main__':
    main()
