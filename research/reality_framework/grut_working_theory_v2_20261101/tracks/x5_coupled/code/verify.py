#!/usr/bin/env python3
"""Verify the published numerical evidence without rerunning expensive evolutions."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from diagnose import compare, handoff
from run import ARCHIVE, SOURCE, TRACK, TRACKS, digest

RES = TRACK / 'results'
REPO = TRACK.parents[4]
PARENT = '83ce5341a3819585d250f6a6bd6f3da6390496f0'


def read(filename):
    return json.loads((RES / filename).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    original, refined = read('h55.json'), read('h55_refined.json')
    runs = {run['case']: run for data in (original, refined) for run in data['cases']}
    ref, finest = runs['reference'], runs['eight_bins']
    diag = read('diagnosis.json')
    ctl = read('controls.json')
    spatial = read('spatial.json')
    phase = read('phase_trace.json')

    check('all eleven coupled and scalar cases', len(runs) == 11)
    check('compiled engine digest', original['engine_sha256'] == refined['engine_sha256'] == digest(SOURCE))
    check('initial H55 input digest', original['parent_archive_sha256'] == refined['parent_archive_sha256'] == digest(ARCHIVE))
    old = json.loads(ARCHIVE.read_text())
    check('archived interface time', all(
        abs(run['terminal']['t']-old['termination']['last_diag_t']) < 1e-12 for run in runs.values()))
    check('total pressure-work budget', max(r['max_energy_budget_error'] for r in runs.values()) < 1e-4)
    check('physical two-helicity states', all(
        r['max_norm'] <= 1 + 1e-9 and r['occupation_range'][0] >= -1e-9 and r['occupation_range'][1] <= 1 + 1e-9
        for r in runs.values()))
    check('fermion feedback changes expansion', ref['terminal']['a'] < 0.95 * runs['scalar_control']['terminal']['a'])
    check('spectator dark-radiation transfer', all(
        r['terminal']['rho_s'] > 0 and r['terminal']['rho_dr'] > 0 for r in runs.values()))
    check('visible decay control transfers energy', runs['decay_control']['terminal']['rho_sm'] > 0)
    check('collisionless visible sector starts empty', finest['terminal']['rho_sm'] == 0)
    check('convergence table recalculated', compare(runs) == diag['comparisons'])
    check('fine time-averaged moments meet gate', diag['numerical_window_pass'])
    check('fine terminal miss remains explicit', not diag['numerical_terminal_pass'])
    check('coarse-grid miss remains explicit', not diag['comparisons'][1]['terminal_pass'])

    for key in ('flat','expanding'):
        rows = ctl['ode'][key]
        check(f'{key} unrotated ODE agreement', rows[-1]['max_relative'] < 1e-3)
        check(f'{key} convergence against ODE', rows[-1]['max_relative'] < rows[0]['max_relative']/10)
    vac = ctl['vacuum']
    check('static vacuum control', vac['static_vacuum_energy'] < 2e-13)
    check('massless control', vac['massless_vs_no_species_max_error'] < 1e-11)
    check('canonical shift charge control', vac['canonical_shift_charge_max_error'] < 1e-11)
    uv = vac['rotating_vacuum_cutoffs']
    check('rotating vacuum quadrature', max(x['quadrature_error'] for x in uv) < 2e-5)
    check('counterterm removes growing rotating vacuum contribution',
          uv[-1]['bare_energy'] > 1.2*uv[0]['bare_energy'] and
          abs(uv[-1]['subtracted_energy']-uv[-2]['subtracted_energy']) < 2e-6)

    cases = spatial['cases']
    wave = cases['wave33']
    check('all occupied spatial states stay orthonormal',
          max(r['orthonormality_error'] for r in cases.values()) < 1e-10)
    check('spatial energy budget', max(r['max_budget_error'] for r in cases.values()) < 2e-4)
    check('spatial production exceeds vacuum artifact',
          wave['final']['energy'] > 1e-3 and abs(cases['mean33']['final']['energy']) < 1e-9)
    check('homogeneous momentum blocks remain separated', cases['uniform17']['momentum_mixing'] < 1e-20)
    check('spatial cutoff and time controls', all(
        abs(cases[x]['final']['energy']/wave['final']['energy']-1) < 0.02
        for x in ('wave65','wave33_half_step')))
    check('spatial inverse dynamics', wave['inversion_error'] < 2e-10)
    check('spatial unrotated matrix ODE', spatial['unrotated_ODE'][-1]['max_absolute_error'] < 2e-5)
    check('archived lattice replay agrees in a and H', max(phase['archive_relative_errors'].values()) < 0.005)
    check('large actual phase variation despite tiny gradient budget', all(
        row['peak_operator_rms'] > 1 and row['gradient_fraction_at_peak'] < 1e-6
        for row in phase['parents'].values()))

    state_file = RES / 'h55_terminal_state.npz'
    check('saved terminal state digest', digest(state_file) == diag['published_terminal_state_sha256'])
    with np.load(state_file) as state:
        q, w, r = state['q'], state['w'], state['bloch']
        check('complete two-helicity spectral restart', r.shape == (4,2,finest['model']['nq'],3)
              and len(q) == len(w) == finest['model']['nq'] and np.all(w>0) and np.all(np.diff(q)>0))
        widths = [s['Grest'] for s in old['run']['channels']]
        recomputed = handoff(finest,q,w,r,widths)
    previous = next(row for row in diag['handoffs'] if row['case']=='eight_bins')
    check('archived coherent stress reconstructed', all(
        abs(recomputed[k]-previous[k]) < 1e-10 for k in (
            'total_pressure_change_if_dephased_over_rho','total_coherent_fermion_power_over_rho')))
    check('resolved Pauli feasibility in saved state', all(
        s['energy_over_Pauli_minimum'] >= 1 for s in recomputed['species']))
    check('dressed splice misses original energy tolerance',
          abs(recomputed['replacing_subtracted_stress_with_dressed_energy_gap']) > 1e-6)
    check('free projection loses pressure and ongoing work',
          abs(recomputed['total_pressure_change_if_dephased_over_rho']) > 1e-6 and
          abs(recomputed['total_coherent_fermion_power_over_rho']) > 1e-3)
    check('scientific status stays conditional', not diag['scientific_acceptance'] and all(
        not r['scientific_acceptance'] for r in runs.values()))

    manifest = read('inputs.json')
    check('published X4 parent pinned', manifest['parent'] == PARENT)
    check('inherited numerical inputs unchanged', all(
        hashlib.sha256((REPO / path).read_bytes()).hexdigest() == sha
        for path,sha in manifest['sha256'].items()))
    verdict = json.loads((TRACKS/'x2_strip/results/blind_verdict.json').read_text())
    check('strict scan still zero of 135', verdict['strict_M2_X_keys'] == [] and all(
        v['grid_points'] == 135 for v in verdict['nodes'].values()))

    result = dict(passed=sum(checks.values()), total=len(checks), checks=checks,
                  numerical_terminal_pass=diag['numerical_terminal_pass'],
                  numerical_window_pass=diag['numerical_window_pass'],
                  scientific_acceptance=False, strict_M2_X_accepted=0, strict_M2_X_tested=135)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'{result["passed"]}/{result["total"]} evidence checks passed')
    if not all(checks.values()):
        raise SystemExit(', '.join(name for name,ok in checks.items() if not ok))


if __name__=='__main__':
    main()
