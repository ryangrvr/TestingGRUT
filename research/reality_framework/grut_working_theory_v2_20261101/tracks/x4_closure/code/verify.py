#!/usr/bin/env python3
"""Check archived X4 evidence and optionally compare a fresh replay."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from trace import TRACK, TRACKS

BASE = '2b0316237a8f8a85b3db22453eb2dfb57231d1bc'
RESULTS = TRACK / 'results'
REPO = TRACK.parents[4]


def read(path):
    return json.loads(path.read_text())


def verify(directory):
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    audit = read(directory / 'source_audit.json')
    species = audit['species']
    n1 = species[0]
    final = n1['samples'][-1]
    check('aligned source ledger', audit['source_ledger_max_relative_error'] < 2e-14)
    check('source scalar evolution replay', max(audit['archive_max_relative_differences'].values()) < 5e-3)
    check('legacy moment recurrences', all(max(s['recurrence_error_over_peak'].values()) < 2e-14 for s in species))
    check('legacy archive column timing', all(max(s['archive_difference_over_peak'].values()) < 5e-3 for s in species))
    check('aligned N1 still fails Pauli', final['same_time']['rho_over_minimum'] < 1)
    check('free streaming cohort control still fails', final['free_streaming_fixed_injections']['rho_over_minimum'] < 1)
    check('removing the floor in cohort control still fails',
          final['free_streaming_without_perturbative_floor']['rho_over_minimum'] < 1)
    check('independent-of-shape propagation obstruction',
          n1['shape_independent_propagation_bound']['energy_maximum_over_Pauli_minimum'] < 1)
    check('stationary spectral support obstruction', all(
        species[i]['perturbative_budget']['requested_over_capacity_bound'] > 1 for i in (0, 3)))
    check('legacy source misses scalar budget gate', audit['source_budget_max'] > 0.01)

    kinetic = read(directory / 'kinetic.json')
    cases = {row['case']: row for row in kinetic['cases']}
    check('eight conditional kinetic cases', len(cases) == 8)
    check('spectator ODE reference controls',
          len(kinetic['spectator_controls']) == 6 and
          max(row['relative_error'] for row in kinetic['spectator_controls']) < 5e-13)
    check('full-source energy with spectator', all(
        abs(r['full_source_energy_relative_gap']) < 1e-6 for r in cases.values() if r['kind'] != 'omitted'))
    check('missing spectator reproduces energy gap',
          abs(cases['matched_omitted_s']['full_source_energy_relative_gap']) > 1e-6)
    check('exact fixed-volume transfer', all(r['maximum_fixed_volume_exchange_error'] < 1e-10 for r in cases.values()))
    check('FRW pressure-work budget', all(r['maximum_FRW_budget_error'] < 0.01 for r in cases.values()))
    check('bounded kinetic occupation', all(
        -1e-14 <= r['occupation_range'][0] <= r['occupation_range'][1] <= 1 + 1e-14 for r in cases.values()))
    check('spectator comoving decay law', all(
        r['maximum_spectator_comoving_decay_error'] < 5e-11 for r in cases.values()))
    check('late spectator depletion', all(
        r['final_spectator_fraction'] < 1e-9 for r in cases.values() if r['kind'] != 'harmonic_stop'))
    check('dark radiation stays separate from visible energy', all(
        r['K_total'] > r['K_visible'] and r['rho_DR_over_SM'] > 0
        for r in (cases['matched'], cases['harmonic'])))
    check('late K momentum and timestep checks', all(
        abs(cases[base + suffix][key] / cases[base][key] - 1) < 0.01
        for base in ('matched', 'harmonic') for suffix in ('_fine_q', '_half_step')
        for key in ('K_total', 'K_visible')))
    check('prepared N1 count explicitly differs', all(
        r['initialized_N_count'][0] < r['native_source_N_count'][0]
        for r in (cases['matched'], cases['harmonic'])))
    check('conditional cases retain scientific decision', all(not r['accepted_cosmological_point'] for r in cases.values()))

    modes = read(directory / 'mean_modes.json')
    probe = modes['cases']
    check('five mean-waveform solver resolutions', len(probe) == 5)
    check('mean-waveform static control', modes['controls']['maximum_static_occupation'] < 1e-23)
    check('independent unrotated integrator', modes['controls']['maximum_ODE_occupation_errors']['4'] < 2e-8)
    check('unitary propagator', max(v['norm_error'] for v in probe.values()) < 1e-8)
    check('bounded mean-waveform modes', max(v['maximum_occupation'] for v in probe.values()) <= 1 + 1e-14)
    check('realizable mean-waveform moments', min(
        s[mode]['rho_over_minimum'] for v in probe.values() for row in v['samples'] for s in row['species']
        for mode in ('instantaneous', 'dressed')) > 1 - 3 / 720**2)
    check('persistent N1 mean-waveform difference', all(
        v['samples'][-1]['species'][0]['dressed']['rho_over_legacy'] < 0.6 for v in probe.values()))
    check('mean-waveform probe remains conditional', not modes['accepted_cosmological_point'])

    restart = read(directory / 'restart.json')
    check('restart at archived source time', restart['interface_time'] == kinetic['final_archived_time'])
    check('restart preserves total energy', abs(restart['restart_energy_relative_gap']) < 1e-9)
    check('restart agrees with uninterrupted history', max(
        abs(restart['K_total_relative_difference']), abs(restart['K_visible_relative_difference'])) < 0.01)
    check('restart remains conditional', not restart['accepted_cosmological_point'])
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks = verify(RESULTS)
    manifest = read(RESULTS / 'inputs.json')
    checks['input parent pinned'] = manifest['parent'] == BASE
    checks['inherited input digests'] = all(
        hashlib.sha256((REPO / filename).read_bytes()).hexdigest() == digest
        for filename, digest in manifest['sha256'].items())
    blind = read(TRACKS / 'x2_strip/results/blind_verdict.json')
    checks['strict scan still zero of 135'] = (
        blind['strict_M2_X_keys'] == [] and
        all(v['grid_points'] == 135 for v in blind['nodes'].values()))
    controls = read(RESULTS / 'mean_modes.json')
    checks['mode source digest'] = controls['source_sha256'] == hashlib.sha256(
        (TRACK / 'code/mean_modes.c').read_bytes()).hexdigest()
    if args.replay:
        checks.update({f'replay {name}': passed for name, passed in verify(args.replay).items()})
        archived = read(RESULTS / 'kinetic.json')
        replayed = read(args.replay / 'kinetic.json')
        checks['fresh replay late K'] = all(
            a['case'] == b['case'] and math.isclose(a['K_total'], b['K_total'], rel_tol=0.01)
            for a, b in zip(archived['cases'], replayed['cases'], strict=True))
    result = dict(passed=sum(checks.values()), total=len(checks), checks=checks,
                  scientific_acceptance=False, strict_M2_X_accepted=0, strict_M2_X_tested=135)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(f'{result["passed"]}/{result["total"]} evidence checks passed')
    if not all(checks.values()):
        raise SystemExit(', '.join(name for name, passed in checks.items() if not passed))


if __name__ == '__main__':
    main()
