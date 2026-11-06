#!/usr/bin/env python3
"""Diagnose the old source without modifying or renormalizing its particle count."""
import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np

from trace import ARCHIVE, COLUMNS, TRACKS

SAMPLE_TIMES = (0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 2.9980937)


def minimum_energy(number, mass):
    if number == 0:
        return 0.0
    k = (3 * math.pi**2 * number)**(1 / 3)
    if k < 1e-3 * mass:
        # Nonrelativistic series, stable as number tends to zero.
        return number * (mass + 0.3 * k**2 / mass - 3 * k**4 / (56 * mass**3))
    integral = (k * math.hypot(k, mass) * (2 * k*k + mass*mass)
                - mass**4 * math.asinh(k / mass)) / 8
    return integral / math.pi**2


def summary_pair(energy, number, mass):
    return {
        'rho': float(energy),
        'n': float(number),
        'rho_over_minimum': float(energy / minimum_energy(number, mass)),
        'mean_energy_over_mass': float(energy / (number * mass)),
    }


def cohorts(t, a, de, dn, mass, width, requested):
    """Two equal-mean-energy cohorts per step, positive free streaming, no cap.

    This checks the effect of the equation of state and decay prescription.
    It is not an occupation closure. Each step's separate channel mean energy
    represents a cohort and does not resolve its native within-cell spectrum.
    """
    ratio = np.divide(de, dn, out=np.full_like(de, mass), where=dn > 0)
    assert np.min(ratio) >= mass * (1 - 5e-9)
    q = np.sqrt(np.maximum((ratio - mass) * (ratio + mass), 0)) * a[1:, None]
    n_comoving = dn * a[1:, None]**3
    history = {}
    for i in range(len(t)):
        view = n_comoving[:i]
        energy = np.hypot(q[:i] / a[i], mass)
        if i in requested:
            history[i] = (np.sum(view * energy) / a[i]**3, np.sum(view) / a[i]**3)
        if i + 1 < len(t):
            # Two midpoint scale factors agree up to quadratic errors in dt.
            midpoint = np.sqrt(a[i] * a[i+1])
            view *= np.exp(-width * mass * (t[i+1] - t[i]) / np.hypot(q[:i] / midpoint, mass))
    return history


def propagation_bound(t, a, de, dn, mass, width, index):
    """Necessary moment bound for any injection shape with the given moments."""
    fraction = a[1:index+1] / a[index]
    survival_minimum = np.exp(-width * (t[index] - t[1:index+1]))
    number_minimum = float((fraction**3) @ (dn[:index] * survival_minimum))
    # E(r*p) <= r*E(p) + (1-r)*M by convexity, even without decays.
    energy_maximum = float((fraction**3) @ (
        fraction * de[:index] + (1 - fraction) * mass * dn[:index]))
    ratio = energy_maximum / minimum_energy(number_minimum, mass)
    return dict(number_minimum=number_minimum, energy_maximum=energy_maximum,
                energy_maximum_over_Pauli_minimum=float(ratio),
                excludes_all_free_streaming_injection_shapes=bool(ratio < 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tape', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    metadata = json.loads(ARCHIVE.read_text())
    params = metadata['run']
    v2 = params['P']['v']**2
    dt = params['dt']
    last = metadata['termination']['last_diag_t']
    data = np.loadtxt(args.tape)
    d = dict(zip(COLUMNS, data.T, strict=True))
    indices = [int(np.argmin(abs(d['t'] - t))) for t in SAMPLE_TIMES]
    i_last = indices[-1]
    assert abs(d['t'][i_last] - last) < 5e-8
    t, a, hubble = d['t'], d['a'], d['H']
    total = d['sigma'] + d['chi'] + d['spectator'] + d['sm']
    total += sum(d[f'rho_{j}'] for j in range(1, 5))
    ledger_err = float(np.max(abs(total / d['rho'] - 1)))
    assert ledger_err < 2e-14

    archive = ARCHIVE.with_name('lower_H_diag.txt.gz')
    with gzip.open(archive, 'rt') as f:
        cols = f.readline().lstrip('# ').split()
        old = dict(zip(cols, np.loadtxt(f).T, strict=True))
    sampled = np.arange(0, len(t), params['ndiag'])
    assert len(sampled) == len(old['t']) and max(abs(t[sampled] - old['t'])) < 5e-8
    comparison = {col: float(max(abs(new[sampled] / old[col] - 1)))
                  for col, new in [('a', a), ('H', hubble), ('rho', d['rho'])]}
    assert max(comparison.values()) < 0.005

    entries = []
    for j, spec in enumerate(params['channels'], 1):
        m, mass, width = spec['m'], spec['M'], spec['Grest']
        rho, n, de, dn, de_pert, dn_pert = (d[f'{key}_{j}'] * v2 for key in
                                            ('rho', 'n', 'de', 'dn', 'de_pert', 'dn_pert'))
        mean = np.divide(rho, n, out=np.full_like(rho, mass), where=n > 1e-299 * v2)
        decay = width * mass / np.maximum(mass, mean)
        rho_next = rho * np.exp(-(4 * hubble + decay) * dt) + de
        n_next = n * np.exp(-(3 * hubble + decay) * dt) + dn
        errors = {
            'energy': float(max(abs(rho_next[:-1] - rho[1:])) / max(rho)),
            'number': float(max(abs(n_next[:-1] - n[1:])) / max(n)),
        }
        assert max(errors.values()) < 2e-14
        sample_energy = old[f'rN{j}'] * v2
        sample_number = old[f'nN{j}'] * v2
        replay_errors = {
            'energy_max_absolute_over_peak': float(max(abs(sample_energy - rho[sampled])) / max(rho)),
            'number_max_absolute_over_peak': float(max(abs(sample_number - n_next[sampled])) / max(n)),
        }
        assert max(replay_errors.values()) < 0.005
        assert min(np.min(de - de_pert), np.min(dn - dn_pert)) >= 0
        e = np.column_stack([de - de_pert, de_pert])[:-1]
        native_counts = np.column_stack([dn - dn_pert, dn_pert])[:-1]
        exact_redshift = cohorts(t, a, e, native_counts, mass, width, set(indices))
        kinematic_counts = native_counts.copy()
        if spec['Gp'] > 0:
            kinematic_counts[:, 1] = e[:, 1] / (m / 2)
        kinematic = cohorts(t, a, e, kinematic_counts, mass, width, set(indices))
        without_floor = cohorts(t, a, e[:, :1], native_counts[:, :1], mass, width, set(indices))
        check_times = [{
            't': float(t[i]),
            'same_time': summary_pair(rho[i], n[i], mass),
            'misaligned_number_relative_change': float(n_next[i] / n[i] - 1),
            'free_streaming_fixed_injections': summary_pair(*exact_redshift[i], mass),
            'free_streaming_and_parent_at_rest': summary_pair(*kinematic[i], mass),
            'free_streaming_without_perturbative_floor': summary_pair(*without_floor[i], mass),
        } for i in indices]

        pert = None
        if spec['Gp'] > 0:
            # Generous support used by the inherited off-shell line at zero
            # thermal mass: E_N <= (m_parent^2 + M_N^2) / (2 m_parent).
            e_max = (m*m + mass*mass) / (2*m)
            p_max = math.sqrt((e_max - mass) * (e_max + mass))
            density_max = p_max**3 / (3 * math.pi**2)  # g=2, f<=1
            ab = a[1:i_last + 1]**3
            request = float(ab @ de_pert[:i_last])
            legacy_count = float(ab @ dn_pert[:i_last])
            # Upper-bound initial/fresh available modes by the entire final
            # comoving ball and every loss by the vacuum rest decay rate.
            max_injected_count = density_max * (a[i_last]**3 + width * dt * np.sum(ab))
            capacity = e_max * max_injected_count
            pert = {
                'E_injected_comoving': request,
                'fraction_of_all_injected_energy': request / float(ab @ de[:i_last]),
                'average_legacy_injection_energy': request / legacy_count,
                'stationary_on_shell_injection_energy': m / 2,
                'generous_stationary_support_E_max': e_max,
                'generous_stationary_support_p_max': p_max,
                'max_number_injected_comoving_with_vacuum_decay': float(max_injected_count),
                'max_energy_injected_comoving_with_vacuum_decay': float(capacity),
                'requested_over_capacity_bound': float(request / capacity),
            }
        entries.append({
            'species': f'N{j}', 'mass': mass, 'parent_mass': m, 'decay_width': width,
            'recurrence_error_over_peak': errors,
            'archive_difference_over_peak': replay_errors,
            'perturbative_budget': pert, 'samples': check_times,
            'shape_independent_propagation_bound': propagation_bound(
                t, a, de, dn, mass, width, i_last),
            'shape_independent_bound_without_floor': propagation_bound(
                t, a, de - de_pert, dn - dn_pert, mass, width, i_last),
        })
    result = {
        'scope': 'same-time legacy source replay and counterfactual free propagation of its requested injections',
        'source_ledger_max_relative_error': ledger_err,
        'archive_max_relative_differences': comparison,
        'source_budget_max': float(np.max(abs(old['budget']))),
        'species': entries,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    for entry in entries:
        row = entry['samples'][-1]
        print(entry['species'], *(row[key]['rho_over_minimum'] for key in
              ('same_time', 'free_streaming_fixed_injections', 'free_streaming_and_parent_at_rest')),
              (entry['perturbative_budget'] or {}).get('requested_over_capacity_bound'), flush=True)


if __name__ == '__main__':
    main()
