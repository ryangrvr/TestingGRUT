#!/usr/bin/env python3
"""Compare existing final field slices; execute no physical evolution.

Run after evolve.py. Writes convergence_slices.json beside this script.
The output is a two-dimensional slice comparison, not a full 3D solution norm.
"""
from pathlib import Path
import json

import numpy as np
from scipy.signal import resample


ROOT = Path(__file__).resolve().parent


def resize(field, count):
    """Periodic Fourier resampling along the two stored spatial axes."""
    return resample(resample(field, count, axis=0), count, axis=1)


def measure(difference, reference, mask, spacing):
    difference_squared = float(np.sum(abs(difference[mask]) ** 2) * spacing ** 2)
    reference_squared = float(np.sum(abs(reference[mask]) ** 2) * spacing ** 2)
    return {
        'max_abs': float(np.max(abs(difference[mask]))),
        'L2_abs_plane': float(np.sqrt(difference_squared)),
        'L2_relative_plane': (
            float(np.sqrt(difference_squared / reference_squared))
            if reference_squared else None
        ),
    }


def region_metrics(A, S, phase, A_reference, S_reference, radius_squared, spacing):
    result = {}
    regions = [
        ('core_r8', radius_squared <= 64),
        ('core_r12', radius_squared <= 144),
        ('whole_plane', np.ones(radius_squared.shape, dtype=bool)),
    ]
    for name, mask in regions:
        result[name] = {
            'A_raw': measure(A - A_reference, A_reference, mask, spacing),
            'A_phase_aligned': measure(
                A * np.exp(-1j * phase) - A_reference,
                A_reference, mask, spacing,
            ),
            'S_raw_organization_normalized': measure(
                S - S_reference, S_reference - 1, mask, spacing,
            ),
        }
    return result


def main():
    base = np.load(ROOT / 'evolution_stress_slice.npz')
    fine = np.load(ROOT / 'evolution_space_slice.npz')
    domain = np.load(ROOT / 'evolution_domain_slice.npz')
    coordinates = base['coords']
    spacing = float(coordinates[1] - coordinates[0])
    radius_squared = coordinates[:, None] ** 2 + coordinates[None, :] ** 2
    A0 = base['A_real'] + 1j * base['A_imag']
    A1 = fine['A_real'] + 1j * fine['A_imag']
    S0, S1 = base['S'], fine['S']
    A1_resampled = resize(A1, len(coordinates))
    S1_resampled = resize(S1, len(coordinates))
    phase = float(np.angle(np.vdot(A0, A1_resampled)))

    space = {
        'fine_grid_N': 128,
        'reference_grid_N': 96,
        'same_box_side': 64.0,
        'same_dt': 0.1,
        'tau': 40.0,
        'resampling': (
            'scipy.signal.resample along both periodic slice axes, fine N128 '
            'to reference N96; field slices have common physical z=0 and '
            'x/y origin=-L/2'
        ),
        'global_A_phase_alignment_radians': phase,
        'metrics': region_metrics(
            A1_resampled, S1_resampled, phase, A0, S0,
            radius_squared, spacing,
        ),
    }
    fine_coordinates = fine['coords']
    fine_spacing = float(fine_coordinates[1] - fine_coordinates[0])
    fine_mask = np.ones(A1.shape, dtype=bool)
    space['fine_to_reference_to_fine_resampling_remainder'] = {
        'A': measure(
            resize(A1_resampled, 128) - A1, A1, fine_mask, fine_spacing,
        ),
        'S_organization_normalized': measure(
            resize(S1_resampled, 128) - S1, S1 - 1, fine_mask, fine_spacing,
        ),
    }

    larger_coordinates = domain['coords']
    indices = np.array([
        np.argmin(abs(larger_coordinates - value)) for value in coordinates
    ])
    A_domain = (domain['A_real'] + 1j * domain['A_imag'])[np.ix_(indices, indices)]
    S_domain = domain['S'][np.ix_(indices, indices)]
    domain_phase = float(np.angle(np.vdot(A0, A_domain)))
    box = {
        'reference_box_side': 64.0,
        'larger_box_side': 256 / 3,
        'same_h': spacing,
        'same_dt': 0.1,
        'tau': 40.0,
        'comparison': 'exact coordinate subset of saved z=0 slice; no resampling',
        'max_coordinate_mismatch': float(np.max(abs(
            larger_coordinates[indices] - coordinates
        ))),
        'global_A_phase_alignment_radians': domain_phase,
        'metrics': region_metrics(
            A_domain, S_domain, domain_phase, A0, S0,
            radius_squared, spacing,
        ),
    }
    output = {
        'scope': (
            'direct saved final two-dimensional field-slice comparisons only; '
            'neither a full three-dimensional solution norm nor a continuum '
            'error enclosure'
        ),
        'space': space,
        'domain': box,
        'caveats': [
            'The time discretization is unchanged in the space comparison.',
            'A single central plane does not constrain every off-plane mode.',
            'Core masks use reference-grid quadrature.',
            'Global phase alignment removes only a common U(1) phase; raw errors are retained.',
            'Spectral resampling is analysis of existing outputs, not a new physical run.',
        ],
    }
    (ROOT / 'convergence_slices.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({
        'space_core_r12': space['metrics']['core_r12'],
        'domain_core_r12': box['metrics']['core_r12'],
        'scope': output['scope'],
    }, indent=2))


if __name__ == '__main__':
    main()
