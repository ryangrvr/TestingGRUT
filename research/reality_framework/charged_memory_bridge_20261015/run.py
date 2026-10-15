"""Fresh CP-12 memory, physical-cutoff and reservoir-reduction controls."""
from __future__ import annotations

from pathlib import Path
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, eye
from scipy.sparse.linalg import expm_multiply

from memory import MemoryBridge, reservoir_krylov

P = Path(__file__).resolve().parent
RF = P.parent
ALLOWANCE = 5e-10


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, RF/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def interval(probability, state_error):
    radius = state_error + ALLOWANCE
    return [max(0., np.sqrt(max(0., probability))-radius)**2,
            min(1., (np.sqrt(max(0., probability))+radius)**2)]


def toy_controls():
    v = .7
    bridge = MemoryBridge(np.zeros((1, 1)), np.array([[v]]),
                          np.zeros((1, 1)), np.ones(1))
    result = []
    for steps in [250, 500, 1000]:
        row = bridge.solve(5., steps, [5.])[0]
        error = float(abs(row['x'][0] + 1j*np.sin(v*5.)))
        result.append(dict(steps=steps, amplitude_error=error))
    assert result[-1]['amplitude_error'] < 1e-5
    assert 3.8 < result[0]['amplitude_error']/result[1]['amplitude_error'] < 4.2
    rng = np.random.default_rng(20261015)
    A = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)); A=(A+A.conj().T)/6
    D = rng.normal(size=(5, 5))+1j*rng.normal(size=(5, 5)); D=(D+D.conj().T)/10
    C = .2*(rng.normal(size=(3, 2))+1j*rng.normal(size=(3, 2))) @ (
        rng.normal(size=(2, 5))+1j*rng.normal(size=(2, 5)))
    z = rng.normal(size=5)+1j*rng.normal(size=5); z /= np.linalg.norm(z)
    bridge = MemoryBridge(A, C, D, z)
    H = np.block([[A, C], [C.conj().T, D]])
    w = expm_multiply(-2j*H, np.r_[np.zeros(3), z])[:3]
    row = bridge.solve(2., 1000, [2.])[0]
    error = float(np.linalg.norm(w-row['x']))
    assert bridge.rank == 2 and error < 1e-5
    return dict(degenerate_two_level=result, complex_rank_two_amplitude_error=error)


def main():
    started = time.monotonic()
    formation = load('charged_memory_formation', 'cp12_unmasked_formation/run.py')
    record = load('charged_memory_record',
        'constructive_sprint_20261015/protected_geometry/neutral_record.py')
    sys.path.insert(0, str(RF/'cp12_distance_history'))
    from evolution import Evolution
    controls = toy_controls()
    times = [2., 5., 10.]
    reference_runs = []
    shift = 2.
    for budget in [128, 256]:
        _, H, outside, states, _, psi = formation.model(4, 4, budget)
        target = np.array([bool(s[1]&15) and record.circulation(s) for s in states])
        ip, iq = np.flatnonzero(target), np.flatnonzero(~target)
        A = H[ip][:, ip].toarray()-shift*np.eye(len(ip))
        C = H[ip][:, iq].toarray()
        D = H[iq][:, iq].toarray()-shift*np.eye(len(iq))
        bridge = MemoryBridge(A, C, D, psi[iq])
        row = dict(electric_square_budget=budget, invariant_dimension=H.shape[0],
            charged_dimension=len(ip), complement_dimension=len(iq),
            numerical_coupling_rank=bridge.rank,
            coupling_factorization_spectral_residual=bridge.coupling_residual,
            direct=[])
        assert np.linalg.norm(psi[ip]) == 0
        for t in times:
            w, _, error = Evolution(H, outside, shift, steps=256)(psi, t)
            g, _ = bridge.integrated_blocks(t)
            probability = float(np.linalg.norm(w[ip])**2)
            first_transfer = float(np.linalg.norm(g)**2)
            row['direct'].append(dict(time=t, charged_probability=probability,
                full_law_state_error_upper_float=float(error),
                probability_interval_with_reporting_allowance=interval(probability, error),
                first_transfer_probability=first_transfer,
                first_transfer_over_full_probability=first_transfer/probability))
        reference_runs.append(row)
        print(json.dumps(dict(stage='direct', **row)), flush=True)
    # The high-budget model is the last one, retained in local variables.
    full_eigenvalues, full_eigenvectors = eigh(H.toarray()-shift*np.eye(H.shape[0]))
    coeff = full_eigenvectors.conj().T @ psi
    references = {t: (full_eigenvectors @ (np.exp(-1j*full_eigenvalues*t)*coeff))[ip]
                  for t in times}
    grids = []
    for steps in [250, 500, 1000, 2000]:
        readouts = bridge.solve(10., steps, times)
        rows = []
        for r in readouts:
            r['amplitude_error_against_finite_direct'] = float(np.linalg.norm(
                r.pop('x')-references[r['time']]))
            r['probability_error_against_finite_direct'] = abs(r['probability']-
                float(np.linalg.norm(references[r['time']])**2))
            rows.append(r)
        grids.append(dict(steps=steps, dt=10./steps, readouts=rows))
        print(json.dumps(dict(stage='volterra', **grids[-1])), flush=True)
    assert max(r['amplitude_error_against_finite_direct'] for r in grids[-1]['readouts']) < 1e-4
    assert all(3.5 < grids[-2]['readouts'][i]['amplitude_error_against_finite_direct'] /
        grids[-1]['readouts'][i]['amplitude_error_against_finite_direct'] < 4.5 for i in range(3))
    full_reference = expm_multiply(-5j*(H-shift*eye(H.shape[0])), psi,
        traceA=-5j*float(H.diagonal().sum()-shift*H.shape[0]))
    reduced = []
    for depth, V in reservoir_krylov(D, C, psi[iq], 70):
        if depth not in [24, 50, 70]:
            continue
        width = len(ip)+V.shape[1]
        E = np.zeros((H.shape[0], width), dtype=V.dtype)
        E[ip, np.arange(len(ip))] = 1.
        E[np.ix_(iq, np.arange(len(ip), width))] = V
        small = E.conj().T @ (H@E)
        small = (small+small.conj().T)/2
        v0 = E.conj().T@psi
        vals, vecs = eigh(small)
        w = vecs @ (np.exp(-5j*(vals-shift))*(vecs.conj().T@v0))
        probability = float(np.linalg.norm(w[:len(ip)])**2)
        r = dict(krylov_depth=depth, charged_dimension=len(ip),
            complement_reservoir_dimension=V.shape[1], total_dimension=width,
            time=5., charged_probability=probability,
            finite_model_state_error=float(np.linalg.norm(E@w-full_reference)),
            isometry_frobenius_error=float(np.linalg.norm(E.conj().T@E-np.eye(width))))
        if depth in [50, 70]:
            # Retained-space and electric-boundary residuals have orthogonal
            # ranges in the original unbounded invariant Hilbert space.
            residual = np.vstack([H@E-E@small, outside@E])
            evolution = Evolution(csr_matrix(small), csr_matrix(residual), shift, steps=256)
            _, _, error = evolution(v0, 5.)
            error += float(np.linalg.norm(psi-E@v0))
            r['full_law_state_error_upper_float'] = float(error)
            r['probability_interval_with_reporting_allowance'] = interval(probability, error)
            assert r['finite_model_state_error'] <= error+ALLOWANCE
        reduced.append(r)
        print(json.dumps(dict(stage='reservoir', **r)), flush=True)
    assert reduced[-1]['full_law_state_error_upper_float'] < 1e-6
    assert reduced[-1]['isometry_frobenius_error'] < 1e-11
    reference_probability = reference_runs[-1]['direct'][1]['charged_probability']
    lower, upper = reduced[-1]['probability_interval_with_reporting_allowance']
    assert lower <= reference_probability <= upper
    kernel0 = np.linalg.norm(bridge.memory_kernel(0.), ord=2)
    kernel_samples = [dict(time=t, spectral_norm=float(np.linalg.norm(
        bridge.memory_kernel(t), ord=2)), norm_relative_to_zero=float(np.linalg.norm(
        bridge.memory_kernel(t), ord=2)/kernel0)) for t in [0., .5, 1., 2., 5., 10.]]
    source_files = [RF/'cp12_ground_connectivity/frozen_model.py',
        RF/'cp12_ground_connectivity/exact_destination.py',
        RF/'cp12_ground_connectivity/invariant_evolution.py',
        RF/'cp12_unmasked_formation/run.py', RF/'cp12_distance_history/evolution.py',
        RF/'constructive_sprint_20261015/protected_geometry/neutral_record.py']
    audit = {str(f.relative_to(RF)): hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files}
    assert audit['cp12_ground_connectivity/frozen_model.py'] == 'e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
    result = dict(status='PASS', parent='3538c9a44afa34a7d8e5ce94d620b00328c3b2df',
        frozen_law_changed=False, original_neutral_ray=True,
        external_record_or_state_filter=False, toy_controls=controls,
        reference_runs=reference_runs, volterra_grid_study=grids,
        reservoir_reductions=reduced, interface_kernel_samples=kernel_samples,
        numerical_reporting_allowance=ALLOWANCE,
        cache_used=bool(os.environ.get('CP12_FORMATION_CACHE')),
        error_scope='Analytic Duhamel upper sums evaluated in ordinary floating point, not outward-rounded intervals. The Volterra discretization is checked by convergence and direct comparison, not certified by a rigorous quadrature bound. Reservoir bounds include full outgoing electric residuals.',
        physical_scope='N4,Q4 validation of a sourced memory reduction; not growing Q=2N evidence, autonomous record formation, dimension selection or gravity.',
        scaling_scope='Interface-history dimension seven in this cutoff control. Full spectral preprocessing still uses the complement; no proved large-N compression.',
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        seconds=time.monotonic()-started)
    (P/'SOURCE_AUDIT.json').write_text(json.dumps(audit, indent=2)+'\n')
    (P/'RESULTS.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], seconds=result['seconds'])), flush=True)


if __name__ == '__main__':
    main()
