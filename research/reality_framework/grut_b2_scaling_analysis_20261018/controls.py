"""B2-SA identity controls only: no trajectories or output probabilities.

All-size statements are analytic, not inferred from these three matrices.
Source ancestry is shared with B1/B2, not an independent implementation.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import eye, kron

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fractions():
    jp, jm, gap = F(1, 64), F(1, 96), F(5, 4)
    a = (jp * jp + jm * jm) / gap
    j = 2 * jp * jm / gap
    omega, tau = 2 - 2 * a, 4 * (jp + jm)
    radius, v = F(1, 4), tau + 2 * j
    delta = radius - 4 * a
    ratio = v / delta
    beta = radius * v / (delta * (delta - v))
    prefactor = radius / (delta - v)
    assert a == F(13, 46080) and j == F(1, 3840)
    assert v == F(67, 640)
    assert ratio < F(3, 7) and beta < F(3, 4)
    assert prefactor < F(7, 4) and delta > v
    high = 4 + F(1, 2) - a - tau - 2 * j
    assert high > 4 + 2 * j and high > 4 + radius
    variance = 2 * j * j
    weight = variance / (high - 2 * omega) ** 2
    assert weight == F(32, 37027225) and weight < F(1, 10**6)
    energy = 2 * omega + 2 * a + tau + 2 * j
    assert energy < F(411, 100)
    bmax = jp + jm
    bulk_gap = gap - 2 * bmax
    assert bulk_gap == F(115, 96)
    a0 = 1 / (1 + 2 * bmax**2 / bulk_gap**2)
    velocity = 4 * jp * jm * a0 / (gap + 2 * bmax)
    assert velocity > F(19, 10) * j
    # The B1 approximate intertwiner's allowance, NOT actual physical error.
    tau1 = 2 * (jp + jm)
    k = 2 * a + 2 * j
    r = tau1 / gap
    u = r**3 + tau1 * k / gap**2
    p = tau1 * u
    q = p + r**2 * k + u * k
    rho2 = p * p + q * q
    e0 = j * tau1 / gap**2
    length, time48 = 48, F(48) / (6 * j)
    assert rho2 * time48**2 > (F(1, 4) - e0) ** 2
    exact = dict(A=a, J=j, Omega=omega, tau2=tau, V_norm_bound=v,
                 contour_radius=radius, diagonal_distance=delta,
                 neumann_ratio=ratio, projector_difference_bound=beta,
                 projector_path_prefactor=prefactor, high_edge=high,
                 variance=variance, high_weight_bound=weight,
                 field_expectation_bound=energy, bulk_gap=bulk_gap,
                 bulk_colocated_band_weight_lower=a0,
                 bulk_velocity_coefficient=velocity,
                 bulk_velocity_over_J=velocity / j,
                 B1_residual_squared=rho2, B1_initial_error=e0)
    return ({k: {'fraction': str(x), 'decimal': float(x)}
             for k, x in exact.items()},
            {'length': length, 'time': str(time48),
             'allowance_exceeds_one_quarter': True,
             'scope': 'Old one-resource error allowance only; not observed error.'})


def matrix_controls(b2):
    result = []
    for length in (4, 6, 8):
        m = b2.b1.PathModel(length, sites=(0, length - 1), resource=2)
        ids = []
        for idx, (plus, minus, fuel) in enumerate(m.states):
            is_p = (plus == minus and
                    (len(plus) < 2 or plus[1] - plus[0] > 1))
            if is_p:
                ids.append(idx)
        count = (length * length + length + 4) // 2
        assert len(ids) == count
        mask = np.zeros(m.dim, dtype=bool)
        mask[ids] = True
        diag = m.D.diagonal()
        a, j = b2.b1.SELF, b2.b1.JPAIR
        assert min(diag[mask]) >= 4 - 4 * a - 1e-12
        assert max(diag[mask]) <= 4 + 1e-12
        assert min(diag[~mask]) >= 4 + .5 - a - 1e-12
        hop_row = float(np.max(np.asarray(abs(m.T).sum(axis=1))))
        conv_row = float(np.max(np.asarray(abs(m.C).sum(axis=1))))
        assert hop_row <= 5 / 48 + 1e-12
        assert conv_row <= 2 * j + 1e-12
        dense = m.H.toarray()
        eigenvalues, eigenvectors = eigh(dense)
        low = abs(eigenvalues - 4) < .25
        assert int(np.sum(low)) == count
        assert max(eigenvalues[low]) <= 4 + 2 * j + 1e-11
        high_bound = 4 + .5 - a - 5 / 48 - 2 * j
        assert min(eigenvalues[~low]) >= high_bound - 1e-11
        vp = eigenvectors[ids][:, low]
        gram_eigenvalues = eigh(vp @ vp.T, eigvals_only=True)
        beta_observed = float(np.sqrt(max(0, 1 - min(gram_eigenvalues))))
        assert beta_observed < .75
        source = int(np.argmax(abs(m.initial())))
        weight = float(np.sum(abs(eigenvectors[source, ~low])**2))
        assert weight < 32 / 37027225
        src = m.initial()
        mean = float(np.vdot(src, m.H @ src).real)
        variance = float(np.linalg.norm(m.H @ src - mean * src)**2)
        assert abs(mean - 2 * b2.b1.OMEGA) < 1e-12
        assert abs(variance - 2 * j * j) < 1e-18
        colored = b2.Colored(m)
        difference = colored.h - kron(m.H, eye(4), format='csr')
        factor_defect = float(np.max(abs(difference.data), initial=0))
        assert factor_defect < 1e-15
        result.append(dict(length=length, scalar_dimension=m.dim,
                           P_dimension=count, colored_dimension=len(colored.states),
                           P_diagonal_min=float(min(diag[mask])),
                           Q_diagonal_min=float(min(diag[~mask])),
                           hopping_absolute_row_max=hop_row,
                           conversion_absolute_row_max=conv_row,
                           highest_low_eigenvalue=float(max(eigenvalues[low])),
                           lowest_high_eigenvalue=float(min(eigenvalues[~low])),
                           projector_difference_norm_observed=beta_observed,
                           high_preparation_weight_observed=weight,
                           source_mean_energy=mean, source_variance=variance,
                           ordered_color_entry_defect=factor_defect))
    return result


def main():
    started = time.perf_counter()
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())['sha256']
    unchanged = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h
                 for p, h in pins.items()}
    assert all(unchanged.values())
    executed = json.loads((HERE / 'EXECUTED_SOURCES.json').read_text())['sha256']
    assert all(hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h
               for p, h in executed.items())
    b2 = load('unchanged_b2', HERE.parent /
              'grut_b2_collision_channel_20261018' / 'checks.py')
    scalar, residual = fractions()
    matrices = matrix_controls(b2)
    data = {'exact_fraction_controls': scalar,
            'old_approximation_scope_check': residual,
            'small_matrix_controls': matrices,
            'pins_unchanged': unchanged,
            'scope': 'No time evolution. Ordinary-double matrix diagnostics, not interval certification.',
            'status': 'PASS: scalar inequalities and small-matrix identities; no scaling operating verdict.'}
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
    dest.mkdir(parents=True, exist_ok=True)
    (dest / 'CONTROLS.json').write_text(json.dumps(data, indent=2) + '\n')
    metrics = {'program_seconds': time.perf_counter() - started,
               'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
               'threads': {x: os.environ.get(x) for x in
                           ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')},
               'python': sys.version, 'numpy': np.__version__}
    (dest / 'EXECUTION.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(json.dumps({'status': data['status'], 'metrics': metrics}))


if __name__ == '__main__':
    main()
