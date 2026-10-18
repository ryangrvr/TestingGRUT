"""B2-SS source spectral/algebra controls. NO time evolution.

Finite source matrices are exact compressions of the half-line law, not
new finite-path candidates or all-size source certifications.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import resource
import sys
import time

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.sparse import diags

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def b1_module():
    path = HERE.parent/'grut_b1_neutral_escape_20261018'/'model.py'
    spec = importlib.util.spec_from_file_location('pinned_source_b1', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scalars():
    jp, jm, delta = F(1, 64), F(1, 96), F(5, 4)
    a, j = F(13, 46080), F(1, 3840)
    b0, bp = jp+jm, jp-jm
    q0 = (2*b0*b0/delta)/(1+2*b0*b0/delta**2)
    d0 = q0-2*a
    dp = 2*a-2*bp*bp/(delta-2*bp)
    assert d0 > 3*j/2 and dp > 3*j/2
    outside = j*j/(j*j+dp*dp)
    assert outside < F(201, 1000)
    assert 1-F(1)/(F(3, 2)**2) == F(5, 9)
    # Positive root residual d/2-J^2/d increases in positive d.
    assert F(3, 2)/2-1/F(3, 2) == F(1, 12)
    speed = F(16928, 8825)
    assert speed**2*F(215, 576) > F(11, 10)**2
    assert F(3, 2)+F(1, 12) == F(19, 12)
    tau = 2*(jp+jm)
    k = 2*a+2*j
    r = tau/delta
    u = r**3+tau*k/delta**2
    p = tau*u
    q = p+r*r*k+u*k
    rho2 = p*p+q*q
    assert rho2 < F(1, 105000)**2
    assert F(1, 105000)/j == F(32, 875)
    e_init = j*tau/delta**2
    assert e_init == F(1, 115200)
    bound2 = F(1, 2)+e_init+F(32, 875)*2
    bound16 = F(1, 16)+e_init+F(32, 875)*16
    assert max(bound2, bound16) < F(13, 20) < F(2, 3)
    exponent16 = -F(3, 20)*(2016//2-2)+F(11, 28)*16
    assert exponent16 < -100
    assert F(101, 100)/2**100 < F(1, 10**6)
    retained = 1-F(1, 10**6)-(F(13, 20)+F(1, 10**6))**2-F(1, 10**12)
    physical_lower = (retained-F(1, 1000))**2
    assert retained > F(1, 1000) and physical_lower > F(33, 100)
    assert (F(17, 50)-F(1, 2000)-F(1, 1000))**2 == F(458329, 4000000)
    assert F(458329, 4000000) > F(1, 10)
    values = dict(lower_band_margin_bound=d0,
                  lower_margin_in_J=d0/j,
                  upper_band_margin_bound=dp,
                  upper_margin_in_J=dp/j,
                  above_ground_band_weight_upper=outside,
                  in_ground_band_weight_lower=1-outside,
                  central_interval_weight_lower=F(5, 9),
                  compatible_velocity_squared_lower=speed**2*F(215, 576),
                  B1_residual_squared=rho2,
                  fuel_bound_at_theta2=bound2,
                  fuel_bound_at_theta16=bound16,
                  retained_matter_source_lower=retained,
                  joint_q2_lower_on_theta2_to16=physical_lower,
                  conditional_joint_output_lower=F(458329, 4000000))
    return {name: {'fraction': str(x), 'decimal': float(x)} for name, x in values.items()}


def finite_controls(b1, out):
    cutoff = 32
    rr = np.arange(-cutoff, cutoff+1)
    dd = 2+abs(rr)+.25*(abs(rr)==1)
    e0v, vec = eigh_tridiagonal(dd, -np.full(len(rr)-1, b1.JP+b1.JM),
                               select='i', select_range=(0, 0))
    epiv = eigh_tridiagonal(dd, -np.full(len(rr)-1, b1.JP-b1.JM),
                           select='i', select_range=(0, 0), eigvals_only=True)
    e0, epi = float(e0v[0]), float(epiv[0])
    ur = vec[:, 0]
    if ur[cutoff] < 0:
        ur = -ur
    d = b1.OMEGA-e0
    assert d > 1.5*b1.JPAIR
    assert epi-b1.OMEGA > 1.5*b1.JPAIR
    bulk = dict(relative_cutoff=cutoff, ground_energy_at_k0=e0,
                ground_energy_at_kpi=epi, relative_coordinates=rr.tolist(),
                computed_ground_vector=ur.tolist(),
                scope='Relative cutoff is a diagnostic only, not the physical half-line source.')
    (out/'BULK_CONTROL.json').write_text(json.dumps(bulk, indent=2)+'\n')
    rows = []
    for length in (8, 12, 16):
        m = b1.PathModel(length, sites=(0,), resource=1)
        # Compress the half-line, retaining only its root counterterm.
        correction = np.zeros(m.dim)
        correction[m.index[((length-1,), (length-1,), 0)]] = b1.SELF
        h = m.H+diags(correction, format='csr')
        source = m.index[((), (), 1)]
        assert sum(f == 1 for p, n, f in m.states) == 1
        f = np.zeros(m.dim); f[source] = 1
        mean = float(f@(h@f))
        variance = float(np.linalg.norm(h@f-mean*f)**2)
        assert abs(mean-b1.OMEGA) < 1e-12
        assert abs(variance-b1.JPAIR**2) < 1e-18
        ev, vv = eigh(h.toarray())
        weight = abs(vv[source])**2
        low = ev < 2.5
        assert int(sum(low)) == length+1
        assert ev[0] >= e0-1e-11
        high_weight = float(sum(weight[~low]))
        assert high_weight < 36/27613**2
        central = abs(ev-b1.OMEGA) <= 1.5*b1.JPAIR
        central_weight = float(sum(weight[central]))
        assert central_weight >= 5/9-1e-10
        band_weight = float(sum(weight[(ev >= e0-1e-11) & (ev <= epi)]))
        assert band_weight > .799-1e-10
        witness = np.zeros(m.dim)
        for idx, (plus, minus, fuel) in enumerate(m.states):
            witness[idx] = (b1.JPAIR*ur[cutoff]/d if fuel else
                            ur[cutoff+minus[0]-plus[0]])
        residual = h@witness-e0*witness
        root = m.index[((0,), (0,), 0)]
        root_predicted = (d/2-b1.JPAIR**2/d)*ur[cutoff]
        assert min(residual) > -1e-10
        assert abs(residual[source]) < 1e-11
        assert abs(residual[root]-root_predicted) < 1e-11
        assert residual[root] > (b1.JPAIR/12)*ur[cutoff]-1e-11
        raw = dict(compression_length=length, source_dimension=m.dim,
                   diagonal_far_end_correction=float(b1.SELF),
                   all_eigenvalues=ev.tolist(), all_fuel_residues=weight.tolist(),
                   scope='Finite compression only. No infinite-volume trapped/outgoing classification or time evolution.')
        (out/f'SPECTRUM_{length}.json').write_text(json.dumps(raw, indent=2)+'\n')
        rows.append(dict(length=length, dimension=m.dim,
                         fuel_subspace_scalar_rank=1,
                         source_mean=mean, source_variance=variance,
                         minimum_eigenvalue=float(ev[0]),
                         ground_band_interval_weight_observed=band_weight,
                         central_energy_weight_observed=central_weight,
                         high_weight_observed=high_weight,
                         minimum_computed_supersolution_residual=float(min(residual)),
                         root_computed_residual=float(residual[root]),
                         root_residual_identity_error=float(abs(residual[root]-root_predicted)),
                         witness_fuel_residual=float(residual[source]),
                         note='Ordinary double controls, not outward-rounded proof.'))
    return rows


def main():
    started = time.perf_counter()
    for name in ('SOURCE_PINS.json', 'EXECUTED_SOURCES.json'):
        for path, expected in json.loads((HERE/name).read_text())['sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected
    out = Path(sys.argv[1]) if len(sys.argv)>1 else HERE
    out.mkdir(parents=True, exist_ok=True)
    exact = scalars()
    finite = finite_controls(b1_module(), out)
    result = dict(exact_scalar_checks=exact, finite_source_controls=finite,
                  status='PASS: explicit scalar/spectral/algebra controls only',
                  full_source_scaling_verdict='MATHEMATICALLY UNRESOLVED',
                  no_source_time_evolution=True, no_operating_probabilities=True,
                  source_and_proof_pins_unchanged=True)
    (out/'CONTROLS.json').write_text(json.dumps(result, indent=2)+'\n')
    metrics = dict(program_seconds=time.perf_counter()-started,
                   peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                   python=sys.version, numpy=np.__version__,
                   threads={x:os.environ.get(x) for x in
                            ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
    (out/'EXECUTION.json').write_text(json.dumps(metrics, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'metrics': metrics}))


if __name__ == '__main__':
    main()
