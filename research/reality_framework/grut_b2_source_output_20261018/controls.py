"""B2-SO exact scalar / small embedding identities; NO propagation.

Reuses the archived B1 scalar catalog. Not an independent implementation.
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
from scipy.sparse import coo_matrix, eye, kron

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def load_b1():
    path = HERE.parent / 'grut_b1_neutral_escape_20261018' / 'model.py'
    spec = importlib.util.spec_from_file_location('unchanged_b1_for_so', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scalar_checks():
    jp, jm, j = F(1, 64), F(1, 96), F(1, 3840)
    a, gap, mu = F(13, 46080), F(5, 4), F(3, 10)
    b, alpha, x = F(35, 576), 2*a+7*j/6, F(1, 19)
    mapping = b/gap+(b+alpha)*x/gap+b*x*x/gap
    contraction = (b+alpha+2*b*x)/gap
    ell = alpha+b*x
    epsilon = b*x*((b/gap)**2+(1+b/gap)*ell/gap)
    complement = (gap-b-2*b*x-alpha*x*x)/(1+x*x)
    assert mapping < x and contraction < F(1, 16)
    assert ell < F(1, 200) and complement > F(1, 2)
    assert epsilon < j/12
    assert (1+x*x)/(1-x*x)**2 < F(101, 100)**2
    exp15_bound = 1+F(3, 20)/(1-F(3, 40))
    assert exp15_bound == F(43, 37) and exp15_bound < F(7, 6)
    assert F(7, 6)-F(6, 7) == F(13, 42)
    c = (F(13, 42)+F(1, 12))*j
    assert c == F(11, 28)*j
    high = 2+gap-2*(jp+jm)
    omega = 2-2*a
    ratio = j/(high-omega)
    assert ratio == F(6, 27613) and ratio < F(1, 4000)
    prefactor = 4*(jp+jm)*F(101, 100)/c*F(10, 7)
    assert prefactor < 1500
    assert -F(3, 40)+F(11, 168) == -F(1, 105)
    length = 2016
    z = F(length, 105)
    lower_exp = sum(z**n/F(factorial(n)) for n in range(17))
    delta_bound = 2*ratio+1500/lower_exp
    tail_bound = F(101, 100)*F(10, 7)/lower_exp
    assert delta_bound < F(1, 1000)
    assert tail_bound < F(1, 1000)
    assert (F(1, 3)-F(1, 1000))**2 > F(1, 10)
    assert (F(5, 9)-F(1, 1000)-F(2, 10**6))**2 > F(1, 4)
    short_exponent = -F(3, 20)*(length//2-2)+c*(4/j)
    assert short_exponent < -100
    assert F(101, 100)/2**100 < F(1, 10**6)
    root_lower = F(3778, 5000)
    assert root_lower**2 < F(22839, 40000)
    assert (root_lower-F(1, 4000)-F(1, 10**6))**2 > F(57, 100)
    short_joint = (F(57, 100)-F(1, 1000))**2
    assert short_joint == F(323761, 10**6) and short_joint > F(1, 4)
    values = dict(mu=mu, graph_radius=x, hopping_weight_bound=b,
                  AP_weight_bound=alpha, fixed_point_ball_image=mapping,
                  contraction_bound=contraction, L1_norm_bound=ell,
                  exact_correction_norm_bound=epsilon,
                  exact_correction_over_J=epsilon/j,
                  complement_gap_lower=complement,
                  weighted_growth_rate=c, initial_high_norm_bound=ratio,
                  delta_prefactor_before_rounding=prefactor,
                  exp_at_L0_lower=lower_exp,
                  delta_at_L0_upper=delta_bound,
                  tail_at_L0_upper=tail_bound,
                  short_joint_retention_lower=short_joint)
    return {key: {'fraction': str(value), 'decimal': float(value)}
            for key, value in values.items()}


def max_entry(matrix):
    return float(np.max(abs(matrix.data), initial=0))


def embeddings(b1):
    output = []
    for length in (4, 6, 8):
        cutoff = length//2-2
        one = b1.PathModel(cutoff+3, sites=(0,), resource=1)
        full = b1.PathModel(length, sites=(0, length-1), resource=2)
        selected = [k for k, (p, m, f) in enumerate(one.states)
                    if f == 1 or max(p[0], m[0]) <= cutoff]
        hcut = one.H[selected][:, selected]
        rows = []
        for il in selected:
            pl, ml, fl = one.states[il]
            for ir in selected:
                pr, mr, fr = one.states[ir]
                p = tuple(sorted(pl+tuple(length-1-y for y in pr)))
                m = tuple(sorted(ml+tuple(length-1-y for y in mr)))
                fields = b1.field(p, m, length)
                assert max(map(abs, fields), default=0) <= 1
                assert fields[length//2-1] == 0
                rows.append(full.index[(p, m, fl+2*fr)])
        assert len(set(rows)) == len(rows)
        dim = len(selected)
        embed = coo_matrix((np.ones(len(rows)), (rows, range(len(rows)))),
                           shape=(full.dim, len(rows))).tocsr()
        hprod = (kron(hcut, eye(dim), format='csr')+
                 kron(eye(dim), hcut, format='csr'))
        compression_defect = max_entry(embed.T@full.H@embed-hprod)
        assert compression_defect < 1e-12
        assert max_entry(embed.T@embed-eye(dim*dim, format='csr')) == 0
        gamma = full.H@embed-embed@hprod
        is_boundary = np.array([f == 0 and max(p[0], m[0]) == cutoff
                                for p, m, f in (one.states[i] for i in selected)])
        boundary_product = np.array([a or b for a in is_boundary for b in is_boundary])
        assert max_entry(gamma[:, ~boundary_product]) < 1e-12
        assert max_entry(embed.T@gamma) < 1e-12
        row_sum = float(np.max(np.asarray(abs(gamma).sum(axis=1))))
        col_sum = float(np.max(np.asarray(abs(gamma).sum(axis=0))))
        assert max(row_sum, col_sum) <= 2*(b1.JP+b1.JM)+1e-12
        # A fixed algebra fixture, NOT time-evolved and NOT an output trial.
        v = np.arange(1, one.dim+1, dtype=float).astype(complex)
        v += 1j*(np.arange(one.dim) % 3-1)
        v /= np.linalg.norm(v)
        cut = v[selected]
        derivative_cut = (one.H@v)[selected]
        phi = embed@np.kron(cut, cut)
        direct_defect = (embed@(np.kron(derivative_cut, cut)+
                               np.kron(cut, derivative_cut))-full.H@phi)
        incoming = derivative_cut-hcut@cut
        decomposed = (embed@(np.kron(incoming, cut)+np.kron(cut, incoming))-
                      gamma@np.kron(cut, cut))
        decomposition_defect = float(np.linalg.norm(direct_defect-decomposed))
        assert decomposition_defect < 1e-12
        shell = np.array([f == 0 and max(p[0], m[0]) in (cutoff, cutoff+1)
                          for p, m, f in one.states])
        shell_norm = float(np.linalg.norm(v[shell]))
        bound = 4*(b1.JP+b1.JM)*shell_norm
        observed = float(np.linalg.norm(direct_defect))
        assert observed <= bound+1e-12
        # Diagnose the tempting but wrong newly terminated path compression.
        naive = b1.PathModel(cutoff+1, sites=(0,), resource=1)
        permutation = [naive.index[one.states[i]] for i in selected]
        naive_defect = max_entry(hcut-naive.H[permutation][:, permutation])
        assert abs(naive_defect-b1.SELF) < 1e-12
        output.append(dict(length=length, cutoff=cutoff, physical_dimension=full.dim,
                           one_body_control_dimension=one.dim,
                           cut_one_body_dimension=dim, product_dimension=dim*dim,
                           compression_entry_defect=compression_defect,
                           gamma_row_max=row_sum, gamma_column_max=col_sum,
                           gamma_away_from_boundary_defect=max_entry(gamma[:, ~boundary_product]),
                           derivative_decomposition_defect=decomposition_defect,
                           algebra_fixture_defect_norm=observed,
                           algebra_fixture_bound=bound,
                           naive_new_endpoint_energy_defect=naive_defect,
                           note='No time evolution or operating probabilities.'))
    return output


def main():
    start = time.perf_counter()
    for name in ('SOURCE_PINS.json', 'EXECUTED_SOURCES.json'):
        for path, value in json.loads((HERE/name).read_text())['sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == value
    scalar = scalar_checks()
    matrix = embeddings(load_b1())
    data = {'exact_scalar_checks': scalar, 'matrix_embedding_controls': matrix,
            'status': 'PASS: rational inequalities and embedding identities only',
            'scaling_operating_gate': 'UNRESOLVED: source bounds not evaluated',
            'source_and_proof_pins_unchanged': True,
            'scope': 'No propagation. Double-precision matrix checks, not interval certification.'}
    out = Path(sys.argv[1]) if len(sys.argv)>1 else HERE
    out.mkdir(parents=True, exist_ok=True)
    (out/'CONTROLS.json').write_text(json.dumps(data, indent=2)+'\n')
    metrics = {'program_seconds': time.perf_counter()-start,
               'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
               'python': sys.version, 'numpy': np.__version__,
               'threads': {k: os.environ.get(k) for k in
                           ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')}}
    (out/'EXECUTION.json').write_text(json.dumps(metrics, indent=2)+'\n')
    print(json.dumps({'status': data['status'], 'metrics': metrics}))


if __name__ == '__main__':
    main()
