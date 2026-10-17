"""Adversarial operator controls, not formation trajectories or proof certification."""
from itertools import product
from math import sqrt, pi
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import hashlib
import json
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from operators import KernelModel, electric, constraint_values, pairs, bonds, field_difference, JA, JB, KAPPA

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def comm_norm(H, diagonals):
    a = H.tocoo()
    return max(float(np.sqrt(np.sum(abs(a.data)**2 * (d[a.row] - d[a.col])**2))) for d in diagonals)


def old_witnesses():
    rows = []
    N, w, kappa = 3, (0, 1, 2), .25
    links = pairs(N)
    for charge in [1, 2]:
        basis = [(p, e) for p in range(N) for e in product(range(-3, 4), repeat=len(links))]
        index = {s: k for k, s in enumerate(basis)}
        G, C = [], []
        for p, fields in basis:
            # A charge-q carrier with charge-q root compensator.
            divergence = [0] * N
            for (i, j), e in zip(links, fields):
                divergence[i] += e; divergence[j] -= e
            divergence[p] -= charge; divergence[0] += charge
            G.append(divergence)
            C.append([e*e if edge not in bonds(w) else 0 for edge, e in zip(links, fields)])
        G, C = np.array(G).T, np.array(C).T
        physical = np.all(G == 0, axis=0) & np.all(C == 0, axis=0)
        assert physical.sum() == N
        def build(kind):
            rr, cc, vv = [], [], []
            for col, (p, e) in enumerate(basis):
                for j in range(N):
                    if j == p:
                        continue
                    isstar = p == 0 or j == 0
                    if kind.startswith('star') and not isstar:
                        continue
                    z = -kappa / sqrt(N - 1) if kind.startswith('star') else -kappa / (N - 1)
                    edge = tuple(sorted((p, j))); link = links.index(edge)
                    new = list(e)
                    if 'dress' in kind or kind == 'local':
                        new[link] += -charge if p < j else charge
                        if abs(new[link]) > 3:
                            continue
                    if kind == 'local':
                        if edge not in bonds(w):
                            continue
                        z = -JA
                    rr.append(index[(j, tuple(new))]); cc.append(col); vv.append(z)
            return coo_matrix((vv, (rr, cc)), shape=(len(basis),)*2).tocsr()
        local = build('local')
        assert comm_norm(local, G) == comm_norm(local, C) == 0
        for name in ['allpairs', 'star', 'allpairs_dress', 'star_dress']:
            H = build(name)
            gcomm, ccomm = comm_norm(H, G), comm_norm(H, C)
            if 'dress' in name:
                assert gcomm == 0 and ccomm > .01
            else:
                assert gcomm > .01
            assert comm_norm(H + local, G) == gcomm
            assert comm_norm(H + local, C) == ccomm
            P = np.flatnonzero(physical)
            block = H[P][:, P]
            if 'dress' not in name:
                assert block.nnz == 0
            rootcolumn = index[(0, (0,)*len(links))]
            z = H[:, rootcolumn].toarray().ravel()
            outside = float(np.linalg.norm(z[~physical]))
            rows.append(dict(charge=charge, witness=name, dimension=len(basis),
                             Gauss_commutator_frobenius=gcomm, flux_support_commutator_frobenius=ccomm,
                             physical_root_leakage_norm=outside,
                             local_addition_does_not_cancel=True,
                             projected_physical_entries=block.nnz))
    return rows


def kernel_checks():
    path = ROOT / 'research/reality_framework/wm2_prepared_disorder_20261016/model.py'
    spec = importlib.util.spec_from_file_location('frozen_s4', path)
    old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    rows = []
    models = {}
    for N in [2, 3, 4]:
        model = KernelModel(N)
        models[N] = model
        reference = old.vacuum(N, 2, model.delta, g=model.g)
        ids = [model.index[(w, sum(1 << v for v in w), o.count(0), None, None)] for w, o in reference['states']]
        block = model.H[ids][:, ids]
        difference = block - reference['H']
        matrix_error = float(np.max(abs(difference.data), initial=0))
        assert matrix_error < 5e-14
        # The whole seeded vacuum component is invariant, not only its diagonal block.
        others = np.ones(len(model.states), bool); others[ids] = False
        assert not model.H[np.flatnonzero(others)][:, ids].nnz
        assert model.audit['maximum_electric_magnitude'] == 3
        assert model.audit['checked_charge_hops'] > 0
        assert model.audit['blocked_nonzero_flux_pops'] > 0
        row = dict(N=N, complete_kernel_dimension=len(model.states), nonzeros=model.H.nnz,
                   delta=list(model.delta), common_delta=.23,
                   formation_block_dimension=len(ids), formation_block_max_difference=matrix_error,
                   formation_component_leakage=0, **model.audit)
        rows.append(row)
    # Integer field covariance, including orientation reversals and nonzero roots.
    rng = np.random.default_rng(48123)
    for _ in range(40):
        N = 6; perm = rng.permutation(N)
        w = (0, 1, 2, 3, 4); pa, pb = 2, 4
        e = electric(N, w, pa, pb)
        neww = tuple(int(perm[p]) for p in w)
        mapped = {}
        for (i, j), x in zip(pairs(N), e):
            u, v = int(perm[i]), int(perm[j])
            mapped[tuple(sorted((u, v)))] = x if u < v else -x
        expected = tuple(mapped[edge] for edge in pairs(N))
        actual = electric(N, neww, int(perm[pa]), int(perm[pb]))
        assert expected == actual
        assert not any(sum((list(z) for z in constraint_values(N, neww, int(perm[pa]), int(perm[pb]), actual)), []))
    return rows, models


def weighted_response(model, w, prefix, receiver, t):
    psi = model.prepared(w, receiver)
    source = model.source(w[0])
    O = model.receiver(receiver, prefix)
    Olocal = model.receiver(receiver)
    tangent = -1j * (source @ psi)
    assert abs(np.vdot(psi, source @ psi)) < 1e-13
    variance = float(np.vdot(source @ psi, source @ psi).real)
    assert abs(variance - .25) < 1e-14
    x = np.column_stack([psi, tangent])
    out = []
    for label, H in [('full', model.H), ('zero_mobility_reference', model.H - model.Hmove)]:
        evolved = expm_multiply(-1j * t * H, x)
        chi = 2 * np.vdot(evolved[:, 0], O @ evolved[:, 1]).real
        pfx = np.array([state[0][:len(prefix)] == prefix for state in model.states])
        lost = float(np.linalg.norm(evolved[~pfx, 0])**2)
        local_chi = 2 * np.vdot(evolved[:, 0], Olocal @ evolved[:, 1]).real
        out.append((label, float(chi), lost, float(local_chi)))
    expected = -np.sin(2 * KAPPA * t)
    assert abs(out[1][1] - expected) < 2e-11
    assert abs(out[1][3] - expected) < 2e-11
    bound = 2 * (2*JA + 2*JB) * abs(t)
    assert abs(out[0][1] - expected) <= bound + 2e-11
    assert abs(out[0][3] - expected) <= bound + 2e-11
    assert out[0][2] <= ((2*JA + 2*JB)*t)**2 + 2e-11
    # Independent +/- phase propagation, not reuse of the tangent formula.
    step = 2e-4
    values, local_values = [], []
    for theta in [-step, step]:
        initial = expm_multiply(-1j * theta * source, psi)
        evolved = expm_multiply(-1j * t * model.H, initial)
        values.append(float(np.vdot(evolved, O @ evolved).real))
        local_values.append(float(np.vdot(evolved, Olocal @ evolved).real))
    finite_difference = (values[1] - values[0]) / (2 * step)
    assert abs(finite_difference - out[0][1]) < 2e-8
    local_difference = (local_values[1] - local_values[0])/(2*step)
    assert abs(local_difference-out[0][3]) < 2e-8
    return dict(N=model.N, receiver_modes=2, initial_word=list(w), retained_prefix=list(prefix),
                receiver_root_distance=w.index(receiver), time=t, source_variance=variance,
                full_success_weighted_response=out[0][1], zero_mobility_response=out[1][1],
                full_unconditioned_local_response=out[0][3],
                exact_reference_formula=float(expected), full_prefix_loss_probability=out[0][2],
                response_difference_bound=bound, finite_phase_difference=finite_difference,
                local_finite_phase_difference=local_difference,
                scope='small full physical-kernel diagnostic; no physical-L formation trajectory')


def flux_word_counts():
    rows = []
    # Disjoint positive charges 1 and 2 cannot cancel each other's cut change
    # with at most one particle of each species.
    for q in [1, 2]:
        other = 3 - q
        for delta_other in [-other, 0, other]:
            assert q + delta_other != 0
    for m in [6, 9, 14]:
        N = m + 2; w = tuple(range(m))
        for species in ['A', 'B']:
            for initial in [0, 1]:
                for final in range(initial + 1, m):
                    for other_start in range(m):
                        for other_end in range(m):
                            pa, pb = (initial, other_start) if species == 'A' else (other_start, initial)
                            qa, qb = (final, other_end) if species == 'A' else (other_end, final)
                            before = electric(N, w, pa, pb); after = electric(N, w, qa, qb)
                            assert len(field_difference(before, after)) >= final - initial
        rows.append(dict(word_length=m, charges=[1, 2], all_two_carrier_endpoint_assignments_checked=True))
    assert Q(2) * (Q(2, 64) + Q(2, 96)) * 4 == Q(5, 12)
    assert ((Q(2, 64) + Q(2, 96)) * 4)**2 == Q(25, 576) < Q(1, 16)
    return rows


def nearby_response(model):
    rows = []
    w, t = (0, 1, 2), 1e-4
    q = sum(1 << v for v in w)
    for species, level, hopping in [('A', 0, JA), ('B', 2, JB)]:
        x = np.zeros(len(model.states), complex)
        source = np.zeros(len(model.states))
        receiver = np.zeros(len(model.states))
        for col, (word, mask, n0, a, b) in enumerate(model.states):
            carrier = a if species == 'A' else b
            if carrier is not None:
                source[col] = int(carrier[0] == 0)
                receiver[col] = int(carrier[0] == 1)
        for p in [0, 1]:
            for n0, phase in [(0, 1), (len(w)-1, 1j)]:
                a, b = ((p, level), None) if species == 'A' else (None, (p, level))
                x[model.index[(w, q, n0, a, b)]] = phase/2
        assert abs(np.linalg.norm(x)-1) < 1e-14
        variance = float(np.vdot(x, source*source*x).real - np.vdot(x, source*x).real**2)
        assert abs(variance-.25) < 1e-14
        evolved = expm_multiply(-1j*t*model.H, np.column_stack([x,-1j*source*x]))
        chi = float(2*np.vdot(evolved[:,0], receiver*evolved[:,1]).real)
        u = Q(17, 50000)  # 3.4 t
        remainder = u*u/(2*(1-u))
        assert remainder < Q(6,100000000)
        lower = Q(1,64 if species=='A' else 96)*Q(1,10000)-remainder
        assert lower > Q(9,10000000)
        assert chi > float(lower) and abs(chi-hopping*t) < 6e-8
        rows.append(dict(species=species, charge=1 if species=='A' else 2,
                         occupied_internal_level=level, time=t, source_variance=variance,
                         full_nearby_response=chi, uniform_analytic_lower=float(lower)))
    return rows


def run():
    start = time.monotonic()
    witnesses = old_witnesses()
    kernels, models = kernel_checks()
    dynamics = [weighted_response(models[3], (0, 1, 2), (0, 1, 2), 2, pi),
                weighted_response(models[4], (0, 1, 2), (0, 1, 2), 2, .7)]
    counts = flux_word_counts()
    nearby = nearby_response(models[3])
    return dict(status='PASS_DIAGNOSTICS_FULL_GATE_NEGATIVE', witnesses=witnesses,
                kernel_isometry=kernels, information_channel=dynamics, cut_change_checks=counts,
                nearby_matter_response=nearby,
                frozen_S4_implementation_sha256=hashlib.sha256((ROOT / 'research/reality_framework/wm2_prepared_disorder_20261016/model.py').read_bytes()).hexdigest(),
                elapsed_seconds=time.monotonic()-start,
                arithmetic='ordinary double precision diagnostics and exact integer/rational checks; not outward rounding or external proof review')


if __name__ == '__main__':
    for name in ['checks.py', 'operators.py']:
        assert (HERE / name).is_file()
    r = run()
    (HERE / 'CHECKS.json').write_text(json.dumps(r, indent=2) + '\n')
    (HERE / 'EXECUTED_SOURCES.json').write_text(json.dumps({name:hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in ['checks.py','operators.py']}, indent=2)+'\n')
    print(json.dumps(r, indent=2))
