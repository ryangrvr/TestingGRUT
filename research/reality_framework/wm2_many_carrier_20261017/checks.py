"""Frozen D1 finite diagnostics; no physical formation trajectories."""
from pathlib import Path
from fractions import Fraction
from itertools import product
from math import ceil, exp, factorial
import hashlib
import importlib.util
import json
import resource
import time
import numpy as np
from scipy.linalg import expm, eigvalsh
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, expm_multiply
from model import (FockChain, KernelModel, electric, constraints, pairs,
                   common_reference, JA, JB, KAPPA)


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def tail(k, x):
    # Positive series; do not subtract two nearly equal exponential values.
    term = x ** k / factorial(k)
    result = term
    for n in range(k + 1, 10000):
        term *= x / n
        result += term
        if term <= max(result * 1e-17, 1e-300):
            return result
    raise RuntimeError('positive tail failed to converge')


def maxentry(a):
    return float(np.max(np.abs(a.data), initial=0))


def delta(t, a, L, R=1, J0=2, mu=1):
    C = 2 / (1 - exp(-mu))
    v = 2 * J0 * (exp(mu * R) - 1)
    return min(2., tail(a + 1, 1.1 * abs(t)) +
               6 * np.sqrt(C) / v * exp(-mu * L / 2) * np.expm1(v * abs(t) / 2))


def capacity():
    result = []
    for N in range(2, 6):
        word = tuple(range(N))
        count = moves = accepted = max_P = 0
        for slots in product(range(4), repeat=N):
            a = tuple(i for i, value in enumerate(slots) if value & 1)
            b = tuple(i for i, value in enumerate(slots) if value & 2)
            field = electric(N, word, a, b)
            g, c = constraints(N, word, a, b, field)
            assert not any(g) and not any(c)
            assert max(map(abs, field)) <= 3 * (N - 1)
            # The last-link flux is precisely minus the head charge.
            last = pairs(N).index((N - 2, N - 1))
            assert (field[last] == 0) == (not slots[-1])
            if max(map(abs, field)) <= 3:
                accepted += 1
                max_P = max(max_P, len(a) + len(b))
                assert len(a) + len(b) <= 5
            for positions, other, charge, is_A in [(a, b, 1, True), (b, a, 2, False)]:
                for old in positions:
                    for new in [old - 1, old + 1]:
                        if not 0 <= new < N or new in positions:
                            continue
                        changed = tuple(sorted(new if x == old else x for x in positions))
                        aa, bb = (changed, other) if is_A else (other, changed)
                        newfield = electric(N, word, aa, bb)
                        altered = tuple(k for k, (x, y) in enumerate(zip(field, newfield)) if x != y)
                        link = pairs(N).index(tuple(sorted((old, new))))
                        assert altered == (link,)
                        assert newfield[link] - field[link] == (-charge if old < new else charge)
                        moves += 1
            count += 1
        result.append({'N': N, 'configurations': count, 'hops_checked': moves,
                       'C3_accepted': accepted, 'C3_max_particles': max_P})
    assert result[-1]['C3_max_particles'] == 5
    prices = []
    for n in [1, 2, 4, 8, 16]:
        N = 32 * n
        a = tuple(4 * q for q in range(n))
        b = tuple(4 * q + 2 for q in range(n))
        charges = [int(x in a) + 2 * int(x in b) for x in range(4 * n)]
        remainder = sum(charges)
        cuts = []
        for q in charges:
            remainder -= q
            cuts.append(remainder)
        L1 = sum(cuts); L2 = sum(x * x for x in cuts)
        assert L1 == 6 * n * n - 2 * n
        assert L2 == 12 * n ** 3 - 6 * n * n + 2 * n
        prices.append({'N': N, 'particles': 2 * n, 'largest_flux': max(cuts),
                       'sum_abs_E': L1, 'sum_E_squared': L2})
    return {'exhaustive': result, 'density_resource_prices': prices}


def local_fock():
    # Independent 12x12 local-site tensor implementation, without matter_actions.
    I3, I4 = np.eye(3), np.eye(4)
    tensor_bond = np.zeros((144, 144))
    for levels, species, j in [(2, 'A', JA), (3, 'B', JB)]:
        for level in range(levels):
            local = np.zeros((levels + 1, levels + 1)); local[0, level + 1] = 1
            annih = np.kron(local, I4) if species == 'A' else np.kron(I3, local)
            hop = np.kron(annih.T, annih)
            tensor_bond -= j * (hop + hop.T)
    ZA = np.kron(np.diag([0., 1., -1.]), I4)
    xb = np.zeros((4, 4)); xb[1, 2] = xb[2, 1] = 1
    XB = np.kron(I3, xb)
    tensor_bond += KAPPA * (np.kron(ZA, XB) + np.kron(XB, ZA))
    # Full two-site hard-core space is the direct sum of these nine sectors.
    bond_max = 0.
    matrix_error = 0.
    for na, nb in product(range(3), repeat=2):
        chain = FockChain(2, na, nb)
        indices = []
        for a, b in chain.states:
            ad, bd = dict(a), dict(b)
            site = [(ad[x] + 1 if x in ad else 0) * 4 +
                    (bd[x] + 1 if x in bd else 0) for x in range(2)]
            indices.append(12 * site[0] + site[1])
        difference = chain.V.toarray() - tensor_bond[np.ix_(indices, indices)]
        matrix_error = max(matrix_error, float(np.max(abs(difference), initial=0)))
        assert matrix_error < 1e-14
        eig = eigvalsh(chain.V.toarray())
        bond_max = max(bond_max, float(np.max(abs(eig), initial=0)))
    analytic_bond = JA + JB + 2 * KAPPA
    assert bond_max <= analytic_bond + 1e-14
    norm_scaling = []
    for n in [1, 2]:
        chain = FockChain(2 * n, n, n)
        if len(chain.states) < 100:
            vals = eigvalsh(chain.V.toarray())
            actual = float(np.max(abs(vals)))
        else:
            actual = float(abs(eigsh(chain.V, k=1, which='LM', return_eigenvectors=False,
                                     v0=np.linspace(1., 2., len(chain.states)))[0]))
        lower = KAPPA * (2 * n - 1)
        upper = 2 * JA * n + 2 * JB * n + 2 * KAPPA * n
        assert lower - 2e-12 <= actual <= upper + 2e-12
        norm_scaling.append({'each_species': n, 'dimension': len(chain.states),
                             'actual_changing_norm': actual, 'linear_lower': lower,
                             'linear_upper': upper})
    chain = FockChain(5, 2, 1)
    H = chain.H.toarray()
    matter = []
    for species in ['A', 'B']:
        psi = np.zeros(len(chain.states), complex)
        if species == 'A':
            columns = [(((p, 0), (3, 0)), ((4, 0),)) for p in [0, 1]]
        else:
            columns = [(((3, 0), (4, 0)), ((p, 2),)) for p in [0, 1]]
        for state in columns:
            psi[chain.index[state]] = 1 / np.sqrt(2)
        source, receiver = chain.number(species, 0), chain.number(species, 1)
        t = 5e-5
        columns = expm_multiply(-1j * t * chain.H,
                               np.column_stack((psi, -1j * source * psi)))
        chi = float(2 * np.vdot(columns[:, 0], receiver * columns[:, 1]).real)
        S, O = np.diag(source), np.diag(receiver)
        first = float(-np.vdot(psi, (S @ (H @ O - O @ H) - (H @ O - O @ H) @ S) @ psi).real)
        j = JA if species == 'A' else JB
        assert abs(first - j) < 1e-14
        lower = j * t - 72 * t * t
        assert chi >= lower - 1e-14 and lower > 3e-7
        variance = float(np.vdot(psi, source * source * psi).real - np.vdot(psi, source * psi).real ** 2)
        assert abs(variance - .25) < 1e-14
        matter.append({'species': species, 'time': t, 'carriers': 3,
                       'source_variance': variance, 'first_response_coefficient': first,
                       'full_Fock_response': chi, 'local_analytic_lower': lower})
    weight_checks = []
    for mu in [.25, 1.]:
        F = np.array([sum(exp(mu * (p - 1)) for p, c in a + b) for a, b in chain.states])
        Cmat = 1j * (H * F[None, :] - F[:, None] * H)
        normalized = Cmat / np.sqrt(F[:, None] * F[None, :])
        actual = float(np.max(abs(eigvalsh(normalized))))
        bound = 4 * (exp(mu) - 1)
        assert actual <= bound + 1e-12
        weight_checks.append({'mu': mu, 'quadratic_form_rate_norm': actual,
                              'density_independent_bound': bound})
    # Local internal source and a distant internal readout on a many-carrier sector.
    X = np.array([[0., 1.], [1., 0.]])
    Y = np.array([[0., -1j, 0.], [1j, 0., 0.], [0., 0., 0.]])
    S = chain.internal('A', 0, X / 2).toarray()
    O = chain.internal('B', 4, Y).toarray()
    nested = O.copy(); zeros = []
    for n in range(4):
        val = float(np.max(abs(S @ nested - nested @ S)))
        zeros.append(val); assert val < 1e-12
        nested = H @ nested - nested @ H
    t = .01
    U = expm(-1j * t * H)
    evolved = U.conj().T @ O @ U
    signal_norm = float(np.linalg.norm(S @ evolved - evolved @ S, 2))
    connected_bound = exp(-4) * np.expm1(8 * exp(1) * t)
    assert signal_norm <= connected_bound + 1e-12
    return {'independent_two_site_matrix_max_error': matrix_error,
            'maximum_full_bond_norm': bond_max, 'analytic_bond_bound': analytic_bond,
            'per_site_upper': 2 * analytic_bond, 'global_norm_scaling': norm_scaling,
            'many_background_nearby': matter, 'weighted_form': weight_checks,
            'connected': {'dimension': len(chain.states), 'distance': 4, 'time': t,
                          'low_order_commutator_entries': zeros,
                          'actual_commutator_norm': signal_norm, 'local_bound': connected_bound}}


def full_geometry():
    spec = importlib.util.spec_from_file_location('frozen_uc1', HERE.parent / 'wm2_charge_constraint_20261017/operators.py')
    old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    vacuum_checks = []
    for N in [2, 3]:
        new = KernelModel(N, 0, 0)
        prior = old.KernelModel(N)
        indices = [prior.index[(word, q, n0, None, None)] for word, q, n0, a, b in new.states]
        difference = new.H - prior.H[indices][:, indices]
        assert maxentry(difference) < 1e-14
        vacuum_checks.append({'N': N, 'full_vacuum_dimension': len(new.states), 'maximum_difference': maxentry(difference)})
    controls = []
    for N in [3, 4]:
        model = KernelModel(N, 2, 1)
        result = {'N': N, 'N_A': 2, 'N_B': 1, 'dimension': len(model.states),
                  'nonzeros': model.H.nnz, 'audit': model.audit}
        if N == 4:
            prefix = (0, 1, 2)
            ref, inside = common_reference(model, prefix)
            psi = np.zeros(len(model.states), complex)
            psi[model.index[((0, 1, 2, 3), 15, 3, ((0, 0), (1, 1)), ((0, 0),))]] = 1
            t = .02
            actual = expm_multiply(-1j * t * model.H, psi)
            reference = expm_multiply(-1j * t * ref, psi)
            error = float(np.linalg.norm(actual - reference))
            account = delta(t, a=1, L=0)
            loss = float(np.linalg.norm((~inside) * actual) ** 2)
            mismatch = np.array([q != sum(1 << x for x in word) for word, q, n0, a, b in model.states])
            storage = float(np.linalg.norm(mismatch * actual) ** 2)
            assert error <= account and loss <= account ** 2 and storage <= account ** 2
            # Nonzero backreaction: the full trajectory is not storage-reflected.
            assert storage > 0 and error > 0
            result['common_reference'] = {'time': t, 'a': 1, 'L': 0,
                'scope': 'small no-macroscopic-buffer diagnostic, not the all-size density sequence',
                'state_error': error, 'state_error_bound': account,
                'prefix_loss_probability': loss, 'storage_error_probability': storage,
                'probability_error_bound': account ** 2}
        controls.append(result)
    return {'vacuum_isometry': vacuum_checks, 'complete_many_carrier_controls': controls}


def constants():
    # Rational coefficient checks in the all-N proof, using 2<e<3.
    assert Fraction(1, 4) - Fraction(33, 40960) > Fraction(1, 5)
    assert (Fraction(63, 1024) - Fraction(1, 17)) * 256 >= Fraction(1, 2)
    assert Fraction(1, 32) - Fraction(3, 512) >= Fraction(1, 40)
    rows = []
    for N in [256, 512, 1024, 2048]:
        a, L = N // 4, N // 8 - 1
        # Use logarithmic/exponential coarse head tail for large k, avoiding
        # underflow of x**k/k! in a nonessential evaluated comparison.
        t = N / 4096
        C = 2 / (1 - exp(-1)); v = 4 * (exp(1) - 1)
        head = exp(exp(1) * 1.1 * t - a - 1)
        evaluated = head + 6 * np.sqrt(C) / v * exp(-L / 2) * np.expm1(v * t / 2)
        coarse = exp(-N / 5) + 3 * exp(-N / 17)
        assert evaluated <= coarse
        nearby_delta = exp(exp(1) * 1.1 * 5e-5 - a - 1) + 6 * np.sqrt(C) / v * exp(-L / 2) * np.expm1(v * 5e-5 / 2)
        assert nearby_delta < 1e-9
        lower = JB * 5e-5 - 72 * (5e-5) ** 2 - 4 * nearby_delta
        assert lower > 3e-7
        rows.append({'N': N, 'particles': N // 16, 'window': t,
                     'evaluated_state_bound': evaluated, 'coarse_state_bound': coarse,
                     'nearby_state_bound': nearby_delta, 'both_species_response_lower': lower})
    return rows


def run():
    start = time.perf_counter()
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    assert hashlib.sha256((HERE / 'CHARTER.md').read_bytes()).hexdigest() == pins['charter_sha256']
    for path, expected in pins['files'].items():
        assert hashlib.sha256((REPO / path).read_bytes()).hexdigest() == expected
    result = {'status': 'RUNNING', 'precision': 'ordinary complex128; no outward rounding',
              'formation_trajectories': 0, 'capacity': capacity(), 'local_Fock': local_fock(),
              'full_geometry': full_geometry(), 'all_size_constants': constants()}
    elapsed = time.perf_counter() - start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    assert elapsed < 60 and peak < 1024
    for path, expected in pins['files'].items():
        assert hashlib.sha256((REPO / path).read_bytes()).hexdigest() == expected
    result.update(status='PASS_FINITE_DIAGNOSTICS', elapsed_seconds=elapsed,
                  peak_resident_MiB=peak, all_size_proof='AUTHOR_DERIVED_INDEPENDENT_REVIEW_PENDING',
                  populated_seed_formation='NOT_DERIVED_NUMBER_CONSERVATION_OBSTRUCTION')
    return result


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, allow_nan=False))
