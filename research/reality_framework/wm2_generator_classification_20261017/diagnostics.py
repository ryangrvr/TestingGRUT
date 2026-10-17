"""Finite IG1 diagnostics, not formation workers or floating-point certificates.

Uses the frozen UC1 physical-kernel builder read-only. Its receiver has two
modes for these controls, NOT the S9 physical receiver. Independent prefix
and qubit controls are implemented below.
"""
from pathlib import Path
from itertools import product
from math import ceil, exp, factorial, pi
import importlib.util
import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

ROOT = Path(__file__).resolve().parent.parent
UPSTREAM = ROOT / 'wm2_charge_constraint_20261017' / 'operators.py'
spec = importlib.util.spec_from_file_location('frozen_uc1_operators', UPSTREAM)
uc1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(uc1)


def norm(x):
    return float(np.linalg.norm(x))


def maxentry(x):
    return float(np.max(np.abs(x.data), initial=0)) if hasattr(x, 'data') else float(np.max(np.abs(x), initial=0))


def comm(a, b):
    return a @ b - b @ a


def response(U, psi, S, O):
    y, tangent = U @ psi, -1j * (U @ (S @ psi))
    return float(2 * np.vdot(y, O @ tangent).real)


def tail(k, x):
    # Positive sum, no exp(x)-partial_sum cancellation.
    term = x ** k / factorial(k)
    answer = term
    for n in range(k + 1, 10000):
        term *= x / n
        answer += term
        if term <= max(1e-300, answer * 1e-17):
            return answer
    raise RuntimeError('tail did not converge')


def qubits():
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]], complex)
    Z = np.diag([1., -1.])
    S = np.kron(X / 2, I)
    O = np.kron(I, Y)
    plus_y = np.array([1, 1j]) / np.sqrt(2)
    zero = np.array([1., 0.])
    plus_x = np.array([1., 1.]) / np.sqrt(2)
    out = {}
    for name, B, initial_B, sign in [('internal_X', X, zero, -1),
                                    ('diagonal_Z', Z, plus_x, 1)]:
        H = .25 * np.kron(Z, B)
        psi = np.kron(plus_y, initial_B)
        errors = []
        for t in [.013, .7, pi]:
            chi = response(expm(-1j * t * H), psi, S, O)
            errors.append(abs(chi - sign * np.sin(.5 * t)))
        double = np.linalg.norm(comm(2 * S, comm(H, O)), 2)
        # A readout commuting with H provides a hypothesis boundary control.
        commuting_readout = np.kron(I, B)
        zero_chi = response(expm(-1j * pi * H), psi, S, commuting_readout)
        assert max(errors) < 2e-14 and abs(double - 1) < 2e-14
        assert abs(zero_chi) < 2e-14
        out[name] = {'maximum_exact_response_error': max(errors),
                     'normalized_double_commutator_norm': float(double),
                     'commuting_readout_response': zero_chi}
    out['source_variance'] = .25
    return out


def restricted_kernel(model):
    """Replace UC1's explicit H_int by its occupied-bond version only.

    All other physical matrix entries remain exactly those of frozen UC1.
    This is a witness for the restricted class, not a new formation result.
    """
    rows, cols, bad_values, bond_values = [], [], [], []
    for col, (word, q, n0, a, b) in enumerate(model.states):
        if a is None or b is None or a[0] == b[0] or b[1] == 2:
            continue
        rows.append(model.index[(word, q, n0, a, (b[0], 1 - b[1]))])
        cols.append(col)
        value = .25 * (1 if a[1] == 0 else -1)
        bad_values.append(value)
        bond_values.append(value if abs(word.index(a[0]) - word.index(b[0])) == 1 else 0.)
    bad = coo_matrix((bad_values, (rows, cols)), shape=model.H.shape).tocsr()
    bond = coo_matrix((bond_values, (rows, cols)), shape=model.H.shape).tocsr()
    H = (model.H - bad + bond).tocsr()
    assert maxentry(H - H.getH()) == 0
    return H, bond, bad


def common_reference(model, H, prefix):
    """On K: fixed-prefix matter, assisted tail, D; no bare geometry.

    On K-perp: the compression of the full H. Keeps the common physical
    Hilbert space. Removes entries explicitly, NOT by renormalizing states.
    """
    ell = len(prefix)
    inside = np.array([w[:ell] == prefix and
                       all(p is None or p[0] in prefix for p in (a, b))
                       for w, q, n0, a, b in model.states])
    h = H.tocoo()
    rr, cc, vv = [], [], []
    removed_bare = 0
    for row, col, value in zip(h.row, h.col, h.data):
        kr, kc = inside[row], inside[col]
        if kr != kc:
            continue
        if not kc:
            rr.append(row); cc.append(col); vv.append(value)
            continue
        wr, qr, nr, ar, br = model.states[row]
        wc, qc, nc, ac, bc = model.states[col]
        if wr != wc and qr == qc:
            removed_bare += 1
            continue
        # Remaining same-word matter entries are entirely within prefix,
        # because both endpoints are in K and the concrete range is one.
        rr.append(row); cc.append(col); vv.append(value)
    ref = coo_matrix((vv, (rr, cc)), shape=H.shape).tocsr()
    assert maxentry(ref - ref.getH()) == 0
    return ref, inside, removed_bare


def nearby(model, H, species):
    word = tuple(range(model.N))
    q = (1 << model.N) - 1
    level = 0 if species == 'A' else 2
    psi = np.zeros(H.shape[0], complex)
    source = np.zeros(H.shape[0])
    receiver = np.zeros(H.shape[0])
    for col, (w, bits, n0, a, b) in enumerate(model.states):
        carrier = a if species == 'A' else b
        if carrier is not None:
            source[col] = carrier[0] == word[0]
            receiver[col] = carrier[0] == word[1]
    for pos in word[:2]:
        a = (pos, level) if species == 'A' else None
        b = (pos, level) if species == 'B' else None
        psi[model.index[(word, q, 0, a, b)]] = 1 / np.sqrt(2)
    t = 1e-4
    columns = expm_multiply(-1j * t * H, np.column_stack((psi, -1j * source * psi)))
    chi = float(2 * np.vdot(columns[:, 0], receiver * columns[:, 1]).real)
    j = uc1.JA if species == 'A' else uc1.JB
    lower = j * t - tail(2, 3.4 * t)
    assert chi >= lower - 2e-14 and lower > 9e-7
    return {'chi': chi, 'analytic_lower': lower, 'source_variance': .25}


def full_kernel_controls():
    out = {}
    for N in [2, 3, 4]:
        model = uc1.KernelModel(N)
        H, bond, bad = restricted_kernel(model)
        vacuum = np.array([a is None and b is None for w, q, n0, a, b in model.states])
        single = np.array([a is None or b is None for w, q, n0, a, b in model.states])
        assert maxentry((H - model.H)[:, vacuum]) == 0
        assert maxentry(bond[:, single]) == 0
        result = {'dimension': H.shape[0], 'nonzeros': H.nnz,
                  'vacuum_difference_from_frozen_UC1': maxentry((H - model.H)[:, vacuum]),
                  'single_species_bond_action': maxentry(bond[:, single])}
        if N == 3:
            result['nearby_A'] = nearby(model, H, 'A')
            result['nearby_B'] = nearby(model, H, 'B')
        if N == 4:
            prefix = (0, 1, 2)
            ref, inside, count = common_reference(model, H, prefix)
            word = (0, 1, 2, 3)
            psi = np.zeros(H.shape[0], complex)
            # Both carriers at root, product internal states, all three
            # allowed quanta at contact to exercise the shortest pop word.
            for a_level, amplitude in [(0, 1), (1, 1j)]:
                psi[model.index[(word, 15, 3, (0, a_level), (0, 0))]] = amplitude / np.sqrt(2)
            S = model.source(0)
            O = model.receiver(2)
            difference = H - ref
            initial_drive = norm(difference @ psi)
            next_drive = norm(difference @ (ref @ psi))
            assert initial_drive < 2e-14 and next_drive > 0
            assert norm((~inside) * psi) == 0
            t = .02
            actual = expm_multiply(-1j * t * H, np.column_stack((psi, -1j * (S @ psi))))
            reference = expm_multiply(-1j * t * ref, np.column_stack((psi, -1j * (S @ psi))))
            err = norm(actual[:, 0] - reference[:, 0])
            chi = float(2 * np.vdot(actual[:, 0], O @ actual[:, 1]).real)
            chi_ref = float(2 * np.vdot(reference[:, 0], O @ reference[:, 1]).real)
            delta = 2.5 * tail(2, 2 * t)  # k=1
            loss = norm((~inside) * actual[:, 0]) ** 2
            bound_chi = tail(2, t) + 2 * delta  # ||S||=1/2, ||O||=1, w<=1/2
            assert err <= delta and loss <= delta ** 2
            assert abs(chi - chi_ref) <= 2 * delta
            assert abs(chi) <= bound_chi
            result['common_space_reference'] = {
                'prefix': prefix, 'm': 4, 'b': 0, 'R': 1, 'k': 1,
                'removed_bare_matrix_entries': count,
                'first_drive_norm': initial_drive,
                'next_drive_norm': next_drive,
                'time': t, 'state_error': err, 'state_error_bound': delta,
                'full_chi': chi, 'reference_chi': chi_ref,
                'response_bound': bound_chi,
                'K_loss_probability': loss, 'K_loss_probability_bound': delta ** 2}
        out[str(N)] = result
    return out


def prefix_control(length, R):
    """Independent fixed-prefix two-species implementation, no UC1 builder.

    Six internal combinations per ordered pair of carrier positions.
    Range-R coherent hopping plus finite-range Z_A X_B interactions.
    """
    states = [(a, la, b, lb) for a, b in product(range(length), repeat=2)
              for la, lb in product(range(2), range(3))]
    index = {s: k for k, s in enumerate(states)}
    rows, cols, values = [], [], []
    clock = np.array([float(lb == 2) for a, la, b, lb in states])
    for col, (a, la, b, lb) in enumerate(states):
        for species, position, j in [('A', a, 1 / 128 / R), ('B', b, 1 / 192 / R)]:
            for shift in range(-R, R + 1):
                if shift == 0 or not 0 <= position + shift < length:
                    continue
                target = (position + shift, la, b, lb) if species == 'A' else (a, la, position + shift, lb)
                rows.append(index[target]); cols.append(col); values.append(-j)
        if 0 < abs(a - b) <= R and lb in [0, 1]:
            rows.append(index[(a, la, b, 1 - lb)]); cols.append(col)
            values.append(.125 * (1 if la == 0 else -1))
    V = coo_matrix((values, (rows, cols)), shape=(len(states),) * 2).toarray()
    D = np.diag(clock)
    H = D + V
    assert np.max(abs(H - H.conj().T)) == 0
    Vnorm = float(np.linalg.norm(V, 2))
    assert Vnorm < .5
    source_event, receiver_event = 0, length - 1
    S = np.zeros_like(H, complex); O = np.zeros_like(H, complex)
    for col, (a, la, b, lb) in enumerate(states):
        if a == source_event:
            S[index[(a, 1 - la, b, lb)], col] = .5
        if b == receiver_event and lb in [0, 1]:
            O[index[(a, la, b, 1 - lb)], col] = 1j if lb == 0 else -1j
    a_order = ceil((length - 1) / R)
    current = O.copy()
    small = []
    # On-site D does not change support; these finite Taylor zeros must hold too.
    for n in range(a_order):
        value = float(np.max(np.abs(comm(S, current))))
        small.append(value)
        assert value < 1e-13
        current = comm(H, current)
    t = .03
    U = expm(-1j * t * H)
    comm_norm = float(np.linalg.norm(comm(S, U.conj().T @ O @ U), 2))
    upper = tail(a_order, 2 * Vnorm * t)
    assert comm_norm <= upper + 1e-13
    # A forbidden, occupation-tagged distant product must fail the same test.
    bad = np.zeros_like(H, complex)
    for col, (aa, la, bb, lb) in enumerate(states):
        if aa == source_event and bb == receiver_event and lb in [0, 1]:
            bad[index[(aa, la, bb, 1 - lb)], col] = .25 * (1 if la == 0 else -1)
    immediate = float(np.linalg.norm(comm(S, comm(bad, O)), 2))
    Ubad = expm(-1j * t * (H + bad))
    bad_norm = float(np.linalg.norm(comm(S, Ubad.conj().T @ O @ Ubad), 2))
    assert immediate > .49 and bad_norm > .01 and bad_norm > upper
    return {'length': length, 'R': R, 'dimension': len(states),
            'changing_norm': Vnorm, 'first_permitted_commutator_order': a_order,
            'low_order_max_entry_norms': small, 'time': t,
            'complete_commutator_norm': comm_norm, 'Dyson_bound': upper,
            'floating_comparison_tolerance': 1e-13,
            'forbidden_product_double_commutator': immediate,
            'forbidden_product_evolved_commutator': bad_norm}


def configuration_count(N, m, ell, b, R):
    """Bounded-depth, label-history-retaining support overapproximation.

    Storage and receiver are allowed every support transition, so omitting
    their amplitudes cannot hide a shorter path. Counts are not probabilities.
    """
    prefix = tuple(range(ell))
    initial = (tuple(range(m)), 0, b)
    c = max(R, 2)
    k = min(m - ell, ceil((ell - c - b) / R))
    frontier = {initial}
    visited = {initial}
    layer_sizes, first_possible = [], None
    for n in range(k + 1):
        layer_sizes.append(len(frontier))
        for word, a, bpos in frontier:
            # Conservative superset of every omitted-action support.
            trigger = len(word) == ell or max(word.index(a), word.index(bpos)) >= ell - c
            if trigger and first_possible is None:
                first_possible = n
            if n < k:
                assert not trigger
            assert word[:ell] == prefix
        if n == k:
            break
        new = set()
        for word, a, bpos in frontier:
            # Reference assisted tail operations, all intermediate label words.
            for v in range(N):
                if v not in word:
                    new.add((word + (v,), a, bpos))
            if len(word) > ell:
                new.add((word[:-1], a, bpos))
            for species, pos in [('A', a), ('B', bpos)]:
                place = word.index(pos)
                for shift in range(-R, R + 1):
                    newplace = place + shift
                    if shift and 0 <= newplace < ell:
                        new.add((word, word[newplace], bpos) if species == 'A' else (word, a, word[newplace]))
            # Simultaneous two-carrier changes within an interval of diameter R.
            for lo in range(max(0, ell - R)):
                interval = prefix[lo:lo + R + 1]
                if a in interval and bpos in interval:
                    for aa, bb in product(interval, repeat=2):
                        new.add((word, aa, bb))
        frontier = new
        visited.update(new)
    return {'N': N, 'm': m, 'ell': ell, 'b': b, 'R': R, 'k': k,
            'layer_sizes': layer_sizes, 'visited_configurations': len(visited),
            'first_possible_omitted_support': first_possible}
