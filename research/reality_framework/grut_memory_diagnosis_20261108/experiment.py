"""E2-01: finite-parent factorial diagnosis and full-auxiliary memory evolution."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

import numpy as np
import scipy
from scipy.linalg import eigh, expm, qr, svd

HERE = Path(__file__).resolve().parent
MATCH = HERE.parent / 'grut_detector_matching_20261107'
E1 = HERE.parent / 'grut_effective_instruments_20261107'
sys.path.insert(0, str(MATCH))
import bosonized
from continuum import key, tv

OLD = [[0., .7, 1.9], [0., .6, 1.5, 2.4]]
NEW = [[0., .43, 1.37, 2.61], [0., .31, 1.07, 1.83, 2.79]]
LEVELS = [32, 128, 512]


def save(path, data):
    with path.open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write('\n')


def fro(x):
    return float(np.linalg.norm(x))


def commutator_upper(ea, ec, b):
    delta = ec[:, None] - ea[None, :]
    d = delta * b
    aa = np.sqrt(2.) * fro(delta * d)
    top = b.conj().T @ d + d.conj().T @ b
    bottom = -b @ d.conj().T - d @ b.conj().T
    vv = np.hypot(fro(top), fro(bottom))
    return {'H0_double_commutator_Frobenius': float(aa),
            'V_double_commutator_Frobenius': float(vv),
            'coefficient': float(aa / 24. + vv / 12.)}


class MemorySplit:
    def __init__(self, ea, ec, b, energy_shift=0.):
        self.k = len(ea)
        self.e = np.concatenate((ea, ec)) - energy_shift
        self.u, self.s, self.vh = svd(b, full_matrices=False)

    def exchange(self, x, step):
        p, q = x[:self.k], x[self.k:]
        pc, qc = self.vh @ p, self.u.conj().T @ q
        c = (np.cos(self.s * step) - 1.)[:, None]
        s = (-1j * np.sin(self.s * step))[:, None]
        return np.vstack((p + self.vh.conj().T @ (c * pc + s * qc),
                          q + self.u @ (c * qc + s * pc)))

    def propagate(self, x, duration, steps):
        if duration == 0.:
            return x.copy()
        h = duration / steps
        phase = np.exp(-.5j * self.e * h)[:, None]
        y = x.copy()
        for _ in range(steps):
            y = phase * self.exchange(phase * y, h)
        return y


def instrument_checks(w):
    identity = np.eye(len(w))
    deficit = (identity - w.conj().T @ w) / 2.
    minimum = float(eigh(deficit, eigvals_only=True, subset_by_index=[0, 0])[0])
    if minimum < -1e-10:
        raise RuntimeError(f'Expansive detector: deficit {minimum}')
    plus, minus = (identity + w) / 2., (identity - w) / 2.
    return {'minimum_overflow_effect_eigenvalue': minimum,
            'completion_residual': fro(plus.conj().T @ plus + minus.conj().T @ minus + deficit - identity)}


def verify_small():
    rng = np.random.default_rng(20261108)
    ea, ec = np.array([-.3, .4]), np.array([.8, 1.2, 2., 2.6])
    b = .2 * (rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2)))
    h0 = np.diag(np.concatenate((ea, ec)))
    v = np.zeros((6, 6), dtype=complex)
    v[2:, :2] = b
    v[:2, 2:] = b.conj().T
    solver = MemorySplit(ea, ec, b)
    identity = np.eye(6, dtype=complex)
    x = rng.normal(size=(6, 3)) + 1j * rng.normal(size=(6, 3))
    exchange_error = fro(solver.exchange(x, .31) - expm(-1j * v * .31) @ x)
    cb = commutator_upper(ea, ec, b)
    first = h0 @ v - v @ h0
    direct_a = fro(h0 @ first - first @ h0)
    direct_v = fro(v @ first - first @ v)
    coefficient_error = abs(cb['coefficient'] - (direct_a / 24. + direct_v / 12.))
    comparisons = []
    for h in [.07, .31]:
        split = solver.propagate(identity, h, 1)
        exact = expm(-1j * (h0 + v) * h)
        error = float(np.linalg.norm(split - exact, ord=2))
        comparisons.append({'step': h, 'error': error, 'bound': cb['coefficient'] * h**3,
                            'unitarity_error': fro(split.conj().T @ split - identity)})
    zero = MemorySplit(ea, ec, np.zeros_like(b))
    commuting_error = fro(zero.propagate(x, .7, 9) - expm(-1j * h0 * .7) @ x)
    w = .8 * expm(1j * (v + h0))
    defect = (identity - w.conj().T @ w) / 2.
    ed, vd = eigh(defect)
    ko = (vd * np.sqrt(ed)) @ vd.conj().T
    kp, km = (identity + w) / 2., (identity - w) / 2.
    dilation_error = fro(kp.conj().T @ kp + km.conj().T @ km + ko.conj().T @ ko - identity)
    expansive_rejected = False
    try:
        instrument_checks(1.1 * identity)
    except RuntimeError:
        expansive_rejected = True
    passed = (exchange_error < 1e-11 and coefficient_error < 1e-11 and commuting_error < 1e-11
              and dilation_error < 1e-11 and expansive_rejected
              and all(z['error'] <= z['bound'] + 1e-12 and z['unitarity_error'] < 1e-11 for z in comparisons))
    return {'pass': bool(passed), 'exchange_vs_expm': exchange_error,
            'commutator_formula_error': coefficient_error, 'commuting_limit_error': commuting_error,
            'overflow_dilation_error': dilation_error, 'expansive_detector_rejected': expansive_rejected,
            'local_bounds': comparisons}


def evolution(e, vectors, shift):
    def propagate(x, duration):
        if duration == 0.:
            return x.copy()
        return vectors @ (np.exp(-1j * (e - shift) * duration)[:, None] * (vectors.conj().T @ x))
    return propagate


def history(propagate, w, initial, times, retained, b=None, reset=False, details=False):
    branches = initial[:, None].copy()
    labels = ['']
    previous = 0.
    losses, resets, stages = [], [], []
    for index, moment in enumerate(times):
        branches = propagate(branches, moment - previous)
        if details:
            rows = {}
            for n, label in enumerate(labels):
                p, q = branches[:retained, n], branches[retained:, n]
                weight = float(np.vdot(branches[:, n], branches[:, n]).real)
                leak = .5 * fro(w[retained:, :retained] @ p)
                returned = .5 * fro(w[:retained, retained:] @ q)
                residual = fro(b @ p) if b is not None else 0.
                rows[label] = {'weight_before': weight, 'Q_weight_before': float(np.vdot(q, q).real),
                               'detector_outgoing_Q_norm': leak, 'detector_return_P_norm': returned,
                               'dynamical_source_norm': residual,
                               'conditional_detector_outgoing_Q_norm': leak / np.sqrt(weight) if weight > 0 else 0.,
                               'conditional_dynamical_source_norm': residual / np.sqrt(weight) if weight > 0 else 0.}
            stages.append({'time': moment, 'prefixes': rows})
        shifted = w @ branches
        loss = (fro(branches)**2 - fro(shifted)**2) / 2.
        if loss < -1e-10:
            raise RuntimeError('Negative physical overflow')
        losses.append(loss)
        branches = np.concatenate(((branches + shifted) / 2., (branches - shifted) / 2.), axis=1)
        labels = [s + '+' for s in labels] + [s + '-' for s in labels]
        reset_loss = 0.
        if reset and index != len(times) - 1:
            reset_loss = fro(branches[retained:])**2
            branches[retained:] = 0.
        resets.append(reset_loss)
        previous = moment
    probabilities = {label: float(np.vdot(branches[:, n], branches[:, n]).real) for n, label in enumerate(labels)}
    probabilities['overflow'] = float(sum(losses) + sum(resets))
    norm_error = abs(sum(probabilities.values()) - 1.)
    if norm_error >= 1e-9 or min(probabilities.values()) < -1e-10:
        raise RuntimeError(f'Invalid history distribution: {norm_error}')
    return {'probabilities': probabilities, 'ordinary_overflow_by_step': losses,
            'reset_overflow_by_step': resets, 'normalization_error': norm_error,
            'prefix_diagnostics': stages}


def autonomous_bound(parent_propagate, diagonal, w, initial, times, retained, preparation_error):
    wr = w[:retained, :retained]
    x = initial[:retained, None].copy()
    labels, previous = [''], 0.
    d = preparation_error
    propagation_budget, instrument_budget = 0., 0.
    records = []
    for moment in times:
        dt = moment - previous
        embedded = np.zeros((len(w), x.shape[1]), dtype=complex)
        embedded[:retained] = x
        exact = parent_propagate(embedded, dt)
        x = np.exp(-1j * diagonal[:retained] * dt)[:, None] * x
        comparison = np.zeros_like(embedded)
        comparison[:retained] = x
        propagation = fro(exact - comparison)
        instrument = fro(w[retained:, :retained] @ x) / np.sqrt(2.)
        records.append({'time': moment, 'propagation_Frobenius_defect': propagation,
                        'instrument_stacked_defect': instrument,
                        'exact_propagation_by_prefix': dict(zip(labels, map(float, np.linalg.norm(exact - comparison, axis=0)))),
                        'detector_leak_by_prefix': dict(zip(labels, map(float, .5 * np.linalg.norm(w[retained:, :retained] @ x, axis=0))))})
        d += propagation + instrument
        propagation_budget += propagation
        instrument_budget += instrument
        shifted = wr @ x
        x = np.concatenate(((x + shifted) / 2., (x - shifted) / 2.), axis=1)
        labels = [s + '+' for s in labels] + [s + '-' for s in labels]
        previous = moment
    return {'bound': min(1., d * (1. + fro(x))), 'initial_vector_error': preparation_error,
            'propagation_budget': propagation_budget, 'instrument_budget': instrument_budget,
            'ordinary_stack_vector_bound': d, 'final_ordinary_norm': fro(x), 'stages': records,
            'scope': 'Exact finite-parent propagator used to evaluate branch-relevant defects; not computational compression.'}


def factor_attribution(rows):
    distributions = {s: r['probabilities'] for s, r in rows.items()}
    p00, p10, p01, p11 = [distributions[s] for s in ['00', '10', '01', '11']]
    terms = {'propagation': {s: p10[s] - p00[s] for s in p00},
             'instrument': {s: p01[s] - p00[s] for s in p00},
             'interaction': {s: p11[s] - p10[s] - p01[s] + p00[s] for s in p00}}
    residual = max(abs(sum(t[s] for t in terms.values()) - (p11[s] - p00[s])) for s in p00)
    norms = {name: .5 * sum(abs(v) for v in term.values()) for name, term in terms.items()}
    return {'signed_terms': terms, 'half_l1_norms': norms, 'total_tv': tv(p11, p00),
            'sum_component_norms': sum(norms.values()), 'reconstruction_residual': residual,
            'interpretation': 'Signed factorial contrast for specified switches; component norms are not additive error shares.'}


def main(out):
    start = time.monotonic()
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=HERE, text=True):
        raise RuntimeError('E2-01 requires a clean committed tree')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=HERE, text=True).strip()
    sources = [HERE / 'experiment.py', HERE / 'DERIVATION.md', HERE / 'PROTOCOL_E2_01.md',
               MATCH / 'bosonized.py', MATCH / 'continuum.py', MATCH / 'protocol_m1_01.json',
               E1 / 'raw/effective_checkpoint.npz', E1 / 'raw/coupling_2.json']
    save(out / 'provenance.json', {'protocol': 'E2-01', 'source_commit': commit,
         'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
         'source_sha256': {str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}})
    check = verify_small()
    save(out / 'verification.json', check)
    if not check['pass']:
        raise RuntimeError('Small independent matrix verification failed')
    card = json.loads((MATCH / 'protocol_m1_01.json').read_text())
    states, omega, h = bosonized.build_hamiltonian(card, 12, 128, .5)
    w_parent, _, _, _ = bosonized.compressed_weyl(card, states, omega)
    checkpoint = np.load(E1 / 'raw/effective_checkpoint.npz')
    j = checkpoint['embedding']
    retained = j.shape[1]
    if (len(h), retained) != (1280, 72):
        raise RuntimeError('Frozen dimensions changed')
    complement = qr(j, mode='full')[0][:, retained:]
    a = j.conj().T @ h @ j
    c = complement.conj().T @ h @ complement
    ea, va = eigh(a)
    ec, vc = eigh(c, driver='evd')
    coordinates = np.column_stack((j @ va, complement @ vc))
    hc = coordinates.conj().T @ h @ coordinates
    b = hc[retained:, :retained]
    diagonal = np.concatenate((ea, ec))
    h0 = np.diag(diagonal)
    exchange = np.zeros_like(hc)
    exchange[retained:, :retained] = b
    exchange[:retained, retained:] = b.conj().T
    w = coordinates.conj().T @ w_parent @ coordinates
    w0 = w.copy()
    w0[:retained, retained:] = 0.
    w0[retained:, :retained] = 0.
    checks = {'parent_dimension': len(h), 'retained_dimension': retained,
              'memory_dimension': len(h) - retained,
              'coordinate_isometry_error': fro(coordinates.conj().T @ coordinates - np.eye(len(h))),
              'block_reconstruction_error': fro(hc - h0 - exchange),
              'archived_A_error': fro(a - checkpoint['effective_hamiltonian']),
              'parent_hermiticity_error': fro(h - h.conj().T),
              'parent_instrument': instrument_checks(w), 'pinched_instrument': instrument_checks(w0)}
    wr = w[:retained, :retained]
    parent_deficit_pp = (np.eye(retained) - (w.conj().T @ w)[:retained, :retained]) / 2.
    compressed_deficit = (np.eye(retained) - wr.conj().T @ wr) / 2.
    extra = w[retained:, :retained].conj().T @ w[retained:, :retained] / 2.
    checks['compressed_overflow_identity_error'] = fro(compressed_deficit - parent_deficit_pp - extra)
    # The full-parent eigensystem is used only for references and diagnostic bounds.
    ep, vp = eigh(h, driver='evd')
    reference, reference_loss = bosonized.histories(ep, vp, w_parent, OLD + NEW, 0.)
    phase = np.vdot(vp[:, 0], coordinates[:, 0])
    phase = phase.conjugate() / abs(phase)
    initial = np.zeros(len(h), dtype=complex)
    initial[0] = phase
    preparation_error = fro(vp[:, 0] - coordinates @ initial)
    checks['preparation_vector_error'] = preparation_error
    checks['prepared_parent_residual'] = fro(h @ (coordinates @ initial) - ea[0] * (coordinates @ initial))
    residual_keys = ['coordinate_isometry_error', 'block_reconstruction_error', 'archived_A_error',
                     'parent_hermiticity_error', 'compressed_overflow_identity_error',
                     'preparation_vector_error', 'prepared_parent_residual']
    checks['pass'] = bool(all(checks[k] < 1e-9 for k in residual_keys)
                          and all(checks[k]['completion_residual'] < 1e-9 for k in ['parent_instrument', 'pinched_instrument']))
    save(out / 'construction.json', checks)
    if not checks['pass']:
        raise RuntimeError('Construction prerequisite failed')
    parent_propagate = evolution(ep, coordinates.conj().T @ vp, ep[0])
    shifted_diagonal = diagonal - ep[0]
    def uncoupled(x, dt):
        return np.exp(-1j * shifted_diagonal * dt)[:, None] * x
    commutators = commutator_upper(ea, ec, b)
    split = MemorySplit(ea, ec, b, energy_shift=ep[0])
    save(out / 'commutators.json', {**commutators, 'exchange_singular_values': split.s.tolist(),
                                   'all_singular_values_retained': True})
    print(json.dumps({'stage': 'construction', 'seconds': time.monotonic() - start, 'checks_pass': True}), flush=True)
    old_data = json.loads((E1 / 'raw/coupling_2.json').read_text())
    records = {}
    rankings = []
    for times in OLD + NEW:
        name = key(times)
        old = times in OLD
        rows = {}
        switches = ['00', '10', '01', '11'] if old else ['00', '11']
        for switch in switches:
            prop = parent_propagate if switch[0] == '1' else uncoupled
            detector = w if switch[1] == '1' else w0
            row = history(prop, detector, initial, times, retained, b=b, details=old and switch == '11')
            row['tv_to_parent'] = tv(row['probabilities'], reference[name])
            rows[switch] = row
        reset = history(parent_propagate, w, initial, times, retained, reset=True)
        reset['tv_to_parent'] = tv(reset['probabilities'], reference[name])
        bound = autonomous_bound(parent_propagate, shifted_diagonal, w, initial, times, retained, preparation_error)
        bound['measured_tv'] = rows['00']['tv_to_parent']
        bound['dominates'] = bool(bound['measured_tv'] <= bound['bound'] + 1e-10)
        bound['accuracy_0.001'] = bool(bound['bound'] <= .001)
        memory = {}
        for steps in LEVELS:
            row = history(lambda x, dt: split.propagate(x, dt, steps), w, initial, times, retained)
            row['tv_to_parent'] = tv(row['probabilities'], reference[name])
            intervals = np.diff([0.] + times)
            row['bound'] = min(1., preparation_error + commutators['coefficient'] * float(np.sum(intervals**3)) / steps**2)
            row['bound_dominates'] = bool(row['tv_to_parent'] <= row['bound'] + 1e-10)
            memory[str(steps)] = row
        replay = None
        attribution = None
        if old:
            replay = {'parent_tv': tv(reference[name], old_data['parent_histories'][name]),
                      'autonomous_tv': tv(rows['00']['probabilities'], old_data['schemes']['joint']['histories'][name]['probabilities'])}
            if max(replay.values()) >= 1e-10:
                raise RuntimeError('Frozen E1 replay failed')
            attribution = factor_attribution(rows)
            if attribution['reconstruction_residual'] >= 1e-12:
                raise RuntimeError('Factorial identity failed')
            # Rank genuine outgoing parent detector sources and autonomous sources separately.
            wr_autonomous = w[:retained, :retained]
            x = initial[:retained, None].copy()
            labels, previous = [''], 0.
            for moment in times:
                x *= np.exp(-1j * shifted_diagonal[:retained] * (moment - previous))[:, None]
                for n, label in enumerate(labels):
                    weight = float(np.vdot(x[:, n], x[:, n]).real)
                    leak = .5 * fro(w[retained:, :retained] @ x[:, n])
                    residual = fro(b @ x[:, n])
                    rankings.append({'history': name, 'time': moment, 'prefix': label, 'weight': weight,
                                     'detector_leak_norm': leak, 'dynamical_residual_norm': residual,
                                     'conditional_detector_leak': leak / np.sqrt(weight) if weight > 0 else 0.,
                                     'conditional_residual': residual / np.sqrt(weight) if weight > 0 else 0.})
                shifted = wr_autonomous @ x
                x = np.concatenate(((x + shifted) / 2., (x - shifted) / 2.), axis=1)
                labels = [s + '+' for s in labels] + [s + '-' for s in labels]
                previous = moment
        record = {'times': times, 'old_calibration': old, 'parent': reference[name],
                  'parent_overflow_by_step': reference_loss[name], 'switches': rows,
                  'factorial': attribution, 'absorbing_reset': reset,
                  'autonomous_diagnostic_bound': bound, 'memory': memory, 'E1_replay': replay}
        save(out / ('history_' + name.replace(',', '_') + '.json'), record)
        records[name] = record
        print(json.dumps({'stage': 'history', 'times': times, 'autonomous_tv': rows['00']['tv_to_parent'],
                          'memory_tv': {n: r['tv_to_parent'] for n, r in memory.items()},
                          'bound_512': memory['512']['bound'], 'seconds': time.monotonic() - start}), flush=True)
    save(out / 'rankings.json', {'by_detector_leak': sorted(rankings, key=lambda r: -r['detector_leak_norm']),
                               'by_dynamical_residual': sorted(rankings, key=lambda r: -r['dynamical_residual_norm']),
                               'scope': 'Unnormalized autonomous prefixes; no low-weight branches omitted.'})
    summaries = {}
    for label, protocols in [('old', OLD), ('new', NEW)]:
        group = [records[key(t)] for t in protocols]
        auto = max(r['switches']['00']['tv_to_parent'] for r in group)
        measured = max(r['memory']['512']['tv_to_parent'] for r in group)
        upper = max(r['memory']['512']['bound'] for r in group)
        summaries[label] = {'autonomous_max_tv': auto, 'memory_max_tv': measured, 'memory_max_bound': upper,
                            'autonomous_diagnostic_max_bound': max(r['autonomous_diagnostic_bound']['bound'] for r in group),
                            'reset_max_tv': max(r['absorbing_reset']['tv_to_parent'] for r in group),
                            'gates': {'memory_measured_accuracy': measured <= .001,
                                      'memory_method_improvement': measured <= .5 * auto,
                                      'memory_analytic_bound_evaluated_in_floating_point': upper <= .001,
                                      'all_memory_bounds_dominate': all(q['bound_dominates'] for r in group for q in r['memory'].values()),
                                      'autonomous_diagnostic_bound_dominates': all(r['autonomous_diagnostic_bound']['dominates'] for r in group),
                                      'autonomous_diagnostic_bound_accuracy': all(r['autonomous_diagnostic_bound']['accuracy_0.001'] for r in group)}}
    summary = {'protocol': 'E2-01', 'groups': summaries, 'construction_pass': checks['pass'],
               'dimensions': {'retained': retained, 'memory': len(h)-retained, 'total': len(h)},
               'state_space_compression': False, 'seconds': time.monotonic() - start,
               'scope': 'Finite E=12 parent only. E1 and continuum verdicts unchanged. No interval-arithmetic certificate.'}
    save(out / 'summary.json', summary)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    try:
        main(output)
    except Exception as error:
        save(output / 'FAILURE.json', {'exception': repr(error), 'traceback': traceback.format_exc(),
                                     'scope': 'Execution failure; do not convert into a physics verdict.'})
        raise
