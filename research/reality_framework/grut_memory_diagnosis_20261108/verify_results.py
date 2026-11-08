"""Verify E2 stored arithmetic and provenance; not an independent physics solver."""
import argparse
import hashlib
import json
from pathlib import Path


def close(actual, expected, tolerance=1e-12):
    if abs(actual - expected) > tolerance:
        raise AssertionError(f'{actual} != {expected} within {tolerance}')


def distance(p, q):
    assert set(p) == set(q)
    return sum(abs(p[s] - q[s]) for s in p) / 2


def verify(raw):
    here = Path(__file__).resolve().parent
    provenance = json.loads((raw / 'provenance.json').read_text())
    frozen = json.loads((here / 'raw/provenance.json').read_text())
    assert frozen['source_commit'] == '3dce73fb6d9e6822281c048f8cc646ebdd8b74ba'
    assert provenance['source_sha256'] == frozen['source_sha256'], 'Frozen scientific sources changed'
    for path, expected in provenance['source_sha256'].items():
        assert hashlib.sha256((here.parent / path).read_bytes()).hexdigest() == expected, path
    c = json.loads((raw / 'commutators.json').read_text())
    coefficient = c['H0_double_commutator_Frobenius'] / 24 + c['V_double_commutator_Frobenius'] / 12
    close(coefficient, c['coefficient'])
    prep = json.loads((raw / 'construction.json').read_text())['preparation_vector_error']
    count, maximum_norm_error, records = 0, 0., []
    for path in sorted(raw.glob('history_*.json')):
        record = json.loads(path.read_text())
        records.append(record)
        parent = record['parent']
        rows = list(record['switches'].values()) + [record['absorbing_reset']] + list(record['memory'].values())
        for p in [parent] + [row['probabilities'] for row in rows]:
            count += 1
            assert min(p.values()) >= -1e-10
            error = abs(sum(p.values()) - 1)
            maximum_norm_error = max(maximum_norm_error, error)
            assert error < 1e-9
        for row in rows:
            close(distance(row['probabilities'], parent), row['tv_to_parent'])
            close(sum(row['ordinary_overflow_by_step']) + sum(row['reset_overflow_by_step']), row['probabilities']['overflow'])
        if record['factorial']:
            f = record['factorial']
            p = {s: r['probabilities'] for s, r in record['switches'].items()}
            for outcome in parent:
                close(f['signed_terms']['propagation'][outcome], p['10'][outcome] - p['00'][outcome])
                close(f['signed_terms']['instrument'][outcome], p['01'][outcome] - p['00'][outcome])
                close(f['signed_terms']['interaction'][outcome], p['11'][outcome] - p['10'][outcome] - p['01'][outcome] + p['00'][outcome])
                close(sum(v[outcome] for v in f['signed_terms'].values()), p['11'][outcome] - p['00'][outcome])
            for name, vector in f['signed_terms'].items():
                close(sum(abs(x) for x in vector.values()) / 2, f['half_l1_norms'][name])
            close(f['total_tv'], distance(p['00'], p['11']))
        b = record['autonomous_diagnostic_bound']
        close(b['propagation_budget'], sum(s['propagation_Frobenius_defect'] for s in b['stages']))
        close(b['instrument_budget'], sum(s['instrument_stacked_defect'] for s in b['stages']))
        delta = b['initial_vector_error'] + b['propagation_budget'] + b['instrument_budget']
        close(b['bound'], min(1, delta * (1 + b['final_ordinary_norm'])))
        assert record['switches']['00']['tv_to_parent'] <= b['bound'] + 1e-10
        previous, cubes = 0., 0.
        for t in record['times']:
            cubes += (t - previous)**3
            previous = t
        for n, row in record['memory'].items():
            close(row['bound'], min(1, prep + coefficient * cubes / int(n)**2))
            assert row['tv_to_parent'] <= row['bound'] + 1e-10
    assert len(records) == 4
    summary = json.loads((raw / 'summary.json').read_text())
    for group, old in [('old', True), ('new', False)]:
        relevant = [r for r in records if r['old_calibration'] == old]
        reported = summary['groups'][group]
        measured = max(r['memory']['512']['tv_to_parent'] for r in relevant)
        auto = max(r['switches']['00']['tv_to_parent'] for r in relevant)
        upper = max(r['memory']['512']['bound'] for r in relevant)
        close(measured, reported['memory_max_tv'])
        close(auto, reported['autonomous_max_tv'])
        close(upper, reported['memory_max_bound'])
        assert reported['gates']['memory_measured_accuracy'] == (measured <= .001)
        assert reported['gates']['memory_method_improvement'] == (measured <= auto / 2)
        assert reported['gates']['memory_analytic_bound_evaluated_in_floating_point'] == (upper <= .001)
        assert reported['gates']['autonomous_diagnostic_bound_accuracy'] == all(r['autonomous_diagnostic_bound']['bound'] <= .001 for r in relevant)
    return {'verified_distributions': count, 'maximum_normalization_error': maximum_norm_error,
            'run_source_commit': provenance['source_commit'],
            'run_at_exact_freeze_commit': provenance['source_commit'] == frozen['source_commit'],
            'source_hashes_match': True, 'factorial_identities_recomputed': True,
            'stored_bounds_and_TV_recomputed': True, 'frozen_gate_decisions_recomputed': True,
            'scope': 'Stored-output integrity checks, not independent dynamics or evaluator validation.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, default=Path(__file__).resolve().parent / 'raw')
    args = parser.parse_args()
    print(json.dumps(verify(args.raw), indent=2))
