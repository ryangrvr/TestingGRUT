"""EX1-01 scorer checks using deterministic SYNTHETIC fixtures, never observations."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

import experiment_verify as ev


def rounded_counts(probabilities, n):
    """Deterministic largest-remainder fixture, with no shot sampling."""
    raw = [n * p for p in probabilities]
    counts = [math.floor(x) for x in raw]
    for i in sorted(range(len(raw)), key=lambda i: (-(raw[i] - counts[i]), i))[:n - sum(counts)]:
        counts[i] += 1
    return {**{format(i, '04b'): count for i, count in enumerate(counts)}, 'erasure': 0}


def fixture():
    data = {'conditions': {}}
    for name, n in ev.COUNTS_REQUIRED.items():
        probabilities = ev.raw_distribution(ev.run('11', dephase=name in ('P11', 'I_P'),
                                                   immediate=name in ('I_C', 'I_P')))
        data['conditions'][name] = rounded_counts(probabilities, n)
    return data


def run(output):
    checks = {}
    baseline = fixture()
    with tempfile.TemporaryDirectory(prefix='EX1_SYNTHETIC_ONLY_') as directory:
        path = Path(directory) / 'SYNTHETIC_ONLY_counts.json'

        def score(data):
            path.write_text(json.dumps(data), encoding='utf-8')
            return ev.score_counts(path)

        def rejects(data):
            try:
                score(data)
            except ValueError:
                return True
            return False

        scored = score(baseline)
        checks['synthetic_ideal_counts_pass_statistical_gates'] = all(scored['gates'].values())
        checks['statistical_success_does_not_certify_apparatus'] = scored['overall_empirical_claim'].startswith('not_evaluated:')
        bad = copy.deepcopy(baseline)
        previous = bad['conditions']['C11']['0000']
        bad['conditions']['C11']['0000'] = -1
        bad['conditions']['C11']['1111'] += previous + 1
        checks['negative_count_rejected_at_unchanged_total'] = rejects(bad)
        bad = copy.deepcopy(baseline)
        del bad['conditions']['I_P']
        checks['missing_condition_rejected'] = rejects(bad)
        bad = copy.deepcopy(baseline)
        del bad['conditions']['C11']['erasure']
        checks['missing_erasure_category_rejected'] = rejects(bad)
        bad = copy.deepcopy(baseline)
        bad['conditions']['C11']['0000'] += 1
        checks['wrong_total_rejected'] = rejects(bad)
        bad = copy.deepcopy(baseline)
        bad['conditions']['C11']['0000'] = True
        checks['boolean_not_accepted_as_integer_count'] = rejects(bad)

        erased = copy.deepcopy(baseline)
        moved = 1000
        erased['conditions']['C11']['0000'] -= moved
        erased['conditions']['C11']['erasure'] += moved
        erasure_scored = score(erased)
        before = scored['records']['C11']['action_agreement']
        after = erasure_scored['records']['C11']['action_agreement']
        checks['erasure_uses_attempted_shot_denominator'] = abs(before - after - moved / 20000) < 1e-14
        checks['erasure_rate_preserved'] = erasure_scored['records']['C11']['erasure_rate'] == moved / 20000
        checks['erasure_does_not_change_registered_total'] = abs(sum(erasure_scored['records']['C11']['joint_S_parity_A_plus_erasure']) - 1) < 1e-14

        all_erased = copy.deepcopy(baseline)
        for name, n in ev.COUNTS_REQUIRED.items():
            all_erased['conditions'][name] = {key: 0 for key in baseline['conditions'][name]}
            all_erased['conditions'][name]['erasure'] = n
        all_erased_scored = score(all_erased)
        checks['all_erasure_has_zero_unconditional_agreement'] = all(record['action_agreement'] == 0 for record in all_erased_scored['records'].values())
        checks['all_erasure_cannot_pass_primary_claim'] = not all_erased_scored['gates']['coherence_contrast'] and all_erased_scored['status'] == 'statistical_gates_fail'

    result = {'status': 'passed' if all(checks.values()) else 'failed',
              'input_kind': 'deterministic_synthetic_fixtures_only', 'actual_experimental_data': False,
              'device_adapter_status': 'pending', 'checks': checks,
              'synthetic_reference_metrics': scored['metrics'],
              'synthetic_erasure_metrics': erasure_scored['metrics'],
              'verifier_sha256': hashlib.sha256(Path(ev.__file__).read_bytes()).hexdigest(),
              'precheck_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'status': result['status'], 'checks': len(checks),
                      'input_kind': result['input_kind'], 'actual_experimental_data': False}))
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    raise SystemExit(run(parser.parse_args().output))
