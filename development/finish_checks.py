"""Verify the exact applied engineering repair; retain the historical snapshot.

This grants engineering integrity only. External review and scientific approval
remain separate. No observed outcome table can authorize a new failure or pass.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from checks import choose_axes, run_check, syntax_check
from integrity import (adjudication_problems, baseline_deltas, classifier_agrees,
                       json_read, locked_manifest, outcome_returncode_agrees,
                       profiles_env, read_junit, source_identity, verify_assets)
from owner_transition import reconciled_manifest

ROOT = Path(__file__).resolve().parents[1]
REVIEW_PIN = 'c32f672288fcad2a9d08772b8ef2bf46ddee70a2'
PATCH_SHA256 = '4c46efcb263d738d1909474cecb5905ecd1c94db5a6da11885ee08c0c2bf8ba1'
PATCH_PATHS = (
    'PUBLIC_NUMBERS.md', 'docs/WHERE_IT_STOPS.md', 'provenance/CALC_MANIFEST.txt',
    'provenance/emit_public_numbers.py', 'provenance/test_expected_red.py',
)
PIN_PATH = 'provenance/test_doc_register_pins.py'
PIN_PATCH_SHA256 = '2ed9a60973e223f4dc2d09d11d300939756e2c750737a870e4f5fc71dc348e90'
REPAIR_PATHS = PATCH_PATHS + (PIN_PATH,)
PASS_CHANGES = (
    'test_expected_red.py::TestExpectedRed::test_a_declaration_citing_a_closed_pass_is_refused',
    'test_expected_red.py::TestExpectedRed::test_a_moot_pass_cannot_be_cited_by_a_declaration',
    'test_public_doc.py::TestPublicDoc::test_the_calc_manifest_matches_the_repository',
    'test_public_numbers.py::TestPublicNumbers::test_no_drift_from_the_register',
    'test_public_numbers.py::TestPublicNumbers::test_the_stamped_test_count_is_true',
)


@contextmanager
def pinned_checkout(root):
    """Use a fixed local commit; never fetch a moving ref or alter the caller."""
    with tempfile.TemporaryDirectory(prefix='grut-pinned-review-') as tmp:
        snapshot = Path(tmp) / 'checkout'
        subprocess.run(['git', 'clone', '--quiet', '--shared', '--no-hardlinks',
                        str(root), str(snapshot)], check=True, capture_output=True)
        subprocess.run(['git', 'checkout', '--quiet', '--detach', REVIEW_PIN],
                       cwd=snapshot, check=True, capture_output=True)
        yield snapshot


def validate_repair_diff(changed_paths, patch, untracked, pin_patch):
    protected = {p for p in changed_paths if p != 'AGENT_COORDINATION.md'
                 and not p.startswith(('development/', '.github/workflows/'))}
    if protected != set(REPAIR_PATHS) or untracked:
        raise ValueError('Repair includes additional or missing protected changes')
    if hashlib.sha256(patch).hexdigest() != PATCH_SHA256:
        raise ValueError('Applied repair differs from the exact reviewed five-file patch')
    if hashlib.sha256(pin_patch).hexdigest() != PIN_PATCH_SHA256:
        raise ValueError('Pin fixture isolation differs from its exact frozen patch')


def repaired_expectation(base, protected_inputs):
    """Only five named FAIL→PASS changes; no declarations/bank/pass refresh."""
    expected = deepcopy(base)
    for node in PASS_CHANGES:
        if expected['pytest_cases'].get(node) != 'FAIL':
            raise ValueError('Repair change is not anchored to the pinned failed test: ' + node)
        expected['pytest_cases'][node] = 'PASS'
    expected['protected_inputs'] = protected_inputs
    return expected


def verify_applied_source(root, snapshot):
    identity = source_identity(root)
    subprocess.run(['git', 'merge-base', '--is-ancestor', REVIEW_PIN, 'HEAD'],
                   cwd=root, check=True, capture_output=True)
    patch_file = root / 'development/review_cycle_02/repair.patch'
    if hashlib.sha256(patch_file.read_bytes()).hexdigest() != PATCH_SHA256:
        raise ValueError('Frozen repair artifact changed')
    diff = subprocess.check_output(['git', 'diff', REVIEW_PIN, '--', *PATCH_PATHS], cwd=root)
    pin_file = root / 'development/finish/pin-fixture.patch'
    pin_diff = subprocess.check_output(['git', 'diff', REVIEW_PIN, '--', PIN_PATH], cwd=root)
    if hashlib.sha256(pin_file.read_bytes()).hexdigest() != PIN_PATCH_SHA256:
        raise ValueError('Frozen pin-fixture repair artifact changed')
    changed = subprocess.check_output(['git', 'diff', REVIEW_PIN, '--name-only'],
                                      cwd=root, text=True).splitlines()
    validate_repair_diff(changed, diff, identity['untracked_protected_paths'], pin_diff)
    base = locked_manifest(snapshot / 'development/expected_red_manifest.json')
    base = reconciled_manifest(base, snapshot / 'development/final_precard/P6_TRANSITION.json', snapshot)
    # Independently derive the resulting protected digest by applying the frozen
    # patch to the pinned tree. Never copy a current digest into the baseline.
    subprocess.run(['git', 'apply', '--check', str(patch_file)], cwd=snapshot,
                   check=True, capture_output=True)
    subprocess.run(['git', 'apply', str(patch_file)], cwd=snapshot,
                   check=True, capture_output=True)
    subprocess.run(['git', 'apply', '--check', str(pin_file)], cwd=snapshot,
                   check=True, capture_output=True)
    subprocess.run(['git', 'apply', str(pin_file)], cwd=snapshot,
                   check=True, capture_output=True)
    protected = source_identity(snapshot)['protected_inputs']
    if identity['protected_inputs'] != protected:
        raise ValueError('Current protected bytes differ from the reconstructed repair')
    subprocess.run(['git', 'reset', '--quiet', '--hard', REVIEW_PIN], cwd=snapshot,
                   check=True, capture_output=True)
    return repaired_expectation(base, protected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'development/results')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix='finished-', dir=args.output.resolve()))
    env = profiles_env(os.environ, 'default')
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    identity = source_identity(ROOT)
    checks, errors, deltas, problems, cases, counts, state = [], [], [], [], {}, None, {}
    agreement = False
    try:
        with pinned_checkout(ROOT) as snapshot:
            manifest = verify_applied_source(ROOT, snapshot)
            checks.append(run_check('pinned_review_integrity',
                [sys.executable, '-m', 'development.final_precard.verify_packet'], snapshot, out, env))
            # These historical owner/packet tests explicitly address the c32
            # input state. Execute every original method in that exact state.
            checks.append(run_check('historical_runner_tests',
                [sys.executable, '-m', 'unittest', 'discover', '-s', 'development', '-p', 'test_*.py'],
                snapshot, out, env))
        checks.append(run_check('applied_repair_tests',
            [sys.executable, '-m', 'unittest', 'discover', '-s', 'development', '-p', 'test_finish.py'],
            ROOT, out, env))
        checks.append(syntax_check(out))
        capture = [sys.executable, 'development/capture_state.py', '--root', str(ROOT)]
        checks.append(run_check('state_before', capture, ROOT, out, env))
        before = json_read(out / 'state_before.log')
        if before['source_commit'] != identity['source_commit']:
            raise ValueError('Captured source SHA differs from the checkout')
        for name, script in (('validate', 'validate.py'), ('validate_scoped', 'validate_scoped.py'),
                             ('bank_gate', 'bankgate.py')):
            checks.append(run_check(name, [sys.executable, 'provenance/' + script], ROOT, out, env))
        junit = out / 'provenance.xml'
        raw = run_check('provenance_suite',
            [sys.executable, '-m', 'pytest', '-q', '--junitxml=' + str(junit)],
            ROOT / 'provenance', out, env, timeout=900, fresh_reports=[junit])
        checks.append(raw)
        cases, counts = read_junit(junit)
        if not outcome_returncode_agrees(cases, raw['returncode']):
            raise ValueError('Raw pytest exit disagrees with the complete case inventory')
        classifier = run_check('expected_red_classifier',
            [sys.executable, 'provenance/expected_red.py'], ROOT, out, env)
        checks.append(classifier)
        checks.append(run_check('state_after', capture, ROOT, out, env))
        state = json_read(out / 'state_after.log')
        if state != before or source_identity(ROOT) != identity:
            raise ValueError('Source/state changed during execution')
        deltas = baseline_deltas(manifest, cases, state)
        problems = adjudication_problems(cases, state)
        agreement = classifier_agrees(cases, problems,
            (out / 'expected_red_classifier.log').read_text(), classifier['returncode'])
        # A bank CLI's zero exit is not scientific approval.
        if state['bank_inventory'] != manifest['bank_inventory']:
            raise ValueError('Bank inventory changed')
        verdicts = re.findall(r'^OVERALL: (\S+)\s*$',
                              (out / 'bank_gate.log').read_text(), re.MULTILINE)
        if verdicts != [state['bank_inventory']['report']['overall']]:
            raise ValueError('Structured bank inventory and CLI verdict disagree')
        assets = verify_assets(out, [c['log'] for c in checks] + ['provenance.xml'])
        (out / 'assets.json').write_text(json.dumps(assets, indent=2, sort_keys=True) + '\n')
    except Exception as exc:
        errors.append(f'{type(exc).__name__}: {exc}')
    critical = {'syntax', 'pinned_review_integrity', 'historical_runner_tests', 'applied_repair_tests',
                'state_before', 'state_after', 'validate', 'validate_scoped'}
    errors.extend(c['name'] + ' did not pass' for c in checks if c['name'] in critical and c['status'] != 'PASS')
    axes = choose_axes(errors, deltas, problems, agreement, state, cases, 'default')
    result = {'scope': 'APPLIED_DEVELOPMENT_REPAIR_ONLY', **identity,
              'review_pin': REVIEW_PIN, 'repair_sha256': PATCH_SHA256,
              'pin_fixture_repair_sha256': PIN_PATCH_SHA256,
              'permitted_outcome_changes': list(PASS_CHANGES), 'checks': checks,
              'provenance_tests': counts, 'case_outcomes': cases, 'axes': axes,
              'baseline_deltas': deltas, 'adjudication_problems': problems,
              'classifier_agreement_with_raw_and_independent_audit': agreement,
              'infrastructure_errors': errors, 'external_approval': False,
              'card_1_authorized': False, 'cards_consumed': 0}
    (out / 'summary.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    (out / 'summary.md').write_text(
        f"Engineering: **{axes['engineering_integrity']}**; science: **{axes['scientific_status']}**.\n\n"
        f"Applied exact five-file repair plus isolated pin fixture; source {identity['source_commit']}; dirty={identity['working_tree_dirty']}.\n\n"
        'External evaluator/CR-5 review pending; Card 1 CLOSED; no bank acceptance.\n')
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
            stream.write(f"engineering_integrity={axes['engineering_integrity']}\nscientific_status={axes['scientific_status']}\n")
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as stream:
            stream.write((out / 'summary.md').read_text())
    print(json.dumps({k: result[k] for k in ('source_commit', 'working_tree_dirty', 'axes',
        'provenance_tests', 'baseline_deltas', 'adjudication_problems', 'infrastructure_errors')}, sort_keys=True))
    print('REPORT=' + str(out / 'summary.json'))
    return 0 if axes['engineering_integrity'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
