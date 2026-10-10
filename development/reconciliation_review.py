"""Issue #2: test a concrete repair in a disposable checkout; never rebaseline.

The proposed patch touches protected files ONLY in that checkout. Applying it to
the tracked record, changing public wave stamps, and disposing of OPEN passes
require a reviewed transition. This program grants none of those authorizations.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from integrity import (adjudication_problems, classifier_agrees, locked_manifest,
                       profiles_env, read_junit, source_identity)

TARGETS = (
    "test_expected_red.py::TestExpectedRed::test_a_declaration_citing_a_closed_pass_is_refused",
    "test_expected_red.py::TestExpectedRed::test_a_moot_pass_cannot_be_cited_by_a_declaration",
    "test_expected_red.py::TestExpectedRed::test_a_symptomless_pass_is_allowed_and_printed",
    "test_expected_red.py::TestExpectedRed::test_the_declared_state_is_accepted",
    "test_public_doc.py::TestPublicDoc::test_the_calc_manifest_matches_the_repository",
    "test_public_numbers.py::TestPublicNumbers::test_no_drift_from_the_register",
    "test_public_numbers.py::TestPublicNumbers::test_the_stamped_test_count_is_true",
)
PATCH_PATHS = (
    "provenance/test_expected_red.py", "provenance/emit_public_numbers.py",
    "provenance/CALC_MANIFEST.txt", "PUBLIC_NUMBERS.md", "docs/WHERE_IT_STOPS.md",
)


def replace_once(body, old, new):
    if body.count(old) != 1:
        raise ValueError("Repair precondition is not unique; review source drift")
    return body.replace(old, new, 1)


def isolated_fixture_source(body):
    body = replace_once(body, "import contextlib\n", "from copy import deepcopy\nimport contextlib\n")
    body = replace_once(body, "        self._passes = X.open_passes\n", """        self._passes = X.open_passes
        self._declared = X.DECLARED
        # Unit inputs must be valid independently of the live integration state.
        X.DECLARED = deepcopy(self._declared)
        for entry in X.DECLARED.values():
            cases = frozenset(entry['cases'])
            entry['enumerate'] = lambda cases=cases: set(cases)
        cited = {pid for entry in X.DECLARED.values() for pid in entry['cases'].values()}
        fixture_passes = {pid: {'status': 'OPEN', 'symptomless': False} for pid in cited}
        fixture_passes['P1B-SHOWN-ON-LEDGER-INPUTS'] = {'status': 'OPEN', 'symptomless': True}
        X.open_passes = lambda: deepcopy(fixture_passes)
""")
    body = replace_once(body, "        X.open_passes = self._passes\n",
                        "        X.open_passes = self._passes\n        X.DECLARED = self._declared\n")
    # This pass is actually cited by the synthetic declarations. P1A is not.
    if body.count('"P1A-EDGE-REPRESENTATION"') != 2:
        raise ValueError("Unexpected closed/moot mutation fixture layout")
    return body.replace('"P1A-EDGE-REPRESENTATION"', '"P4-TERMINATION-KAPPA-RESULT"')


def run(argv, cwd, output, name, env, expected=0, timeout=900):
    result = subprocess.run(argv, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)
    (output / (name + ".log")).write_text(result.stdout + result.stderr)
    if result.returncode != expected:
        raise RuntimeError(f"{name}: exit {result.returncode}; expected {expected}")
    return result.stdout


def review(root, output, stamp, reconciliation=None):
    root, output = root.resolve(), output.resolve()
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", stamp):
        raise ValueError("Explicit fixed review stamp required")
    manifest = locked_manifest(root / "development/expected_red_manifest.json")
    if reconciliation:
        from owner_transition import reconciled_manifest
        manifest = reconciled_manifest(manifest, reconciliation, root)
    original = source_identity(root)
    if original['protected_inputs'] != manifest['protected_inputs'] or original['untracked_protected_paths']:
        raise ValueError("Protected source differs from frozen review input")
    output.mkdir(parents=True, exist_ok=False)
    # Retain OS/runtime settings; neutralize pytest profile switches explicitly.
    import os
    env = profiles_env(os.environ, "default")
    # Same-length stamp edits within one timestamp tick can otherwise reuse the
    # pre-repair bytecode created by the before-run, concealing the source edit.
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    preconditions = {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in PATCH_PATHS}
    with tempfile.TemporaryDirectory(prefix="grut-reviewed-overlay-") as tmp:
        overlay = Path(tmp) / "checkout"
        run(["git", "clone", "--quiet", "--shared", "--no-hardlinks", str(root), str(overlay)],
            root, output, "clone", env)
        if reconciliation:
            # A local pre-push review can contain the authorized pass/declaration
            # edits not yet in HEAD. Copy only these hash-checked two paths.
            from owner_transition import PATHS as OWNER_PATHS
            for p in OWNER_PATHS:
                (overlay / p).write_bytes((root / p).read_bytes())
        for p, digest in preconditions.items():
            if hashlib.sha256((overlay / p).read_bytes()).hexdigest() != digest:
                raise ValueError("Clone differs from frozen protected input")
        before = output / "target-before.xml"
        run([sys.executable, "-m", "pytest", "-q", *TARGETS, "--junitxml=" + str(before)],
            overlay / "provenance", output, "target-before", env, expected=1)
        old, _ = read_junit(before)
        if old != {n: 'FAIL' for n in TARGETS}:
            raise ValueError("The seven claimed defects were not independently reproduced")
        path = overlay / PATCH_PATHS[0]
        path.write_text(isolated_fixture_source(path.read_text()))
        count = len(manifest['pytest_cases'])
        path = overlay / PATCH_PATHS[1]
        body = replace_once(path.read_text(), "STAMPED_TEST_COUNT = 240", f"STAMPED_TEST_COUNT = {count}")
        body = replace_once(body, 'STAMP_DATE = "2026-08-12"', f'STAMP_DATE = "{stamp}"')
        body = replace_once(body, 'WAVE_DATE = "2026-08-18"', f'WAVE_DATE = "{stamp}"')
        path.write_text(body)
        for script in ("build_calc_manifest.py", "emit_public_numbers.py", "build_public_doc.py"):
            run([sys.executable, "provenance/" + script, "--write"], overlay, output, script, env)
        patch = run(["git", "diff", "--", *PATCH_PATHS], overlay, output, "proposal-diff", env)
        (output / "review-only.patch").write_text(patch)
        changed = run(["git", "diff", "--name-only"], overlay, output, "changed-paths", env).splitlines()
        extra = set(changed) - set(PATCH_PATHS)
        allowed_extra = OWNER_PATHS if reconciliation else set()
        if not set(PATCH_PATHS) <= set(changed) or not extra <= allowed_extra:
            raise ValueError("Repair touched an unexpected protected input")
        full = output / "overlay.xml"
        run([sys.executable, "-m", "pytest", "-q", "--junitxml=" + str(full)],
            overlay / "provenance", output, "overlay-suite", env, expected=1)
        cases, counts = read_junit(full)
        expected = {**manifest['pytest_cases'], **{n: 'PASS' for n in TARGETS}}
        if cases != expected:
            raise ValueError("Overlay introduced an unapproved outcome change or missing test")
        # This runs the live classifier and real enumerators, never the synthetic
        # unit fixture. Owner scoping changes only the recorded obligations.
        classifier = run([sys.executable, "provenance/expected_red.py"], overlay,
                         output, "overlay-live-classifier", env, expected=1)
        new_cases = re.findall(r"\*\*\* NEW RED \(case\): [^\n]+\n[ \t]+([^\n]+)", classifier)
        orphans = re.findall(r"\*\*\* ORPHANED OPEN PASS '([^']+)'", classifier)
        expected_new = 0 if reconciliation else 4
        expected_orphans = ['P6-STALE-NETS-IN-STANDING-DOCS'] if reconciliation else ['P1A-EDGE-REPRESENTATION', 'P6-STALE-NETS-IN-STANDING-DOCS']
        if len(new_cases) != expected_new or sorted(orphans) != expected_orphans:
            raise ValueError("Live guard obligations changed in the fixture repair")
        if re.search(r"\*\*\* NEW RED: ", classifier):
            raise ValueError("An undeclared failing test remains in the proposed overlay")
        captured = json.loads(run([sys.executable, 'development/capture_state.py', '--root', str(overlay)],
                                  overlay, output, 'overlay-state', env))
        for key in ('declarations', 'open_passes', 'bank_inventory'):
            if captured[key] != manifest[key]:
                raise ValueError('Overlay altered owner/scientific observations: ' + key)
        audit = adjudication_problems(cases, captured)
        if not classifier_agrees(cases, audit, classifier, 1):
            raise ValueError('Overlay classifier disagrees with independent audit/raw results')
        # Exact protected declaration/pass/seal byte invariance, including unused records.
        invariant_paths = ['provenance/expected_red.py', 'provenance/OPEN_PASSES.txt',
                           *subprocess.check_output(['git', 'ls-files', 'provenance/prereg/'], cwd=root, text=True).splitlines()]
        for p in invariant_paths:
            if (root / p).read_bytes() != (overlay / p).read_bytes():
                raise ValueError(f"Guard/pass/seal mutation in proposal: {p}")
    if source_identity(root)['protected_inputs'] != original['protected_inputs']:
        raise ValueError("Original protected inputs changed during review")
    result = {'scope': 'REVIEW_ONLY_NOT_APPLIED', 'source_commit': original['source_commit'],
              'owner_reconciliation': str(reconciliation) if reconciliation else None,
              'protected_inputs': original['protected_inputs'], 'precondition_sha256': preconditions,
              'proposed_wave_stamp': stamp, 'historical_correction_stamp_unchanged': True,
              'before_target_cases': old, 'overlay_counts': dict(Counter(cases.values())),
              'pytest_counts': counts, 'changed_outcomes': list(TARGETS), 'changed_paths': changed,
              'repair_paths': list(PATCH_PATHS),
              'owner_paths_copied_from_review_input': sorted(OWNER_PATHS) if reconciliation else [],
              'classifier_agreement_with_raw_and_independent_audit': True,
              'overlay_adjudication_problems': audit,
              'remaining_new_pointer_cases': sorted(new_cases), 'remaining_orphan_passes': sorted(orphans),
              'allowlist_passes_and_all_seals_byte_identical_to_review_input': True,
              'patch_sha256': hashlib.sha256(patch.encode()).hexdigest(),
              'engineering_clearance': False, 'scientific_approval': False,
              'external_review': 'NOT_PERFORMED_BY_THIS_RUN', 'card_1_authorized': False}
    (output / 'review.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stamp', required=True)
    parser.add_argument('--reconciliation', type=Path)
    args = parser.parse_args()
    review(args.root, args.output, args.stamp, args.reconciliation)
