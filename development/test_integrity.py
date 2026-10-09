"""Hostile checks for false greens, case masking, metadata, and report freshness."""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from checks import choose_axes, command_status, run_check
from integrity import (adjudication_problems, baseline_deltas, classifier_agrees,
    json_read, locked_manifest, outcome_returncode_agrees, profiles_env,
    profile_baseline, read_junit, source_identity, verify_assets)


def state_fixture():
    return {"declarations": {}, "open_passes": {},
        "bank_inventory": {"report": {"overall": "CLEAN", "flags": []}},
        "protected_inputs": {"digest": "fixture"}, "tool_versions": {"python": "fixture"}}


def manifest_fixture(cases, state):
    return {"pytest_cases": dict(cases), **deepcopy(state)}


class HostileIntegrityTests(unittest.TestCase):
    def test_stale_junit_is_deleted_before_child_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report = root / "old.xml"
            report.write_text('<testsuites><testsuite tests="1"/></testsuites>')
            run_check("fixture", [sys.executable, "-c", "print('no report emitted')"],
                      root, root, os.environ, fresh_reports=[report])
            self.assertFalse(report.exists())
            with self.assertRaises(FileNotFoundError):
                read_junit(report)

    def test_missing_empty_or_nonlocal_asset_blocks_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "empty").write_text("")
            for name in ("absent", "empty", "../outside"):
                with self.assertRaises(ValueError):
                    verify_assets(root, [name])

    def test_conflicting_validator_lines_cannot_produce_green(self):
        self.assertEqual(command_status("validate", 0, "PASS: ok\nFAIL: bad\n"), "FAIL")
        self.assertEqual(command_status("validate", 0, "PASS: a\nPASS: b\n"), "FAIL")
        self.assertEqual(command_status("validate_scoped", 0, "PASS: ok\n"), "FAIL")

    def test_pytest_collection_or_process_error_cannot_look_like_known_red(self):
        cases = {"test_fixture.py::test_a": "FAIL"}
        for rc in (0, 2, 3, 4, 5, 124, 127):
            self.assertFalse(outcome_returncode_agrees(cases, rc))
        self.assertTrue(outcome_returncode_agrees(cases, 1))

    def test_duplicate_json_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            path.write_text('{"status":"FAIL","status":"PASS"}')
            with self.assertRaises(ValueError):
                json_read(path)

    def test_manifest_cannot_change_without_hash_lock_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            payload = json.dumps({"authorization": "OBSERVATIONS_ONLY_NO_ADJUDICATION"})
            path.write_text(payload)
            path.with_suffix(".json.sha256").write_text(hashlib.sha256(payload.encode()).hexdigest())
            self.assertEqual(locked_manifest(path)["authorization"], "OBSERVATIONS_ONLY_NO_ADJUDICATION")
            path.write_text(payload + " ")
            with self.assertRaises(ValueError):
                locked_manifest(path)

    def test_empty_or_duplicate_or_conflicting_junit_is_rejected(self):
        case = '<testcase classname="test_fixture" name="test_a"/>'
        bad = [('<testsuite tests="0" failures="0" errors="0" skipped="0"/>'),
               ('<testsuite tests="2" failures="0" errors="0" skipped="0">' + case * 2 + '</testsuite>'),
               ('<testsuite tests="1" failures="1" errors="0" skipped="1"><testcase classname="test_fixture" name="test_a"><failure/><skipped/></testcase></testsuite>'),
               ('<testsuite tests="99" failures="0" errors="0" skipped="0">' + case + '</testsuite>')]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "suite.xml"
            for xml in bad:
                path.write_text(xml)
                with self.assertRaises(ValueError):
                    read_junit(path)

    def test_new_declared_set_member_is_detected_even_when_node_stays_red(self):
        state = state_fixture()
        node = "test_fixture.py::test_a"
        state["declarations"][node] = {"cases": {"old": "P1"}, "live_cases": ["old"]}
        state["open_passes"]["P1"] = {"status": "OPEN", "symptomless": False}
        cases = {node: "FAIL"}
        manifest = manifest_fixture(cases, state)
        state["declarations"][node]["live_cases"].append("new")
        self.assertTrue(baseline_deltas(manifest, cases, state))
        self.assertIn("NEW_DECLARED_SET_MEMBER", {p["kind"] for p in adjudication_problems(cases, state)})

    def test_removed_declared_case_is_detected_as_stale(self):
        state = state_fixture()
        node = "test_fixture.py::test_a"
        state["declarations"][node] = {"cases": {"old": "P1"}, "live_cases": []}
        state["open_passes"]["P1"] = {"status": "OPEN", "symptomless": False}
        self.assertIn("STALE_DECLARED_CASE", {p["kind"] for p in adjudication_problems({node: "FAIL"}, state)})

    def test_closed_unknown_or_orphan_pass_cannot_authorize_red(self):
        state = state_fixture()
        node = "test_fixture.py::test_a"
        state["declarations"][node] = {"cases": {"old": "P1"}, "live_cases": ["old"]}
        for passes in ({}, {"P1": {"status": "CLOSED", "symptomless": False}}):
            state["open_passes"] = passes
            self.assertIn("NON_OPEN_OR_UNKNOWN_ADJUDICATION", {p["kind"] for p in adjudication_problems({node: "FAIL"}, state)})
        state["open_passes"] = {"P2": {"status": "OPEN", "symptomless": False}}
        self.assertIn("ORPHANED_OPEN_PASS", {p["kind"] for p in adjudication_problems({node: "FAIL"}, state)})

    def test_bank_inventory_delta_is_not_silenced_by_exit_zero(self):
        state = state_fixture()
        manifest = manifest_fixture({"test_fixture.py::test_a": "PASS"}, state)
        state["bank_inventory"]["report"]["flags"].append("new")
        self.assertIn("BANK_INVENTORY", {d["kind"] for d in baseline_deltas(manifest, manifest["pytest_cases"], state)})

    def test_classifier_must_match_raw_failures_and_case_audit(self):
        cases = {"test_fixture.py::test_a": "FAIL"}
        problems = [{"kind": "UNDECLARED_FAILING_TEST", "node": "test_fixture.py::test_a"}]
        output = "  *** NEW RED: test_fixture.py::test_a\nFAIL: the failing set is not the declared set.\n"
        self.assertTrue(classifier_agrees(cases, problems, output, 1))
        self.assertFalse(classifier_agrees(cases, problems, output, 0))
        self.assertFalse(classifier_agrees(cases, problems, "No new red.\n", 0))
        problems.append({"kind": "NEW_DECLARED_SET_MEMBER", "node": "test_fixture.py::test_a", "case": "hidden"})
        self.assertFalse(classifier_agrees(cases, problems, output, 1))

    def test_matching_observed_red_is_not_automatically_authorized(self):
        state = state_fixture()
        cases = {"test_fixture.py::test_a": "FAIL"}
        problems = adjudication_problems(cases, state)
        axes = choose_axes([], [], problems, True, state, cases, "default")
        self.assertEqual(axes["engineering_integrity"], "BLOCKED")
        self.assertEqual(axes["scientific_status"], "REVIEW_REQUIRED")

    def test_engineering_green_does_not_make_open_science_green(self):
        state = state_fixture()
        state["open_passes"] = {"P": {"status": "OPEN", "symptomless": True}}
        axes = choose_axes([], [], [], True, state, {"test_fixture.py::test_a": "PASS"}, "default")
        self.assertEqual(axes["engineering_integrity"], "PASS")
        self.assertEqual(axes["scientific_status"], "REVIEW_REQUIRED")
        self.assertFalse(axes["scientific_approval"])

    def test_profiles_are_separate_and_explicit(self):
        original = {"GRUT_FULL_MUTATION": "1", "GRUT_RUN_SLOW": "1"}
        self.assertNotIn("GRUT_FULL_MUTATION", profiles_env(original, "default"))
        self.assertNotIn("GRUT_RUN_SLOW", profiles_env(original, "full-mutation"))
        self.assertNotIn("GRUT_FULL_MUTATION", profiles_env(original, "slow-falsifiers"))
        self.assertEqual(profiles_env(original, "slow-falsifiers")["GRUT_RUN_SLOW"], "1")
        node = "test_layer5_overturning.py::TestLayer5::test_every_cited_falsifier_exits_zero"
        manifest = {"pytest_cases": {node: "SKIP", "test_mutation_battery.py::TestBattery::test_a": "PASS"}}
        self.assertEqual(profile_baseline(manifest, "slow-falsifiers")[node], "PASS")
        self.assertEqual(len(profile_baseline(manifest, "full-mutation")), 1)

    def test_source_sha_and_dirty_state_follow_actual_git_checkout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / "input.txt").write_text("original")
            subprocess.run(["git", "add", "input.txt"], cwd=root, check=True)
            subprocess.run(["git", "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "fixture"], cwd=root, check=True)
            first = source_identity(root)
            sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
            self.assertEqual(first["source_commit"], sha)
            self.assertFalse(first["working_tree_dirty"])
            (root / "input.txt").write_text("changed")
            second = source_identity(root)
            self.assertTrue(second["working_tree_dirty"])
            self.assertEqual(second["source_commit"], sha)
            self.assertNotEqual(first["protected_inputs"], second["protected_inputs"])

    def test_engineering_code_changes_are_detected_even_if_head_and_dirty_bit_stay_same(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / "input.txt").write_text("original")
            subprocess.run(["git", "add", "input.txt"], cwd=root, check=True)
            subprocess.run(["git", "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "fixture"], cwd=root, check=True)
            (root / "development").mkdir()
            code = root / "development/fixture.py"
            code.write_text("before")
            first = source_identity(root)
            code.write_text("after")
            second = source_identity(root)
            self.assertEqual(first["source_commit"], second["source_commit"])
            self.assertTrue(first["working_tree_dirty"] and second["working_tree_dirty"])
            self.assertEqual(first["protected_inputs"], second["protected_inputs"])
            self.assertNotEqual(first["engineering_inputs"], second["engineering_inputs"])

    def test_actual_pytest_new_red_and_removed_declared_red_are_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            test = root / "test_fixture.py"
            report = root / "suite.xml"
            def execute(text):
                test.write_text(text)
                result = subprocess.run([sys.executable, "-m", "pytest", "-q", "--junitxml=" + str(report)],
                    cwd=root, env=profiles_env(os.environ, "default"), capture_output=True)
                cases, _ = read_junit(report)
                self.assertTrue(outcome_returncode_agrees(cases, result.returncode))
                return cases
            good = execute('def test_a():\n    assert True\n\ndef test_b():\n    assert True\n')
            state = state_fixture()
            manifest = manifest_fixture(good, state)
            bad = execute('def test_a():\n    assert True\n\ndef test_b():\n    assert False\n')
            self.assertEqual(len(baseline_deltas(manifest, bad, state)), 1)
            self.assertEqual(adjudication_problems(bad, state)[0]["kind"], "UNDECLARED_FAILING_TEST")
            state["declarations"]["test_fixture.py::test_b"] = {"cases": {"leak": "P"}, "live_cases": ["leak"]}
            state["open_passes"]["P"] = {"status": "OPEN", "symptomless": False}
            good_again = execute('def test_a():\n    assert True\n\ndef test_b():\n    assert True\n')
            self.assertEqual(adjudication_problems(good_again, state)[0]["kind"], "STALE_DECLARATION")


if __name__ == "__main__":
    unittest.main()
