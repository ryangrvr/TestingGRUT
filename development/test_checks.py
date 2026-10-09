"""Regression checks for status reporting; no scientific claims are exercised here."""
from pathlib import Path
import tempfile
import unittest

from checks import command_status, junit_counts


class ReportingTests(unittest.TestCase):
    def test_exit_zero_cannot_hide_bank_flags(self):
        self.assertEqual(command_status("bank_gate", 0,
                         "OVERALL: FLAG-FOR-FIREWALL\n"), "REVIEW_REQUIRED")
        self.assertEqual(command_status("bank_gate", 0, "OVERALL: CLEAN\n"), "PASS")

    def test_missing_or_conflicting_gate_verdict_is_not_clean(self):
        for output in ("", "CLEAN\n", "OVERALL: CLEAN\nOVERALL: BLOCK\n"):
            self.assertEqual(command_status("bank_gate", 0, output), "REVIEW_REQUIRED")

    def test_failed_command_cannot_be_promoted_by_clean_output(self):
        self.assertEqual(command_status("bank_gate", 1, "OVERALL: CLEAN\n"), "FAIL")

    def test_junit_preserves_failures_errors_and_skips(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "suite.xml"
            self.assertIsNone(junit_counts(path))
            path.write_text('<testsuites><testsuite tests="5" failures="1" errors="1" '
                            'skipped="2"/></testsuites>', encoding="utf-8")
            self.assertEqual(junit_counts(path),
                             {"tests": 5, "failures": 1, "errors": 1, "skipped": 2})


if __name__ == "__main__":
    unittest.main()
