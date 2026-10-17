"""Control-plane tests; no scientific campaign execution."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import runner


class DiscoveryTests(unittest.TestCase):
    def test_manifest_detects_mutation(self):
        frozen = runner.tracked_manifest()
        altered = dict(frozen)
        altered[next(iter(altered))] = "0" * 64
        with patch.object(Path, "read_text", return_value=json.dumps({"sha256": altered})):
            with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                runner.verify()

    def test_unsafe_id_rejected_before_tools(self):
        for value in ("../escape", "/tmp/escape", "a/b", "", "x" * 65):
            with self.assertRaises(ValueError):
                runner.run(value)

    def test_existing_run_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            here = Path(temp)
            (here / "runs" / "same").mkdir(parents=True)
            (here / "PROTOCOL.json").write_text(json.dumps({"budgets": {}}))
            with patch.object(runner, "HERE", here), patch.object(runner, "verify", return_value={}):
                with self.assertRaises(FileExistsError):
                    runner.run("same")

    def test_ledger_has_no_fake_review(self):
        ledger = json.loads((runner.HERE / "claims.json").read_text())
        self.assertFalse(ledger["independent_review"])
        self.assertEqual(ledger["experimental_observations"], 0)
        self.assertEqual(sum(c["ledger_delta"] for c in ledger["claims"]), 2)

    def test_protocol_prevents_sampled_proof(self):
        protocol = json.loads((runner.HERE / "PROTOCOL.json").read_text())
        self.assertIn("continuous-window", protocol["acceptance"].lower())
        self.assertEqual(protocol["budgets"]["max_variants"], 1)


if __name__ == "__main__":
    unittest.main()
