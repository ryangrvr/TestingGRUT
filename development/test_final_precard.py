"""Hostile tests of packet provenance and failure reporting, not physics."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from final_precard import run_review as R


class FinalPacketIntegrity(unittest.TestCase):
    def test_wrong_input_bytes_reject_before_subject_import(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'bad.zip';path.write_bytes(b'not the pinned input archive')
            with patch('development.cr5_review.check.load') as loader:
                with self.assertRaises(ValueError): R.verified_inputs(path)
                loader.assert_not_called()

    def test_relocked_oracle_file_cannot_pass_original_freeze(self):
        read=Path.read_bytes
        def hostile(path):
            data=read(path)
            return data+b'\n# injected comparison-dependent rule\n' if path==R.HERE/'price_oracle.py' else data
        with patch.object(Path,'read_bytes',hostile):
            with self.assertRaises(ValueError): R.verify_oracle_freeze()

    def test_failed_command_never_leaves_READY_result(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t);bad=t/'bad.zip';bad.write_bytes(b'wrong')
            out=t/'fresh'
            result=subprocess.run([sys.executable,'-m','development.final_precard.run_review',
                '--inputs',str(bad),'--output',str(out)],cwd=R.ROOT,capture_output=True,text=True)
            self.assertEqual(result.returncode,1)
            self.assertFalse((out/'RESULTS.json').exists())
            self.assertTrue(json.loads((out/'BLOCKED.json').read_text())['status'].startswith('BLOCKED'))


if __name__=='__main__': unittest.main()
