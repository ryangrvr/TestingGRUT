"""Reproduce or verify checked-in control results; fail on stale code or data."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
import unittest
from controls import run_controls

ROOT = Path(__file__).resolve().parent


def encode(value):
    if isinstance(value, Fraction):
        return {"rational": [value.numerator, value.denominator]}
    raise TypeError(type(value))


def result():
    return {"schema": "FOUNDATIONS_STANDARD_CONTROLS_V1",
            "scope": "ISSUE_3_STANDARD_CONTROLS_ONLY_NO_CANDIDATE_NO_CARD",
            "review": "SAME_AUTHOR_CHECKS_EXTERNAL_REVIEW_PENDING",
            "python": platform.python_version(),
            "code_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                            for p in ["controls.py", "test_controls.py", "run.py"]},
            "controls": run_controls()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Verify checked-in output, do not overwrite")
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.discover(str(ROOT), pattern="test_controls.py")
    ran = unittest.TextTestRunner(verbosity=2).run(suite)
    if not ran.wasSuccessful():
        return 1
    output = json.dumps(result(), default=encode, indent=2, sort_keys=True, allow_nan=False)+"\n"
    path = ROOT/"RESULTS.json"
    if args.check:
        if not path.is_file() or path.read_text() != output:
            print("FAIL: missing/stale/different RESULTS.json (including Python version or code hashes)", file=sys.stderr)
            return 1
        print("13 standard/control groups reproduced; checked-in evidence matches. NO GRUT BANK DECISION.")
    else:
        path.write_text(output)
        print("Wrote RESULTS.json; same-author reproduction, no scientific approval.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
