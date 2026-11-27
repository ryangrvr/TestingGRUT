"""Verify byte-for-byte replay and exact agreement on frozen F5 observables."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def replay(program, archive):
    raw = subprocess.check_output([sys.executable, str(HERE / program)])
    if raw != (HERE / archive).read_bytes():
        raise SystemExit(f"Replay mismatch for {archive}")
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def main():
    params = json.loads((HERE / "PARAMETERS.json").read_text())
    production, prod_hash = replay("run_f5.py", "RESULTS.json")
    independent, ind_hash = replay("independent_check.py", "INDEPENDENT_RESULTS.json")
    expected = params["expected"]
    for label, results in [("production", production), ("independent", independent)]:
        actual = results["observables"]
        if actual.keys() != expected.keys():
            raise SystemExit(f"Wrong observable coverage in {label}")
        wrong = [name for name, value in actual.items() if value != expected[name]]
        if wrong:
            raise SystemExit(f"{label} differs from the frozen values: {', '.join(wrong)}")
    print(f"{len(expected)} registered quantities agree exactly in both implementations.")
    print(f"Production replay  {prod_hash}")
    print(f"Independent replay {ind_hash}")


if __name__ == "__main__":
    main()
