"""Explicit first capture of issue #2's observational baseline. Never auto-refresh."""
import argparse
from datetime import datetime, timezone
import json
import hashlib
from pathlib import Path
import subprocess
import sys

from integrity import read_junit

PINNED_SOURCE = "114af03f5a9851ccc7ae5d681e2457203d192e80"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--junit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    capture = Path(__file__).with_name("capture_state.py")
    state = json.loads(subprocess.check_output([sys.executable, str(capture), "--root", str(args.root)]))
    if state["source_commit"] != PINNED_SOURCE or state["working_tree_dirty"]:
        raise ValueError("Capture requires the pinned source in a pristine reference worktree")
    cases, counts = read_junit(args.junit)
    if counts != {"tests": 242, "failures": 9, "errors": 0, "skipped": 1}:
        raise ValueError("Raw suite does not reproduce the owner's pinned 242-case baseline")
    out = {"schema_version": 1, "authorization": "OBSERVATIONS_ONLY_NO_ADJUDICATION",
           "source_commit": PINNED_SOURCE, "utc": datetime.now(timezone.utc).isoformat(),
           "pytest_cases": cases, "counts": counts,
           "prior_hosted_run": {"run": 37924406643, "artifact": 11613965534,
               "artifact_digest": "6a5f12772c48e95884797c0b3ccfbf9ba14c9577c3f09b7254634f702368bd2c",
               "python_reported_by_owner": "3.12.15"},
           **{key: state[key] for key in ("declarations", "open_passes", "bank_inventory",
                                          "protected_inputs", "tool_versions")}}
    lock = args.output.with_suffix(args.output.suffix + ".sha256")
    if args.output.exists() or lock.exists():
        raise FileExistsError("A frozen manifest/lock already exists; refresh requires explicit review")
    payload = json.dumps(out, indent=2, sort_keys=True) + "\n"
    with args.output.open("x", encoding="utf-8") as target:
        target.write(payload)
    with lock.open("x", encoding="utf-8") as target:
        target.write(hashlib.sha256(payload.encode()).hexdigest() + "\n")
    print(f"Frozen observational manifest: {args.output}")


if __name__ == "__main__":
    main()
