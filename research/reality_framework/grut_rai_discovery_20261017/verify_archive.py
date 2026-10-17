"""Read-only verification of archived results, frozen source and preservation."""
import hashlib
import json
from pathlib import Path
import subprocess

import runner


def main():
    frozen = runner.verify()
    base = frozen["base_commit"]
    changed = subprocess.check_output(
        ["git", "-C", str(runner.ROOT), "diff", "--name-only", base, "HEAD"],
        text=True).splitlines()
    prefix = str(runner.HERE.relative_to(runner.ROOT)) + "/"
    if any(not name.startswith(prefix) for name in changed):
        raise ValueError("Historical files changed relative to the P1 base")
    count = 0
    for directory in sorted((runner.HERE / "runs").iterdir()):
        if not directory.is_dir():
            continue
        start = json.loads((directory / "START.json").read_text())
        end = json.loads((directory / "END.json").read_text())
        if start["run_id"] != directory.name or end["run_id"] != directory.name:
            raise ValueError("Archive ID mismatch")
        if start["inputs"] != frozen["sha256"] or end["inputs"] != frozen["sha256"]:
            raise ValueError("Archive input manifest mismatch")
        if start["freeze_sha256"] != runner.sha(runner.HERE / "FREEZE.json"):
            raise ValueError("Archive freeze mismatch")
        if end["result_sha256"] != runner.sha(directory / "RESULTS.json"):
            raise ValueError("Archived result mutated")
        if end["status"] == "completed":
            result = json.loads((directory / "RESULTS.json").read_text())
            if result["independent_validation"] or result["experimental_observations"]:
                raise ValueError("Archive claims unauthorized validation")
        count += 1
    print(json.dumps({"frozen_inputs": "pass", "historical_preservation": "pass",
                      "archived_runs": count, "certifies_physical_truth": False}))


if __name__ == "__main__":
    main()
