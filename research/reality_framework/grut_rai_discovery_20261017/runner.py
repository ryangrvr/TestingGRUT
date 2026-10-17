"""Hash-frozen, append-only, resource-bounded GRUT-RAI discovery runner."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OWNED = ["README.md", "CHARTER.md", "PROTOCOL.json", "claims.json",
         "MATHEMATICS.md", "experiment.py", "runner.py", "test_discovery.py"]
PINNED = ["provenance/auditor.py", "provenance/claims.json", "CURRENT_STATE.json",
          "research/reality_framework/grut_p1_preparation_20261017/model.py",
          "research/reality_framework/grut_p1_preparation_20261017/CHARTER.md",
          "research/reality_framework/grut_p1_preparation_20261017/MANIFEST.json"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tracked_manifest():
    return {str(path.relative_to(ROOT)): sha(path)
            for path in [*(HERE / f for f in OWNED), *(ROOT / f for f in PINNED)]}


def verify():
    frozen = json.loads((HERE / "FREEZE.json").read_text())
    current = tracked_manifest()
    if frozen["sha256"] != current:
        raise ValueError("Frozen input changed; create a separately named successor, do not refreeze")
    # Require the freeze and all executable inputs to exist in HEAD unmodified.
    frozen_paths = [*current, str((HERE / "FREEZE.json").relative_to(ROOT))]
    for path in frozen_paths:
        committed = subprocess.run(["git", "-C", str(ROOT), "show", f"HEAD:{path}"],
                                   capture_output=True, check=True).stdout
        if hashlib.sha256(committed).hexdigest() != sha(ROOT / path):
            raise ValueError(f"Uncommitted frozen input: {path}")
    return frozen


def freeze():
    path = HERE / "FREEZE.json"
    with path.open("x") as handle:
        json.dump({"schema": 1, "base_commit": json.loads((HERE / "PROTOCOL.json").read_text())["base_commit"],
                   "sha256": tracked_manifest()}, handle, indent=2)
        handle.write("\n")
    print("Frozen. Commit the packet before running. No experiment executed.")


def limits(budget):
    resource.setrlimit(resource.RLIMIT_CPU, (budget["cpu_seconds"], budget["cpu_seconds"]))
    cap = budget["address_space_MiB"] * 1024**2
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    cap = budget["output_file_MiB"] * 1024**2
    resource.setrlimit(resource.RLIMIT_FSIZE, (cap, cap))


def run(run_id):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", run_id):
        raise ValueError("Run ID must be 1-64 safe letters, digits, hyphens or underscores")
    frozen = verify()
    budget = json.loads((HERE / "PROTOCOL.json").read_text())["budgets"]
    output = HERE / "runs" / run_id
    output.mkdir(parents=True, exist_ok=False)
    metadata = {"run_id": run_id, "status": "started", "freeze_sha256": sha(HERE / "FREEZE.json"),
                "revision": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
                "budgets": budget, "inputs": frozen["sha256"]}
    (output / "START.json").write_text(json.dumps(metadata, indent=2) + "\n")
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1")
    started = time.monotonic()
    status, exitcode = "failed", None
    try:
        with (output / "RESULTS.json").open("x") as stdout, (output / "stderr.txt").open("x") as stderr:
            child = subprocess.run([sys.executable, str(HERE / "experiment.py")], cwd=ROOT, env=env,
                                   stdout=stdout, stderr=stderr, timeout=budget["wall_seconds"],
                                   preexec_fn=lambda: limits(budget))
        exitcode = child.returncode
        if exitcode == 0:
            result = json.loads((output / "RESULTS.json").read_text())
            status = "completed"
            print(result["decision"])
        else:
            if exitcode in (-signal.SIGXCPU, -signal.SIGKILL, -signal.SIGXFSZ):
                status = "resource-stopped"
            print(f"Experiment failed; preserved in {output}", file=sys.stderr)
    except subprocess.TimeoutExpired:
        status = "resource-stopped"
    except Exception as error:
        metadata["runner_error"] = repr(error)
    finally:
        try:
            verify()
        except Exception as error:
            status = "integrity-failed"
            metadata["integrity_error"] = repr(error)
        metadata.update(status=status, exit_code=exitcode, wall_seconds=time.monotonic() - started,
                        result_sha256=sha(output / "RESULTS.json") if (output / "RESULTS.json").exists() else None)
        (output / "END.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return 0 if status == "completed" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["freeze", "verify", "run"])
    parser.add_argument("--id")
    args = parser.parse_args()
    if args.command == "freeze":
        freeze()
    elif args.command == "verify":
        verify()
        print("Frozen committed inputs intact; integrity is not mathematical validation.")
    else:
        if not args.id:
            parser.error("run requires --id")
        return run(args.id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
