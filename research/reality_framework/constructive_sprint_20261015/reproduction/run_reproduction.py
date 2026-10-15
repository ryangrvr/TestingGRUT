"""Execute immutable accessibility packet without cache and retain exact outputs."""
from pathlib import Path
import argparse
import json
import os
import platform
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1] / "cp12_charged_accessibility_gate"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("which", choices=["node", "python"])
    args = parser.parse_args()
    command = ["node", str(PACKET / "verify.mjs")] if args.which == "node" else [sys.executable, str(PACKET / "audit.py")]
    env = os.environ.copy()
    env.pop("CP12_FORMATION_CACHE", None)
    for key in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]:
        env[key] = "1"
    started = time.time()
    clock = time.monotonic()
    with (HERE / f"{args.which}.stdout.log").open("w") as stdout, (HERE / f"{args.which}.stderr.log").open("w") as stderr:
        result = subprocess.run(command, cwd=PACKET, env=env, stdout=stdout, stderr=stderr, check=False)
    metadata = dict(command=command, exit_code=result.returncode,
        runtime_seconds=time.monotonic()-clock, started_unix=started,
        finished_unix=time.time(), python=sys.version, platform=platform.platform(),
        CP12_FORMATION_CACHE="unset", threads=1,
        immutable_source_commit="1801afee27096c03c2a3b64cc811ac0be0d2a4a7")
    (HERE / f"{args.which}.execution.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata, indent=2), flush=True)
    raise SystemExit(result.returncode)

if __name__ == "__main__":
    main()
