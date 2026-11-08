"""Run P1 from a clean committed checkout, preserving resource stops."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
from analyze import analyze
from point import write

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/"grut_detector_matching_20261107"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(subprocess.check_output(["git","rev-parse","--show-toplevel"],cwd=HERE,text=True).strip())
    out = args.output_dir.resolve()
    if out.is_relative_to(root):
        raise ValueError("Production outputs must remain outside the frozen checkout")
    if subprocess.check_output(["git","status","--porcelain"],cwd=HERE,text=True):
        raise ValueError("Clean committed source required")
    out.mkdir(parents=True, exist_ok=False)
    card = json.loads((HERE/"protocol_p1_01.json").read_text())
    sources = sorted(HERE.glob("*.py"))+sorted(HERE.glob("*.md"))+[HERE/"protocol_p1_01.json"]
    sources += [OLD/n for n in ("bosonized.py","continuum.py","protocol_m1_01.json")]
    sources += [OLD/"raw_m1_01"/n for n in (
        "boson_m05_E10_q128.json","boson_m05_E12_q128.json","boson_m05_E12_q128.npz","lattice_N20.json")]
    write(out/"provenance.json", {
        "card":card, "source_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=HERE,text=True).strip(),
        "clean_source":True, "python":platform.python_version(), "numpy":np.__version__, "scipy":scipy.__version__,
        "source_and_input_sha256":{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    })
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    began, events = time.monotonic(), []
    for point in card["points"]:
        label = point["label"]
        remaining = card["total_wall_seconds"]-(time.monotonic()-began)
        event = {"label":label}
        start = time.monotonic()
        if remaining <= 0:
            event["status"] = "total_wall_cap"
        else:
            command = [sys.executable,str(HERE/"point.py"),"--label",label,"--output",str(out/(label+".json"))]
            try:
                r = subprocess.run(command,env=env,capture_output=True,text=True,
                    timeout=min(remaining,card["point_wall_seconds"]))
                event.update(status="complete" if r.returncode==0 else "exception",
                             returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
            except subprocess.TimeoutExpired as exc:
                event.update(status="point_wall_cap",stdout=str(exc.stdout),stderr=str(exc.stderr))
        event["elapsed_seconds"] = time.monotonic()-start
        events.append(event)
        write(out/("event_"+label+".json"),event)
        print(json.dumps(event),flush=True)
    write(out/"events.json",events)
    result = analyze(out)
    write(out/"analysis.json",result)
    print(json.dumps(result),flush=True)


if __name__ == "__main__":
    main()
