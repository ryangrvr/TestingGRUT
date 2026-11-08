"""Execute M1-01 from a clean checkout, retaining individual failures and caps."""
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

from analyze_m1 import analyze

HERE=Path(__file__).resolve().parent


def write(path,value):
    with path.open("x") as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write("\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    card=json.loads((HERE/"protocol_m1_01.json").read_text())
    commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=HERE,text=True).strip()
    dirty=subprocess.check_output(["git","status","--porcelain"],cwd=HERE,text=True)
    if dirty:
        raise RuntimeError("M1-01 execution requires a clean committed tree")
    write(out/"provenance.json",{"card":card,"source_commit":commit,"source_tree_clean":True,
          "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
          "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted(HERE.iterdir()) if p.suffix in [".py",".md",".json"]}})
    tasks=[("verification",["bosonized.py","--verify"])]
    for e in card["free_checks"]:
        tasks.append((f"boson_m0_E{e}_q128",["bosonized.py","--mass","0","--cutoff",str(e)]))
    for e in card["energy_cutoffs"]:
        tasks.append((f"boson_m05_E{e}_q128",["bosonized.py","--mass","0.5","--cutoff",str(e)]))
    tasks.append(("boson_m05_E10_q256",["bosonized.py","--mass","0.5","--cutoff","10","--quadrature","256"]))
    for n in card["lattice_N"]:
        tasks.append((f"lattice_N{n}",["lattice_massive.py","--N",str(n)]))
    tasks.append(("lattice_N14_refined",["lattice_massive.py","--N","14","--refined"]))
    began=time.monotonic()
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1")
    events=[]
    for label,arguments in tasks:
        remaining=card["resources"]["total_wall_seconds"]-(time.monotonic()-began)
        event={"label":label}
        if remaining<=0:
            event["status"]="total_wall_cap"
        else:
            command=[sys.executable,str(HERE/arguments[0]),*arguments[1:],"--output",str(out/(label+".json"))]
            try:
                p=subprocess.run(command,env=env,capture_output=True,text=True,
                    timeout=min(remaining,card["resources"]["point_wall_seconds"]))
                event.update(status="complete" if p.returncode==0 else "exception",returncode=p.returncode,
                             stdout=p.stdout,stderr=p.stderr)
            except subprocess.TimeoutExpired as exc:
                event.update(status="point_wall_cap",stdout=str(exc.stdout),stderr=str(exc.stderr))
        events.append(event)
        write(out/("event_"+label+".json"),event)
        print(json.dumps(event),flush=True)
        if label=="verification":
            if event["status"]!="complete":
                break
            v=json.loads((out/"verification.json").read_text())
            if v["vertex_max_error"]>=card["gates"]["vertex_error"] or max(v["normal_factor_images_vs_integral"],v["normal_factor_image_doubling"])>=card["gates"]["normal_factor_error"]:
                break
    write(out/"events.json",events)
    result=analyze(out,card)
    write(out/"analysis.json",result)
    print(json.dumps({"verdict":result["verdict"],"gates":result.get("gates")}),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
