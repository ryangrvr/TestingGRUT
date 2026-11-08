"""Execute the separately frozen M0-02 boundary-hypothesis comparison, preserving failures."""
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

from continuum import all_histories, interval_spectrum, tv, verify_reference
from analyze import analyze

HERE = Path(__file__).resolve().parent


def write(path, data):
    with path.open("x") as f:
        json.dump(data,f,indent=2,allow_nan=False)
        f.write("\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    card_path=HERE/"protocol_m0_02.json"
    card=json.loads(card_path.read_text())
    try:
        commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=HERE,text=True).strip()
    except subprocess.CalledProcessError:
        commit=None
    metadata={"protocol":card,"source_commit":commit,"python":platform.python_version(),
              "numpy":np.__version__,"scipy":scipy.__version__,
              "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(HERE.iterdir()) if p.suffix in (".py",".json",".md")}}
    write(out/"provenance.json",metadata)
    verification=verify_reference(card)
    write(out/"oscillator_verification.json",verification)
    if max(verification["tv_gaussian_vs_high"],verification["tv_cutoff_change"])>=card["gates"]["oscillator_tv"]:
        print("STOP: independent oscillator verification failed",flush=True)
        return 1
    ref_start=time.monotonic()
    omega,h2,fn,fnerr,spectrum=interval_spectrum(card,card["modes"],card["quadrature_epsabs"])
    histories,imaginary=all_histories(card,omega,h2,spectrum["mean"])
    wrong_mean,_=all_histories(card,omega,h2,0.0)
    spectrum["covariance"]={f"{t:g}":[float(z.real),float(z.imag)]
        for t in (0.0,1.0,2.0) for z in [np.sum(h2*np.exp(-1j*omega*t))]}
    controls={}
    for label,modes,eps in (("mode_doubling_tv",card["refined_modes"],card["quadrature_epsabs"]),
                            ("quadrature_tightening_tv",card["modes"],card["quadrature_refined_epsabs"])):
        o,h,_,_,s=interval_spectrum(card,modes,eps)
        p,_=all_histories(card,o,h,s["mean"])
        controls[label]=max(tv(histories[k],p[k]) for k in histories)
    nmax=max(map(len,card["protocols"]))
    reference={"spectrum":spectrum,"histories":histories,"wrong_zero_mean_histories":wrong_mean,
               "wrong_mean_tv_by_protocol":{k:tv(p,wrong_mean[k]) for k,p in histories.items()},
               "controls":controls,"oscillator_verification":verification,
               "maximum_history_tail_bound":2**nmax*nmax**2*card["coupling"]**2*spectrum["tail_variance_bound"],
               "imaginary_residual":imaginary,"seconds":time.monotonic()-ref_start}
    write(out/"reference.json",reference)
    np.savez_compressed(out/"continuum_modes.npz",omega=omega,h2=h2,coefficients=fn,
                        quadrature_error_estimates=fnerr)
    print(json.dumps({"stage":"reference_complete","mean":spectrum["mean"],
                      "variance":spectrum["variance"],"seconds":reference["seconds"]}),flush=True)
    tasks=[(n,"primary",[]) for n in card["primary_N"]]
    tasks += [(n,"smoothed",["--smoothed"]) for n in card["smoothed_control_N"]]
    tasks += [(card["refined_N"],"refined",["--refined"])]
    began=time.monotonic()
    events=[]
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1")
    for n,label,flags in tasks:
        remaining=card["resources"]["total_lattice_wall_seconds"]-(time.monotonic()-began)
        event={"N":n,"label":label}
        if remaining<=0:
            event["status"]="total_wall_cap"
            events.append(event)
            continue
        command=[sys.executable,str(HERE/"lattice_boundary.py"),"--N",str(n),"--output",
                 str(out/f"point_N{n}_{label}.json")]+flags
        try:
            completed=subprocess.run(command,env=env,capture_output=True,text=True,
                timeout=min(remaining,card["resources"]["point_wall_seconds"]))
            event.update(status="complete" if completed.returncode==0 else "exception",
                         returncode=completed.returncode,stdout=completed.stdout,stderr=completed.stderr)
        except subprocess.TimeoutExpired as exc:
            event.update(status="point_wall_cap",stdout=str(exc.stdout),stderr=str(exc.stderr))
        events.append(event)
        write(out/f"event_N{n}_{label}.json",event)
        print(json.dumps(event),flush=True)
    write(out/"execution_events.json",events)
    analysis=analyze(out,card)
    write(out/"analysis.json",analysis)
    print(json.dumps({"verdict":analysis["verdict"],"gates":analysis["gates"]}),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
