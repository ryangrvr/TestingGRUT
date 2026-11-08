"""Frozen D1-01 finite-parent reduction, with constructive branch-error bounds."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.linalg import eigh

from bosonized import build_hamiltonian, compressed_weyl
from continuum import key
from analyze_m1 import distance

HERE=Path(__file__).resolve().parent
PROTOCOLS=[[0.,.7,1.9],[0.,.6,1.5,2.4]]
COUNTS=[4,8,16,32]


def coherent(energies,w,times,k=None):
    dim=len(energies)
    k=dim if k is None else k
    states=np.zeros((k,1),complex);states[0,0]=1
    bounds=np.zeros(1)
    labels=[""]
    extra=[];baseline=[];previous=0.
    e=energies[:k]-energies[0]
    for t in times:
        states*=np.exp(-1j*e*(t-previous))[:,None]
        all_shifted=w[:,:k]@states
        shifted=all_shifted[:k,:]
        q=.5*np.linalg.norm(all_shifted[k:,:],axis=0)
        extra.append(float(2*np.dot(q,q)))
        baseline.append(float((np.linalg.norm(states)**2-np.linalg.norm(all_shifted)**2)/2))
        bounds=np.concatenate((bounds+q,bounds+q))
        states=np.concatenate(((states+shifted)/2,(states-shifted)/2),axis=1)
        labels=[s+"+" for s in labels]+[s+"-" for s in labels]
        previous=t
    norms=np.linalg.norm(states,axis=0)
    p={s:float(norms[j]**2) for j,s in enumerate(labels)}
    p["overflow"]=1-sum(p.values())
    terms=2*norms*bounds+bounds*bounds
    return {"probabilities":p,"tv_bound":min(1.,float(terms.sum())),
            "branch_vector_bounds":dict(zip(labels,map(float,bounds))),
            "additional_projection_loss_by_step":extra,
            "parent_instrument_loss_on_reduced_paths_by_step":baseline,
            "overflow_accounting_residual":float(p["overflow"]-sum(extra)-sum(baseline))}


def dephased(w,times):
    eye=np.eye(len(w))
    plus=np.abs((eye+w)/2)**2
    minus=np.abs((eye-w)/2)**2
    populations=np.zeros((len(w),1));populations[0,0]=1
    labels=[""]
    for _ in times:
        populations=np.concatenate((plus@populations,minus@populations),axis=1)
        labels=[s+"+" for s in labels]+[s+"-" for s in labels]
    p={s:float(populations[:,j].sum()) for j,s in enumerate(labels)}
    p["overflow"]=1-sum(p.values())
    return p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    if subprocess.check_output(["git","status","--porcelain"],cwd=HERE,text=True):
        raise RuntimeError("Execution requires a clean committed tree")
    source=subprocess.check_output(["git","rev-parse","HEAD"],cwd=HERE,text=True).strip()
    began=time.monotonic()
    card=json.loads((HERE/"protocol_m1_01.json").read_text())
    states,omega,h=build_hamiltonian(card,12,128,.5)
    energies,vectors=eigh(h,driver="evd")
    w,_,_,_=compressed_weyl(card,states,omega)
    w=vectors.T@w@vectors
    archive=json.loads((HERE/"raw_m1_01/boson_m05_E12_q128.json").read_text())
    replay={key(t):coherent(energies,w,t) for t in card["protocols"]}
    replay_error=max(distance(v["probabilities"],archive["histories"][k]) for k,v in replay.items())
    if replay_error>=1e-10:
        raise RuntimeError(f"Parent replay fails: {replay_error}")
    records={}
    for times in PROTOCOLS:
        parent=coherent(energies,w,times)
        pop=dephased(w,times)
        reduced={}
        for k in COUNTS:
            v=coherent(energies,w,times,k)
            v["tv_from_parent"]=distance(v["probabilities"],parent["probabilities"])
            v["bound_dominates"]=v["tv_from_parent"]<=v["tv_bound"]+1e-10
            reduced[str(k)]=v
        records[key(times)]={"parent":parent,"dephased":pop,
            "dephased_tv_from_parent":distance(pop,parent["probabilities"]),"reduced":reduced}
    calibration=distance(dephased(w,[0.]),replay["0"]["probabilities"])
    gates={"parent_replay":replay_error<1e-10,"dephased_calibration":calibration<1e-10,
           "coherence_contrast":max(v["dephased_tv_from_parent"] for v in records.values())>=.001,
           "K32_accuracy":all(v["reduced"]["32"]["tv_from_parent"]<=1e-5 for v in records.values()),
           "K32_bound":all(v["reduced"]["32"]["tv_bound"]<=.001 for v in records.values()),
           "all_bounds_dominate":all(r["bound_dominates"] for v in records.values() for r in v["reduced"].values()),
           "overflow_accounting":max(abs(r["overflow_accounting_residual"]) for v in records.values() for r in v["reduced"].values())<1e-10}
    result={"protocol":"D1-01","source_commit":source,"parent_M1_01_cutoff":12,
            "source_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                 ["reduce_d1.py","bosonized.py","protocol_m1_01.json","PROTOCOL_D1_01.md"]},
            "parent_replay_max_tv":replay_error,"dephased_calibration_tv":calibration,
            "records":records,"gates":gates,"seconds":time.monotonic()-began,
            "verdict":"finite_parent_reduction_pass" if all(gates.values()) else "finite_parent_reduction_gate_miss",
            "scope":"No bound on the parent continuum error; M1-01 remains unresolved."}
    with (out/"result.json").open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    print(json.dumps({k:result[k] for k in ["verdict","gates","seconds"]}),flush=True)


if __name__=="__main__":
    raise SystemExit(main())
