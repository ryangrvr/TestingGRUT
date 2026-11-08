"""E1-01: a fixed joint Hamiltonian/instrument correction in the E=12 parent."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import eigh

HERE=Path(__file__).resolve().parent
MATCH=HERE.parent/"grut_detector_matching_20261107"
sys.path.insert(0,str(MATCH))
import bosonized
from analyze_m1 import distance
from continuum import key

TIMES=[[0.,.7,1.9],[0.,.6,1.5,2.4]]
STRENGTHS=[2.,.5]


def save(path,data):
    with path.open("x") as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write("\n")


def integral_gram(energies,matrix,t):
    delta=energies[:,None]-energies[None,:]
    return matrix*t*np.exp(.5j*delta*t)*np.sinc(delta*t/(2*np.pi))


def verify_integral():
    e=np.array([.2,1.3,2.7])
    r=np.array([[.3,.2j,-.7],[.1j,.5,-.2j]])
    z=np.array([.4,.3+.2j,-.1j]);z/=np.linalg.norm(z)
    t=.7
    b=integral_gram(e,r.conj().T@r,t)
    analytic=float(np.vdot(z,b@z).real)
    integral=quad(lambda u:np.linalg.norm(r@(np.exp(-1j*e*u)*z))**2,0,t,epsabs=1e-13)[0]
    norm_integral=quad(lambda u:np.linalg.norm(r@(np.exp(-1j*e*u)*z)),0,t,epsabs=1e-13)[0]
    upper=np.sqrt(t*analytic)
    return {"squared_residual_integral_error":abs(analytic-integral),
            "minimum_gram_eigenvalue":float(eigh(b,eigvals_only=True)[0]),
            "integrated_norm":norm_integral,"cauchy_schwarz_upper":float(upper),
            "pass":bool(abs(analytic-integral)<1e-11 and norm_integral<=upper+1e-12)}


def construct(h,states,omega):
    free=states@omega
    p=np.flatnonzero(free<=6+1e-11)
    q=np.flatnonzero(free>6+1e-11)
    if len(p)!=72:raise RuntimeError("Retained basis dimension is not 72")
    a=h[np.ix_(p,p)]
    b=h[np.ix_(q,p)]
    c=h[np.ix_(q,q)]
    ec,vc=eigh(c,driver="evd")
    cb=vc.T@b
    anchor=float(eigh(a,eigvals_only=True)[0])
    iterations=[]
    converged=False
    for i in range(12):
        denominator=anchor-ec
        gap=float(np.min(np.abs(denominator)))
        if gap<=1e-8:raise RuntimeError("Anchor enters the eliminated-block spectrum")
        t=vc@(cb/denominator[:,None])
        g=np.eye(len(p))+t.T@t
        eg,vg=eigh(g)
        gi=(vg*(eg**-.5))@vg.T
        j=np.zeros((len(h),len(p)))
        j[p]=gi;j[q]=t@gi
        he=j.T@h@j
        new=float(eigh(he,eigvals_only=True)[0])
        iterations.append({"iteration":i,"anchor_in":anchor,"anchor_out":new,"resolvent_distance":gap})
        if abs(new-anchor)<1e-12:
            converged=True;break
        anchor=new
    normt=float(np.sqrt(eigh(t.T@t,eigvals_only=True)[-1]))
    barej=np.eye(len(h))[:,p]
    return {"bare":(barej,a),"joint":(j,he)},p,{
        "iterations":iterations,"converged":converged,"T_operator_norm":normt,
        "isometry_error":float(np.linalg.norm(j.T@j-np.eye(len(p)),ord=2)),
        "retained_dimension":len(p),"eliminated_dimension":len(q),
        "construction_inputs":"parent matrix blocks only; full-parent eigenvectors not used"}


def bounded_histories(h,omega_parent,j,hr,w,wr,protocols):
    er,vr=eigh(hr)
    vacuum=vr[:,0].astype(complex)
    embedded=j@vacuum
    overlap=np.vdot(omega_parent,embedded)
    phase=overlap.conjugate()/abs(overlap) if abs(overlap) else 1.
    vacuum*=phase
    preparation_error=float(np.linalg.norm(omega_parent-j@vacuum))
    residual=h@j-j@hr
    re=residual@vr
    gram=re.conj().T@re
    detector=w@j-j@wr
    defect=np.eye(len(wr))-wr.conj().T@wr
    mindef=float(eigh(defect,eigvals_only=True)[0])
    records={}
    smallest_gram=0.;smallest_quadratic=0.
    for times in protocols:
        branches=vacuum[:,None]
        errors=np.array([preparation_error])
        labels=[""]
        previous=0.
        dynamics=[];instruments=[]
        for moment in times:
            dt=moment-previous
            coefficients=vr.T@branches
            b=integral_gram(er,gram,dt)
            smallest_gram=min(smallest_gram,float(eigh(b,eigvals_only=True)[0]))
            squared=np.real(np.sum(coefficients.conj()*(b@coefficients),axis=0))
            smallest_quadratic=min(smallest_quadratic,float(squared.min()))
            if smallest_gram < -1e-9 or smallest_quadratic < -1e-10:
                raise RuntimeError("Residual integral positivity failed")
            propagation=np.sqrt(dt*np.maximum(0.,squared))
            errors+=propagation
            branches=vr@(np.exp(-1j*(er-er[0])*dt)[:,None]*coefficients)
            q=.5*np.linalg.norm(detector@branches,axis=0)
            dynamics.append(dict(zip(labels,map(float,propagation))))
            instruments.append(dict(zip(labels,map(float,q))))
            shifted=wr@branches
            branches=np.concatenate(((branches+shifted)/2,(branches-shifted)/2),axis=1)
            errors=np.concatenate((errors+q,errors+q))
            labels=[s+"+" for s in labels]+[s+"-" for s in labels]
            previous=moment
        norms=np.linalg.norm(branches,axis=0)
        probabilities={s:float(norms[n]**2) for n,s in enumerate(labels)}
        probabilities["overflow"]=1-sum(probabilities.values())
        records[key(times)]={"probabilities":probabilities,
            "bound":min(1.,float(np.sum(2*norms*errors+errors**2))),
            "branch_vector_bounds":dict(zip(labels,map(float,errors))),
            "propagation_bounds_by_prefix":dynamics,"detector_errors_by_prefix":instruments}
    return {"histories":records,"preparation_vector_error":preparation_error,
            "ground_energy":float(er[0]),"gap":float(er[1]-er[0]),
            "ground_parent_residual":float(np.linalg.norm(h@(j@vacuum)-er[0]*(j@vacuum))),
            "dynamical_residual_norm":float(np.sqrt(max(0.,eigh(gram,eigvals_only=True)[-1]))),
            "hermiticity_residual":float(np.linalg.norm(hr-hr.T,ord=np.inf)),
            "minimum_instrument_defect":mindef,"minimum_integral_gram_eigenvalue":smallest_gram,
            "minimum_integral_quadratic_form":smallest_quadratic}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    if subprocess.check_output(["git","status","--porcelain"],cwd=HERE,text=True):
        raise RuntimeError("E1-01 requires a clean committed tree")
    start=time.monotonic()
    commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=HERE,text=True).strip()
    provenance={"protocol":"E1-01","source_commit":commit,"python":platform.python_version(),
        "numpy":np.__version__,"scipy":scipy.__version__,"source_sha256":{
            str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [HERE/"experiment.py",HERE/"DERIVATION.md",HERE/"PROTOCOL_E1_01.md",MATCH/"bosonized.py",MATCH/"continuum.py",MATCH/"protocol_m1_01.json"]}}
    save(out/"provenance.json",provenance)
    verification=verify_integral();save(out/"integral_verification.json",verification)
    if not verification["pass"]:raise RuntimeError("Integral implementation verification failed")
    card=json.loads((MATCH/"protocol_m1_01.json").read_text())
    states,omega,h=bosonized.build_hamiltonian(card,12,128,.5)
    schemes,p,construction=construct(h,states,omega)
    save(out/"construction.json",construction)
    if not construction["converged"]:raise RuntimeError("Anchor did not converge within 12 iterations")
    # Only now construct the full-parent spectral reference.
    ep,vp=eigh(h,driver="evd")
    parent_ground=vp[:,0].astype(complex)
    old=json.loads((MATCH/"raw_d1_01/result.json").read_text())
    all_results={}
    for coupling in STRENGTHS:
        workcard=dict(card,coupling=coupling)
        w,_,_,_=bosonized.compressed_weyl(workcard,states,omega)
        parent,loss=bosonized.histories(ep,vp,w,TIMES,0.)
        replay=None
        if coupling==2:
            replay=max(distance(parent[k],old["records"][k]["parent"]["probabilities"]) for k in parent)
            if replay>=1e-10:raise RuntimeError("Strong parent does not reproduce D1")
        rows={}
        barew=w[np.ix_(p,p)]
        for name in ["bare","joint","hamiltonian_only"]:
            j,hr=schemes["bare" if name=="bare" else "joint"]
            wr=j.T@w@j if name=="joint" else barew
            data=bounded_histories(h,parent_ground,j,hr,w,wr,TIMES)
            for k,record in data["histories"].items():
                record["tv"]=distance(record["probabilities"],parent[k])
                record["bound_dominates"]=record["tv"]<=record["bound"]+1e-10
            data["ground_energy_error"]=data["ground_energy"]-float(ep[0])
            data["gap_error"]=data["gap"]-float(ep[1]-ep[0])
            data["max_tv"]=max(v["tv"] for v in data["histories"].values())
            data["max_bound"]=max(v["bound"] for v in data["histories"].values())
            rows[name]=data
        gates={"construction":construction["isometry_error"]<1e-10 and rows["joint"]["preparation_vector_error"]<1e-8 and rows["joint"]["ground_parent_residual"]<1e-8,
            "numerics":all(r["hermiticity_residual"]<1e-10 and r["minimum_instrument_defect"]>=-1e-10 and all(min(x["probabilities"].values())>=-1e-10 and abs(sum(x["probabilities"].values())-1)<1e-9 for x in r["histories"].values()) for r in rows.values()),
            "method_improvement":rows["joint"]["max_tv"]<=.5*rows["bare"]["max_tv"] and rows["joint"]["max_tv"]<=rows["hamiltonian_only"]["max_tv"],
            "finite_parent_accuracy":rows["joint"]["max_tv"]<=.001,
            "finite_parent_certification":rows["joint"]["max_bound"]<=.001,
            "bounds_dominate":all(x["bound_dominates"] for r in rows.values() for x in r["histories"].values())}
        result={"coupling":coupling,"parent_histories":parent,"parent_overflow_by_step":loss,
                "strong_parent_replay_max_tv":replay,"schemes":rows,"gates":gates}
        save(out/f"coupling_{coupling:g}.json",result)
        all_results[f"{coupling:g}"]={"gates":gates,"max_tv":{n:r["max_tv"] for n,r in rows.items()},
             "max_bound":{n:r["max_bound"] for n,r in rows.items()}}
        print(json.dumps({"coupling":coupling,**all_results[f"{coupling:g}"]}),flush=True)
    np.savez_compressed(out/"effective_checkpoint.npz",parent_energies=ep,retained_indices=p,
                        embedding=schemes["joint"][0],effective_hamiltonian=schemes["joint"][1])
    summary={"protocol":"E1-01","results":all_results,"construction":construction,
             "seconds":time.monotonic()-start,"scope":"Only the finite E=12 parent; M1 continuum gate remains unresolved."}
    save(out/"summary.json",summary)
    print(json.dumps({"complete":True,"seconds":summary["seconds"]}),flush=True)


if __name__=="__main__":main()
