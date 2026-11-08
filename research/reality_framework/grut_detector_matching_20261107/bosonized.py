"""M1-01 independent interval Hamiltonian truncation.

No lattice Hamiltonian, state, eigenvalue, or measured covariance is imported.
Vertex operators use infinite-line normal ordering converted to the Dirichlet
vacuum. Detector Weyl operators are compressed, with explicit overflow.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh, expm
from scipy.special import eval_genlaguerre, gammaln, k0, roots_legendre

from continuum import interval_spectrum, all_histories, key

HERE = Path(__file__).resolve().parent


def basis(length, cutoff):
    modes = int(np.floor(length * np.sqrt(max(0, cutoff**2 - 1/np.pi)) / np.pi))
    omega = np.sqrt((np.arange(1, modes+1)*np.pi/length)**2 + 1/np.pi)
    states, energies = [], []
    def visit(j, left, occupations):
        if j == modes:
            states.append(occupations)
            energies.append(cutoff-left)
            return
        for n in range(int((left+1e-12)/omega[j])+1):
            visit(j+1, left-n*omega[j], occupations+[n])
    visit(0, cutoff, [])
    order = np.argsort(energies, kind="stable")
    return np.asarray(states, dtype=int)[order], np.asarray(energies)[order], omega


def normal_vertex_one(lam, maximum):
    """<n|exp(i lam a^dagger) exp(i lam a)|m>, independently of truncation."""
    n, m = np.indices((maximum+1, maximum+1))
    lo, hi = np.minimum(n,m), np.maximum(n,m)
    d = hi-lo
    return np.exp((gammaln(lo+1)-gammaln(hi+1))/2) * (1j*lam)**d * eval_genlaguerre(lo,d,lam*lam)


def vertex(states, lambdas):
    result = np.ones((len(states),len(states)),complex)
    for j, lam in enumerate(lambdas):
        n = states[:,j]
        if n.max():
            v = normal_vertex_one(lam, int(n.max()))
            result *= v[n[:,None],n[None,:]]
    return result


def log_normal_factor(z, length, images):
    mu = 1/np.sqrt(np.pi)
    js = np.arange(-images,images+1)
    return float(np.sum(k0(2*mu*np.abs(z+js*length)))
                 -2*np.sum(k0(2*mu*np.arange(1,images+1)*length)))


def log_normal_integral(z, length):
    """Independent integral form, u=cosh(t), stable near the endpoints."""
    mu=1/np.sqrt(np.pi)
    def fun(t):
        u=np.cosh(t)
        q=np.exp(-2*mu*length*u)
        return (-2*q+np.exp(-2*mu*z*u)+np.exp(-2*mu*(length-z)*u))/(1-q)
    return quad(fun,0,30,epsabs=1e-11,epsrel=1e-11)[0]


def build_hamiltonian(card, cutoff, quadrature, mass):
    states, energies, omega = basis(card["length"],cutoff)
    if len(states)>card["resources"]["dimension_cap"]:
        raise RuntimeError("Registered dimension cap exceeded")
    h = np.diag(energies)
    if mass:
        nodes, weights = roots_legendre(quadrature)
        length=card["length"]
        z=(nodes+1)*length/2
        weights*=length/2
        mu=1/np.sqrt(np.pi)
        beta=2*np.sqrt(np.pi)
        profile=np.sqrt(np.pi)*np.cosh(mu*(z-length/2))/(4*np.cosh(mu*length/2))
        for x,w,cl in zip(z,weights,profile):
            u=np.sqrt(2/length)*np.sin(np.arange(1,len(omega)+1)*np.pi*x/length)
            v=vertex(states,beta*u/np.sqrt(2*omega))
            cosine=np.real(np.exp(1j*beta*cl)*v)
            h-=w*mass*np.exp(np.euler_gamma)/(2*np.pi**1.5)*np.exp(
                log_normal_factor(x,length,card["normal_order_images"]))*cosine
    return states, omega, h


def compressed_weyl(card, states, cutoff_omega):
    omega,h2,fn,_,meta=interval_spectrum(card,card["reference_modes"],card["quadrature_epsabs"])
    # The unoccupied higher modes contribute their vacuum factor. They are
    # not removed by normalizing the projected operator to be unitary.
    hs=fn[:len(cutoff_omega)]/np.sqrt(2*np.pi*cutoff_omega)
    c=card["coupling"]
    w=vertex(states,c*hs)*np.exp(1j*c*meta["mean"]-c*c*meta["variance"]/2)
    return w, omega, h2, meta


def histories(energies, vectors, w, protocols, phase):
    vacuum=vectors[:,0].astype(complex)
    e=energies-energies[0]
    result={}
    loss_by_step={}
    for times in protocols:
        branches=vacuum[:,None]
        labels=[""]
        previous=0.
        losses=[]
        for t in times:
            if t!=previous:
                branches=vectors@(np.exp(-1j*e*(t-previous))[:,None]*(vectors.T@branches))
            shifted=np.exp(1j*phase)*(w@branches)
            loss=float((np.linalg.norm(branches)**2-np.linalg.norm(shifted)**2)/2)
            if loss < -1e-10:
                raise RuntimeError("Negative instrument overflow beyond tolerance")
            losses.append(loss)
            branches=np.concatenate(((branches+shifted)/2,(branches-shifted)/2),axis=1)
            labels=[s+"+" for s in labels]+[s+"-" for s in labels]
            previous=t
        p={s:float(np.vdot(branches[:,j],branches[:,j]).real) for j,s in enumerate(labels)}
        p["overflow"]=sum(losses)
        result[key(times)]=p
        loss_by_step[key(times)]=losses
    return result,loss_by_step


def run(card,cutoff,quadrature,mass):
    start=time.monotonic()
    states,omega,h=build_hamiltonian(card,cutoff,quadrature,mass)
    herm=float(np.linalg.norm(h-h.T,ord=np.inf))
    energies,vectors=eigh(h,driver="evd")
    residual=max(float(np.linalg.norm(h@vectors[:,j]-energies[j]*vectors[:,j])) for j in range(2))
    w,all_omega,h2,meta=compressed_weyl(card,states,omega)
    defect=np.eye(len(w))-w.conj().T@w
    mineig=float(eigh(defect,eigvals_only=True,subset_by_index=[0,0])[0])
    if mineig < -card["gates"]["negative_defect"]:
        raise RuntimeError(f"Compressed instrument is not a contraction: {mineig}")
    p,loss=histories(energies,vectors,w,card["protocols"],card["phase"])
    reference,_=all_histories(card,all_omega,h2,meta["mean"])
    free_tv={k:.5*(sum(abs(v[s]-reference[k][s]) for s in reference[k])+abs(v["overflow"]))
             for k,v in p.items()}
    return {"protocol":card["id"],"cutoff":cutoff,"quadrature":quadrature,"mass_R_over_g":mass,
            "dimension":len(states),"occupied_modes":len(omega),"energies_low":energies[:6].tolist(),
            "gap":float(energies[1]-energies[0]),"hermiticity_residual":herm,
            "eigen_residual":residual,"minimum_defect_eigenvalue":mineig,
            "histories":p,"overflow_by_step":loss,"free_reference_tv":free_tv,
            "free_reference_tv_role":"verification for m=0; physical contrast only for m>0",
            "probability_norm_error":max(abs(sum(q.values())-1) for q in p.values()),
            "detector_reference_metadata":meta,"seconds":time.monotonic()-start},states,energies,vectors[:,0]


def verify():
    # Matrix elements are checked against a much larger direct oscillator.
    a=np.diag(np.sqrt(np.arange(1,80)),1)
    lam=.37
    direct=expm(1j*lam*(a+a.T))*np.exp(lam*lam/2)
    small=normal_vertex_one(lam,8)
    sites=[.02,.2,1.,2.5,4.8]
    return {"vertex_max_error":float(np.max(np.abs(direct[:9,:9]-small))),
            "normal_factor_images_vs_integral":max(abs(log_normal_factor(z,5,16)-log_normal_integral(z,5)) for z in sites),
            "normal_factor_image_doubling":max(abs(log_normal_factor(z,5,8)-log_normal_factor(z,5,16)) for z in sites),
            "basis_dimensions":{str(c):len(basis(5,c)[0]) for c in (4,6,8,10,12)}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--verify",action="store_true")
    parser.add_argument("--cutoff",type=float,default=4)
    parser.add_argument("--mass",type=float,default=.5)
    parser.add_argument("--quadrature",type=int,default=128)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix(".npz").exists():
        raise FileExistsError("Refusing to replace a preserved result")
    card=json.loads((HERE/"protocol_m1_01.json").read_text())
    if args.verify:
        data=verify()
    else:
        data,states,energies,ground=run(card,args.cutoff,args.quadrature,args.mass)
        np.savez_compressed(args.output.with_suffix(".npz"),basis=states,energies=energies,ground=ground)
    data["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open("x") as f:
        json.dump(data,f,indent=2,allow_nan=False)
        f.write("\n")
    print(json.dumps({k:v for k,v in data.items() if k in ["cutoff","mass_R_over_g","dimension","gap","seconds","vertex_max_error","normal_factor_images_vs_integral"]}),flush=True)
