"""B1: fixed lattice Hamiltonian, newly derived detector, explicit backaction."""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.sparse.linalg import eigsh,expm_multiply
from dictionary import weights,references,gaussian_from_cov,disturbance,key,OLD

sys.path.insert(0,str(OLD))
from lattice import schwinger

HERE=Path(__file__).resolve().parent


def write(path,data):
    with path.open('x') as f: json.dump(data,f,indent=2,allow_nan=False);f.write('\n')


def run(card,n,mass,refined):
    start=time.monotonic(); a,edges,v,b,shift=weights(card,n)
    w,electric,states=schwinger.hamiltonian(n,a**-2,mass-a/8)
    h=w*a/2
    tol=card['eigensolver_refined_tol'] if refined else card['eigensolver_tol']
    initial=np.random.default_rng(card['eigensolver_seed_base']+n).standard_normal(h.shape[0])
    e,vec=eigsh(h,k=2,which='SA',tol=tol,v0=initial,maxiter=10000)
    order=np.argsort(e);e=e[order];vec=vec[:,order];ground=vec[:,0].astype(complex)
    residual=max(float(np.linalg.norm(h@vec[:,j]-e[j]*vec[:,j])) for j in range(2))
    z=(edges[:-1]+edges[1:])/2
    oldweights=np.exp(-.5*((z-card['length']/2)/card['sigma'])**2);oldweights/=oldweights.sum()
    oldF=electric@oldweights
    F=electric@v+shift
    # Verify the exact Gauss-law rewrite using independently extracted spin bits.
    Z=2*schwinger.bits(states,n)-1
    S=np.cumsum(Z,axis=1)[:,:-1]/2
    identity=float(np.max(abs(F-(.25+S@v))))
    trace=float(h.diagonal().sum()); steps=2 if refined else 1
    def evolve(x,dt):
        for _ in range(steps):
            if dt: x=expm_multiply(-1j*(dt/steps)*h,x,traceA=-1j*(dt/steps)*trace)
        return x
    norm_error=0.
    def histories(field):
        nonlocal norm_error
        phase=np.exp(1j*(card['phase']+card['coupling']*field))[:,None]
        result={}
        for times in card['protocols']:
            branch=ground[:,None];labels=[''];previous=0.
            for t in times:
                before=np.linalg.norm(branch)**2
                branch=evolve(branch,t-previous)
                norm_error=max(norm_error,abs(float(np.linalg.norm(branch)**2-before)))
                shifted=phase*branch
                branch=np.concatenate(((branch+shifted)/2,(branch-shifted)/2),axis=1)
                labels=[s+'+' for s in labels]+[s+'-' for s in labels];previous=t
            result[key(times)]={s:float(np.vdot(branch[:,j],branch[:,j]).real) for j,s in enumerate(labels)}
        return result
    new,old=histories(F),histories(oldF)
    mean=float(np.vdot(ground,F*ground).real);delta=(F-mean)*ground
    cov={str(t):complex(np.vdot(delta,evolve(delta,t))*np.exp(1j*e[0]*t)) for t in (0,1,2)}
    result={'N':n,'a':a,'dimension':len(ground),'mass_R_over_g':mass,'mass_lat_over_g':mass-a/8,
        'refined':refined,'ground_energy':float(e[0]),'gap':float(e[1]-e[0]),
        'new_histories':new,'raw_replay_histories':old,'mean':mean,'covariance':{k:[z.real,z.imag] for k,z in cov.items()},
        'disturbance':disturbance(new),'weights':v.tolist(),'boundary_weight':b,'constant_shift':shift,
        'gauss_detector_identity_error':identity,'eigen_residual':residual,'free_norm_error':norm_error,
        'probability_norm_error':max(abs(sum(p.values())-1) for qs in (new,old) for p in qs.values())}
    if mass==0:
        reference,cc=references(card,n)
        glat=gaussian_from_cov(card,mean,cov)
        static={k:z*cov['0'].real/cc['0'].real for k,z in cc.items()}
        gstatic=gaussian_from_cov(card,mean,static)
        result.update(reference=reference,gaussian_lattice_summary=glat,gaussian_static_summary=gstatic)
    result['seconds']=time.monotonic()-start
    return result,{'ground':ground,'eigenvalues':e,'states':states,'cell_weights':v}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--N',type=int,required=True)
    p.add_argument('--mass',type=float,default=0);p.add_argument('--refined',action='store_true')
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    card=json.loads((HERE/'protocol_b1_01.json').read_text())
    if args.N not in card['N'] or args.mass not in (0,.5): raise ValueError('Unregistered point')
    if args.output.exists() or args.output.with_suffix('.npz').exists():raise FileExistsError('Preserved result')
    data,checkpoint=run(card,args.N,args.mass,args.refined)
    np.savez_compressed(args.output.with_suffix('.npz'),**checkpoint)
    write(args.output,data)
    print(json.dumps({k:data[k] for k in ('N','mass_R_over_g','refined','seconds')}),flush=True)
