"""Supplied common phase with small fluctuations; no phase-selector added."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import numpy as np
from model import minimum,initial_z,YQ,yield_factor
from run import evolve,write,ROOT
from scan import worker,compact
from build import WORK


def ensemble(H,cH,delta,sigma,N=64,seed=20261025,workers=8,rtol=2e-12,kind='source_phase'):
    rng=np.random.default_rng(seed)
    if kind=='source_phase':
        seeds=rng.normal(0,sigma,N)
        jobs=[dict(H=H,cH=cH,delta=delta+float(d),rtol=rtol,atol=2e-15) for d in seeds]
    else:
        seeds=rng.normal(0,sigma,(N,3))
        jobs=[]
        base=initial_z(minimum(H,cH),delta)
        for d in seeds:
            z=base*np.exp(1j*d)
            jobs.append(dict(H=H,cH=cH,delta=delta,rtol=rtol,atol=2e-15,z0=z.tolist()))
    # Complex initial values cannot be directly JSON cache keys.
    if kind=='source_phase':
        with ProcessPoolExecutor(max_workers=workers) as pool:r=list(pool.map(worker,jobs))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:r=list(pool.map(_direct,jobs))
    rr=[compact(x) for x in r]
    nx=np.array([x['NX'] for x in rr]);nm=np.array([x['NcircX'] for x in rr]);eps=nx/nm
    mean=nx.mean()/nm.mean()
    # Delta-method SEM for ratio of means, preserving correlated denominator.
    sem=np.std(nx-mean*nm,ddof=1)/(np.sqrt(N)*nm.mean())
    qmean=nx.mean()*yield_factor(H)/YQ
    qsem=nx.std(ddof=1)/np.sqrt(N)*yield_factor(H)/YQ
    return dict(H_GeV=H,cH=cH,delta=delta,seed_sigma_rad=sigma,N=N,rng_seed=seed,seed_kind=kind,rtol=rtol,
                eps_net=float(mean),eps_net_sem=float(sem),local_eps_rms=float(np.sqrt(np.mean(eps**2))),
                local_eps_mean=float(eps.mean()),local_eps_std=float(eps.std(ddof=1)),
                fraction_positive=float(np.mean(nx>0)),YX_over_YQ_mean=float(qmean),YX_over_YQ_sem=float(qsem),
                n_greater_than_019=int(np.sum(abs(eps)>=.19)),n_greater_than_099=int(np.sum(abs(eps)>=.99)),
                rows=rr)


def _direct(job):
    return evolve(**job)


def main():
    H=8718444157511.949
    configs=[(4.7e13,1.,0.,1e-3,64,'three_phases'),
             (H,1.,np.pi/6,1e-3,64,'source_phase'),
             (H,1.,7*np.pi/12,1e-3,64,'source_phase'),
             (H,1.,np.pi/6,.03,64,'source_phase'),
             (H,1.,7*np.pi/12,.03,64,'source_phase')]
    results=[]
    for i,(h,c,d,s,n,k) in enumerate(configs):
        r=ensemble(h,c,d,s,n,20261025+i,kind=k);results.append(r)
        print('ENSEMBLE',h,c,d,s,k,'eps',r['eps_net'],'sem',r['eps_net_sem'],'rms',r['local_eps_rms'],flush=True)
        write(WORK/'initial_ensembles_raw.json',dict(scope='Homogeneous prescribed phase-noise ensembles; not spatial-volume weights or the inflationary stochastic solution.',rows=results))


if __name__=='__main__':main()
