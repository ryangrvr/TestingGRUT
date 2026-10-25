"""Test every refinement candidate with phase seeds; scrutinize raw eps>=.19.

The ensemble is an explicit sensitivity experiment, not a probability law
derived from inflation. Passing its threshold would remain conditional.
"""
import json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import minimize
from run import ROOT,evolve,write
from model import minimum,parameters,potential,force
from ensemble import ensemble
from scan import compact,worker
from build import WORK


def conditional_preparation(H,cH,delta):
    point=minimum(H,cH);p=parameters(H)
    phase=np.array([1.,-1.,-np.exp(1j*delta)])
    def v(r):return potential(r*phase,p,cH)
    def g(r):return 2*np.real(np.conj(phase)*force(r*phase,p,cH))
    best=None
    for scale in [1.,.3,3.]:
        r=minimize(v,np.array(point['amplitudes_over_H'])*scale,jac=g,
                   bounds=[(0.,None)]*3,method='L-BFGS-B',options=dict(ftol=1e-15,gtol=1e-10,maxiter=10000))
        if best is None or r.fun<best.fun:best=r
    return best.x*phase,dict(amplitudes_over_H=best.x.tolist(),
        interpretation='Radial minimum conditional on externally held phase, including possible amplitude boundary. Not an unforced inflationary minimum.')


def main():
    candidates=json.loads((ROOT/'REFINED_CANDIDATES.json').read_text())['rows']
    out=dict(scope='Conditional homogeneous seed robustness. Not an inflationary distribution, spatial domain calculation, or final relic census.',
             candidate_ensembles=[],high_fraction_controls=[])
    for i,r in enumerate(candidates):
        H,cH,delta=r['H_GeV'],r['cH'],r['delta']
        high=abs(r['epsX'])>=.19
        n=128 if high else 48
        e=ensemble(H,cH,delta,1e-3,N=n,seed=2026102500+i,rtol=2e-12)
        e['nominal_epsX']=r['epsX'];e['candidate_index']=i
        out['candidate_ensembles'].append(e)
        print('CANDIDATE',i,'/',len(candidates),'cH',cH,'H',f'{H:.6g}','raw',r['epsX'],
              'net',e['eps_net'],'sem',e['eps_net_sem'],flush=True)
        write(WORK/'candidate_validation_raw.json',out)
        if high:
            larger=ensemble(H,cH,delta,.03,N=128,seed=2026102600+i,rtol=2e-12)
            jobs=[dict(H=H,cH=cH,delta=delta,rtol=t,atol=t/100,ns=120) for t in [2e-12,2e-14]]
            jobs.append(dict(H=H,cH=cH,delta=delta+1e-9,rtol=2e-14,atol=2e-16,ns=120))
            jobs.append(dict(H=H,cH=cH,delta=delta,rtol=2e-12,atol=2e-14,end=3000,ns=240))
            with ProcessPoolExecutor(max_workers=4) as pool:controls=list(pool.map(worker,jobs))
            z,prep=conditional_preparation(H,cH,delta)
            radial=evolve(H,cH,delta,z0=z,rtol=2e-12,atol=2e-14)
            ctrl=dict(H_GeV=H,cH=cH,delta=delta,candidate_index=i,larger_seed_ensemble=larger,
                      tighter_and_longer_controls=controls,conditional_preparation=prep,
                      conditional_radial_release=radial)
            out['high_fraction_controls'].append(ctrl)
            write(WORK/'candidate_validation_raw.json',out)
    print('CANDIDATE_VALIDATION_COMPLETE',len(candidates),flush=True)


if __name__=='__main__':main()
