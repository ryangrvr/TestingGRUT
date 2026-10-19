#!/usr/bin/env python3
"""Independent evolution and predeclared attacks on the actual candidate."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from calculate import BASE, product, state, hamiltonian, kron_hamiltonian, evolve

HERE = Path(__file__).resolve().parent


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def bit_success(E,V,z,L,kind):
    ts = np.linspace(0,40,201)
    means = []
    for r in [1,-1]:
        if kind == 'local':
            psi = state(L,'retention',1)
            if r == -1:
                psi = psi[np.arange(len(psi)) ^ (1 << L)]
        elif kind == 'write':
            plusx = np.array([1.,1.])/np.sqrt(2)
            matter = [plusx]*L
            matter[0] = [1.,0.] if r == 1 else [0.,1.]
            psi = product(matter + [plusx]*L)
        else:
            psi = state(L,'retention',r)
        evolved = evolve(E,V,psi,ts)
        means.append(np.sum(np.abs(evolved)**2*z[L,:,None],axis=0))
    success = .5+np.abs(means[0]-means[1])/4
    answer = dict(times=ts.tolist(),means=[x.tolist() for x in means],success=success.tolist(),
                  minimum_t10=float(np.min(success[ts<=10])),p10=float(success[50]),p40=float(success[-1]),
                  maximum=float(np.max(success)),maximum_time=float(ts[np.argmax(success)]))
    if kind == 'write':
        witnesses=[]
        for j,t in enumerate(ts):
            if t<=10 and success[j]>=.8 and np.min(success[j:j+51])>=.8:
                witnesses.append(float(t))
        answer['write_and_hold_pass']=bool(witnesses)
        answer['witness_times']=witnesses
    return answer


def main():
    output = HERE/'challenge_results.json'
    if output.exists():
        raise SystemExit('Refusing to overwrite challenge record')
    start=time.monotonic()
    initial=json.loads((HERE/'results_initial.json').read_text())
    result=dict(initial_results_sha256=digest(HERE/'results_initial.json'),
                challenge_protocol_sha256=digest(HERE/'challenge_protocol.json'),
                challenge_source_sha256=digest(Path(__file__)),checks={})
    L=4
    H=kron_hamiltonian(L,BASE)
    direct,z=hamiltonian(L,BASE)
    result['independent_full_H_max_error']=float(np.max(np.abs(H-direct)))
    E,V=eigh(H)
    sparse=csr_matrix(-1j*H)
    independent=[]
    for r in [1,-1]:
        psi=state(L,'retention',r)
        trajectory=expm_multiply(sparse,psi,start=0,stop=40,num=201,endpoint=True)
        m=np.sum(np.abs(trajectory)**2*z[L,None,:],axis=1)
        saved=np.array(initial['runs']['candidate_L4']['retention']['means'][0 if r==1 else 1])
        independent.append(dict(bit=r,max_mean_difference=float(np.max(np.abs(m-saved)))))
    result['krylov_retention']=independent
    t=initial['runs']['candidate_L4']['readout']['best_registered_time']
    measurement=[]
    ids=np.arange(len(H))
    Y=np.zeros_like(H,dtype=complex)
    Y[ids,ids^2]=-1j*z[1]
    for r in [1,-1]:
        psi=expm_multiply(t*sparse,state(L,'readout',r))
        energy=float(np.vdot(psi,H@psi).real)
        ypsi=Y@psi
        energy_after=.5*(energy+float(np.vdot(ypsi,H@ypsi).real))
        meanY=float(np.vdot(psi,ypsi).real)
        original=initial['runs']['candidate_L4']['readout']['means'][0 if r==1 else 1][int(round(t/.02))]
        measurement.append(dict(bit=r,time=t,energy_before=energy,
                                energy_after_nonselective_Y=energy_after,
                                mean_energy_change=energy_after-energy,mean_Y=meanY,
                                spectral_mean_difference=abs(meanY-original)))
    result['measurement']=measurement
    result['local_bit']={}
    result['autonomous_write']={}
    for L in [4,5]:
        h,zl=hamiltonian(L,BASE)
        en,vec=eigh(h)
        result['local_bit'][str(L)]=bit_success(en,vec,zl,L,'local')
        result['autonomous_write'][str(L)]=bit_success(en,vec,zl,L,'write')
    result['unfacilitated_perturbations']=[]
    for kappa in [.1,.3]:
        pert=H.copy()
        for i in range(4):
            pert[ids^(1<<(4+i)),ids]-=kappa
        en,vec=eigh(pert)
        attacked=bit_success(en,vec,z,4,'repetition')
        attacked['kappa']=kappa
        result['unfacilitated_perturbations'].append(attacked)
    p=np.exp(-(E-E[0]));p/=p.sum();y=p**.25;sp=np.sqrt(p)
    A=V.T@V[ids^1]
    B=V.T@V[ids^2]
    ma=float(p@np.diag(A));mb=float(p@np.diag(B))
    A-=ma*np.eye(len(E));B-=mb*np.eye(len(E))
    fd=float(np.sum(sp[:,None]*sp[None,:]*A*A)*np.sum(sp[:,None]*sp[None,:]*B*B))
    F=[];ts=np.linspace(0,4,41)
    for tt in ts:
        at=A*np.exp(1j*(E[:,None]-E[None,:])*tt)
        d=(y[:,None]*at*y[None,:])@B
        F.append(float(np.einsum('nm,mn->',d,d).real))
    result['centered_regulated_otoc']=dict(times=ts.tolist(),regulated_F=F,
            removed_mean_A=ma,removed_mean_B=mb,factorized_reference=fd,
            initial_relative_factorization_defect=abs(F[0]-fd)/abs(fd),
            lyapunov_exponent=None)
    result['checks']={
        'independent_full_H_match':result['independent_full_H_max_error']<1e-13,
        'krylov_retention_matches_spectral':max(x['max_mean_difference'] for x in independent)<1e-9,
        'krylov_readout_matches_spectral':max(x['spectral_mean_difference'] for x in measurement)<1e-9,
        'local_independent_bit_retains_to_t10_L4':result['local_bit']['4']['minimum_t10']>=.9,
        'local_independent_bit_retains_to_t10_L5':result['local_bit']['5']['minimum_t10']>=.9,
        'autonomous_write_and_hold_L4':result['autonomous_write']['4']['write_and_hold_pass'],
        'autonomous_write_and_hold_L5':result['autonomous_write']['5']['write_and_hold_pass'],
        'perturbed_kappa0_1_retains_to_t10':result['unfacilitated_perturbations'][0]['minimum_t10']>=.9,
        'perturbed_kappa0_3_retains_to_t10':result['unfacilitated_perturbations'][1]['minimum_t10']>=.9
    }
    result['wall_seconds']=time.monotonic()-start
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=result['checks'],local_p10={k:v['p10'] for k,v in result['local_bit'].items()},
                         writes={k:v['maximum'] for k,v in result['autonomous_write'].items()},
                         measurements=measurement,seconds=result['wall_seconds'])),flush=True)


if __name__=='__main__':
    main()
