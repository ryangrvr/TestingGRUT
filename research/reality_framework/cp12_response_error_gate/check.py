"""Independent dense and complete physical-matrix tests of the error bound."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix,coo_matrix
from scipy.sparse.linalg import expm_multiply
from control import response,error_evolve,parent_error_evolve
from run import model
from frozen_model import build
from probe import detectors

P=Path(__file__).resolve().parent

def dense_checks():
    rows=[]
    for seed in range(24):
        rng=np.random.default_rng(9010+seed);d=9;p=4
        a=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));H=(a+a.conj().T)/5
        HP=csr_matrix(H[:p,:p]);B=csr_matrix(H[p:,:p])
        psi=rng.normal(size=p)+1j*rng.normal(size=p);psi/=np.linalg.norm(psi)
        source=rng.uniform(size=d);O=rng.uniform(-1,1,size=d)
        # Norm-bounded nonunitary insertions, and detector signs of both kinds.
        w=.03+.01*(seed%5);lag=.04+.02*(seed%4)
        result,_,_=response(HP,B,psi,source[:p],[O[:p]],w,lag,128,0.)
        v=expm(-1j*w*H)@np.r_[psi,np.zeros(d-p)]
        x=expm(-1j*lag*H)@v;y=expm(-1j*lag*H)@(source*v)
        actual=float(2*np.vdot(x,O*y).imag);row=result[0]
        error=abs(actual-row['projected_response'])
        assert error<=row['observable_absolute_error_upper']+1e-12
        rows.append(dict(seed=seed,actual_error=error,observable_error_upper=row['observable_absolute_error_upper'],global_error_upper=row['global_absolute_error_upper']))
    return rows


def physical(n,q,budget):
    sub,H,B,states,sizes,psi=model(n,q,budget,True)
    full,fs,_,_=build(n,q,2);index={s:i for i,s in enumerate(states)}
    rr=[];cc=[];vv=[]
    for i,s in enumerate(fs):
        if sum(x*x for x in s[2])>budget:continue
        key,sign,size=sub.canonical(s)
        if sign and key in index:rr.append(i);cc.append(index[key]);vv.append(sign/np.sqrt(size))
    T=coo_matrix((vv,(rr,cc)),shape=(len(fs),len(states))).tocsr()
    pulse=np.array([b&1 for g,b,E in states],float)
    O=detectors(n,states,sub.es)
    result,_,_=response(H,B,psi,pulse,O,.1,.1,128,q/2)
    uf=expm_multiply(-.1j*full,T@psi)
    pulsef=np.array([b&1 for g,b,E in fs],float)
    Z=expm_multiply(-.1j*full,np.column_stack([uf,pulsef*uf]))
    Of=detectors(n,fs,sub.es)
    rows=[]
    for j,row in enumerate(result):
        actual=float(2*np.vdot(Z[:,0],Of[j]*Z[:,1]).imag)
        error=abs(actual-row['projected_response'])
        assert error<=row['observable_absolute_error_upper']+1e-12
        assert row['observable_absolute_error_upper']<row['global_absolute_error_upper']
        rows.append(dict(distance=j+1,direct_regulated_response=actual,actual_error=error,**row))
    return dict(N=n,Q=q,budget=budget,full_cutoff_dimension=len(fs),rows=rows)

if __name__=='__main__':
    # Multi-time exponential optimization must preserve the parent method.
    rng=np.random.default_rng(84);a=rng.normal(size=(8,8));h=csr_matrix((a+a.T)/4)
    leak=csr_matrix(rng.normal(size=(3,8)));v=rng.normal(size=(8,3))+1j*rng.normal(size=(8,3))
    vf,ef=error_evolve(h,leak,v,.2,32,0.)
    vs,es=parent_error_evolve(h,leak,v,.2,32,0.)
    state_delta=float(np.max(abs(vf-vs)));bound_delta=float(np.max(abs(ef-es)))
    assert state_delta<1e-10 and bound_delta<1e-10
    # Hostile case: an observable connecting retained and omitted spaces.
    h=np.array([[0.,1.],[1.,0.]])
    oy=np.array([[0.,-1j],[1j,0.]])
    u=expm(-.1j*h)@np.array([1.,0.])
    error=float(abs(np.vdot(u,oy@u).real))
    invalid_bound=.1**2
    assert error>invalid_bound
    hostile=dict(observable_commutes_with_projection=False,actual_expectation_error=error,invalid_block_diagonal_bound=invalid_bound,applicability_rejection_required=True)
    out=dict(status='CP12_OBSERVABLE_LEAKAGE_CONTROLS_PASS',dense=dense_checks(),physical=[physical(3,3,1),physical(4,4,1),physical(4,4,3)],hostile_control=hostile,time_grid_control=dict(maximum_state_difference=state_delta,maximum_bound_difference=bound_delta))
    (P/'CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
