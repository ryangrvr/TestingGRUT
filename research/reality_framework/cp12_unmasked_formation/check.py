"""Independent full-basis/intertwiner checks for budget and probe solver."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,time
import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import expm_multiply
from run import model,error_evolve,graphstats
from frozen_model import build

P=Path(__file__).resolve().parent

def check(n,q,budget,source):
    sub,H,B,states,sizes,psi=model(n,q,budget,source)
    full,fs,lookup,_=build(n,q,2)
    rr=[];cc=[];vv=[];idx={s:i for i,s in enumerate(states)}
    for i,s in enumerate(fs):
        if sum(x*x for x in s[2])>budget:continue
        key,sign,size=sub.canonical(s)
        if sign and key in idx:rr.append(i);cc.append(idx[key]);vv.append(sign/np.sqrt(size))
    T=coo_matrix((vv,(rr,cc)),shape=(len(fs),len(states))).tocsr()
    assert np.linalg.norm((T.T@T-eye(H.shape[0])).data)<1e-10
    projected=T.T@full@T
    herr=float(np.linalg.norm((projected-H).data));assert herr<1e-10
    residual=full@T-T@H
    gram_error=float(np.linalg.norm((residual.T@residual-B.T@B).data));assert gram_error<1e-9
    evolved,err=error_evolve(H,B,psi,.2,128,q/2)
    actual=expm_multiply(-.2j*full,T@psi)
    diff=float(np.linalg.norm(actual-np.exp(-0j)*T@evolved*np.exp(-.2j*q/2)))
    # error_evolve returns evolution under H-shift, hence restore global phase.
    assert diff<=err+1e-11
    # Gauge-invariant local pulse: independent central finite difference.
    pulse=None
    if source:
        density=np.array([b&1 for g,b,E in states],float)
        detector=np.array([sum((b>>j)&1 for j in range(1,n))/(n-1) for g,b,E in states])
        wait=expm_multiply(-.1j*H,psi)
        Z=expm_multiply(-.1j*H,np.column_stack([wait,density*wait]))
        chi=float(2*np.vdot(Z[:,0],detector*Z[:,1]).imag)
        eta=1e-4
        W=expm_multiply(-.1j*H,np.column_stack([np.exp(-1j*eta*density)*wait,np.exp(1j*eta*density)*wait]))
        fd=float((detector@abs(W[:,0])**2-detector@abs(W[:,1])**2)/(2*eta))
        assert abs(chi-fd)<1e-9
        pulse=dict(analytic_response=chi,central_pulse_response=fd,absolute_error=abs(chi-fd))
    # Exact invariance of volume B_1=1+2Q/N for EVERY graph at fixed Q.
    for g in {s[0] for s in states}:
        stats=graphstats(n,g,sub.es)
        assert abs(stats[6]-(1+2*q/n))<1e-12
    return dict(N=n,Q=q,budget=budget,source_fixed=source,full_cutoff_dimension=len(fs),quotient_dimension=len(states),projected_H_error=herr,leakage_Gram_error=gram_error,direct_full_state_error=diff,analytic_error_upper=err,pulse_control=pulse)

if __name__=='__main__':
    assert 3*F(181,63000000)==F(2353,273000000)
    rows=[check(3,3,3,False),check(4,4,3,False),check(4,4,3,True)]
    out=dict(status='CP12_UNMASKED_SOLVER_CHECK_PASS',rows=rows,coefficient_reconciliation=True)
    (P/'CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
