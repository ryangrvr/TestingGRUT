"""Auxiliary exact-arithmetic lemma diagnostics, never used by primary MPS."""
from pathlib import Path
from math import cos,sin,sqrt
import json
import numpy as np


def monotone(in_range):
    rng=np.random.default_rng(80311);n=4;size=2**n
    x=rng.normal(size=size)+1j*rng.normal(size=size);x/=np.linalg.norm(x)
    initial=x.copy();unitaries=[];projectors=[];errors=[];inside=[]
    nesting=commuting=0.;previous=np.eye(size)
    for j in range(n-1):
        local=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))[0]
        unitary=np.kron(np.eye(2**j),np.kron(local,np.eye(2**(n-j-2))))
        unitaries.append(unitary)
        for old in projectors:commuting=max(commuting,float(np.linalg.norm(unitary@old-old@unitary)))
        y=unitary@x;left=2**(j+1);right=2**(n-j-1)
        u,s,v=np.linalg.svd(y.reshape(left,right),full_matrices=False)
        prefix=u[:,:1];p=np.kron(prefix@prefix.conj().T,np.eye(right))
        nesting=max(nesting,float(np.linalg.norm(previous@p-p)))
        projected=p@y;errors.append(y-projected)
        f=(.001*(j+1) if in_range else 0)*projected
        inside.append(f);x=projected+f
        projectors.append(p);previous=p
    reference=initial.copy()
    for u in unitaries:reference=u@reference
    transported=[]
    for j,error in enumerate(errors):
        e=error.copy()
        for u in unitaries[j+1:]:e=u@e
        transported.append(e)
    overlap=np.array([[np.vdot(a,b) for b in transported] for a in transported])
    off=overlap-np.diag(np.diag(overlap))
    actual=float(np.linalg.norm(reference-x))
    orthogonal=sqrt(sum(float(np.vdot(e,e).real) for e in errors))
    fee=sum(float(np.linalg.norm(f)) for f in inside)
    assert nesting<1e-12 and commuting<1e-12 and np.max(abs(off))<1e-12
    assert actual<=orthogonal+fee+2e-13
    if not in_range:assert abs(actual-orthogonal)<2e-13
    return dict(dimension=size,nesting_error=nesting,later_gate_commutator=commuting,
        maximum_transported_error_overlap=float(np.max(abs(off))),
        actual_wavefunction_error=actual,orthogonal_discard_bound=orthogonal,
        in_range_error_sum=fee,sweep_bound=orthogonal+fee)


def counterexample():
    theta=.1;u=np.array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
    p=np.diag([1.,0.]);initial=np.array([1.,0.]);x=initial.copy();loss=[]
    for _ in range(2):
        y=u@x;loss.append(float(np.linalg.norm(y-p@y)));x=p@y
    actual=float(np.linalg.norm(u@u@initial-x));wrong=sqrt(sum(e*e for e in loss))
    assert actual>wrong and actual<=sum(loss)
    return dict(actual_error=actual,invalid_all_time_square_sum=wrong,
        valid_linear_bound=sum(loss),future_gate_commutes_with_projection=False)


def metric_control():
    rng=np.random.default_rng(81093)
    l=np.linalg.qr(rng.normal(size=(5,3))+1j*rng.normal(size=(5,3)))[0]
    l=l@(np.eye(3)+1e-5*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))))
    t=np.kron(l,np.eye(2));h=t.conj().T@t
    k=np.linalg.qr(rng.normal(size=(6,2))+1j*rng.normal(size=(6,2)))[0]
    hn=k.conj().T@h@k
    x=rng.normal(size=(6,4))+1j*rng.normal(size=(6,4))
    y=np.linalg.solve(hn,k.conj().T@h@x)
    z=y+1e-6*(rng.normal(size=y.shape)+1j*rng.normal(size=y.shape))
    w=hn@z-k.conj().T@h@x
    actual=float(np.linalg.norm(t@k@(z-y)))
    eigen=np.linalg.eigvalsh(hn);bound=float(sqrt(eigen[-1])/eigen[0]*np.linalg.norm(w))
    p=t@k@np.linalg.solve(hn,k.conj().T@t.conj().T)
    pe=float(np.linalg.norm(p@p-p));he=float(np.linalg.norm(p-p.conj().T))
    assert pe<1e-12 and he<1e-12 and actual<=bound+1e-14
    return dict(prefix_Gram_error=float(np.linalg.norm(h-np.eye(6),2)),
        retained_Gram_minimum=float(eigen[0]),retained_Gram_maximum=float(eigen[-1]),
        in_range_defect=actual,direct_metric_residual_bound=bound,
        projector_idempotence_error=pe,projector_Hermiticity_error=he)


def run():
    return dict(status='PASS',used_in_S8_or_G1=False,new_primary_run=False,
        arithmetic='ordinary double diagnostics of a separate analytical lemma',
        exact_projection=monotone(False),in_range_error=monotone(True),
        coherent_return_counterexample=counterexample(),nonorthogonal_prefix_metric=metric_control())


if __name__=='__main__':
    result=run();Path('SWEEP_CONTROLS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
