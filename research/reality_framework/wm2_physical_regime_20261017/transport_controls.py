"""Reference-walk and exact support controls of the analytic transport lemma.

These are not new formation trajectories, a changed physical model, or a
numerical replacement for the full-H tail comparison.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
import time
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply


def primitive_words(N,m,ell):
    initial=tuple(range(m));prefix=initial[:ell];q=ell-2
    current={(initial,0),(initial,1)};visited=set(current);first=None;counts=[]
    def successors(w,p,reference):
        out=set()
        for v in range(N):
            if v not in w:
                out.add(('assist+',w+(v,),p))
                if not reference and p in [w[-1],v]:out.add(('bare+',w+(v,),p))
        if len(w)>1:
            if not reference or len(w)>ell:out.add(('assist-',w[:-1],p))
            if not reference and p in w[-2:]:out.add(('bare-',w[:-1],p))
        if reference:
            x=prefix.index(p)
            if x:out.add(('hop',w,prefix[x-1]))
            if x+1<ell:out.add(('hop',w,prefix[x+1]))
            else:out.add(('dummy',w,ell))
        elif p in w:
            x=w.index(p)
            for y in [x-1,x+1]:
                if 0<=y<len(w):out.add(('hop',w,w[y]))
        return out
    for depth in range(q+1):
        future=set();counts.append(len(current))
        for w,p in current:
            a=successors(w,p,False);b=successors(w,p,True)
            if a!=b:
                first=depth+1 if first is None else min(first,depth+1)
                if depth<q:raise AssertionError('short physical/reference word differs')
            if depth<q:
                for kind,ww,pp in a:
                    assert not kind.startswith('bare') and pp in prefix
                    assert ww[:ell]==prefix
                    if (ww,pp) not in visited:future.add((ww,pp))
        visited.update(future);current=future
    return dict(N=N,initial_length=m,prefix_length=ell,q=q,states_checked=len(visited),
        layer_sizes=counts,first_difference_order=first,
        scope='exact support overapproximation; contact occupations not suppressed')


def reference_walk():
    J=.25;L=2048;x=np.arange(L,dtype=float)
    H=diags([-J*np.ones(L-1),-J*np.ones(L-1)],[-1,1],format='csr')
    X=diags(x,format='csr');V=1j*(H@X-X@H)
    acceleration=1j*(H@V-V@H)
    boundary=diags(np.r_[2*J*J,np.zeros(L-2),-2*J*J],format='csr')
    difference=acceleration-boundary
    assert not difference.nnz or max(abs(difference.data))==0
    initial=np.zeros(L,complex);initial[:2]=1/np.sqrt(2)
    rows=[]
    for t in [4.,16.,64.,256.,1024.]:
        psi=expm_multiply(-1j*t*H,initial,traceA=0.)
        p=abs(psi)**2;mean=float(np.dot(x,p));second=float(np.dot(x*x,p))
        threshold=int(np.ceil(J*t/128));prob=float(sum(p[threshold:]))
        assert mean>=J*t/64 and second<=9*J*J*t*t+1e-9
        assert prob>=1/147456
        rows.append(dict(t=t,comparison_sites=L,mean=mean,second_moment=second,
            distance_threshold=threshold,arrival_probability=prob,
            upper_end_probability=float(p[-1]),physical_formation_run=False))
    assert Q(7,2**24)<Q(1,294912)
    return dict(finite_matrix_acceleration_identity_error=0.,walks=rows,
        conditional_arrival_lower='1/294912',original_trial_joint_lower='1/589824')


def run():
    start=time.monotonic()
    return dict(status='PASS',primitive_support=[primitive_words(N,N//2,N//4) for N in [12,16,20]],
        reference=reference_walk(),elapsed_seconds=time.monotonic()-start,
        scope='analytic-comparison controls, not a full probe simulation or external validation')


if __name__=='__main__':
    r=run();Path(__file__).with_name('TRANSPORT_CONTROLS.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
