"""Algebra and independent small dense checks; no registered production point."""
import itertools,json
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from model import HERE,algebra_checks,reconstruction,tv
from continuum import verify_reference
from lattice import schwinger

def precheck(card):
    checks=algebra_checks(card)
    toy=dict(card,oscillator_verification={'mean':.17,'h':.21,'omega':.83,
        'times':[0,.23,.91],'phases':[.13,-.21,.34],'cutoffs':[16,24]})
    oscillator=verify_reference(toy)
    checks['oscillator_vs_ordered_gaussian']=oscillator['tv_gaussian_vs_high']
    checks['oscillator_cutoff_step']=oscillator['tv_cutoff_change']
    n=6;a=card['length']/(n+1);mass=-a/8
    sparse,electric,states=schwinger.hamiltonian(n,a**-2,mass);h=sparse*a/2
    # Independent element assembly from spin configurations and Gauss law.
    dense=np.zeros(h.shape);lookup={int(s):i for i,s in enumerate(states)}
    for col,state in enumerate(states):
        spin=np.array([1 if int(state)&(1<<j) else -1 for j in range(n)])
        flux=0.;energy=0.
        for j in range(n):
            energy+=mass/2*(1+(-1)**j*spin[j])
            flux+=(spin[j]+(-1)**j)/2
            if j<n-1:
                energy+=a*flux*flux/2
                if spin[j]!=spin[j+1]:dense[lookup[int(state)^((1<<j)|(1<<(j+1)))],col]=1/(2*a)
        dense[col,col]+=energy
    checks['independent_small_hamiltonian']=float(np.max(abs(h.toarray()-dense)))
    energies,vectors=np.linalg.eigh(dense);ground=vectors[:,0].astype(complex)
    _,_,weights,boundary,_=reconstruction(card,n)
    field=electric@weights+(boundary.sum()-weights@((-1.)**np.arange(n-1)))/4
    W=np.diag(np.exp(1j*card['coupling']*field));times=[0,.23,.91]
    direct={}
    for signs in itertools.product((1,-1),repeat=3):
        v=ground.copy();previous=0.
        for t,s in zip(times,signs):
            v=expm(-1j*dense*(t-previous))@v;v=(v+s*(W@v))/2;previous=t
        direct[''.join('+' if s==1 else '-' for s in signs)]=float(np.vdot(v,v).real)
    branches=ground[:,None];labels=[''];previous=0.
    for t in times:
        if t>previous:branches=expm_multiply(-1j*h*(t-previous),branches,traceA=-1j*(t-previous)*h.diagonal().sum())
        shifted=np.diag(W)[:,None]*branches
        branches=np.concatenate(((branches+shifted)/2,(branches-shifted)/2),axis=1)
        labels=[s+'+' for s in labels]+[s+'-' for s in labels];previous=t
    batched={s:float(np.vdot(branches[:,j],branches[:,j]).real) for j,s in enumerate(labels)}
    checks['dense_vs_sparse_complete_history']=tv(direct,batched)
    checks['small_instrument_completeness']=float(np.max(abs(((np.eye(len(W))+W)/2).conj().T@((np.eye(len(W))+W)/2)+
      ((np.eye(len(W))-W)/2).conj().T@((np.eye(len(W))-W)/2)-np.eye(len(W)))))
    # A deliberately wrong detector normalization is rejected by this identity.
    bad=(np.eye(len(W))+W)/np.sqrt(2)
    wrong=float(np.linalg.norm(bad.conj().T@bad+((np.eye(len(W))-W)/np.sqrt(2)).conj().T@((np.eye(len(W))-W)/np.sqrt(2))-np.eye(len(W))))
    return {'checks':checks,'negative_control_wrong_normalization_rejected':wrong>1,
      'passed':bool(max(checks.values())<card['gates']['algebra'] and wrong>1),
      'scope':'Synthetic oscillator and N6 at times 0,.23,.91 only; no B2 production probability examined.'}

if __name__=='__main__':
    print(json.dumps(precheck(json.loads((HERE/'protocol_b2_01.json').read_text())),indent=2))
