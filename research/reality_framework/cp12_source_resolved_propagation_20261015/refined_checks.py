"""Full labelled reference for the refined instrument and endpoint source."""
from math import sqrt
import json,pickle,time
import numpy as np
from shared import *
from response_run import response


def main():
    start=time.monotonic();n,q,budget=4,5,256
    m=quotient_model(n,q,budget);full=physical_model(n,q,budget)
    z=np.load(cache_root()/'N4_Q5_B256_W5.npz');w=z['w'];ew=float(z['state_error'])
    H,B,states=(full[k] for k in ['H','B','states']);by={s:i for i,s in enumerate(m['states'])}
    Cs=cycle_masks(n);labels=sorted({(c,s[0]) for s in states for c in [flow(n,s)[0]] if c in Cs})
    assert len(labels)==6
    physical=np.zeros(len(states),complex);assign=np.full(len(states),6,int)
    for j,s in enumerate(states):
        key,sign,size=Invariant.canonical(m['sub'],s)
        if sign:physical[j]=sign*w[by[key]]/sqrt(size)
        label=(flow(n,s)[0],s[0])
        if label in labels:assign[j]=labels.index(label)
    vectors=np.column_stack([(assign==k)*physical for k in range(7)]+[physical])
    probabilities=np.linalg.norm(vectors[:,:7],axis=0)**2
    work=float(np.vdot(vectors[:,:7],H@vectors[:,:7]).real-np.vdot(physical,H@physical).real)
    entropy=float(-sum(p*np.log2(p) for p in probabilities if p>0))
    published=json.loads((P/'REFINED_RECORD_N4_Q5_B256_W5.json').read_text())
    assert abs(work-published['complete_record_work'])<1e-11
    assert abs(entropy-published['complete_record_entropy_bits'])<1e-11
    C=cycle_mask(n);host=C|(1<<1);k=labels.index((C,host));a=sqrt(probabilities[k])
    with (cache_root()/'N4_Q5_B256_W5_G200000_F2_L0.pkl').open('rb') as f:seed=pickle.load(f)
    ids={s:i for i,s in enumerate(seed['states'])};embedding=np.zeros(len(states),complex)
    for j,s in enumerate(states):
        key,sign,size=Invariant.canonical(seed['sub'],s)
        if sign and key in ids:embedding[j]=seed['X'][ids[key]]*sign/sqrt(size)
    lift_error=float(np.linalg.norm(vectors[:,k]/a-embedding));assert lift_error<3e-11
    n0=np.array([s[1]&1 for s in states],float);width=vectors.shape[1]
    trajectory=Propagator(H,{'B':B},q/2,dt=1/64,batch=8).run(
        np.column_stack([vectors,n0[:,None]*vectors]),[.25,.5])
    adaptive=json.loads((P/'RESPONSE_N4_Q5_B256_W5_G200000_F2_L12.json').read_text())
    checks=[];effects=[]
    for lag,(z,er) in trajectory.items():
        rho=ew/sqrt(probabilities[:6].sum());e0=sqrt(2)*rho/sqrt(1+sqrt(1-rho*rho))
        x,y=z[:,k]/a,z[:,k+width]/a
        ex=e0+float(er['B'][k])/a;ey=e0+float(er['B'][k+width])/a
        old=next(r for r in adaptive['rows'] if r['lag']==lag)
        for receiver in [1,2]:
            O=np.array([(s[1]>>receiver)&1 for s in states],float)
            ref=response(x,y,O,ex,ey);other=old['fixed_receiver_pairs'][str(receiver)]['density_response']
            assert max(ref['interval'][0],other['interval'][0])<=min(ref['interval'][1],other['interval'][1])
            checks.append(dict(lag=lag,receiver=receiver,labelled_interval=ref['interval'],
                rooted_interval=other['interval'],central_difference=abs(ref['value']-other['value'])))
        O=np.array([(s[1]>>1)&1 for s in states],float);effect={}
        for name,columns in [('recorded',list(range(7))),('unrecorded',[7])]:
            yy=np.array(columns)+width
            effect[name]=response(z[:,columns],z[:,yy],O[:,None],
                ew+float(np.linalg.norm(er['B'][columns])),ew+float(np.linalg.norm(er['B'][yy])))
        effect['lag']=lag;aa,bb=effect['recorded'],effect['unrecorded']
        effect['response_change_interval']=[aa['interval'][0]-bb['interval'][1],aa['interval'][1]-bb['interval'][0]]
        effects.append(effect)
    result=dict(status='PASS',physical_dimension=len(states),complete_record_outcomes=7,
        record_work=work,work_identity_difference=abs(work-published['complete_record_work']),
        record_entropy_bits=entropy,entropy_identity_difference=abs(entropy-published['complete_record_entropy_bits']),
        endpoint_lift_difference=lift_error,labelled_response_checks=checks,
        complete_record_effect_on_fixed_receiver=effects,
        interpretation='All refined record outcomes including failure are retained for the nonselective comparison; no autonomous record formation or zero-cost measurement is assumed.',
        seconds=time.monotonic()-start)
    (P/'REFINED_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
