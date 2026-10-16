"""Full labelled N4 check of host-only lift, source and record disturbance."""
from math import sqrt
import json,pickle,time
import numpy as np
from common import *


def main():
    start=time.monotonic();n,q,b=4,5,256;m=quotient_model(n,q,b);full=physical_model(n,q,b)
    states,H,B=(full[k] for k in ['states','H','B']);by={s:i for i,s in enumerate(m['states'])}
    zz=np.load(neutral_path(n,q,b,5));w=zz['w'];ew=float(zz['state_error'])
    physical=np.zeros(len(states),complex)
    for j,s in enumerate(states):
        key,sign,size=Invariant.canonical(m['sub'],s)
        if sign:physical[j]=w[by[key]]*sign/sqrt(size)
    hosts=sorted({s[0] for s in states});vectors=np.column_stack([physical*np.array([s[0]==g for s in states]) for g in hosts]+[physical])
    probabilities=np.linalg.norm(vectors[:,:len(hosts)],axis=0)**2
    work=float(np.vdot(vectors[:,:len(hosts)],H@vectors[:,:len(hosts)]).real-np.vdot(physical,H@physical).real)
    entropy=float(-sum(p*np.log2(p) for p in probabilities if p>0))
    neutral=json.loads((P/'HOST_N4_Q5_B256.json').read_text());last=neutral['rows'][-1]
    assert abs(work-last['complete_record_work'])<1e-11 and abs(entropy-last['complete_record_entropy_bits'])<1e-11
    with (cache_root()/'HOST_N4_Q5_B256_W5_G100000_L0.pkl').open('rb') as f:seed=pickle.load(f)
    root=seed['sub'];lookup={s:i for i,s in enumerate(seed['states'])};embedding=np.zeros(len(states),complex)
    host=seed['initial']['host'];k=hosts.index(host);a=sqrt(probabilities[k])
    for j,s in enumerate(states):
        key,sign,size=Invariant.canonical(root,s)
        if sign and key in lookup:embedding[j]=seed['X'][lookup[key]]*sign/sqrt(size)
    lift=float(np.linalg.norm(embedding-vectors[:,k]/a));assert lift<2e-10
    bits=np.array([s[1] for s in states]);n0=(bits&1).astype(float);width=vectors.shape[1]
    trajectory=Propagator(H,{'electric':B},q/2,dt=1/64,batch=4).run(np.column_stack([vectors,n0[:,None]*vectors]),[.25,.5])
    adaptive=json.loads((P/'HOST_RESPONSE_N4_Q5_B256_L6.json').read_text());checks=[];effects=[]
    for t,(z,er) in trajectory.items():
        old=next(r for r in adaptive['rows'] if r['lag']==t);rho=ew/sqrt(6)/a
        e0=sqrt(2)*rho/sqrt(1+sqrt(1-rho*rho))
        x,y=z[:,k]/a,z[:,k+width]/a;ex=e0+float(er['electric'][k])/a;ey=e0+float(er['electric'][k+width])/a
        for receiver in [1,3]:
            O=((bits>>receiver)&1).astype(float);r=response(x,y,O,ex,ey);s=old['receivers'][str(receiver)]['response']
            assert max(r['interval'][0],s['interval'][0])<=min(r['interval'][1],s['interval'][1])
            checks.append(dict(lag=t,receiver=receiver,labelled_interval=r['interval'],rooted_interval=s['interval'],central_difference=abs(r['value']-s['value'])))
        O=((bits>>1)&1).astype(float);effect={}
        for name,ids in [('recorded',list(range(len(hosts)))),('unrecorded',[len(hosts)])]:
            sy=np.array(ids)+width
            effect[name]=response(z[:,ids],z[:,sy],O[:,None],ew+float(np.linalg.norm(er['electric'][ids])),ew+float(np.linalg.norm(er['electric'][sy])))
        aa,bb=effect['recorded'],effect['unrecorded'];effect.update(lag=t,response_change_interval=[aa['interval'][0]-bb['interval'][1],aa['interval'][1]-bb['interval'][0]])
        effects.append(effect)
    result=dict(status='PASS',labelled_dimension=len(states),host_only_initial_lift_difference=lift,
        record_work_identity_difference=abs(work-last['complete_record_work']),
        record_entropy_identity_difference=abs(entropy-last['complete_record_entropy_bits']),
        fixed_receiver_checks=checks,complete_record_effects=effects,seconds=time.monotonic()-start)
    (P/'RESPONSE_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
