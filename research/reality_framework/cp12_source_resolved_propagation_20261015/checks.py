"""Independent labelled N4 reference, signed lift and record-effect controls."""
from pathlib import Path
from math import sqrt,factorial,lcm
import json,pickle,time
import numpy as np
from shared import *
from response_run import response,ratio


def main():
    start=time.monotonic();n,q,budget,wait=4,5,256,5.
    full=physical_model(n,q,budget);H,B,states=(full[k] for k in ['H','B','states'])
    neutral=quotient_model(n,q,budget);sub=neutral['sub']
    w=np.load(cache_root()/'N4_Q5_B256_W5.npz')['w']
    ew=float(np.load(cache_root()/'N4_Q5_B256_W5.npz')['state_error'])
    by={s:i for i,s in enumerate(neutral['states'])};physical=np.zeros(len(states),complex)
    cores=np.zeros(len(states),np.int64);nv=np.zeros(len(states),int)
    sub.cache.clear()
    for j,s in enumerate(states):
        key,sign,size=Invariant.canonical(sub,s)
        if sign:physical[j]=w[by[key]]*sign/sqrt(size)
        cores[j]=flow(n,s)[0];nv[j]=(s[1]&15).bit_count()
    assert abs(np.linalg.norm(physical)-1)<1e-11
    Cs=cycle_masks(n);C=cycle_mask(n);which=Cs.index(C)
    success=np.isin(cores,Cs)
    records=np.column_stack([(cores==c)*physical for c in Cs]+[(~success)*physical])
    probabilities=np.linalg.norm(records,axis=0)**2;p=float(probabilities[:-1].sum())
    expected=json.loads((P/'NEUTRAL_N4_Q5_B256.json').read_text())['rows'][0]
    assert abs(p-expected['formation_probability'])<1e-12
    energy=float(np.vdot(records,H@records).real-np.vdot(physical,H@physical).real)
    assert abs(energy-expected['record_energy_change'])<1e-12
    entropy=float(-sum(a*np.log2(a) for a in probabilities if a>0))
    assert abs(entropy-expected['record_entropy_bits'])<1e-11
    crossings=0;checked=0
    for s,c in zip(states,cores):
        if c not in Cs:continue
        for t,a in transitions(n,s,120):
            tc=flow(n,t)[0]
            crossings+=tc in Cs and tc!=c;checked+=1
    assert crossings==0
    with (cache_root()/'N4_Q5_B256_W5_G200000_L0.pkl').open('rb') as f:seed=pickle.load(f)
    rooted=seed['sub'];root_index={s:i for i,s in enumerate(seed['states'])}
    embedding=np.zeros(len(states),complex)
    for j,s in enumerate(states):
        key,sign,size=Invariant.canonical(rooted,s)
        if sign and key in root_index:embedding[j]=seed['X'][root_index[key]]*sign/sqrt(size)
    vector=records[:,which]/sqrt(probabilities[which])
    lift_difference=float(np.linalg.norm(vector-embedding))
    assert lift_difference<3e-11
    bits=np.array([s[1] for s in states]);n0=(bits&1).astype(float)
    vectors=np.column_stack([records,physical]);width=vectors.shape[1]
    trajectory=MatrixFreePropagator(H,{'B':B},q/2,dt=1/64,batch=2).run(
        np.column_stack([vectors,n0[:,None]*vectors]),[.25,.5])
    reference=[];effects=[];d=distances(n,C)
    adaptive=json.loads((P/'RESPONSE_N4_Q5_B256_W5_G200000_L12.json').read_text())
    for lag,(z,errors) in trajectory.items():
        a=sqrt(probabilities[which]);rho=ew/sqrt(p)
        e0=sqrt(2)*rho/sqrt(1+sqrt(1-rho*rho))
        x,y=z[:,which]/a,z[:,which+width]/a
        ex=e0+float(errors['B'][which])/a;ey=e0+float(errors['B'][which+width])/a
        oldrow=next(r for r in adaptive['rows'] if r['lag']==lag);values={}
        S=np.array([ring_present(n,s,C) for s in states],float)
        for radius in [1,2]:
            O=np.mean(np.array([(bits>>j)&1 for j in range(n) if d[j]==radius]),axis=0)
            for joint in [False,True]:
                r=response(x,y,O*S if joint else O,ex,ey)
                name='density_and_same_core_response' if joint else 'density_response'
                old=oldrow['shells'][str(radius)][name]
                assert max(r['interval'][0],old['interval'][0])<=min(r['interval'][1],old['interval'][1])
                values[f'd{radius}'+('_joint_core' if joint else '')]=dict(
                    full_labelled=r['interval'],adaptive=old['interval'],
                    central_difference=abs(r['value']-old['value']))
        reference.append(dict(lag=lag,comparisons=values))
        O=((bits>>1)&1).astype(float)
        stored={}
        for name,ids in [('recorded',list(range(len(Cs)+1))),('unrecorded',[len(Cs)+1])]:
            xx,yy=z[:,ids],z[:,np.array(ids)+width]
            xe=ew+float(np.linalg.norm(errors['B'][ids]));ye=ew+float(np.linalg.norm(errors['B'][np.array(ids)+width]))
            stored[name]=response(xx,yy,O[:,None],xe,ye)
        a=stored['recorded'];b=stored['unrecorded']
        effects.append(dict(lag=lag,**stored,
            response_change_interval=[a['interval'][0]-b['interval'][1],a['interval'][1]-b['interval'][0]],
            response_change=a['value']-b['value']))
    # Direct integer-minor proof is inherited unchanged; the new cone bound
    # is checked on random complex vectors, allowing arbitrary physical phase.
    rng=np.random.default_rng(417);cone_checks=0
    for _ in range(1000):
        x=rng.normal(size=5)+1j*rng.normal(size=5);delta=(rng.normal(size=5)+1j*rng.normal(size=5))*.02
        y=x+delta;a=np.linalg.norm(x);rho=np.linalg.norm(delta)/a
        phase=np.vdot(x,y);distance=sqrt(max(0.,2-2*abs(phase)/(a*np.linalg.norm(y))))
        bound=sqrt(2)*rho/sqrt(1+sqrt(1-rho*rho))
        assert distance<=bound+1e-12;cone_checks+=1
    result=dict(status='PASS',labelled_physical_dimension=len(states),rooted_initial_lift_difference=lift_difference,
        raw_record_transition_checks=int(checked),cross_record_H_transitions=int(crossings),
        record_work_difference=abs(energy-expected['record_energy_change']),
        record_entropy_difference=abs(entropy-expected['record_entropy_bits']),
        complete_record_effect_on_fixed_receiver=effects,independent_response_controls=reference,
        conditional_ray_cone_checks=cone_checks,
        interpretation='Full labelled physical propagation independently checks the rooted source calculation. Nonselective complete-core-record response is compared with the original unrecorded state under the same full H; all record outcomes including failure are retained.',
        seconds=time.monotonic()-start)
    (P/'CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',lift_difference=lift_difference,
        transition_checks=checked,record_effects=effects,seconds=result['seconds'])),flush=True)


if __name__=='__main__':main()
