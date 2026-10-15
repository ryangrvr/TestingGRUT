"""Full-law density response at fixed recorded-ring distances (auxiliary N4).

This is a neutral-origin recorded probe, not a Q=2N scaling experiment.
Both neighbours are averaged, so no arbitrary neighbour label is preferred.
No term is removed from the evolution Hamiltonian.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from boundary import RF
from operators import physical_model
from integrals import Propagator
from recorded_transport import ES,ring_for,ring_present,selected,opposite,response,ratio
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
from stream import probability_interval

P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--budget',type=int,default=256)
    ap.add_argument('--wait',type=float,default=5.);ap.add_argument('--lags',default='.125,.25,.5')
    ap.add_argument('--out',default=str(P/'N4_CORE_RESPONSE.json'));args=ap.parse_args()
    start=time.monotonic();m=physical_model(4,5,args.budget)
    H,B,states,psi=(m[k] for k in ['H','B','states','psi'])
    w,er=Propagator(H,{'B':B},2.5).run(psi,[args.wait])[args.wait];ew=float(er['B'][0])
    hosts=np.array([s[0] for s in states]);gs=sorted(set(hosts));G=np.array([selected(s) for s in states])
    occ=np.array([[(s[1]>>j)&1 for j in range(4)] for s in states],float)
    X=np.column_stack([(hosts==g)*G*w for g in gs]);Y=occ[:,0,None]*X
    records=np.column_stack([(hosts==g)*(G if flag else ~G)*w for g in gs for flag in [0,1]])
    probs=np.linalg.norm(records,axis=0)**2
    record_energy=float(np.vdot(records,H@records).real-np.vdot(w,H@w).real)
    workerr=2*(np.sqrt(40)+.1*5/6)*(2*ew+ew*ew)+5e-10
    entropy=float(-sum(p*np.log2(p) for p in probs if p>0))
    tv=min(1.,2*ew+ew*ew)
    h2=lambda x:0. if x<=0 or x>=1 else -x*np.log2(x)-(1-x)*np.log2(1-x)
    entropyerr=min(np.log2(12),h2(tv)+tv*np.log2(11)) if tv<=11/12 else np.log2(12)
    labels=[];neighbors=[];far=[];survival=[]
    for g in gs:
        C,chord=ring_for(int(g));j=opposite(C)
        ns=sorted(b if a==0 else a for e,(a,b) in enumerate(ES) if C>>e&1 and 0 in (a,b))
        assert len(ns)==2
        neighbors.append(occ[:,ns].mean(axis=1));far.append(occ[:,j])
        survival.append(np.array([ring_present(s,C) for s in states],float))
        labels.append(dict(graph_mask=int(g),ring_mask=C,neighbors=ns,opposite=j,
            opposite_has_active_shortcut=0 in ES[chord]))
    neighbors=np.column_stack(neighbors);far=np.column_stack(far);survival=np.column_stack(survival)
    lags=sorted(float(x) for x in args.lags.split(','));width=len(gs)
    trajectory=Propagator(H,{'B':B},2.5).run(np.column_stack([X,Y]),lags);rows=[]
    for lag,(Z,err) in trajectory.items():
        for group in ['all','opposite_shortcut','opposite_no_shortcut']:
            ids=[i for i,d in enumerate(labels) if group=='all' or d['opposite_has_active_shortcut']==(group=='opposite_shortcut')]
            x,y=Z[:,ids],Z[:,np.array(ids)+width]
            ex=ew+float(np.linalg.norm(err['B'][ids]));ey=ew+float(np.linalg.norm(err['B'][np.array(ids)+width]))
            p=float(np.linalg.norm(X[:,ids])**2);pI=probability_interval(p,ew)
            row=dict(lag=lag,group=group,initial_probability=p,initial_probability_interval=pI,
                state_errors=[ex,ey],responses={})
            S=survival[:,ids];sp=float(np.linalg.norm(S*x)**2)
            row['same_ring_survival_interval']=ratio(probability_interval(sp,ex),pI)
            for name,allO in [('core_neighbors_distance1',neighbors),('opposite_distance2',far)]:
                O=allO[:,ids]
                for joint in [False,True]:
                    r=response(x,y,O*S if joint else O,ex,ey)
                    r['conditional_value']=r['joint_value']/p
                    r['conditional_interval']=ratio(r['joint_interval'],pI)
                    row['responses'][name+('_and_same_ring' if joint else '')]=r
            rows.append(row)
            print(json.dumps(dict(stage='core_response',lag=lag,group=group,
                neighbor=row['responses']['core_neighbors_distance1']['conditional_interval'],
                neighbor_and_ring=row['responses']['core_neighbors_distance1_and_same_ring']['conditional_interval'],
                same_ring=row['same_ring_survival_interval'])),flush=True)
    result=dict(status='PASS',N=4,Q=5,budget=args.budget,wait=args.wait,dimension=H.shape[0],
        waiting_state_error=ew,record_labels=labels,record_probabilities=probs.tolist(),
        complete_record_energy_change=record_energy,
        complete_record_energy_interval=[record_energy-workerr,record_energy+workerr],
        complete_record_entropy_bits=entropy,complete_record_entropy_error_bits=float(entropyerr),
        source='exp(-i theta n0) immediately after the complete twelve-outcome host/ring record.',
        receivers='Average density at the two fixed recorded-ring neighbours (distance1), and density at the fixed opposite vertex (distance2). Receiver labels and the ring pointer are fixed before evolution.',
        evolution='Full frozen H only, including all active off-ring hopping and graph rewiring.',
        scope='Auxiliary N4,Q5 recorded-ring control, not a Q=2N scaling point. Nonzero distance1 response with ring retention demonstrates finite conditional matter response, not that the ring alone controls it or that a local metric has emerged.',
        error_method='All outgoing full-H residuals; direct-sum Frobenius state error and bilinear response intervals. Complete record work/entropy included. Ordinary floating point plus5e-10 reporting allowance, not outward rounding.',
        rows=rows,seconds=time.monotonic()-start)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
