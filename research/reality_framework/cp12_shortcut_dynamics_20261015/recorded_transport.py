"""Neutral-origin CP-12 transport on a recorded ring with an active diagonal.

N=4,Q=5 is the smallest spanning ring with an extra active connection. Every
labelled host is K4 minus one edge and contains a unique four-cycle. No flux
or matter state is supplied: the original neutral ray first evolves under H.
The diagnostic record and the subsequent comparison intervention are explicit.
"""
from pathlib import Path
from itertools import combinations
from math import sqrt
import argparse,hashlib,json,sys,time
import numpy as np
from scipy.sparse import vstack
from operators import physical_model,flow,RF
from integrals import Propagator
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
from stream import probability_interval

ALLOWANCE=5e-10
ES=list(combinations(range(4),2))


def ring_for(g):
    missing=next(e for e in range(6) if not (g>>e&1))
    a,b=ES[missing];other=tuple(v for v in range(4) if v not in (a,b))
    chord=ES.index(other);C=g^(1<<chord)
    assert C.bit_count()==4
    return C,chord


def ring_present(s,C):
    g,bits,E=s;out=[0]*4;inc=[0]*4
    if g&C!=C:return False
    for e,(a,b) in enumerate(ES):
        if not(C>>e&1):continue
        x,y=E[2*e:2*e+2]
        if x*y>=0:return False
        u,v=(a,b) if x>0 else (b,a);out[u]+=1;inc[v]+=1
    return out==inc==[1]*4


def selected(s):
    g,bits,E=s;C,chord=ring_for(g)
    return bool(bits&15) and ring_present(s,C) and E[2*chord]==E[2*chord+1]==0


def opposite(C,source=0):
    neighbors={b if a==source else a for e,(a,b) in enumerate(ES) if C>>e&1 and source in (a,b)}
    return next(v for v in range(4) if v!=source and v not in neighbors)


def ratio(I,J):
    if J[0]<=0:return None
    return [min(x/y for x in I for y in J),max(x/y for x in I for y in J)]


def response(x,y,O,ex,ey):
    val=float(2*np.vdot(x,O*y).imag)
    radius=2*(ex*np.linalg.norm(O*y)+ey*np.linalg.norm(O*x)+ex*ey)+ALLOWANCE
    plus=x+(-1j-1)*y;minus=x+(1j-1)*y
    phase=float((np.vdot(plus,O*plus).real-np.vdot(minus,O*minus).real)/2)
    assert abs(phase-val)<1e-12
    theta=.001
    plus=x+np.expm1(-1j*theta)*y;minus=x+np.expm1(1j*theta)*y
    finite=float((np.vdot(plus,O*plus).real-np.vdot(minus,O*minus).real)/(2*theta))
    return dict(joint_value=val,joint_interval=[val-radius,val+radius],error_bound=radius,
        zero_excluded=bool(abs(val)>radius),phase_shift_identity_error=abs(phase-val),
        finite_theta=theta,finite_pulse_susceptibility=finite,finite_pulse_difference=finite-val)


def main():
    p=argparse.ArgumentParser();p.add_argument('--budget',type=int,required=True)
    p.add_argument('--wait',type=float,default=5.);p.add_argument('--lags',default='.125,.25,.5,1')
    p.add_argument('--dt',type=float,default=1/32);p.add_argument('--out',required=True)
    args=p.parse_args();started=time.monotonic();m=physical_model(4,5,args.budget)
    H,B,states,psi=(m[k] for k in ['H','B','states','psi']);Ks=m['Ks'];BKs=m['BKs']
    w,er=Propagator(H,{'electric':B},2.5,dt=args.dt).run(psi,[args.wait])[args.wait]
    ew=float(er['electric'][0]);gs=sorted({s[0] for s in states})
    host=np.array([s[0] for s in states]);G=np.array([selected(s) for s in states])
    n0=np.array([s[1]&1 for s in states],float)
    nv=np.array([(s[1]&15).bit_count() for s in states])
    all_record=np.column_stack([(host==g)*(G if flag else ~G)*w for g in gs for flag in [0,1]])
    probabilities=np.linalg.norm(all_record,axis=0)**2
    X=np.column_stack([(host==g)*G*w for g in gs]);Y=n0[:,None]*X
    eventp=float(np.linalg.norm(X)**2);pI=probability_interval(eventp,ew)
    deltaE=float(np.vdot(all_record,H@all_record).real-np.vdot(w,H@w).real)
    b=sqrt(40)+.1*5/6
    work_error=2*b*(2*ew+ew*ew)+ALLOWANCE
    entropy=float(-sum(a*np.log2(a) for a in probabilities if a>0))
    tv=min(1.,2*ew+ew*ew)
    h2=lambda a:0. if a<=0 or a>=1 else -a*np.log2(a)-(1-a)*np.log2(1-a)
    entropy_error=min(np.log2(12),h2(tv)+tv*np.log2(11)) if tv<=11/12 else np.log2(12)
    labels=[]
    for g in gs:
        C,chord=ring_for(g);j=opposite(C)
        labels.append(dict(graph_mask=g,ring_mask=C,chord_edge=chord,receiver=j,
            ring_distance=2,active_graph_distance=1 if 0 in ES[chord] else 2,
            shortcut=0 in ES[chord],initial_probability=float(np.linalg.norm(X[:,gs.index(g)])**2)))
    lags=sorted(float(s) for s in args.lags.split(','));width=len(gs)
    full=Propagator(H,{'electric':B},2.5,dt=args.dt).run(np.column_stack([X,Y]),lags)
    core={t:dict(x=np.empty_like(X),y=np.empty_like(Y),ex=np.empty(width),ey=np.empty(width),
        outside=np.empty(width)) for t in lags}
    for C in sorted({d['ring_mask'] for d in labels}):
        indices=[i for i,d in enumerate(labels) if d['ring_mask']==C]
        V=sum((Ks[e] for e in range(6) if not(C>>e&1)),start=0*H)
        BV=sum((BKs[e] for e in range(6) if not(C>>e&1)),start=0*B)
        control=Propagator(H-V,{'electric':B-BV,'outside_hopping':vstack([V,BV],format='csr')},
            2.5,dt=args.dt).run(np.column_stack([X[:,indices],Y[:,indices]]),lags)
        k=len(indices)
        for lag,(z,res) in control.items():
            core[lag]['x'][:,indices]=z[:,:k];core[lag]['y'][:,indices]=z[:,k:]
            core[lag]['ex'][indices]=res['electric'][:k];core[lag]['ey'][indices]=res['electric'][k:]
            core[lag]['outside'][indices]=res['outside_hopping'][:k]
    rows=[]
    for lag in lags:
        z,res=full[lag];xf,yf=z[:,:width],z[:,width:]
        for shortcut in [True,False]:
            ids=[i for i,d in enumerate(labels) if d['shortcut']==shortcut]
            p0=float(np.linalg.norm(X[:,ids])**2);p0I=probability_interval(p0,ew)
            O=np.column_stack([np.array([(s[1]>>labels[i]['receiver'])&1 for s in states],float) for i in ids])
            S=np.column_stack([np.array([ring_present(s,labels[i]['ring_mask']) for s in states],float) for i in ids])
            SH=np.column_stack([(host==labels[i]['graph_mask']).astype(float) for i in ids])
            new_chord=np.column_stack([np.array([s[2][2*labels[i]['chord_edge']]*s[2][2*labels[i]['chord_edge']+1]<0
                for s in states],float) for i in ids])
            unequal_ends=np.column_stack([np.array([((s[1]>>ES[labels[i]['chord_edge']][0])&1)!=
                ((s[1]>>ES[labels[i]['chord_edge']][1])&1) for s in states],float) for i in ids])
            onset4=float(np.linalg.norm(unequal_ends*X[:,ids])**2)/4
            onsetI=[a/4 for a in probability_interval(4*onset4,ew)]
            vacancies=np.column_stack([np.array([2-((s[1]>>ES[labels[i]['chord_edge']][0])&1)-
                ((s[1]>>ES[labels[i]['chord_edge']][1])&1) for s in states],float) for i in ids])
            initial_Vnorm=float(np.sqrt(np.sum(vacancies*abs(X[:,ids])**2)))
            out=dict(lag=lag,shortcut_present_at_record=shortcut,ring_distance=2,
                active_graph_distance=1 if shortcut else 2,initial_joint_probability=p0,
                initial_joint_probability_interval=p0I,
                exact_initial_outside_action_norm=initial_Vnorm,
                exact_initial_outside_action_norm_interval=[max(0.,initial_Vnorm-4*ew-ALLOWANCE),initial_Vnorm+4*ew+ALLOWANCE],
                initial_outside_action_norm_per_initial_norm=initial_Vnorm/sqrt(p0),
                new_chord_and_ring_probability_T4_coefficient_joint=onset4,
                new_chord_and_ring_probability_T4_coefficient_joint_interval=onsetI,
                new_chord_and_ring_probability_T4_coefficient_conditional=onset4/p0,
                new_chord_and_ring_probability_T4_coefficient_conditional_interval=ratio(onsetI,p0I),models={})
            stored={}
            for label in ['full_frozen','recorded_ring_hopping_only']:
                if label=='full_frozen':
                    x,y=xf[:,ids],yf[:,ids]
                    ex=ew+float(np.linalg.norm(res['electric'][ids]));ey=ew+float(np.linalg.norm(res['electric'][np.array(ids)+width]))
                else:
                    x,y=core[lag]['x'][:,ids],core[lag]['y'][:,ids]
                    ex=ew+float(np.linalg.norm(core[lag]['ex'][ids]));ey=ew+float(np.linalg.norm(core[lag]['ey'][ids]))
                r=response(x,y,O,ex,ey);r['conditional_value']=r['joint_value']/p0
                r['conditional_interval']=ratio(r['joint_interval'],p0I)
                rs=response(x,y,O*S,ex,ey);rs['conditional_value']=rs['joint_value']/p0
                rs['conditional_interval']=ratio(rs['joint_interval'],p0I)
                surv=float(np.linalg.norm(S*x)**2);survI=probability_interval(surv,ex)
                same=float(np.linalg.norm(SH*x)**2);sameI=probability_interval(same,ex)
                chordp=float(np.linalg.norm(S*new_chord*x)**2);chordI=probability_interval(chordp,ex)
                values=dict(receiver_response=r,receiver_and_initial_ring_joint_response=rs,
                    initial_ring_survival_joint=surv,initial_ring_survival_conditional=surv/p0,
                    initial_ring_survival_conditional_interval=ratio(survI,p0I),
                    host_survival_joint=same,host_survival_conditional=same/p0,
                    host_survival_conditional_interval=ratio(sameI,p0I),
                    new_chord_and_ring_joint_probability=chordp,
                    new_chord_and_ring_joint_interval=chordI,
                    new_chord_and_ring_conditional_probability=chordp/p0,
                    new_chord_and_ring_conditional_interval=ratio(chordI,p0I),
                    full_X_state_error=ex,full_Y_state_error=ey,
                    final_receiver_density_joint=float(np.vdot(x,O*x).real))
                out['models'][label]=values;stored[label]=(x,y,ex,ey)
            af=out['models']['full_frozen']['receiver_response'];ac=out['models']['recorded_ring_hopping_only']['receiver_response']
            diff=af['joint_value']-ac['joint_value'];differr=af['error_bound']+ac['error_bound']
            out['response_change_joint']=diff;out['response_change_joint_interval']=[diff-differr,diff+differr]
            out['response_change_conditional']=diff/p0;out['response_change_conditional_interval']=ratio(out['response_change_joint_interval'],p0I)
            xx,yy,ex,ey=stored['full_frozen'];xc,yc,exc,eyc=stored['recorded_ring_hopping_only']
            dist=float(np.linalg.norm(xx-xc));radius=ex+exc
            out['state_distance']=dist;out['state_distance_interval']=[max(0.,dist-radius-ALLOWANCE),dist+radius+ALLOWANCE]
            # Each classical record may carry its own physically irrelevant
            # global phase; optimize those phases separately before combining.
            ray2=sum(max(0.,np.linalg.norm(xx[:,k])**2+np.linalg.norm(xc[:,k])**2-
                2*abs(np.vdot(xx[:,k],xc[:,k]))) for k in range(len(ids)))
            ray=float(np.sqrt(ray2));out['phase_optimized_state_distance']=ray
            out['phase_optimized_state_distance_interval']=[max(0.,ray-radius-ALLOWANCE),ray+radius+ALLOWANCE]
            out['normalized_state_distance']=dist/sqrt(p0)
            out['normalized_state_distance_interval']=[out['state_distance_interval'][0]/(sqrt(p0)+ew),
                min(2.,out['state_distance_interval'][1]/max(sqrt(p0)-ew,1e-300))]
            Iv=float(np.linalg.norm(core[lag]['outside'][ids]));ec=float(np.linalg.norm(core[lag]['ex'][ids]))
            out['outside_hopping_integral_upper_sum']=Iv
            out['full_law_Duhamel_upper']=Iv+4*lag*(ew+ec)
            rows.append(out)
            print(json.dumps(dict(stage='recorded_response',budget=args.budget,lag=lag,shortcut=shortcut,
                probability=p0,full_chi=af['conditional_value'],core_chi=ac['conditional_value'],
                difference=out['response_change_conditional_interval'],ring_survival=out['models']['full_frozen']['initial_ring_survival_conditional_interval'])),flush=True)
    result=dict(status='PASS',N=4,Q=5,electric_budget=args.budget,dimension=H.shape[0],
        outgoing_dimension=B.shape[0],H_nonzeros=H.nnz,wait=args.wait,wait_state_error=ew,
        initial_energy=2.5,initial_variance=10.,neutral_energy_check=float(np.vdot(psi,H@psi).real),
        neutral_variance_check=float(np.linalg.norm(H@psi-2.5*psi)**2),
        initial_preparation='Original coherent neutral ray over all SIX labelled Q5 hosts. No supplied loop, flux, matter or host.',
        complete_record='Twelve projectors {P_g G_g, P_g (1-G_g)}. G_g requires a directed four-cycle, nv>=1 and a neutral fifth edge. Every Q5 host has a unique spanning four-cycle.',
        record_labels=labels,record_outcome_probabilities=probabilities.tolist(),
        event_probability=eventp,event_probability_interval=pI,
        complete_record_energy_change=deltaE,complete_record_energy_interval=[deltaE-work_error,deltaE+work_error],
        complete_record_entropy_bits=entropy,complete_record_entropy_error_bits=float(entropy_error),
        comparison='For each recorded ring C, V_C=sum_{e not in C} K_e on the FULL physical Hilbert space; H_C=H-V_C. The pointer C stays fixed, all graph rewiring and all electric/mass terms remain. This is an external diagnostic intervention, not a changed GRUT law.',
        direct_charge_transfer='For e=(i,j), T_ji=a_j^dagger U_je U_ie^dagger a_i is gauge invariant. K_e^2 contains (2n_b,e-1)(T_ji+T_ij). On the recorded neutral fifth edge n_b=1,E=0, the two-hop transfer amplitude has magnitude lag^2/2 with no gauge cancellation. Its final fields are (-1,+1) or (+1,-1), while the original ring fields are unchanged.',
        exact_shortcut_onset='For each recorded branch, Pr(new fifth-edge through-flow AND original ring) = lag^4 <n_i(1-n_j)+n_j(1-n_i)>/4 + O(lag^5). The measured finite-wait state fixes the expectation; its coefficient interval includes waiting-state error. This finite derivative does not assert a long-time rate. The corresponding probability is exactly zero under H_C because all non-C half-leg fields stay zero.',
        source='Gauge-invariant exp(-i theta n_0) immediately after the record. It commutes with all twelve projectors, so the initial event probability has zero derivative.',
        receiver='The vertex opposite source0 on the initial recorded ring, fixed before the lag, independent of the final graph. A joint classical record stores its label.',
        density_contrast_limit='Shortcut and no-shortcut categories have source/receiver degree3 and degree2 respectively. Their direct comparison is not a degree-matched scaling test. Within each full-versus-core comparison the initial state and recorded host are identical.',
        hopping_off_response='Exactly zero for all receivers: [n_j,H_diagonal+H_rewire]=0 for every j.',
        results=rows,
        error_method='Full physical electric-budget basis and ALL outgoing rows. Analytic midpoint residual upper sums, Duhamel state bounds, and bilinear response intervals. Errors for orthogonal record branches combine in Frobenius norm; waiting-state error enters each direct-sum group once. Ordinary floating point plus5e-10 reporting allowance, not outward rounding.',
        frozen_sources={name:hashlib.sha256((RF/'cp12_ground_connectivity'/name).read_bytes()).hexdigest()
            for name in ['frozen_model.py','exact_destination.py']},seconds=time.monotonic()-started)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='done_record',out=args.out,event_probability=eventp,event_interval=pI,
        wait_error=ew,dimension=H.shape[0],seconds=result['seconds'])),flush=True)


if __name__=='__main__':main()
