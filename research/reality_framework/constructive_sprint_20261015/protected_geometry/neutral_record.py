"""Neutral-origin dynamics followed by an explicitly costed geometry record.

Individual graph records break the source-fixing permutation sector. The wait
evolves in that exact quotient; recording and subsequent evolution use the full
physical budget basis and a verified signed isometric embedding.
"""
from __future__ import annotations
from itertools import combinations
from math import comb, log2, sqrt, lcm
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys
import time
import numpy as np
from scipy.sparse import coo_matrix,eye

P=Path(__file__).resolve().parent
RF=P.parents[1]
FROZEN=RF/'cp12_ground_connectivity'
sys.path.insert(0,str(FROZEN))
sys.path.insert(0,str(RF/'cp12_distance_history'))
from frozen_model import graph_basis
from exact_destination import transitions
from evolution import Evolution


def load_formation():
    spec=importlib.util.spec_from_file_location('grut_neutral_formation',RF/'cp12_unmasked_formation'/'run.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def positive(p,e):
    e+=5e-12
    return [min(1.,max(0.,sqrt(max(p,0.))-e)**2),min(1.,(sqrt(max(p,0.))+e)**2)]


def ratio(a,b):
    assert b[0]>0
    z=[x/y for x in a for y in b]
    return [min(z),max(z)]


def full_model(budget,quotient):
    sub,Hq,Bq,sq,sizes,psi=quotient
    gs=[sum(1<<e for e in c) for c in combinations(range(6),4)]
    states=[s for g in gs for s in graph_basis(4,g,int(sqrt(budget)))
            if sum(x*x for x in s[2])<=budget]
    lookup={s:i for i,s in enumerate(states)}
    outer={};rr=[];cc=[];vv=[];br=[];bc=[];bv=[]
    for col,s in enumerate(states):
        for t,a in transitions(4,s,120):
            if t in lookup:rr.append(lookup[t]);cc.append(col);vv.append(a/120)
            else:
                row=outer.setdefault(t,len(outer));br.append(row);bc.append(col);bv.append(a/120)
    H=coo_matrix((vv,(rr,cc)),shape=(len(states),len(states))).tocsr()
    B=coo_matrix((bv,(br,bc)),shape=(len(outer),len(states))).tocsr()
    assert np.linalg.norm((H-H.T).data)<1e-11
    index={s:i for i,s in enumerate(sq)};rr=[];cc=[];vv=[]
    for row,s in enumerate(states):
        key,sign,size=sub.canonical(s)
        if sign and key in index:rr.append(row);cc.append(index[key]);vv.append(sign/sqrt(size))
    J=coo_matrix((vv,(rr,cc)),shape=(len(states),len(sq))).tocsr()
    iso=float(np.linalg.norm((J.T@J-eye(len(sq))).data))
    intertwiner=float(np.linalg.norm((H@J-J@Hq).data))
    assert iso<1e-11 and intertwiner<1e-10
    return H,B,states,J,dict(quotient_dimension=len(sq),physical_dimension=len(states),
        physical_outgoing_dimension=len(outer),isometry_error=iso,H_intertwiner_error=intertwiner)


def circulation(s):
    g,b,E=s
    es=list(combinations(range(4),2))
    deg=[sum(g>>e&1 for e,pair in enumerate(es) if v in pair) for v in range(4)]
    if deg!=[2]*4:return False
    order=[0];previous=None
    while len(order)<4:
        current=order[-1]
        neighbors=[v for v in range(4) if v!=current and v!=previous
                   and g>>es.index(tuple(sorted((v,current))))&1]
        nxt=min(neighbors) if previous is None else neighbors[0]
        previous=current;order.append(nxt)
    signs=[]
    for a,c in zip(order,order[1:]+order[:1]):
        e=es.index(tuple(sorted((a,c))));leg=0 if a<c else 1
        if E[2*e]*E[2*e+1]>=0:return False
        signs.append(E[2*e+leg]>0)
    return all(signs) or not any(signs)


def main():
    started=time.monotonic()
    expected='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
    assert hashlib.sha256((FROZEN/'frozen_model.py').read_bytes()).hexdigest()==expected
    # A cache is optional only when the caller explicitly supplies
    # CP12_FORMATION_CACHE.  Fresh reproduction needs no workspace path.
    formation=load_formation();q=formation.model(4,4,256,source=True)
    sub,Hq,Bq,sq,sz,psi=q
    wait=.5*10;lag=.5
    Eq=Evolution(Hq,Bq,2.,steps=128)
    wq,_,ew=Eq(psi,wait)
    H,B,states,J,checks=full_model(256,q)
    w=J@wq
    n0=np.array([b&1 for g,b,E in states],float)
    Gcirc=np.array([circulation(s) for s in states],float)
    G=Gcirc*np.array([bool(b&15) for g,b,E in states],float)
    gs=sorted({g for g,b,E in states});mask=np.array([g for g,b,E in states])
    all_columns=[];outcome_labels=[]
    for g in gs:
        for flag in [0,1]:
            f=(mask==g)*(G if flag else 1-G)
            all_columns.append(f*w);outcome_labels.append((g,flag))
    complete=np.column_stack(all_columns)
    probs=np.linalg.norm(complete,axis=0)**2
    assert abs(probs.sum()-np.linalg.norm(w)**2)<1e-12
    entropy=float(-sum(p*log2(p) for p in probs if p>0))
    before=float(np.vdot(w,H@w).real)
    after=float(np.vdot(complete,H@complete).real)
    selected=[i for i,(g,flag) in enumerate(outcome_labels) if flag and probs[i]>1e-25]
    X=complete[:,selected];Y=n0[:,None]*X
    labels=[outcome_labels[i][0] for i in selected]
    pG=float(np.linalg.norm(X)**2);pI=positive(pG,float(ew))
    E=Evolution(H,B,2.,steps=32)
    k=X.shape[1]
    Z,individual,grouped=E(np.column_stack([X,Y]),lag,groups=[range(k),range(k,2*k)])
    x,y=Z[:,:k],Z[:,k:]
    ex=float(ew+grouped[0]);ey=float(ew+grouped[1])
    exn=float(np.linalg.norm(x));eyn=float(np.linalg.norm(y))
    response_error=2*(ex*eyn+ey*exn+ex*ey)+5e-12
    joint=hostsurv=circ=both=0.
    joint_host=joint_circ=joint_both=0.
    per_record=[]
    es=list(combinations(range(4),2))
    for i,g in enumerate(labels):
        near={b if a==0 else a for e,(a,b) in enumerate(es) if g>>e&1 and 0 in (a,b)}
        receiver=next(v for v in range(1,4) if v not in near)
        D=np.array([(bits>>receiver)&1 for gg,bits,EE in states],float)
        S=(mask==g).astype(float);C=Gcirc
        a=x[:,i];b=y[:,i]
        chi=float(2*np.vdot(a,D*b).imag)
        chiS=float(2*np.vdot(a,S*D*b).imag)
        chiC=float(2*np.vdot(a,C*D*b).imag)
        chiB=float(2*np.vdot(a,S*C*D*b).imag)
        pr=abs(a)**2
        hs=float(S@pr);cs=float(C@pr);bs=float((S*C)@pr)
        joint+=chi;hostsurv+=hs;circ+=cs;both+=bs
        joint_host+=chiS;joint_circ+=chiC;joint_both+=chiB
        per_record.append(dict(graph_mask=int(g),receiver_at_record_distance_two=receiver,
            initial_probability=float(np.linalg.norm(X[:,i])**2),
            joint_density_response=chi,same_host_joint_density_response=chiS,
            final_circulation_joint_density_response=chiC,
            same_host_and_circulation_joint_density_response=chiB))
    def response(v):
        I=[v-response_error,v+response_error]
        return dict(unnormalized_joint_value=v,unnormalized_interval=I,
            conditional_on_initial_record_event=v/pG,
            conditional_interval=ratio(I,pI),zero_excluded=I[0]*I[1]>0)
    def survival(v):
        I=positive(v,ex)
        conditional=ratio(I,pI)
        return dict(joint_probability=v,joint_interval=I,
                    conditional_probability=v/pG,
                    conditional_interval=[max(0.,conditional[0]),min(1.,conditional[1])])
    tv=min(1.,2*float(ew)+float(ew)**2)
    h2=lambda p:0. if p in (0.,1.) else -p*log2(p)-(1-p)*log2(1-p)
    entropy_error=h2(tv)+tv*log2(len(outcome_labels)-1)
    kinetic_norm=sqrt(2*4*4)
    rewire_offdiag_norm=.1/6*4*(6-4)
    work_error=2*(kinetic_norm+rewire_offdiag_norm)*(2*float(ew)+float(ew)**2)+5e-12
    result=dict(status='PASS',N=4,Q=4,electric_budget=256,wait=wait,lag=lag,
        initial_preparation='Original coherent neutral ray, equal amplitude on all15 Q4 graphs, no graph or loop mask.',
        initial_energy=2.,initial_variance=8.,model_checks=checks,
        event='At waiting time, all four active connections belong to one directed electric circulation and nv>=1. This forces a spanning ring at N4,Q4.',
        event_probability_at_wait=pG,event_probability_interval=pI,
        complete_record='Orthogonal projective channel {P_g G, P_g(1-G)} for every labelled graph g; all outcomes retained. G1 branches are analyzed conditionally afterward.',
        complete_record_outcome_probabilities=[dict(graph_mask=int(g),event_flag=flag,probability=float(p)) for (g,flag),p in zip(outcome_labels,probs)],
        complete_record_entropy_bits=entropy,entropy_error_estimate_bits=entropy_error,
        energy_before_complete_record=before,energy_after_complete_record=after,
        complete_record_energy_cost=after-before,complete_record_energy_cost_error_estimate=work_error,
        record_energy_error_method='A=Phi(H)-H=Phi(K+Roff)-(K+Roff) is bounded because electric, mass and diagonal rewiring terms commute with the complete record. ||A||<=2(sqrt(2NQ)+Gamma Q(M-Q)/M); expectation uncertainty <=||A||(2 ew+ew^2), plus5e-12 reporting tolerance. No unbounded-energy expectation is bounded merely from wavefunction norm.',
        record_entropy_error_method='The full pure waiting state gives orthogonal complete-record branches, so pinched measured-system entropy equals Shannon entropy of30 outcome probabilities. Conservative total variation <=2 ew+ew^2; continuity error h2(TV)+TV log2(29).',
        wait_state_error=float(ew),lag_X_Frobenius_error=float(grouped[0]),lag_Y_Frobenius_error=float(grouped[1]),
        full_X_Frobenius_error=ex,full_Y_Frobenius_error=ey,
        source='Gauge-invariant exp(-i theta n0) after recording; commutes with every record projector, so initial G probability has zero source derivative.',
        detector='Receiver is opposite source0 on the recorded ring. Its vertex density is measured independently of the final graph; no final-distance reclassification.',
        density_response=response(joint),same_host_joint_density_response=response(joint_host),
        final_circulation_joint_density_response=response(joint_circ),
        same_host_and_circulation_joint_density_response=response(joint_both),
        host_survival=survival(hostsurv),circulation_survival=survival(circ),
        joint_host_and_circulation_survival=survival(both),per_record=per_record,
        recorded_state_entropy_scope='Shannon entropy equals measured-system von Neumann entropy because this complete projective channel acts on a pure state with mutually orthogonal branch vectors. A coherent record apparatus can retain a globally pure enlarged state; apparatus preparation/reset cost is not derived.',
        error_scope='Full unbounded-law electric truncation controlled by all outgoing residuals and analytic Duhamel upper sums; ordinary floating-point evaluation with separate5e-12 reporting tolerance, no formal interval arithmetic.',
        interpretation='Neutral dynamics naturally reaches the rare G sector; an explicit external diagnostic record selects graph-relative detector pointers and changes coherence/energy. This is an operational recorded-history experiment, not autonomous selection of a geometric phase, macroscopic locality or three-dimensional emergence.',
        seconds=time.monotonic()-started)
    (P/'NEUTRAL_RECORD.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
