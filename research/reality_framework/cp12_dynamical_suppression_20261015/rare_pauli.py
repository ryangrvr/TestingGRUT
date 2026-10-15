"""A neutral-origin, rare Pauli-filled ring; full frozen dynamics afterwards.

The diagonal event has a one-dimensional range inside the exact S5 invariant
sector. Its normalized conditional shape is therefore exact, even though the
formation probability retains the separate full-law waiting error interval.
This is an exploratory mechanism control, not an appreciable geometric phase.
"""
from pathlib import Path
from itertools import combinations
import json,sys,time
import numpy as np
from boundary import quotient_model,RF
from mechanisms import statistics,action_operators
from action_integrals import ActionEvolution,ALLOWANCE
from operators import flow
from integrals import Propagator
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
from stream import probability_interval

P=Path(__file__).resolve().parent


def ratio(I,J):
    return [I[0]/J[1],I[1]/J[0]] if J[0]>0 else [0.,None]


def event(s):
    g,b,E=s;cm=flow(5,s)[0];degree=[0]*5
    for e,(i,j) in enumerate(combinations(range(5),2)):
        if cm>>e&1:degree[i]+=1;degree[j]+=1
    return ((b&31)==31 and cm.bit_count()==5 and degree==[2]*5 and sum(a*a for a in E)==25 and
        all(E[2*e]==E[2*e+1]==0 for e in range(10) if not(cm>>e&1)))


def main():
    start=time.monotonic();m=quotient_model(5,10,50);H,B=m['H'],m['B']
    mask=np.array([event(s) for s in m['states']]);ids=np.flatnonzero(mask)
    assert len(ids)==1 and m['sizes'][ids[0]]==24
    phi=np.zeros(len(mask),complex);phi[ids[0]]=1
    w,er=Propagator(H,{'B':B},5.).run(m['psi'],[2.])[2.];ew=float(er['B'][0])
    p=float(np.linalg.norm(mask*w)**2);pI=probability_interval(p,ew)
    assert pI[0]>0
    actions=action_operators(m,5);ds,db=statistics(m,5,10,'N5_B50')
    ds['filled_ring_event']=mask.astype(float);db['filled_ring_event']=1.
    nv=np.array([(s[1]&31).bit_count() for s in m['states']]);ds['full_vertex_filling']=(nv==5).astype(float);db['full_vertex_filling']=1.
    out0=sum(np.linalg.norm(L@phi)**2 for L in actions['outside'])
    in0=sum(np.linalg.norm(L@phi)**2 for L in actions['inside'])
    vh=sum(np.linalg.norm(L@(H@phi))**2 for L in actions['outside'])
    assert out0==0 and abs(in0-10)<1e-12 and abs(vh-20)<1e-11
    assert np.linalg.norm(B@phi)==0
    results=ActionEvolution(H,B,actions,ds,db,10,5.,dt=1/64).run(phi,[.0625,.125,.25,.5,1.])
    rows=[]
    for t,d in results.items():
        ints=d['integrals'];eta=ints['outside']['value']/ints['inside']['value']
        rows.append(dict(lag=t,integrals=ints,eta=eta,
            eta_interval=ratio(ints['outside']['full_law_interval'],ints['inside']['full_law_interval']),
            final=d['final'],state_error=d['state_error'],phase_pulse_fixed_receiver_response=0.))
        print(json.dumps(dict(stage='rare_Pauli',lag=t,eta=eta,interval=rows[-1]['eta_interval'],
            outside=ints['outside']['full_law_interval'],inside=ints['inside']['full_law_interval'])),flush=True)
    x=mask*w;z=(1-mask)*w;de=-2*float(np.vdot(x,H@z).real)
    Aw=mask*(H@z)+(1-mask)*(H@x)
    work_error=2*ew*float(np.linalg.norm(Aw))+10*ew*ew+ALLOWANCE
    # All selected fields have electric square25; their one-hop neighbours have
    # square at most28. P_out=0 and B P=0, so the full A w_hat has no exit part.
    h2=lambda a:0. if a<=0 or a>=1 else -a*np.log2(a)-(1-a)*np.log2(1-a)
    entropy=[h2(pI[0]),h2(pI[1])]
    result=dict(status='PASS',N=5,Q=10,budget=50,wait=2.,
        formation_probability=p,formation_interval=pI,waiting_state_error=ew,
        conditional_shape='Exactly one S5 invariant orbit, size24: all vertices filled, a spanning charged five-cycle at unit circulation, all five non-core edges neutral. The shape follows from symmetry and the event; it is not supplied at time0. Its global phase is immaterial.',
        complete_binary_record_energy_change=de,complete_binary_record_energy_interval=[de-work_error,de+work_error],
        complete_binary_record_entropy_interval_bits=entropy,
        conditional_initial_state_error=0.,
        initial_outside_action_square=out0,initial_inside_action_square=in0,
        outside_H_action_square=vh,eta_T2_coefficient=vh/(3*in0),
        phase_response_identity='n0 phi=phi for every fixed source vertex0. An immediate local density phase pulse is a global phase, so every fixed-density receiver susceptibility is exactly zero at every subsequent lag. This does not assert absence of all motion or of responses to other source operators.',
        rare_event_scope='Exploratory neutral-origin postselection with complete record resources. Formation weight is tiny and not an appreciable phase; it cannot satisfy the requested locality gate merely because eta is initially small.',
        rows=rows,error_scope='All outgoing full-H electric residuals, analytic squared-action quadrature and ordinary floating-point reporting allowance; no hopping deletion.',
        seconds=time.monotonic()-start)
    (P/'N5_RARE_PAULI.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='rare_Pauli_done',formation_interval=pI,record_energy=result['complete_binary_record_energy_interval'])),flush=True)


if __name__=='__main__':main()
