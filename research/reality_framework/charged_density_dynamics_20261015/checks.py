"""Independent controls for symmetry, circulation and record accounting."""
from __future__ import annotations
from pathlib import Path
from itertools import combinations,permutations
import importlib.util,json,sys
import numpy as np
from formation import model
from diagnostics import arrays,event,cycle_edges
from stream import StreamEvolution
from fast_symmetry import FastInvariant

P=Path(__file__).resolve().parent
RF=P.parent
spec=importlib.util.spec_from_file_location('density_control_reference',RF/'cp12_unmasked_formation/run.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from exact_destination import transitions


def explicit_cycles(n,state,balanced=False):
    g,b,E=state;es=list(combinations(range(n),2));arcs={}
    for e,(i,j) in enumerate(es):
        a,c=E[2*e:2*e+2]
        if g>>e&1 and a*c<0 and (not balanced or a==-c):arcs[(i,j) if a>0 else (j,i)]=e
    answer=set()
    for length in range(3,n+1):
        for vs in combinations(range(n),length):
            a=min(vs)
            for rest in permutations(v for v in vs if v!=a):
                walk=(a,)+rest;links=list(zip(walk,walk[1:]+walk[:1]))
                if all(e in arcs for e in links):answer|={arcs[e] for e in links}
    return answer


def main():
    matrix=[];count=0
    for n,q,b,source in [(3,2,6,False),(4,4,12,False),(4,4,12,True),(5,10,8,False)]:
        fast=model(n,q,b,source);slow=old.model(n,q,b,source)
        sub,H,B,states,sizes,psi=fast;ref,Ho,Bo,so,szo,po=slow
        assert states==so and np.array_equal(sizes,szo)
        h=float(np.linalg.norm((H-Ho).data));leak=float(np.linalg.norm((B.T@B-Bo.T@Bo).data))
        assert h<1e-11 and leak<1e-11
        for state in states:
            assert sub.canonical(state)==ref.canonical(state)
            for balanced in [False,True]:
                assert {e for e,u,v in cycle_edges(n,state,balanced)}==explicit_cycles(n,state,balanced)
            count+=1
        matrix.append(dict(N=n,Q=q,budget=b,source=source,dimension=len(states),hamiltonian_error=h,outgoing_gram_error=leak))
    # Exact equivalence to the published charged N4 prototype; and comparison
    # of streamwise leakage upper sums with the inherited full-grid routine.
    sub,H,B,states,sizes,psi=model(4,4,256)
    record_spec=importlib.util.spec_from_file_location('density_legacy_record',RF/'constructive_sprint_20261015/protected_geometry/neutral_record.py')
    record=importlib.util.module_from_spec(record_spec);record_spec.loader.exec_module(record)
    d=arrays(4,4,states);mask=event(d,4,4,1,.25)
    expected=np.asarray([record.circulation(s) and bool(s[1]&15) for s in states])
    assert np.array_equal(mask,expected)
    assert not np.any((d[2]==4)&(d[1]>0))
    sys.path.insert(0,str(RF/'cp12_distance_history'));from evolution import Evolution
    w,_,e=Evolution(H,B,2.,steps=64)(psi,1.)
    ws,es=StreamEvolution(H,B,2.,max_dt=1/64).run(psi,[1.])[1.]
    stream_error=float(np.linalg.norm(w-ws));assert stream_error<1e-11
    # Non-hopping terms preserve the corrected circulation and matter event.
    commute_checks=0
    for state in states:
        flag=record.circulation(state) and bool(state[1]&15)
        for target,c in transitions(4,state,120,hopping=False):
            assert flag==(record.circulation(target) and bool(target[1]&15));commute_checks+=1
    saved=json.loads((P/'N4_Q4_B256_A1_V0.25.json').read_text())
    old_result=json.loads((RF/'constructive_sprint_20261015/protected_geometry/NEUTRAL_RECORD.json').read_text())
    rc=saved['complete_graph_event_record']
    energy=abs(rc['energy_change']-old_result['complete_record_energy_cost'])
    entropy=abs(rc['entropy_bits']-old_result['complete_record_entropy_bits'])
    assert energy<1e-10 and entropy<1e-10
    result=dict(status='PASS',matrix_controls=matrix,explicit_cycle_controls=count,
        corrected_N4_event_matches=True,strict_balanced_charged_full_event_empty=True,
        nonhopping_event_commutation_checks=commute_checks,
        streaming_state_difference=stream_error,
        graph_record_energy_difference_from_full_physical_calculation=energy,
        graph_record_entropy_difference_from_full_physical_calculation=entropy,
        scope='Independent legacy Python orbit/matrix/outgoing-Gram comparison, explicit directed-cycle enumeration, full-grid/stream propagation, and prior full-physical record-channel comparison. Internal computational checks, not external peer review.')
    (P/'CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
