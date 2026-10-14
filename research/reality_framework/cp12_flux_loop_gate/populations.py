"""Decompose the already published states; no new main-state evolution."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from loops import P,REL,UP,model,resources,observables

CASES=[(4,4,256,[0.,1.,3.,5.,10.]),(5,5,128,[0.,1.,3.,5.,10.]),
       (6,6,128,[0.,1.,3.,5.,10.]),(6,12,10,[0.,.4,.6,1.])]

SCOPE='Existing neutral-sea states at identical published times, decomposed without a physical intervention. Hole-supported means at least one hole on a cycle. Occupied-edge-flux cycles are undirected supports; their directed subset is separately reported. Electric circulation uses both nonzero opposite-sign half-leg fields and can include hole edges. The two hole/occupied channel events can coexist; disjoint hole-only, flux-only and both events sum to protected-cycle probability. Positive-score intervals use [max(0,sqrt(p)-e)^2,min(1,(sqrt(p)+e)^2)] for every positive contraction. Analytic cutoff bounds evaluated in ordinary floating point; no interval roundoff certificate.'


def run(n,q,budget,times,snapshot_dir):
    started=time.monotonic();sub,H,B,states,sizes,psi=model(n,q,budget)
    names,O,mm=observables(n,q,states)
    r,k,nv,e2=resources(n,q,sub,states)
    source=REL/f'RELATIONAL_N{n}_Q{q}_B{budget}.json'
    old=json.loads(source.read_text());trajectory=json.loads((UP/old['input_trajectory']).read_text())
    rows=[];difference=0.
    for t in times:
        path=Path(snapshot_dir)/f'N{n}_Q{q}_B{budget}_T{t:g}.npz'
        v=np.load(path)['wave'];p=abs(v)**2
        ref=min(old['rows'],key=lambda a:abs(a['T']-t));assert abs(ref['T']-t)<1e-9
        resource=min(trajectory['rows'],key=lambda a:abs(a['T']-t))
        vals=p@O;data=dict(zip(names,map(float,vals)))
        error=max(abs(data['protected_cycle_any']-ref['means']['has_cycle']),
                  abs(data['vacuum_closed_flux_cycle']-resource['vacuum_flux_probability']),
                  abs(float(p@k/q)-ref['flux_protected_fraction']))
        difference=max(difference,error);assert error<1e-9
        e=ref['state_error_upper']
        intervals={name:[max(0.,float(np.sqrt(max(0.,value))-e))**2,
                         min(1.,float(np.sqrt(max(0.,value))+e)**2)] for name,value in data.items()}
        rows.append(dict(T=t,values=data,state_error_upper=ref['state_error_upper'],
            bounded_score_error_upper=ref['bounded_operator_error_upper'],
            positive_score_intervals=intervals,
            main_state_full_mean_energy=q/2,main_state_global_entropy_bits=0,
            snapshot_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            projected_norm_drift=float(abs(p.sum()-1))))
        print(json.dumps(dict(N=n,Q=q,T=t,protected=data['protected_cycle_any'],
            flux=data['occupied_flux_cycle_any'],hole_only=data['hole_only_cycle_event'],
            directed=data['directed_flux_cycle_any'],spanning=data['directed_cycle_network_spanning'])),flush=True)
    out=dict(N=n,Q=q,budget=budget,dimension=len(states),metrics=names,rows=rows,
        previous_relational_record_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        maximum_previous_record_difference=difference,seconds=time.monotonic()-started,
        scope=SCOPE)
    (P/f'CHANNELS_N{n}_Q{q}_B{budget}.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    a=ap.parse_args()
    for case in CASES:run(*case,a.snapshot_dir)
