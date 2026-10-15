"""No initial loop filter: finite CP12 response and record-interference audit.

The frozen Hamiltonian and supplied neutral ray are unchanged. The only
physical intervention on the main path is a gauge-invariant density phase.
Pinched/selected paths are explicit counterfactual diagnostics, not free
preparations. All electric truncations retain the full outgoing residual.
"""
from pathlib import Path
from itertools import combinations
import argparse, hashlib, json, sys, time
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
FLUX=P.parent/'cp12_flux_loop_gate'
sys.path.insert(0,str(FLUX))
from loops import model, circuits, masks, incidence_masks, UP, REL
from references import inclusion, energy_basis
from persistence import positive_interval, h2
from evolution import Evolution


def subtract(a,b):
    return [float(a[0]-b[1]),float(a[1]-b[0])]


def multiply(a,b):
    x=[u*v for u in a for v in b]
    return [float(min(x)),float(max(x))]


def divide(a,b):
    if b[0]<=0:return None
    return multiply(a,[1/b[1],1/b[0]])


def detectors(n,states):
    """All are invariant under permutations fixing source vertex zero."""
    es=list(combinations(range(n),2));outside=[c for c in circuits(n) if 0 not in c[3]]
    source=np.array([s[1]&1 for s in states],float)
    density=np.array([sum(s[1]>>j&1 for j in range(1,n))/(n-1) for s in states])
    incident=[e for e,(i,j) in enumerate(es) if i==0]
    nonincident=[e for e,(i,j) in enumerate(es) if i!=0]
    h_inc=np.array([sum(1-(s[1]>>(n+e)&1) for e in incident)/(n-1) for s in states])
    h_non=np.array([sum(1-(s[1]>>(n+e)&1) for e in nonincident)/len(nonincident) for s in states])
    ward=source+(n-1)*density-(n-1)*h_inc-len(nonincident)*h_non
    assert np.max(abs(ward))<1e-12
    occupied=[];circulation=[]
    for s in states:
        protected,flux,holes,positive=masks(n,s)
        through,sign=incidence_masks(n,s)
        occupied.append(any(flux&c==c for _,c,_,_ in outside))
        circulation.append(any(through&c==c and sign&c in (a,c^a) for _,c,a,_ in outside))
    F=np.array(occupied,float);C=np.array(circulation,float)
    return dict(density=density,circulation=C,circulation_density=C*density,
                source_density=source,incident_edge_holes=h_inc,other_edge_holes=h_non),dict(occupied_flux=F,electric_circulation=C)


def propagate(E,V,lag,step=.1):
    steps=max(1,int(np.ceil(lag/step))) if lag else 0
    errors=np.zeros(V.shape[1] if V.ndim==2 else 1)
    z=V.copy()
    for _ in range(steps):
        z,de,_=E(z,lag/steps)
        errors+=de
    return z,errors


def summarize_paths(paths):
    plus,minus=paths['plus'],paths['minus'];zero=paths['zero']
    slopes={}
    for name in zero['values']:
        a=plus['values'][name];b=minus['values'][name]
        I=subtract(plus['intervals'][name],minus['intervals'][name])
        slopes[name]=dict(value=(a-b)/2,interval=[I[0]/2,I[1]/2])
    f=zero['values']['circulation_density'];g=zero['values']['circulation']
    fp=slopes['circulation_density'];gp=slopes['circulation']
    fI=zero['intervals']['circulation_density'];gI=zero['intervals']['circulation']
    numerator=subtract(multiply(fp['interval'],gI),multiply(fI,gp['interval']))
    conditional=dict(value=(fp['value']*g-f*gp['value'])/g**2 if g else None,
        interval=divide(numerator,[gI[0]**2,gI[1]**2]))
    return dict(paths=paths,susceptibilities=slopes,
                density_response_inside_final_circulation=conditional,
                zero_source_density_inside_final_circulation=f/g if g else None)


def evolve_paths(X,ewait,source,H,B,O,lag,save_waves=True):
    """X is an optionally subnormalized vector or a purification matrix."""
    if X.ndim==1:X=X[:,None]
    E=Evolution(H,B,0.,4);paths={};waves={}
    for label,phase in [('zero',0.),('plus',np.pi/2),('minus',-np.pi/2)]:
        kick=np.exp(-1j*phase*source)
        initial=kick[:,None]*X
        Z,de=propagate(E,initial,lag)
        e=float(ewait+np.linalg.norm(de))
        p=np.sum(abs(Z)**2,axis=1)
        values={k:float(p@d) for k,d in O.items()}
        intervals={k:positive_interval(v,e) for k,v in values.items()}
        paths[label]=dict(values=values,intervals=intervals,full_path_error_upper=e,
            projected_norm=float(p.sum()),phase=phase,
            projected_source_work=float(np.vdot(initial,H@initial).real-np.vdot(X,H@X).real))
        if save_waves:waves[label]=Z
    return summarize_paths(paths),waves


def interference(full,a,b,wa,wb,O):
    paths={};complex_cross={}
    for phase in ['zero','plus','minus']:
        vals={k:a['paths'][phase]['values'][k]+b['paths'][phase]['values'][k] for k in O}
        ints={k:[max(0.,a['paths'][phase]['intervals'][k][0]+b['paths'][phase]['intervals'][k][0]),
                 min(1.,a['paths'][phase]['intervals'][k][1]+b['paths'][phase]['intervals'][k][1])] for k in O}
        paths[phase]=dict(values=vals,intervals=ints)
        ea=a['paths'][phase]['full_path_error_upper'];eb=b['paths'][phase]['full_path_error_upper']
        cross={}
        for k,d in O.items():
            za=wa[phase];zb=wb[phase]
            z=np.vdot(za,d[:,None]*zb)
            error=float(ea*np.linalg.norm(d[:,None]*zb)+eb*np.linalg.norm(d[:,None]*za)+ea*eb)
            difference=full['paths'][phase]['values'][k]-vals[k]
            assert abs(2*z.real-difference)<1e-10
            cross[k]=dict(real=float(z.real),imaginary=float(z.imag),absolute_error_upper=error,
                coherent_minus_pinched=float(2*z.real),
                coherent_minus_pinched_interval=[float(2*z.real-2*error),float(2*z.real+2*error)])
        complex_cross[phase]=cross
    pinched=summarize_paths(paths)
    slopes={}
    for k in O:
        a=complex_cross['plus'][k];b=complex_cross['minus'][k]
        I=subtract(a['coherent_minus_pinched_interval'],b['coherent_minus_pinched_interval'])
        v=(a['coherent_minus_pinched']-b['coherent_minus_pinched'])/2
        assert abs(v-full['susceptibilities'][k]['value']+pinched['susceptibilities'][k]['value'])<1e-10
        slopes[k]=dict(value=v,interval=[I[0]/2,I[1]/2])
    return dict(nonselective_record=pinched,offdiagonal_history_terms=complex_cross,
                response_interference=slopes)


def main_case(n,q,budget,wait,lag,snapshot_dir,records=True):
    start=time.monotonic();full=model(n,q,budget)
    source_model=model(n,q,budget,source=True)
    J,checks=inclusion(full,source_model)
    sub,H,B,states,sizes,psi=source_model
    O,events=detectors(n,states);S=O['source_density']
    snap=Path(snapshot_dir)/f'N{n}_Q{q}_B{budget}_T{wait:g}.npz'
    v=np.load(snap)['wave'];X=(J@v)[:,None]
    old=json.loads((REL/f'RELATIONAL_N{n}_Q{q}_B{budget}.json').read_text())
    row=min(old['rows'],key=lambda r:abs(r['T']-wait));assert abs(row['T']-wait)<1e-10
    ew=row['state_error_upper']
    result,waves=evolve_paths(X,ew,S,H,B,O,lag)
    record_rows=[]
    if records:
        hop_norm=np.sqrt(2*n*q)
        for label,event in events.items():
            A=event[:,None]*X;Z=(1-event)[:,None]*X
            p=float(np.linalg.norm(A)**2)
            a,wa=evolve_paths(A,ew,S,H,B,O,lag)
            b,wb=evolve_paths(Z,ew,S,H,B,O,lag)
            tests=interference(result,a,b,wa,wb,O)
            pI=positive_interval(p,ew)
            work=float(-2*np.vdot(A,H@Z).real)
            selection={k:dict(value=z['value']/p if p else None,
                interval=divide(z['interval'],pI)) for k,z in a['susceptibilities'].items()}
            record_rows.append(dict(record=label,initial_probability=p,initial_probability_interval=pI,
                selected_branch=a,complement_branch=b,selected_branch_normalized_susceptibilities=selection,
                nonselective_recording_work=work,
                nonselective_recording_work_error_upper=float(2*hop_norm*(2*ew+ew**2)),
                measured_system_pinching_entropy_bits=h2(p),
                trace_distance_to_pinched_state=float(np.sqrt(p*(1-p))),**tests))
    slopes=result['susceptibilities']
    ward=slopes['source_density']['value']+(n-1)*slopes['density']['value']-(n-1)*slopes['incident_edge_holes']['value']-((n-1)*(n-2)/2)*slopes['other_edge_holes']['value']
    assert abs(ward)<1e-10
    out=dict(N=n,Q=q,budget=budget,wait=wait,lag=lag,dimension=len(states),
        invariant_dimension=len(full[3]),embedding_checks=checks,initial_wait_error_upper=ew,
        snapshot_sha256=hashlib.sha256(snap.read_bytes()).hexdigest(),
        source_vertex=0,detector_vertices=list(range(1,n)),
        no_initial_filter=result,record_comparisons=record_rows,
        exact_total_fermion_response_ward_residual=float(ward),
        full_initial_energy=q/2,full_initial_variance=2*q,
        source_kick_work_absolute_error_upper=float(2*np.sqrt(2*n*q)*(2*ew+ew**2)),
        main_preparation_entropy_bits=0,seconds=time.monotonic()-start,
        scope='The initial evolved main state is unfiltered. The sole physical intervention on its response path is exp(-i eta n0). Final circulation conditioning is readout postselection, not an initial selection or autonomous geometry phase. Circulation event means at least one directed incidence cycle avoiding source zero. At N4 it is the fixed triangle (1,2,3). Pinched/selected paths are explicitly counterfactual recording diagnostics. Density is graph-distance independent. Source phase +/-pi/2 gives an exact derivative. Number conservation includes edge modes, so vertex-density response is not by itself transfer of a tagged particle. Analytic electric-leakage bounds in ordinary floating point, without formal roundoff intervals.')
    path=P/f'UNFILTERED_N{n}_Q{q}_B{budget}_W{wait:g}_L{lag:g}.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(saved=path.name,dimension=len(states),chi=slopes['density'],
        conditional=result['density_response_inside_final_circulation'],seconds=out['seconds'])),flush=True)
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=4);ap.add_argument('--q',type=int)
    ap.add_argument('--budget',type=int,default=256);ap.add_argument('--wait',type=float,default=10.)
    ap.add_argument('--lag',type=float,default=.5);ap.add_argument('--snapshot-dir',required=True)
    ap.add_argument('--no-records',action='store_true');ap.add_argument('--prepare-only',action='store_true')
    a=ap.parse_args();q=a.q or a.n
    if a.prepare_only:model(a.n,q,a.budget,source=True)
    else:main_case(a.n,q,a.budget,a.wait,a.lag,a.snapshot_dir,not a.no_records)
