"""Bounded S6 tests, original root weights, complete stop disclosures."""
import hashlib
import json
from math import factorial,sqrt
from pathlib import Path
import sys
import numpy as np
from scipy.fft import dst
from scipy.sparse import diags

import conditioned as C
import controls
from fock import space

BASE=Path(__file__).resolve().parent
S5DIR=BASE.parent/'wm2_history_formation_20261016'


def save(path,obj):
    temporary=Path(str(path)+'.partial')
    temporary.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def physical_reference():
    e=C.Evolution([1,1],[-.25,.25],576)
    r=e.run(36,max_seconds=180)
    # The complete reference is at the actual completed endpoint. If a budget
    # stops this solver, it cannot establish a T=36 formation result.
    t=r['completed_steps_time']
    exact=controls.M.vacuum(3,576,[0,-.25,.25],.9,classes=True)
    y=controls.M.evolve(exact['H'],controls.M.unit(len(exact['states']),exact['seed']),t)
    projected={}
    for u in e.words:
        k=len(u); w=dst(e.frames[u],type=1,norm='ortho',axis=0)
        if k==0:
            projected[u]=np.array([y[exact['index'][u,()]]])
        elif k==1:
            field=np.array([y[exact['index'][u,(a,)]] for a in range(576)])
            projected[u]=w.conj().T@field
        else:
            tensor=np.zeros((576,576),complex)
            for a in range(576):
                for b in range(a,576):
                    z=y[exact['index'][u,(a,b)]]/(1 if a==b else sqrt(2))
                    tensor[a,b]=tensor[b,a]=z
            reduced=w.conj().T@tensor@w.conj()
            projected[u]=np.array([reduced[a,b]*(1 if a==b else sqrt(2)) for a,b in space(w.shape[1],2).occ])
    missing=max(0.0,float(np.vdot(y,y).real)-C.norm2(projected))
    error=sqrt(sum(float(np.linalg.norm(e.x[u]-projected[u])**2) for u in e.words)+missing)
    assert error<=r['analytic_state_error_account']+3e-9
    return dict(trajectory=r,complete_physical_dimension=len(y),reference_time=t,
                full_wavefunction_error=error,complete_reference_norm_error=abs(float(np.vdot(y,y).real)-1),
                target_time_reached=(r['status']=='COMPLETE'),comparison_uses_entire_two_photon_tensor=True)


def uniform_control():
    n,modes,t=4,1024,64
    e=C.Evolution([3],[0.0],modes)
    r=e.run(t,max_seconds=180)
    actual=r['completed_steps_time']
    off=np.r_[-.9/n,-np.ones(modes-1)]
    h=diags([off,np.zeros(modes+1),off],[-1,0,1],format='csr')
    z=controls.M.evolve(h,controls.M.unit(modes+1,0),actual)
    spectral=dst(z[1:],type=1,norm='ortho')
    projected={}
    for u in e.words:
        k=len(u); field=e.frames[u].conj().T@spectral
        entries=[]
        for o in space(e.frames[u].shape[1],k).occ:
            a=sqrt(factorial(n-1)/factorial(n-1-k))*z[0]**(n-1-k)
            for j in o: a*=field[j]
            for j in set(o): a/=sqrt(factorial(o.count(j)))
            entries.append(a)
        projected[u]=np.array(entries)
    missing=max(0.0,1-C.norm2(projected))
    error=sqrt(sum(float(np.linalg.norm(e.x[u]-projected[u])**2) for u in e.words)+missing)
    assert error<=r['analytic_state_error_account']+3e-9
    return dict(trajectory=r,quadratic_full_state_error=error,reference_time=actual,
                exact_pump_survival=float(abs(z[0])**2),target_time_reached=r['status']=='COMPLETE')


def families():
    return [dict(N=n,bar_delta=bar,values=[bar-.25,bar+.25],
                 global_counts=[(n+1)//2,n//2],g=.9)
            for n in (4,5,6) for bar in (0.0,.5)]


def component(family,root):
    counts=family['global_counts'].copy(); counts[root]-=1
    values=family['values'].copy()
    return dict(counts=counts,values=values,g=family['g'],root_type=root,
                weight=family['global_counts'][root]/family['N'])


def key(counts,values,g):
    return json.dumps([list(counts),list(values),g])


def mirror(counts,values,g):
    pairs=sorted((-v,n) for n,v in zip(counts,values))
    return key([n for _,n in pairs],[v for v,_ in pairs],g)


def run(skip_validation=False):
    original_sha=hashlib.sha256((S5DIR/'manifest.json').read_bytes()).hexdigest()
    plan=dict(parent='98f1e861c3509c67a06032b7742c295069a9a273',
       charter_commit='e627bb5f8287efb47d2ab4c79efbad6fd2984132',physical_law_changed=False,
       observation='4*N^2',families=families(),rank_ceiling=10,coefficient_ceiling=2_000_000,
       cached_overlap_entry_ceiling=8_000_000,per_instance_step_boundary_seconds=180,
       strang_budget=.15,projection_budget=.20,residual_budget=.20,total_error_envelope=.60,
       polynomial_degree=8,residual_terms=3,dense_overlap_entry_threshold=65536,
       original_root_weights=True,roundoff_certified=False)
    save(BASE/'RUN_PLAN.json',plan)
    result=dict(plan=plan,checks={},attempts=[],families=[],all_size_formation_proved=False,
                independent_external_validation=False,S5_manifest_sha256=original_sha)
    if not skip_validation:
        for name in ('fock_checks','cross_history_checks','emission_step_checks','dynamics_checks'):
            result['checks'][name]=getattr(controls,name)()
            print('CHECK',name,'PASS',flush=True)
        result['checks']['physical_receiver']=physical_reference()
        print('PHYSICAL',result['checks']['physical_receiver']['trajectory']['status'],flush=True)
        save(BASE/'RESULTS.json',result)
        result['checks']['uniform_receiver']=uniform_control()
        print('UNIFORM',result['checks']['uniform_receiver']['trajectory']['status'],flush=True)
        save(BASE/'RESULTS.json',result)
    seen={}
    for family in families():
        row=dict(**family,components=[])
        for root in (0,1):
            c=component(family,root)
            direct=key(c['counts'],c['values'],c['g']); mirrored=mirror(c['counts'],c['values'],c['g'])
            if direct in seen:
                index=seen[direct]; relation='identical component'
            elif mirrored in seen:
                index=seen[mirrored]; relation='S5 proved sign identity'
            else:
                index=len(result['attempts']); relation='new trajectory'
                print('START',index,family['N'],family['bar_delta'],c['counts'],flush=True)
                e=C.Evolution(c['counts'],c['values'],64*family['N']**2)
                r=e.run(4*family['N']**2,progress=lambda v:print('PROGRESS',index,json.dumps(v),flush=True))
                result['attempts'].append(r); seen[direct]=index
                print('FINISH',index,r['status'],r['stop_reason'],r['completed_steps_time'],flush=True)
                save(BASE/'RESULTS.json',result)
            row['components'].append(dict(**c,index=index,relation=relation))
        lower=full=missing=0.0
        for c in row['components']:
            r=result['attempts'][c['index']]
            if r['status']=='COMPLETE':
                lower+=c['weight']*r['formation_lower']; full+=c['weight']*r['full_assembly_lower']
            else:
                missing+=c['weight']
        row.update(status='COMPLETE' if missing==0 else 'UNRESOLVED',formation_lower=lower,
                   full_assembly_lower=full,unresolved_root_weight=missing,
                   gate_pass=bool(missing==0 and lower>=.5 and
                                  all(result['attempts'][c['index']]['gate_pass'] for c in row['components'])))
        result['families'].append(row); save(BASE/'RESULTS.json',result)
    assert original_sha==hashlib.sha256((S5DIR/'manifest.json').read_bytes()).hexdigest()
    save(BASE/'RESULTS.json',result)
    return result


if __name__=='__main__':
    run()
