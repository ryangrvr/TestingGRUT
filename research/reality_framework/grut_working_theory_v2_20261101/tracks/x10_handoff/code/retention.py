#!/usr/bin/env python3
"""Carry bath directions selected by their actual energy at a movable boundary."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from handoff import (BOUNDARIES, CASES, OBSERVABLES, adjoint, boundary_spread, comparison, evolve,
                     ledger, make_box, measure, pinched_frame, project, reaction_basis)

CASES_TESTED=tuple(CASES)


def adaptive(slab,state,extra,return_layout=False):
    D=slab.dim
    n=slab.geometry.bath_bins
    basis=np.broadcast_to(np.eye(D),(slab.l,D,D)).copy()
    terms=slab.l*slab.n*2*2*2
    limit=0.5e-6/terms
    rho=measure(slab,state,slab.observe(state,source=True))['rho']
    groups=[]
    orders=[]
    layout=np.zeros((slab.l,slab.n,2,2,2),int)
    bound=0.
    for ell in range(slab.l):
        active=list(slab.sterile_indices)
        for i in range(slab.n):
            for h in (0,1):
                for b in (0,1):
                    begin=(2*i+h)*slab.block+b*slab.nb+1
                    for sign in (0,1):
                        first=begin+sign*n
                        sec=slice(first,first+n)
                        f=slab.frequency[ell,i,h,sign*n:(sign+1)*n]
                        g=slab.g[ell,i,h,b,sign*n:(sign+1)*n]
                        frame=state.frame[ell,sec,:]

                        def assess(order):
                            u=reaction_basis(f,g,order)
                            retained,tail=u[:,:order],u[:,order:]
                            mixed=tail.T @ frame @ adjoint(retained.T @ frame)
                            hcross=(retained.T*f) @ tail
                            gap=abs(np.einsum('ij,ji->',hcross,mixed).real)
                            return u,float(slab.box.weight[ell]*gap/(state.a**4*rho))

                        for order in range(1,n+1):
                            u,epsilon=assess(order)
                            if epsilon<limit:
                                break
                        if extra and order<n:
                            order=min(n,order+int(extra))
                            u,epsilon=assess(order)
                            while epsilon>=limit and order<n:
                                order+=1
                                u,epsilon=assess(order)
                        basis[ell,sec,sec]=u
                        active.extend(range(first,first+order))
                        orders.append(order)
                        layout[ell,i,h,b,sign]=order
                        bound+=epsilon
        groups.append([np.asarray(sorted(active))])
    new,clip,cost=pinched_frame(state.frame,basis,groups)
    result=(replace(state,frame=new),dict(pauli_roundoff_clip=clip,stored_covariance_fraction=cost,
                                         normalized_energy_jump_bound=bound,
                                         orders_min=min(orders),orders_max=max(orders),
                                         orders_mean=float(np.mean(orders))))
    return (*result,basis,layout) if return_layout else result


def derivative_checks(slab,state):
    h=slab.full_hamiltonian(state.a,state.phi)
    S=np.zeros_like(h)
    for i in range(slab.n):
        for helicity in (0,1):
            start=(2*i+helicity)*slab.block
            idx=np.array([start,start+slab.nb])
            k=slab.k_signed[:,i,helicity]/state.a
            mass=slab.cfg.mass
            e=np.hypot(k,mass)
            S[:,idx[0],idx[0]]=(1+k/e)/2
            S[:,idx[1],idx[1]]=(1-k/e)/2
            S[:,idx[0],idx[1]]=S[:,idx[1],idx[0]]=mass/(2*e)
    C=np.zeros(slab.dim)
    view=C.reshape(slab.n,2,2,slab.nb)
    view[:,:,0,1:1+slab.geometry.bath_bins]=1.
    channel=np.broadcast_to(np.diag(C),h.shape).astype(complex)
    result={}
    for r in (0,1,2,4,6):
        out,_,_=project(slab,state,f'r{r}')
        delta=out.frame @ adjoint(out.frame)-state.frame @ adjoint(state.frame)
        checks={}
        for name,original in (('sterile',S),('one_daughter_channel',channel)):
            op=original.copy()
            errors=[]
            for _ in range(r+2):
                scale=max(1.,float(np.max(np.linalg.norm(op,axis=(-1,-2)))))
                error=np.einsum('lij,lji->l',delta,op).real
                errors.append(float(np.max(abs(error))/scale))
                op=1j*(h @ op-op @ h)
            checks[name]=errors
        result[str(r)]=checks
    return result


def run(args):
    name,folder=args
    slab=make_box(CASES[name])
    dt=CASES[name][-1]
    base,checkpoints=evolve(slab,slab.prepare(dt),0.,dt,checkpoints=BOUNDARIES)
    tests={}
    for t,state in sorted(checkpoints.items()):
        row={}
        for extra in (False,True):
            label='energy_plus_one' if extra else 'energy'
            prepared,info=adaptive(slab,state,extra)
            trace,_=evolve(slab,prepared,t,dt)
            row[label]=dict(**info,**comparison(base,trace),**ledger(trace),
                            terminal={k:float(trace[k][-1]) for k in OBSERVABLES})
        if name=='bath_spacing':
            for k in ('r6','r7'):
                prepared,clip,cost=project(slab,state,k)
                trace,_=evolve(slab,prepared,t,dt)
                row[k]=dict(pauli_roundoff_clip=clip,stored_covariance_fraction=cost,
                            **comparison(base,trace),**ledger(trace))
        tests[str(t)]=row
        print(name,t,{k:(v['jump_pass'],v['terminal_pass']) for k,v in row.items()},flush=True)
    result=dict(boundaries=tests,
                boundary_spreads={key:boundary_spread(tests,key) for key in ('energy','energy_plus_one')},
                passes_all_boundaries={key:all(v[key]['jump_pass'] and v[key]['terminal_pass'] and
                                              v[key]['energy_pressure_error']<1e-5 for v in tests.values())
                                       for key in next(iter(tests.values()))})
    if name=='reference':
        result['derivative_hierarchy']=derivative_checks(slab,checkpoints[2.25])
    (folder/f'{name}.json').write_text(json.dumps(result,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    with ProcessPoolExecutor(max_workers=3) as pool:
        list(pool.map(run,[(name,args.work) for name in CASES_TESTED]))
    result=dict(scope='state retention using actual omitted bath energy and no occupation adjustment',
                cases={k:json.loads((args.work/f'{k}.json').read_text()) for k in CASES_TESTED},
                scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
