#!/usr/bin/env python3
"""Move the handoff boundary and measure the state information needed downstream."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from coupled import Configuration, CoupledSlab, Geometry, State, adjoint

T=4.5
BOUNDARIES=(1.25,2.25,3.25)
CASES={
    'reference':(1,2,8,5.5,0.01),
    'half_step':(1,2,8,5.5,0.005),
    'bath_spacing':(1,2,12,5.5,0.01),
    'bath_cutoff':(1,2,11,7.2,0.01),
    'longitudinal':(2,2,8,5.5,0.01),
    'transverse':(1,3,8,5.5,0.01),
}
NAMES=('bare','dressed','separated','dressed_half','r0','r1','r2','r4','r6','r7')
OBSERVABLES=('sterile_count','bath_energy','daughter_plus','daughter_minus','daughter_total','daughter_difference')


def make_box(settings):
    longitudinal,radial,bins,cutoff,_=settings
    cfg=Configuration(longitudinal=longitudinal,transverse_cut=1.6,transverse_bins=radial,
                      points=25,scalar_modes=2)
    return CoupledSlab(Geometry(cfg,bath_bins=bins,bath_cutoff=cutoff))


def reaction_basis(frequency,coupling,order):
    """Keep the Krylov directions. Diagonalize the discarded bath energy."""
    length=len(frequency)
    if order==0:
        return np.eye(length)
    v=coupling/np.linalg.norm(coupling)
    columns=[v]
    for _ in range(1,order):
        candidate=frequency*columns[-1]
        u=np.column_stack(columns)
        for _ in range(2):
            candidate-=u @ (u.T @ candidate)
        candidate/=np.linalg.norm(candidate)
        columns.append(candidate)
    q,_=np.linalg.qr(np.column_stack(columns),mode='complete')
    if order<length:
        rest=q[:,order:]
        _,axes=np.linalg.eigh((rest.T*frequency) @ rest)
        q[:,order:]=rest @ axes
    return q


def reaction_partition(slab,order):
    n=slab.geometry.bath_bins
    basis=np.broadcast_to(np.eye(slab.dim),(slab.l,slab.dim,slab.dim)).copy()
    retain=list(slab.sterile_indices)
    for i in range(slab.n):
        for h in (0,1):
            for b in (0,1):
                start=(2*i+h)*slab.block+b*slab.nb+1
                for sign in (0,1):
                    sl=slice(start+sign*n,start+(sign+1)*n)
                    target=slice(sign*n,(sign+1)*n)
                    for ell in range(slab.l):
                        basis[ell,sl,sl]=reaction_basis(slab.frequency[ell,i,h,target],
                                                       slab.g[ell,i,h,b,target],order)
                    retain.extend(range(sl.start,sl.start+order))
    return basis,np.asarray(sorted(retain))


def pinched_frame(frame,basis,groups):
    """A block pinching preserves positivity and the fermionic Pauli interval."""
    coords=adjoint(basis) @ frame
    diag=np.sum(abs(coords)**2,axis=-1)
    factors=np.zeros((frame.shape[0],frame.shape[1],frame.shape[1]),complex)
    idx=np.arange(frame.shape[1])
    factors[:,idx,idx]=np.sqrt(np.clip(diag,0,1))
    clip=float(np.max(abs(diag-np.clip(diag,0,1))))
    for ell,blocks in enumerate(groups):
        for ids in blocks:
            keep=coords[ell,ids,:]
            vals,u=np.linalg.eigh(keep @ adjoint(keep))
            corrected=np.clip(vals,0,1)
            clip=max(clip,float(np.max(abs(vals-corrected))))
            factors[ell,ids[:,None],ids[None,:]]=u*np.sqrt(corrected)
    entries=[frame.shape[1]-sum(map(len,b))+sum(len(k)**2 for k in b) for b in groups]
    return basis @ factors,clip,float(np.mean(entries)/frame.shape[1]**2)


def project(slab,state,name):
    D=slab.dim
    if name in ('bare','dressed','dressed_half'):
        _,basis=np.linalg.eigh(slab.full_hamiltonian(state.a,state.phi,free=name=='bare'))
        groups=[[] for _ in range(slab.l)]
        if name=='dressed_half':
            score=np.sum(abs(basis[:,slab.sterile_indices,:])**2,axis=1)
            groups=[[np.sort(np.argsort(s)[-D//2:])] for s in score]
    elif name=='separated':
        basis=np.broadcast_to(np.eye(D),(slab.l,D,D))
        bath=np.setdiff1d(np.arange(D),slab.sterile_indices)
        groups=[[slab.sterile_indices,bath] for _ in range(slab.l)]
    else:
        basis,active=reaction_partition(slab,int(name.removeprefix('r')))
        groups=[[active] for _ in range(slab.l)]
    frame,clip,cost=pinched_frame(state.frame,basis,groups)
    return replace(state,frame=frame),clip,cost


def measure(slab,state,obs):
    d=slab.measure(state,obs,source=True)
    d['daughter_total']=d['daughter_plus']+d['daughter_minus']
    d['daughter_difference']=d['daughter_plus']-d['daughter_minus']
    return d


def evolve(slab,start,t0,dt,end=T,checkpoints=()):
    state=start
    obs=slab.observe(state,source=True)
    keep={}
    history=[]
    budget=0.
    steps=round((end-t0)/dt)
    for i in range(steps+1):
        t=t0+i*dt
        d=measure(slab,state,obs)
        history.append(dict(t=t,pressure_work=budget,**d))
        if any(abs(t-u)<1e-11 for u in checkpoints):
            keep[t]=state
        if i<steps:
            state,obs,amount=slab.step(state,obs,dt)
            budget+=amount
    trace={k:np.array([row[k] for row in history]) for k in history[0]}
    return trace,keep


def ledger(trace):
    vol=trace['a']**3
    E=vol*trace['rho']
    err=float(np.max(abs((E+trace['pressure_work'])/E[0]-1)))
    mom=float(np.max(abs(trace['momentum']-trace['momentum'][0]))/E[0])
    qwork=cumulative_trapezoid(vol*trace['work'],trace['t'],initial=0.)
    qstress=vol*trace['rho_q']-(vol*trace['rho_q'])[0]+cumulative_trapezoid(
        3*trace['H']*vol*trace['P_q'],trace['t'],initial=0.)
    return dict(energy_pressure_error=err,momentum_error=mom,
                reciprocal_quantum_work_error=float(np.max(abs(qwork-qstress))/E[0]),
                integrated_scalar_to_quantum_work=float(qwork[-1]))


def comparison(reference,current):
    at=int(np.searchsorted(reference['t'],current['t'][0]-1e-10))
    scale=reference['rho'][at]
    jumps={k:float((current[k][0]-reference[k][at])/scale)
           for k in ('rho','P','work')}
    jumps['momentum']=float((current['momentum'][0]-reference['momentum'][at])/
                            (scale*reference['a'][at]**3))
    daughter_floor=1e-3*max(abs(reference['daughter_total']))
    floors={key:(1e-3*max(abs(reference[key])) if key=='bath_energy' else daughter_floor)
            for key in OBSERVABLES}
    scales={key:max(abs(reference[key][-1]),floors[key],1e-15) for key in OBSERVABLES}
    errors={key:float(abs(current[key][-1]-reference[key][-1])/scales[key]) for key in OBSERVABLES}
    query=current['t']
    ref_density=np.interp(query,reference['t'],reference['rho'])
    stress={k:float(np.max(abs(current[k]-np.interp(query,reference['t'],reference[k]))/ref_density))
            for k in ('rho','P','work')}
    return dict(jump_over_reference_density=jumps,terminal_relative_errors=errors,
                terminal_scales=scales,downstream_normalized_stress_error=stress,
                jump_pass=bool(max(map(abs,jumps.values()))<1e-6),
                terminal_pass=bool(max(errors.values())<0.01))


def boundary_spread(rows,method):
    data=[v[method] for v in rows.values()]
    return {key:float(np.ptp([d['terminal'][key] for d in data])/data[0]['terminal_scales'][key])
            for key in OBSERVABLES}


def static_projection(slab):
    zero=np.zeros_like(slab.box.x)
    initial=State(1.,zero,zero,slab.vacuum())
    ref=measure(slab,initial,slab.observe(initial,source=True))
    result={}
    for name in NAMES:
        changed,clip,_=project(slab,initial,name)
        out=measure(slab,changed,slab.observe(changed,source=True))
        result[name]=dict(rho_change=out['rho']-ref['rho'],count_change=out['daughter_total']-
                          ref['daughter_total'],pauli_roundoff_clip=clip)
    return result


def run(args):
    name,folder,replay=args
    slab=make_box(CASES[name])
    dt=CASES[name][-1]
    start=slab.prepare(dt)
    data,points=evolve(slab,start,0.,dt,checkpoints=BOUNDARIES if replay else ())
    budget=ledger(data)
    np.savez_compressed(folder/f'{name}.npz',**data)
    base=dict(parameters=dict(longitudinal=slab.cfg.longitudinal,transverse=slab.cfg.transverse_bins,
                              cutoff=slab.geometry.bath_cutoff,bath_bins=slab.geometry.bath_bins,dt=dt),
              dimension=slab.dim,recurrence=slab.recurrence,ledger=budget,
              initial={k:float(v[0]) for k,v in data.items()},
              terminal={k:float(v[-1]) for k,v in data.items()},
              max_abs_work=float(max(abs(data['work']))))
    print(name,'reference budget',budget['energy_pressure_error'],'daughter counts',
          base['terminal']['daughter_plus'],base['terminal']['daughter_minus'],flush=True)
    if replay:
        out={}
        for t,state in sorted(points.items()):
            result={}
            for method in NAMES:
                changed,clip,cost=project(slab,state,method)
                trace,_=evolve(slab,changed,t,dt)
                result[method]=dict(**comparison(data,trace),**ledger(trace),
                                    pauli_roundoff_clip=clip,stored_covariance_fraction=cost,
                                    terminal={k:float(trace[k][-1]) for k in OBSERVABLES})
            out[str(t)]=result
            print(name,t,'minimum jump',min(max(map(abs,v['jump_over_reference_density'].values()))
                                             for v in result.values()),flush=True)
        base['handoffs']=out
        base['static_projection']=static_projection(slab)
        base['boundary_spreads']={name:boundary_spread(out,name) for name in NAMES}
        base['passes_all_boundaries']={
            method:all(rows[method]['jump_pass'] and rows[method]['terminal_pass']
                       and rows[method]['energy_pressure_error']<1e-5 for rows in out.values())
            for method in NAMES}
    (folder/f'{name}.json').write_text(json.dumps(base,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--only',nargs='+',choices=CASES)
    p.add_argument('--gather',action='store_true')
    args=p.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    names=args.only or CASES
    if not args.gather:
        jobs=[(name,args.work,name in ('reference','half_step')) for name in names]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(run,jobs))
        if args.only:
            return
    cases={name:json.loads((args.work/f'{name}.json').read_text()) for name in CASES}
    ref=cases['reference']['terminal']
    comparison_keys=('rho_q','bath_energy','daughter_plus','daughter_minus','daughter_total','daughter_difference')
    refinements={name:{k:float(case['terminal'][k]/ref[k]-1) if abs(ref[k])>1e-15 else None
                       for k in comparison_keys} for name,case in cases.items() if name!='reference'}
    summary=dict(scope='finite one-axis Gaussian vacuum bath and scalar, averaged FRW stress',
                 boundaries=BOUNDARIES,cases=cases,terminal_resolution_differences=refinements,
                 scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
