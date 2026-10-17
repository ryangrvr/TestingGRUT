"""One bounded A1 cycle: finite controls, not large-N proof validation."""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import resource
import time
import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal, expm
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from model import (JP,JM,PathModel,exact_cell_matrix,geometry_bound,neutral_constants)

HERE=Path(__file__).resolve().parent


def expectation(v,h): return float(np.vdot(v,h@v).real)


def cell_controls():
    mod=PathModel(2,matching=((0,1),),resource=1)
    order=[((),(),1),((0,),(1,),0),((1,),(0,),0),((0,),(0,),0),((1,),(1,),0)]
    ix=[mod.index[x] for x in order]
    manual=exact_cell_matrix()
    entry=float(np.max(abs(mod.H[ix,:][:,ix].toarray()-manual)))
    assert entry<1e-12
    z=2.3+.17j; a,b=JP+JM,JP-JM
    gd=(z-2)/2*(1/((z-13/4)*(z-2)-a*a)+1/((z-13/4)*(z-2)-b*b))
    gf=1/(z-13/4-gd)
    error=abs(np.linalg.inv(z*np.eye(5)-manual)[0,0]-gf)
    assert error<1e-12
    v=mod.initial(); results=[]
    for t in [math.pi/2,math.pi]:
        u=mod.evolve(v,t)
        prob=float(np.vdot(u,(mod.particles>=1)*u).real)
        results.append({'time':t,'matter_probability':prob,
                        'energy':expectation(u,mod.H)})
    assert results[0]['matter_probability']>.99
    assert results[1]['matter_probability']<.01
    assert all(abs(x['energy']-13/4)<1e-12 for x in results)
    return {'dimension':5,'independent_matrix_entry_error':entry,
            'fuel_green_function_error':float(error),'observations':results}


def coupled_cells():
    mod=PathModel(8,matching=((1,2),(5,6)),resource=2)
    v=mod.initial(); e0=expectation(v,mod.H)
    variance=expectation(v,mod.H@mod.H)-e0*e0
    assert abs(e0-6.5)<1e-12 and abs(variance-2)<1e-12
    eta=5/48
    lower=1-eta*eta*math.pi*math.pi/2
    upper=2*eta*eta*math.pi*math.pi
    results=[]
    for t in [math.pi/2,math.pi]:
        u=mod.evolve(v,t)
        p=float(np.vdot(u,(mod.particles>=2)*u).real)
        results.append({'time':t,'P_at_least_K_probability':p,
                        'mean_particle_number':float(np.vdot(u,mod.particles*u).real),
                        'energy_error':abs(expectation(u,mod.H)-e0)})
    assert results[0]['P_at_least_K_probability']>=lower-1e-10
    assert results[1]['P_at_least_K_probability']<=upper+1e-10
    assert max(x['energy_error'] for x in results)<1e-10
    # Independent six-site ordinary hard-core hopping norm, all sectors.
    h=np.zeros((64,64))
    for col in range(64):
        for x in range(5):
            if ((col>>x)&1)!=((col>>(x+1))&1):
                h[col^(1<<x)^(1<<(x+1)),col]=-1
    norm=float(max(abs(np.linalg.eigvalsh(h))))
    assert norm<math.sqrt(15)<4
    return {'events':8,'cells':2,'dimension':mod.dim,
            'all_physical_states_in_exact_resource_sector':True,
            'initial_energy':e0,'initial_energy_variance':variance,
            'six_site_unit_hopping_norm':norm,
            'many_cell_formation_lower':lower,'many_cell_return_upper':upper,
            'observations':results}


def neutral_controls():
    c=neutral_constants()
    assert c['ground_dispersion_error']<8.1e-6
    assert c['excited_dispersion_error']<4.2e-5
    assert c['collision_diagonal_difference_lower']>.444
    r=np.arange(-32,33); diag=2+abs(r)+.25*(abs(r)==1)
    rows=[]; e0max=0.; e1max=0.
    for k in np.linspace(0,math.pi,17):
        alpha2=JP*JP+JM*JM+2*JP*JM*math.cos(k)
        vals,vec=eigh_tridiagonal(diag,np.full(64,-math.sqrt(alpha2)),select='i',select_range=(0,2))
        parity=np.array([float(np.vdot(vec[:,i],vec[::-1,i]).real) for i in range(3)])
        assert parity[0]>.99999 and parity[1]<-.99999 and parity[2]>.99999
        pred=[2-2*alpha2/(5/4),13/4-4*alpha2/3,13/4+4*alpha2/15]
        errors=abs(vals-np.array(pred))
        e0max=max(e0max,float(errors[0]));e1max=max(e1max,float(max(errors[1:])))
        assert errors[0]<=c['ground_dispersion_error']+1e-11
        assert max(errors[1:])<=c['excited_dispersion_error']+1e-11
        rows.append({'k':float(k),'finite_R32_energies':vals.tolist(),'second_order':pred})
    # Exact finite-path one-pair matrix; all separations, no cutoff.
    m=PathModel(12,resource=1)
    vals,vec=eigh(m.H.toarray())
    pp=np.array([i for i,(p,n,f) in enumerate(m.states) if p==n])
    qq=np.array([i for i in range(m.dim) if i not in set(pp)])
    h=m.H.toarray(); pvq=h[np.ix_(pp,qq)]
    sigma=-pvq@np.linalg.solve(h[np.ix_(qq,qq)]-2*np.eye(len(qq)),pvq.conj().T)
    leading=np.zeros((12,12)); aa=13/46080; tt=1/3840
    for x in range(12): leading[x,x]=-aa*(1 if x in (0,11) else 2)
    for x in range(11): leading[x,x+1]=leading[x+1,x]=-tt
    remainder=float(max(abs(np.linalg.eigvalsh(sigma-leading))))
    assert remainder<=c['ground_remainder_at_E2']+1e-12
    low=vec[:,:12]
    low_q_norm=float(np.linalg.svd(low[qq,:],compute_uv=False)[0])
    assert low_q_norm<=1/23+1e-12
    adjacent=m.index[((5,),(6,),0)]
    low_weight=float(np.sum(abs(low[adjacent,:])**2))
    assert low_weight<=1/529+1e-12
    # Complete two-pair physical field space. Self-energy is NOT exponentiated.
    two=PathModel(8,resource=2)
    cp=np.array([i for i,(p,n,f) in enumerate(two.states) if p==n])
    cq=np.array([i for i in range(two.dim) if i not in set(cp)])
    hh=two.H.toarray(); v=hh[np.ix_(cp,cq)]
    se=v@np.linalg.solve(4.5*np.eye(len(cq))-hh[np.ix_(cq,cq)],v.conj().T)
    se_norm=float(max(abs(np.linalg.eigvalsh(se))))
    assert se_norm<=c['two_pair_self_energy_upper']+1e-12
    eff=hh[np.ix_(cp,cp)]+se
    a=list(cp).index(two.index[((2,3),(2,3),0)])
    b=list(cp).index(two.index[((1,6),(1,6),0)])
    difference=float((eff[a,a]-eff[b,b]).real)
    assert difference>=c['collision_diagonal_difference_lower']-1e-12
    return {'constants':c,'relative_cutoff':32,'finite_relative_blocks':rows,
            'max_ground_dispersion_error_in_controls':e0max,
            'max_excited_dispersion_error_in_controls':e1max,
            'full_one_pair_path':{'events':12,'dimension':m.dim,
                'complete_Schur_remainder':remainder,'low_band_separated_operator_norm':low_q_norm,
                'adjacent_preparation_low_band_weight':low_weight},
            'two_pair_collision':{'events':8,'dimension':two.dim,
                'self_energy_norm':se_norm,'near_minus_far_diagonal':difference},
            'scattering_matrix_derived':False,'finite_density_gas_derived':False}


def geometry_control():
    source=HERE.parent/'grut_p1_preparation_20261017'/'model.py'
    spec=importlib.util.spec_from_file_location('p1_frozen',source)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    mod=old.Model(4,matching=((1,2),),modes=2,
                  delta=.23+np.linspace(-.25,.25,4),fuel_delta=[.25],lam=[1.])
    initial0=mod.initial(); init=np.zeros(mod.dim,complex)
    for i in np.flatnonzero(initial0):
        w,q,r,p,m,f=mod.states[i]
        init[mod.index[(w,q,r,p,m,1)]]=initial0[i]
    h=mod.Hop+mod.C
    resource_count=np.array([len(p)+f.bit_count() for w,q,r,p,m,f in mod.states])
    selected=np.flatnonzero(resource_count==1)
    outside=np.flatnonzero(resource_count!=1)
    assert h[outside,:][:,selected].nnz==0
    ref=mod.reference(h,(0,1,2))
    hs=h[selected,:][:,selected]; rs=ref[selected,:][:,selected]; psi=init[selected]
    t=.02
    actual=expm_multiply((-1j*t)*hs,psi); reference=expm_multiply((-1j*t)*rs,psi)
    error=float(np.linalg.norm(actual-reference))
    upper=geometry_bound(4,t,a=1,L=-1)
    assert error<=upper+1e-10
    prefix=mod.retained((0,1,2))[selected]
    loss=float(np.vdot(actual,(~prefix)*actual).real)
    inconsistent=np.array([q!=sum(1<<x for x in w) for w,q,r,p,m,f in mod.states])[selected]
    storage=float(np.vdot(actual,inconsistent*actual).real)
    assert storage>0 and storage<=upper*upper+1e-10
    vac=np.flatnonzero(resource_count==0)
    difference=h[vac,:][:,vac]-mod.H0[vac,:][:,vac]
    entry=float(np.max(abs(difference.data),initial=0))
    assert entry==0
    return {'events':4,'receiver_modes':2,'physical_receiver_modes_would_be':1024,
            'full_builder_dimension':mod.dim,'exact_resource_dimension':len(selected),
            'time':t,'state_error':error,'state_bound':upper,
            'prefix_loss':loss,'storage_inconsistency':storage,
            'vacuum_matrix_difference':entry,'conversion_present':True,
            'complete_receiver_and_geometry_states_in_resource_sector':True}


def run():
    start=time.perf_counter()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for p,d in pins['sources'].items():
        assert hashlib.sha256((HERE.parents[2]/p).read_bytes()).hexdigest()==d
    result={'arithmetic':'ordinary double precision; no outward rounding or independent proof review',
            'cell':cell_controls(),'coupled_cells':coupled_cells(),
            'neutral':neutral_controls(),'moving_geometry':geometry_control()}
    for p,d in pins['sources'].items():
        assert hashlib.sha256((HERE.parents[2]/p).read_bytes()).hexdigest()==d
    eta=5/48
    result['analytic_formula_evaluation_not_simulation']={
        'formation_lower':1-eta*eta*math.pi*math.pi/2,
        'return_upper':2*eta*eta*math.pi*math.pi,
        'moving_geometry_error_at_N1024_t_pi':geometry_bound(1024,math.pi),
        'primary_min_N':2**18,'primary_use_window_contains_pi':True,
        'scalable_retention_gate':'refuted at permitted resonant point'}
    result['source_pins_unchanged']=len(pins['sources'])
    result['resources']={'wall_seconds':time.perf_counter()-start,
        'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        'written_bytes':sum(p.stat().st_size for p in HERE.iterdir() if p.is_file())}
    assert result['resources']['wall_seconds']<60
    assert result['resources']['peak_rss_MiB']<1024
    assert result['resources']['written_bytes']<100*1024**2
    return result


if __name__=='__main__':
    r=run()
    (HERE/'CHECKS.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k!='neutral'},indent=2))
    print(json.dumps({'neutral_controls':{k:v for k,v in r['neutral'].items()
                                         if k!='finite_relative_blocks'}},indent=2))
