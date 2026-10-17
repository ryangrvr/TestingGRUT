"""Independent full physical-state controls for charge-sector MPS evolution."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
from math import sqrt
import json
import time
import numpy as np
from scipy.linalg import expm
import bounds as E
import gates as G
from mps import Receiver

ROOT=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('s8_frozen_s4',ROOT/'wm2_prepared_disorder_20261016/model.py')
S4=importlib.util.module_from_spec(spec);spec.loader.exec_module(S4)


def vector(e,model):
    return np.array([e.amplitude(word,occupation) for word,occupation in model['states']],complex)


def arithmetic():
    rng=np.random.default_rng(820171);rows=[]
    for n in [1,17,200]:
        x=rng.normal(size=n)+1j*rng.normal(size=n)
        true=sum((Q.from_float(float(z.real))**2+Q.from_float(float(z.imag))**2 for z in x),Q(0))
        lo,hi=E.norm2(x);assert Q.from_float(lo)<=true<=Q.from_float(hi)
        q=np.linalg.qr(rng.normal(size=(n+3,min(n,6)))+1j*rng.normal(size=(n+3,min(n,6))))[0]
        gram=float(np.linalg.norm(q.conj().T@q-np.eye(q.shape[1]),2))
        bound=E.gram(q);assert gram<=bound
        rows.append(dict(entries=n,norm_interval=[lo,hi],observed_Gram_error=gram,Gram_upper=bound))
    return rows


def field_and_commutators():
    rows=[]
    for R,dt,M in [(4,.03,2),(8,.1,3),(16,.2,5)]:
        h=np.diag(-np.ones(R-1),1)+np.diag(-np.ones(R-1),-1)
        exact=expm(-1j*dt*h);p=np.eye(R,dtype=complex)
        for a in [G.A,G.B,G.A]:
            _,_,t,_=G.field_gates(M,a*dt/2)
            for seq in [range(R-1),range(R-2,-1,-1)]:
                for i in seq:p[[i,i+1]]=t@p[[i,i+1]]
        observed=float(np.linalg.norm(exact-p,2));certificate=G.one_particle_audit(R,dt)
        assert observed<=certificate['one_particle_error_upper']
        rows.append(dict(R=R,dt=dt,one_particle_error=observed,certificate=certificate))
    # The complete four-mode sparse model independently tests the first-three
    # local-algebra/spectator maximum, including all modes in the reference.
    comm=[]
    for n,delta in [(3,[0,-.25,.25]),(4,[0,-.25,.25,-.25])]:
        model=S4.vacuum(n,4,delta,.9,classes=True)
        f=(model['H']-model['H_emission']).toarray()
        d=np.diag(np.array([sum(model['values'][a] for a in u) for u,o in model['states']]))
        f-=d;g=model['H_emission'].toarray()+d
        c=g@f-f@g;first=g@c-c@g;second=-f@c+c@f
        values=sorted(set(delta[1:]));counts=[delta[1:].count(v) for v in values]
        certified=G.commutator_certificate(counts,values,.9)
        observed=[float(np.linalg.norm(x,2)) for x in [first,second]]
        assert observed[0]<=certified['GGF_upper'] and observed[1]<=certified['FFG_upper']
        comm.append(dict(N=n,physical_modes=4,full_dimension=len(f),GGF=observed[0],FFG=observed[1],
            GGF_upper=certified['GGF_upper'],FFG_upper=certified['FFG_upper']))
    return dict(field=rows,complete_physical_commutators=comm)


def wavefunctions():
    rows=[]
    for n,modes,delta in [(3,5,[0,-.25,.25]),(4,4,[0,-.25,.25,-.25])]:
        values=sorted(set(delta[1:]));counts=[delta[1:].count(v) for v in values]
        model=S4.vacuum(n,modes,delta,.9,classes=True)
        for T in [.5,2,8]:
            e=Receiver(counts,values,modes);r=e.run(T,strang_budget=.001,cut=False,max_seconds=90)
            assert r['status']=='COMPLETE',r
            reference=S4.evolve(model['H'],S4.unit(len(model['states']),model['seed']),T)
            approximate=vector(e,model);error=float(np.linalg.norm(reference-approximate))
            assert error<=r['analytic_state_error_account']+2e-11
            macro=sum(abs(approximate[j])**2 for j,(u,o) in enumerate(model['states']) if len(u)+1>=max(2,(n+1)//2))
            norm=float(np.vdot(approximate,approximate).real)
            assert r['computed_normalized_event_intervals']['macro'][0]<=macro/norm<=r['computed_normalized_event_intervals']['macro'][1]
            rows.append(dict(N=n,modes=modes,T=T,physical_dimension=len(reference),
                full_wavefunction_error=error,complete_error_account=r['analytic_state_error_account'],
                error_components=r['error_components'],maximum_charge_rank=r['maximum_charge_rank'],
                maximum_bond_rank=r['maximum_bond_rank'],gates=r['two_site_gates']))
    return rows


def physical_reference():
    started=time.monotonic();n,L,T=3,576,36
    model=S4.vacuum(n,L,[0,-.25,.25],.9,classes=True)
    reference=S4.evolve(model['H'],S4.unit(len(model['states']),model['seed']),T)
    e=Receiver([1,1],[-.25,.25],L)
    r=e.run(T,max_seconds=180)
    assert r['status']=='COMPLETE',r
    approximate=vector(e,model)
    error=float(np.linalg.norm(reference-approximate))
    assert error<=r['analytic_state_error_account']+2e-9
    return dict(N=n,L=L,T=T,physical_dimension=len(reference),
        computational_sites=r['computational_receiver_sites'],
        full_wavefunction_error=error,complete_error_account=r['analytic_state_error_account'],
        error_components=r['error_components'],maximum_charge_rank=r['maximum_charge_rank'],
        maximum_bond_rank=r['maximum_bond_rank'],two_site_gates=r['two_site_gates'],
        elapsed_seconds=time.monotonic()-started,
        vector_comparison='all physical amplitudes, including outside the computational prefix')


def run():
    started=time.monotonic()
    result=dict(arithmetic=arithmetic(),operator_controls=field_and_commutators(),
        full_wavefunctions=wavefunctions(),physical_576_reference=physical_reference(),
        elapsed_seconds=time.monotonic()-started,status='PASS')
    return result


if __name__=='__main__':
    result=run();Path('controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
