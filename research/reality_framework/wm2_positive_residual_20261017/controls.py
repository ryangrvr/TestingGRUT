"""Independent arithmetic, archived-Python and full physical-matrix checks."""
from fractions import Fraction as Q
from math import factorial,sqrt
from pathlib import Path
import importlib.util
import sys
import time
import numpy as np
from scipy.fft import dst
import enclosures as E
from layout import Layout,S5
import native
import rotation
from solver import Evolution

OLD=Path(__file__).resolve().parent.parent/'wm2_history_formation_20261016'
sys.path.append(str(OLD))
spec=importlib.util.spec_from_file_location('s7_reference_controls',OLD/'controls.py')
REF=importlib.util.module_from_spec(spec);spec.loader.exec_module(REF)
M=REF.M


def arithmetic():
    rng=np.random.default_rng(71017)
    cases=0
    for size in [1,3,41,300]:
        x=(rng.normal(size=size)+1j*rng.normal(size=size))*2.**rng.integers(-30,30,size=size)
        exact=sum((Q.from_float(float(v.real))**2+Q.from_float(float(v.imag))**2 for v in x),Q(0))
        lo,hi=E.norm2_bounds(x)
        assert Q.from_float(lo)<=exact<=Q.from_float(hi)
        cases+=1
    trig=0
    for x in [-3.14,-.5,-1e-10,0,.23,1,3.14]:
        for co in [False,True]:
            q=Q.from_float(x);degree=100+(0 if co else 1)
            exact=sum(((-1)**j*q**(2*j+(0 if co else 1))/factorial(2*j+(0 if co else 1)) for j in range(50)),Q(0))
            rem=abs(q)**(100 if co else 101)/factorial(100 if co else 101)
            a,b=E.trig_interval(E.point(x),co)
            assert Q.from_float(float(a))<=exact-rem<=exact+rem<=Q.from_float(float(b))
            trig+=1
    # A positive directly formed residual remains enclosed even if a squared
    # norm subtraction loses it entirely. No negative clipping is performed.
    w=np.array([[1.0],[0.0]],complex);c=np.array([1.0,2.**-40]);z=w.conj().T@c
    sigma,coord=E.contact_bounds(w,z,c,0,E.gram_upper(w))
    actual=2.**-40;subtraction=float(c@c-np.vdot(z,z).real)
    assert subtraction==0 and sigma>=actual>0
    pi=E.pi_interval()
    return dict(exact_rational_positive_norm_cases=cases,rational_Taylor_trig_cases=trig,
        pi_interval=pi,cancellation_example=dict(true_residual=actual,
        squared_norm_subtraction=subtraction,positive_residual_upper=sigma),
        long_double_significand_bits=native.LONG_DIGITS)


def algebra():
    rng=np.random.default_rng(20261017);rows=[]
    for counts,rank in [([1,1],3),([2,1],4),([2,2],4),([2,3],3)]:
        values=[-.25,.25];g=.9
        a=Layout(counts,values,rank,g);b=S5.Layout(counts,values,rank)
        x=rng.normal(size=a.size)+1j*rng.normal(size=a.size);x/=np.linalg.norm(x)
        z=rng.normal(size=rank)+1j*rng.normal(size=rank);z/=np.linalg.norm(z)
        e=float(np.linalg.norm(a.hint(x,z)-b.hint(x,z,g)))
        herm=abs(np.vdot(x,a.hint(x,z)).imag)
        de=float(np.linalg.norm(a.density(x)-b.density(x)))
        expected=[]
        for j in range(rank):
            expected.append(sum(float(np.sum(abs(a.block(x,k)[:,[i for i,o in enumerate(a.occ[k]) if j in o]])**2)) for k in range(1,a.M+1)))
        ve=float(np.max(abs(a.vacancy(x)-expected)))
        t=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0]
        matrices,delta,rowsum,eta=rotation.fast_matrices(t,a.M)
        interval_matrices,slow_delta,_,_=rotation.matrices(t,a.M)
        y=x.copy();a.beam(y,0,rank-1,matrices)
        expected=b.beam(x.copy(),0,rank-1,t)
        be=float(np.linalg.norm(y-expected))
        assert e<3e-14 and herm<3e-14 and de<3e-14 and ve<3e-14 and be<3e-14
        assert slow_delta<=delta
        rows.append(dict(counts=counts,rank=rank,dimension=a.size,emission_action_error=e,
            Hermitian_quadratic_form_imaginary=float(herm),density_error=de,vacancy_error=ve,
            Fock_rotation_error=be,fast_coefficient_error_upper=delta,
            independent_interval_coefficient_error_upper=slow_delta))
    return rows


def residual_control():
    rng=np.random.default_rng(719)
    from itertools import combinations_with_replacement
    rows=[]
    for counts in [[1,1],[2,1]]:
        modes=5;g=.9
        full=Layout(counts,[-.25,.25],modes,g)
        rank=2;small=Layout(counts,[-.25,.25],rank,g)
        frame=np.linalg.qr(rng.normal(size=(modes,rank))+1j*rng.normal(size=(modes,rank)))[0]
        contact=rng.normal(size=modes)+1j*rng.normal(size=modes);contact/=np.linalg.norm(contact)
        z=frame.conj().T@contact
        c=rng.normal(size=small.size)+1j*rng.normal(size=small.size);c/=np.linalg.norm(c)
        def lift(x):
            out=np.zeros(full.size,complex)
            for k in range(small.M+1):full.block(out,k)[:]=small.block(x,k)@REF.field_isometry(frame,k).T
            return out
        residual=full.hint(lift(c),contact)-lift(small.hint(c,z))
        actual=float(np.linalg.norm(residual))
        sigma=float(np.linalg.norm(contact-frame@z))
        bound=g*sqrt(sum(counts))/(sum(counts)+1)*sigma
        assert actual<=bound+2e-14
        rows.append(dict(counts=counts,full_dimension=full.size,full_action_residual=actual,
            positive_structural_operator_bound=bound))
    return rows


def full_wavefunctions():
    rows=[]
    for n,modes,delta in [(3,5,[0,-.25,.25]),(4,4,[0,-.25,.25,-.25])]:
        values=sorted(set(delta[1:]));counts=[delta[1:].count(v) for v in values]
        model=M.vacuum(n,modes,delta,.9,classes=True)
        for t in [.5,2,8]:
            e=Evolution(counts,values,modes,rank_cap=modes)
            result=e.run(t,strang_budget=.001,projection_budget=.01,max_seconds=90)
            assert result['status']=='COMPLETE',result
            exact=M.evolve(model['H'],M.unit(len(model['states']),model['seed']),t)
            error=float(np.linalg.norm(exact-REF.lift(e,model)))
            assert error<=result['analytic_state_error_account']+2e-11
            rows.append(dict(N=n,modes=modes,time=t,full_dimension=len(exact),
                full_wavefunction_error=error,complete_error_account=result['analytic_state_error_account'],
                error_components=result['error_components']))
    return rows


def physical_reference():
    started=time.monotonic();n,modes,t=3,576,36
    model=M.vacuum(n,modes,[0,-.25,.25],.9,classes=True)
    exact=M.evolve(model['H'],M.unit(len(model['states']),model['seed']),t)
    e=Evolution([1,1],[-.25,.25],modes)
    result=e.run(t,max_seconds=120)
    assert result['status']=='COMPLETE',result
    frame=dst(e.w,type=1,norm='ortho',axis=0)
    projection=np.zeros_like(e.c)
    for wi,word in enumerate(e.layout.words):
        k=len(word);start=e.layout.offset[wi]
        if not k:projection[start]=exact[model['index'][word,()]]
        elif k==1:
            field=np.array([exact[model['index'][word,(a,)]] for a in range(modes)])
            projection[start:start+e.layout.rank]=frame.conj().T@field
        else:
            tensor=np.zeros((modes,modes),complex)
            for a in range(modes):
                for b in range(a,modes):
                    value=exact[model['index'][word,(a,b)]]
                    tensor[a,b]=tensor[b,a]=value if a==b else value/sqrt(2)
            reduced=frame.conj().T@tensor@frame.conj()
            for j,(a,b) in enumerate(e.layout.occ[2]):projection[start+j]=reduced[a,b]*(1 if a==b else sqrt(2))
    # Used only as a comparison diagnostic, never as the production residual
    # certificate. Pad the missing-projection subtraction conservatively.
    missing=max(0.,float(np.vdot(exact,exact).real-np.vdot(projection,projection).real))
    error=sqrt(missing+float(np.linalg.norm(e.c-projection)**2))
    assert error<=result['analytic_state_error_account']+2e-9
    return dict(N=n,physical_receiver_modes=modes,time=t,physical_dimension=len(exact),
        full_wavefunction_error=error,complete_error_account=result['analytic_state_error_account'],
        error_components=result['error_components'],rank=e.layout.rank,
        maximum_coefficient_dimension=result['maximum_coefficient_dimension'],
        reference_norm_error=abs(float(np.vdot(exact,exact).real)-1),
        elapsed_seconds=time.monotonic()-started,uses_full_two_photon_tensor=True)


def all_checks():
    started=time.monotonic()
    return dict(arithmetic=arithmetic(),native_against_S5=algebra(),
        positive_complete_action_residual=residual_control(),
        physical_full_wavefunctions=full_wavefunctions(),physical_576_reference=physical_reference(),
        elapsed_seconds=time.monotonic()-started,status='PASS')
