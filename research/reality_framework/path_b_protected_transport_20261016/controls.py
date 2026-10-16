"""Reproducible all-terms mediator and geometry-change controls."""
from __future__ import annotations
import json
import math
from pathlib import Path
import platform
import sys
import numpy as np
import scipy
from scipy import sparse
from scipy.linalg import eigvalsh, eigh
from scipy.sparse.linalg import splu, expm_multiply
import module as wm
from cp12_activation import run as activation

P=Path(__file__).resolve().parent


def module_control(kind,size):
    a=wm.cycle(size) if kind=='cycle' else wm.torus(size)
    n=a.shape[0]; z=int(a.sum(axis=1).max()); d=wm.distances(a)
    h=wm.microscopic(a); f,w,spectrum=wm.exact_low_band(a)
    eye=np.eye(n)
    schur=-.25**2*splu((8*sparse.eye(n)-a).tocsc()).solve(eye)
    const=wm.constants(z=4)
    q=4/8
    bound=.25**2/(8-4)*q**d
    # Ordinary floating-point checks have an absolute roundoff floor.
    schur_violation=float(np.max(np.abs(schur[:,0])-bound))
    dressed_bound=const['dressed_C']*np.exp(-const['mu']*d)
    dressed_violation=float(np.max(np.abs(f[:,0])-dressed_bound))
    nearest=np.flatnonzero(d==1)
    near=float(np.min(np.abs(f[nearest,0])))
    residual=float(np.linalg.norm(h@w-w@f,ord='fro'))
    ortho=float(np.linalg.norm(w.T@w-eye,ord='fro'))
    far=int(np.argmax(d)); far_schur=abs(float(schur[far,0]))
    green_error=None
    if kind=='cycle':
        green=np.array([-.25**2*wm.ring_green(n,int(di)) for di in d])
        green_error=float(np.max(np.abs(green-schur[:,0])))
    assert schur_violation<2e-13 and dressed_violation<2e-13
    assert residual<2e-10 and ortho<2e-10
    assert near>const['dressed_near_lower']>0
    assert spectrum['high_min']-spectrum['low_max']>=const['spectral_gap_lower']-1e-12
    # Check localized dressing in the same independently defined metric.
    weights=np.exp(const['mu']*d)
    dl=float(np.sum(np.abs(w[:n,0]-eye[:,0])*weights))
    dh=float(np.sum(np.abs(w[n:,0])*weights))
    # Weighting a numerical roundoff floor at large d is not a proof check.
    # Limit direct weighted sums to diameters <= 16; the analytic bound is all N.
    dressing_check=None
    if d.max()<=16:
        assert dl<const['dressing_L_minus_identity_weighted_upper']+1e-8
        assert dh<const['dressing_H_weighted_upper']+1e-8
        dressing_check=dict(L_minus_I=dl,H=dh)
    return dict(kind=kind,size=size,N=n,max_degree=z,diameter=int(d.max()),
        spectrum=spectrum,near_dressed_coupling=near,far_receiver=far,
        far_E0_Feshbach_coupling=far_schur,
        far_E0_Feshbach_bound=float(bound[far]),
        far_dressed_entry_float=float(f[far,0]),
        roundoff_notice='Dressed entries at the numerical noise floor are not resolved signals.',
        intertwining_Frobenius_error=residual,orthogonality_Frobenius_error=ortho,
        max_Feshbach_bound_excess=schur_violation,
        max_dressed_bound_excess=dressed_violation,
        cycle_Green_function_max_error=green_error,
        localized_dressing_check=dressing_check)


def latch_control(n, gap=8.,delta=8.,lam=.25):
    h=wm.latch(wm.cycle(n),gap=gap,delta=delta,lam=lam)
    q=h[n:,n:].tocsc(); b=h[n:,:n].tocsc(); j=n//2
    c=wm.latch_constants(n,gap=gap,delta=delta,lam=lam)
    rhs=b[:,0].toarray().ravel()
    rows=[]
    for energy in [0.,c['low_bottom_lower']]:
        sol=splu(q-energy*sparse.eye(3*n,format='csc')).solve(rhs)
        result=-np.asarray(b.T@sol).ravel()
        far=abs(float(result[j])); near=abs(float(result[1]))
        assert far+1e-16>=c['shortcut_floor_all_low_energies']
        assert near>=c['nearest_lower_all_low_energies']-1e-13
        residual=float(np.linalg.norm((q-energy*sparse.eye(3*n))@sol-rhs))
        assert residual<1e-12
        rows.append(dict(E=energy,near=near,far=far,far_near=far/near,
                         linear_solve_residual=residual))
    spectrum=None
    if n<=64:
        eigen,u=eigh(h.toarray()); qeigen=eigvalsh(q.toarray())
        assert eigen[n-1]<=1e-12 and eigen[n]>=c['full_spectral_gap_lower']-1e-12
        assert eigen[0]>=c['low_bottom_lower']-1e-12
        assert qeigen[0]>=c['Q_gap_lower']-1e-12
        spectrum=dict(low_min=float(eigen[0]),low_max=float(eigen[n-1]),
            high_min=float(eigen[n]),Q_min=float(qeigen[0]),
            actual_gap=float(eigen[n]-eigen[n-1]))
        low=u[:,:n]; r=low[:n,:]; gram=r@r.T
        gv,gu=eigh(gram); inv=(gu/np.sqrt(gv))@gu.T
        w=low@r.T@inv
        dressed=w.T@(h@w)
        spectrum['canonical_dressed_far']=abs(float(dressed[j,0]))
        spectrum['canonical_dressed_near']=abs(float(dressed[1,0]))
        spectrum['canonical_isometry_error']=float(np.linalg.norm(w.T@w-np.eye(n)))
        assert spectrum['canonical_isometry_error']<1e-9
        initial=np.zeros(4*n); initial[0]=1.
        probabilities=[]
        for t in [.25,1.,5.,32.]:
            state=expm_multiply(-1j*t*h,initial)
            total_q=float(np.vdot(state[n:],state[n:]).real)
            geom_one=float(np.vdot(state[2*n:],state[2*n:]).real)
            assert total_q<c['bare_P0_total_Q_probability_upper']+1e-10
            probabilities.append(dict(t=t,Q_probability=total_q,
                geometry_one_probability=geom_one,
                receiver_L0_probability=float(abs(state[j])**2)))
        spectrum['full_dynamics_checks']=probabilities
    return dict(N=n,distance=j,Delta=gap,delta=delta,lam=lam,
        constants=c,Schur=rows,spectrum=spectrum)


def dense_control(n):
    gap=8.; rate=1.; g=.25
    # Full normalized complete graph, not a deleted-edge comparison law.
    off=g*g*rate/((n-1)*(gap-rate)*(gap+rate/(n-1)))
    lower=(n-1)*off
    return dict(N=n,reference_ring_distance=n//2,normalized_full_degree_rate=rate,
        Q_gap=gap-rate,far_Feshbach_E0=off,
        N_minus_1_times_far=lower,
        meaning='A size-stable gap and bounded total hopping still give a 1/N tail in the independently fixed ring metric.')


def source_control(n):
    a=wm.cycle(n); h=wm.microscopic(a); f,w,_=wm.exact_low_band(a)
    packet=np.zeros(n); packet[0]=packet[1]=1/math.sqrt(2)
    psi=w@packet
    number_psi=np.zeros(2*n); number_psi[0]=psi[0]
    p=float(psi[0]**2); var=p*(1-p)
    assert var>.24
    rows=[]; d=wm.distances(a)
    for t in [.25,1.,32.,128.]:
        x=expm_multiply(-1j*t*h,psi)
        y=expm_multiply(-1j*t*h,number_psi)
        chi=2*np.imag(x[:n].conj()*y[:n])
        # Independent finite pulse identity for a projection-valued source.
        plus=psi+(np.exp(-1j*np.pi/2)-1)*number_psi
        minus=psi+(np.exp(1j*np.pi/2)-1)*number_psi
        xp=expm_multiply(-1j*t*h,plus)
        xm=expm_multiply(-1j*t*h,minus)
        fd=(np.abs(xp[:n])**2-np.abs(xm[:n])**2)/2
        err=float(np.max(np.abs(fd-chi)))
        assert err<2e-10
        bare=np.zeros(2*n); bare[0]=1.
        bare=expm_multiply(-1j*t*h,bare)
        high=float(np.vdot(bare[n:],bare[n:]).real)
        assert high<=wm.constants()['bare_high_probability_upper']+1e-10
        propagation_bound=np.exp(wm.constants()['microscopic_propagation_exponent_rate']*t
            -wm.constants()['mu']*d)
        assert np.all(np.abs(bare[:n])<=propagation_bound+2e-11)
        rows.append(dict(t=t,norm_error=abs(float(np.linalg.norm(x))-1),
            finite_pulse_error=err,
            bare_source_high_probability=high,
            receiver_response=[dict(j=j,distance=int(d[j]),chi=float(chi[j]))
                for j in range(n//2+1)]))
    return dict(N=n,variance=var,preparation='Supplied dressed (|0>+|1>)/sqrt(2) fixture. NOT the CP-12 neutral preparation.',
        physical_source='Bare local neutral-probe L-number phase; full two-band H retained.',rows=rows)


def main():
    data=dict(status='PASS',environment=dict(python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__),
        constants=wm.constants(),
        modules=[module_control('cycle',n) for n in [8,16,32,64,128,256]]+
            [module_control('torus',l) for l in [4,6,8,12,16]],
        latch=[latch_control(n) for n in [8,16,32,64,128,256,512]],
        robustness=[latch_control(64,gap=d,delta=d,lam=l)
                    for d in [8.,12.,16.] for l in [0.,.125,.25,.5]],
        dense=[dense_control(n) for n in [8,16,32,64,128,256,512,1024,4096]],
        source=[source_control(n) for n in [16,32,64]],
        CP12_activation=activation(),
        interpretation='Conditional all-orders transport module: proved. Complete WM1 architecture: not accepted. Coupled original-neutral spatial phase and universal protected matter: not established.',
        arithmetic='Analytic all-size bounds and exact rational CP-12 words; ordinary floating-point matrix controls, not outward-rounded numerical certificates.')
    (P/'RESULTS.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=data['status'],module_rows=len(data['modules']),
        latch_rows=len(data['latch']),robustness_rows=len(data['robustness']),
        activation_rows=len(data['CP12_activation']['rows']),constants=data['constants'])),flush=True)


if __name__=='__main__': main()
