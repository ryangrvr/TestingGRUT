#!/usr/bin/env python3
"""Finite 3D nonlinear Candidate 1 tests with charge-exact Hamiltonian splitting.

No filtering, damping, external controller, gravity or quantum dynamics.
See PROTOCOL.md and its pre-trajectory AMENDMENT_01 for the actual matrix.
"""
from pathlib import Path
import argparse, hashlib, json, platform, sys, time
import numpy as np
import scipy
from scipy.fft import rfftn, irfftn
from scipy.interpolate import CubicSpline

ROOT=Path(__file__).resolve().parent
BACKGROUND=ROOT.parent/'grut_candidate_1_20261018/code/results_profile.npz'
REFERENCE=ROOT.parent/'grut_candidate_1_20261018/code/results.json'
CONFIG={
    'a':1.,'b':.1,'Qbar_target':399.4384624223012,
    'T':40.,'diagnostic_interval':1.,'fft_workers':1,
    'integrator':'Strang: exact linear KG + local potential kicks',
    'boundary':'periodic Fourier collocation; no filter or damping',
    'preparation':'volume-preserving deformation M=diag(exp(epsilon),exp(-epsilon),1); on-grid fixed Q',
    'run_matrix':[
        {'name':'control','L':64.,'N':96,'dt':.1,'epsilon':0.},
        {'name':'small','L':64.,'N':96,'dt':.1,'epsilon':.05},
        {'name':'stress','L':64.,'N':96,'dt':.1,'epsilon':.2},
        {'name':'time','L':64.,'N':96,'dt':.05,'epsilon':.2},
        {'name':'coarse','L':64.,'N':64,'dt':.1,'epsilon':.2},
        {'name':'space','L':64.,'N':128,'dt':.1,'epsilon':.2},
        {'name':'domain','L':256/3,'N':128,'dt':.1,'epsilon':.2},
    ],
    'validity_relative_Q_drift':1e-9,'validity_relative_E_drift':2e-3,
    'operating_sampled_core_Q_fraction':.9,
    'scope':'finite classical scalar trajectories; not a global/quantum stability proof',
}


def js(o):
    if isinstance(o,np.generic): return o.item()
    if isinstance(o,np.ndarray): return o.tolist()
    raise TypeError(type(o).__name__)


def kick(q,p,step):
    S=1+q[2]; amp2=q[0]*q[0]+q[1]*q[1]
    common=(S*S-1)+CONFIG['b']*amp2
    p[0]-=step*common*q[0]
    p[1]-=step*common*q[1]
    p[2]-=step*((CONFIG['a']*(S*S-1)+amp2)*S-2*CONFIG['a']*q[2])


def run_case(case, reference, data):
    start=time.perf_counter();n=case['N'];L=case['L'];h=L/n;dt=case['dt']
    coords=(np.arange(n)-n//2)*h
    x,y,z=np.meshgrid(coords,coords,coords,indexing='ij',sparse=True)
    radius=np.sqrt(x*x+y*y+z*z)
    spline=CubicSpline(data['rho'],data['y'],axis=1)
    def profile(r):
        limited=np.minimum(r,data['rho'][-1]); result=spline(limited)
        result[:,r>data['rho'][-1]]=np.array([0.,0.,1.,0.])[:,None]
        return result[0],result[2]
    baseF,baseS=profile(radius)
    eps=case['epsilon'];deformed=np.sqrt(np.exp(2*eps)*x*x+np.exp(-2*eps)*y*y+z*z)
    F,S=profile(deformed)
    q=np.array([F,np.zeros_like(F),S-1]);p=np.zeros_like(q)
    Omega=CONFIG['Qbar_target']/(np.sum(F*F)*h**3)
    p[1]=-Omega*F
    k=2*np.pi*np.fft.fftfreq(n,d=h);kr=2*np.pi*np.fft.rfftfreq(n,d=h)
    k2=k[:,None,None]**2+k[None,:,None]**2+kr[None,None,:]**2
    w=np.sqrt(k2[None,:,:,:]+np.array([1.,1.,2.*CONFIG['a']])[:,None,None,None])
    cosine=np.cos(w*dt);sin_over_w=np.sin(w*dt)/w;w_sin=w*np.sin(w*dt)
    packed_weight=np.ones(n//2+1)*2;packed_weight[0]=packed_weight[-1]=1
    masks={str(r):radius<=r for r in [8,12]}
    background_norm=np.sum(baseF**2)*h**3
    organization_norm=np.sum((baseS-1)**2)*h**3
    rows=[]
    initial_expected=reference['Ebar']+(4/3)*reference['Egrad_bar']*np.sinh(eps)**2
    boundary=np.zeros((n,n,n),dtype=bool)
    for axis in range(3):
        sl=[slice(None)]*3;sl[axis]=0;boundary[tuple(sl)]=True
    initial_boundary_F=float(np.max(abs(F[boundary])))

    def diagnose(t):
        A=q[0]+1j*q[1];At=p[0]+1j*p[1];s=1+q[2]
        amp2=abs(A)**2;qd=-np.imag(np.conj(A)*At)
        kinetic=.5*np.sum(p*p,axis=0)
        U=.25*CONFIG['a']*(s*s-1)**2+.5*s*s*amp2+.25*CONFIG['b']*amp2**2
        qhat=rfftn(q,axes=(1,2,3),workers=CONFIG['fft_workers'])
        minus_lap=irfftn(k2[None,:,:,:]*qhat,s=(n,n,n),axes=(1,2,3),workers=CONFIG['fft_workers'])
        grad=.5*np.sum(q*minus_lap,axis=0)
        ed=kinetic+U+grad
        Q=float(np.sum(qd)*h**3); E=float(np.sum(ed)*h**3)
        # q*(-lap q) is correct globally but differs from local |grad q|² by a divergence.
        # Regional physical gradient energies are computed separately using spectral derivatives.
        local_grad=np.zeros_like(radius)
        for kval in [k[:,None,None],k[None,:,None],kr[None,None,:]]:
            derivative=irfftn(1j*kval[None,:,:,:]*qhat,s=(n,n,n),axes=(1,2,3),workers=CONFIG['fft_workers'])
            local_grad+=.5*np.sum(derivative*derivative,axis=0)
        local_ed=kinetic+U+local_grad
        phase=np.angle(np.sum(A*baseF))
        shape_error=float(np.sqrt(np.sum(abs(A*np.exp(-1j*phase)-baseF)**2)*h**3/background_norm))
        organization_error=float(np.sqrt(np.sum((s-baseS)**2)*h**3/organization_norm))
        positive=np.maximum(qd,0);positive_total=np.sum(positive)
        quadrupoles=np.array([np.sum(amp2*(x*x-y*y)),np.sum(amp2*(2*z*z-x*x-y*y)),
                              2*np.sum(amp2*x*y),2*np.sum(amp2*x*z),2*np.sum(amp2*y*z)])
        second_moment=float(np.sum(amp2*radius**2))
        rr=radius.ravel();ww=positive.ravel();bins=np.linspace(0,L*.9,601)
        counts,_=np.histogram(rr,bins=bins,weights=ww)
        csum=np.cumsum(counts)
        radius90=float(np.interp(.9*positive_total,csum,bins[1:]))
        row={'tau':t,'Qbar':Q,'Ebar':E,'core_Qbar':{},'core_Ebar':{},
             'core_Q_fraction':{},'core_E_fraction':{},'central_S':float(s[n//2,n//2,n//2]),
             'peak_A':float(np.max(abs(A))), 'positive_charge_radius90':radius90,
             'negative_charge_fraction':float(np.sum(np.maximum(-qd,0))/(positive_total+1e-300)),
             'phase_aligned_profile_error':shape_error,'organization_profile_error':organization_error,
             'quadrupole_components_over_second_moment':quadrupoles/second_moment,
             'quadrupole_norm':float(np.linalg.norm(quadrupoles)/second_moment),
             'spectral_local_gradient_total_difference':float((np.sum(local_ed)*h**3-E)/E)}
        for rad,mask in masks.items():
            cq=float(np.sum(qd[mask])*h**3);ce=float(np.sum(local_ed[mask])*h**3)
            row['core_Qbar'][rad]=cq;row['core_Ebar'][rad]=ce
            row['core_Q_fraction'][rad]=cq/CONFIG['Qbar_target'];row['core_E_fraction'][rad]=ce/E
        rows.append(row)
        if t==0 or abs(t%10)<1e-8:
            print(json.dumps({'case':case['name'],'tau':t,'E':E,'Q':Q,
                              'Q_in_r8':row['core_Q_fraction']['8'],
                              'shape_error':shape_error},default=js),flush=True)

    diagnose(0.)
    steps=round(CONFIG['T']/dt);every=round(CONFIG['diagnostic_interval']/dt)
    for step in range(1,steps+1):
        kick(q,p,.5*dt)
        qh=rfftn(q,axes=(1,2,3),workers=CONFIG['fft_workers'])
        ph=rfftn(p,axes=(1,2,3),workers=CONFIG['fft_workers'])
        newq=cosine*qh+sin_over_w*ph;newp=cosine*ph-w_sin*qh
        q=irfftn(newq,s=(n,n,n),axes=(1,2,3),workers=CONFIG['fft_workers'])
        p=irfftn(newp,s=(n,n,n),axes=(1,2,3),workers=CONFIG['fft_workers'])
        kick(q,p,.5*dt)
        if step%every==0: diagnose(float(step*dt))
    Q0=rows[0]['Qbar'];E0=rows[0]['Ebar']
    summary={'case':case,'h':h,'Omega_grid_initial':Omega,
             'initial_E_continuum_analytic':initial_expected,
             'initial_E_grid':E0,'initial_E_grid_minus_analytic':E0-initial_expected,
             'initial_boundary_abs_A_max':initial_boundary_F,
             'max_Q_relative_drift':max(abs(x['Qbar']/Q0-1) for x in rows),
             'max_E_relative_drift':max(abs(x['Ebar']/E0-1) for x in rows),
             'min_core_Q_fraction_r8':min(x['core_Q_fraction']['8'] for x in rows),
             'min_core_Q_fraction_r12':min(x['core_Q_fraction']['12'] for x in rows),
             'max_phase_aligned_profile_error':max(x['phase_aligned_profile_error'] for x in rows),
             'seconds':time.perf_counter()-start}
    summary['numerical_validity_pass']=summary['max_Q_relative_drift']<CONFIG['validity_relative_Q_drift'] and summary['max_E_relative_drift']<CONFIG['validity_relative_E_drift']
    summary['operating_diagnostic_pass']=summary['min_core_Q_fraction_r8']>CONFIG['operating_sampled_core_Q_fraction']
    report={'summary':summary,'trajectory':rows}
    (ROOT/f"evolution_{case['name']}.json").write_text(json.dumps(report,indent=2,default=js)+'\n')
    mid=n//2
    np.savez_compressed(ROOT/f"evolution_{case['name']}_slice.npz",coords=coords,
                        A_real=q[0,:,:,mid],A_imag=q[1,:,:,mid],S=1+q[2,:,:,mid])
    print('COMPLETE '+json.dumps(summary,default=js),flush=True)
    return summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',choices=[x['name'] for x in CONFIG['run_matrix']])
    args=parser.parse_args()
    prescription={'configuration':CONFIG,'versions':{'python':sys.version,'numpy':np.__version__,
                  'scipy':scipy.__version__,'platform':platform.platform()},
                  'background_sha256':hashlib.sha256(BACKGROUND.read_bytes()).hexdigest(),
                  'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'evolution_prescription.json').write_text(json.dumps(prescription,indent=2)+'\n')
    data=np.load(BACKGROUND);reference=json.loads(REFERENCE.read_text())['observables']
    cases=[x for x in CONFIG['run_matrix'] if args.case is None or x['name']==args.case]
    summaries=[run_case(x,reference,data) for x in cases]
    (ROOT/'evolution_run_summaries.json').write_text(json.dumps(summaries,indent=2,default=js)+'\n')


if __name__=='__main__':main()
