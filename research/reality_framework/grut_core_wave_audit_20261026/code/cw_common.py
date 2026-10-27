"""Exact local Hamiltonian ledger on the pinned periodic FD4 lattice."""
from pathlib import Path
import sys,json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[2]
PILOT=ROOT.parent/'grut_nonlinear_pilot_20261026'
sys.path.insert(0,str(PILOT/'code'))
import pilot as P
import analyze as A
from scipy.ndimage import gaussian_filter
C=P.C

def save(name,obj):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(obj,indent=2,default=convert,allow_nan=False)+'\n')

def read(name):return json.loads((ROOT/name).read_text())

def ledger(q,p,L):
    """Densities and exact instantaneous local derivatives; q,p are real."""
    d=L/q.shape[-1];r2=q[0]**2+q[1]**2
    j=q[0]*p[1]-q[1]*p[0]
    kX=.5*np.sum(p[:2]**2,axis=0);kphi=.5*p[2]**2
    mix=.5*C.G**2*q[2]**2*r2;selfX=C.LX/4*r2**2
    vphi=C.LP/4*(q[2]**2-C.NU**2)**2
    gradX=np.zeros_like(r2);gradphi=np.zeros_like(r2)
    dEX_flux=np.zeros_like(r2);dEphi_flux=np.zeros_like(r2);dQ=np.zeros_like(r2)
    force=-C.gradient(q)
    for axis in range(3):
        # Axis below refers only to spatial dimensions, not the field axis.
        for step,weight in [(1,4/3),(2,-1/12)]:
            coeff=weight/d**2
            qp=np.roll(q,-step,axis+1);qm=np.roll(q,step,axis+1)
            pp=np.roll(p,-step,axis+1);dq=qp-q
            force+=coeff*(qp+qm-2*q)
            gx=.25*coeff*np.sum(dq[:2]**2,axis=0)
            gs=.25*coeff*dq[2]**2
            gradX+=gx+np.roll(gx,step,axis)
            gradphi+=gs+np.roll(gs,step,axis)
            fx=.5*coeff*np.sum(dq[:2]*(p[:2]+pp[:2]),axis=0)
            fs=.5*coeff*dq[2]*(p[2]+pp[2])
            fq=coeff*(q[0]*qp[1]-q[1]*qp[0])
            dEX_flux+=fx-np.roll(fx,step,axis)
            dEphi_flux+=fs-np.roll(fs,step,axis)
            dQ+=fq-np.roll(fq,step,axis)
    work=C.G**2*q[2]*p[2]*r2
    return dict(j=j,r2=r2,kX=kX,kphi=kphi,gradX=gradX,gradphi=gradphi,
                mix=mix,selfX=selfX,vphi=vphi,
                eX=kX+gradX+mix+selfX,ephi=kphi+gradphi+vphi,
                dEX_flux=dEX_flux,dEphi_flux=dEphi_flux,dQ=dQ,
                work_to_X=work,force=force)

def region(q,p,fields,mask,L):
    """A frozen region at one instant, never an isolated soliton by definition."""
    w=np.asarray(mask,float);dv=(L/q.shape[-1])**3
    def integral(x):return float(np.sum(w*x)*dv)
    def mean(x):return float(np.sum(w*x)/w.sum())
    Q=integral(fields['j']);I=integral(fields['r2']);omega=Q/I if I>0 else 0.
    r2=fields['r2'];safe=np.maximum(r2,1e-300)
    radial=.5*(q[0]*p[0]+q[1]*p[1])**2/safe
    shear=.5*(fields['j']-omega*r2)**2/safe
    rot=Q**2/(2*I) if I>0 else 0.
    K=integral(fields['kX']);Kr=integral(radial);Ks=integral(shear)
    Kn=integral(np.where(r2==0,fields['kX'],0.))
    positive=np.maximum(fields['j'],0)*w
    psi=(q[0]+1j*q[1])/np.sqrt(safe)
    coherence=float(abs(np.sum(positive*psi))/positive.sum()) if positive.sum()>0 else None
    kinetic_res=.5*((p[0]+omega*q[1])**2+(p[1]-omega*q[0])**2)
    stationary_res=np.sum((fields['force'][:2]+omega**2*q[:2])**2,axis=0)
    EX=integral(fields['eX']);Ephi=integral(fields['ephi']);E=EX+Ephi
    fluxX=-integral(fields['dEX_flux']);fluxphi=-integral(fields['dEphi_flux'])
    fluxQ=-integral(fields['dQ']);work=integral(fields['work_to_X'])
    return dict(cells=int(np.count_nonzero(w)),volume=integral(np.ones_like(w)),Q=Q,
        E_X=EX,E_phi=Ephi,E_total=E,E_over_Q=E/Q if Q!=0 else None,
        E_X_over_Q=EX/Q if Q!=0 else None,omega_fit=omega,
        kinetic_X=K,kinetic_common_rotation=rot,kinetic_radial=Kr,kinetic_angular_dispersion=Ks,kinetic_at_exact_nodes=Kn,
        kinetic_excess_fraction=(Kr+Ks+Kn)/K if K>0 else None,
        kinetic_identity_relative=abs(K-rot-Kr-Ks-Kn)/max(K,1),
        rotation_fit_identity_relative=abs(integral(kinetic_res)-Kr-Ks-Kn)/max(K,1),
        phase_coherence_weighted_positive_current=coherence,
        stationary_X_force_residual_over_m2=float(np.sqrt(integral(stationary_res)/max(I,1e-300))),
        phi_velocity_RMS_over_mv=float(np.sqrt(mean(p[2]**2))/C.NU),
        phi_acceleration_RMS_over_m2v=float(np.sqrt(mean(fields['force'][2]**2))/C.NU),
        outward_energy_flux_X=fluxX,outward_energy_flux_phi=fluxphi,outward_charge_flux=fluxQ,
        reversible_work_phi_to_X=work,
        fixed_region_dE_X=-fluxX+work,fixed_region_dE_phi=-fluxphi-work,
        fixed_region_dE_total=-fluxX-fluxphi,fixed_region_dQ=-fluxQ,
        fixed_region_fractional_dE_total=(-fluxX-fluxphi)/E if E!=0 else None,
        fixed_region_fractional_dQ=-fluxQ/Q if Q!=0 else None,
        components={k:integral(fields[k]) for k in ['kX','kphi','gradX','gradphi','mix','selfX','vphi']},
        scope='Fixed-region instantaneous exact ledger with shared edge energies. No ambient subtraction, mass certificate, escape spectrum or assigned particle count.')

def distances(N,L,centre):
    axis=np.arange(N)*L/N
    diff=[(axis-c+L/2)%L-L/2 for c in centre]
    return np.sqrt(diff[0][:,None,None]**2+diff[1][None,:,None]**2+diff[2][None,None,:]**2)

def sign_shell(smoothed_phi,rr,r90):
    shell=(rr>=r90+1)&(rr<r90+3)
    if not shell.any():return dict(cells=0)
    x=smoothed_phi[shell]/C.NU
    return dict(cells=int(shell.sum()),r_inner=r90+1,r_outer=r90+3,
                strong_positive_volume_fraction=float(np.mean(x>.5)),
                strong_negative_volume_fraction=float(np.mean(x<-.5)),
                weak_sign_volume_fraction=float(np.mean(abs(x)<=.5)),
                scope='Instantaneous smoothed two-sign neighborhood, not a persistent/topological wall certificate.')
