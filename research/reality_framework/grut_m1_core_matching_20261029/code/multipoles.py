"""Time-resolved surface multipoles in the recreated periodic 3D test.

We project the fields and radial derivatives, including a matched undisturbed
control. An angular-gradient fraction is not a particle collision kernel.
Surfaces are inside the 64-wide periodic box; near-field, interpolation,
finite-duration and periodic-return effects remain explicit.
"""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.special import sph_harm_y,roots_legendre

ROOT=Path(__file__).resolve().parents[1]
DYN=ROOT.parent/'grut_m1_dynamics_20261027'
sys.path.insert(0,str(DYN/'code'))
from common import snapshot

def sphere(ntheta,nphi,lmax):
    mu,w=roots_legendre(ntheta);phi=2*np.pi*np.arange(nphi)/nphi
    theta=np.arccos(mu);st=np.sqrt(1-mu*mu)
    directions=np.stack(np.broadcast_arrays(st[:,None]*np.cos(phi),st[:,None]*np.sin(phi),mu[:,None]),axis=0).reshape(3,-1)
    weights=np.broadcast_to(w[:,None]*(2*np.pi/nphi),(ntheta,nphi)).ravel()
    modes=[(l,m) for l in range(lmax+1) for m in range(-l,l+1)]
    Y=np.stack([sph_harm_y(l,m,theta[:,None],phi[None,:]).ravel() for l,m in modes])
    P=np.conj(Y)*weights[None,:]
    gram=P@Y.T;orth_err=float(np.max(np.abs(gram-np.eye(len(modes)))))
    return directions,modes,P,orth_err

def sample(a,N,L,r,directions,mode='full'):
    dx=L/N;coords=N/2+(r/dx)*directions
    def get(i):return map_coordinates(a[i],coords,order=1,mode='nearest',prefilter=False)
    return np.stack((get(0)+1j*get(1),get(2)-10+0j,get(3)+1j*get(4),get(5)+0j))

def integrate(times,vals):return np.trapezoid(vals,x=times,axis=0)

def analyze(N,radii,max_time,ntheta,nphi,lmax,dr,suffix=''):
    directions,modes,P,orth=sphere(ntheta,nphi,lmax)
    base=ROOT/'.work'/f'object{N}_control{suffix}';dist=ROOT/'.work'/f'object{N}_perturbed{suffix}'
    names=sorted(set(x.name for x in base.glob('state_s*.bin')) & set(x.name for x in dist.glob('state_s*.bin')))
    data=[];times=[]
    for name in names:
        n,L,t,a=snapshot(base/name);n2,L2,t2,b=snapshot(dist/name)
        if t>max_time+1e-8:continue
        if (n,L,t)!=(n2,L2,t2) or n!=N:raise ValueError('Mismatched 3D controls')
        if max(radii)+dr>=L/2:raise ValueError('Surface exceeds half the periodic box')
        frame=[]
        for r in radii:
            work=[]
            for arr in (a,b):
                center=sample(arr,n,L,r,directions)
                inner=sample(arr,n,L,r-dr,directions)[:2];outer=sample(arr,n,L,r+dr,directions)[:2]
                q=center[:2]@P.T;mom=center[2:]@P.T;dq=((outer-inner)/(2*dr))@P.T
                work.append(np.stack((q,mom,dq)))
            # full perturbation and its difference against the identical-grid control
            frame.append(np.stack((work[1],work[1]-work[0])))
        data.append(frame);times.append(t)
    if len(times)<3:raise RuntimeError('Not enough saved snapshots')
    times=np.array(times);data=np.array(data) # t,r,case,q-p-dq,channel,lm
    rows=[];ell=np.array([l for l,m in modes]);ll=ell*(ell+1)
    for ir,r in enumerate(radii):
        for ic,case in enumerate(('full_disturbed','difference_to_control')):
            for ch,label in enumerate(('charged_complex','neutral_real')):
                Q=data[:,ir,ic,0,ch,:];Pmom=data[:,ir,ic,1,ch,:];DQ=data[:,ir,ic,2,ch,:]
                amp=integrate(times,abs(Q)**2);rad=integrate(times,abs(DQ)**2)
                by_l=[float(amp[ell==l].sum()) for l in range(lmax+1)]
                angular=float(np.dot(amp,ll)/r**2);radial=float(rad.sum());total=sum(by_l)
                flux=integrate(times,-np.real(np.conj(Pmom)*DQ)).sum()
                rows.append({'radius_mX_inverse':r,'case':case,'channel':label,
                             'field_power_by_l':by_l,'field_power_fraction_by_l':[x/total for x in by_l] if total>0 else None,
                             'integrated_angular_gradient_over_total_gradient':angular/(angular+radial) if angular+radial>0 else None,
                             'radial_surface_flux_integral_unscaled':float(flux),
                             'warning':'Field-energy projection; not a phase-space or collision-angle distribution.'})
    return {'N':N,'L':L,'times_mX_inverse':times.tolist(),'angular_modes':modes,'orthonormality_max_error':orth,
            'rows':rows,'last_time':float(times[-1]),'dr':dr},data

def main():
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True);p.add_argument('--max-time',type=float,default=30.)
    p.add_argument('--ntheta',type=int,default=24);p.add_argument('--nphi',type=int,default=48);p.add_argument('--lmax',type=int,default=4)
    p.add_argument('--dr',type=float,default=.75);p.add_argument('--suffix',default='');p.add_argument('--tag',required=True);a=p.parse_args()
    result,data=analyze(a.N,[16.,20.,24.],a.max_time,a.ntheta,a.nphi,a.lmax,a.dr,a.suffix)
    np.savez_compressed(ROOT/(a.tag+'.npz'),coefficients=data,times=result['times_mX_inverse'],modes=np.array(result['angular_modes']))
    out={'configuration':vars(a),'result':result,
         'scope':'Recreated periodic 64-wide box, 5% radial+quadrupole disturbance, first 30/m_X, N=64 or 96; interpolation and finite-radius near field. No 60/m_X surface, 10% disturbance, unique quantum occupation, true outgoing phase-space distribution, or pair-angle correlation.'}
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
