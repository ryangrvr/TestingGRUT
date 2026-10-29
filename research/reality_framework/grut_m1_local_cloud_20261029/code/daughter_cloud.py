"""First-order direct h-daughter vs primary-beam scattering in a local plume.

Sample an explicit, outward-moving h; propagate its two decay daughters as
straight worldlines in flat space and integrate their local s-channel NWA
collisions with the opposite-charge primary radial beam. Spatially/time
resolved emissions use positive-projected windowed field fluxes. No produced
h, secondary scattering, core interaction, finite-width, or mixing feedback
is evolved. NWA caustics and source assumptions require convergence checks.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import simpson
from scipy.stats import qmc
from local_rates import pulse_series,GH,MH,GAMMA

ROOT=Path(__file__).resolve().parents[1]
BETA_STAR=np.sqrt(1-4/MH**2)
ESTAR=MH/2


def reconstruct_times(centers,width,uniform,stop):
    """Inverse CDF of the squared cosine window, truncated to [0,stop]."""
    def cdf(x):return x+.5+np.sin(2*np.pi*x)/(2*np.pi)
    lo=np.maximum(-.5,-centers/width);hi=np.minimum(.5,(stop-centers)/width)
    target=cdf(lo)+uniform*(cdf(hi)-cdf(lo))
    a=lo.copy();b=hi.copy()
    for _ in range(40):
        mid=(a+b)/2;less=cdf(mid)<target;a=np.where(less,mid,a);b=np.where(less,b,mid)
    return centers+width*(a+b)/2


def source_interp(t,E,tgrid,egrid,flux):
    """Bilinear source interpolation with explicit physical support."""
    dt=tgrid[1]-tgrid[0];de=egrid[1]-egrid[0]
    valid=(t>=tgrid[0])&(t<=tgrid[-1])&(E>=1.)&(E<=egrid[-1]+de/2)
    it=np.floor((t-tgrid[0])/dt).astype(int);ft=(t-tgrid[0])/dt-it
    it=np.clip(it,0,len(tgrid)-2);ft=np.clip(ft,0,1)
    # Hold the nearest nodal density across the half-bin endcaps.
    xe=(E-egrid[0])/de;ie=np.floor(xe).astype(int);fe=xe-ie
    fe=np.where(xe<=0,0,np.where(xe>=len(egrid)-1,1,fe));ie=np.clip(ie,0,len(egrid)-2)
    val=((1-ft)*((1-fe)*flux[it,ie]+fe*flux[it,ie+1])+
         ft*((1-fe)*flux[it+1,ie]+fe*flux[it+1,ie+1]))
    return np.where(valid,val,0.)


def root_rates(Ed,pd,mu,t,r,r0,tgrid,egrid,flux):
    """Integrate the on-shell delta over the primary beam energy."""
    A=Ed[:,None];B=pd[:,None]*mu;D=A*A-B*B;disc=4.-D
    good=disc>=0;root=np.sqrt(np.maximum(disc,0));total=np.zeros_like(mu)
    diagnostics={'near_tangent':0,'two_roots':0,'max_inverse_jacobian':0.}
    validity=[]
    for sign in (1.,-1.):
        Ep=(2*A+sign*B*root)/D;pp=np.sqrt(np.maximum(Ep*Ep-1,0))
        residual=A*Ep-B*pp-2
        valid=good&(Ep>1+1e-12)&(abs(residual)<1e-7)
        if sign<0:valid &= abs(B)>1e-10 # remove duplicate root at B=0
        vp=np.divide(pp,Ep,out=np.ones_like(Ep),where=Ep>0)
        deriv=2*(A-np.divide(B*Ep,pp,out=np.zeros_like(Ep),where=pp>0))
        inv=np.divide(1.,abs(deriv),out=np.zeros_like(deriv),where=abs(deriv)>1e-12)
        valid &= abs(deriv)>1e-12
        diagnostics['near_tangent']+=int(np.count_nonzero(good&(abs(deriv)<=1e-12)))
        safe_vp=np.where(valid & (vp>1e-12),vp,1.)
        safe_Ep=np.where(valid,Ep,1.)
        ts=t-(r-r0)/safe_vp
        src=source_interp(ts,safe_Ep,tgrid,egrid,flux)
        rate=src/(4*np.pi*r*r*safe_vp)*(np.pi*GH**2/(2*A*safe_Ep))*inv
        total+=np.where(valid,rate,0.)
        validity.append(valid)
        if np.any(valid):diagnostics['max_inverse_jacobian']=max(diagnostics['max_inverse_jacobian'],float(np.max(inv[valid])))
    diagnostics['two_roots']=int(np.count_nonzero(validity[0]&validity[1]))
    return total,diagnostics


def daughter_integrals(Ed,prad,ptrans,td,rd,r0,Rmax,tend,nsteps,tgrid,egrid,target):
    """Integrate along straight paths in z, using impact-parameter geometry."""
    p=np.sqrt(np.maximum(Ed*Ed-1,0));v=np.divide(p,Ed,out=np.zeros_like(p),where=Ed>0)
    mu=np.divide(prad,p,out=np.zeros_like(p),where=p>0)
    impact=rd*np.sqrt(np.maximum(0,1-mu*mu));z0=rd*mu
    active=(rd<Rmax)&(td<tend)&(v>1e-12)
    zout=np.sqrt(np.maximum(0,Rmax*Rmax-impact*impact))
    would_hit=active&(mu<0)&(impact<r0)
    core_z=-np.sqrt(np.maximum(0,r0*r0-impact*impact))
    hits=would_hit&(td+np.divide(core_z-z0,v,out=np.full_like(v,np.inf),where=v>0)<=tend)
    stop=np.where(would_hit,core_z,zout)
    stop=np.minimum(stop,z0+v*np.maximum(tend-td,0));stop=np.where(active,np.maximum(stop,z0),z0)
    result=np.zeros(len(Ed));diags={'near_tangent':0,'two_roots':0,'max_inverse_jacobian':0.}
    uu=np.linspace(0,1,nsteps)
    for negative in (True,False):
        if negative:
            a=np.minimum(z0,0);b=np.minimum(stop,0);select=active&(z0<0)&(b>a+1e-14)
        else:
            a=np.maximum(z0,0);b=np.maximum(stop,0);select=active&(stop>0)&(b>a+1e-14)
        ii=np.flatnonzero(select)
        if not len(ii):continue
        if negative:
            L0=np.log1p(abs(a[ii])/r0);L1=np.log1p(abs(b[ii])/r0);L=L0[:,None]+(L1-L0)[:,None]*uu
            zz=-r0*np.expm1(L);jac=r0*(L0-L1)[:,None]*np.exp(L)
        else:
            L0=np.log1p(a[ii]/r0);L1=np.log1p(b[ii]/r0);L=L0[:,None]+(L1-L0)[:,None]*uu
            zz=r0*np.expm1(L);jac=r0*(L1-L0)[:,None]*np.exp(L)
        rr=np.sqrt(impact[ii,None]**2+zz**2)
        tt=td[ii,None]+(zz-z0[ii,None])/v[ii,None]
        murel=np.divide(zz,rr,out=np.zeros_like(zz),where=rr>0)
        rates,d=root_rates(Ed[ii],p[ii],murel,tt,rr,r0,tgrid,egrid,target)
        result[ii]+=simpson(rates*jac/v[ii,None],x=uu,axis=1)
        for k in ('near_tangent','two_roots'):diags[k]+=d[k]
        diags['max_inverse_jacobian']=max(diags['max_inverse_jacobian'],d['max_inverse_jacobian'])
    return result,hits,active,diags


def run(z,power,seed,nsteps,dt,Rmax,tend,chunk_size):
    tgrid,egrid,plus,minus=pulse_series(z,dt);de=float(np.diff(z['energy_edges'])[0]);minus=minus/de
    h=z['neutral_positive'];weights=h.ravel();normal=weights.sum();cum=np.cumsum(weights)/normal
    sampler=qmc.Sobol(d=5,scramble=True,seed=seed);U=sampler.random_base2(power);n=len(U)
    width=float(z['window_width_mX_t']);edges=z['energy_edges']
    r0=float(z['probe_radius_mX_inverse']) if 'probe_radius_mX_inverse' in z else (2400+.5)*.025
    acc=[{'sum':0.,'sum2':0.,'max':0.,'hits':0,'active':0,'nonzero':0,
          'diag':{'near_tangent':0,'two_roots':0,'max_inverse_jacobian':0.}} for _ in range(2)]
    for start in range(0,n,chunk_size):
        V=U[start:start+chunk_size]
        chosen=np.minimum(np.searchsorted(cum,V[:,0]),len(cum)-1);win=chosen//h.shape[1];j=chosen%h.shape[1]
        centers=z['window_center_mX_t'][win]
        temit=reconstruct_times(centers,width,V[:,1],float(z['window_center_mX_t'][-1]))
        lo=np.maximum(MH,edges[j]);hi=edges[j+1]
        if np.any(hi<=lo):raise ValueError('neutral source bin wholly below physical threshold')
        Eh=lo+(hi-lo)*V[:,2];gamma=Eh/MH;beta=np.sqrt(np.maximum(0,1-1/gamma**2))
        delay=-np.log(np.maximum(1-V[:,3],1e-300))*gamma/GAMMA
        td=temit+delay;rd=r0+beta*delay;x=2*V[:,4]-1
        for which,(target,xx) in enumerate(((minus,x),(plus,-x))):
            energy=gamma*ESTAR*(1+beta*BETA_STAR*xx)
            radial=gamma*ESTAR*(beta+BETA_STAR*xx)
            transverse=ESTAR*BETA_STAR*np.sqrt(np.maximum(0,1-xx*xx))
            exposure,hits,active,diag=daughter_integrals(energy,radial,transverse,td,rd,r0,Rmax,tend,nsteps,tgrid,egrid,target)
            a=acc[which];a['sum']+=float(exposure.sum());a['sum2']+=float(exposure@exposure)
            a['max']=max(a['max'],float(exposure.max()));a['hits']+=int(hits.sum());a['active']+=int(active.sum());a['nonzero']+=int(np.count_nonzero(exposure>0))
            for key in ('near_tangent','two_roots'):a['diag'][key]+=diag[key]
            a['diag']['max_inverse_jacobian']=max(a['diag']['max_inverse_jacobian'],diag['max_inverse_jacobian'])
    rows=[]
    for particle,a in zip(('positive','antiparticle'),acc):
        mean=a['sum']/n;variance=max(0.,(a['sum2']-n*mean**2)/(n-1))
        rows.append({'daughter':particle,'mean_exposure_per_parent':mean,
                     'events_per_reference_core':normal*mean,'max_sample_exposure':a['max'],
                     'standard_error_iid_diagnostic':normal*np.sqrt(variance/n),
                     'direct_matching_surface_hits_per_core':normal*a['hits']/n,
                     'fraction_decayed_inside_domain_before_horizon':a['active']/n,
                     'nonzero_exposure_samples':a['nonzero'],'diagnostics':a['diag']})
    return {'neutral_source_count':float(normal),'samples':n,'seed':seed,'time_horizon_mX_inverse':tend,
            'outer_radius_mX_inverse':Rmax,'time_step':dt,'trajectory_steps_per_segment':nsteps,
            'rows':rows,
            'scope':'Unattenuated NWA scattering of direct h-decay daughters against primary radial beams; point-source radial targets, finite local horizon. No daughter-daughter collisions, rescattering, core boundary physics, finite-width effects, or interactions before r0.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='source_t400');p.add_argument('--power',type=int,default=13)
    p.add_argument('--seed',type=int,default=19);p.add_argument('--nsteps',type=int,default=65);p.add_argument('--dt',type=float,default=10.)
    p.add_argument('--rmax',type=float,default=46868.85);p.add_argument('--tend',type=float,default=100000.)
    p.add_argument('--chunk-size',type=int,default=8192);p.add_argument('--tag',required=True)
    a=p.parse_args();z=np.load(ROOT/(a.source+'.npz'))
    out={'configuration':vars(a),'result':run(z,a.power,a.seed,a.nsteps,a.dt,a.rmax,a.tend,a.chunk_size)}
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
