#!/usr/bin/env python3
"""Controlled quantum-mode tests. Neither is the H55 cosmological Dirac sea.

1. Inherited 1D line solver, homogeneous pulse, against independent DOP853.
2. Full 3D four-component Dirac *one-mode* propagation, resolved transverse
   background, spatial/time refinement, unitary norm and vacuum control.
"""
import math
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp
from moments import TRACKS, RES, dump

sys.path.insert(0,str(TRACKS/'lf_fermion/code'))
import lf1d


def line_run(dt,dk=.08,A=10.,M=1.6,static=False):
    L=2*np.pi/dk;N=1<<int(math.ceil(math.log2(2*(1.3*A+3*M+4)/dk)))
    line=lf1d.Line(N,L,[M],n_r=1,theta0=2*A,phases='ones')
    nstep=math.ceil(math.pi/dt);h=math.pi/nstep
    for i in range(nstep):
        theta=2*A if static else 2*A*math.cos((i+.5)*h)
        line.step(i*h,h,1.,theta)
    theta_end=2*A if static else -2*A
    occ=line.occupation(1.,theta_end)[0]
    norm=np.sum(abs(line.psik[0,0])**2,axis=0)
    indices=np.array([int(round(k/dk))%N for k in [-18,-10,-7,-3,0,3,7,10,18]])
    k=line.k[indices]
    # Independent unrotated two-level ODE, no split propagation or inherited spinor helper.
    vec=[]
    for p in k:
        H=np.array([[p,M*np.exp(2j*A)],[M*np.exp(-2j*A),-p]])
        _,u=np.linalg.eigh(H);vec.append(u[:,0])
    vec=np.array(vec).T
    def rhs(t,z):
        z=z.reshape(2,-1);th=2*A if static else 2*A*math.cos(t)
        c=M*np.exp(1j*th)
        return (-1j*np.array([k*z[0]+c*z[1],c.conjugate()*z[0]-k*z[1]])).ravel()
    ode=solve_ivp(rhs,[0,math.pi],vec.ravel(),rtol=2e-11,atol=2e-13,method='DOP853')
    final=ode.y[:,-1].reshape(2,-1);nref=[]
    for i,p in enumerate(k):
        H=np.array([[p,M*np.exp(1j*theta_end)],[M*np.exp(-1j*theta_end),-p]])
        _,u=np.linalg.eigh(H);nref.append(abs(np.vdot(u[:,1],final[:,i]))**2)
    return dict(N=N,dt=h,dk=dk,A=A,M=M,static=static,
        energy=float(line.eps3d(occ[None,:],1.)[0]),
        max_norm_error=float(abs(norm-1).max()),max_occupation=float(occ.max()),
        min_occupation=float(occ.min()),sample_k=k.tolist(),sample_occupation=occ[indices].tolist(),
        independent_occupation=nref,max_sample_absolute_error=float(abs(occ[indices]-nref).max()))


def sigma(kx,ky,kz,u):
    return np.stack([kz*u[0]+(kx-1j*ky)*u[1],(kx+1j*ky)*u[0]-kz*u[1]])


def dirac3d(N,dt,profile='full',static=False):
    """4-component H=[[sigma.k,M exp(i theta)],[M exp(-i theta),-sigma.k]]."""
    M=1.;T=2.;L=2*math.pi
    axis=2*np.pi*np.fft.fftfreq(N,d=L/N)
    kx,ky,kz=np.meshgrid(axis,axis,axis,indexing='ij');km=np.sqrt(kx*kx+ky*ky+kz*kz)
    xx=np.arange(N)*L/N;x,y,z=np.meshgrid(xx,xx,xx,indexing='ij')
    shape=(1.+.3*np.cos(x)+.4*np.cos(y)+.2*np.cos(z)) if profile=='full' else np.ones_like(x)
    k0=np.array([1,0,1]);sk=np.array([[1,1],[1,-1]],complex)
    H0=np.block([[sk,M*np.eye(2)],[M*np.eye(2),-sk]])
    _,V=np.linalg.eigh(H0)
    psi=np.zeros((4,N,N,N),complex);psi[:,1,0,1]=V[:,0]
    steps=math.ceil(T/dt);h=T/steps
    cf=np.cos(km*h/2);sf=np.sinc(km*h/(2*np.pi))*h/2
    cm,sm=math.cos(M*h),math.sin(M*h)
    def free(p):
        return np.concatenate([cf*p[:2]-1j*sf*sigma(kx,ky,kz,p[:2]),
                               cf*p[2:]+1j*sf*sigma(kx,ky,kz,p[2:])])
    for j in range(steps):
        psi=free(psi);xp=np.fft.ifftn(psi,axes=(1,2,3),norm='ortho')
        theta=0. if static else 1.2*math.sin(math.pi*(j+.5)*h/T)**2*shape
        phase=np.exp(1j*theta)
        rp=cm*xp[:2]-1j*sm*phase*xp[2:]
        lp=cm*xp[2:]-1j*sm*np.conj(phase)*xp[:2]
        psi=free(np.fft.fftn(np.concatenate([rp,lp]),axes=(1,2,3),norm='ortho'))
    norm=float(np.sum(abs(psi)**2))
    hp=np.concatenate([sigma(kx,ky,kz,psi[:2])+M*psi[2:],M*psi[:2]-sigma(kx,ky,kz,psi[2:])])
    positive=.5*(psi+hp/np.sqrt(km*km+M*M))
    nk=np.sum(abs(positive)**2,axis=0)
    occupied=np.sum(abs(psi)**2,axis=0)
    result=dict(N=N,L=L,dt=h,steps=steps,profile=profile,static=static,
        norm_error=abs(norm-1),positive_probability=float(nk.sum()),
        positive_energy=float(np.sum(nk*np.sqrt(km*km+M*M))),
        transverse_sideband_probability=float(np.sum(occupied[(kx!=k0[0])|(ky!=k0[1])])),
        uv_edge_probability=float(np.sum(occupied[km>.7*np.max(abs(axis))])),
        full_cosmological_fermion_solution=False,filled_sea=False)
    if profile=='hom':
        def rhs(t,p):
            th=0. if static else 1.2*math.sin(math.pi*t/T)**2
            H=np.block([[sk,M*np.exp(1j*th)*np.eye(2)],[M*np.exp(-1j*th)*np.eye(2),-sk]])
            return -1j*H@p
        sol=solve_ivp(rhs,[0,T],V[:,0],rtol=2e-12,atol=2e-14,method='DOP853')
        p=sol.y[:,-1];proj=.5*(np.eye(4)+H0/math.sqrt(3))
        result['independent_positive_probability']=float(np.linalg.norm(proj@p)**2)
        result['independent_absolute_error']=abs(result['positive_probability']-result['independent_positive_probability'])
    return result


def main():
    start=time.time();line=[];three=[]
    for dt,dk,A,stat in [(0.004,.08,10.,False),(.002,.08,10.,False),(.001,.08,10.,False),
                         (.002,.04,10.,False),(.002,.08,42.3,False),(.001,.08,42.3,False),
                         (.004,.08,10.,True),(.002,.08,10.,True),(.001,.08,10.,True)]:
        r=line_run(dt,dk,A,static=stat);line.append(r)
        print('line',A,dt,dk,stat,'energy',r['energy'],'ODE error',r['max_sample_absolute_error'],flush=True)
    for N,dt,prof,stat in [(16,.02,'full',False),(16,.01,'full',False),(16,.005,'full',False),
                           (24,.01,'full',False),(32,.01,'full',False),
                           (16,.01,'hom',False),(16,.005,'hom',False),
                           (16,.01,'hom',True),(16,.005,'hom',True)]:
        r=dirac3d(N,dt,prof,stat);three.append(r)
        print('3D',N,dt,prof,stat,'P+',r['positive_probability'],'norm',r['norm_error'],flush=True)
    need=math.ceil(1.2e5*.42/math.pi)
    out=dict(line=line,three_dimensional_mode=three,runtime_seconds=time.time()-start,
        cosmological_resolution=dict(reported_H84_hotspot_comoving_k=1.2e5,L=.42,
          minimum_cells_per_axis=need,sites=need**3,complex128_bytes_one_four_spinor=need**3*4*16,
          notes='Nyquist lower bound only; no convergence margin, replicas or work arrays included'),
        verdict='algorithm controls only; no H55 3D sea, backreaction or drain recalibration')
    assert max(r['max_norm_error'] for r in line)<1e-9
    assert max(r['norm_error'] for r in three)<1e-9
    dump(RES/'quantum_controls.json',out)


if __name__=='__main__':main()
