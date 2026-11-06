"""One-axis scalar and sterile field with a continuously coupled signed vacuum bath."""
from dataclasses import dataclass, replace
import math
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline

TRACKS=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(TRACKS/'x8_unitary/code'))
from unitary import Configuration, UnitaryMass, adjoint
from slab import SIGMA_X, SIGMA_Z
sys.path.insert(0,str(TRACKS/'x9_memory/code'))
from bath import build, mass_matched


@dataclass(frozen=True)
class State:
    a: float
    phi: np.ndarray
    pi: np.ndarray
    frame: np.ndarray


@dataclass(frozen=True)
class Geometry:
    cfg: Configuration
    bath_bins: int = 8
    bath_cutoff: float = 5.5
    gamma: float = 0.012
    reference_nodes: int = 29


class CoupledSlab:
    def __init__(self,geometry):
        self.geometry=geometry
        self.cfg=geometry.cfg
        box=UnitaryMass(replace(self.cfg,mass=1.,expand=False,subtract=False),None)
        self.box=box
        self.l,self.n=len(box.perp),box.n
        self.nb=1+2*geometry.bath_bins
        self.block=2*self.nb
        self.dim=4*self.n*self.nb
        self.view_shape=(self.l,self.n,2,2,self.nb)
        self.weight=box.weight[:,None,None]
        self.q=np.sqrt(box.q2)
        _,spin=np.linalg.eigh(box.k[None,:,None,None]*SIGMA_Z+box.perp[:,None,None,None]*SIGMA_X)
        self.transform=np.zeros((self.l,4*self.n,4*self.n))
        for j in range(self.n):
            for h in (0,1):
                for b in (0,1):
                    self.transform[:,2*(b*self.n+j):2*(b*self.n+j+1),(2*j+h)*2+b]=spin[:,j,:,h]
        signed=self.q[:,:,None]*np.array([-1.,1.])
        self.k_signed=signed
        self.kinetic=signed[...,None]*np.array([1.,-1.])
        freq=np.empty((*signed.shape,self.nb-1))
        g=np.empty((*signed.shape,2,self.nb-1))
        max_spacing=0.
        for index,k in np.ndenumerate(signed):
            span=np.hypot(k,geometry.bath_cutoff)-abs(k)
            star,spacing=build(float(k),geometry.gamma,geometry.bath_cutoff,
                               span/(geometry.bath_bins-0.25))
            if len(star.frequency)!=self.nb-1:
                raise ValueError('Unexpected bath size')
            freq[index]=star.frequency
            g[index]=star.coupling
            max_spacing=max(max_spacing,spacing)
        self.frequency=freq
        self.bath_frequency=np.broadcast_to(freq[...,None,:],g.shape)
        self.g=g
        self.recurrence=2*np.pi/max_spacing
        self.sterile_indices=np.arange(self.dim).reshape(self.n,2,2,self.nb)[...,0].ravel()
        self.kinetic_blocks=np.zeros((self.l,self.n,2,self.block,self.block))
        for b in (0,1):
            off=b*self.nb
            bath=np.arange(off+1,off+self.nb)
            self.kinetic_blocks[...,off,off]=self.kinetic[...,b]
            self.kinetic_blocks[...,bath,bath]=freq
            self.kinetic_blocks[...,off,bath]=g[...,b,:]
            self.kinetic_blocks[...,bath,off]=g[...,b,:]
        self.mass_unit=np.zeros((self.block,self.block))
        self.mass_unit[0,self.nb]=self.mass_unit[self.nb,0]=1.
        self.phase_unit=np.zeros((self.block,self.block),complex)
        self.phase_unit[0,self.nb]=1j/self.cfg.phase_scale
        self.phase_unit[self.nb,0]=-1j/self.cfg.phase_scale
        self.prepare_reference()

    def mass(self,a):
        m=self.cfg.mass
        if m==0:
            return 0.
        return float(mass_matched(a*m,self.geometry.gamma,self.geometry.bath_cutoff)/a)

    def mass_slope(self,a):
        eps=1e-5
        return (self.mass(a*math.exp(eps))-self.mass(a*math.exp(-eps)))/(2*eps)

    def eigensystem(self,a):
        return np.linalg.eigh(self.kinetic_blocks/a+self.mass(a)*self.mass_unit)

    def sector_measure(self,rows):
        """Unweighted expectation by transverse momentum, longitudinal momentum and helicity."""
        s=rows[...,0,:]
        n=np.sum(abs(s)**2,axis=-1)
        mass_x=2*np.sum((s[...,0,:].conj()*s[...,1,:]).real,axis=-1)
        z=n[...,0]-n[...,1]
        bath=rows[...,1:,:]
        diag=np.sum(abs(bath)**2,axis=-1)
        bath_e=np.sum(diag*self.bath_frequency,axis=(-1,-2))
        cross=2*np.sum((s.conj()*np.einsum('lqhbr,lqhbrc->lqhbc',self.g,bath,optimize=True)).real,axis=(-1,-2))
        kinetic=np.sum(self.kinetic*n,axis=-1)
        plus=np.sum(diag[...,:self.geometry.bath_bins],axis=-1)
        minus=np.sum(diag[...,self.geometry.bath_bins:],axis=-1)
        number=n.sum(axis=-1)
        return n,z,mass_x,kinetic,bath_e,cross,plus,minus,number

    def vacuum_row(self,a):
        e,u=self.eigensystem(a)
        half=self.nb
        if np.any(e[...,half-1]>=0) or np.any(e[...,half]<=0):
            raise ValueError('Vacuum split is not gapped')
        occupied=u[...,:half]
        empty=u[...,half:]
        rows=occupied.reshape(*self.view_shape,half)
        n,z,mx,kin,bath,cross,plus,minus,count=self.sector_measure(rows)
        gap=e[...,half:,None]-e[...,None,:half]
        dm=self.mass_slope(a)
        v_phase=self.mass(a)*self.phase_unit
        v_metric=-self.kinetic_blocks/a+dm*self.mass_unit
        phase=adjoint(empty) @ v_phase @ occupied
        metric=adjoint(empty) @ v_metric @ occupied
        Zt=np.sum(self.weight*np.sum(abs(phase)**2/gap**3,axis=(-1,-2)))/a**3
        Zg=np.sum(self.weight*np.sum(abs(metric)**2/gap**3,axis=(-1,-2)))/a**3
        potential=np.sum(self.weight*np.sum(e[...,:half],axis=-1))
        slope=np.sum(self.weight*(-(kin+bath+cross)/a+dm*mx))
        aux=self.auxiliary(a,z,mx,kin,bath,cross,plus,minus,count)
        return potential,slope,np.r_[Zt,Zg,aux]

    def prepare_reference(self):
        grid=np.linspace(0,np.log(1.55),self.geometry.reference_nodes)
        rows=[self.vacuum_row(float(np.exp(x))) for x in grid]
        self.vac_energy=CubicHermiteSpline(grid,[r[0] for r in rows],[r[1] for r in rows],extrapolate=False)
        self.vac_aux=CubicSpline(grid,np.array([r[2] for r in rows]),extrapolate=False)

    def matching(self,a):
        x=math.log(a)
        if x<-1e-14 or x>math.log(1.55):
            raise ValueError(f'Scale outside matched table: {a}')
        aux=self.vac_aux(max(0.,x))
        A=1-aux[0]
        if not A>0:
            raise ValueError('Nonpositive scalar kinetic coefficient')
        return float(A),float(-self.vac_aux(max(0.,x),1)[0]),float(aux[1]),aux[2:]

    def vacuum(self,a=1.):
        _,u=self.eigensystem(a)
        frame=np.zeros((self.l,self.n,2,self.block,self.dim//2),complex)
        for j in range(self.n):
            for h in (0,1):
                span=slice((2*j+h)*self.nb,(2*j+h+1)*self.nb)
                frame[:,j,h,:,span]=u[:,j,h,:,:self.nb]
        return frame.reshape(self.l,self.dim,self.dim//2)

    def full_hamiltonian(self,a,phi,free=False):
        blocks=self.kinetic_blocks/a+self.mass(a)*self.mass_unit
        out=np.zeros((self.l,self.dim,self.dim),complex)
        for j in range(self.n):
            for h in (0,1):
                section=slice((2*j+h)*self.block,(2*j+h+1)*self.block)
                out[:,section,section]=blocks[:,j,h]
        if free:
            for h in range(self.l):
                ii=self.sterile_indices
                out[h,ii,:]=0.
                out[h,:,ii]=0.
                for j in range(self.n):
                    for k in (0,1):
                        start=(2*j+k)*self.block
                        iq=[start,start+self.nb]
                        free_h=np.diag(self.kinetic[h,j,k]/a)+self.cfg.mass*SIGMA_X
                        out[h,np.ix_(iq,iq)[0],np.ix_(iq,iq)[1]]=free_h
            return out
        rm=np.kron(self.box.residual(phi/self.cfg.phase_scale),np.eye(2))*self.mass(a)
        idx=self.sterile_indices
        out[:,idx[:,None],idx[None,:]]+=adjoint(self.transform) @ rm @ self.transform
        return out

    def quantum_free(self,a,frame,dt):
        rows=frame.reshape(*self.view_shape,-1)
        s=rows[...,0,:]
        bath=rows[...,1:,:]
        z=2j/dt
        g=self.g/a
        d=1/(z-self.bath_frequency/a)
        gd=g*d
        diag=z-self.kinetic/a-np.sum(g*gd,axis=-1)
        target=s+np.einsum('lqhbr,lqhbrc->lqhbc',gd,bath,optimize=True)
        mass=self.mass(a)
        inv=1/(diag[...,0]*diag[...,1]-mass**2)
        sterile=np.stack([diag[...,1,None]*target[...,0,:]+mass*target[...,1,:],
                          mass*target[...,0,:]+diag[...,0,None]*target[...,1,:]],axis=-2)*inv[...,None,None]
        out=np.empty_like(rows)
        out[...,0,:]=2*z*sterile-s
        out[...,1:,:]=2*z*d[...,None]*(bath+g[...,None]*sterile[...,None,:])-bath
        return out.reshape(frame.shape)

    def lab(self,frame):
        lab=self.transform @ frame[:,self.sterile_indices,:]
        return lab.reshape(self.l,2*self.n,2,-1)

    def local(self,a,phi,frame,dt,force):
        m=self.mass(a)
        phase=phi/self.cfg.phase_scale
        basis_phase=self.box.phase_basis(phase)
        energy,basis=np.linalg.eigh(m*self.box.residual(phase,basis_phase))
        propagator=(basis*np.exp(-1j*energy*dt)) @ adjoint(basis)
        lab=self.lab(frame)
        evolved=np.einsum('ij,ljsm->lism',propagator,lab,optimize=True)
        out=frame.copy()
        out[:,self.sterile_indices,:]=adjoint(self.transform) @ evolved.reshape(self.l,4*self.n,-1)
        if not force:
            return out,None
        cov=adjoint(basis) @ self.box.cov(lab) @ basis
        gap=energy[:,None]-energy[None,:]
        integral=dt*np.exp(-0.5j*gap*dt)*np.sinc(gap*dt/(2*np.pi))
        integrated=basis @ (cov*integral) @ adjoint(basis)
        impulse=-m*self.box.project_scalar(self.box.mass_derivative(phase,integrated,basis_phase))
        return out,impulse

    def auxiliary(self,a,z,mx,kin,bath,cross,plus,minus,count):
        E=np.hypot(self.k_signed,a*self.cfg.mass)
        polar=np.divide(self.k_signed*z+a*self.cfg.mass*mx,E,out=np.zeros_like(E),where=E>0)
        occupation=(count+polar)/2
        w=self.weight
        return np.r_[np.sum(w*kin),np.sum(w*bath),np.sum(w*cross),
                     np.sum(w[...,None]*plus,axis=(0,1,2)),
                     np.sum(w[...,None]*minus,axis=(0,1,2)),np.sum(w*occupation)]

    def observe(self,state,source=False):
        rows=state.frame.reshape(*self.view_shape,-1)
        _,z,mx,kin,bath,cross,plus,minus,count=self.sector_measure(rows)
        cov=self.box.cov(self.lab(state.frame))
        phase=state.phi/self.cfg.phase_scale
        phase_basis=self.box.phase_basis(phase)
        mass_unit=self.box.chiral(np.eye(self.n))+self.box.residual(phase,phase_basis)
        mass_value=float(np.einsum('ij,ji->',mass_unit,cov).real)
        force=self.box.project_scalar(self.box.mass_derivative(phase,cov,phase_basis)) if source else None
        momentum=float(np.sum(self.weight*self.box.k[None,:,None]*(count+(plus+minus).sum(axis=-1))))
        return z,mx,kin,bath,cross,plus,minus,count,mass_value,momentum,force

    def measure(self,state,observation=None,source=False):
        a,phi=state.a,state.phi
        A,As,Zg,vac=self.matching(a)
        obs=self.observe(state,source=source) if observation is None else observation
        z,mx,kin,bath,cross,plus,minus,count,mass_value,momentum,force=obs
        aux=self.auxiliary(a,z,mx,kin,bath,cross,plus,minus,count)
        kval=sum(aux[:3])
        m=self.mass(a)
        loga=math.log(a)
        energy_q=(kval/a+m*mass_value-self.vac_energy(loga))/(2*a**3)
        pressure_q=(kval/a-self.mass_slope(a)*mass_value+self.vac_energy(loga,1))/(6*a**3)
        velocity=state.pi/(a**3*A)
        dx=self.box.derivative(phi)
        grad=np.mean(dx**2)/a**2
        v2=np.mean(velocity**2)
        cfg=self.cfg
        V=cfg.vacuum_density+0.5*cfg.scalar_mass**2*np.mean(phi**2)
        density_scalar=V+0.5*A*(v2+grad)
        pressure_scalar=-V+A*v2/2+As*v2/6-(A+As)*grad/6
        rho=float(density_scalar+energy_q)
        if rho<=0:
            raise ValueError('Nonpositive total density')
        H=math.sqrt(rho/(3*cfg.planck**2+Zg/2)) if cfg.expand else 0.
        delta=(aux-vac)/2
        result=dict(a=a,rho=rho,P=float(pressure_scalar+pressure_q),H=H,
                    rho_q=float(energy_q),P_q=float(pressure_q),rho_scalar=float(density_scalar),
                    gradient=float(A*grad/2),bath_energy=float(delta[1]/a),
                    interaction_energy=float(delta[2]/a),daughter_plus=float(delta[3]),
                    daughter_minus=float(delta[4]),bath_holes=float(-sum(delta[5:7])),
                    sterile_count=float(delta[7]),momentum=float(momentum/2-np.mean(state.pi*dx)))
        if source:
            if force is None:
                force=self.observe(state,source=True)[-1]
            result['work']=float(m*np.mean(force*velocity)/a**3)
        return result

    def linear(self,state,dt):
        a=state.a
        A,*_=self.matching(a)
        inertia=a**3*A
        omega=np.sqrt(self.cfg.scalar_mass**2/A+(self.box.scalar_k/a)**2)
        c,s=np.cos(omega*dt),np.sin(omega*dt)
        phi,v=np.fft.fft(state.phi),np.fft.fft(state.pi)/inertia
        return State(a,np.fft.ifft(c*phi+s*v/omega).real,
                     np.fft.ifft(inertia*(c*v-s*omega*phi)).real,
                     self.quantum_free(a,state.frame,dt))

    def expand(self,state,observation,dt):
        if not self.cfg.expand:
            return state,0.
        a=state.a
        H=self.measure(state,observation)['H']
        middle=self.measure(replace(state,a=a*math.exp(H*dt/2)),observation)
        rate=middle['H']*dt
        return replace(state,a=a*math.exp(rate)),a**3*math.expm1(3*rate)*middle['P']

    def step(self,state,observation,dt):
        before,first=self.expand(state,observation,dt/2)
        mid=self.linear(before,dt/2)
        frame,kick=self.local(mid.a,mid.phi,mid.frame,dt,True)
        after=self.linear(replace(mid,frame=frame,pi=mid.pi+kick),dt/2)
        obs=self.observe(after,source=True)
        end,last=self.expand(after,obs,dt/2)
        return end,obs,first+last

    def prepare(self,dt):
        frame=self.vacuum()
        x=2*np.pi*self.box.x/self.cfg.length
        duration=1.5
        steps=round(duration/dt)
        delta=duration/steps
        for i in range(steps):
            t=-duration+(i+0.5)*delta
            u=(t+duration)/duration
            ramp=u**3*(10+u*(-15+6*u))
            phase=ramp*(0.15*np.sin(x)+t*(0.9*np.cos(x)+0.30))
            frame=self.quantum_free(1.,frame,delta/2)
            frame,_=self.local(1.,phase,frame,delta,False)
            frame=self.quantum_free(1.,frame,delta/2)
        A,*_=self.matching(1.)
        return State(1.,0.15*np.sin(x),A*(0.9*np.cos(x)+0.30),frame)
