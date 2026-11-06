#!/usr/bin/env python3
"""Finite-regulator four-component fermions with one resolved spatial coordinate."""
from dataclasses import dataclass, replace
import math

import numpy as np

SIGMA_X = np.array([[0, 1], [1, 0]])
SIGMA_Z = np.diag([1, -1])
AX = np.kron(SIGMA_Z, SIGMA_Z)
AP = np.kron(SIGMA_Z, SIGMA_X)
BETA = np.kron(SIGMA_X, np.eye(2))
IDENTITY4 = np.eye(4)


def adjoint(a):
    return a.swapaxes(-1, -2).conj()


@dataclass(frozen=True)
class Configuration:
    longitudinal: int = 6
    transverse_cut: float = 8.0
    transverse_bins: int = 16
    points: int = 49
    scalar_modes: int = 3
    mass: float = 0.7
    phase_scale: float = 1.15
    scalar_mass: float = 2.6
    planck: float = 16.0
    vacuum_density: float = 1.5
    expand: bool = True
    subtract: bool = True
    length: float = 2 * math.pi


@dataclass
class State:
    a: float
    phi: np.ndarray
    pi: np.ndarray
    occupied: np.ndarray


class Slab:
    def __init__(self, cfg):
        self.cfg = cfg
        if cfg.longitudinal < 1 or cfg.points <= 4*cfg.longitudinal or cfg.points <= 2*cfg.scalar_modes:
            raise ValueError('Quadrature must resolve every Galerkin mode difference')
        self.k = 2*math.pi/cfg.length * np.arange(-cfg.longitudinal, cfg.longitudinal+1)
        self.n = len(self.k)
        self.d = 2*self.n
        edges = np.linspace(0, cfg.transverse_cut, cfg.transverse_bins+1)
        self.perp = (edges[:-1]+edges[1:])/2
        self.weight = np.diff(edges**2)/(4*math.pi*cfg.length)
        self.x = cfg.length*np.arange(cfg.points)/cfg.points
        self.plane = np.exp(1j*self.x[:,None]*self.k[None,:])/math.sqrt(cfg.points)
        index = np.fft.fftfreq(cfg.points)*cfg.points
        self.scalar_k = 2*math.pi/cfg.length*index
        self.scalar_mask = abs(index) <= cfg.scalar_modes+1e-8
        self.q2 = self.k[None,:]**2+self.perp[:,None]**2
        self.kx = self.k[None,:,None,None]*AX
        self.kt = self.perp[:,None,None,None]*AP

    def project_scalar(self, f):
        return np.fft.ifft(np.fft.fft(f)*self.scalar_mask).real

    def derivative(self, f):
        return np.fft.ifft(1j*self.scalar_k*np.fft.fft(f)).real

    def blocks(self, occupied):
        return occupied.reshape(len(self.perp),2,self.n,2,self.d).transpose(0,2,1,3,4).reshape(
            len(self.perp),self.n,4,self.d)

    def orbitals(self, blocks):
        return blocks.reshape(len(self.perp),self.n,2,2,self.d).transpose(0,2,1,3,4).reshape(
            len(self.perp),self.d,2,self.d)

    def energy(self, a):
        return np.sqrt(self.cfg.mass**2+self.q2/a**2)

    def free_blocks(self, a):
        return (self.kx+self.kt)/a+self.cfg.mass*BETA

    def matrix(self, profile):
        """Local multiplication with the same projection on both sides."""
        return adjoint(self.plane) @ (profile[:,None]*self.plane)

    def mass_matrix(self, profile):
        off = self.matrix(profile)
        out = np.zeros((self.d,self.d),complex)
        out[:self.n,self.n:] = off
        out[self.n:,:self.n] = adjoint(off)
        return out

    def residual(self, phase):
        return self.mass_matrix(self.cfg.mass*np.expm1(1j*phase))

    def full_hamiltonian(self, a, phase, kinetic_only=False):
        out = np.zeros((len(self.perp),2*self.d,2*self.d),complex)
        free = (self.kx+self.kt)/a if kinetic_only else self.free_blocks(a)
        for j in range(self.n):
            ind = np.array([2*j,2*j+1,2*(j+self.n),2*(j+self.n)+1])
            out[:,ind[:,None],ind[None,:]] = free[:,j]
        if not kinetic_only:
            out += np.kron(self.residual(phase),np.eye(2))
        return out

    def cov(self, occupied):
        weighted = (occupied*np.sqrt(self.weight)[:,None,None,None]).transpose(1,0,2,3)
        amplitudes = weighted.reshape(self.d,-1)
        return amplitudes @ adjoint(amplitudes)

    def mass_derivative(self, phase, covariance):
        cross = self.plane @ covariance[self.n:,:self.n]
        local = np.sum(cross*self.plane.conj(),axis=1)
        return self.cfg.points*(1j*self.cfg.mass/self.cfg.phase_scale*np.exp(1j*phase)*local).real

    def spatial_match(self, a):
        cfg = self.cfg
        if not cfg.subtract or cfg.mass == 0:
            return 0.0
        energy = self.energy(a)
        dot = cfg.mass**2 + (self.perp[:,None]**2+self.k[:-1]*self.k[1:])/a**2
        lo, hi = energy[:,:-1], energy[:,1:]
        virtual = (1+dot/(lo*hi))/(lo+hi)
        curvature = (1/energy).sum(axis=1)-virtual.sum(axis=1)
        k_ref = 2*math.pi/cfg.length
        return cfg.mass**2*(self.weight @ curvature)/(cfg.phase_scale**2*a*k_ref**2)

    def matching(self, a, slopes=True):
        cfg = self.cfg
        if not cfg.subtract or cfg.mass == 0:
            return 1.0, 1.0, 0.0, 0.0, 0.0
        e = self.energy(a)
        factor = cfg.mass**2/(4*a**3)
        zt = factor*(self.weight @ (e**-3).sum(axis=1))/cfg.phase_scale**2
        zg = factor*(self.weight @ (self.q2/a**2*e**-5).sum(axis=1))
        zx = self.spatial_match(a)
        at = 3*factor*cfg.mass**2*(self.weight @ (e**-5).sum(axis=1))/cfg.phase_scale**2
        da = 2e-5
        bx = ((self.spatial_match(a*math.exp(-da))-self.spatial_match(a*math.exp(da)))/(2*da)
              if slopes else 0.0)
        if zt>=1 or zx>=1:
            raise ValueError(f'Nonpositive matched kinetic coefficient at a={a}: {zt}, {zx}')
        return 1-zt, 1-zx, at, bx, zg

    def first_graph(self, a, phi, speed, hubble, energy, basis):
        cfg=self.cfg
        angle=phi/cfg.phase_scale
        vertex=self.mass_matrix(1j*cfg.mass*np.exp(1j*angle)*speed/cfg.phase_scale)
        derivative=np.kron(vertex,np.eye(2))-hubble*self.full_hamiltonian(a,angle,kinetic_only=True)
        gap=energy[:,self.d:,None]-energy[:,None,:self.d]
        mix=-1j*(adjoint(basis[:,:,self.d:]) @ derivative @ basis[:,:,:self.d])/gap**2
        return gap,mix

    def graph_orbitals(self, basis, mix):
        vals,axes=np.linalg.eigh(np.eye(self.d)+adjoint(mix) @ mix)
        inverse_root=(axes*vals[:,None,:]**-0.5) @ adjoint(axes)
        occupied=(basis[:,:,:self.d]+basis[:,:,self.d:] @ mix) @ inverse_root
        return occupied.reshape(len(self.perp),self.d,2,self.d)

    def initial(self, phi, speed, a=1.0, adiabatic=True):
        phi, speed = self.project_scalar(phi), self.project_scalar(speed)
        cfg = self.cfg
        angle = phi/cfg.phase_scale
        h = self.full_hamiltonian(a,angle)
        energy, basis = np.linalg.eigh(h)
        if np.any(energy[:,self.d-1]>=0) or np.any(energy[:,self.d]<=0):
            raise ValueError('Prepared state has no occupied/unoccupied gap')
        A,B,At,_,zg = self.matching(a)
        shape=(len(self.perp),self.d,2,self.d)
        pi=a**3*A*speed
        occupied=basis[:,:,:self.d].reshape(shape)
        if adiabatic:
            vclass = cfg.vacuum_density+0.5*np.mean(speed**2+cfg.scalar_mass**2*phi**2
                                                   +(self.derivative(phi)/a)**2)
            hubble = math.sqrt(vclass/(3*cfg.planck**2)) if cfg.expand else 0.0
            gap,mix=self.first_graph(a,phi,speed,hubble,energy,basis)
            occupied=self.graph_orbitals(basis,mix)
            if adiabatic==2:
                state=State(a,phi,pi,occupied)
                measures=self.measure(state)
                hubble=measures['H']
                force=self.project_scalar(self.mass_derivative(angle,self.cov(occupied)))/a**3
                laplace=self.derivative(self.derivative(phi))/a**2
                acceleration=(B*laplace-cfg.scalar_mass**2*phi-force-hubble*(3*A+At)*speed)/A
                zg_slope=-5*cfg.mass**4*(self.weight @ (self.q2/a**2*self.energy(a)**-7).sum(axis=1))/(4*a**3)
                if not cfg.subtract:
                    zg_slope=0.0
                hdot=-(3*(measures['rho']+measures['P'])+zg_slope*hubble**2/2)/(6*cfg.planck**2+zg)
                if not cfg.expand:
                    hdot=0.0
                gap,mix=self.first_graph(a,phi,speed,hubble,energy,basis)

                def correction(t):
                    scale=a*math.exp(hubble*t+hdot*t*t/2)
                    field=phi+speed*t+acceleration*t*t/2
                    vel=speed+acceleration*t
                    e,b=np.linalg.eigh(self.full_hamiltonian(scale,field/cfg.phase_scale))
                    _,x=self.first_graph(scale,field,vel,hubble+hdot*t,e,b)
                    block=b[:,:,self.d:] @ x @ adjoint(b[:,:,:self.d])
                    return block+adjoint(block)

                dt=2e-5
                derivative=(correction(dt)-correction(-dt))/(2*dt)
                mix+=1j*(adjoint(basis[:,:,self.d:]) @ derivative @ basis[:,:,:self.d])/gap
                occupied=self.graph_orbitals(basis,mix)
        return State(a,phi,pi,occupied)

    def measure(self, state, work=False):
        cfg, a, phi, pi = self.cfg, state.a, state.phi, state.pi
        A, B, At, Bt, zg = self.matching(a)
        v = pi/(a**3*A)
        gradient = self.derivative(phi)
        potential = cfg.vacuum_density+0.5*cfg.scalar_mass**2*np.mean(phi**2)
        grad_energy = 0.5*B*np.mean(gradient**2)/a**2
        kinetic = 0.5*A*np.mean(v**2)
        scalar_energy = kinetic+potential+grad_energy
        scalar_p = kinetic-potential+At*np.mean(v**2)/6-(B+Bt)*np.mean(gradient**2)/(6*a**2)
        scalar_momentum = -A*np.mean(v*gradient)/a

        blocks = self.blocks(state.occupied)
        free = self.free_blocks(a)
        e = self.energy(a)
        vacuum = 0.5*(IDENTITY4-free/e[...,None,None])
        delta = blocks @ adjoint(blocks) - vacuum

        def trace(matrix):
            return self.weight @ np.einsum('...ij,...ji->...',matrix,delta).real.sum(axis=1)

        cov = self.cov(state.occupied)
        phase = phi/cfg.phase_scale
        residual_energy = np.trace(self.residual(phase) @ cov).real
        fermion_energy = (trace(free)+residual_energy)/(2*a**3)
        px = trace(self.kx/a)/(2*a**3)
        pt = trace(self.kt/a)/(4*a**3)
        fermion_p = (px+2*pt)/3
        delta_count = np.trace(delta,axis1=-2,axis2=-1).real
        fermion_momentum = (self.weight @ (delta_count @ self.k))/(2*a**4)
        rho = scalar_energy+fermion_energy
        pressure = scalar_p+fermion_p
        if rho<=0:
            raise ValueError('Nonpositive total density in Friedmann equation')
        hubble = math.sqrt(rho/(3*cfg.planck**2+zg/2)) if cfg.expand else 0.0
        zt,zx=1-A,1-B
        speed2=np.mean(v**2)
        grad2=np.mean(gradient**2)/a**2
        rho_ct=(zt*speed2+zx*grad2)/2
        p_ct=(zt/2-At/6)*speed2-(zx-Bt)*grad2/6
        power=power_matched=0.0
        if work:
            force=self.project_scalar(self.mass_derivative(phase,cov))/a**3
            laplace=self.derivative(gradient)/a**2
            acceleration=(B*laplace-cfg.scalar_mass**2*phi-force-hubble*(3*A+At)*v)/A
            matched_force=force-zt*acceleration-hubble*(3*zt-At)*v+zx*laplace
            power=np.mean(v*force)
            power_matched=np.mean(v*matched_force)
        # Directional scalar terms treat coefficients as functions of the volume scale.
        return dict(rho=float(rho), H=float(hubble), P=float(pressure), rho_F=float(fermion_energy),
                    P_F=float(fermion_p), P_F_x=float(px), P_F_perp=float(pt), Jv=float(power),
                    rho_F_matched=float(fermion_energy-rho_ct), P_F_matched=float(fermion_p-p_ct),
                    Jv_matched=float(power_matched),
                    P_x=float(scalar_p+4*grad_energy/3+px), P_perp=float(scalar_p-2*grad_energy/3+pt),
                    momentum=float(scalar_momentum+fermion_momentum), fermion_momentum=float(fermion_momentum),
                    rho_scalar=float(scalar_energy), scalar_gradient=float(grad_energy),
                    A=float(A), B=float(B), Z_g=float(zg))

    def free_step(self, a, occupied, dt):
        energy = self.energy(a)[...,None,None]
        phase = energy*dt
        u = np.cos(phase)*IDENTITY4-1j*np.sin(phase)/energy*self.free_blocks(a)
        return self.orbitals(u @ self.blocks(occupied))

    def linear(self, state, dt):
        cfg, a = self.cfg, state.a
        A, B, _, _, _ = self.matching(a,slopes=False)
        omega = np.sqrt((cfg.scalar_mass**2+B*(self.scalar_k/a)**2)/A)
        cos = np.cos(omega*dt)
        sin = np.sin(omega*dt)
        q = np.fft.fft(state.phi)
        v = np.fft.fft(state.pi)/(a**3*A)
        phi = np.fft.ifft(cos*q+sin*v/omega).real
        pi = np.fft.ifft(cos*v-sin*omega*q).real*(a**3*A)
        return State(a,phi,pi,self.free_step(a,state.occupied,dt))

    def local(self, phase, occupied, dt, force):
        residual = self.residual(phase)
        energy, basis = np.linalg.eigh(residual)
        phase_u = np.exp(-1j*energy*dt)
        u = (basis*phase_u[None,:]) @ adjoint(basis)
        result = np.einsum('ij,ljsm->lism',u,occupied,optimize=True)
        if not force:
            return result, None
        cov = adjoint(basis) @ self.cov(occupied) @ basis
        gap = energy[:,None]-energy[None,:]
        integral = dt*np.exp(-0.5j*gap*dt)*np.sinc(gap*dt/(2*math.pi))
        integrated_cov = basis @ (cov*integral) @ adjoint(basis)
        impulse = -self.project_scalar(self.mass_derivative(phase,integrated_cov))
        return result, impulse

    def matter(self, state, dt):
        middle = self.linear(state,dt/2)
        occupied, impulse = self.local(middle.phi/self.cfg.phase_scale,middle.occupied,dt,True)
        return self.linear(State(middle.a,middle.phi,middle.pi+impulse,occupied),dt/2)

    def expand(self, state, dt):
        if not self.cfg.expand:
            return state, np.zeros(3)
        a0 = state.a
        start = self.measure(state)
        mid = replace(state,a=a0*math.exp(start['H']*dt/2))
        middle = self.measure(mid)
        a1 = a0*math.exp(middle['H']*dt)
        volume_change = a0**3*math.expm1(3*middle['H']*dt)
        return replace(state,a=a1), volume_change*np.array([middle['P'],middle['P_F'],middle['P_F_matched']])

    def step(self, state, dt):
        before, first = self.expand(state,dt/2)
        after = self.matter(before,dt)
        state, second = self.expand(after,dt/2)
        return state, first+second
