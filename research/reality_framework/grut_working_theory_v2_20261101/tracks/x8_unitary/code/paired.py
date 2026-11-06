"""Common action subtraction on the same momentum pairs as the unitary bubble."""
from dataclasses import replace
import math

import numpy as np

from unitary import State, UnitarySlab, adjoint


def apply(diagonal,profile):
    return np.fft.ifft(diagonal*np.fft.fft(profile)).real


class PairedSlab(UnitarySlab):
    def paired_coefficients(self,a):
        cfg,box=self.cfg,self.spatial
        p2=self.q2/a**2
        e=np.sqrt(cfg.mass**2+p2)
        base=cfg.mass**2/(4*a**3)
        z=np.zeros(cfg.points)
        slope=np.zeros(cfg.points)
        e3=box.weight @ e**-3
        e5=box.weight @ e**-5
        orders=np.rint(abs(box.scalar_k)*cfg.length/(2*np.pi)).astype(int)
        for r in range(cfg.scalar_modes+1):
            # Reflective free energies make the upper and lower sums equal.
            last=box.n-r
            z[orders==r]=base*np.sum(e3[:last])/cfg.phase_scale**2
            slope[orders==r]=3*base*cfg.mass**2*np.sum(e5[:last])/cfg.phase_scale**2
        zg=base*np.sum(box.weight[:,None]*p2/e**5)
        zg_slope=-5*base*cfg.mass**2*np.sum(box.weight[:,None]*p2/e**7)
        if np.any(z>=1):
            raise ValueError('Nonpositive scalar kinetic coefficient')
        return 1-z,slope,float(zg),float(zg_slope)

    def measure(self,state,expectation,work=False):
        a,phi,cfg=state.a,state.phi,self.cfg
        box=self.spatial
        A,slope,zg,zg_slope=self.paired_coefficients(a)
        v=apply(1/A,state.pi)/a**3
        uniform_a,*_=self.matching(a)
        raw=super().measure(replace(state,pi=a**3*uniform_a*v),expectation)
        dx=box.derivative(phi)
        gradient=dx/a
        potential=cfg.vacuum_density+0.5*cfg.scalar_mass**2*np.mean(phi**2)
        kinetic=0.5*np.mean(v*apply(A,v))
        grad_energy=0.5*np.mean(gradient*apply(A,gradient))
        pressure_scalar=kinetic-potential+np.mean(v*apply(slope,v)-gradient*apply(A+slope,gradient))/6
        density_scalar=kinetic+potential+grad_energy
        rho=density_scalar+raw['rho_F']
        P=pressure_scalar+raw['P_F']
        if rho<=0 or 3*cfg.planck**2+zg/2<=0:
            raise ValueError('Invalid density or gravitational matching')
        H=math.sqrt(rho/(3*cfg.planck**2+zg/2)) if cfg.expand else 0.
        Z=1-A
        rho_ct=0.5*np.mean(v*apply(Z,v)+gradient*apply(Z,gradient))
        pressure_ct=np.mean(v*apply(Z/2-slope/6,v)-gradient*apply((Z-slope)/6,gradient))
        raw.update(rho=rho,P=P,H=H,A=float(np.min(A)),Z_g=zg,
                   rho_F_matched=raw['rho_F']-rho_ct,P_F_matched=raw['P_F']-pressure_ct,
                   P_x=pressure_scalar+4*grad_energy/3+raw['P_F_x'],
                   P_perp=pressure_scalar-2*grad_energy/3+raw['P_F_perp'],
                   rho_scalar=density_scalar,scalar_gradient=grad_energy,
                   momentum=(expectation.comoving_momentum-np.mean(state.pi*dx))/a**4)
        if work:
            J=expectation.force/a**3
            laplace=box.derivative(dx)/a**2
            acceleration=laplace+apply(1/A,-cfg.scalar_mass**2*phi-J-H*apply(3*A+slope,v))
            matched=J-apply(Z,acceleration)-H*apply(3*Z-slope,v)+apply(Z,laplace)
            hdot=-(3*(rho+P)+zg_slope*H*H/2)/(6*cfg.planck**2+zg) if cfg.expand else 0.
            raw.update(Jv=np.mean(v*J),Jv_matched=np.mean(v*matched),
                       acceleration=acceleration,Hdot=hdot)
        return raw

    def linear(self,state,dt):
        cfg,a=self.cfg,state.a
        A,*_=self.paired_coefficients(a)
        inertia=a**3*A
        omega=np.sqrt(cfg.scalar_mass**2/A+(self.spatial.scalar_k/a)**2)
        cos,sin=np.cos(omega*dt),np.sin(omega*dt)
        phi,v=np.fft.fft(state.phi),np.fft.fft(state.pi)/inertia
        occupied=self.spatial.free_step(a,state.orbitals[0],dt)
        return State(a,np.fft.ifft(cos*phi+sin*v/omega).real,
                     np.fft.ifft(inertia*(cos*v-sin*omega*phi)).real,(occupied,))

    def initial(self,phi,speed):
        cfg,box=self.cfg,self.spatial
        phi,speed=box.project_scalar(phi),box.project_scalar(speed)
        A,*_=self.paired_coefficients(1.)
        pi=apply(A,speed)
        h=box.full_hamiltonian(1.,phi/cfg.phase_scale)
        energy,basis=np.linalg.eigh(h)
        if np.any(energy[:,box.d-1]>=0) or np.any(energy[:,box.d]<=0):
            raise ValueError('Prepared state has no occupied/unoccupied gap')
        grad=box.derivative(phi)
        classical=cfg.vacuum_density+0.5*np.mean(speed**2+cfg.scalar_mass**2*phi**2+grad**2)
        hubble=math.sqrt(classical/(3*cfg.planck**2)) if cfg.expand else 0.
        _,mix=box.first_graph(1.,phi,speed,hubble,energy,basis)
        state=State(1.,phi,pi,(box.graph_orbitals(basis,mix),))
        result=self.measure(state,self.observe(state),work=True)
        H,Hdot,acceleration=(result[k] for k in ('H','Hdot','acceleration'))
        gap,mix=box.first_graph(1.,phi,speed,H,energy,basis)

        def embedded(t):
            a=math.exp(H*t+0.5*Hdot*t*t)
            field=phi+speed*t+0.5*acceleration*t*t
            energy,basis=np.linalg.eigh(box.full_hamiltonian(a,field/cfg.phase_scale))
            _,mix=box.first_graph(a,field,speed+acceleration*t,H+Hdot*t,energy,basis)
            block=basis[:,:,box.d:] @ mix @ adjoint(basis[:,:,:box.d])
            return block+adjoint(block)

        eps=2e-5
        tangent=(embedded(eps)-embedded(-eps))/(2*eps)
        mix+=1j*(adjoint(basis[:,:,box.d:]) @ tangent @ basis[:,:,:box.d])/gap
        return replace(state,orbitals=(box.graph_orbitals(basis,mix),))
