"""A band-projected mass phase with an exactly unitary finite vertex."""
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2] / 'x7_regulator/code'))
from evolution import RegulatedSlab, State
from regulator import Configuration, RadialSlab
from slab import adjoint


class UnitaryMass(RadialSlab):
    def phase_basis(self,phase):
        return np.linalg.eigh(self.matrix(phase))

    def chiral(self,off):
        block=np.zeros((self.d,self.d),complex)
        block[:self.n,self.n:]=off
        block[self.n:,:self.n]=adjoint(off)
        return block

    def residual(self,phase,basis=None):
        if self.cfg.mass==0:
            return np.zeros((self.d,self.d),complex)
        angle,u=self.phase_basis(phase) if basis is None else basis
        return self.chiral(self.cfg.mass*(u*np.expm1(1j*angle)) @ adjoint(u))

    @staticmethod
    def frechet(angle):
        left,right=angle[:,None],angle[None,:]
        return 1j*np.exp(0.5j*(left+right))*np.sinc((left-right)/(2*np.pi))

    def mass_derivative(self,phase,covariance,basis=None):
        if self.cfg.mass==0:
            return np.zeros_like(phase)
        angle,u=self.phase_basis(phase) if basis is None else basis
        crossing=adjoint(u) @ covariance[self.n:,:self.n] @ u
        adjoint_vertex=u @ (self.frechet(angle).T*crossing) @ adjoint(u)
        local=np.sum((self.plane @ adjoint_vertex)*self.plane.conj(),axis=1)
        return self.cfg.points*self.cfg.mass/self.cfg.phase_scale*local.real

    def first_graph(self,a,phi,speed,hubble,energy,basis):
        phase=phi/self.cfg.phase_scale
        angle,u=self.phase_basis(phase)
        variation=adjoint(u) @ self.matrix(speed/self.cfg.phase_scale) @ u
        tangent=self.cfg.mass*u @ (self.frechet(angle)*variation) @ adjoint(u)
        velocity=np.kron(self.chiral(tangent),np.eye(2))
        derivative=velocity-hubble*self.full_hamiltonian(a,phase,kinetic_only=True)
        gap=energy[:,self.d:,None]-energy[:,None,:self.d]
        mix=-1j*(adjoint(basis[:,:,self.d:]) @ derivative @ basis[:,:,:self.d])/gap**2
        return gap,mix

    def local(self,phase,occupied,dt,force):
        if self.cfg.mass==0:
            return occupied,np.zeros_like(phase) if force else None
        phase_basis=self.phase_basis(phase)
        energy,basis=np.linalg.eigh(self.residual(phase,phase_basis))
        propagator=(basis*np.exp(-1j*energy*dt)) @ adjoint(basis)
        result=np.einsum('ij,ljsm->lism',propagator,occupied,optimize=True)
        if not force:
            return result,None
        cov=adjoint(basis) @ self.cov(occupied) @ basis
        gap=energy[:,None]-energy[None,:]
        integral=dt*np.exp(-0.5j*gap*dt)*np.sinc(gap*dt/(2*np.pi))
        integrated=basis @ (cov*integral) @ adjoint(basis)
        impulse=-self.project_scalar(self.mass_derivative(phase,integrated,phase_basis))
        return result,impulse


class UnitarySlab(RegulatedSlab):
    """Reuse the X7 kinetic, stress, preparation and expansion ledger with one mass."""
    def __init__(self,cfg,order):
        self.cfg=cfg
        self.masses=np.array([cfg.mass])
        self.signed=np.ones(1)
        self.spatial=UnitaryMass(replace(cfg,expand=False,subtract=False),order)
        self.sectors=[self.spatial]
        self.weight=self.spatial.weight
        self.qx2=self.spatial.k**2
        self.qperp2=self.spatial.perp[:,None]**2
        self.q2=self.spatial.q2
        self.mass2=self.masses[:,None,None]**2
        self.signed_weight=self.weight[None,:,None]
