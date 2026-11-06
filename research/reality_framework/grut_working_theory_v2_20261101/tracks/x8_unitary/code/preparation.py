"""Third-order ultraviolet graph initialization for the paired slab."""
from dataclasses import replace
import math

import numpy as np

from paired import PairedSlab, apply
from unitary import adjoint


class ThirdOrderSlab(PairedSlab):
    def __init__(self,cfg,order,transition):
        super().__init__(cfg,order)
        self.transition=transition

    def tangent(self,state,obs):
        A,*_=self.paired_coefficients(state.a)
        v=apply(1/A,state.pi)/state.a**3
        d=self.measure(state,obs,work=True)
        return v,d['H'],d['Hdot'],d['acceleration']

    def eigensystem(self,a,phi):
        return np.linalg.eigh(self.spatial.full_hamiltonian(a,phi/self.cfg.phase_scale))

    def second_projector(self,state,obs):
        box=self.spatial
        d=box.d
        v,H,Hdot,acceleration=self.tangent(state,obs)
        e,basis=self.eigensystem(state.a,state.phi)
        plus,minus=basis[:,:,d:],basis[:,:,:d]
        gap,x1=box.first_graph(state.a,state.phi,v,H,e,basis)

        def first(t):
            a=state.a*math.exp(H*t+0.5*Hdot*t*t)
            phi=state.phi+v*t+0.5*acceleration*t*t
            e,basis=self.eigensystem(a,phi)
            _,x=box.first_graph(a,phi,v+acceleration*t,H+Hdot*t,e,basis)
            part=basis[:,:,d:] @ x @ adjoint(basis[:,:,:d])
            return part+adjoint(part)

        eps=3e-5
        d1=(first(eps)-first(-eps))/(2*eps)
        x2=1j*(adjoint(plus) @ d1 @ minus)/gap
        crossing=plus @ x2 @ adjoint(minus)
        occupied=minus @ adjoint(x1)
        empty=plus @ x1
        diagonal=empty @ adjoint(empty)-occupied @ adjoint(occupied)
        return crossing+adjoint(crossing)+diagonal

    def initial(self,phi,speed):
        state=super().initial(phi,speed)
        box=self.spatial
        d=box.d
        obs=self.observe(state)
        v,H,*_=self.tangent(state,obs)
        e,basis=self.eigensystem(state.a,state.phi)
        plus,minus=basis[:,:,d:],basis[:,:,:d]
        gap,x1=box.first_graph(state.a,state.phi,v,H,e,basis)
        eps=1e-4
        nearby=[]
        for sign in (-1,1):
            point,expectation,_=self.step(state,obs,sign*eps)
            nearby.append(self.second_projector(point,expectation))
        dx2=(nearby[1]-nearby[0])/(2*eps)
        x3=1j*(adjoint(plus) @ dx2 @ minus)/gap
        cubic=x3+x1 @ (adjoint(x1) @ x1)
        o=state.orbitals[0].reshape(len(box.perp),2*d,d)
        lower=adjoint(minus) @ o
        upper=adjoint(plus) @ o
        graph=np.linalg.solve(lower.swapaxes(-1,-2),upper.swapaxes(-1,-2)).swapaxes(-1,-2)
        lo,hi=(p*self.cfg.scalar_mass for p in self.transition)
        x=np.clip((gap-lo)/(hi-lo),0,1)
        weight=x*x*(3-2*x)
        final=box.graph_orbitals(basis,graph+weight*cubic)
        return replace(state,orbitals=(final,))
