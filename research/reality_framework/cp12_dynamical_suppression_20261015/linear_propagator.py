"""Reuse real sparse H rather than storing real-shifted and complex copies."""
import numpy as np
from scipy.sparse.linalg import LinearOperator
from boundary import RF
from integrals import Propagator,schur_norm


def shifted_operators(H,shift):
    def sparse_action(matrix,x):
        if np.iscomplexobj(x):
            # Avoid an implicit complex copy of every real sparse entry.
            return (matrix@np.ascontiguousarray(x.real)+
                1j*(matrix@np.ascontiguousarray(x.imag)))
        return matrix@x
    def real_action(x):return sparse_action(H,x)-shift*x
    def real_adjoint(x):return sparse_action(H.T,x)-shift*x
    def generator(x):return -1j*real_action(x)
    def adjoint_generator(x):return 1j*real_adjoint(x)
    A=LinearOperator(H.shape,matvec=real_action,matmat=real_action,
        rmatvec=real_adjoint,rmatmat=real_adjoint,dtype=np.float64)
    G=LinearOperator(H.shape,matvec=generator,matmat=generator,
        rmatvec=adjoint_generator,rmatmat=adjoint_generator,dtype=np.complex128)
    trace=-1j*complex(H.diagonal().sum()-shift*H.shape[0])
    return A,G,trace


class MatrixFreePropagator(Propagator):
    def __init__(self,H,actions,shift,dt=1/32,batch=2):
        self.A,self.generator,self.trace=shifted_operators(H,shift)
        self.actions=actions;self.norms={name:schur_norm(L) for name,L in actions.items()}
        self.dt=dt;self.batch=batch
