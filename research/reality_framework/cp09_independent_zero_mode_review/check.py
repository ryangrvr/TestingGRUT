"""Independent N=2 CP-09 zero-mode and empty-state spectral-moment controls.

Does not import the frozen CP-09 implementation. Dependencies: numpy.
"""
import numpy as np
from numpy.linalg import eigvalsh, norm
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],dtype=complex)
Z=np.diag([1.,-1.]).astype(complex)
K=np.kron
J=1.
def build(gamma,Omega):
    P01=K(Z,X); P10=K(X,Z)
    U01=np.cos(gamma)*np.eye(4)-1j*np.sin(gamma)*P01
    U10=np.cos(gamma)*np.eye(4)-1j*np.sin(gamma)*P10
    H=np.zeros((12,12),complex)
    onsite=Omega*(2*np.eye(4)-K(X,I)-K(I,X))
    # Geometry basis: [empty, 0<1, 1<0]. Each link kinetic is L^dagger L.
    for q,U in ((1,U01),(2,U10)):
        sl=slice(4*q,4*q+4)
        H[:4,:4]+=J*np.eye(4)
        H[sl,sl]+=J*np.eye(4)
        H[sl,:4]+=-J*U
        H[:4,sl]+=-J*U.conj().T
    for q in range(3):
        H[4*q:4*q+4,4*q:4*q+4]+=onsite
    assert norm(H-H.conj().T)<1e-12
    return H,U01,U10
def test():
    gam=np.pi/4;om=.25
    for g,o in ((0.,om),(gam,0.),(gam,om)):
        H,U,V=build(g,o)
        ev=eigvalsh(H)
        print(f"gamma={g:.8f} Omega={o:.4f}: E0={ev[0]:.12f}, nullity={sum(abs(ev)<1e-9)}")
        assert ev[0]>-1e-12
        if g==0 or o==0:
            assert abs(ev[0])<1e-12
        else:
            assert ev[0]>1e-5
            assert abs(ev[0]-.09475528369770218)<1e-11
    H,U,V=build(gam,om)
    v=np.zeros(12,complex);v[0]=1
    mean=np.vdot(v,H@v).real
    variance=np.vdot(H@v,H@v).real-mean**2
    assert abs(mean-(2*J+2*om))<1e-12
    assert abs(variance-(2*J**2+2*om**2))<1e-12
    print(f"empty |00> N=2: mean={mean:.12f}, variance={variance:.12f}")
    # Exact zero-mode impossibility: unique onsite kernel |++>,
    # but U_ij |++> has orthogonal |-+> component when sin gamma != 0.
    plus=np.array([1,1],complex)/np.sqrt(2)
    plusplus=K(plus,plus)
    for U in (U,V):
        overlap=np.vdot(plusplus,U@plusplus)
        off=norm(U@plusplus-overlap*plusplus)
        assert abs(abs(overlap)-abs(np.cos(gam)))<1e-12
        assert abs(off-abs(np.sin(gam)))<1e-12
    print("PASS: Ω=0 N=2 flat connection has zero modes; Ω>0 γ!=kπ has none")
    print("PASS: exact empty-state mean/variance for N=2")
if __name__=="__main__":
    test()
