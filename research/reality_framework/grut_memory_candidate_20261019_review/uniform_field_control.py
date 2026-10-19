# Independent re-derivation of GRUT-M1 core numbers (L=4) plus one control the packet lacks:
# uniform small transverse field (F_i == eps for all i, no facilitation).
import numpy as np
from scipy.linalg import eigh
I=np.eye(2); X=np.array([[0,1],[1,0.]]); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1.,-1])
def op(L,ops):
    out=np.array([[1.]])
    for k in range(2*L): out=np.kron(out, ops.get(k,I))
    return out
def H_M1(L,J=1,a=.5,hx=1.05,hz=.5,K=1,D=1,eps=.05,facil=True):
    H=np.zeros((4**L,4**L))
    for i in range(L):
        j=(i+1)%L; p=(i-1)%L
        H+=J*op(L,{i:Z,j:Z})+a*J*op(L,{i:Z,j:Z,L+i:Z})+hx*op(L,{i:X})+hz*op(L,{i:Z})-K*op(L,{L+i:Z,L+j:Z})
        if facil:
            H-=D*(1+eps)/2*op(L,{L+i:X}); H+=D*(1-eps)/2*op(L,{L+p:Z,L+i:X,L+j:Z})
        else:
            H-=D*eps*op(L,{L+i:X})
    return H
def prod(vs):
    out=np.array([1.+0j])
    for v in vs: out=np.kron(out,v)
    return out
px=np.array([1,1])/np.sqrt(2); up=np.array([1.,0]); dn=np.array([0,1.])
def run(L,facil,label):
    H=H_M1(L,facil=facil); E,V=eigh(H)
    tr=np.linspace(0,40,201); tw=np.linspace(0,2,101)
    tz0=op(L,{L:Z}); sy1=op(L,{1:Y})
    def ev(psi,ts):
        c=V.conj().T@psi; return V@(c[:,None]*np.exp(-1j*E[:,None]*ts[None,:]))
    m={}
    for r,zb in ((1,up),(-1,dn)):
        psi=prod([px]*L+[zb]*L); ps=ev(psi,tr); m[r]=np.einsum('nt,nm,mt->t',ps.conj(),tz0,ps).real
    succ=.5+np.abs(m[1]-m[-1])/4
    y={}
    for r,zb in ((1,up),(-1,dn)):
        mat=[up]*L; mat[1]=px; psi=prod(mat+[zb]*L); ps=ev(psi,tw); y[r]=np.einsum('nt,nm,mt->t',ps.conj(),sy1,ps).real
    rs=.5+np.abs(y[1]-y[-1])/4; k=int(np.argmax(rs))
    print(f"{label:28s} L={L}: min retention t<=10 = {succ[tr<=10].min():.6f}; at t=40 = {succ[-1]:.6f}; min over [0,40] = {succ.min():.6f}; best readout = {rs[k]:.6f} at t={tw[k]:.2f}")
for L in (3,4,5):
    run(L,True,"M1 candidate (facilitated)")
    run(L,False,"control: uniform field 0.05")
