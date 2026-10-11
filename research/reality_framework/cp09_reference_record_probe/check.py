"""CP-09 independent reference-record implementation and back-action control.

No replacement Hamiltonian. Original CP-09 J=1, Omega=1/4, gamma=pi/4.
Requires NumPy/SciPy. Runs on all posets for N=2,3.
"""
import itertools
import math
import numpy as np
from scipy.linalg import expm

J=1.0
OMEGA=0.25
GAMMA=math.pi/4
I2=np.eye(2, dtype=complex)
X=np.array([[0,1],[1,0]], dtype=complex)
Z=np.diag([1.,-1.]).astype(complex)

def at_site(n, idx, pauli):
    out=np.array([[1.]], dtype=complex)
    for k in range(n):
        out=np.kron(out, pauli if k==idx else I2)
    return out

def enumerate_posets(n):
    pairs=[(i,j) for i in range(n) for j in range(n) if i!=j]
    masks=[]
    for mask in range(1<<len(pairs)):
        C=np.zeros((n,n),dtype=int)
        for p,(i,j) in enumerate(pairs): C[i,j]=(mask>>p)&1
        if np.any(C&C.T):continue
        if np.all(((C@C)>0)<=C): masks.append(mask)
    return masks,pairs

def build(n):
    masks,pairs=enumerate_posets(n)
    idx={m:k for k,m in enumerate(masks)}
    d=1<<n
    Xs=[at_site(n,i,X) for i in range(n)]
    Zs=[at_site(n,i,Z) for i in range(n)]
    H=np.zeros((len(masks)*d,len(masks)*d),complex)
    onsite=OMEGA*(n*np.eye(d)-sum(Xs))
    for k,m in enumerate(masks):
        H[k*d:(k+1)*d,k*d:(k+1)*d]+=onsite
        for p,(i,j) in enumerate(pairs):
            add=m|(1<<p)
            if add==m or add not in idx:continue
            l=idx[add]
            U=np.cos(GAMMA)*np.eye(d)-1j*np.sin(GAMMA)*(Zs[i]@Xs[j])
            H[k*d:(k+1)*d,k*d:(k+1)*d]+=J*np.eye(d)
            H[l*d:(l+1)*d,l*d:(l+1)*d]+=J*np.eye(d)
            H[l*d:(l+1)*d,k*d:(k+1)*d]+=-J*U
            H[k*d:(k+1)*d,l*d:(l+1)*d]+=-J*U.conj().T
    assert np.linalg.norm(H-H.conj().T)<1e-11
    return H,masks,pairs,Xs,Zs

def role_operators(n,masks,pairs,Xs,Zs):
    d=1<<n; groups=[]; A=[];B=[]
    for mask in masks:
        C=np.zeros((n,n),dtype=bool)
        for p,(i,j) in enumerate(pairs):C[i,j]=bool((mask>>p)&1)
        lo=set(np.flatnonzero(~C.any(axis=0)).tolist())
        up=set(np.flatnonzero(~C.any(axis=1)).tolist())
        isolated=lo&up
        lo=tuple(sorted(lo-isolated));up=tuple(sorted(up-isolated))
        groups.append((lo,up))
        B.append(sum((Xs[i] for i in lo),np.zeros((d,d),complex))/max(1,len(lo)))
        A.append(sum((Zs[i] for i in up),np.zeros((d,d),complex))/max(1,len(up)))
    tags={role:i for i,role in enumerate(sorted(set(groups)))}
    tag_index=np.array([tags[g] for g in groups],dtype=int)
    refdim=1<<(len(tags)-1).bit_length()
    Bmat=np.zeros((len(masks)*d,len(masks)*d),complex)
    Amat=np.zeros_like(Bmat)
    for k in range(len(masks)):
        Bmat[k*d:(k+1)*d,k*d:(k+1)*d]=B[k]
        Amat[k*d:(k+1)*d,k*d:(k+1)*d]=A[k]
    assert np.linalg.norm(Bmat@Amat-Amat@Bmat)<1e-12
    return A,B,Amat,Bmat,tags,tag_index,refdim

def simulate(n):
    H,masks,pairs,Xs,Zs=build(n)
    A,B,Amat,Bmat,tags,tagidx,rd=role_operators(n,masks,pairs,Xs,Zs)
    ener,vec=np.linalg.eigh(H)
    g=vec[:,0];E0=ener[0]
    G=len(masks);m=1<<n
    groups=[]
    for tag in range(len(tags)):
        v=np.zeros_like(g)
        indices=np.where(tagidx==tag)[0]
        for i in indices:v[i*m:(i+1)*m]=g[i*m:(i+1)*m]
        groups.append(v)
    # Exact reversible oracle V: |C,m,0> --> |C,m,f(C)>.
    eps=0.123
    copied=np.zeros((G,m,rd),complex)
    gg=g.reshape(G,m)
    for k in range(G):copied[k,:,tagidx[k]]=gg[k]
    # Controlled pulse, using record as coherent control; uncompute XOR.
    kicked=np.zeros_like(copied)
    for k in range(G):
        gate=expm(1j*eps*B[k])
        kicked[k,:,tagidx[k]]=gate@copied[k,:,tagidx[k]]
    restored=np.zeros_like(copied)
    for k in range(G):restored[k,:,0]=kicked[k,:,tagidx[k]]
    direct=np.concatenate([expm(1j*eps*B[k])@gg[k] for k in range(G)])
    oracle_err=np.linalg.norm(restored[:,:,0].reshape(-1)-direct)
    assert oracle_err<1e-12
    assert np.linalg.norm(restored[:,:,1:])<1e-12

    # Persistent record -> dephasing between distinct *role* classes.
    mixture_energy=sum(np.vdot(v,H@v).real for v in groups)
    deltaE=mixture_energy-E0
    role_prob=np.array([np.vdot(v,v).real for v in groups])
    record_fidelity=float(np.dot(role_prob,role_prob))
    gap=float(ener[1]-ener[0])
    lower_bound=gap*(1-record_fidelity)
    assert deltaE >= lower_bound-1e-11
    assert deltaE>1e-5
    time=.4; Ut=expm(-1j*time*H)
    base=Ut@g
    orig_geom=np.sum(abs(g.reshape(G,m))**2,axis=1)
    recorded_geom=sum((np.sum(abs((Ut@v).reshape(G,m))**2,axis=1) for v in groups),np.zeros(G))
    tv=.5*np.sum(abs(recorded_geom-orig_geom))
    origR=-2*np.vdot(base,Amat@(Ut@(Bmat@g))).imag
    recR=sum((-2*np.vdot(Ut@v,Amat@(Ut@(Bmat@v))).imag for v in groups))
    impulse=1e-5
    pulseplus=expm(1j*impulse*Bmat)
    pulseminus=expm(-1j*impulse*Bmat)
    def expected_A(v):
        v=Ut@v
        return np.vdot(v,Amat@v).real
    impulse_original=(expected_A(pulseplus@g)-expected_A(pulseminus@g))/(2*impulse)
    impulse_recorded=sum((expected_A(pulseplus@v)-expected_A(pulseminus@v))/(2*impulse) for v in groups)
    assert abs(impulse_original-origR)<1e-8
    assert abs(impulse_recorded-recR)<1e-8
    assert tv>1e-6
    assert abs(origR-recR)>1e-6
    assert abs(np.linalg.norm(base)-1)<1e-12
    print(f'N={n}: posets={G}, different role records={len(tags)}, ancilla bits={int(math.log2(rd))}')
    print(f' E0={E0:.10f}, gap={gap:.10f}, record_fidelity={record_fidelity:.10f}, gap_lower_bound={lower_bound:.10f}')
    print(f' record_energy_increase={deltaE:.10f}, geometry_TV_at_t.4={tv:.10f}')
    print(f' unrecorded_relational_response={origR:.10f}, persistent_record_response={recR:.10f}, oracle_error={oracle_err:.3g}, recorded_impulse_error={abs(impulse_recorded-recR):.3g}')
    assert abs(origR - {2:0.05258851,3:0.10919014}[n])<5e-7
    return dict(N=n, posets=G, roles=len(tags), bits=int(math.log2(rd)),E0=E0, record_energy=deltaE, recorded_geom_tv=tv, response_no_record=origR, response_record=recR)

if __name__=='__main__':
    for n in [2,3]:simulate(n)
    print('CP09_REFERENCE_RECORD_AND_BACKACTION_PASS')
