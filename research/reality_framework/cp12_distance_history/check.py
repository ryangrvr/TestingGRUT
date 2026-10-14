"""Independent full-basis, recording, and high-order quadrature controls."""
from pathlib import Path
from math import log2
import json,numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix,coo_matrix,eye
from scipy.sparse.linalg import expm_multiply
from basis import model,distance
from evolution import Evolution,dual
from frozen_model import build

P=Path(__file__).resolve().parent


def dense():
    records=[]
    for seed in range(16):
        rng=np.random.default_rng(489+seed);dim=9;p=4
        a=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim));H=(a+a.conj().T)/3
        HP=csr_matrix(H[:p,:p]);B=csr_matrix(H[p:,:p]);E=Evolution(HP,B,0.,8)
        V=rng.normal(size=(p,3))+1j*rng.normal(size=(p,3));V/=np.linalg.norm(V)
        out,ec,ef=E(V,.12,groups=[range(2),[2]])
        full=expm(-.12j*H)@np.vstack([V,np.zeros((dim-p,3))])
        actual=np.linalg.norm(full-np.vstack([out,np.zeros((dim-p,3))]),axis=0)
        assert np.all(actual<=ec+1e-12)
        assert np.linalg.norm(actual[:2])<=ef[0]+1e-12
        ref=expm(-.12j*HP.toarray())@V
        assert np.linalg.norm(out-ref)<1e-11
        records.append(dict(seed=seed,max_state_error=float(max(actual)),maximum_column_bound=float(max(ec)),group_bounds=ef.tolist()))
    return records


def physical():
    n=4;q=3;budget=3;sub,H,B,states,sizes,psi=model(n,q,budget)
    full,fs,_,_=build(n,q,2);index={s:i for i,s in enumerate(states)}
    for g,bits,field in fs:
        removable=sum(bool(g>>e&1) and bool(bits>>(n+e)&1) and field[2*e]==field[2*e+1]==0 for e in range(sub.M))
        assert sum(x*x for x in field)>=q-removable
    rr=[];cc=[];vv=[]
    for i,s in enumerate(fs):
        if sum(x*x for x in s[2])>budget:continue
        key,sign,size=sub.canonical(s)
        if sign and key in index:rr.append(i);cc.append(index[key]);vv.append(sign/np.sqrt(size))
    Tmap=coo_matrix((vv,(rr,cc)),shape=(len(fs),len(states))).tocsr()
    iso=float(np.linalg.norm((Tmap.T@Tmap-eye(H.shape[0])).data));assert iso<1e-10
    ph=float(np.linalg.norm((Tmap.T@full@Tmap-H).data));assert ph<1e-10
    L=full@Tmap-Tmap@H
    gram=float(np.linalg.norm((L.T@L-B.T@B).data));assert gram<1e-9
    E=Evolution(H,B,q/2,16);D=distance(n,states,sub.es);DF=distance(n,fs,sub.es)
    S=np.array([b&1 for g,b,e in states]);Sf=np.array([b&1 for g,b,e in fs])
    O=np.array([(b>>1)&1 for g,b,e in states]);Of=np.array([(b>>1)&1 for g,b,e in fs])
    v,_,ew=E(psi,.1);vf=expm_multiply(-.1j*full,Tmap@psi)
    cs=sorted(set(D.tolist()));F=np.column_stack([D==r for r in cs]);FF=np.column_stack([DF==r for r in cs])
    W=F*v[:,None];Wf=FF*vf[:,None]
    X,ex,gx=E(W,.1);Y,ey,gy=E(S[:,None]*W,.1)
    XF=expm_multiply(-.1j*full,Wf);YF=expm_multiply(-.1j*full,Sf[:,None]*Wf)
    group_record_gram=W.conj().T@W
    assert np.max(abs(group_record_gram-np.diag(np.diag(group_record_gram))))<1e-13
    rows=[]
    for j,r in enumerate(cs):
        f=F[:,j]
        bound=dual(E,X[:,j],Y[:,j],O,lambda z,f=f:f*z,lambda z,f=f:f*S*z,.1,.1,float(ew+ex[j]),float(ew+ey[j]))
        actual=float(2*np.vdot(XF[:,j],Of*YF[:,j]).imag)
        assert abs(actual-bound['value'])<=bound['error']+1e-12
        post=dual(E,X,Y,O*f,lambda z:np.sum(F*z,axis=1),lambda z:S*np.sum(F*z,axis=1),.1,.1,float(ew+gx),float(ew+gy))
        af=float(2*np.vdot(XF,Of[:,None]*FF[:,j,None]*YF).imag)
        assert abs(af-post['value'])<=post['error']+1e-12
        rows.append(dict(distance=r,past_actual=actual,past_projected=bound['value'],past_error_upper=bound['error'],post_actual=af,post_projected=post['value'],post_error_upper=post['error']))
    # Exact classical-history sum rules and coherence/backaction identity.
    C=np.array([[2*np.vdot(X[:,i],O*F[:,j]*Y[:,i]).imag for j in range(len(cs))] for i in range(len(cs))])
    recorded=2*np.vdot(X,O[:,None]*Y).imag
    coherent=2*np.vdot(X.sum(axis=1),O*Y.sum(axis=1)).imag
    assert abs(C.sum()-recorded)<1e-12
    cross=2*sum(np.vdot(X[:,i],O*Y[:,j]) for i in range(len(cs)) for j in range(len(cs)) if i!=j).imag
    assert abs(coherent-recorded-cross)<1e-12
    # Weighted cross-sector inequality, including charged/rotor configurations.
    absent=np.flatnonzero(D!=1);present=np.flatnonzero(D==1)
    Cedge=H[present,:][:,absent]
    movable=[]
    for k in absent:
        g,b,field=states[k]
        movable.append(sum(bool(g>>e&1) and bool(b>>(n+e)&1) and field[2*e]==field[2*e+1]==0 for e in range(sub.M)))
    weighted=np.diag(np.array(movable)*(sub.M-q)/(100*sub.M**2))-(Cedge.T@Cedge).toarray()
    minimum=float(np.linalg.eigvalsh(weighted)[0]);assert minimum>=-1e-12
    return dict(N=n,Q=q,budget=budget,full_dimension=len(fs),ordered_dimension=len(states),isometry_error=iso,projected_H_error=ph,leakage_Gram_error=gram,rows=rows,recording_orthogonality_and_sum_rules=True,weighted_pinning_inequality_min_eigenvalue=minimum,electric_energy_pinning_inequality_exact_basis_check=True)

if __name__=='__main__':
    out=dict(status='CP12_DISTANCE_HISTORY_CONTROLS_PASS',dense=dense(),physical=physical())
    (P/'CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
