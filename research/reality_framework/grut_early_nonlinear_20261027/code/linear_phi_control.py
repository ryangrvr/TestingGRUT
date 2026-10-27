"""Same finite-band phi realization driven by the unperturbed reduced orbit."""
from en_common import *

def main():
    tag='N48';params=json.loads((WORK/tag/'parameters.json').read_text());N=params['N'];L=params['L']
    raw=np.memmap(WORK/tag/'initial.bin',offset=8+4+8+32,dtype='<f8',shape=(6,N,N,N))
    q=np.array(raw[:3]);p=np.array(raw[3:]);kh=symbol(N,L);qh=rfftn(q[2],norm='forward');ph=rfftn(p[2],norm='forward')
    weights=np.full(N//2+1,2.);weights[[0,-1]]=1;weights=np.broadcast_to(weights[None,None,:],kh.shape)
    use=(abs(qh)+abs(ph))>1e-12
    values,ids=np.unique(np.round(kh[use],12),return_inverse=True);nk=len(values)
    qq=np.bincount(ids,weights=(weights*abs(qh)**2)[use],minlength=nk)
    pp=np.bincount(ids,weights=(weights*abs(ph)**2)[use],minlength=nk)
    qp=np.bincount(ids,weights=(weights*np.real(qh.conj()*ph))[use],minlength=nk)
    trace=np.genfromtxt(WORK/tag/'state.tsv',names=True,delimiter='\t');elapsed=trace['eta'][::10]
    initial=np.r_[q[:2].mean(axis=(1,2,3)),p[:2].mean(axis=(1,2,3)),np.tile(np.eye(2).ravel(),nk)]
    def rhs(s,y):
        x=y[:2];v=y[2:4];r=x@x;Tphi=y[4:].reshape(nk,2,2)
        out=np.empty_like(y);out[:4]=np.r_[v,-LX*r*x]
        out[4:]=np.stack([Tphi[:,1,:],-(values+G2*r)[:,None]*Tphi[:,0,:]],axis=1).ravel()
        return out
    sol=solve_ivp(rhs,(0,elapsed[-1]),initial,t_eval=elapsed,method='DOP853',rtol=2e-11,atol=2e-13);assert sol.success
    rows=[];deterr=0.
    for j,s in enumerate(elapsed):
        y=sol.y[:,j];r=y[:2]@y[:2];t=y[4:].reshape(nk,2,2)
        cqq=t[:,0,0]**2*qq+2*t[:,0,0]*t[:,0,1]*qp+t[:,0,1]**2*pp
        cpp=t[:,1,0]**2*qq+2*t[:,1,0]*t[:,1,1]*qp+t[:,1,1]**2*pp
        e=float(.5*np.sum(cpp+(values+G2*r)*cqq))
        deterr=max(deterr,float(abs(np.linalg.det(t)-1).max()))
        A=params['clock'];tau=A['tau0']+A['a0']*s+.5*A['C']*s*s
        rows.append(dict(tau=tau,phi_plus_interaction_raw_fraction=e/trace['rho'][0]))
    values=np.array([r['phi_plus_interaction_raw_fraction'] for r in rows]);first=np.flatnonzero(values>=1)
    save('LINEAR_PHI_CONTROL.json',dict(rows=rows,first_sampled_unit_fraction=float(rows[first[0]]['tau']) if len(first) else None,
        maximum_fraction=float(values.max()),transfer_determinant_error=deterr,
        scope='Same N48 finite-band phi seed, exact discrete momenta, homogeneous reduced X orbit and linear phi. Includes interaction energy; no vacuum subtraction. An extrapolation/control, not a nonlinear population or the full review neutral calculation.'))
    print('PHI CONTROL',nk,values[-1],deterr,flush=True)

if __name__=='__main__':main()
