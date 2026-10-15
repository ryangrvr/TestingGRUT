"""Neutral full-H Galerkin/Krylov evolution with explicit residual bounds.

Real orthonormal vectors save memory. Projection is purely numerical:
both the last Arnoldi residual and every omitted physical row B are kept.
"""
from math import sqrt
import argparse,json,time
import numpy as np
from scipy.sparse import csr_matrix
from shared import *
from neutral_prepare import save_results


def reduced(m,k):
    H,B=m['H'],m['B'];psi=m['psi'];dimension=len(psi)
    assert np.linalg.norm(psi.imag)<1e-14
    V=np.empty((dimension,k),float,order='F');V[:,0]=psi.real
    start=time.monotonic();last=start
    for j in range(k-1):
        u=H@V[:,j]
        for _ in range(2):u-=V[:,:j+1]@(V[:,:j+1].T@u)
        beta=float(np.linalg.norm(u));assert beta>1e-12
        V[:,j+1]=u/beta
        if time.monotonic()-last>15:
            print(json.dumps(dict(stage='krylov_basis',columns=j+2,total=k,seconds=time.monotonic()-start)),flush=True)
            last=time.monotonic()
    orth=float(np.linalg.norm(V.T@V-np.eye(k)));assert orth<1e-10
    T=np.empty((k,k))
    for j in range(k):T[:,j]=V.T@(H@V[:,j])
    asym=float(np.linalg.norm(T-T.T));assert asym<1e-10;T=(T+T.T)/2
    # A last-column residual plus an explicitly evaluated full numerical
    # relation defect covers reorthogonalization and symmetrization effects.
    tail=H@V[:,-1]-V@T[:,-1];beta=float(np.linalg.norm(tail))
    defect2=0.
    for j in range(k-1):
        residual=H@V[:,j]-V@T[:,j];defect2+=float(np.vdot(residual,residual).real)
    defect=sqrt(defect2)
    BV=B@V
    # Thin QR preserves ||B V c|| while reducing a tall boundary to k rows.
    R=np.linalg.qr(BV,mode='r')
    rng=np.random.default_rng(128);qr_error=0.
    for _ in range(4):
        c=rng.normal(size=k)+1j*rng.normal(size=k);c/=np.linalg.norm(c)
        reference=np.linalg.norm(BV@c);actual=np.linalg.norm(R@c)
        qr_error=max(qr_error,abs(reference-actual)/max(1.,reference))
    assert qr_error<1e-11;del BV
    last_action=np.zeros((1,k));last_action[0,-1]=beta
    info=dict(krylov_dimension=k,orthogonality_error=orth,projected_H_asymmetry=asym,
        full_finite_H_relation_defect_frobenius=defect,last_column_residual_norm=beta,
        boundary_QR_relative_norm_check=qr_error,seconds=time.monotonic()-start)
    print(json.dumps(dict(stage='krylov_ready',**info)),flush=True)
    return V,T,{'electric_boundary':csr_matrix(R),'last_column':csr_matrix(last_action)},info


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True)
    ap.add_argument('--q',type=int,required=True);ap.add_argument('--budget',type=int,required=True)
    ap.add_argument('--columns',type=int,default=128);ap.add_argument('--waits',default='5')
    a=ap.parse_args();start=time.monotonic();m=quotient_model(a.n,a.q,a.budget)
    V,T,actions,info=reduced(m,a.columns);initial=np.zeros(a.columns,complex);initial[0]=1
    times=sorted(float(v) for v in a.waits.split(','))
    result=Propagator(csr_matrix(T),actions,a.q/2,dt=1/64,batch=16).run(initial,times)
    physical={};bounds=[]
    for t,(c,er) in result.items():
        # The two boundary channels are added by the triangle inequality;
        # this is conservative and does not assume their orthogonality.
        error=float(er['electric_boundary'][0]+er['last_column'][0]+t*info['full_finite_H_relation_defect_frobenius']+ALLOWANCE)
        physical[t]=(V@c,{'B':np.array([error])})
        bounds.append(dict(wait=t,electric_boundary=float(er['electric_boundary'][0]),
            last_column=float(er['last_column'][0]),relation_defect=t*info['full_finite_H_relation_defect_frobenius'],
            total=error))
    meta=dict(backend='real full-H Galerkin/Krylov, two-pass reorthogonalization',
        krylov=info,krylov_bounds=bounds,dt=1/64,
        bound='Duhamel bound from complete electric boundary B V, last Arnoldi column, and full finite-H relation defect; thin QR preserves the boundary action norm. Projection changes no physical Hamiltonian term.')
    save_results(m,a.n,a.q,a.budget,physical,meta,start)


if __name__=='__main__':main()
