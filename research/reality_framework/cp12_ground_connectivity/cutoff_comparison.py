"""Numerical Schur-comparison checks in EVERY S_N sector.

This extends the analytic unbounded-rotor test in cutoff_certificate.py to
larger multiplicity blocks. Eigenvalues here use floating Lanczos, so the
reported positive minimum is numerical evidence, not rational certification.
"""
import argparse
import json
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from representation_ground import get_orbits,Reduced,partitions,P
from cutoff_certificate import boundary_data,tail_lower

def run(n,L,z,resume=True):
    sub=get_orbits(n,L);count,_=boundary_data(n,L,sub.states);count=np.array(count)
    c=float(tail_lower(n,L));assert c>z
    output=P/f'cutoff_comparison_N{n}_L{L}.json'
    rows=[]
    if resume and output.exists():
        old=json.loads(output.read_text())
        if old['comparison_threshold']==z:rows=old['sectors']
    previous={tuple(r['shape']) for r in rows}
    for shape in partitions(n):
        if shape in previous:continue
        start=time.monotonic();red=Reduced(sub,shape)
        if not red.multiplicity:continue
        diagonal=np.empty(red.multiplicity)
        for ids,V,indices in red.patterns:diagonal[indices]=count[ids,None]
        A=LinearOperator(red.H.shape,matvec=lambda v:red.H@v-z*v-diagonal*v/(c-z),dtype=float)
        initial=np.random.default_rng(9107).standard_normal(red.multiplicity)
        E,v=eigsh(A,k=1,which='SA',v0=initial,tol=2e-9,maxiter=20000,ncv=30)
        res=float(np.linalg.norm(A@v[:,0]-E[0]*v[:,0]));assert res<2e-8
        row=dict(shape=shape,multiplicity=red.multiplicity,comparison_minimum=float(E[0]),
                 residual=res,seconds=time.monotonic()-start)
        rows.append(row)
        result=dict(status='NUMERICAL_ALL_SECTOR_COMPARISON_PASS'
                    if len(rows)==len(list(partitions(n))) and min(r['comparison_minimum'] for r in rows)>0
                    else 'PARTIAL_OR_NOT_POSITIVE',N=n,Q=n,L=L,comparison_threshold=z,
                    tail_lower=c,tail_lower_rational=str(tail_lower(n,L)),
                    full_dimension=sub.full_dimension,sectors=rows,
                    exact_rational_certificate=False,
                    scope='All-sector floating smallest-eigenvalue tests; NOT interval-certified positivity.')
        output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(row),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True)
    p.add_argument('--L',type=int,required=True);p.add_argument('--z',type=float,required=True)
    p.add_argument('--fresh',action='store_true',help='Recompute all sectors.')
    a=p.parse_args();run(a.N,a.L,a.z,not a.fresh)
