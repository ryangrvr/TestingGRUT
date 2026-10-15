"""Compare the neutral Krylov reduction with independent full-basis evolution."""
import json,time
import numpy as np
from scipy.sparse import csr_matrix
from shared import *
from neutral_krylov import reduced


def main():
    start=time.monotonic();rows=[]
    for n,q,b,k in [(4,5,256,64),(5,6,256,96),(6,7,192,144)]:
        m=quotient_model(n,q,b);V,T,actions,info=reduced(m,k)
        c=np.zeros(k,complex);c[0]=1
        reduced_result=Propagator(csr_matrix(T),actions,q/2,dt=1/64,batch=16).run(c,[5.])
        zz=np.load(cache_root()/f'N{n}_Q{q}_B{b}_W5.npz')
        coef,er=reduced_result[5.];w=V@coef
        bound=sum(float(v[0]) for v in er.values())+5*info['full_finite_H_relation_defect_frobenius']+ALLOWANCE
        difference=float(np.linalg.norm(w-zz['w']))
        # Both approximations have bounds against the same exact full state.
        # N6 compares an enlarged Krylov space to the original 128-column
        # computation, while N4/5 use independent full-basis exponentials.
        assert difference<=bound+float(zz['state_error'])+ALLOWANCE
        rows.append(dict(N=n,Q=q,budget=b,columns=k,full_basis_difference=difference,
            krylov_full_state_bound=bound,independent_full_basis_bound=float(zz['state_error']),**info))
    data=dict(status='PASS',rows=rows,seconds=time.monotonic()-start)
    (P/'KRYLOV_CHECK.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data),flush=True)


if __name__=='__main__':main()
