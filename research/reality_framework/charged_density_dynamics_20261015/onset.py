"""Exact first-order charged-density coefficient in weighted orbit integers."""
from __future__ import annotations
from fractions import Fraction
from math import comb,factorial
from pathlib import Path
import argparse,json
import numpy as np
from scipy.sparse import coo_matrix
from formation import model
from diagnostics import arrays,event

P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--budget',type=int,required=True)
    args=ap.parse_args();n=args.n;q=2*n
    _,H,B,states,sizes,psi=model(n,q,args.budget)
    diag=arrays(n,q,states);mask=event(diag,q,n)
    coo=H.tocoo();keep=diag[1,coo.row]!=diag[1,coo.col]
    row,col=coo.row[keep],coo.col[keep]
    raw=coo.data[keep]*np.sqrt(sizes[row]/sizes[col]);integer=np.rint(raw).astype(np.int64)
    rounding=float(np.max(abs(raw-integer)))
    assert rounding<1e-10
    K=coo_matrix((integer,(row,col)),shape=H.shape).tocsr()
    del coo,row,col,raw
    M=n*(n-1)//2;sea=sum(1<<(n+e) for e in range(M))
    y=np.asarray([int(sizes[i]) if b==sea and not any(E) else 0 for i,(g,b,E) in enumerate(states)],dtype=np.int64)
    rows=[];d=2*int(np.ceil(q/4))+int(np.ceil(n/4))
    assert d==8 and args.budget>=12
    bound=max(np.asarray(abs(K).sum(axis=1)).reshape(-1))
    for k in range(d+1):
        assert np.all(y%sizes==0)
        physical=y//sizes
        norm=sum(int(sizes[i])*int(physical[i])**2 for i in np.flatnonzero(mask& (y!=0)))
        rows.append(dict(order=k,target_squared_integer_norm=str(norm)))
        if k<d:
            assert norm==0
            assert int(np.max(abs(y)))*int(bound)<2**63
            y=K@y
    coefficient=Fraction(norm,comb(M,q)*factorial(d)**2)
    assert coefficient>0
    result=dict(status='PASS',N=n,Q=q,c_min=int(np.ceil(q/4)),nv_min=int(np.ceil(n/4)),
        first_amplitude_order=d,leading_probability_order=2*d,
        leading_probability_coefficient=str(coefficient),coefficient_float=float(coefficient),
        integer_squared_norm=str(norm),neutral_graphs=comb(M,q),rows=rows,
        exact_integer_action_after_matrix_entry_recovery=True,entry_rounding_residual=rounding,
        scope='Minimum electric L1 is 2c+nv=8. Every contributing path is monotone in L1; all target fields then have electric square <=12, so budget>=12 contains every leading path. Integer weighted orbit action is normalized by the original neutral graph count. No long-time extrapolation.')
    (P/f'ONSET_N{n}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
