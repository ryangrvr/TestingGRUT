"""Leading-order formulas, checked against all existing main matrices."""
from fractions import Fraction as F
from math import comb,factorial
import json,time
import numpy as np
from loops import P,model,resources,observables


def formulas(n,q):
    M=comb(n,2)
    triangles=F(comb(n,3)*comb(q,3),comb(M,3)) if q>=3 else F(0)
    return dict(protected_triangle_T6=2*triangles,
                occupied_flux_triangle_T12=F(2,225)*triangles,
                occupied_flux_fraction_T6=F(2*(n-2)*(q-1),9*(M-1)))


def run():
    start=time.monotonic();rows=[]
    for n,q,budget in [(4,4,256),(5,5,128),(6,6,128),(6,12,10)]:
        sub,H,B,states,sizes,psi=model(n,q,budget)
        names,O,_=observables(n,q,states);r,k,nv,e2=resources(n,q,sub,states)
        expected=formulas(n,q);vectors={};v=psi.copy()
        for j in range(1,7):
            v=H@v
            if j==3:vectors[3]=v.copy()
        vectors[6]=v
        p3=abs(vectors[3])**2/factorial(3)**2
        p6=abs(vectors[6])**2/factorial(6)**2
        actual=dict(protected_triangle_T6=float(p3@O[:,names.index('protected_cycle_3_any')]),
            occupied_flux_triangle_T12=float(p6@O[:,names.index('occupied_flux_cycle_3_any')]),
            occupied_flux_fraction_T6=float(p3@k/q))
        error=max(abs(actual[key]-float(value)) for key,value in expected.items())
        assert error<1e-9
        rows.append(dict(N=n,Q=q,budget=budget,exact_formulas={key:str(v) for key,v in expected.items()},
            numerical_matrix_coefficients=actual,maximum_formula_difference=error))
        print(json.dumps(dict(stage='uniform_onset',N=n,Q=q,maximum_error=error)),flush=True)
    out=dict(rows=rows,seconds=time.monotonic()-start,
        scope='Exact leading formulas derived from minimal hopping paths; matrix checks use existing projected Hamiltonians. Minimal triangle-flux paths have six distinct unit-field changes, so their electric-square support never exceeds six and all retained budgets include them. No long-time Taylor extrapolation.')
    (P/'UNIFORM_ONSET.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':run()
