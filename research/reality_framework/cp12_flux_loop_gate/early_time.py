"""Resolve the onset channels at early times in the unchanged N4 state."""
import json,time
import numpy as np
from loops import P,model,resources,observables
from onset_uniform import formulas
from evolution import Evolution


def run():
    start=time.monotonic();n=q=4;budget=256
    sub,H,B,states,sizes,psi=model(n,q,budget)
    names,O,_=observables(n,q,states);r,k,nv,e2=resources(n,q,sub,states)
    detectors=dict(protected_triangle=O[:,names.index('protected_cycle_3_any')],
        occupied_flux_triangle=O[:,names.index('occupied_flux_cycle_3_any')],
        occupied_flux_fraction=k/q,protected_edge_fraction=(q-r)/q)
    out=[]
    for t in [.02,.04,.08,.16,.32]:
        E=Evolution(H,B,2.,max(8,round(t/.005)))
        v,_,e=E(psi,t);p=abs(v)**2
        vals={key:float(p@value) for key,value in detectors.items()}
        intervals={key:[max(0.,np.sqrt(value)-e)**2,min(1.,(np.sqrt(value)+e)**2)] for key,value in vals.items()}
        out.append(dict(T=t,state_error_upper=e,values=vals,positive_score_intervals=intervals))
    f=formulas(n,q)
    result=dict(N=n,Q=q,budget=budget,rows=out,exact_leading_coefficients={key:str(v) for key,v in f.items()},
        protected_edge_fraction_T2_coefficient=2,seconds=time.monotonic()-start,
        scope='Additional earlier sample times in the same supplied ray and unchanged Hamiltonian, no new size or longer evolution. Full-law outgoing bounds plus positive-contraction amplitude intervals. Leading orders are fixed-system statements, not macroscopic time scales. Floating-point roundoff is not formally interval-certified.')
    (P/'EARLY_TIME.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True);return result


if __name__=='__main__':run()
