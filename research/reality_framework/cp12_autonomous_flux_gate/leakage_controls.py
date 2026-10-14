"""Fresh parent residual controls without rewriting the parent packet."""
from pathlib import Path
import importlib.util,json,sys
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix

P=Path(__file__).resolve().parent
parent=P.parent/'cp12_distance_history'
sys.path.insert(0,str(parent))
spec=importlib.util.spec_from_file_location('cp12_previous_leakage_checks',parent/'check.py')
checks=importlib.util.module_from_spec(spec);spec.loader.exec_module(checks)
from evolution import Evolution


def expectation_dense():
    out=[]
    for seed in range(8):
        rng=np.random.default_rng(801+seed);dim=8;p=4
        a=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim));H=(a+a.conj().T)/3
        steps=2 if seed%2 else 8
        HP=csr_matrix(H[:p,:p]);B=csr_matrix(H[p:,:p]);E=Evolution(HP,B,0.,steps)
        v=rng.normal(size=p)+1j*rng.normal(size=p);v/=np.linalg.norm(v)
        O=rng.uniform(size=dim);z,_,e=E(v,.2)
        _,_,d=E(O[:p]*z,.2,back=True)
        ref=expm(-.2j*H)@np.r_[v,np.zeros(dim-p)]
        actual=abs(float(abs(ref)**2@O-abs(z)**2@O[:p]))
        bound=2*e*d+e*e
        assert actual<=bound+1e-12
        out.append(dict(seed=seed,leakage_subintervals=steps,actual_error=actual,detector_bound=bound))
    return out


def main():
    out=dict(status='CP12_AUTONOMOUS_FLUX_LEAKAGE_CONTROLS_PASS',
             dense=checks.dense(),expectation_dense=expectation_dense(),ordered_physical=checks.physical())
    (P/'LEAKAGE_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],dense_tests=len(out['dense']),
                         physical_dimension=out['ordered_physical']['full_dimension'])),flush=True)


if __name__=='__main__':main()
