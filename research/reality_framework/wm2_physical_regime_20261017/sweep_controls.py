"""Production residual checks plus complete independent physical references."""
from pathlib import Path
import json
import time
import numpy as np
import solver
from certificate import residual_enclosure
import bounds as E
import controls as S8_CONTROLS
import sweep_control as S8_LEMMA


def direct_metrics():
    rng=np.random.default_rng(918172);rows=[]
    for left,physical,right,rank in [(3,2,4,2),(4,3,5,3),(2,2,3,1)]:
        for perturb in [0.,2e-9]:
            L=np.linalg.qr(rng.normal(size=(left+2,left))+1j*rng.normal(size=(left+2,left)))[0]
            L=L@np.diag(np.sqrt(1+np.linspace(-2e-8,2e-8,left)))
            T=np.kron(L,np.eye(physical));H=T.conj().T@T
            K=np.linalg.qr(rng.normal(size=(left*physical,rank))+1j*rng.normal(size=(left*physical,rank)))[0]
            K=K@np.diag(1+np.linspace(0,perturb,rank))
            X=(rng.normal(size=(left*physical,right))+1j*rng.normal(size=(left*physical,right)))/10
            Z=K.conj().T@X+1e-8*(rng.normal(size=(rank,right))+1j*rng.normal(size=(rank,right)))
            rightmap=np.linalg.qr(rng.normal(size=(right+2,right))+1j*rng.normal(size=(right+2,right)))[0]
            mapped=T@K;P=mapped@np.linalg.solve(mapped.conj().T@mapped,mapped.conj().T)
            y=T@X@rightmap.conj().T;x=mapped@Z@rightmap.conj().T
            loss=float(np.linalg.norm((np.eye(len(P))-P)@y))
            defect=float(np.linalg.norm(P@y-x))
            cert=residual_enclosure(X,K,Z,E.gram(K),2e-7,1/(1-2e-7))
            assert loss<=cert['orthogonal_loss_upper']+1e-14
            assert defect<=cert['in_range_defect_upper']+1e-14
            rows.append(dict(shape=list(X.shape),rank=rank,nonorthogonality=perturb,
                actual_loss=loss,actual_in_range_defect=defect,certificate=cert,
                observed_prefix_Gram=float(np.linalg.norm(H-np.eye(len(H)),2))))
    return rows


def structural_guards():
    e=solver.Receiver([1,1],[-.25,.25],5);e.initialize(5)
    e.begin_batch('right',4)
    try:e.begin_batch('left',4)
    except RuntimeError:overlap_rejected=True
    else:overlap_rejected=False
    e.end_batch()
    assert overlap_rejected
    return dict(overlapping_batch_rejected=True,production_checks='direction and consecutive bond indices on every eligible gate')


def run():
    start=time.monotonic()
    result=dict(direct_nonorthogonal_prefix=direct_metrics(),structural=structural_guards(),
        nested_projectors=S8_LEMMA.run())
    S8_CONTROLS.Receiver=solver.Receiver
    result['independent_physical_and_arithmetic']=S8_CONTROLS.run()
    result.update(status='PASS',elapsed_seconds=time.monotonic()-start,
        arithmetic='production enclosure under declared IEEE/BLAS model; dense controls ordinary double')
    return result


if __name__=='__main__':
    r=run();Path(__file__).with_name('SWEEP_CONTROLS.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
