"""Post-campaign independent small-matrix test; not used by S7 workers."""
import json
from pathlib import Path
from math import sqrt
import numpy as np
from scipy.linalg import eigh
from layout import Layout
from controls import REF


def run():
    rng=np.random.default_rng(73491);rows=[]
    for counts in [[1,1],[2,1]]:
        L,rank,g=5,2,.9;M=sum(counts);N=M+1
        full=Layout(counts,[-.25,.25],L,g);small=Layout(counts,[-.25,.25],rank,g)
        independent=REF.A.Layout(counts,[-.25,.25],L)
        w=np.linalg.qr(rng.normal(size=(L,rank))+1j*rng.normal(size=(L,rank)))[0]
        c=rng.normal(size=L)+1j*rng.normal(size=L);c/=np.linalg.norm(c)
        z=w.conj().T@c;sigma=float(np.linalg.norm(c-w@z))
        v=np.zeros((full.size,small.size),complex)
        field=[REF.field_isometry(w,k) for k in range(M+1)]
        for wi,u in enumerate(small.words):
            k=len(u);v[full.offset[wi]:full.offset[wi]+full.nocc[k],
                small.offset[wi]:small.offset[wi]+small.nocc[k]]=field[k]
        eye=np.eye(full.size,dtype=complex)
        b=np.column_stack([independent.hint(np.ascontiguousarray(eye[:,j]),c,g) for j in range(full.size)])
        compressed=v.conj().T@b@v;d=b@v-v@compressed
        weight=np.sqrt(M-np.repeat(small.lengths,small.nocc[small.lengths]))
        rhs=np.diag((g*sigma/N)**2*weight**2)
        identity_error=float(np.linalg.norm(d.conj().T@d-rhs,2))
        comm=(weight[:,None]-weight[None,:])*compressed
        comm_norm=float(np.linalg.norm(comm,2));G=g*float(np.linalg.norm(z))
        assert identity_error<2e-14 and comm_norm<=G+2e-14
        eb,ub=eigh(b);ec,uc=eigh(compressed)
        packets=[]
        for kind in ['random','fully_assembled']:
            x=rng.normal(size=small.size)+1j*rng.normal(size=small.size)
            if kind=='fully_assembled':x[:small.koff[M]]=0
            x/=np.linalg.norm(x)
            for h in [.02,.1]:
                exact=ub@(np.exp(-1j*h*eb)*(ub.conj().T@(v@x)))
                approx=v@(uc@(np.exp(-1j*h*ec)*(uc.conj().T@x)))
                error=float(np.linalg.norm(exact-approx))
                bound=(g*sigma/N)*(h*np.linalg.norm(weight*x)+h*h*G*np.linalg.norm(x)/2)
                assert error<=bound+3e-14
                packets.append(dict(packet=kind,h=h,full_emission_wavefunction_error=error,
                    weighted_bound=float(bound),S7_operator_bound=g*sqrt(M)/N*sigma*h))
        rows.append(dict(counts=counts,physical_modes=L,physical_dimension=full.size,
            protected_dimension=small.size,weighted_residual_identity_error=identity_error,
            commutator_norm=comm_norm,commutator_upper=G,packets=packets))
    return dict(status='PASS',post_campaign=True,used_in_S7_trajectories=False,
        primary_formation_runs=False,physical_B_implementation='immutable S5 Python full-Fock operator',
        arithmetic='ordinary double diagnostics',checks=rows)


if __name__=='__main__':
    result=run();Path('WEIGHTED_CONTROLS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
