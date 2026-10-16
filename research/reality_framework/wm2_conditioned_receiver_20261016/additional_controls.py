"""Common-frame limit, larger matrix-free maps and exact storage counts."""
from math import comb,sqrt
import json
from pathlib import Path
import numpy as np
import fock as F
import conditioned as C
import controls


def run():
    rng=np.random.default_rng(623019)
    cache_checks=[]
    for rank,k in [(4,3),(9,4)]:
        u=controls.frame(rng,rank,rank)
        x=rng.normal(size=F.space(rank,k).size)+1j*rng.normal(size=F.space(rank,k).size)
        cached=F.apply_unitary(F.unitary_plan(u,k,True),x,rank,k)
        uncached=F.apply_unitary(F.unitary_plan(u,k,False),x,rank,k)
        assert np.array_equal(cached,uncached)
        cache_checks.append(dict(rank=rank,photons=k,bit_identical=True))
    stress=[]
    for rd,rs,k in [(9,10,4),(10,9,4),(10,10,3)]:
        wd,ws=controls.frame(rng,14,rd),controls.frame(rng,14,rs)
        g=wd.conj().T@ws
        x=rng.normal(size=F.space(rs,k).size)+1j*rng.normal(size=F.space(rs,k).size)
        x/=np.linalg.norm(x)
        dense=F.symmetric_power(g,k)@x
        other=F.Map(g,k,dense_limit=0)(x)
        error=float(np.linalg.norm(dense-other))
        assert error<5e-11
        stress.append(dict(source_rank=rs,destination_rank=rd,photons=k,error=error))
    common=[]
    for counts,values,rank in [([1,1],[-.25,.25],3),([2,1],[-.25,.25],4),([2,2],[-.25,.25],3)]:
        e=C.Evolution(counts,values,8)
        q=e.contact[:,None].astype(complex)
        z=controls.frame(rng,8,rank-1)
        w=np.linalg.qr(np.column_stack((q,z)))[0]
        e.frames={u:w.copy() for u in e.words}
        e.x={u:rng.normal(size=F.space(rank,len(u)).size)+1j*rng.normal(size=F.space(rank,len(u)).size)
             for u in e.words}
        e.normalize()
        layout=C.S5.Layout(counts,values,rank)
        flat=np.zeros(layout.size,complex)
        for i,u in enumerate(layout.words):flat[layout.offset[i]:layout.offset[i]+len(e.x[u])]=e.x[u]
        s5=layout.hint(flat,w.conj().T@e.contact,.9)
        op=C.Operator(e.words,e.frames,e.edges,e.contact,.9)
        actual=op.hint(e.x)
        error=sqrt(sum(float(np.linalg.norm(actual[u]-s5[layout.offset[i]:layout.offset[i]+len(e.x[u])])**2)
                       for i,u in enumerate(layout.words)))
        residual=op.residual(e.x,actual)
        assert error<3e-11 and residual<5e-8
        common.append(dict(N=e.N,rank=rank,coefficient_dimension=layout.size,
                           common_frame_S5_operator_error=error,full_residual=residual))
    cost=[]
    for n in (4,5,6,7,8):
        counts=[(n+1)//2-1,n//2]
        words=C.S5.type_words(counts)
        cost.append(dict(N=n,nonroot_counts=counts,histories=len(words),
            coefficient_ceiling_if_all_ranks_10=sum(comb(10+len(u)-1,len(u)) for u in words),
            frame_entries_if_all_nonempty_ranks_10=64*n*n*10*(len(words)-1),
            conditional_rank_is_not_proved_sufficient=True))
    out=dict(rotation_cache_identity=cache_checks,matrix_free_stress=stress,common_frame_limit=common,exact_storage_counts=cost,
             external_validation=False,roundoff_certified=False)
    Path(__file__).with_name('ADDITIONAL_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    print(json.dumps(run(),indent=2))
