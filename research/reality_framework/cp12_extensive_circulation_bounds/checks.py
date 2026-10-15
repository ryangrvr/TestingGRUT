"""Independent finite-matrix and rational checks of the new analytic bounds."""
from fractions import Fraction as F
from math import sqrt
import json,time
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import svds
from bounds import P,residual_probability
from prepared import FROZEN
from frozen_model import build
from exact_destination import transitions


def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def add(x,y):return x[0]+y[0],x[1]+y[1]
def times(x,c):return x[0]*c,x[1]*c


def rational_control():
    one=(F(1),F(0));eta=(F(60,291),F(-32,291));a=F(19,8)
    polynomial=add(add(times(mul(eta,eta),F(873)),times(eta,F(-360))),(F(16),F(0)))
    assert polynomial==(0,0)
    bright=times(add((F(1,2),F(0)),times(eta,a)),F(1,2))
    assert mul(bright,bright)==times(mul(eta,add(one,times(eta,-1))),2)
    local_mean=add(times(eta,a),times(bright,-2));assert local_mean==(F(-1,2),0)
    local_second=add(add((F(2),F(0)),times(eta,a*a)),times(bright,-2*a))
    assert local_second==(F(13,16),0)
    assert add(local_second,times(mul(local_mean,local_mean),-1))==(F(9,16),0)
    assert F(31,16)+F(33,16)==4
    assert F(33,16*127)<F(51,400)**2
    assert F(12247,10000)**2<F(3,2)
    assert F(35,32)+F(51,400)<F(12247,10000)
    return dict(quadratic_field='Q(sqrt(2))',neutral_mean_energy_identity_exact=True,
        local_second_moment_exact='13/16',local_variance_exact='9/16',
        polynomial_identity_exact=True,
        two_orientation_variance_match_exact='31/16 + (N-1) * 33/[16(N-1)] = 4',
        phase_cosines_valid_for_all_N_at_least_128_by_rational_inequalities=True)


def physical_blocks():
    H,states,lookup,degrees=build(4,4,3)
    nv=np.array([(s[1]&15).bit_count() for s in states]);co=H.tocoo()
    hopmask=nv[co.row]!=nv[co.col]
    from scipy.sparse import coo_matrix
    K=coo_matrix((co.data[hopmask],(co.row[hopmask],co.col[hopmask])),shape=H.shape).tocsr()
    rows=[]
    for m in range(4):
        block=K[np.flatnonzero(nv==m+1)][:,np.flatnonzero(nv==m)]
        norm=float(svds(block,k=1,which='LM',return_singular_vectors=False)[0]) if block.nnz else 0.
        bound=sqrt(8*(m+1));assert norm<=bound+1e-10
        rows.append(dict(vertex_fermions=m,raising_block_norm=norm,global_unregulated_norm_bound=bound))
    # Direct unregulated hopping images: every vacuum parent has exactly 2Q
    # distinct outgoing targets, and no two vacuum parents share a target.
    targets={};vacuum=0
    for s in states:
        if s[1]&15:continue
        vacuum+=1;local=[]
        for target,c in transitions(4,s,120):
            if target[1]&15:
                assert abs(c)==120 and target not in targets
                targets[target]=s;local.append(target)
        assert len(local)==8
    return dict(raw_physical_dimension=len(states),raising_blocks=rows,
        unregulated_vacuum_parent_columns=vacuum,unregulated_distinct_hopping_targets=len(targets),
        exact_vacuum_hopping_Gram='8 I',vacuum_columns_have_unique_targets=True)


def residual_dilations():
    rng=np.random.default_rng(941237);cases=0;nonvacuous=0;worst=0.
    for j in range(24):
        k=3;d=8;A=rng.normal(size=(k,k));pp=(4+j/8)*np.eye(k)+A.T@A/5
        B=rng.normal(size=(k,d-k))/5;C=rng.normal(size=(d-k,d-k));cc=(C+C.T)/5
        H=np.block([[pp,B],[B.T,cc]])
        x=np.r_[np.zeros(k),rng.normal(size=d-k)+1j*rng.normal(size=d-k)];x/=np.linalg.norm(x)
        E=float(np.vdot(x,H@x).real);sigma=float(np.linalg.norm((H-E*np.eye(d))@x))
        delta=float(np.linalg.eigvalsh(pp)[0]-E);b=float(np.linalg.norm(B,2))
        bound=residual_probability(delta,sigma,b)
        if bound<1:nonvacuous+=1
        for t in [0.,.1,1.,10.,100.]:
            y=expm(-1j*t*H)@x;p=float(np.linalg.norm(y[:k])**2)
            worst=max(worst,p-bound);assert p<=bound+1e-10;cases+=1
    return dict(independent_matrices=24,time_evolution_cases=cases,
        matrices_with_nonvacuous_bound=nonvacuous,maximum_probability_bound_violation=worst)


def run():
    start=time.monotonic();out=dict(status='CP12_EXTENSIVE_BOUND_CHECKS_PASS',
        exact_algebra=rational_control(),physical=physical_blocks(),residual_probability=residual_dilations(),
        seconds=time.monotonic()-start,
        scope='Exact rational identities and unregulated vacuum hopping columns; finite-matrix numerical norm and residual-bound checks independently test analytic steps, not the large-N theorem by numerical extrapolation.')
    (P/'CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True);return out


if __name__=='__main__':run()
