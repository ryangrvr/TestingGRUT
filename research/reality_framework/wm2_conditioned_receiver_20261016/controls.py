"""Independent physical lifts for cross-history algebra and dynamics."""
from pathlib import Path
import importlib.util
from math import sqrt
import numpy as np
from scipy.fft import dst
from scipy.linalg import expm
from scipy.sparse import block_diag

import fock as F
import conditioned as C

_path=Path(__file__).resolve().parent.parent/'wm2_history_formation_20261016/controls.py'
_spec=importlib.util.spec_from_file_location('frozen_s5_controls',_path)
S5=importlib.util.module_from_spec(_spec)
# Its explicit `import active` must resolve to the immutable S5 source.
import sys
sys.path.insert(0,str(_path.parent))
_spec.loader.exec_module(S5)
sys.path.pop(0)
M=S5.M


def frame(rng,modes,rank):
    return np.linalg.qr(rng.normal(size=(modes,rank))+1j*rng.normal(size=(modes,rank)))[0]


def fock_checks():
    rng=np.random.default_rng(623016)
    rows=[]
    for rd,rs,k in [(2,3,1),(3,2,2),(3,4,3),(4,3,3),(4,4,4)]:
        wd,ws=frame(rng,6,rd),frame(rng,6,rs)
        exact=S5.field_isometry(wd,k).conj().T@S5.field_isometry(ws,k)
        direct=F.symmetric_power(wd.conj().T@ws,k)
        x=rng.normal(size=exact.shape[1])+1j*rng.normal(size=exact.shape[1])
        x/=np.linalg.norm(x)
        matrixfree=F.Map(wd.conj().T@ws,k,dense_limit=0)(x)
        errors=[float(np.linalg.norm(direct-exact)),float(np.linalg.norm(matrixfree-exact@x))]
        assert max(errors)<3e-11,(rd,rs,k,errors)
        rows.append(dict(destination_rank=rd,source_rank=rs,photons=k,
                         lifted_dense_error=errors[0],matrix_free_error=errors[1]))
    return rows


def lift(e,exact):
    out=np.zeros(len(exact['states']),complex)
    for u in e.words:
        w=dst(e.frames[u],type=1,norm='ortho',axis=0)
        field=S5.field_isometry(w,len(u))@e.x[u]
        for i,o in enumerate(F.space(e.modes,len(u)).occ):
            out[exact['index'][u,o]]=field[i]
    return out


def cross_history_checks():
    rng=np.random.default_rng(623017)
    rows=[]
    for counts,values,modes in [([1,1],[-.25,.25],4),([2,1],[-.25,.25],4),([1,1,1],[-.25,0,.25],3)]:
        e=C.Evolution(counts,values,modes,rank_cap=modes)
        for j,u in enumerate(e.words):
            r=1+(j%(modes-1)) if u else 1
            e.frames[u]=frame(rng,modes,r)
            e.x[u]=rng.normal(size=F.space(r,len(u)).size)+1j*rng.normal(size=F.space(r,len(u)).size)
        e.normalize()
        deltas=[0]+[v for n,v in zip(counts,values) for _ in range(n)]
        exact=M.vacuum(e.N,modes,deltas,.9,classes=True)
        x=lift(e,exact)
        op=C.Operator(e.words,e.frames,e.edges,e.contact,.9)
        compressed=op.hint(e.x)
        old=e.x; e.x=compressed; lifted=lift(e,exact); e.x=old
        full=exact['H_emission']@x
        residual_direct=float(np.linalg.norm(full-lifted))
        residual_formula=op.residual(e.x,compressed)
        action_norm=float(np.linalg.norm(full)**2)
        norm_formula=op.full_action_norm2(e.x)
        # Project each complete physical action into its actual history frame.
        projected={}
        for u in e.words:
            values_physical=np.array([full[exact['index'][u,o]] for o in F.space(modes,len(u)).occ])
            w=dst(e.frames[u],type=1,norm='ortho',axis=0)
            projected[u]=S5.field_isometry(w,len(u)).conj().T@values_physical
        block_error=sqrt(sum(float(np.linalg.norm(projected[u]-compressed[u])**2) for u in e.words))
        assert block_error<3e-11
        assert abs(action_norm-norm_formula)<3e-11
        assert abs(residual_direct-residual_formula)<3e-11
        # Deliberately remove coherent cross terms from the full norm.
        diagonal_only=0.0
        for u in e.words:
            one={v:np.zeros_like(a) for v,a in e.x.items()}; one[u]=e.x[u]
            diagonal_only+=op.full_action_norm2(one)
        # Verify Hermiticity independently on random coefficient packets.
        z={u:rng.normal(size=len(a))+1j*rng.normal(size=len(a)) for u,a in e.x.items()}
        bz=op.hint(z)
        adjoint_error=abs(sum(np.vdot(z[u],compressed[u])-np.vdot(bz[u],e.x[u]) for u in e.words))
        assert adjoint_error<3e-11
        rows.append(dict(N=e.N,modes=modes,histories=len(e.words),coefficient_dimension=sum(len(a) for a in e.x.values()),
            projected_append_pop_error=block_error,full_action_norm_squared_error=abs(action_norm-norm_formula),
            full_residual_direct=residual_direct,full_residual_formula=residual_formula,
            adjoint_error=float(adjoint_error),interference_change=norm_formula-diagonal_only))
    assert max(abs(r['interference_change']) for r in rows)>1e-4
    return rows


def emission_step_checks():
    rng=np.random.default_rng(623018)
    rows=[]
    for h in [.02,.1,.3]:
        e=C.Evolution([1,1],[-.25,.25],4,rank_cap=4)
        for j,u in enumerate(e.words):
            r=1+j%3 if u else 1
            e.frames[u]=frame(rng,4,r)
            e.x[u]=rng.normal(size=F.space(r,len(u)).size)+1j*rng.normal(size=F.space(r,len(u)).size)
        e.normalize()
        exact=M.vacuum(3,4,[0,-.25,.25],.9,classes=True)
        y0=lift(e,exact)
        op=C.Operator(e.words,e.frames,e.edges,e.contact,.9)
        out,rho,tail=op.polynomial(e.x,h)
        e.x=out
        approximated=lift(e,exact)
        y=expm(-1j*h*exact['H_emission'].toarray())@y0
        error=float(np.linalg.norm(y-approximated))
        assert error<=rho+tail+3e-12
        rows.append(dict(time=h,state_error=error,residual_integral_upper=rho,polynomial_tail=tail))
    return rows


def dynamics_checks():
    rows=[]
    for counts,values,modes in [([1,1],[-.25,.25],5),([2,1],[-.25,.25],4)]:
        for duration in [.5,2,8]:
            e=C.Evolution(counts,values,modes,rank_cap=modes)
            r=e.run(duration,strang_budget=.002,projection_budget=.005,residual_budget=.02,max_seconds=60)
            assert r['status']=='COMPLETE',r['stop_reason']
            deltas=[0]+[v for n,v in zip(counts,values) for _ in range(n)]
            exact=M.vacuum(e.N,modes,deltas,.9,classes=True)
            y=M.evolve(exact['H'],M.unit(len(exact['states']),exact['seed']),duration)
            error=float(np.linalg.norm(y-lift(e,exact)))
            assert error<=r['analytic_state_error_account']+3e-10
            rows.append(dict(N=e.N,modes=modes,time=duration,full_dimension=len(y),
                             full_wavefunction_error=error,analytic_state_error_account=r['analytic_state_error_account'],
                             maximum_word_rank=r['maximum_word_rank'],projections=r['projections']))
    return rows


if __name__=='__main__':
    import json
    result={}
    for name in ['fock_checks','cross_history_checks','emission_step_checks','dynamics_checks']:
        result[name]=globals()[name]()
        print(name,'PASS',flush=True)
    Path(__file__).with_name('CONTROLS.json').write_text(json.dumps(result,indent=2)+'\n')
