"""Internal graph records in the unrecorded, unforced pure state.

Trace graph variables against fermion/electric variables in their common
kinematic embedding. Every expanded amplitude is a physical Gauss state.
No measuring apparatus or graph-distance record is inserted in evolution.
"""
from pathlib import Path
from math import comb,log2,sqrt
import argparse,json,sys,time
import numpy as np
from scipy.sparse import coo_matrix

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cp12_unmasked_formation'))
sys.path.insert(0,str(P.parent/'cp12_distance_history'))
from run import model
from evolution import Evolution
from dynamics import resources


def entropy(vals):
    x=np.asarray(vals,float);x=x[x>0]
    return float(-sum(x*np.log2(x)))


def expansion(sub,states,sizes,r):
    graphs={};env={};rows=[];cols=[];orbit=[];weights=[];er=[]
    for j,s in enumerate(states):
        seen={}
        for p in sub.group:
            t,sign=sub.transform(s,p)
            if t in seen:assert seen[t]==sign
            else:seen[t]=sign
        assert len(seen)==sizes[j]
        for (g,b,e),sign in seen.items():
            if g not in graphs:graphs[g]=len(graphs)
            key=(b,e)
            if key not in env:env[key]=len(env);er.append(int(r[j]))
            assert er[env[key]]==r[j]
            rows.append(graphs[g]);cols.append(env[key]);orbit.append(j);weights.append(sign/sqrt(sizes[j]))
    rows=np.asarray(rows);cols=np.asarray(cols);orbit=np.asarray(orbit);weights=np.asarray(weights)
    assert np.max(abs(np.bincount(orbit,weights=weights**2,minlength=len(states))-1))<1e-12
    return rows,cols,orbit,weights,np.asarray(er),len(graphs),len(env)


def run(n,q,budget,times):
    started=time.monotonic();sub,H,B,states,sizes,psi=model(n,q,budget)
    r,k,nv,e2=resources(n,q,sub,states)
    rr,cc,jj,weights,er,G,C=expansion(sub,states,sizes,r)
    assert G==comb(sub.M,q)
    capacity=np.array([log2(comb(sub.M-q+int(v),int(v))) for v in er])
    E=Evolution(H,B,q/2,8);v=psi.copy();err=0.;last=0.;out=[]
    for t in times:
        for _ in range(round((t-last)/.1)):
            v,_,de=E(v,.1);err+=de
        last=t
        A=coo_matrix((weights*v[jj],(rr,cc)),shape=(G,C)).tocsr()
        rho=(A@A.conj().T).toarray()
        assert abs(np.trace(rho)-1)<1e-8
        vals=np.linalg.eigvalsh(rho).clip(0,None)
        pg=rho.diagonal().real
        penv=np.asarray(abs(A).power(2).sum(axis=0)).ravel()
        entries=abs(weights*v[jj])**2
        positive=entries>0
        Hcond=float(-sum(entries[positive]*np.log2(entries[positive]/penv[cc[positive]])))
        Hbound=float(penv@capacity)
        assert Hcond<=Hbound+1e-8
        A0=A@__import__('scipy').sparse.diags((er==0).astype(float))
        rho0=(A0@A0.conj().T).toarray()
        p0=float(np.trace(rho0).real)
        off0=rho0-np.diag(rho0.diagonal())
        assert np.max(abs(off0))<1e-13
        d=min(1.,err)
        S_error=(d*log2(G-1)+entropy([d,1-d])) if d<=1-1/G else log2(G)
        out.append(dict(T=t,graph_entropy_bits=entropy(vals),
            graph_entropy_error_upper=S_error,graph_purity=float(np.trace(rho@rho).real),
            graph_classical_entropy_bits=entropy(pg),
            graph_coherence_l1_normalized=float((sum(abs(rho).ravel())-np.trace(rho).real)/(G-1)),
            classical_graph_uncertainty_given_full_configuration_bits=Hcond,
            classical_conditional_entropy_support_bound_bits=Hbound,
            zero_removability_probability=p0,
            zero_removability_graph_offdiagonal_max=float(np.max(abs(off0))),
            state_error_upper=err,total_system_entropy_bits=0.,
            graph_environment_mutual_information_bits=2*entropy(vals)))
    result=dict(N=n,Q=q,budget=budget,graph_count=G,
                expanded_physical_states=len(rr),environment_configurations=C,
                rows=out,seconds=time.monotonic()-started,
                scope='Graph/environment partial trace of a pure physical state in the specified kinematic embedding; no external record inserted. The r=0 sector is analyzed as a diagnostic, not selected in the main evolution. Entropy errors use trace-distance continuity and analytic cutoff bounds, evaluated in ordinary floating point.')
    (P/f'RECORDS_N{n}_Q{q}_B{budget}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True);return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True)
    ap.add_argument('--q',type=int,required=True);ap.add_argument('--budget',type=int,required=True)
    a=ap.parse_args();run(a.n,a.q,a.budget,[0.,1.,3.,5.,10.])
