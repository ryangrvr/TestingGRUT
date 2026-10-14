"""Unforced, topology-neutral CP12 flux-protection dynamics.

The Hamiltonian and initial ray are exactly the previous frozen law and
uniform neutral sea. Electric-square budgets are projections, with ALL
outgoing transitions retained for a full-law Duhamel error estimate.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time
from math import lcm, sqrt
import numpy as np
from scipy.sparse import coo_matrix, diags

P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent / 'cp12_unmasked_formation'))
sys.path.insert(0, str(P.parent / 'cp12_distance_history'))
from run import model, graphstats
from evolution import Evolution
from exact_destination import transitions


def resources(n, q, sub, states):
    rows = []
    for g, bits, fields in states:
        nv = (bits & ((1 << n)-1)).bit_count()
        r = sum(bool(g >> e & 1) and bool(bits >> (n+e) & 1)
                and fields[2*e] == fields[2*e+1] == 0
                for e in range(sub.M))
        k = sum(bool(g >> e & 1) and bool(bits >> (n+e) & 1)
                and (fields[2*e] != 0 or fields[2*e+1] != 0)
                for e in range(sub.M))
        e2 = sum(x*x for x in fields)
        assert r+k+nv == q
        assert e2 >= 2*k+nv
        rows.append((r, k, nv, e2))
    return np.asarray(rows, float).T


def rewire_matrix(n, q, sub, states, sizes, H=None):
    if H is not None:
        # Hopping changes total vertex number by exactly one; rewiring and
        # diagonal terms preserve it. This decomposition remains exact even
        # when graph masks fold under the signed permutation quotient.
        nv=np.array([(s[1]&((1<<n)-1)).bit_count() for s in states])
        e2=np.array([sum(x*x for x in s[2]) for s in states])
        coo=H.tocoo();keep=nv[coo.row]==nv[coo.col]
        return coo_matrix((coo.data[keep],(coo.row[keep],coo.col[keep])),shape=H.shape).tocsr()-diags(q/2+2*nv+e2/8)
    # Independent slow transition assembly retained as a small-system control.
    index = {s:i for i,s in enumerate(states)}
    rr=[]; cc=[]; vv=[]; scale=lcm(8, 10*sub.M)
    for j,s in enumerate(states):
        for target,c in transitions(n,s,scale,hopping=False):
            key, sign, size=sub.canonical(target)
            if sign:
                assert key in index
                rr.append(index[key]);cc.append(j)
                vv.append(c/scale*sign*sqrt(sizes[j]/size))
    no_hop=coo_matrix((vv,(rr,cc)),shape=(len(states),)*2).tocsr()
    nv=np.array([(s[1]&((1<<n)-1)).bit_count() for s in states])
    e2=np.array([sum(x*x for x in s[2]) for s in states])
    return no_hop-diags(q/2+2*nv+e2/8)


def run(n,q,budget,end,step,steps,build_only=False):
    started=time.monotonic()
    sub,H,B,states,sizes,psi=model(n,q,budget)
    metadata=dict(N=n,Q=q,budget=budget,dimension=len(states),
                  outside_orbits=B.shape[0],nnz=H.nnz)
    if build_only:
        print(json.dumps(metadata),flush=True);return metadata
    r,k,nv,e2=resources(n,q,sub,states)
    R=rewire_matrix(n,q,sub,states,sizes,H=H)
    comm=(R@diags(r)-diags(r)@R).tocsr()
    comm_norm=float(np.linalg.norm(comm.data))
    assert comm_norm<1e-10
    stats={g:graphstats(n,g,sub.es) for g in {s[0] for s in states}}
    gs=np.array([stats[g] for g,_,_ in states],float).T
    shape_names=['connected','triangles','squares','short_cycle_packing',
                 'degree_square_sum']+[f'ball_{j}' for j in range(n)]+[f'shell_{j}' for j in range(n)]
    graph_order=sorted(stats)
    graph_id={g:i for i,g in enumerate(graph_order)}
    graph_ids=np.array([graph_id[s[0]] for s in states])
    initial_graph=np.bincount(graph_ids,weights=abs(psi)**2,minlength=len(graph_order))
    low=(r <= .1*q).astype(float)
    vacuum_flux=((nv == 0)&(k>0)).astype(float)
    evo=Evolution(H,B,q/2,steps)
    v=psi.copy();error=0.;rows=[]
    count=round(end/step)
    assert abs(count*step-end)<1e-10
    for j in range(count+1):
        p=abs(v)**2
        frac_error=min(1.,2*error+error**2)
        rate=float(-2*np.vdot(H@v,r*v).imag)
        pg=np.bincount(graph_ids,weights=p,minlength=len(graph_order))
        row=dict(T=j*step,removable_fraction=float(p@r/q),
                 flux_protected_fraction=float(p@k/q),
                 hole_protected_fraction=float(p@nv/q),
                 low_removability_probability=float(p@low),
                 vacuum_flux_probability=float(p@vacuum_flux),
                 removable_derivative=rate,
                 graph_TV_from_initial=float(abs(pg-initial_graph).sum()/2),
                 graph_TV_error_upper=min(1.,error+error**2/2),
                 electric_square_projected=float(p@e2),
                 mass_energy_projected=float(2*p@nv),
                 hopping_energy_projected=float(np.vdot(v,(H-R)@v).real-q/2-2*p@nv-p@e2/8),
                 rewiring_energy_projected=float(np.vdot(v,R@v).real),
                 state_error_upper=error,bounded_fraction_error_upper=frac_error,
                 norm_error=float(abs(p.sum()-1)),
                 energy_error_projected=float(abs(np.vdot(v,H@v).real-q/2)),
                 geometry=dict(zip(shape_names,map(float,gs@p))))
        row['removable_fraction_interval']=[max(0.,row['removable_fraction']-frac_error),
                                            min(1.,row['removable_fraction']+frac_error)]
        rows.append(row)
        if j%max(1,round(1/step))==0 or j==count:
            print(json.dumps(dict(stage='time',N=n,Q=q,budget=budget,**{a:row[a] for a in ['T','removable_fraction','flux_protected_fraction','low_removability_probability','state_error_upper']})),flush=True)
        if j<count:
            v,_,local=evo(v,step)
            error+=local
    # On the full physical basis, hopping has at most 2Q unit-amplitude
    # transitions per row and each changes r or k by at most one. Thus
    # |d<r/Q>/dt| and |d<k/Q>/dt| <= 2|t_hop| = 2. Midpoint distance
    # to a sampled time is <= step/2. This certifies whole intervals.
    late=[row for row in rows if row['T']>=min(3.,end)]
    continuous=dict(normalized_count_Lipschitz_bound=2.,
        removable_fraction_lower_over_entire_window=max(0.,min(row['removable_fraction']-row['bounded_fraction_error_upper'] for row in rows)-step),
        persistent_flux_window=[late[0]['T'],end],
        flux_protected_fraction_lower_over_late_window=max(0.,min(row['flux_protected_fraction']-row['bounded_fraction_error_upper'] for row in late)-step))
    result=dict(**metadata,initial_energy=q/2,initial_variance=2*q,
                initial_system_entropy_bits=0,initial_graph_count=int(__import__('math').comb(sub.M,q)),
                time_step=step,leakage_subintervals=steps,
                rewiring_removability_commutator_norm=comm_norm,rows=rows,
                continuous_time_bounds=continuous,
                seconds=time.monotonic()-started,
                scope='Unrecorded topology-neutral supplied pure ray, unchanged full law. Bounded-observable cutoff errors are analytically justified and numerically evaluated, without interval arithmetic. Projected unbounded electric-energy diagnostics are not certified full-law expectations.')
    out=P/f'N{n}_Q{q}_B{budget}_T{end:g}.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='saved',path=out.name,seconds=result['seconds'])),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--n',type=int,required=True);ap.add_argument('--q',type=int,required=True)
    ap.add_argument('--budget',type=int,required=True);ap.add_argument('--end',type=float,default=10)
    ap.add_argument('--step',type=float,default=.1);ap.add_argument('--steps',type=int,default=8)
    ap.add_argument('--build-only',action='store_true')
    a=ap.parse_args();run(a.n,a.q,a.budget,a.end,a.step,a.steps,a.build_only)
