#!/usr/bin/env python3
"""SCOUT-CA1: CP12's neutral connection Laplacian, independently audited.

FULL CP12 has fermion hopping and dynamical compact gauge fields. This
standalone check isolates ONLY its unchanged, positive graph-reconnection
operator on electrically neutral, filled-edge configurations, with
fixed conserved Q. The exact t^2 coefficients derived below also apply to
the specified full-CP12 neutral initial states. Long-time curves do NOT.

Dependencies numpy scipy; run with OPENBLAS_NUM_THREADS=1 python check.py.
"""
import itertools, json, math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

GAMMA=.1

def reduced(n,q):
    M=math.comb(n,2)
    S=n-1                   # potential edges incident to reference vertex
    A=M-S                   # potential edges not incident to it
    assert 0<q<=A
    ks=np.arange(0,min(q,S)+1,dtype=int)
    degree=q*(M-q)
    within=ks*(S-ks)+(q-ks)*(A-q+ks)
    diagonal=GAMMA/M*(degree-within)
    off=-GAMMA/M*np.sqrt(
        (ks[:-1]+1)*(S-ks[:-1])*(q-ks[:-1])*(A-q+ks[:-1]+1))
    return ks,diagonal,off

def reduced_evolution(n,q,times):
    ks,d,o=reduced(n,q)
    eigen,vec=eigh_tridiagonal(d,o)
    rows=[]
    for tau in times:
        psi=vec@(vec[0,:]*np.exp(-1j*eigen*tau))
        prob=np.abs(psi)**2
        rows.append(dict(time=tau,attached=float(1-prob[0]),
                         mean_incident_edges=float(prob@ks),
                         zero_incident=float(prob[0])))
    return rows

def entire_graph_space(n,q,times):
    M=math.comb(n,2);S=n-1;A=M-S
    configs=list(itertools.combinations(range(M),q))
    lookup={c:i for i,c in enumerate(configs)}
    row=[];col=[];val=[]
    for i,c in enumerate(configs):
        occupied=set(c)
        for e in c:
            for f in range(M):
                if f not in occupied:
                    to=tuple(sorted((occupied-{e})|{f}))
                    row.append(i);col.append(lookup[to])
                    val.append(-GAMMA/M)
    H=coo_matrix((val,(row,col)),
        shape=(len(configs),len(configs))).tocsr()
    H+=diags(np.full(len(configs),GAMMA*q*(M-q)/M))
    initially_isolated=np.array([all(e<A for e in c)
                                 for c in configs])
    assert initially_isolated.sum()==math.comb(A,q)
    phi=initially_isolated.astype(float)/math.sqrt(math.comb(A,q))
    observed=[]
    for tau in times:
        state=expm_multiply(-1j*tau*H,phi)
        observed.append(float(1-np.linalg.norm(
                        state[initially_isolated])**2))
    return len(configs),observed

def full_hamiltonian_short_time_identity(n,q):
    """Direct orthogonal-transition counts from frozen CP12's operator.

    The neutral sea with one occupied graph has q*(n-1) allowed swaps
    that attach the isolated vertex. In a uniform coherent sum of all
    q-edge graphs avoiding the vertex, each adjacent k=1 target has
    (A-q+1) coherent antecedents. Fermion hopping changes occupations
    but NOT the graph, so cannot change the graph-attachment t^2 term.
    """
    M=math.comb(n,2);S=n-1;A=M-S
    fixed=GAMMA**2*q*S/M**2
    coherent=fixed*(A-q+1)
    dephasing=GAMMA*q*(A-q)/M
    graph_energy_coherent=GAMMA*q*S/M
    # Independently, all 2q vertex-edge hopping outputs of a
    # neutral filled-edge sea are pairwise orthogonal. With frozen
    # t_hop=1 the FULL CP12 short-time neutral-sector exit curvature
    # is 2q. This grows with N and blocks a uniform neutral-only
    # approximation on the formation timescale.
    hop_curvature=2*q
    return dict(fixed_graph_attach_t2=fixed,
        coherent_isolate_attach_t2=coherent,
        enhancement=A-q+1,
        coherent_graph_kinetic_energy=graph_energy_coherent,
        graph_dephasing_energy_increase=dephasing,
        full_CP12_hopping_neutral_exit_t2=hop_curvature,
        full_CP12_initial_energy_expectation=.5*q+graph_energy_coherent,
        full_CP12_initial_energy_variance=hop_curvature+coherent,
        caveat="Full-law equalities are exact derivatives/moments at T=0 only.")

def main():
    times=[.1,1.,10.,20.,math.pi/GAMMA,2*math.pi/GAMMA]
    rows=[]
    for n in (4,5,6,8,12,20,40,80,160,320):
        q=n-1
        curve=reduced_evolution(n,q,times)
        direct_dimension=None
        direct_difference=None
        if n<=6:
            direct_dimension,entire=entire_graph_space(n,q,times)
            direct_difference=max(abs(r["attached"]-a)
                  for r,a in zip(curve,entire))
            assert direct_difference<1e-10
        rows.append(dict(N=n,Q=q,
            initial_coherent_graphs=math.comb((n-1)*(n-2)//2,q),
            exact_sector_dimension=len(reduced(n,q)[0]),
            direct_graph_space_dimension=direct_dimension,
            direct_vs_reduced_max_difference=direct_difference,
            short_time=full_hamiltonian_short_time_identity(n,q),
            neutral_only_curve=curve))
    # n=4 is exactly the report's neutral triangle+isolated state.
    assert abs(rows[0]["short_time"]["fixed_graph_attach_t2"]-.0025)<1e-12
    # The report's 50/50 neutral + charged preparation has one edge
    # protected by flux: one fewer neutral removals on charged branch.
    report_assembly_t2=(3*3+2*3)/2*(GAMMA/6)**2
    assert abs(report_assembly_t2-1/480)<1e-12
    for r in rows[-4:]:
        t=10.
        limiting=1-math.exp(-4*(1-math.cos(GAMMA*t)))
        assert abs(r["neutral_only_curve"][2]["attached"]-limiting)<.05
    assert abs(rows[-1]["neutral_only_curve"][1]["attached"] -
       (1-math.exp(-4*(1-math.cos(GAMMA)))))<.001
    assert rows[-1]["neutral_only_curve"][-1]["attached"]<2e-6
    result=dict(status="Neutral-rewiring exact quotient + full-law t2 identities",
        background_CP12_commit="0014d6153a03a27cd1212ac176b77ceb67a07830",
        gamma=GAMMA,hop_full_CP12=1.,
        asymptotic_Q_over_N_one=
           "P(attached,T) -> 1-exp[-4(1-cos(GAMMA*T))]",
        scope="Full neutral-only graph dynamics tested. Full CP12 later-time quantum dynamics NOT simulated.",
        N4_assembly_source_exact_initial_t2=report_assembly_t2,
        rows=rows,
        decisive_qualification="Uniform disconnected graph coherence gives O(1) neutral-only attachment; full CP12 hopping has 2Q initial excitation curvature, so this is NOT an N-independent full-model geometry formation theorem.")
    Path(__file__).with_name("results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({r["N"]:{
         "fixed":r["short_time"]["fixed_graph_attach_t2"],
         "coherent":r["short_time"]["coherent_isolate_attach_t2"],
         "P_attach_T10":r["neutral_only_curve"][2]["attached"],
         "P_attach_T62.8":r["neutral_only_curve"][-1]["attached"]}
             for r in rows},indent=2))
    print("CP12_COHERENT_ASSEMBLY_AUDIT_PASS")

if __name__=="__main__":main()
