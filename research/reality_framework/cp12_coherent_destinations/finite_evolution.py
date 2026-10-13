"""Full frozen CP-12 finite-time coherent destination control.

Primary sector N=5,Q=4. All graphs, all Gauss-law fermion states, electric
cutoffs L=2,3. Initial vacancy coherence and its square-tail reference are
supplied. No postselection, modified Hamiltonian, or inferred metric.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
import time

import numpy as np
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
EXPECTED='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
assert hashlib.sha256((P/'frozen_model.py').read_bytes()).hexdigest()==EXPECTED
sys.path.insert(0,str(P))
from frozen_model import build


def main(L,end):
    start=time.monotonic()
    n=5;Q=4
    es=list(itertools.combinations(range(n),2));M=len(es)
    pool=sum(1 << es.index(tuple(sorted(e))) for e in [(0,1),(1,2),(2,3),(3,0),(3,4)])
    initial=[pool ^ (1 << e) for e in range(M) if pool >> e & 1]
    absent=[e for e in range(M) if not (pool >> e & 1)]
    H,states,lookup,_=build(n,Q,L)
    print(json.dumps(dict(stage='built',L=L,dimension=len(states))),flush=True)
    sea=sum(1 << (n+e) for e in range(M));zero=(0,)*(2*M)
    ci=len(initial);ui=ci+1
    psi=np.zeros((len(states),len(initial)+2),complex)
    for j,g in enumerate(initial):
        psi[lookup[g,sea,zero],j]=1
    psi[:,ci]=psi[:,:ci].sum(axis=1)/np.sqrt(len(initial))
    neutral_ids=[i for i,s in enumerate(states) if s[1]==sea and s[2]==zero]
    psi[neutral_ids,ui]=1/np.sqrt(len(neutral_ids))
    Hpsi=H@psi
    energy=np.einsum('ij,ij->j',psi.conj(),Hpsi).real
    variance=np.sum(abs(Hpsi)**2,axis=0)-energy**2
    assert max(abs(energy[:ci]-2.24))<1e-12 and abs(energy[ci]-2.2)<1e-12
    assert abs(variance[ci]-8.004)<1e-12
    assert abs(energy[ui]-2.)<1e-12 and abs(variance[ui]-8.)<1e-12
    masks=np.array([s[0] for s in states],np.int64)
    readouts=np.array([((masks >> e) & 1) for e in absent],float)
    boundary=np.array([max(map(abs,s[2]))==L for s in states],float)
    nv=np.array([(s[1]&((1 << n)-1)).bit_count() for s in states],float)
    E2=np.array([sum(x*x for x in s[2]) for s in states],float)
    graphs=sorted(set(masks.tolist()))
    shape={}
    for g in graphs:
        adj=[set() for _ in range(n)]
        for e,(a,b) in enumerate(es):
            if g >> e & 1:adj[a].add(b);adj[b].add(a)
        reached={0};todo=[0]
        while todo:
            a=todo.pop()
            for b in adj[a]-reached:reached.add(b);todo.append(b)
        connected=len(reached)==n
        degrees=[len(x) for x in adj]
        shape[g]=[float(connected),float(connected and max(degrees)<=2),
                  sum(d*(d-1)/2 for d in degrees),max(degrees),sum(d==0 for d in degrees)]
    shape_read=np.array([shape[g] for g in masks],float).T
    labels=['connected','spanning_path','degree_wedges','maximum_degree','isolated_vertices']
    uniform_baseline=np.mean(np.array([shape[g] for g in graphs]),axis=0)
    assert abs(uniform_baseline[0]-125/210)<1e-12
    assert abs(uniform_baseline[1]-60/210)<1e-12
    assert abs(uniform_baseline[2]-4.)<1e-12
    # Neutral-only projection is a declared comparator, not the full law.
    Hn=H[neutral_ids][:,neutral_ids]
    vn=psi[neutral_ids,ci].copy()
    neutral_read=readouts[:,neutral_ids]
    del states,lookup,Hpsi
    rows=[];last=0.
    for T in [.02,.1,.25,.5,1.,2.,4.,8.,16.]:
        if T>end+1e-10:continue
        steps=max(1,int(np.ceil(T-last)))
        dt=(T-last)/steps
        for _ in range(steps):
            psi=expm_multiply(-1j*dt*H,psi)
        vn=expm_multiply(-1j*(T-last)*Hn,vn);last=T
        prob=abs(psi)**2
        observations=readouts@prob
        coherent=observations[:,ci]
        dephased=observations[:,:ci].mean(axis=1)
        neutral=neutral_read@abs(vn)**2
        assert np.ptp(neutral)<2e-12
        norm_error=float(np.max(abs(prob.sum(axis=0)-1)))
        current_energy=np.einsum('ij,ij->j',psi.conj(),H@psi).real
        energy_error=float(np.max(abs(current_energy-energy)))
        assert norm_error<2e-9 and energy_error<2e-8
        which=[]
        for eta in [0.,.5,1.]:
            values=eta*coherent+(1-eta)*dephased
            which.append(dict(record_overlap=eta,mean_initial_system_energy=2.24-.04*eta,
                              occupations=values.tolist(),
                              normalized_destination_occupation_distribution=(values/values.sum()).tolist()))
        row=dict(T=T,coherent_occupations=coherent.tolist(),dephased_occupations=dephased.tolist(),
                 neutral_only_occupations=neutral.tolist(),record_visibility_controls=which,
                 coherent_cutoff_boundary=float(boundary@prob[:,ci]),
                 dephased_cutoff_boundary=float((boundary@prob[:,:ci]).mean()),
                 coherent_charge_density=float(nv@prob[:,ci]/n),
                 coherent_electric_square_per_vertex=float(E2@prob[:,ci]/n),
                 norm_error=norm_error,energy_error=energy_error)
        shape_values=shape_read@prob[:,ui]
        all_edge_means=np.array([((masks >> e)&1)@prob[:,ui] for e in range(M)])
        assert np.ptp(all_edge_means)<2e-12
        assert max(abs(all_edge_means-Q/M))<2e-12
        row['unordered_pregeometry_control']=dict(
            **dict(zip(labels,shape_values.tolist())),
            all_labelled_edges_mean_occupation=float(all_edge_means.mean()),
            edge_mean_spread=float(np.ptp(all_edge_means)),
            cutoff_boundary=float(boundary@prob[:,ui]),
            charge_density=float(nv@prob[:,ui]/n),
            electric_square_per_vertex=float(E2@prob[:,ui]/n))
        # Specified equal-degree destinations: 04 is distance two, 14 distance three in F.
        ia=absent.index(es.index((0,4)));ib=absent.index(es.index((1,4)))
        row['near_minus_far_coherent']=float(coherent[ia]-coherent[ib])
        row['near_minus_far_dephased']=float(dephased[ia]-dephased[ib])
        row['near_minus_far_over_T10']=row['near_minus_far_coherent']/T**10
        rows.append(row)
        print(json.dumps(dict(L=L,T=T,near_minus_far=row['near_minus_far_coherent'],
                              boundary=row['coherent_cutoff_boundary'])),flush=True)
    out=dict(status='CP12_COHERENT_DESTINATION_FINITE_PASS',N=n,Q=Q,L=L,
             dimension=H.shape[0],pool_edges=[es[e] for e in range(M) if pool >> e & 1],
             common_absent_edges=[es[e] for e in absent],initial_graph_count=len(initial),
             initial_energies=energy.tolist(),initial_variances=variance.tolist(),
             unordered_graph_count=len(graphs),unordered_initial_shape=dict(zip(labels,uniform_baseline.tolist())),
             rows=rows,seconds=time.monotonic()-start,
             scope='Full finite frozen Hamiltonian at one connection sector and supplied preparations; no continuum or assembly claim.')
    (P/f'finite_L{L}.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--cutoff',type=int,default=2)
    ap.add_argument('--end',type=float,default=16.)
    a=ap.parse_args();main(a.cutoff,a.end)
