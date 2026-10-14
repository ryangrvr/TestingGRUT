"""Derive physical summaries from the separately recorded eigensolves.

This script does NOT rerun those eigensolves. Uniform graph counts are freshly
enumerated. Projector bounds refer to the invariant blocks unless stated.
"""
from pathlib import Path
import itertools
import json
import math
from ground import topology

P=Path(__file__).resolve().parent

def uniform_counts(n):
    es=list(itertools.combinations(range(n),2));connected=regular=0
    for occupied in itertools.combinations(range(len(es)),n):
        parent=list(range(n));degrees=[0]*n
        def root(i):
            while parent[i]!=i:i=parent[i]
            return i
        for e in occupied:
            a,b=es[e];degrees[a]+=1;degrees[b]+=1;parent[root(a)]=root(b)
        connected+=len({root(i) for i in range(n)})==1
        regular+=all(d==2 for d in degrees)
    total=math.comb(len(es),n)
    return dict(total=total,connected=connected,degree_two=regular,
                connected_probability=connected/total,degree_two_probability=regular/total)

def main():
    rows=[]
    for n in [4,5,6,7]:
        L=2 if n==7 else 3
        raw=json.loads((P/f'ground_N{n}_L{L}_invariant.json').read_text())
        es=list(itertools.combinations(range(n),2));cycles={};largest={};components={}
        for r in raw['graph_class_probability']:
            flags=topology(n,r['graph_mask'],es);prob=r['probability']
            if flags[0]:
                ell=n-int(flags[6]);cycles[ell]=cycles.get(ell,0)+prob
            size=round(n*flags[2]);largest[size]=largest.get(size,0)+prob
            nc=int(flags[1]);components[nc]=components.get(nc,0)+prob
        pc=raw['observables']['connected'];assert abs(sum(cycles.values())-pc)<1e-10
        row=dict(N=n,Q=n,L=L,scope='all-sector finite-cutoff ground checked' if n<7 else
                 'invariant-sector provisional; neither all-irrep nor cutoff-converged',
            energy=raw['energies'][0],connected_probability=pc,
            largest_component_fraction=raw['observables']['largest_component_fraction'],
            degree_two_probability=raw['observables']['degree_two'],
            single_ring_probability=cycles.get(n,0.),
            connected_cycle_length_probability={str(k):v/pc for k,v in sorted(cycles.items())},
            mean_cycle_length_given_connected=sum(k*v for k,v in cycles.items())/pc,
            largest_component_size_probability=largest,component_count_probability=components,
            uniform_labelled_graph_count=uniform_counts(n))
        if 'sector_comparison' in raw:
            s=raw['sector_comparison'];E=raw['energies'][0]
            ec=s['connected_energy'];ed=s['disconnected_energy'];b=s['cross_norm']
            dc=ec-E;dd=ed-E;assert dc>0 and dd>0
            lower=dd*dd/(dd*dd+b*b);upper=b*b/(dc*dc+b*b)
            scalar_energy_lower=(ec+ed-math.sqrt((ec-ed)**2+4*b*b))/2
            assert lower-1e-9<=pc<=upper+1e-9
            assert scalar_energy_lower-1e-9<=E<=min(ec,ed)+1e-9
            row['block_comparison']=dict(scope='invariant connected/disconnected blocks; their projected ground components stay invariant',
                connected_energy=ec,disconnected_energy=ed,cross_norm=b,
                connected_probability_lower=lower,connected_probability_upper=upper,
                scalar_energy_lower=scalar_energy_lower,
                graph_mixing_energy=s['graph_mixing_energy'],
                connectedness_dephasing_energy_cost=-s['graph_mixing_energy'],
                conditional_connected_energy=s['conditional_connected_energy'],
                conditional_disconnected_energy=s['conditional_disconnected_energy'])
        rows.append(row)
    result=dict(status='CP12_GROUND_CONNECTIVITY_ANALYSIS_PASS',rows=rows,
        exact_uniform_counts=True,eigensolves_rerun_by_this_script=False,
        probabilities_scope='Born probabilities of graph projectors in computed finite-system states; not experimental data')
    (P/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
