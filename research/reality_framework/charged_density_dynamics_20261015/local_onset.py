"""Exact host-summed local onset for the quarter-density event at N=5..8.

At minimum L1, N5/N6 require one directed triangle and N7/N8 one directed
square, with exactly two holes at distinct positive endpoints. Unused host
edges cannot participate at this order. Each target's SIGNED hopping sum
is computed before squaring; host multiplicities are exact combinations.
"""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations,permutations,product
from math import comb,factorial
from pathlib import Path
import json,sys,time

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'constructive_sprint_20261015/charged_access'))
from cycle_cover import direct_amplitude
sys.path.insert(0,str(P.parent/'cp12_ground_connectivity'))
from exact_destination import gauss


def calculate(n,length):
    start=time.monotonic();q=2*n;es=list(combinations(range(n),2));M=len(es)
    cycle=[(v,(v+1)%length) for v in range(length)]
    cycle_edges={es.index(tuple(sorted(e))) for e in cycle}
    positive={es.index(tuple(sorted((a,b)))):a for a,b in cycle}
    group=[]
    for shift in range(length):
        for outside in permutations(range(length,n)):
            group.append(tuple((v+shift)%length for v in range(length))+outside)
    def key(holes):
        return min(tuple(sorted((min(p[es[e][0]],p[es[e][1]]),
            max(p[es[e][0]],p[es[e][1]]),p[v]) for e,v in holes)) for p in group)
    sea=sum(1<<(n+e) for e in range(M));fields=[0]*(2*M)
    for a,b in cycle:
        e=es.index(tuple(sorted((a,b))));fields[2*e+(a>b)]=1;fields[2*e+(a<b)]=-1
    cases={};total=0
    for pair in combinations(range(M),2):
        options=[(positive[e],) if e in positive else es[e] for e in pair]
        for vertices in product(*options):
            if vertices[0]==vertices[1]:continue
            holes=tuple(zip(pair,vertices));signature=key(holes)
            k=len(cycle_edges|set(pair));parents=comb(M-k,q-k)
            if signature not in cases:
                E=fields.copy();bits=sea
                for e,v in holes:
                    E[2*e+int(v==es[e][1])]+=1
                    bits^=1<<(n+e);bits|=1<<v
                g=sum(1<<e for e in cycle_edges|set(pair))
                gauss(n,(g,bits,tuple(E)))
                amplitude,peak=direct_amplitude(n,(g,bits,tuple(E)))
                cases[signature]=dict(signed_integer_amplitude=str(amplitude),
                    squared_integer_amplitude=amplitude**2,
                    changed_edges=k,compatible_hosts=parents,multiplicity=0,
                    maximum_DP_layer=peak)
            row=cases[signature];assert row['changed_edges']==k and row['compatible_hosts']==parents
            row['multiplicity']+=1
            total+=parents*row['squared_integer_amplitude']
    directed_cycles=factorial(n)//factorial(n-length)//length
    squared_norm=directed_cycles*total;order=2*length+2
    coefficient=Fraction(squared_norm,comb(M,q)*factorial(order)**2)
    rows=[dict(signature=[list(x) for x in s],**r) for s,r in sorted(cases.items())]
    result=dict(N=n,Q=q,c_min=length,nv_min=2,first_amplitude_order=order,
        leading_probability_order=2*order,coefficient=str(coefficient),
        coefficient_float=float(coefficient),integer_squared_norm=str(squared_norm),
        directed_cycles=directed_cycles,neutral_graphs=comb(M,q),motif_classes=len(cases),
        fixed_cycle_target_configurations=sum(x['multiplicity'] for x in rows),rows=rows,
        seconds=time.monotonic()-start)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    return result


def main():
    data=[calculate(n,3 if n<7 else 4) for n in [5,6,7,8]]
    for row in data[:2]:
        saved=json.loads((P/f'ONSET_N{row["N"]}.json').read_text())
        assert row['coefficient']==saved['leading_probability_coefficient']
        assert row['integer_squared_norm']==saved['integer_squared_norm']
    result=dict(status='PASS',cases=data,
        independent_N5_N6_orbit_and_matrix_match=True,
        scope='Exact first-permitted-order coefficient for c/Q>=1/4,nv/N>=1/4 at N5..N8. These integer thresholds are c=3,3,4,4 and nv=2. A triangle/square is the only possible minimal circulating core. All compatible Q-edge hosts are counted, with no supplied degree or dimension. This is not finite-time formation weight or a thermodynamic phase.')
    (P/'LOCAL_ONSET.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
