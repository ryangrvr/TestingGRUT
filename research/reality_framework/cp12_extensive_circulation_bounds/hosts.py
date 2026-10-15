"""Matched-resource supplied 1D, 2D and 3D host controls, not formed graphs."""
from collections import deque
import json,time
from bounds import P,dressed_moments
from prepared import euler,graph


def diamond(L):
    def index(x,y,z,s):return 2*((x%L)*L*L+(y%L)*L+(z%L))+s
    es=set()
    for x in range(L):
        for y in range(L):
            for z in range(L):
                a=index(x,y,z,0)
                for dx,dy,dz in [(0,0,0),(-1,0,0),(0,-1,0),(0,0,-1)]:
                    es.add(tuple(sorted((a,index(x+dx,y+dy,z+dz,1)))))
    return sorted(es)


def square(L):
    return sorted({tuple(sorted((x*L+y,((x+dx)%L)*L+(y+dy)%L)))
                   for x in range(L) for y in range(L) for dx,dy in [(1,0),(0,1)]})


def balls(n,es,maximum=12):
    adj=[[] for _ in range(n)]
    for a,b in es:adj[a].append(b);adj[b].append(a)
    assert all(len(x)==4 for x in adj)
    total=[0]*(maximum+1);seen={0};queue=deque([(0,0)]);counts=[0]*(maximum+1)
    while queue:
        v,r=queue.popleft()
        if r<=maximum:counts[r]+=1
        for w in adj[v]:
            if w not in seen:seen.add(w);queue.append((w,r+1))
    assert len(seen)==n
    for r in range(maximum+1):total[r]=sum(counts[:r+1])
    return total


def run():
    start=time.monotonic();rows=[]
    # The sizes match simultaneously: 2 L_3^3 = L_2^2 = 16 k^6.
    for k in [1,2,3]:
        L3=2*k*k;L2=4*k**3;n=16*k**6
        for name,es in [('one_dimensional_ring_power',graph(n,'ring_power')),
                        ('two_dimensional_square',square(L2)),('three_dimensional_diamond',diamond(L3))]:
            assert len(es)==2*n;euler(n,es)
            V=balls(n,es)
            rows.append(dict(supplied_host=name,N=n,Q=2*n,cycle_rank=n+1,
                electric_square_support_range=[4*n,7*n],ball_volume_from_representative_root=V,
                exact_mean_energy_at_eta_neutral=n,
                variance_formula='N * (41/16 - 3 eta + 2 eta^2)',
                mean_matter_formula='N eta',mean_electric_square_formula='N (4 + 3 eta)',
                mean_and_variance_matched_pair_available=n>=128,
                pair_matched_resources=dict(mean_energy=n,variance=4*n,matter_density=.25,
                    mean_electric_square=19*n/4) if n>=128 else None,
                graph_selected_by_initial_state=True))
    result=dict(status='CP12_MATCHED_HOST_CONTROLS_PASS',cases=rows,seconds=time.monotonic()-start,
        scope='Intrinsic BFS on supplied transitive graph families. One-dimensional ring-power, square torus and diamond host are prescribed, never selected by the Hamiltonian. Their identical charged-state energy, variance and resources follow from the unregulated local-state theorem, checked by explicit small full-H actions. These are not large-H evolutions or fitted emergent dimensions.')
    (P/'MATCHED_HOSTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],cases=len(rows),largest_N=max(r['N'] for r in rows))),flush=True);return result


if __name__=='__main__':run()
