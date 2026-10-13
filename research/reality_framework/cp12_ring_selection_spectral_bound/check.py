#!/usr/bin/env python3
"""CP12 ring-sector cut: exact full-graph enumeration, no GRUT imports.

Counts ALL spanning simple rings and ALL graphs accessible by a single
occupied/empty edge swap. The graph-only adjacency norm bounds the CP12
graph-changing block including fermion, U(1), and Gauss neutrality.
Run: python check.py (stdlib only).
"""
import itertools,collections,math,json
from pathlib import Path
GAMMA=.1

def test(N):
    pairs=list(itertools.combinations(range(N),2))
    edgeid={e:i for i,e in enumerate(pairs)}
    M=len(pairs);Q=N;R=set()
    for perm in itertools.permutations(range(1,N)):
        v=(0,)+perm+(0,)
        R.add(frozenset(edgeid[tuple(sorted(edge))] for edge in
                        zip(v,v[1:])))
    assert len(R)==math.factorial(N-1)//2
    count=collections.Counter()
    for G in R:
        for old in G:
            for new in range(M):
                if new in G:continue
                H=frozenset((G-{old})|{new})
                assert H not in R
                count[H]+=1
    hist=collections.Counter(count.values())
    assert set(hist).issubset({1,2}),hist
    assert sum(k*v for k,v in hist.items())==len(R)*N*(M-N)
    assert hist[2]==len(R)*N*(N-3)
    scalar_sq=sum(k*k*v for k,v in hist.items())/len(R)
    assert scalar_sq==N*(M+N-6)
    b=GAMMA/M*math.sqrt(scalar_sq)
    assert b <= GAMMA/M*math.sqrt(2*N*(M-N))+1e-14
    return dict(N=N,M=M,spanning_rings=len(R),
        distinct_external_graphs=len(count),
        incoming_preimage_distribution=dict(hist),
        exact_scalar_rewire_cut_norm=b,
        universal_physical_cut_upper_bound=b)

def main():
    rows=[test(N) for N in range(4,9)]
    d=(-.635203537847)-(-.637262028038)
    b=rows[2]["universal_physical_cut_upper_bound"]
    upper=b*b/(b*b+d*d)
    assert abs(upper-.9989417755586539)<1e-11
    assert all(x["exact_scalar_rewire_cut_norm"]>0 for x in rows)
    result=dict(model="Frozen CP12 Q=N neutral bond swap",
        exact_scalar_cut="Gamma/M * sqrt[N(M+N-6)]",
        full_physical_cut_operator_norm_bound="<= scalar cut",
        conditional_ground_state_theorem="If E_R-E_comp >= delta>0 then P_R <= b^2/(b^2+delta^2)",
        regulated_N6=dict(delta=d,ring_probability_upper_bound=upper),
        status="N6 bound valid but too weak to reject ring dominance",
        rows=rows,
        warnings=[
          "CP12 matter-block N6 energies are taken from published 7be83fc; not recalculated",
          "No large-N positive gap lower bound established",
          "Does not simulate full many-body Hilbert space",
          "Proof applies to projected ring class, not arbitrary spatial phases"])
    Path(__file__).with_name("results.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    print("CP12_RING_SECTOR_COUPLING_CHECK_PASS")

if __name__=="__main__":main()
