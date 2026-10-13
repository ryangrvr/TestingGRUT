"""Exact coherent sparse-pool destination and record-visibility counts.

The supplied ring is a reference pool. Its Q=N-1 coherent graphs are paths.
No physical dynamics or formation inference is replaced by these counts.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import itertools
import json

P=Path(__file__).resolve().parent


def case(n,R,locked=False):
    es=list(itertools.combinations(range(n),2));M=len(es)
    ids={e:i for i,e in enumerate(es)}
    pool={ids[tuple(sorted((i,(i+1)%n)))] for i in range(n)}
    g=sum(1 << e for e in pool)
    lock=min(pool) if locked else None
    cuts=sorted(pool-({lock} if locked else set()))
    graphs=[g ^ (1 << e) for e in cuts]
    lower=Counter()
    for h in graphs:
        for e in pool:
            if h >> e & 1 and e!=lock:
                lower[h ^ (1 << e)]+=1
    Q=n-1;r=Q-int(locked)
    coherent=F(sum(a*a for a in lower.values()),len(graphs))
    assert coherent==2*r
    common_absent=set(range(M))-pool
    local={e for e in common_absent if min(abs(es[e][0]-es[e][1]),n-abs(es[e][0]-es[e][1]))<=R}
    assert len(local)==n*(R-1)
    fraction=F(len(local),len(common_absent))
    assert fraction==F(2*(R-1),n-3)
    visibility=[]
    for eta in [F(0),F(1,2),F(1)]:
        effective=(1-eta)*r+eta*coherent
        per_edge=effective/F(100*M*M)
        visibility.append(dict(record_overlap=str(eta),effective_removal_norm_square=str(effective),
                               per_common_absent_creation_coefficient=str(per_edge),
                               fraction_with_destination_within_R=str(fraction)))
    return dict(N=n,Q=Q,R=R,locked_connections=int(locked),initial_coherent_graphs=len(graphs),
                deletion_cores=len(lower),common_absent=len(common_absent),local_common_absent=len(local),
                local_fraction=str(fraction),visibility=visibility)


def main():
    out=dict(status='CP12_COHERENT_RING_RECORD_COUNT_PASS',
             cases=[case(n,R,lock) for n,R in [(7,2),(9,3),(12,3),(25,3),(101,3)]
                    for lock in [False,True]],
             scope='Exact first-step amplitude counts in supplied coherent spanning-path pools; no large-time or emergent metric claim.')
    (P/'ring_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],local_fractions=[{'N':x['N'],'fraction':x['local_fraction']}
                       for x in out['cases'] if not x['locked_connections']]),indent=2))


if __name__=='__main__':main()
