"""Exact leaf migration in a NEW two-child branching grammar, not S2 path law."""
from collections import defaultdict


def neighbors(n,state):
    edges,q,p=state
    children=[0]*n;bound={0}
    for u,v in edges:children[u]+=1;bound.add(v)
    free=set(range(n))-bound
    # Local counter at each bound parent has values 0,1,2. Removal tests
    # the child's local zero-child counter; these are finite-support terms.
    for u in bound:
        if children[u]<2:
            for v in free:
                assert not(q>>v&1)
                yield (edges|{(u,v)},q|(1<<v),p),-1
    for u,v in edges:
        if children[v]==0:
            assert q>>v&1
            yield (edges-{(u,v)},q^(1<<v),p),-1
        if p==u:yield (edges,q,v),-1
        elif p==v:yield (edges,q,u),-1


def apply(n,vector):
    out=defaultdict(int)
    for c,a in vector.items():
        for dest,v in neighbors(n,c):out[dest]+=a*v
    return {c:a for c,a in out.items() if a}


def arms(length):
    n=2*length+1;edges=set()
    for arm in [list(range(1,length+1)),list(range(length+1,2*length+1))]:
        for u,v in zip([0]+arm,arm):edges.add((u,v))
    return n,frozenset(edges),length,2*length


def check(length):
    n,edges,i,j=arms(length);q=(1<<n)-1
    source=(edges,q,i);target=(edges,q,j)
    a=[{source:1}];b=[{target:1}]
    for _ in range(2):a.append(apply(n,a[-1]));b.append(apply(n,b[-1]))
    dot=lambda x,y:sum(v*y.get(c,0) for c,v in x.items())
    lower=[dot(a[0],b[0]),dot(a[1],b[0]),dot(a[2],b[0]),dot(a[2],b[1]),dot(a[2],b[2])]
    assert lower==[0]*5
    fifth=sum(amp*v*b[2].get(dest,0) for c,amp in a[2].items() for dest,v in neighbors(n,c))
    assert fifth==-2
    parent=next(u for u,v in edges if v==j)
    chain=[source,
        (edges-{(parent,j)},q^(1<<j),i),
        ((edges-{(parent,j)})|{(i,j)},q,i),
        ((edges-{(parent,j)})|{(i,j)},q,j),
        (edges-{(parent,j)},q^(1<<j),j),target]
    for old,new in zip(chain,chain[1:]):assert (new,-1) in list(neighbors(n,old))
    for ee,qq,pp in chain:
        bound={0}|{v for u,v in ee};assert qq==sum(1<<v for v in bound)
    return dict(N=n,initial_tree_distance=2*length,
        protected_primitive_word_length=5,full_fifth_power_unit_coefficient=fifth,
        normalized_coefficient='-2*g^4*J/N^2 if each assembly move is g/sqrt(N)',
        lower_powers_zero_through=4,source_H2_states=len(a[2]),target_H2_states=len(b[2]),
        original_tree_and_storage_fully_restored=True,
        all_five_transitions_stay_fuel_consistent=True,
        meaning='Protected operational distance is <=5; the old growing tree distance cannot be retained as that metric.')


def run():
    rows=[check(length) for length in [3,6,12,24]]
    binary=[dict(depth=h,N=2**(h+1)-1,leaves=2**h,
        leaf_operational_ball_radius=5,leaf_ball_volume_lower=2**h)
        for h in [3,4,5,6,7,8]]
    return dict(status='PASS',law='Separate finite-support two-child/leaf-removal proposal; not inserted into single-head S2.',
        rows=rows,binary_tree_ball_obstruction=binary,
        conclusion='The consistent branching class contains bounded-radius operational balls with extensive leaf count. This rejects a uniform finite-dimensional ball-growth claim for this naive class, not every branching architecture.',
        limits='Exact primitive words and full fifth coefficients. Not a finite-time large-N remainder bound or an original-neutral tree-formation result.')
