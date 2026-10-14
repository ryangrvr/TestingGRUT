"""Unregulated CP-12 destination Taylor controls, exact rational arithmetic.

Independent reachable-state builder: no frozen_model import or electric cutoff.
The physical graph source is supplied; these are not geometry-selection tests.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import math

P = Path(__file__).resolve().parent


def edges(n):
    return list(itertools.combinations(range(n), 2))


def mask(n, pairs):
    es = edges(n)
    return sum(1 << es.index(tuple(sorted(e))) for e in pairs)


def sea(n, g):
    M = len(edges(n))
    return g, sum(1 << (n + e) for e in range(M)), (0,) * (2 * M)


def gauss(n, state):
    g, bits, E = state
    es = edges(n)
    gv = [-((bits >> i) & 1) for i in range(n)]
    assert bits.bit_count() == len(es)
    for e, (i, j) in enumerate(es):
        a, b = E[2*e:2*e+2]
        gv[i] += a
        gv[j] += b
        assert -a-b-((bits >> (n+e)) & 1)+1 == 0
        if not (g >> e & 1):
            assert a == b == 0 and bits >> (n+e) & 1
    assert not any(gv)


def hop(bits, dst, src):
    if not (bits >> src & 1) or bits >> dst & 1:
        return None
    new = bits ^ (1 << src)
    parity = ((bits & ((1 << src)-1)).bit_count()
              + (new & ((1 << dst)-1)).bit_count())
    return new | (1 << dst), (-1 if parity % 2 else 1)


def transitions(n, state, scale, hopping=True):
    es = edges(n)
    M = len(es)
    g, bits, E = state
    Q = g.bit_count()
    occupied = [e for e in range(M) if g >> e & 1]
    removable = []
    for e in occupied:
        if E[2*e] == E[2*e+1] == 0 and bits >> (n+e) & 1:
            removable.append(e)
        if hopping:
            i, j = es[e]
            for leg, vertex, orient in [(0, i, 1), (1, j, -1)]:
                for dst, src, shift in [(vertex, n+e, 1), (n+e, vertex, -1)]:
                    z = hop(bits, dst, src)
                    if z is None:
                        continue
                    newbits, sign = z
                    field = list(E)
                    field[2*e+leg] += shift
                    yield (g, newbits, tuple(field)), scale*orient*sign
    for e in removable:
        for f in range(M):
            if not (g >> f & 1):
                yield (g ^ (1 << e) ^ (1 << f), bits, E), -scale//(10*M)
    diagonal = (scale*Q//2 + 2*scale*(bits & ((1 << n)-1)).bit_count()
                + scale//8*sum(x*x for x in E)
                + scale//(10*M)*len(removable)*(M-Q))
    yield state, diagonal


def apply(n, v, scale, hopping=True):
    out = defaultdict(int)
    for s, a in v.items():
        for target, c in transitions(n, s, scale, hopping):
            out[target] += c*a
    return {s: a for s, a in out.items() if a}


def powers(n, v, order, hopping=True):
    scale = math.lcm(8, 10*len(edges(n)))
    out = [v]
    for _ in range(order):
        out.append(apply(n, out[-1], scale, hopping))
    for row in out:
        for s in row:
            gauss(n, s)
    return scale, out


def coefficient(vs, scale, edge, order, norm):
    total = F(0)
    for a in range(order+1):
        # i**a * (-i)**(order-a) is an exact fourth root of unity.
        phase = (2*a-order) % 4
        sign = 1 if phase == 0 else (-1 if phase == 2 else 0)
        if not sign:
            continue
        dot = sum(x*vs[order-a].get(s, 0) for s, x in vs[a].items()
                  if s[0] >> edge & 1)
        total += F(sign*dot, norm*scale**order*math.factorial(a)*math.factorial(order-a))
    return total


def deletion(n, v):
    out = defaultdict(int)
    for (g, bits, E), a in v.items():
        for e in range(len(edges(n))):
            if g >> e & 1 and bits >> (n+e) & 1 and E[2*e] == E[2*e+1] == 0:
                out[g ^ (1 << e), bits, E] += a
    return {s: a for s, a in out.items() if a}


def first_step(n, v, name):
    M = len(edges(n))
    Q = next(iter(v))[0].bit_count()
    assert all(s[0].bit_count() == Q for s in v)
    union = 0
    for s in v:
        gauss(n, s)
        union |= s[0]
    absent = [e for e in range(M) if not (union >> e & 1)]
    assert absent
    norm = sum(x*x for x in v.values())
    lowered = deletion(n, v)
    effective = F(sum(x*x for x in lowered.values()), norm)
    scale, vs = powers(n, v, 1)
    expected = effective/F(100*M*M)
    rows = []
    for e in absent:
        actual = F(sum(x*x for s, x in vs[1].items() if s[0] >> e & 1), norm*scale**2)
        assert actual == expected
        rows.append(dict(edge=edges(n)[e],coefficient=str(actual)))
    return dict(name=name,N=n,Q=Q,initial_states=len(v),effective_removal_norm_square=str(effective),
                common_absent_count=len(absent),rows=rows)


def taylor(n, v, order, name):
    union = 0
    for s in v:
        union |= s[0]
    norm = sum(x*x for x in v.values())
    sc, vs = powers(n, v, order)
    sc0, neutral = powers(n, v, order, hopping=False)
    rows=[]
    for e, pair in enumerate(edges(n)):
        if union >> e & 1:
            continue
        row={str(k):str(coefficient(vs,sc,e,k,norm)) for k in range(2,order+1,2)}
        rows.append(dict(edge=pair,coefficients=row))
    # Without hopping all common-vacant edges are indistinguishable at every tested order.
    for k in range(2,order+1,2):
        values={coefficient(neutral,sc0,edges(n).index(tuple(r['edge'])),k,norm) for r in rows}
        assert len(values)==1
    return dict(name=name,N=n,Q=next(iter(v))[0].bit_count(),norm_square=norm,
                reachable_counts=[len(x) for x in vs],maximum_reachable_electric_field=
                max(abs(e) for x in vs for s in x for e in s[2]),rows=rows,
                no_hopping_destination_equality_through_order=order)


def main():
    # Two histories share two bonds; two destination bonds are absent in both.
    a=mask(4,[(0,1),(1,2),(2,3)])
    b=mask(4,[(0,2),(1,2),(2,3)])
    plus={sea(4,a):1,sea(4,b):1}
    minus={sea(4,a):1,sea(4,b):-1}
    charged={}
    for (g,bits,E),amp in plus.items():
        e=edges(4).index((1,2))
        field=list(E);field[2*e]=1
        charged[g,bits^(1 << (4+e))^(1 << 1),tuple(field)]=amp
    tests=[first_step(4,plus,'two_histories_plus'),first_step(4,minus,'two_histories_minus'),
           first_step(4,charged,'two_histories_shared_charged_bond')]
    # Orthogonal real/imaginary quadratures implement an exact relative phase i.
    vr={sea(4,a):1};vi={sea(4,b):1};sc,xr=powers(4,vr,1);_,xi=powers(4,vi,1)
    dr=deletion(4,vr);di=deletion(4,vi)
    effective=F(sum(z*z for z in dr.values())+sum(z*z for z in di.values()),2)
    values=[]
    for e in range(6):
        if (a|b) >> e & 1:continue
        actual=F(sum(z*z for s,z in xr[1].items() if s[0] >> e & 1)
                 +sum(z*z for s,z in xi[1].items() if s[0] >> e & 1),2*sc*sc)
        assert actual==effective/F(3600)
        values.append(dict(edge=edges(4)[e],coefficient=str(actual)))
    tests.append(dict(name='two_histories_relative_phase_i',N=4,Q=3,initial_states=2,
                      effective_removal_norm_square=str(effective),common_absent_count=2,rows=values))
    # A supplied square with a tail; Q=N-1. Equal coherent vacancy state on its five bonds.
    pool=mask(5,[(0,1),(1,2),(2,3),(3,0),(3,4)])
    v={sea(5,pool^(1 << e)):1 for e in range(10) if pool >> e & 1}
    tests.append(first_step(5,v,'square_tail_coherent_vacancy'))
    coherent=taylor(5,v,10,'square_tail_coherent_vacancy')
    singles=[taylor(5,{s:1},10,'recorded_'+str(i)) for i,s in enumerate(v)]
    # Every single history has additional absent bonds. Compare only common-vacant destinations.
    common=[tuple(r['edge']) for r in coherent['rows']]
    mixture={pair:{str(k):sum((F(next(r for r in data['rows'] if tuple(r['edge'])==pair)
                         ['coefficients'][str(k)]) for data in singles),F(0))/len(singles)
                       for k in range(2,11,2)} for pair in common}
    rows={tuple(r['edge']):r['coefficients'] for r in coherent['rows']}
    difference={str(k):F(rows[(0,4)][str(k)])-F(rows[(1,4)][str(k)]) for k in range(2,11,2)}
    assert all(difference[str(k)]==0 for k in [2,4,6,8])
    assert difference['10']==-F(1103,3240000000)
    small=taylor(4,plus,6,'two_histories_plus')
    srows={tuple(r['edge']):r['coefficients'] for r in small['rows']}
    assert F(srows[(1,3)]['6'])-F(srows[(0,3)]['6'])==F(1,72000)
    out=dict(status='CP12_COHERENT_DESTINATION_EXACT_PASS',first_step_controls=tests,
             two_history_taylor=small,primary_coherent_taylor=coherent,
             primary_dephased_taylor=[dict(edge=e,coefficients={k:str(x) for k,x in row.items()})
                                      for e,row in mixture.items()],
             primary_distance_two_minus_distance_three={k:str(x) for k,x in difference.items()},
             scope='Full unregulated frozen-law Taylor coefficients at fixed finite N; supplied graph pools; no continuum inference.')
    (P/'exact_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],distance_difference=out['primary_distance_two_minus_distance_three'],
                         reachable_counts=coherent['reachable_counts']),indent=2))


if __name__=='__main__':
    main()
