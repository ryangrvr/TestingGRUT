"""Exact graph-choice accessibility via signed cycle/path covers."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,permutations
from math import comb,factorial
import argparse,json,time
from cycle_cover import cover_coefficient,euler_orientations,target_state,direct_amplitude

P=Path(__file__).resolve().parent


def automorphisms(n,es):
    edge_set=set(es)
    return [p for p in permutations(range(n))
            if {tuple(sorted((p[a],p[b]))) for a,b in es}==edge_set]


def marked_orbits(n,es):
    orientations=list(euler_orientations(n,es));autos=automorphisms(n,es)
    maps=[]
    for p in autos:
        edge_map=[]
        for a,b in es:edge_map.append((es.index(tuple(sorted((p[a],p[b])))),p[a]<p[b]))
        maps.append((p,edge_map))
    groups={}
    for arcs in orientations:
        pattern=sum(1<<e for e,(a,b) in enumerate(arcs) if a<b)
        transformed=[]
        for p,edge_map in maps:
            pp=sum(1<<new for e,(new,same) in enumerate(edge_map) if bool(pattern>>e&1)==same)
            transformed.append((p,pp))
        for marked in combinations(range(n),2):
            key=min((pp,sum(1<<p[v] for v in marked)) for p,pp in transformed)
            groups[key]=groups.get(key,0)+1
    assert sum(groups.values())==len(orientations)*comb(n,2)
    return orientations,autos,groups


def one_host(n,es,name):
    start=time.monotonic();orientations,autos,groups=marked_orbits(n,es)
    rows=[];total=F(0)
    for (pattern,marked_mask),multiplicity in sorted(groups.items()):
        arcs=[(a,b) if pattern>>e&1 else (b,a) for e,(a,b) in enumerate(es)]
        marked=tuple(v for v in range(n) if marked_mask>>v&1)
        coefficient=cover_coefficient(n,arcs,marked)
        # Each distinct marked vertex has two outgoing choices of its hole.
        # Their magnitudes are equal; the Slater row-order sign can differ.
        physical_count=4*multiplicity
        total+=physical_count*coefficient*coefficient
        rows.append(dict(orientation_pattern=pattern,marked_vertices=list(marked),
            orientation_marked_orbit_size=multiplicity,physical_hole_configurations=physical_count,
            amplitude_taylor_coefficient=str(coefficient)))
    labelled_hosts=factorial(n)//len(autos)
    neutral_graphs=comb(comb(n,2),2*n)
    initial_coefficient=F(labelled_hosts,neutral_graphs)*total
    triangles=sum(all(tuple(sorted(e)) in set(es) for e in combinations(c,2))
                  for c in combinations(range(n),3))
    out=dict(N=n,Q=2*n,host=name,edges=[list(e) for e in es],
        automorphism_group_order=len(autos),labelled_host_count=labelled_hosts,
        neutral_graph_count=neutral_graphs,euler_orientations=len(orientations),
        marked_orbits=len(rows),minimum_target_configurations_per_host=4*len(orientations)*comb(n,2),
        triangles=triangles,minimum_hops=4*n+2,probability_leading_order=8*n+4,
        electric_square_minimum=4*n+6,
        conditional_on_this_host_coefficient=str(total),conditional_on_this_host_coefficient_float=float(total),
        actual_neutral_ray_joint_coefficient=str(initial_coefficient),actual_neutral_ray_joint_coefficient_float=float(initial_coefficient),
        rows=rows,seconds=time.monotonic()-start)
    print(json.dumps({k:v for k,v in out.items() if k not in ['rows','edges']}),flush=True)
    return out


def checks(hosts):
    results=[]
    # Independent integer half-leg DP checks extremal nonzero and first
    # orbit coefficients, including the largest graph-choice hosts.
    for host in hosts:
        n=host['N'];es=[tuple(e) for e in host['edges']]
        candidates={0,len(host['rows'])-1,
                    min(range(len(host['rows'])),key=lambda j:abs(F(host['rows'][j]['amplitude_taylor_coefficient']))),
                    max(range(len(host['rows'])),key=lambda j:abs(F(host['rows'][j]['amplitude_taylor_coefficient'])))}
        for j in sorted(candidates):
            row=host['rows'][j];pattern=row['orientation_pattern']
            arcs=[(a,b) if pattern>>e&1 else (b,a) for e,(a,b) in enumerate(es)]
            marked=tuple(row['marked_vertices']);state=target_state(n,es,arcs,marked)
            actual,layer=direct_amplitude(n,state);order=4*n+2
            expected=F(row['amplitude_taylor_coefficient'])
            assert abs(F(actual,factorial(order)))==abs(expected),(host['host'],j)
            results.append(dict(host=host['host'],orbit_index=j,N=n,integer_amplitude=str(actual),
                cover_coefficient=str(expected),maximum_DP_layer=layer))
    # Check a simple one-in/one-out loop, independent of 4-regular graphs.
    for n in [3,4]:
        es=sorted(tuple(sorted((j,(j+1)%n))) for j in range(n))
        arcs=[(a,b) if (b-a)%n==1 else (b,a) for a,b in es]
        c=cover_coefficient(n,arcs,())
        s=target_state(n,es,arcs,())
        a,layer=direct_amplitude(n,s)
        assert abs(F(a,factorial(2*n)))==abs(c)
        results.append(dict(host=f'cycle_{n}',N=n,integer_amplitude=str(a),cover_coefficient=str(c),maximum_DP_layer=layer))
    return results


def run(output=None,include_seven=True):
    start=time.monotonic();data=[]
    cases=[(5,list(combinations(range(5),2)),'K5'),
           (6,[e for e in combinations(range(6),2) if e not in {(0,1),(2,3),(4,5)}],'octahedron')]
    if include_seven:
        cycle={tuple(sorted((j,(j+1)%7))) for j in range(7)}
        split={(0,1),(1,2),(0,2),(3,4),(4,5),(5,6),(3,6)}
        cases.extend([(7,[e for e in combinations(range(7),2) if e not in cycle],'complement_C7'),
                      (7,[e for e in combinations(range(7),2) if e not in split],'complement_C3_plus_C4')])
    for n,es,name in cases:data.append(one_host(n,es,name))
    assert data[0]['actual_neutral_ray_joint_coefficient']=='1326415312/848779513441154588027109375'
    checked=checks(data)
    seven=[x for x in data if x['N']==7]
    if seven:
        # Every 4-regular simple graph on seven vertices has a 2-regular
        # complement, hence C7 or C3 disjoint C4. Thus this is exhaustive.
        total=sum(F(x['actual_neutral_ray_joint_coefficient']) for x in seven)
        bias=dict(full_target_leading_coefficient=str(total),
            conditional_C7_complement_fraction=str(F(seven[0]['actual_neutral_ray_joint_coefficient'])/total),
            prior_C7_complement_fraction=str(F(seven[0]['labelled_host_count'],sum(x['labelled_host_count'] for x in seven))),
            rate_ratio_per_host=str(F(seven[0]['conditional_on_this_host_coefficient'])/F(seven[1]['conditional_on_this_host_coefficient'])))
    else:bias=None
    out=dict(status='EXACT_GRAPH_CHOICE_ACCESS_PASS',hosts=data,independent_DP_checks=checked,
        N7_joint_kernel=bias,seconds=time.monotonic()-start,
        scope='Exact first-permitted-order probabilities from the original uniform coherent neutral graph ray. No diagonal/rewiring term contributes at this order. Final charged-sector readout is not an initial mask. Conditional rare-sector host weights do not establish dominant geometry, long-time persistence, dimension, or a causal cone.')
    dest=Path(output) if output else P/'RESULTS.json';dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],seconds=out['seconds'],N7_joint_kernel=bias)),flush=True)
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output');ap.add_argument('--skip-seven',action='store_true')
    a=ap.parse_args();run(a.output,not a.skip_seven)
