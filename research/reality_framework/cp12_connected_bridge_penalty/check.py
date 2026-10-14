#!/usr/bin/env python3
"""Exact CP12 Q=N connected-graph rewiring counts, including Gauss-locked bridges."""
import itertools,json
from collections import deque
from pathlib import Path

GAMMA=.1
def is_connected(n,edges):
    adj=[[] for _ in range(n)]
    for a,b in edges:adj[a].append(b);adj[b].append(a)
    todo=deque([0]);seen={0}
    while todo:
        for v in adj[todo.popleft()]:
            if v not in seen:seen.add(v);todo.append(v)
    return len(seen)==n

def components(n,edges):
    adj=[[] for _ in range(n)]
    for a,b in edges:adj[a].append(b);adj[b].append(a)
    seen=set();parts=[]
    for i in range(n):
        if i in seen:continue
        todo=[i];seen.add(i);part={i}
        while todo:
            for v in adj[todo.pop()]:
                if v not in seen:seen.add(v);part.add(v);todo.append(v)
        parts.append(part)
    return parts

def gauss(n,present,holes,vertex,fields):
    g=[-int(i in vertex) for i in range(n)]
    for a,b in itertools.combinations(range(n),2):
        x,y=fields.get((a,b),(0,0));g[a]+=x;g[b]+=y
        if x+y!=int((a,b) in holes):return False
        if (a,b) not in present and (x or y or (a,b) in holes):return False
    return all(v==0 for v in g)

def check(n,kind,lock_one=False):
    if kind=='cycle':pairs=[(i,(i+1)%n) for i in range(n)]
    elif kind=='triangle_tail':pairs=[(0,1),(1,2),(2,0)]+[(i,i+1) for i in range(2,n-1)]
    elif kind=='square_tail':pairs=[(0,1),(1,2),(2,3),(3,0)]+[(i,i+1) for i in range(3,n-1)]
    else:raise ValueError(kind)
    present={tuple(sorted(e)) for e in pairs}
    allpairs=set(itertools.combinations(range(n),2));M=len(allpairs)
    assert len(present)==n and is_connected(n,present)
    absent=allpairs-present
    bridge_data={}
    for e in present:
        parts=components(n,present-{e})
        if len(parts)==2:bridge_data[e]=len(parts[0])
    lock=min(bridge_data) if lock_one and bridge_data else None
    holes=set();vertex=set();fields={}
    if lock:
        a,b=lock;holes.add(lock);vertex.add(a);fields[lock]=(1,0)
    assert gauss(n,present,holes,vertex,fields)
    removable=present-({lock} if lock else set())
    expected=sum(M-n+1-s*(n-s) for e,s in bridge_data.items() if e in removable)
    actual=sum(not is_connected(n,(present-{e})|{f}) for e in removable for f in absent)
    assert expected==actual
    lower=max(0,M-n+1-n*n//4)*(len(bridge_data)-int(bool(lock)))
    assert actual>=lower
    degree=[0]*n
    for a,b in present:degree[a]+=1;degree[b]+=1
    D2=sum(abs(d-2) for d in degree)
    assert D2<=2*len(bridge_data)
    return dict(N=n,kind=kind,locked_bridge=bool(lock),bridges=len(bridge_data),
       degree_defect=D2,exit_swaps=actual,exit_lower=lower,
       short_time_T2_probability_coefficient=GAMMA**2*actual/M**2)

def main():
    data=[]
    for n in range(5,11):
        for kind in ('cycle','triangle_tail','square_tail'):
            data.append(check(n,kind))
            if kind!='cycle':data.append(check(n,kind,True))
    out=dict(status='CP12_CONNECTED_BRIDGE_BOUNDARY_EXACT_PASS',
             scope='exact Q=N one-swap connectivity, Gauss-valid bridge locks; not ground-state connectivity',
             data=data)
    Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
    for x in data:
        if x['N'] in (6,9) and not x['locked_bridge']:print(x)
if __name__=='__main__':main()
