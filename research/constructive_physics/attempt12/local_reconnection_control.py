"""Kinematic compatibility test on all port matchings, not a new Hamiltonian."""
from pathlib import Path
import itertools,json

def matchings(items):
 if not items:yield ();return
 a=items[0]
 for j in range(1,len(items)):
  for tail in matchings(items[1:j]+items[j+1:]):yield tuple(sorted(((a,items[j]),)+tail))
def canonical(p):return tuple(sorted(tuple(sorted(e)) for e in p))
def components(p,n):
 par=list(range(n))
 def root(i):
  while par[i]!=i:i=par[i]
  return i
 for a,b in p:
  x,y=root(a//2),root(b//2)
  if x!=y:par[x]=y
 return tuple(sorted(tuple(i for i in range(n) if root(i)==r) for r in set(root(i) for i in range(n))))
rows=[]
for n in [3,4,5]:
 ms=list(matchings(tuple(range(2*n))));parts={p:components(p,n) for p in ms};counts=dict(all_moves=0,same_component_both_directions=0,partition_changes=0,forbidden_partition_changes=0)
 for p in ms:
  for i,j in itertools.combinations(range(n),2):
   a,b=p[i];c,d=p[j];support={a//2,b//2,c//2,d//2}
   for pair in [((a,c),(b,d)),((a,d),(b,c))]:
    nxt=list(p);nxt[i],nxt[j]=pair;nxt=canonical(nxt)
    if nxt<=p:continue
    counts['all_moves']+=1;old,new=parts[p],parts[nxt];change=old!=new
    local_old=any(support<=set(comp) for comp in old)
    local_new=any(support<=set(comp) for comp in new)
    both=local_old and local_new
    counts['same_component_both_directions']+=both;counts['partition_changes']+=change
    if change:
     assert not both;counts['forbidden_partition_changes']+=1
 rows.append(dict(N=n,matchings=len(ms),**counts))
Path(__file__).with_name('local_reconnection_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2));print('CP12_RG02_LOCALITY_COMPONENT_COMPATIBILITY_PASS')
