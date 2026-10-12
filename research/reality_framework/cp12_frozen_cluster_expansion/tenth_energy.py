"""Evaluate exact inclusion-exclusion combination for edge order ten at t=1.
Each energy is a finite Gauss-law graph block with common host a=.1.
"""
from pathlib import Path
from fractions import Fraction
import sys,json,itertools,gc,hashlib,math
import networkx as nx
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build

def known_energy(g,L):
 for f in P.glob(f'energy_*_L{L}.json'):
  d=json.loads(f.read_text());h=nx.cycle_graph(d['N']) if d['name'].startswith('cycle') else nx.path_graph(d['N'])
  if len(h)==len(g) and h.number_of_edges()==g.number_of_edges() and nx.is_isomorphic(h,g):return d
 return None

def solve(g,L,idx):
 cache=P/f'tenth_energy_{idx}_L{L}.json'
 if cache.exists():
  d=json.loads(cache.read_text());h=nx.Graph();h.add_nodes_from(range(d['N']));h.add_edges_from(d['pairs'])
  assert d['L']==L and nx.is_isomorphic(g,h)
  return d
 known=known_energy(g,L)
 if known is not None:
  known=dict(known,reused=True,index=idx,pairs=list(g.edges()));cache.write_text(json.dumps(known,indent=2)+'\n');return known
 n=len(g);q=g.number_of_edges();edges=list(itertools.combinations(range(n),2));M=len(edges);active=[edges.index(tuple(sorted(e))) for e in g.edges()];mask=sum(1<<e for e in active)
 H,states,lookup,_=build(n,q,L,masks=[mask],rewire=False);del lookup
 rem=np.array([sum((bits>>(n+e)&1) and E[2*e]==E[2*e+1]==0 for e in active) for _,bits,E in states]);boundary=np.array([max(map(abs,s[2]))==L for s in states]);dim=len(states);del states;gc.collect()
 H=H+diags(.1*rem);del rem;gc.collect()
 E,U=eigsh(H,k=1,which='SA',tol=1e-13,v0=np.random.default_rng(20261012+idx).normal(size=dim));v=U[:,0]
 d=dict(index=idx,pairs=list(g.edges()),N=n,Q=q,L=L,a=.1,dimension=dim,energy=float(E[0]),residual=float(np.linalg.norm(H@v-E[0]*v)),cutoff_boundary=float(abs(v)**2@boundary),reused=False)
 assert d['residual']<1e-10;cache.write_text(json.dumps(d,indent=2)+'\n');return d

def main(L):
 expression=json.loads((P/'tenth_energy_expression.json').read_text());rows=[]
 for idx,r in enumerate(expression):
  g=nx.Graph();g.add_nodes_from(range(r['N']));g.add_edges_from(r['pairs']);d=solve(g,L,idx)
  row=dict(**r,energy=d['energy'],residual=d['residual'],cutoff_boundary=d['cutoff_boundary'],dimension=d['dimension']);rows.append(row);print(json.dumps(row),flush=True)
 value=math.fsum(float(Fraction(r['coefficient']))*r['energy'] for r in rows)
 residual=sum(abs(float(Fraction(r['coefficient'])))*r['residual'] for r in rows)
 previous=json.loads((P/f'winding_weight_l9_L{L}.json').read_text())['thin_minus_balanced_first_sensitive_density']
 out=dict(L=L,a=.1,hopping=1,rows=rows,tenth_increment=value,cumulative_difference_through_ten=previous+value,ninth_difference=previous,residual_indicator=residual,scope='Truncated edge-cluster density difference only, without thermodynamic remainder bound.')
 (P/f'tenth_results_L{L}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='rows'}),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--L',type=int,default=3);main(p.parse_args().L)
