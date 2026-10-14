"""Unregulated full-Hamiltonian one-action controls on Q=2N and Q=3N hosts.
All rewiring destinations are counted; no fixed-graph diagonal is dropped.
The prepared hosts and chosen loop are supplied benchmarks, not emergent space.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,lcm
import hashlib,json,sys,time,os

P=Path(__file__).resolve().parent;S=P.parent/'cp12_ground_connectivity'
sys.path.insert(0,str(S));from exact_destination import transitions,sea,gauss
sys.path.insert(0,str(P));from check import torus

def shift(n,state,cycle,sign):
    from itertools import combinations
    all_edges=list(combinations(range(n),2));ids={e:i for i,e in enumerate(all_edges)}
    g,b,E=state;E=list(E)
    for i,j in zip(cycle,cycle[1:]+cycle[:1]):
        e=ids[tuple(sorted((i,j)))];assert g>>e&1
        z=sign*(1 if i<j else -1);E[2*e]+=z;E[2*e+1]-=z
    s=(g,b,tuple(E));gauss(n,s);return s

def moment(n,s):
    scale=lcm(8,10*comb(n,2));diagonal=0;off_square=0;hops=0;rewires=0
    for target,c in transitions(n,s,scale):
        if target==s:diagonal+=c
        else:
            off_square+=c*c
            if target[0]!=s[0]:rewires+=1
            else:hops+=1
    return dict(mean=F(diagonal,scale),variance=F(off_square,scale**2),hopping_outputs=hops,rewiring_outputs=rewires)

def main():
    from itertools import combinations
    out=[]
    for dims in [(6,6),(5,5,5)]:
        start=time.monotonic();n,edges,faces=torus(dims);q=len(edges);M=comb(n,2)
        ids={e:i for i,e in enumerate(combinations(range(n),2))};g=sum(1<<ids[e] for e in edges)
        s=sea(n,g);gauss(n,s);cycle=list(faces[0]);ell=len(cycle)
        rows=[]
        for sign in [0,1,-1]:
            a=moment(n,s if not sign else shift(n,s,cycle,sign));removable=q-(ell if sign else 0)
            expected=F(q,2)+F(ell,4)*bool(sign)+F(removable*(M-q),10*M)
            assert a['mean']==expected and a['hopping_outputs']==2*q
            assert a['rewiring_outputs']==removable*(M-q)
            assert a['variance']==2*q+F(removable*(M-q),100*M*M)
            rows.append(dict(sign=sign,**{k:str(v) if isinstance(v,F) else v for k,v in a.items()}))
        slope=F(rows[1]['mean'])+F(rows[2]['mean'])-2*F(rows[0]['mean'])
        expected_slope=2*ell*(F(1,4)-F(M-q,10*M));assert slope==expected_slope
        row=dict(shape=list(dims),N=n,Q=q,cycle_length=ell,cycle_rank=q-n+1,
                 normalized_Wilson_retarded_slope=str(slope),slope_float=float(slope),
                 physical_one_action_records=rows,seconds=time.monotonic()-start)
        out.append(row);print(json.dumps(row),flush=True)
    result=dict(status='CP12_HIGHER_DENSITY_WILSON_MOMENTS_PASS',hosts=out,
                scope='Exact rational full-law moments on prepared neutral graphs; not a ground state, wave dispersion, or geometry-selection calculation.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'wilson_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
if __name__=='__main__':main()
