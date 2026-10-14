"""Independent exact unregulated full-H moments, N=6,Q=12.

For cubic graph moments only H psi and the projection of H^2 psi onto the
one-action support are required. This projection is an exact contraction,
not a truncation of time evolution or of allowed electric fields.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb,lcm
from pathlib import Path
import hashlib,json,os,sys,time
P=Path(__file__).resolve().parent;S=P.parent/'cp12_ground_connectivity'
assert hashlib.sha256((S/'exact_destination.py').read_bytes()).hexdigest()=='1a1ccb4f06602179f14c2834dea3566fada022ca6a23f58370abe0d4f1e5b35b'
sys.path.insert(0,str(S));from exact_destination import transitions,sea,gauss
sys.path.insert(0,str(P));from topology import A6,B6,coefficients

ZERO=(0,0)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def times(c,x):return c*x[0],c*x[1]
def mul(x,y):return x[0]*y[0]-15*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def conj(x):return x[0],-x[1]
def power(x,k):
    out=(1,0)
    for _ in range(k):out=mul(out,x)
    return out
def dot(v,w,O=lambda s:1):
    out=ZERO
    for s,x in v.items():out=add(out,times(O(s),mul(conj(x),w.get(s,ZERO))))
    return out
def action(v,scale,hopping=True,allowed=None):
    out={}
    for s,x in v.items():
        for dest,c in transitions(6,s,scale,hopping=hopping):
            if allowed is not None and dest not in allowed:continue
            out[dest]=add(out.get(dest,ZERO),times(c,x))
    return {s:x for s,x in out.items() if x!=ZERO}


def graph_stats(g,edges):
    A=[[0]*6 for _ in range(6)]
    for e,(i,j) in enumerate(edges):
        if g>>e&1:A[i][j]=A[j][i]=1
    degrees=list(map(sum,A));tri=sum(A[i][j]*A[j][k]*A[i][k] for i,j,k in combinations(range(6),3))
    sq=0
    for i,j,k,l in combinations(range(6),4):
        sq+=A[i][j]*A[j][k]*A[k][l]*A[l][i]
        sq+=A[i][j]*A[j][l]*A[l][k]*A[k][i]
        sq+=A[i][k]*A[k][j]*A[j][l]*A[l][i]
    return dict(triangles=tri,squares=sq,wedges=sum(comb(d,2) for d in degrees))


def run(mask,name,stats,graphs,edges):
    start=time.monotonic();Q=12;M=15;D=comb(M,Q);den=8**7;norm=D*den**2;scale=lcm(8,150)
    sm=sum(1<<edges.index(e) for e in mask)
    v0={sea(6,g):times(8**(7-(g&sm).bit_count()),power((-7,1),(g&sm).bit_count())) for g in graphs}
    assert dot(v0,v0)==(norm,0)
    v1=action(v0,scale);support=set(v1)|set(v0);v2=action(v1,scale,allowed=support)
    reverse={s:conj(x) for s,x in v0.items()}
    assert action(reverse,scale)=={s:conj(x) for s,x in v1.items()}
    for row in [v0,v1,v2]:
        for s in row:gauss(6,s)
    mean=F(dot(v0,v1)[0],norm*scale);var=F(dot(v1,v1)[0],norm*scale**2)-mean**2
    co={}
    for obs in ['triangles','squares','wedges']:
        O=lambda s:stats[s[0]][obs]
        ov0={s:times(O(s),x) for s,x in v0.items() if O(s)}
        u1=action(ov0,scale)
        assert set(u1)<=support
        c0=F(dot(v0,v0,O)[0],norm)
        c1=F(2*dot(v0,v1,O)[1],norm*scale)
        c2=F(dot(v1,v1,O)[0]-dot(u1,v1)[0],norm*scale**2)
        c3=F(3*dot(v1,v2,O)[1]-dot(u1,v2)[1],3*norm*scale**3)
        # c1,c3 are coefficients divided by sqrt(15).
        co[obs]=[str(c0),str(c1),str(c2),str(c3)]
    # Independently enumerate the conditional motif counts in the complete graph.
    colors={x:[0]*(r+1) for x,r in [('wedges',2),('triangles',3),('squares',4)]}
    for center in range(6):
        for i,j in combinations([x for x in range(6) if x!=center],2):
            colored=sum(tuple(sorted(e)) in mask for e in [(center,i),(center,j)])
            colors['wedges'][colored]+=1
    for i,j,k in combinations(range(6),3):
        colored=sum(tuple(sorted(e)) in mask for e in [(i,j),(j,k),(i,k)])
        colors['triangles'][colored]+=1
    for i,j,k,l in combinations(range(6),4):
        for path in [(i,j,k,l),(i,j,l,k),(i,k,j,l)]:
            colored=sum(tuple(sorted(e)) in mask for e in zip(path,path[1:]+path[:1]))
            colors['squares'][colored]+=1
    for obs,counts in colors.items():
        coef,c,ss=coefficients(M,Q,7,counts)
        full=list(map(F,co[obs]));full[1]*=8;full[3]*=8
        assert full==coef,(name,obs,full,coef)
    # Non-invariant control: full cubic coefficient differs from neutral-only
    # by -(2/3)t_hop^2 times the linear coefficient. This tests the correction
    # which vanishes for every invariant observable in the main comparison.
    O=lambda s:(s[0]&sm).bit_count()
    u1=action({s:times(O(s),x) for s,x in v0.items() if O(s)},scale)
    full1=F(2*dot(v0,v1,O)[1],norm*scale)
    full3=F(3*dot(v1,v2,O)[1]-dot(u1,v2)[1],3*norm*scale**3)
    r1=action(v0,scale,hopping=False)
    r2=action(r1,scale,hopping=False)
    ru1=action({s:times(O(s),x) for s,x in v0.items() if O(s)},scale,hopping=False)
    neutral1=F(2*dot(v0,r1,O)[1],norm*scale)
    neutral3=F(3*dot(r1,r2,O)[1]-dot(ru1,r2)[1],3*norm*scale**3)
    assert full1==neutral1 and full1!=0
    assert full3-neutral3==-F(2,3)*full1
    out=dict(mask=name,mean_energy=str(mean),variance=str(var),graph_count=D,
             one_action_states=len(v1),projected_two_action_states=len(v2),
             full_coefficients=co,coefficient_convention='Odd-order entries multiply sqrt(15); even entries are rational.',
             noninvariant_cubic_control=dict(linear=str(full1),full_cubic=str(full3),neutral_cubic=str(neutral3),
                                             full_minus_neutral=str(full3-neutral3),expected=str(-F(2,3)*full1)),
             gauss_verified=True,seconds=time.monotonic()-start)
    out['phase_conjugation']=dict(real_H_action_verified=True,
        reversed_cubic_over_sqrt15={obs:str(-F(x[3])) for obs,x in co.items()},
        energy_spectral_measure='Identical by K H K = H; this is a theorem, not a numerical spectral scan.')
    print(json.dumps(out,indent=2),flush=True)
    return out


def main():
    edges=list(combinations(range(6),2));graphs=[sum(1<<e for e in c) for c in combinations(range(15),12)]
    stats={g:graph_stats(g,edges) for g in graphs}
    rows=[run(mask,name,stats,graphs,edges) for mask,name in [(A6,'square_rich'),(B6,'triangle_rich')]]
    assert rows[0]['mean_energy']==rows[1]['mean_energy']
    assert rows[0]['variance']==rows[1]['variance']
    for obs in rows[0]['full_coefficients']:
        assert rows[0]['full_coefficients'][obs][:3]==rows[1]['full_coefficients'][obs][:3]
    assert rows[0]['full_coefficients']['triangles'][3]!=rows[1]['full_coefficients']['triangles'][3]
    assert rows[0]['full_coefficients']['squares'][3]!=rows[1]['full_coefficients']['squares'][3]
    out=dict(status='CP12_TOPOLOGY_FULL_H_UNREGULATED_CUBIC_PASS',N=6,Q=12,rows=rows,
             scope='Auxiliary Q=2N independent full-H exact moments. Not a full time evolution or main Q=3N phase result.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'full_moments_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':main()
