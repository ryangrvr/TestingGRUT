"""Full frozen CP12 propagation for an auxiliary N=4,Q=4 sector.
Large-N/high-density results live in preparation.py and are initial derivatives,
not finite-time extrapolations of this small-system calculation.
"""
from pathlib import Path
from itertools import combinations
from math import comb,factorial
import hashlib,json,os,sys,time
import numpy as np
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
S=P.parent/'cp12_ground_connectivity'
assert hashlib.sha256((S/'frozen_model.py').read_bytes()).hexdigest()=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
sys.path.insert(0,str(S));from frozen_model import build
sys.path.insert(0,str(P));from preparation import pattern,quotient,motif_colors,conditional


def coefficients(H,V,O,order=6):
    powers=[V]
    for j in range(1,order+1):powers.append((-1j)*(H@powers[-1])/j)
    return [float(sum(np.sum(powers[a].conj()*O[:,None]*powers[j-a]) for a in range(j+1)).real) for j in range(order+1)]


def run(L):
    start=time.monotonic();n=4;Q=4;M=6;q=quotient(n,Q)
    H,states,lookup,_=build(n,Q,L);dim=len(states);es=list(combinations(range(n),2))
    bits=sum(1<<(n+e) for e in range(M));seas=[i for i,(g,f,E) in enumerate(states) if f==bits and not any(E)]
    D=comb(M,Q);assert len(seas)==D
    masks=[states[i][0] for i in seas];W=np.zeros((dim,D),complex);W[seas,np.arange(D)]=1/np.sqrt(D)
    graphstats={}
    for g in masks:
        G=np.zeros((n,n),int)
        for e,(i,j) in enumerate(es):
            if g>>e&1:G[i,j]=G[j,i]=1
        d=G.sum(axis=1);wedges=int(sum(d*(d-1)//2));tri=int(np.trace(G@G@G)//6)
        square=int((np.trace(G@G@G@G)-2*Q-4*wedges)//8)
        graphstats[g]=[wedges,tri,square]
    diagstats=np.array([graphstats[g]+[(f&15).bit_count(),sum(x*x for x in E)/8,int(any(abs(x)==L for x in E))] for g,f,E in states],float).T
    names=['wedges','triangles','squares','vertex_fermions','electric_energy','cutoff_boundary_probability']
    coeffs={};ps=[]
    for kind in ['random','concentrated','spread']:
        A=pattern(n,kind);ks=np.array([sum(A[i,j] for e,(i,j) in enumerate(es) if g>>e&1) for g in masks])
        phase=np.exp(1j*q['theta']*ks);v=W@phase;ps.append(v)
        expect=float(np.vdot(v,H@v).real);variance=float(np.vdot(H@v,H@v).real-expect**2)
        assert abs(expect-(Q/2+q['rewiring_energy']))<1e-12
        assert abs(variance-(2*Q+q['pure_rewiring_variance']))<2e-12
        stats={}
        Hr=H[seas,:][:,seas];vr=phase[:,None]/np.sqrt(D)
        for j,name in enumerate(names[:3]):
            full=coefficients(H,v[:,None],diagstats[j]);neutral=coefficients(Hr,vr,diagstats[j,seas])
            exact=conditional(q,motif_colors(A)[name])
            assert abs(full[0]-q['w']@exact)<1e-12
            assert abs(full[1]-q['p1']@exact)<1e-12
            assert abs(full[2]-q['p2']@exact)<2e-12
            assert max(abs(full[z]-neutral[z]) for z in range(3))<2e-12
            stats[name]=dict(full_T_coefficients=full,neutral_only_T_coefficients=neutral)
        coeffs[kind]=dict(energy=expect,variance=variance,observables=stats)
    mixture_energy=float(np.sum(W.conj()*(H@W)).real)
    mixture_variance=float(np.sum(abs(H@W)**2).real-mixture_energy**2)
    assert abs(mixture_energy-(Q/2+q['rewiring_energy']))<2e-12
    assert abs(mixture_variance-(2*Q+q['mixed_rewiring_variance']))<2e-12
    mixedcoef={name:coefficients(H,W,diagstats[j]) for j,name in enumerate(names[:3])}
    assert all(abs(a)<2e-12 for v in mixedcoef.values() for a in v[1:3])
    initial=np.column_stack([*ps,W]);energies=np.sum(initial.conj()*(H@initial),axis=0).real
    results=[];maxnorm=0.;maxenergy=0.
    for T in [0.,.4,2.,8.,20.]:
        V=initial if T==0 else expm_multiply(-1j*T*H,initial)
        norm=np.sum(abs(V)**2,axis=0);maxnorm=max(maxnorm,float(np.max(abs(norm-np.r_[np.ones(3),np.full(D,1/D)]))))
        ev=np.sum(V.conj()*(H@V),axis=0).real;maxenergy=max(maxenergy,float(np.max(abs(ev-energies))))
        probs=np.column_stack([abs(V[:,:3])**2,np.sum(abs(V[:,3:])**2,axis=1)])
        val=diagstats@probs
        row=dict(T=T,preparations={kind:dict(zip(names,map(float,val[:,j]))) for j,kind in enumerate(['random','concentrated','spread','dephased'])})
        results.append(row)
    assert maxnorm<1e-10 and maxenergy<1e-9
    out=dict(N=n,Q=Q,L=L,dimension=dim,mixture_energy=mixture_energy,mixture_variance=mixture_variance,
             pure_initial_coefficients=coeffs,mixture_initial_coefficients=mixedcoef,finite_time=results,
             max_norm_error=maxnorm,max_energy_error=maxenergy,seconds=time.monotonic()-start)
    print(json.dumps(dict(L=L,dimension=dim,seconds=out['seconds'],finite_time=results)),flush=True)
    return out


def main():
    rows=[run(L) for L in [1,2,3,4,5]]
    comparisons=[]
    for a,b in zip(rows,rows[1:]):
        deltas=[]
        for ar,br in zip(a['finite_time'],b['finite_time']):
            differences={k:{name:abs(ar['preparations'][k][name]-br['preparations'][k][name]) for name in ['wedges','triangles','squares','vertex_fermions','electric_energy']} for k in ar['preparations']}
            deltas.append(dict(T=ar['T'],absolute_differences=differences))
        comparisons.append(dict(cutoffs=[a['L'],b['L']],rows=deltas))
    folder=Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))
    index=[]
    for row in rows:
        name=f"quantum_L{row['L']}.json"
        (folder/name).write_text(json.dumps(row,indent=2)+'\n')
        index.append(dict(file=name,L=row['L'],dimension=row['dimension'],max_norm_error=row['max_norm_error'],max_energy_error=row['max_energy_error']))
    out=dict(status='CP12_EQUAL_ENERGY_FULL_LAW_PROPAGATION_PASS',models=index,cutoff_comparisons=comparisons,
             scope='Complete finite-cutoff CP12 propagation in auxiliary N=4,Q=4. Mean energy matched; variance differs between pure and mixed preparations. No higher-density finite-time or infinite-cutoff certification.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'quantum_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':main()
