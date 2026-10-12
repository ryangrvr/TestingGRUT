#!/usr/bin/env python3
"""Independent CP11 fixed-bond entropy and positive swap-kinetic controls.

A declared counterfactual K-conserving model, NOT the frozen CP11 law.
Requires numpy/scipy. Outputs results.json. No GRUT modules imported.
"""
import itertools, json, math
from pathlib import Path
import numpy as np
from scipy.special import gammaln

def logchoose(n,k):
    return float(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1))

def expected_cycles(n,k,ell=4):
    m=n*(n-1)//2
    if min(n,k)<ell: return 0.0
    lf=lambda x: gammaln(x+1)-gammaln(x-ell+1)
    return float(np.exp(lf(n)+lf(k)-lf(m)-math.log(2*ell)))

def fermion_extrema(n,chosen):
    """Exact canonical M-fermion spectral min/max from one-body Dirac D."""
    pairs=list(itertools.combinations(range(n),2)); M=len(pairs)
    B=np.zeros((n,M),complex)
    for j in chosen:
        u,v=pairs[j]; B[u,j]=1; B[v,j]=-1
    D=np.block([[np.eye(n),B],[B.conj().T,-np.eye(M)]])
    E=np.linalg.eigvalsh(D)
    assert abs(np.trace(B@B.conj().T).real-2*len(chosen))<1e-10
    return float(sum(E[:M])+M),float(sum(E[-M:])+M)

def fermion_controls(n,k):
    m=n*(n-1)//2; worst=0
    for chosen in itertools.combinations(range(m),k):
        lo,hi=fermion_extrema(n,chosen)
        worst=max(worst,abs(lo),abs(hi))
    bound=2*n+math.sqrt(2*n*k)
    assert worst<=bound+1e-10
    return dict(N=n,K=k,graphs=math.comb(m,k),max_abs_fermion_energy=worst,
                exact_norm_upper_bound=bound)

def swap_controls(n,k,Gamma=0.1):
    """Conserves K: positive Laplacian on all K-edge configurations."""
    m=n*(n-1)//2
    graphs=list(itertools.combinations(range(m),k))
    lookup={g:i for i,g in enumerate(graphs)}
    L=np.zeros((len(graphs),len(graphs)))
    for i,g in enumerate(graphs):
        S=set(g)
        for x in g:
            for y in range(m):
                if y not in S:
                    h=tuple(sorted((S-{x})|{y}))
                    L[i,i]+=1;L[i,lookup[h]]-=1
    eig=np.linalg.eigvalsh(L)
    degree=k*(m-k)
    assert abs(eig.min())<1e-10
    assert np.linalg.norm(L@np.ones(len(graphs)))<1e-9
    assert eig.max()<=2*degree+1e-9
    return dict(N=n,K=k,dim=len(graphs),degree=degree,
                kinetic_norm=float(Gamma*eig.max()/n**2),
                proven_norm_upper_bound=float(2*Gamma*degree/n**2))

def count_controls():
    rows=[]
    for n in (16,25,36,64,100,200,500,1000):
        k=2*n; m=n*(n-1)//2
        log_bound=float(gammaln(n+1)-logchoose(m,k))
        C=2+math.sqrt(2*k/n)+0.2*k*(m-k)/n**3
        rows.append(dict(N=n,K=k,
            log10_uniform_template_probability_upper_bound=log_bound/math.log(10),
            uniform_expected_four_cycles=expected_cycles(n,k),
            operator_norm_coefficient_C=C,
            sufficient_beta_scale_from_bound=-log_bound/(2*C*n)))
    assert rows[-1]['log10_uniform_template_probability_upper_bound']<-3000
    return rows

def main():
    n,k=4,4
    edges=list(itertools.combinations(range(n),2))
    successes=0
    for idx in itertools.combinations(range(len(edges)),k):
        selected=[edges[j] for j in idx]
        if all(sum(v in e for e in selected)==2 for v in range(n)):
            successes+=1
    p=successes/math.comb(len(edges),k)
    assert abs(p-expected_cycles(n,k))<1e-12
    result=dict(
        theorem="Quantum fixed-inventory Gibbs class-probability bound",
        finite_bound="Pr(A)<=|A|/choose(N(N-1)/2,K)*exp(2 beta C N)",
        asymptotic="log Pr(A)<=-(alpha-1) N log N+O_beta(N)",
        conditions="K=alpha N+O(1), alpha>1, uniform graph-sector dimension, ||H-Eref||<=CN, beta fixed, |A|<=N! exp(O(N))",
        c4_probability_n4_k4=p,
        canonical_fermion=[fermion_controls(4,4),fermion_controls(5,6)],
        conserved_graph_swap=[swap_controls(4,4),swap_controls(5,6)],
        counting=count_controls(),
        cautions=[
            "No assumption that CP11 already conserves edge inventory.",
            "No restriction on noncommutation of quantum graph and matter terms.",
            "Not a no-go against zero-temperature geometry selection.",
            "Unlabelled gauge-state counting and physical clock/reference sectors need separate definition.",
            "Norm and counting are mathematically derived; originality and real-world physics unverified."
        ])
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('CP11_FIXED_BOND_ENTROPY_THEOREM_CHECK_PASS')

if __name__=="__main__": main()
